#!/usr/bin/env ruby

# Copyright 2018 Ryan Moore
# Contact: moorer@udel.edu
#
# This file is part of iterate_mmesqs.
#
# iterate_mmseqs is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# iterate_mmseqs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with iterate_mmseqs.  If not, see
# <http://www.gnu.org/licenses/>.

# def easy_search! queries, subjects, outdir, basename
#   tmpdir = File.join outdir, "tmp"
#   FileUtils.rm_r tmpdir if File.exist?(tmpdir)

#   outf = "#{File.join outdir, basename}.txt"
#   FileUtils.rm outf if File.exist?(outf)

#   cmd = "#{MMSEQS} easy-search " \
#         "#{queries} " \
#         "#{subjects} " \
#         "#{outf} " \
#         "#{tmpdir} " \
#         "--num-iterations 1 " \
#         "-s 1 " \
#         "--format-mode 2 " \
#         "--threads #{THREADS} " \
#         ">> #{log} 2>&1"

#   Process.run_and_time_it! "Running search", cmd

#   outf
# end

# def iterate_easy_search! queries, subjects, outdir, basename
#   iter = 0
#   new_queries = queries
#   all_hit_names = Set.new


#   loop do
#     iter += 1

#     btab = easy_search! new_queries, subjects, outdir, basename

#     new_query_info = make_new_queries queries,
#                                       subjects,
#                                       btab,
#                                       outdir,
#                                       basename,
#                                       all_hit_names,
#                                       iter

#     new_queries = new_query_info[:new_queries]

#     increase =
#       new_query_info[:new_hit_count] /
#       new_query_info[:all_hits_count].to_f

#     Rya::AbortIf.logger.info do
#       "Iter: #{iter}, " \
#       "new hits: #{new_query_info[:new_hit_count]}, " \
#       "all hits: #{new_query_info[:all_hits_count]}, " \
#       "increase: #{(increase * 100).round(2)}%"
#     end

#     if iter > MAX_ITERS ||
#        new_query_info[:new_hit_count].zero? ||
#        increase <= STOP

#       break
#     end
#   end
# end

Signal.trap("PIPE", "EXIT")

require "rya"
require "fileutils"
require "optimist"
require "tempfile"
require "set"

Process.extend Rya::CoreExtensions::Process
include Rya::AbortIf

RNR_REFS = File.join __dir__, "assets", "ClassIandII_ref_subset5.fasta"

GOOD_HITS_FILES = Set.new

DB_SUFFIX = ".mmseqs_db"

VERSION   = "v0.3.0"
COPYRIGHT = "2018 - 2019 Ryan Moore"
CONTACT   = "moorer@udel.edu"
WEBSITE   = "https://github.com/mooreryan/iterate_mmseqs"
LICENSE   = "GPLv3"


VERSION_BANNER = "# Version:   #{VERSION}
# Copyright: #{COPYRIGHT}
# Contact:   #{CONTACT}
# License:   #{LICENSE}"

def pasv!(exe:,
          outdir:,
          refs:,
          queries:,
          start:,
          stop:,
          threads:,
          posns:)

  cmd = "#{exe} " \
  "--outdir #{outdir} " \
  "--aligner clustalo " \
  "--alignment-parameters '\\--threads 1' " \
  "--refs #{refs} " \
  "--queries #{queries} " \
  "--start #{start} " \
  "--end #{stop} " \
  "--threads #{threads} " \
  "#{posns.join(" ")}"

  Process.run_and_time_it! "Running PASV", cmd
end

def get_good_pasv_seqs_path pasv_outdir, good_file
  # pasv_outdir/
  # -- pasv.partition_CF_Yes.fa
  # -- pasv.partition_ED_No.fa
  # -- pasv.partition_ED_Yes.fa
  # -- pasv_counts.txt
  # -- pasv_types.txt

  good_file_path = File.join pasv_outdir,
                             "pasv.partition_#{good_file}.fa"

  # we only want to send the path if the file exists and actually has sequences.
  if File.exist?(good_file_path) && count_seqs(good_file_path) > 0
    good_file_path
  else
    nil
  end
end

def create_db! seqs, outdir
  basename = File.basename seqs
  db_name  = File.join outdir, "#{basename}#{DB_SUFFIX}"

  cmd = "#{MMSEQS} createdb #{seqs} #{db_name} >> #{MMSEQS_LOG} 2>&1"
  Process.run_and_time_it! "Making DB", cmd

  db_name
end

def index_db! db, outdir
  tmpdir = File.join outdir, "tmp_index"
  FileUtils.rm_r tmpdir if File.exist? tmpdir

  cmd = "#{MMSEQS} createindex " \
        "#{db} " \
        "#{tmpdir} " \
        "--include-headers " \
        "--threads #{THREADS} " \
        ">> #{MMSEQS_LOG} 2>&1"

  Process.run_and_time_it! "Making DB", cmd
end

def count_seqs seqs
  cmd = "grep -c '^>' #{seqs}"

  `#{cmd}`.chomp.to_i
end

def seq_names fasta
  Set.new `grep '^>' #{fasta} | sed 's/^>//'`.
    chomp.
    split("\n").
    map { |header| header.split(" ").first }
end

# @note that if there are duplicates in the queries and the subjects,
#   they will not be treated any differently than if they were not
#   duplicates.
#
# @note all_hit_names will be modified
def make_new_queries(subjects:,
                     btab:,
                     outdir:,
                     basename:,
                     all_hit_names:,
                     iter:,
                     pasv_use: nil,
                     pasv_exe: nil,
                     pasv_refs: nil,
                     pasv_start: nil,
                     pasv_stop: nil,
                     pasv_threads: nil,
                     pasv_posns: nil,
                     pasv_good_file: nil)
  new_queries_fname =
    "#{File.join(outdir, basename)}.new_queries_iter_#{iter}.faa"
  FileUtils.rm new_queries_fname if File.exist?(new_queries_fname)

  new_subjects_fname =
    "#{File.join(outdir, basename)}.new_subjects_iter_#{iter}.faa"
  FileUtils.rm new_subjects_fname if File.exist?(new_subjects_fname)

  new_hits = nil
  Tempfile.open do |ids_f|
    cmd = "cut -f2 #{btab} > #{ids_f.path}"

    Process.run_and_time_it! "Getting ids from homologous subject seqs",
                             cmd

    # ids_f is the file with subject hits for the current iteration.
    #
    # This will remove any seq from the subject DB that was a hit this round.
    #
    # TODO Note that the new subjects will have all hits from this round removed, even those that did not pass the PASV filter if it is run.  Thus, there may be a weird thing in the logs re. the total hits and database size not adding up quite right
    Process.run_and_time_it! "Making new subject seqs for next round",
                             "#{ANTI_GREP_IDS} " \
                             "#{ids_f.path} " \
                             "#{subjects} " \
                             "> #{new_subjects_fname}"


    Tempfile.open do |seqs_f|
      cmd = "#{GREP_IDS} " \
            "#{ids_f.path} " \
            "#{subjects} " \
            "> #{seqs_f.path}"

      Process.run_and_time_it! "Pulling seqs from btab", cmd

      # Get the hit names for this round of searching.  Including seqs that may later fail PASV.
      hit_names = seq_names seqs_f.path

      # Get a list of new hits
      new_hits = hit_names - all_hit_names

      Rya::AbortIf.logger.info do
        "There were #{new_hits.count} total new hits"
      end

      # Add in the new hits to all hits for the next possible
      # iteration.
      unless USE_PASV
        new_hits.each do |hit|
          all_hit_names << hit
        end
      end


      Tempfile.open do |new_names_f|
        new_names_f.puts new_hits.to_a

        new_names_f.fsync

        Process.run_and_time_it! "Grepping new hits for next iteration",
                                 "#{GREP_IDS} " \
                                 "#{new_names_f.path} " \
                                 "#{seqs_f.path} " \
                                 "> #{new_queries_fname}"

      end
    end
  end

  # Now that I have the new queries and the new subjects.  I need to run the new_queries through pasv if that is what the user asked for.
  new_hit_count = new_hits.count
  if pasv_use
    pasv_outdir = File.join outdir, "PASV_iter_#{iter}"
    # Then do pasv!
    # TODO somewhere we need to clean up the PASV output.
    #
    # We want to run this on the file currently marked as the new queries file so that we filter any of these through PASV.
    pasv!(exe:     pasv_exe,
          outdir:  pasv_outdir,
          refs:    pasv_refs,
          queries: new_queries_fname,
          start:   pasv_start,
          stop:    pasv_stop,
          threads: pasv_threads,
          posns:   pasv_posns)

    # And pull the good seqs
    good_pasv_seqs_fname = get_good_pasv_seqs_path(pasv_outdir, pasv_good_file)

    Rya::AbortIf.logger.debug do
      "good_pasv_seqs_fname: #{good_pasv_seqs_fname}"
    end

    if good_pasv_seqs_fname
      # Then we have new queries that passed PASV filtering.
      new_queries_fname = good_pasv_seqs_fname

      # We need to adjust the numbers of hits and new hits to reflect what passed PASV.
      hit_names = seq_names new_queries_fname

      # Get a list of new hits
      new_hits = hit_names - all_hit_names
      new_hit_count = new_hits.count

      Rya::AbortIf.logger.info do
        "There were #{new_hit_count} hits that also made it through PASV filtering"
      end

      # Add in the new hits to all hits for the next possible
      # iteration.
      new_hits.each do |hit|
        all_hit_names << hit
      end
    else
      new_queries_fname = nil
      new_hit_count     = nil
    end
  end

  { new_queries:    new_queries_fname,
    new_subjects:   new_subjects_fname,
    new_hit_count:  new_hit_count,
    all_hits_count: all_hit_names.count }
end

def search! queries, subject_db, outdir, basename, iter
  tmpdir = File.join outdir, "tmp"
  FileUtils.rm_r tmpdir if File.exist?(tmpdir)

  out_db = "#{File.join outdir, basename}.iter_#{iter}#{DB_SUFFIX}"
  outf   = out_db + ".btab.txt"
  FileUtils.rm outf if File.exist?(out_db)
  FileUtils.rm outf if File.exist?(outf)

  # First you need to convert queries to a query db
  query_db = create_db! queries, outdir

  cmd = "#{MMSEQS} search " \
        "#{query_db} " \
        "#{subject_db} " \
        "#{out_db} " \
        "#{tmpdir} " \
        "--num-iterations #{NUM_ITERS} " \
        "-s #{SENS} " \
        "--threads #{THREADS} " \
        ">> #{MMSEQS_LOG} 2>&1"

  Process.run_and_time_it! "Running search", cmd

  # Then you have to convert to btab
  cmd = "#{MMSEQS} convertalis " \
        "#{query_db} " \
        "#{subject_db} " \
        "#{out_db} " \
        "#{outf} " \
        "--threads #{THREADS} " \
        "--format-mode 2 " \
        ">> #{MMSEQS_LOG} 2>&1"

  Process.run_and_time_it! "Running convertalis", cmd

  # Then remove the query_db
  FileUtils.rm query_db

  outf
end

def iterate_search! queries, subjects, outdir, basename
  iter          = 0
  new_queries   = queries
  all_hit_names = Set.new

  # First make and index the subject database
  new_subject_db = create_db! subjects, outdir
  index_db! new_subject_db, outdir

  new_subjects = subjects

  loop do
    iter += 1

    Rya::AbortIf.logger.info do
      "Just started iteration #{iter}"
    end

    # Try and remove last iterations DBs as you go to save space.
    # Minus 2 rather than minus 1 is correct.
    # TODO this still doesn't get all the files that could go after each iteration.
    previous_dbs = Dir.glob File.join outdir, "*iter_#{iter - 2}*mmseqs_db*"
    previous_dbs.each do |fname|
      unless fname.include? ".btab.txt"
        FileUtils.rm fname
      end
    end

    btab = search! new_queries, new_subject_db, outdir, basename, iter

    new_query_info = make_new_queries(subjects:       new_subjects,
                                      btab:           btab,
                                      outdir:         outdir,
                                      basename:       basename,
                                      all_hit_names:  all_hit_names,
                                      iter:           iter,
                                      pasv_use:       USE_PASV,
                                      pasv_exe:       PASV_EXE,
                                      pasv_refs:      PASV_REFS,
                                      pasv_start:     PASV_ROI_START,
                                      pasv_stop:      PASV_ROI_END,
                                      pasv_threads:   THREADS,
                                      pasv_posns:     PASV_KEY_POSITIONS,
                                      pasv_good_file: PASV_GOOD_FILE)

    new_queries = new_query_info[:new_queries]
    GOOD_HITS_FILES << new_queries
    new_subjects = new_query_info[:new_subjects]

    if new_queries.nil?
      # Then were using PASV but got no new seqs.
      Rya::AbortIf.logger.info do
        "No seqs passed PASV filtering Iter: #{iter}.  Stopping iteration."
      end

      break
    end

    # Make the new subject DB (it won't have any sequence already hit)
    # TODO this only speeds things up once you get pretty far into the
    # test and you're collecting most of the things in the database.
    new_subject_db = create_db! new_subjects, outdir

    # TODO if using PASV this might not be quite right as we haven't updated all_hits_count to reflect only the hits that passed PASV.
    increase =
      new_query_info[:new_hit_count] /
        new_query_info[:all_hits_count].to_f

    Rya::AbortIf.logger.info do
      "Iter: #{iter}, " \
      "new hits: #{new_query_info[:new_hit_count]}, " \
      "all hits: #{new_query_info[:all_hits_count]}, " \
      "increase: #{(increase * 100).round(2)}%, " \
      "new db size: #{count_seqs new_subjects}" # TODO this could potentially be pretty slow if the db is big enough.
    end

    # TODO pretty sure this will do one extra iter than user asked for
    if iter > MAX_ITERS ||
      new_query_info[:new_hit_count].zero? ||
      increase <= STOP

      break
    end
  end

  # Remove the rest of the DBs
  previous_dbs = Dir.glob File.join outdir, "*mmseqs_db*"
  previous_dbs.each do |fname|
    unless fname.include? ".btab.txt"
      FileUtils.rm fname
    end
  end
end

opts = Optimist.options do
  version VERSION_BANNER

  banner <<-EOS

#{VERSION_BANNER}

  Run mmseqs iteratively.

  Use this if you want to pull out everything in the subject sequences
  even remotely similar to your query sequences.

  After each iteration, any new hits will become the queries for the
  next iteration.

  I will stop if I get past --max-iters or if I'm not giving you
  enough new hits (--min-percent-increase).

  Options:
  EOS

  opt(:queries, "Queries", type: :string)
  opt(:subject, "Subject", type: :string)
  opt(:basename, "Basename of outfiles", default: "search")
  opt(:outdir, "Outdir", default: ".")

  opt(:threads, "Number of cores (mmseqs option)", default: 1)
  opt(:num_iters, "Number of iterations (mmseqs option)", default: 2)
  opt(:sensitivity, "Sensitivity (mmseqs option)", default: 5.7)

  opt(:max_iters, "Max number of iterations", default: 10)
  opt(:min_percent_increase,
      "Minimum percent increase in new hits to continue",
      default: 10)

  opt(:mmseqs, "/path/to/mmseqs", default: "mmseqs")
  opt(:grep_ids, "/path/to/grep_ids", default: "grep_ids")
  opt(:anti_grep_ids, "/path/to/anti_grep_ids", default: "anti_grep_ids")

  # PASV opts
  opt(:pasv_use, "Set this flag to use PASV to filter hits", default: false)
  opt(:pasv, "/path/to/pasv", default: "pasv")
  opt(:pasv_refs, "Fasta with refs", default: RNR_REFS)
  opt(:pasv_roi_start, "Start of ROI", default: 437)
  opt(:pasv_roi_end, "End of ROI", default: 625)
  opt(:pasv_key_positions, "List of key positions", default: [437, 439, 441, 462])
  opt(:pasv_good_file, "The file with seqs to keep", default: "NCEC_YES")
end

queries  = opts[:queries]
subjects = opts[:subject]
basename = opts[:basename]
outdir   = opts[:outdir]

# PASV opts
# TODO check these opts for good input
USE_PASV           = opts[:pasv_use]
PASV_EXE           = opts[:pasv]
PASV_REFS          = opts[:pasv_refs]
PASV_ROI_START     = opts[:pasv_roi_start]
PASV_ROI_END       = opts[:pasv_roi_end]
PASV_KEY_POSITIONS = opts[:pasv_key_positions]
PASV_GOOD_FILE     = opts[:pasv_good_file]

# MMseqs opts
THREADS   = opts[:threads]
NUM_ITERS = opts[:num_iters]
SENS      = opts[:sensitivity]

MAX_ITERS = opts[:max_iters]
STOP      = opts[:min_percent_increase] / 100.0

MMSEQS        = opts[:mmseqs]
GREP_IDS      = opts[:grep_ids]
ANTI_GREP_IDS = opts[:anti_grep_ids]

MMSEQS_LOG = File.join outdir, "mmseqs_log.txt"

abort_unless queries && File.exist?(queries), "#{queries} does not exist."
abort_unless subjects && File.exist?(subjects), "#{subjects} does not exist."

FileUtils.mkdir_p outdir

work_dir = File.join outdir, "work"
FileUtils.mkdir_p work_dir

final_outdir = File.join outdir, "hits"
FileUtils.mkdir_p final_outdir

iterate_search! queries,
                subjects,
                work_dir,
                basename

# hits_glob = File.join work_dir, "*.new_queries_iter_*.faa"

all_hits = File.join final_outdir, "#{basename}.hits.faa"
Process.run_and_time_it! "Making hits files",
                         "cat #{GOOD_HITS_FILES.to_a.join(" ")} > #{all_hits}"

Rya::AbortIf.logger.info { "Final output: #{all_hits}" }
