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
require "trollop"
require "tempfile"
require "set"

Process.extend Rya::CoreExtensions::Process
include Rya::AbortIf

DB_SUFFIX = ".mmseqs_db"

VERSION   = "v0.1.0"
COPYRIGHT = "2018 Ryan Moore"
CONTACT   = "moorer@udel.edu"
WEBSITE   = "https://github.com/mooreryan/iterate_mmseqs"
LICENSE   = "GPLv3"


VERSION_BANNER = "# Version:   #{VERSION}
# Copyright: #{COPYRIGHT}
# Contact:   #{CONTACT}
# License:   #{LICENSE}"

def create_db! seqs, outdir
  basename = File.basename seqs
  db_name = File.join outdir, "#{basename}#{DB_SUFFIX}"

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
def make_new_queries queries,
                     subjects,
                     btab,
                     outdir,
                     basename,
                     all_hit_names,
                     iter
  new_seqs_fname =
    "#{File.join(outdir, basename)}.new_queries_iter_#{iter}.faa"
  FileUtils.rm new_seqs_fname if File.exist?(new_seqs_fname)

  new_hits = nil
  Tempfile.open do |ids_f|
    cmd = "cut -f2 #{btab} > #{ids_f.path}"

    Process.run_and_time_it! "Getting ids", cmd

    Tempfile.open do |seqs_f|
      cmd = "#{GREP_IDS} " \
            "#{ids_f.path} " \
            "#{subjects} " \
            "> #{seqs_f.path}"

      Process.run_and_time_it! "Pulling seqs from btab", cmd

      # Get the hit names for this round of searching.
      hit_names = seq_names seqs_f.path

      # Get a list of new hits
      new_hits = hit_names - all_hit_names

      # Add in the new hits to all hits for the next possible
      # iteration.
      new_hits.each do |hit|
        all_hit_names << hit
      end

      Tempfile.open do |new_names_f|
        new_names_f.puts new_hits.to_a

        new_names_f.fsync

        Process.run_and_time_it! "Grepping new hits for next iteration",
                                 "#{GREP_IDS} " \
                                 "#{new_names_f.path} " \
                                 "#{seqs_f.path} " \
                                 "> #{new_seqs_fname}"

      end
    end
  end

  { new_queries: new_seqs_fname,
    new_hit_count: new_hits.count,
    all_hits_count: all_hit_names.count }
end

def search! queries, subject_db, outdir, basename, iter
  tmpdir = File.join outdir, "tmp"
  FileUtils.rm_r tmpdir if File.exist?(tmpdir)

  out_db = "#{File.join outdir, basename}.iter_#{iter}#{DB_SUFFIX}"
  outf = out_db + ".btab.txt"
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
  iter = 0
  new_queries = queries
  all_hit_names = Set.new

  # First make and index the subject database
  subject_db = create_db! subjects, outdir
  index_db! subject_db, outdir

  loop do
    iter += 1

    btab = search! new_queries, subject_db, outdir, basename, iter

    new_query_info = make_new_queries queries,
                                      subjects,
                                      btab,
                                      outdir,
                                      basename,
                                      all_hit_names,
                                      iter

    new_queries = new_query_info[:new_queries]

    increase =
      new_query_info[:new_hit_count] /
      new_query_info[:all_hits_count].to_f

    Rya::AbortIf.logger.info do
      "Iter: #{iter}, " \
      "new hits: #{new_query_info[:new_hit_count]}, " \
      "all hits: #{new_query_info[:all_hits_count]}, " \
      "increase: #{(increase * 100).round(2)}%"
    end

    if iter > MAX_ITERS ||
       new_query_info[:new_hit_count].zero? ||
       increase <= STOP

      break
    end
  end
end

opts = Trollop.options do
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

end

queries = opts[:queries]
subjects = opts[:subject]
basename = opts[:basename]
outdir = opts[:outdir]

# MMseqs opts
THREADS = opts[:threads]
NUM_ITERS = opts[:num_iters]
SENS = opts[:sensitivity]

MAX_ITERS = opts[:max_iters]
STOP = opts[:min_percent_increase] / 100.0

MMSEQS = opts[:mmseqs]
GREP_IDS = opts[:grep_ids]

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

hits_glob = File.join work_dir, "*.new_queries_iter_*.faa"

all_hits = File.join final_outdir, "#{basename}.hits.faa"
Process.run_and_time_it! "Making hits files",
                         "cat #{hits_glob} > #{all_hits}"

Rya::AbortIf.logger.info { "Final output: #{all_hits}" }
