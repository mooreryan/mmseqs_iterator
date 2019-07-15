TEST_FILES = test_files
TEST_OUTDIR = test_outdir

QUERIES = $(TEST_FILES)/queries.fa
SUBJECT = $(TEST_FILES)/subject.fa
THREADS = 3

MMSEQS_ITERATOR = mmseqs_iterator.rb

test_all: test_mmseqs_iterator test_rnr test_rnr_no_pasv test_random_seqs_with_rnr_filter

test_mmseqs_iterator:
	rm -r $(TEST_OUTDIR); ./$(MMSEQS_ITERATOR) -q $(QUERIES) -s $(SUBJECT) -o $(TEST_OUTDIR) --mmseqs-threads $(THREADS) --mmseqs-num-iterations 1 --mmseqs-sensitivity 1 --pipeline-max-iterations 3

test_rnr:
	rm -r $(TEST_OUTDIR).rnr; ./$(MMSEQS_ITERATOR) -q assets/ClassIandII_ref_subset5.fasta -s $(SUBJECT) -o $(TEST_OUTDIR).rnr --mmseqs-threads $(THREADS) --mmseqs-num-iterations 1 --mmseqs-sensitivity 1 --use-pasv --pasv /Users/moorer/projects/pasv/pasv

test_rnr_multi_good_files:
	rm -r $(TEST_OUTDIR).rnr_multi_good_file; ./$(MMSEQS_ITERATOR) -q assets/ClassIandII_ref_subset5.fasta -s $(SUBJECT) -o $(TEST_OUTDIR).rnr_multi_good_file --mmseqs-threads $(THREADS) --mmseqs-num-iterations 1 --mmseqs-sensitivity 1 --use-pasv --pasv /Users/moorer/projects/pasv/pasv --pasv-good-file NCEC_Yes NCEC_No

test_rnr_no_pasv:
	rm -r $(TEST_OUTDIR).rnr_no_pasv; ./$(MMSEQS_ITERATOR) -q assets/ClassIandII_ref_subset5.fasta -s $(SUBJECT) -o $(TEST_OUTDIR).rnr_no_pasv --mmseqs-threads $(THREADS) --mmseqs-num-iterations 1 --mmseqs-sensitivity 1

test_random_seqs_with_rnr_filter:
	rm -r $(TEST_OUTDIR).random_with_rnr_filter; ./$(MMSEQS_ITERATOR) -q $(TEST_FILES)/small.fa -s $(SUBJECT) -o $(TEST_OUTDIR).random_with_rnr_filter --mmseqs-threads $(THREADS) --mmseqs-num-iterations 1 --mmseqs-sensitivity 1 --use-pasv --pasv /Users/moorer/projects/pasv/pasv
