TEST_FILES = test_files
TEST_OUTDIR = test_outdir

QUERIES = $(TEST_FILES)/queries.fa
SUBJECT = $(TEST_FILES)/subject.fa
THREADS = 3

ITERATE_MMSEQS = iterate_mmseqs.rb

all: test_iterate_mmseqs

test_iterate_mmseqs:
	rm -r $(TEST_OUTDIR); ./$(ITERATE_MMSEQS) -q $(QUERIES) -s $(SUBJECT) -o $(TEST_OUTDIR) -t $(THREADS) --num-iters 1 --sensitivity 1 --max-iters 3

test_rnr:
	rm -r $(TEST_OUTDIR).rnr; ./$(ITERATE_MMSEQS) -q assets/ClassIandII_ref_subset5.fasta -s $(SUBJECT) -o $(TEST_OUTDIR).rnr -t $(THREADS) --num-iters 1 --sensitivity 1 --pasv-use --pasv /Users/moorer/projects/pasv/pasv

test_rnr_no_pasv:
	rm -r $(TEST_OUTDIR).rnr_no_pasv; ./$(ITERATE_MMSEQS) -q assets/ClassIandII_ref_subset5.fasta -s $(SUBJECT) -o $(TEST_OUTDIR).rnr_no_pasv -t $(THREADS) --num-iters 1 --sensitivity 1

test_random_seqs_with_rnr_filter:
	rm -r $(TEST_OUTDIR).random_with_rnr_filter; ./$(ITERATE_MMSEQS) -q $(TEST_FILES)/small.fa -s $(SUBJECT) -o $(TEST_OUTDIR).random_with_rnr_filter -t $(THREADS) --num-iters 1 --sensitivity 1 --pasv-use --pasv /Users/moorer/projects/pasv/pasv
