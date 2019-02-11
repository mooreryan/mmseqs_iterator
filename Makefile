TEST_FILES = test_files
TEST_OUTDIR = test_outdir

QUERIES = $(TEST_FILES)/queries.fa
SUBJECT = $(TEST_FILES)/subject.fa
THREADS = 3

ITERATE_MMSEQS = iterate_mmseqs.rb

all: test_iterate_mmseqs

test_iterate_mmseqs:
	rm -r $(TEST_OUTDIR); ./$(ITERATE_MMSEQS) -q $(QUERIES) -s $(SUBJECT) -o $(TEST_OUTDIR) -t $(THREADS) --num-iters 1 --sensitivity 1
