# MMseqs2 Iterator

Run mmseqs iteratively.

Use this if you want to pull out everything in the subject sequences
even remotely similar to your query sequences.  After each
iteration, any new hits will become the queries for the next
iteration.

I will stop if I get past --pipeline-max-iterations or if I'm not
giving you enough new hits (as specified by --min-percent-increase).

## False positives

You can get a lot of "false positives" this way.  To avoid this, you
can filter sequences with [PASV](https://github.com/mooreryan/pasv)
after each round of searching.  To do this, you need to pass the
--use-pasv flag, and then set all theother --pasv-* options as well.

--pasv-good-file should be the basename of the sequences that you
want to keep.  E.g., 'NCEC_Yes' would keep the sequences that had
NCEC at the key positions and spanned the region of interest.

Specify key positions like this: --pasv-key-positions 437 439
... separated by a space (not by a comma).

## Running with Docker

Note the use of `-v $(pwd):$(pwd)`.  This means that the contents of the folder in which you are running the Docker run command will be available inside the running container. 

Also note that you must prefix all of the filenames you pass in to options with `$(pwd)`, i.e., `$(pwd)/test_files/subject.fa`.  This is because we mounted the current working directory in the docker container (`$(pwd)`), and we need to pass in a full path so it knows where those files are located.

Note that mounting the current working directory in this way allows the Docker container to write to that directory.

```
time docker run --rm -v $(pwd):$(pwd) \
  mooreryan/mmseqs_iterator:0.4.2 \
  mmseqs_iterator \
  -q $(pwd)/assets/ClassIandII_ref_subset5.fasta \
  -s $(pwd)/test_files/subject.fa \
  -o $(pwd)/output_dir \
  --mmseqs-threads 3 \
  --mmseqs-num-iterations 1 \
  --mmseqs-sensitivity 1 \
  --use-pasv \
  --pasv-good-file NCEC_Yes NCEC_No
  --pasv-key-positions 437 439 441 462
```