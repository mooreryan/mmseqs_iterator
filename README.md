# iterate_mmesqs

Run mmseqs iteratively.

Use this if you want to pull out everything in the subject sequences even remotely similar to your query sequences.

After each iteration, any new hits will become the queries for the next iteration.

I will stop if I get past --max-iters or if I'm not giving you enough new hits (--min-percent-increase).
