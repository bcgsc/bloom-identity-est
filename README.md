# Description

These scripts provide a fast, memory-efficient method for estimating the percent sequence identity between two genomes using a probabilistic data structure called a Bloom filter ([wikipedia/Bloom_filter](https://en.wikipedia.org/wiki/Bloom_filter)).

# Method

1. The sequence(s) of the two genomes are cut up into k-mers and loaded into two separate Bloom filters using `abyss-bloom build`.
2. The number of overlapping k-mers between the two genomes is estimated from the bitwise intersection of the two Bloom filters using `abyss-bloom intersect`.
3. The percent sequence identity is estimated from the number of overlapping k-mers, according to the following formula.
   ```
   O/G = I^k
   ```
   where O is the number of overlapping kmers between the two genomes, G is the number of kmers in the smaller of the two genomes, I is the percent sequence identity between the two genomes, and k is the k-mer
   size.

The main assumption of the method is that single-nucleotide differences between the two genomes are randomly and independently distributed.

# Usage

The main for script for calculating the percent identity is `monte-carlo-experiment/real-genomes/monte-carlo.mk`:

```
$ monte-carlo.mk genome1=ecoli-strain1.fasta genome1_name=strain1 \
	genome2=ecoli-strain2.fasta genome2_name=strain2 \
	k=20 b=10G out=results.txt
```

# Parameters for `monte-carlo.mk`

* `genome1`: FASTA file for genome 1 [required]
* `genome1_name`: label for genome 1 (used to generate names of temp files) [required]
* `genome2`: FASTA file for genome 2 [required]
* `genome2_name`: label for genome 2 (used to generate names of temp files) [required]
* `k`: k-mer size [required]
* `b`: Bloom filter size. Allowable units are kilobytes ('k'), megabytes ('M'), gigabytes ('G') [required]
* `j`: number of threads [1]
* `s`: exclude genome sequences shorter than this length [0]

# Authors

* Rene Warren
* Ben Vandervalk - [GitHub/benvvalk](https://github.com/benvvalk)
