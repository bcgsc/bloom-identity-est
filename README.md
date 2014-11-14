# APPLICATION/BACKGROUND

Rene has invented a new method for estimating percent sequence identity between two genomes, using a probabilistic data structure called a Bloom filter.  Two Bloom filters are used to represent the set of kmers in each genome, and the idea of the method is to estimate the percent sequence identity based on the number of overlapping kmers.

The main assumption of the method is that single-nucleotide differences between the two genomes are randomly and independently distributed.  Under this assumption, the equation that relates percent sequence identity to number of overlapping kmers is:

X/Y = Z^k

where X is the number of overlapping kmers, Y is the number of kmers in the smaller of the two genomes, and Z is the percent sequence identity between the two genomes.

# THE PROBLEM

We have a random variable Z that is a function of two other random
variables X and Y:

Z = (X/Y)^(1/k)

where k is a constant.

We are given functions for:

1. MLE(X)
2. MLE(Y)
3. f\_x(d) such that p(MLE(X) - d <= X <= MLE(X) + d) >= f\_x(d)
4. f\_y(d) such that p(MLE(Y) - d <= Y <= MLE(Y) + d) >= f\_y(d)

Using those functions, is it possible to determine d such that:

p(MLE(Z) - d <= Z <= MLE(Z) + d) >= 0.95 ?

# EQUATIONS

MLE(X) is given by Equation (6) in attached bloom\_set\_ops.pdf (Section 5.2)
MLE(Y) is given by Equation (3) in attached bloom\_set\_ops.pdf (Section 5.1)
f\_x(d) is given by Theorem 2 in attached bloom\_set\_ops.pdf (Section 5.2)
f\_y(d) is given by Theorem 1 in attached bloom\_set\_ops.pdf (Section 5.1)

# IDEAS

I came across the idea of "Probability Bounds Analysis" which sounds like it might help, but thus far I don't understand the technique.

Wikipedia: [http://en.wikipedia.org/wiki/Probability\_bounds\_analysis](http://en.wikipedia.org/wiki/Probability\_bounds\_analysis)

Tutorial: attached pba\_tutorial.pdf
