#!/usr/bin/env r

PROG_NAME='calc-k'
DESCRIPTION = '
	Description: Given the size of two genomes, choose
	a kmer size that makes the probability of kmers
	occurring in both genomes by chance small.
'

#------------------------------------------------------------
# parse command line opts
#------------------------------------------------------------

library('getopt')
spec = matrix(c(
	'help', 'h', 0, 'logical', 'print this help message',
	'length', 'l', 1, 'double', 'genome length [required]'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (!is.null(opt$help) || is.null(opt$length)) {
	cat(strwrap(DESCRIPTION, simplify=TRUE))
	cat('\n')
	cat(getopt(spec, command=PROG_NAME, usage=TRUE))
	q(status=1)
}

#------------------------------------------------------------
# do the calculation
#------------------------------------------------------------

l = opt$length

# NOTE: The probability that a given kmer will be observed
# once or more in both genomes by chance is:
#
#   (1 - (1 - 1/4^k)^l1) * (1 - (1 - 1/4^k)^l2) 
#
# where l1 and l2 are the sizes of the two genomes.
#
# It follows that the expected number of kmers observed by
# chance in both genomes is:
#
#  min(l1,l2) * (1 - (1 - 1/4^k)^l1) * (1 - (1 - 1/4^k)^l2)
#
# To simplify things, we assume that l1 == l2 (call this
# value l). Then the equation becomes:
#
#     l * (1 - (1 - 1/4^k)^l)^2
#
# We want to choose k such that this number is small, so 
# we set
# 
#     l * (1 - (1 - 1/4^K)^l)^2 < 1
#
# and solve for k:
#
#     k > log4(1/(1-(1-1/sqrt(l))))

k = logb(1/(1-(1-1/sqrt(l))^(1/l)),base=4)
cat(ceiling(k))
