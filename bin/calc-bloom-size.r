#!/usr/bin/env r

PROG_NAME='calc-bloom-size'
DESCRIPTION = '
	Description: Estimate Bloom filter size in bytes
	for a given number of hash functions, inserted
	elements, and target false positive rate.
'

#------------------------------------------------------------
# parse command line opts
#------------------------------------------------------------

library('getopt')
spec = matrix(c(
	'help', 'h', 0, 'logical', 'print this help message',
	'size', 'n', 1, 'integer', 'number of elements [required]',
	'fpr', 'f', 1, 'double', 'desired false positive rate [required]',
	'num_hash', 'N', 0, 'integer', 'number of hash functions [1]'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (!is.null(opt$help) || is.null(opt$size) || is.null(opt$fpr)) {
	cat(strwrap(DESCRIPTION, simplify=TRUE))
	cat('\n')
	cat(getopt(spec, command=PROG_NAME, usage=TRUE))
	q(status=1)
}

# default value
if (is.null(opt$num_hash)) { opt$num_hash = 1 }

#------------------------------------------------------------
# do the calculation
#------------------------------------------------------------

n = opt$size
fpr = opt$fpr
h = opt$num_hash

# NOTE: The equation for estimating false postive rate is:
#
#    fpr = (1 - e^(hn/m))^h
#
# Rearranging this equation for m (size in bits) gives us:
#
#    m = hn/ln(1-fpr^(1/h))
#
# We divide by 8 to convert the size from bits to bytes.

bytes = -h*n/(8*log(1-fpr^(1/h)))
cat(ceiling(bytes))
