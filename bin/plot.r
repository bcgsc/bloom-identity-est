#!/usr/bin/env r

PROG_NAME = 'plot.r'

#------------------------------------------------------------
# parse command line args
#------------------------------------------------------------

if (length(argv) != 1) {
	cat(sprintf('Usage: %s <data.txt>\n', PROG_NAME))
	q(status=1)
}

data_file=argv[1]

#------------------------------------------------------------
# input data
#------------------------------------------------------------

data=read.delim(data_file, header=TRUE)

#------------------------------------------------------------
# output file
#------------------------------------------------------------

pdf('plot.pdf', useDingbats=FALSE)

#------------------------------------------------------------
# layout for multi-part plot
#------------------------------------------------------------

layout_matrix=matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
layout(layout_matrix, heights=c(1,1))

#------------------------------------------------------------
# plot: X (overlapping kmers)
#------------------------------------------------------------

hist(data$intersect_kmers_est, freq=FALSE,
	main='X\n(|G1 intersect G2|)', xlab='kmers',
	col='grey')

#------------------------------------------------------------
# plot: Y (kmers in genome)
#------------------------------------------------------------

hist(data$min_kmers_est, freq=FALSE,
	main='Y\n(min(|G1|,|G2|)', xlab='kmers',
	col='grey')

#------------------------------------------------------------
# plot: Z (percent identity)
#------------------------------------------------------------

hist(data$percent_id_est, freq=FALSE,
	main='Z=(X/Y)^(1/k)\n(percent_identity(G1,G2))', xlab='kmers',
	col='grey')

abline(v=0.80, lwd=2, lty=2)

#------------------------------------------------------------
# legend
#------------------------------------------------------------

savemar=par()$mar
par(mar=c(1,1,1,1))
legend_text = c(
	'kmer size: 13 bp',
	'simulated genome sizes: 100,000 bp',
	'true percent seq identity: 0.8',
	'Bloom filter false positive rate: 0.1'
)
plot(1, type = 'n', axes=FALSE, xlab='', ylab='')
legend('center', title='Parameters', pch=19, legend=legend_text)

#------------------------------------------------------------
# write plot output file
#------------------------------------------------------------

dev.off()
