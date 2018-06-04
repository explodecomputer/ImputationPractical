dat <- read.table("results/data_chr16.snp-stats.gz", header=TRUE)

pdf(file="output/maf.pdf")
hist(dat$MAF, breaks=100, main="Histogram of MAF", xlab="MAF")
dev.off()

pdf(file="output/info.pdf")
hist(dat$information, breaks=100, main="Histogram of info scores", xlab="Info scores")
dev.off()

library(ggplot2)
p <- ggplot(dat, aes(x=MAF, y=information)) +
geom_smooth()
ggsave(p, file="output/maf_vs_info.pdf")
