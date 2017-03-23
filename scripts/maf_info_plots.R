dat <- read.table("dat.snp-stats", header=TRUE)

pdf(file="maf.pdf")
hist(dat$MAF, breaks=100, main="Histogram of MAF", xlab="MAF")
dev.off()

pdf(file="info.pdf")
hist(dat$information, breaks=100, main="Histogram of info scores", xlab="Info scores")
dev.off()

library(ggplot2)
ggplot(dat, aes(x=MAF, y=information)) +
geom_smooth()
ggsave(file="maf_vs_info.pdf")
