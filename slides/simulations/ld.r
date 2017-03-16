library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(reshape2)

n <- 1000
nsnp <- 10

r <- matrix(0.8, nsnp, nsnp)
diag(r) <- 1

x1 <- rmvnorm(n, rep(0, nsnp), sigma=r)
x2 <- rmvnorm(n, rep(0, nsnp), sigma=r)

x1[1:500, nsnp] <- x2[1:500, 1]

for(i in 1:nsnp)
{
	x2[,i] <- x2[,i] + rnorm(n, sd=i/(nsnp*2))
	x1[,i] <- x1[,i] + rnorm(n, sd=(nsnp-i+1)/(nsnp*2))
}

x <- cbind(x1, x2)


cv <- 10
g1 <- 5
g2 <- 15


make_plots <- function(x, cv, g1, g2, vg)
{
	melted_correlation_matrix <- melt(cor(x))
dot_dat <- data.frame(Var1=c(cv, g1, g2), Var2=c(cv, g1, g2), value=1, val=c("Causal variant", "Genotyped", "Genotyped"))
# dot_dat <- data.frame(Var1=nsnp, Var2=nsnp, value=1)
p1 <- ggplot(melted_correlation_matrix, aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() + 
	scale_fill_gradient(low="white", high="grey") +
	theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
	geom_point(data=dot_dat, aes(colour=val), size=3) +
	labs(x="SNP number", y="SNP number", fill="R-square", colour="") +
	scale_colour_manual(values=c("white", "black"))


	y <- x[,10] * sqrt(vg) + rnorm(n, sd=sqrt(1-vg))
	cor(y, x[,10])^2


mod1 <- summary(lm(y ~ x[,10]))
mod2 <- summary(lm(y ~ x[,5]))
mod3 <- summary(lm(y ~ x[,15]))
mod4 <- summary(lm(y ~ x[,5] + x[,15]))

pval <- rep(0, ncol(x))
for(i in 1:ncol(x))
{
	pval[i] <- summary(lm(y ~ x[,i]))$coefficients[2,4]
}
dat <- data.frame(pval=-log10(pval), snp=1:ncol(x), lab="Untyped SNP", stringsAsFactors=FALSE)
dat$lab[cv] <- "Causal variant"
dat$lab[g1] <- "Genotyped variant"
dat$lab[g2] <- "Genotyped variant"

p2 <- ggplot(dat, aes(x=snp, y=pval)) + geom_point(size=6) + geom_point(size=4, aes(colour=lab)) + scale_colour_manual(values=c("white", "grey", "black"))

t1 <- data.frame(
	model = c("1. y ~ causal", "2. y ~ g1", "3. y ~ g2", "4. y ~ g1 + g2"),
	rsq = round(c(mod1$r.squared, mod2$r.squared, mod3$r.squared, mod4$r.squared), 3)
)

t1$fac <- "R-square"

p3 <- ggplot(t1, aes(y=model, x=1)) +
geom_text(aes(label=rsq)) +
facet_grid(. ~ fac) +
theme_bw() +
scale_y_discrete(limits=rev(levels(t1$model))) +
theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), rect=element_blank(), axis.line=element_blank(), panel.border=element_blank())

p <- arrangeGrob(p1, p3, p2, nrow=2)
print(p)

grid.arrange(p1, p2, p3, nrow=2)
return(list(q))





y <- x[,5] * sqrt(0.1) + rnorm(n, sd=sqrt(0.9))
cor(y, x[,5])^2

summary(lm(y ~ x[,5]))
summary(lm(y ~ x[,4]))
summary(lm(y ~ x[,15]))
summary(lm(y ~ x[,4] + x[,15]))

pval <- rep(0, ncol(x))
for(i in 1:ncol(x))
{
	pval[i] <- summary(lm(y ~ x[,i]))$coefficients[2,4]
}

plot(-log10(pval), ylim=c(0,max(-log10(pval))))

