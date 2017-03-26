load("phen.RData")
# a <- read.table("scratch/chr16.fam.orig", stringsAsFactors=FALSE)
a <- read.table("scratch/chr16.fam", stringsAsFactors=FALSE)
b <- read.table("~/data/alspac_1kg/data/data.sample", stringsAsFactors=FALSE, skip=2)
b1 <- subset(b, V1 %in% a$V1)
phen$V1 <- a$V1
phen$V2 <- a$V2


g <- scan("~/repo/ImputationPractical/data/rsid.gen", what="numeric")
g <- as.numeric(g[-c(1:6)])
i1 <- seq(1, length(g), by=3)
i2 <- i1 + 1
i3 <- i1 + 2

G <- g[i1] * 2 + g[i2]


# index <- match(phen$V1, b1$V1)
index <- scan("chr16.index", what=numeric())
phen2 <- phen[order(index),]
stopifnot(all(phen2$V1 == b1$V1))
covars2 <- covars[order(index),]
stopifnot(all(covars2$V1 == phen2$fid))
covars2 <- covars2[,-c(1,2)]
names(covars2)[1:10] <- paste0("PC", 1:10)
phen2 <- subset(phen2, select=-c(V1, V2))
f <- data.frame(ID_1=phen2$fid, ID_2=phen2$iid, missing=0, mother=0, father=0, sex=covars2$sex)
phen2 <- phen2[,-c(1:2)]
covars2 <- subset(covars2, select=-c(sex))
dat <- cbind(f, phen2, covars2)

summary(lm(dat$bmi ~ G))


write.table(dat, file="~/repo/ImputationPractical/data.sample", row=FALSE, col=TRUE, qu=FALSE)

