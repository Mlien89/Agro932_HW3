library("BGLR")
ped <- read_ped("largedata/sativas413.ped") 

p=ped$p
n=ped$n
out=ped$x
#Recode snp to 0,1,2 format using allele 1
# 0 --> 0
# 1 --> 1
# 2 --> NA
# 3 --> 2
out[out==2]=NA
out[out==3]=2
Z <- matrix(out, nrow=p, ncol=n, byrow=TRUE)
Z <- t(Z) 
dim(Z) # # 413 x 36901

fam <- read.table("largedata/sativas413.fam", header = FALSE, stringsAsFactors = FALSE)  
head(fam)
rownames(Z) <- paste0("NSFTV_", fam$V2) # 413 x 36901

for (j in 1:ncol(Z)){
  Z[,j] <- ifelse(is.na(Z[,j]), mean(Z[,j], na.rm=TRUE), Z[,j])
}

map <- read.table("largedata/sativas413.map", header = FALSE, stringsAsFactors = FALSE)
mygeno <- data.frame(marker=map[,2], chrom=map[,1], pos=map[,4], t(Z-1), check.names = FALSE) # Z = \in{-1, 0, 1}
pheno <- read.delim("http://ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", header=TRUE)
pheno$NSFTVID <- paste0("NSFTV_", pheno$NSFTVID)
mypheno <- data.frame(NSFTV_ID=pheno$NSFTVID, y=pheno$Plant.height) 
res2 <- GWAS(mypheno, mygeno, fixed=NULL, n.PC=3, min.MAF=0.05, P3D=TRUE, plot=FALSE)

library(qqman)
pdf("graphs/mht_res2.pdf", width=10, height=5)
manhattan(x = res2, chr = "chrom", bp = "pos", p = "y", snp = "marker", col = c("blue4", "orange3"), logp = FALSE)
dev.off()