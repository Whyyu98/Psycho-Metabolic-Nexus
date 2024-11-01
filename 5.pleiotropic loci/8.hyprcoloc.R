setwd("E:/Dsektop/ND-PD/HyPrColoc/ALS")
library(data.table)
library(TwoSampleMR)
library(hyprcoloc)

SNP=fread("SNP.csv")
data1=fread("ALS_noP.txt")
data1_b=data1[, c("SNP", "b")]
data2=fread("ADHD_noP.txt")
data2_b=data2[, c("SNP", "b")]
data3<-extract_outcome_data(snps = SNP$SNP,outcomes ="ebi-a-GCST90001391")
data3_b=data3[, c("SNP", "beta.outcome")]
df1 <- merge(SNP, data1_b, by="SNP")
df2<- merge(df1, data2_b, by="SNP")
df_b<- merge(df2, data3_b, by="SNP")
betas=data.frame(SNPs=df_b$SNP,T1=df_b$b.x,T2=df_b$b.y,T3=df_b$beta.outcome)
rownames(betas)<-betas$SNPs
betas=subset(betas,select=-SNPs)


data1_se=data1[, c("SNP", "se")]
data2_se=data2[, c("SNP", "se")]
data3_se=data3[, c("SNP", "se.outcome")]
df3<- merge(SNP, data1_se, by="SNP")
df4<- merge(df3, data2_se, by="SNP")
df_se<- merge(df4, data3_se, by="SNP")
ses=data.frame(SNPs=df_se$SNP,T1=df_se$se.x,T2=df_se$se.y,T3=df_se$se.outcome)
rownames(ses)<-ses$SNPs
ses=subset(ses,select=-SNPs)


traits <- paste0("T", 1:dim(betas)[2])
rsid <- rownames(betas)
betas=data.matrix(betas)
ses=data.matrix(ses)
res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
a <- data.frame(do.call("cbind", res))
