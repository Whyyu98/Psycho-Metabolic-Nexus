data1=fread("finngen.txt")
data1<- data1[, c("rsids", "beta", "se")]
colnames(data1) <- c("SNP", "beta1", "se1")

data2=fread("PGC.txt")
data2<- data2[, c("MarkerName", "LogOR", "StdErrLogOR")]
colnames(data2) <- c("SNP", "beta2", "se2")
f=merge(data1,data2,by="SNP")

T1 <- c(f$beta1)
T2 <- c(f$beta2)
V1 <- c(f$se1)^2
V2 <- c(f$se2)^2

T =(T1-T2)^2/(V1+V2)

lambda_meta=median(T)/median(0.455)
print(lambda_meta)
