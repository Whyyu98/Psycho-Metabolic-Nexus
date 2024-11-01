#! /usr/bin/env Rscript
#library(xlsx)
library(bigreadr)
library(plyr)
library(dplyr)
library(stringr)
# read ukb file 
PHENO_FILE <- "/public/home/Datasets/ukb/ukb47503.csv.gz"
ICD10_main <- paste0("41202-0.", 0:74)
ocd_main <- c(paste0('20002-0.',0:33),paste0('20002-1.',0:33),paste0('20002-2.',0:33),paste0('20002-3.',0:33))
df_eid <- fread2(PHENO_FILE, select = c("eid",'31-0.0','21022-0.0','738-0.0'))
df_ICD10 <- fread2(PHENO_FILE, select = ICD10_main)
df_ocd <- fread2(PHENO_FILE, select = ocd_main)

df_index <- fread2(PHENO_FILE, select = c('21001-0.0'))


colnames(df_index) <- c('BMI')
colnames(df_eid)[2:4] <- c('sex','age','es2')
# psychiatry 
ind_mdd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F32"))
))))
mdd <- rep(0, nrow(df_eid))
mdd[ind_mdd] <- 1

ind_bd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F31"))
))))
bd <- rep(0, nrow(df_eid))
bd[ind_bd] <- 1

ind_scz <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F20"))
))))
scz <- rep(0, nrow(df_eid))
scz[ind_scz] <- 1

ind_ad <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F41"))
))))
ad <- rep(0, nrow(df_eid))
ad[ind_ad] <- 1


psy_dis <- data.frame(mdd,bd,scz,ad)
colnames(psy_dis) <- c('MDD','BD','SCZ','ANX')

# metabolic diseases
ind_t2d <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
t2d <- rep(0, nrow(df_eid))
t2d[ind_t2d] <- 1

ind_cad <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I25"))
))))
cad <- rep(0, nrow(df_eid))
cad[ind_cad] <- 1

metabolic_dis <- data.frame(t2d,cad)

colnames(metabolic_dis) <- c('T2D','CAD')

# total variates
eur_pop <- fread2('/public/home/Datasets/ukb/pheno/eid_EUR.txt',header=F)
total_var <- cbind(df_eid,df_index,psy_dis,metabolic_dis)

# cov 
cov <- fread2('/public/home/Datasets/ukb/pheno/sqc.txt')

# match and keep eur population
total_var <- total_var %>% filter(eid %in% eur_pop[,1])
cov <- cov %>% filter(eid %in% eur_pop[,1])

total_var <- total_var[match(eur_pop[,1],total_var$eid),]
cov <- cov[match(eur_pop[,1],cov$eid),]


# logit regression
total_var_cov <- cbind(total_var,cov[,c(paste0('PC',1:10))])
cov_use <- paste0(c('age','sex','es2',paste0('PC',1:10)),collapse='+')

library(pROC)
library(pscl)
# get three-models auc
pairs <- fread2('/public/home/biostat07/project/mr/new_PRS/pairs.txt',header=F)
thr_models1 <- plyr::alply(1:nrow(pairs),1,function(xvar){
    var1 <- pairs[xvar,1]
    var2 <- pairs[xvar,2]
    pgs1 <- fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var1,'/pred.profile'))[,6]
    pgs2 <- fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var1,'-',var2,'/pred.profile'))[,6]
    pgs3 <-  fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var2,'/pred.profile'))[,6]
    list_var <- cbind(total_var_cov,pgs1,pgs2,pgs3)
    if ( var1 == 'BMI' ) {
        # null model
        null_model <- lm(paste0(var1,'~',cov_use),data=list_var)
        beta_null <- coef(summary(null_model))[,1]
        pred_null <- beta_null[1] + as.matrix(list_var[,c('age','sex','es2',paste0('PC',1:10))]) %*% beta_null[-1]
        frame_null <- na.omit(data.frame(list_var$BMI,as.vector(pred_null)))
        ci_null <- plyr::alply(1:1000,1,function(xvar){
            set.seed(2023110801+xvar)
            data <- frame_null[sample(nrow(frame_null),10000,replace=F),]
            return(cor(data[,1],data[,2])**2)
        }) %>% do.call('rbind',.) %>% as.vector()
        R2_null <- mean(ci_null)
        # uni-pgs model
        uni_model <- lm(paste0(var1,'~','pgs1+',cov_use),data=list_var)
        beta_uni <- coef(summary(uni_model))[,1]
        pred_uni <- beta_uni[1] + as.matrix(list_var[,c('pgs1','age','sex','es2',paste0('PC',1:10))]) %*% beta_uni[-1]
        frame_uni <- na.omit(data.frame(list_var$BMI,as.vector(pred_uni)))
        ci_uni <- plyr::alply(1:1000,1,function(xvar){
            set.seed(2023110801+xvar)
            data <- frame_uni[sample(nrow(frame_uni),10000,replace=F),]
            return(cor(data[,1],data[,2])**2)
        }) %>% do.call('rbind',.) %>% as.vector()
        R2_uni <- cor(frame_uni[,1],frame_uni[,2])**2
        rise_uni <- (R2_uni-R2_null)/R2_null
        # uni-cov model
        cov_model <- lm(paste0(var1,'~','pgs3+',cov_use),data=list_var)
        beta_cov <- coef(summary(cov_model))[,1]
        pred_cov <- beta_cov[1] + as.matrix(list_var[,c('pgs3','age','sex','es2',paste0('PC',1:10))]) %*% beta_cov[-1]
        frame_cov <- na.omit(data.frame(list_var$BMI,as.vector(pred_cov)))
        ci_cov <- plyr::alply(1:1000,1,function(xvar){
            set.seed(2023110801+xvar)
            data <- frame_cov[sample(nrow(frame_cov),10000,replace=F),]
            return(cor(data[,1],data[,2])**2)
        }) %>% do.call('rbind',.) %>% as.vector()
        R2_cov <- cor(frame_cov[,1],frame_cov[,2])**2
        rise_cov <- (R2_cov-R2_null)/R2_null
        # mtpgs model
        mt_model <- lm(paste0(var1,'~','pgs2+',cov_use),data=list_var)
        beta_mt <- coef(summary(mt_model))[,1]
        pred_mt <- beta_mt[1] + as.matrix(list_var[,c('pgs2','age','sex','es2',paste0('PC',1:10))]) %*% beta_mt[-1]
        frame_mt <- na.omit(data.frame(list_var$BMI,as.vector(pred_mt)))
        ci_mt <- plyr::alply(1:1000,1,function(xvar){
            set.seed(2023110801+xvar)
            data <- frame_mt[sample(nrow(frame_mt),10000,replace=F),]
            return(cor(data[,1],data[,2])**2)
        }) %>% do.call('rbind',.) %>% as.vector()
        R2_mt <- cor(frame_mt[,1],frame_mt[,2])**2
        rise_mt <- (R2_mt-R2_null)/R2_null
        # two dbslmm model
        bi_model <- lm(paste0(var1,'~','pgs1+pgs3+',cov_use),data=list_var)
        beta_bi <- coef(summary(bi_model))[,1]
        pred_bi <- beta_bi[1] + as.matrix(list_var[,c('pgs1','pgs3','age','sex','es2',paste0('PC',1:10))]) %*% beta_bi[-1]
        frame_bi <- na.omit(data.frame(list_var$BMI,as.vector(pred_bi)))
        ci_bi <- plyr::alply(1:1000,1,function(xvar){
            set.seed(2023110801+xvar)
            data <- frame_bi[sample(nrow(frame_bi),10000,replace=F),]
            return(cor(data[,1],data[,2])**2)
        }) %>% do.call('rbind',.) %>% as.vector()
        R2_bi <- cor(frame_bi[,1],frame_bi[,2])**2
        rise_bi <- (R2_bi-R2_null)/R2_null
        # results
        outcome <- rep(var1,4)
        exposure <- rep(var2,4)
        type <- c('DBSLMM','var2-pgs-only','MTPGS','TWO-MODELS')
        AUC_OR_R2 <- c(mean(ci_uni),mean(ci_cov),mean(ci_mt),mean(ci_bi))
        lowCI <- c(quantile(ci_uni,0.025),quantile(ci_cov,0.025),quantile(ci_mt,0.025),quantile(ci_bi,0.025))
        highCI <- c(quantile(ci_uni,0.975),quantile(ci_cov,0.975),quantile(ci_mt,0.975),quantile(ci_bi,0.975))
        null_model_R2 <- rep(R2_null,4)
        full_model_R2 <- c(R2_uni,R2_cov,R2_mt,R2_bi)
        R2_rise_ratio <- c(rise_uni,rise_cov,rise_mt,rise_bi)
        results <- data.frame(outcome,exposure,type,AUC_OR_R2,lowCI,highCI,null_model_R2,full_model_R2,R2_rise_ratio)
        return(results)
    } else {
        # null model
        model0 <- glm(paste0(var1,'~',1),family='binomial',data=list_var)
        null_model <- glm(paste0(var1,'~',cov_use),family='binomial',data=list_var)
        R2_null <- pR2(null_model,model0)['McFadden']
        beta_null <- coef(summary(null_model))[,1]
        pred_null <- beta_null[1] + as.matrix(list_var[,c('age','sex','es2',paste0('PC',1:10))]) %*% beta_null[-1]
        frame_null <- na.omit(data.frame(list_var$BMI,as.vector(pred_null)))
        #R2_null <- cor(frame_null[,1],frame_null[,2])**2
        # uni-pgs model
        uni_model <- glm(paste0(var1,'~','pgs1+',cov_use),family='binomial',data=list_var)
        R2_uni <- pR2(uni_model,model0)['McFadden']
        beta_uni <- coef(summary(uni_model))[,1]
        pred_uni <- beta_uni[1] + as.matrix(list_var[,c('pgs1','age','sex','es2',paste0('PC',1:10))]) %*% beta_uni[-1]
        roc_uni <- roc(list_var[,var1],as.vector(pred_uni))
        auc_uni <- ci.auc(roc_uni)
        rise_uni <- (R2_uni-R2_null)/R2_null
        # uni-cov model
        cov_model <- glm(paste0(var1,'~','pgs3+',cov_use),family='binomial',data=list_var)
        R2_cov <- pR2(cov_model,model0)['McFadden']
        beta_cov <- coef(summary(cov_model))[,1]
        pred_cov <- beta_cov[1] + as.matrix(list_var[,c('pgs3','age','sex','es2',paste0('PC',1:10))]) %*% beta_cov[-1]
        roc_cov <- roc(list_var[,var1],as.vector(pred_cov))
        auc_cov <- ci.auc(roc_cov)
        rise_cov <- (R2_cov-R2_null)/R2_null
        # mtpgs model
        mt_model <- glm(paste0(var1,'~','pgs2+',cov_use),family='binomial',data=list_var)
        R2_mt <- pR2(mt_model,model0)['McFadden']
        beta_mt <- coef(summary(mt_model))[,1]
        pred_mt <- beta_mt[1] + as.matrix(list_var[,c('pgs2','age','sex','es2',paste0('PC',1:10))]) %*% beta_mt[-1]
        roc_mt <- roc(list_var[,var1],as.vector(pred_mt))
        auc_mt <- ci.auc(roc_mt)
        rise_mt <- (R2_mt-R2_null)/R2_null
        # two dbslmm model
        bi_model <- glm(paste0(var1,'~','pgs1+pgs3+',cov_use),family='binomial',data=list_var)
        R2_bi <- pR2(bi_model,model0)['McFadden']
        beta_bi <- coef(summary(bi_model))[,1]
        pred_bi <- beta_bi[1] + as.matrix(list_var[,c('pgs1','pgs3','age','sex','es2',paste0('PC',1:10))]) %*% beta_bi[-1]
        roc_bi <- roc(list_var[,var1],as.vector(pred_bi))
        auc_bi <- ci.auc(roc_bi)
        rise_bi <- (R2_bi-R2_null)/R2_null
        # results
        outcome <- rep(var1,4)
        exposure <- rep(var2,4)
        type <- c('DBSLMM','var2-pgs-only','MTPGS','TWO-MODELS')
        AUC_OR_R2 <- c(auc_uni[2],auc_cov[2],auc_mt[2],auc_bi[2])
        lowCI <- c(auc_uni[1],auc_cov[1],auc_mt[1],auc_bi[1])
        highCI <- c(auc_uni[3],auc_cov[1],auc_mt[3],auc_bi[3])
        null_model_R2 <- rep(R2_null,4)
        full_model_R2 <- c(R2_uni,R2_cov,R2_mt,R2_bi)
        R2_rise_ratio <- c(rise_uni,rise_cov,rise_mt,rise_bi)
        results <- data.frame(outcome,exposure,type,AUC_OR_R2,lowCI,highCI,null_model_R2,full_model_R2,R2_rise_ratio)
        return(results)
    } 
}) %>% do.call('rbind',.)
write.csv(thr_models1,file='/public/home/biostat07/project/mr/new_PRS/r2_rise_ratio.csv',row.names=F,quote=F)

thr_models2 <- plyr::alply(1:nrow(pairs),1,function(xvar){
    var1 <- pairs[xvar,1]
    var2 <- pairs[xvar,2]
    pgs1 <- fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var1,'/pred.profile'))[,6]
    pgs2 <- fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var1,'-',var2,'/pred.profile'))[,6]
    pgs3 <-  fread2(paste0('/public/home/biostat07/project/mr/new_PRS/',var2,'/pred.profile'))[,6]
    list_var <- cbind(total_var_cov,pgs1,pgs2,pgs3)
    if ( var1 == 'BMI' ) {
        # uni-pgs model
        uni_model <- lm(paste0(var1,'~','pgs1+',cov_use),data=list_var)
        beta_uni <- coef(summary(uni_model))[2,1]
        se_uni <- coef(summary(uni_model))[2,2]
        z_uni <- coef(summary(uni_model))[2,3]
        p_uni <- coef(summary(uni_model))[2,4]
        # uni-cov model
        cov_model <- lm(paste0(var1,'~','pgs3+',cov_use),data=list_var)
        beta_cov <- coef(summary(cov_model))[2,1]
        se_cov <- coef(summary(cov_model))[2,2]
        z_cov <- coef(summary(cov_model))[2,3]
        p_cov <- coef(summary(cov_model))[2,4]
        # two dbslmm model
        bi_model <- lm(paste0(var1,'~','pgs1+pgs3+',cov_use),data=list_var)
        anova_model <- anova(uni_model,bi_model)
        f_bi <- anova_model$F[2]
        p_bi <- anova_model$P[2]
        # results
        outcome <- rep(var1,3)
        exposure <- rep(var2,3)
        type <- c('DBSLMM','var2-pgs-only','TWO-MODELS')
        beta <- c(beta_uni,beta_cov,NA)
        se <- c(se_uni,se_cov,NA)
        statistics <- c('Z-SCORE','Z-SCORE','F')
        statistic_value <- c(z_uni,z_cov,f_bi)
        P <- c(p_uni,p_cov,p_bi)
        results <- data.frame(outcome,exposure,type,beta,se,statistics,statistic_value,P)
        return(results)
    } else {
        # uni-pgs model
        uni_model <- glm(paste0(var1,'~','pgs1+',cov_use),family='binomial',data=list_var)
        beta_uni <- coef(summary(uni_model))[2,1]
        se_uni <- coef(summary(uni_model))[2,2]
        z_uni <- coef(summary(uni_model))[2,3]
        p_uni <- coef(summary(uni_model))[2,4]
        # uni-cov model
        cov_model <- glm(paste0(var1,'~','pgs3+',cov_use),family='binomial',data=list_var)
        beta_cov <- coef(summary(cov_model))[2,1]
        se_cov <- coef(summary(cov_model))[2,2]
        z_cov <- coef(summary(cov_model))[2,3]
        p_cov <- coef(summary(cov_model))[2,4]
        # two dbslmm model
        bi_model <- glm(paste0(var1,'~','pgs1+pgs3+',cov_use),family='binomial',data=list_var)
        null_loglik <- logLik(uni_model)
        full_loglik <- logLik(bi_model)
        LR_bi <- 2*(as.numeric(full_loglik) - as.numeric(null_loglik))
        df <- 1
        p_bi <- pchisq(LR_bi, df, lower.tail = FALSE)
        # results
        outcome <- rep(var1,3)
        exposure <- rep(var2,3)
        type <- c('DBSLMM','var2-pgs-only','TWO-MODELS')
        beta <- c(beta_uni,beta_cov,NA)
        se <- c(se_uni,se_cov,NA)
        statistics <- c('Z-SCORE','Z-SCORE','Likelihood-ratio')
        statistic_value <- c(z_uni,z_cov,LR_bi)
        P <- c(p_uni,p_cov,p_bi)
        results <- data.frame(outcome,exposure,type,beta,se,statistics,statistic_value,P)
        return(results)
    } 
}) %>% do.call('rbind',.)
thr_models2$p.adjusted <- p.adjust(thr_models2$P,method='BH')
write.csv(thr_models2,file='/public/home/biostat07/project/mr/new_PRS/models_pgs_pvalue.csv',row.names=F,quote=F)

# pgs
pgs_tg <- fread2('/public/home/biostat07/project/mr/prs/tg/pred.profile')[,6]
pgs_ocd <- fread2('/public/home/biostat07/project/mr/prs/ocd/pred.profile')[,6]
pgs_ob <- fread2('/public/home/biostat07/project/mr/prs/ob/pred.profile')[,6]
pgs_ob_adhd <- fread2('/public/home/biostat07/project/mr/prs/ob-adhd/pred.profile')[,6]
pgs_mdd_ob <- fread2('/public/home/biostat07/project/mr/prs/mdd-ob/pred.profile')[,6]
pgs_bmi_ocd <- fread2('/public/home/biostat07/project/mr/prs/bmi-ocd/pred.profile')[,6]
pgs_bmi_adhd <- fread2('/public/home/biostat07/project/mr/prs/bmi-adhd/pred.profile')[,6]
pgs_adhd_tg <- fread2('/public/home/biostat07/project/mr/prs/adhd-tg/pred.profile')[,6]
pgs_adhd <- fread2('/public/home/biostat07/project/mr/prs/adhd/pred.profile')[,6]
pgs_bmi <- fread2('/public/home/biostat07/project/mr/prs/bmi/pred.profile')[,6]

total_pgs <- cbind(total_var_cov,pgs_tg,pgs_ocd,pgs_ob,
                  pgs_ob_adhd,pgs_mdd_ob,pgs_bmi_ocd,
                  pgs_bmi_adhd,pgs_adhd_tg,pgs_adhd,pgs_bmi)
library(pROC)
##### 1. bmi - ocd
uni_ocd <- glm(paste0('`Obsessive-compulsive disorder`~pgs_ocd+',cov_use),family='binomial',data=total_pgs)
beta_ocd <- coef(summary(uni_ocd))[,1]
pred_ocd <- beta_ocd[1] + as.matrix(total_pgs[,c('pgs_ocd','age','sex','es2',paste0('PC',1:10))]) %*% beta_ocd[-1]
auc_ocd <- roc(total_pgs$`Obsessive-compulsive disorder`,as.vector(pred_ocd))


uni_bmi_ocd <- glm(paste0('`Obsessive-compulsive disorder`~pgs_bmi_ocd+',cov_use),family='binomial',data=total_pgs)
beta_bmi_ocd <- coef(summary(uni_bmi_ocd))[,1]
pred_bmi_ocd <- beta_bmi_ocd[1] + as.matrix(total_pgs[,c('pgs_bmi_ocd','age','sex','es2',paste0('PC',1:10))]) %*% beta_bmi_ocd[-1]
auc_bmi_ocd <- roc(total_pgs$`Obsessive-compulsive disorder`,as.vector(pred_bmi_ocd))

bi_bmi_ocd <- glm(paste0('`Obsessive-compulsive disorder`~pgs_ocd+pgs_bmi+',cov_use),family='binomial',data=total_pgs)
beta_bi_ocd <- coef(summary(bi_bmi_ocd))[,1]
pred_bi_ocd <- beta_bi_ocd[1] + as.matrix(total_pgs[,c('pgs_ocd','pgs_bmi','age','sex','es2',paste0('PC',1:10))]) %*% beta_bi_ocd[-1]
auc_bi_ocd <- roc(total_pgs$`Obsessive-compulsive disorder`,as.vector(pred_bi_ocd))

##### 2. mdd - ob
uni_ob <- glm(paste0('`Obesity`~pgs_ob+',cov_use),family='binomial',data=total_pgs)
beta_ob <- coef(summary(uni_ob))[,1]
pred_ob <- beta_ob[1] + as.matrix(total_pgs[,c('pgs_ob','age','sex','es2',paste0('PC',1:10))]) %*% beta_ob[-1]
auc_ob <- roc(total_pgs$`Obesity`,as.vector(pred_ob))

uni_mdd_ob <- glm(paste0('`Obesity`~pgs_mdd_ob+',cov_use),family='binomial',data=total_pgs)
beta_mdd_ob <- coef(summary(uni_mdd_ob))[,1]
pred_mdd_ob <- beta_mdd_ob[1] + as.matrix(total_pgs[,c('pgs_mdd_ob','age','sex','es2',paste0('PC',1:10))]) %*% beta_mdd_ob[-1]
auc_mdd_ob <- roc(total_pgs$`Obesity`,as.vector(pred_mdd_ob))

bi_mdd_ob <- glm(paste0('`Obesity`~pgs_ob+',cov_use),family='binomial',data=total_pgs)
beta_mdd_ob <- coef(summary(uni_mdd_ob))[,1]
pred_mdd_ob <- beta_mdd_ob[1] + as.matrix(total_pgs[,c('pgs_mdd_ob','age','sex','es2',paste0('PC',1:10))]) %*% beta_mdd_ob[-1]
auc_mdd_ob <- roc(total_pgs$`Obesity`,as.vector(pred_mdd_ob))