
#library(xlsx)
library(bigreadr)
library(plyr)
library(dplyr)
library(stringr)
# read ukb file 
PHENO_FILE <- "/public/home/Datasets/ukb/ukb47503.csv.gz"
ICD10_main <- paste0("41202-0.", 0:74)
df_eid <- fread2(PHENO_FILE, select = c("eid",'31-0.0','21022-0.0','738-0.0'))
df_ICD10 <- fread2(PHENO_FILE, select = ICD10_main)

df_index <- fread2(PHENO_FILE, select = c('21001-0.0','48-0.0','30740-0.0',
                                        '30750-0.0','30760-0.0','30780-0.0',
                                        '30870-0.0','30690-0.0','4080-0.0',
                                        '4079-0.0'))
colnames(df_index) <- c('BMI','Waist to hip ratio','Fasting glucose',
                        'Glycated hemoglobin','High-density lipoprotein cholesterol',
                        'Low-density lipoprotein cholesterol','Triglycerides',
                        'Total cholesterol','Systolic blood pressure','Diastolic blood pressure')
colnames(df_eid)[2:4] <- c('sex','age','es2')
# psychiatry 
ind_mdd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F32"))
))))
mdd <- rep(0, nrow(df_eid))
mdd[ind_mdd] <- 1

ind_bp <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F31"))
))))
bp <- rep(0, nrow(df_eid))
bp[ind_bp] <- 1

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

ind_ptsd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F43"))
))))
ptsd <- rep(0, nrow(df_eid))
ptsd[ind_ptsd] <- 1

ind_adhd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F90"))
))))
adhd <- rep(0, nrow(df_eid))
adhd[ind_adhd] <- 1

ind_asd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F84"))
))))
asd <- rep(0, nrow(df_eid))
asd[ind_asd] <- 1

ind_ocd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F42"))
))))
ocd <- rep(0, nrow(df_eid))
ocd[ind_ocd] <- 1

psy_dis <- data.frame(mdd,bp,scz,ad,ptsd,asd,ocd)
colnames(psy_dis) <- c('Major depressive disorder','Bipolar disorder','Schizophrenia',
                       'Anxiety disorders','Post traumatic stress disorder',
                       'Autism spectrum disorder','Obsessive-compulsive disorder')
# metabolic diseases
ind_t2d <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
t2d <- rep(0, nrow(df_eid))
t2d[ind_t2d] <- 1

ind_hyt <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I10"))
))))
hyt <- rep(0, nrow(df_eid))
hyt[ind_hyt] <- 1

ind_hyc <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E78"))
))))
hyc <- rep(0, nrow(df_eid))
hyc[ind_hyc] <- 1


ind_cad <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I25"))
))))
cad <- rep(0, nrow(df_eid))
cad[ind_cad] <- 1

ind_hf <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I50"))
))))
hf <- rep(0, nrow(df_eid))
hf[ind_hf] <- 1

ind_stroke <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I63"))
))))
stroke <- rep(0, nrow(df_eid))
stroke[ind_stroke] <- 1

ind_mi <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I21"))
))))
mi <- rep(0, nrow(df_eid))
mi[ind_mi] <- 1

ind_af <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I48"))
))))
af <- rep(0, nrow(df_eid))
af[ind_af] <- 1

ind_tia <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "G45"))
))))
tia <- rep(0, nrow(df_eid))
tia[ind_tia] <- 1

ind_vt <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I82"))
))))
vt <- rep(0, nrow(df_eid))
vt[ind_vt] <- 1

ind_ih <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I61"))
))))
ih <- rep(0, nrow(df_eid))
ih[ind_ih] <- 1

ind_ckd <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "N18"))
))))
ckd <- rep(0, nrow(df_eid))
ckd[ind_ckd] <- 1

ind_cardi <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I42"))
))))
cardi <- rep(0, nrow(df_eid))
cardi[ind_cardi] <- 1

ind_as <- sort(unique(unlist(c(
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "I71"))
))))
as <- rep(0, nrow(df_eid))
as[ind_as] <- 1

metabolic_dis <- data.frame(t2d,hyt,hyc,cad,hf,stroke,
                            mi,af,tia,vt,ih,ckd,cardi,as)

colnames(metabolic_dis) <- c('Type 2 diabetes','Hypertension','Hypercholesterolaemia',
                            'Coronary artery disease','Heart failure','Stroke',
                            'Myocardial infarction','Atrial fibrillation','Transient ischemic attack',
                            'Venous thromboembolism','Intracerebral hemorrhage',
                            'Chronic kidney disease','Cardiomyopathy','Aortic aneurysm')

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

case_control <- plyr::alply(15:35,1,function(xvar){
  case <- length(which(total_var[,xvar]==1))
  return(data.frame(disease=colnames(total_var)[xvar],case))
}) %>% do.call('rbind',.)
# chisq test
chisq <- plyr::alply(15:21,1,function(mm){
    metab <- plyr::alply(22:35,1,function(kk){
        four_table <- table(total_var[,mm],total_var[,kk])
        if (four_table[1,2]+four_table[2,1] >= 40) {
            results <- mcnemar.test(four_table,correct=F)
        } else {
            results <- mcnemar.test(four_table,correct=T)
        }
        p <- results$p.value
        chi <- results$statistic
        b <- four_table[1,2]
        c <- four_table[2,1]
        var1 <- colnames(total_var)[mm]
        var2 <- colnames(total_var)[kk]
        last_res <- data.frame(var1,var2,b,c,chi,p)
        print(kk)
        return(last_res)
    }) %>% do.call('rbind',.)
    return(metab)
}) %>% do.call('rbind',.)

write.csv(chisq,file = '/public/home/biostat07/project/mr/chisq_test.csv',
        row.names=F,quote=F)

# logit regression
total_var_cov <- cbind(total_var,cov[,c(paste0('PC',1:10))])
cov_use <- paste0(c('age','sex','es2',paste0('PC',1:10)),collapse='+')
## univariate logistic regression
logit <- plyr::alply(15:21,1,function(mm){
    logit1 <- plyr::alply(c(5:14,22:35),1,function(kk){
        formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+',cov_use)
        model0 <- glm(formula,family = 'binomial',data=total_var_cov)
        beta <- coef(summary(model0))[2,1]
        se <- coef(summary(model0))[2,2]
        Z <- coef(summary(model0))[2,3]
        P <- coef(summary(model0))[2,4]
        outcome <- colnames(total_var_cov)[mm]
        exposure <- rownames(coef(summary(model0)))[2]
        OR <- exp(beta)
        lowCI <- exp(beta-1.96*se)
        highCI <- exp(beta+1.96*se)
        results <- data.frame(outcome,exposure,beta,se,OR,lowCI,highCI,Z,P)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)

write.csv(logit,file = '/public/home/biostat07/project/mr/univar_logistic.csv',
        row.names=F,quote=F)

## interaction analysis
logit <- plyr::alply(15:21,1,function(mm){
    logit1 <- plyr::alply(c(5:14),1,function(kk){
        logit2 <- plyr::alply(c(22:35),1,function(vv){
            formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+','`',
            colnames(total_var_cov)[vv],'`+','`',colnames(total_var_cov)[kk],'`*`',colnames(total_var_cov)[vv],'`+',cov_use)
            model0 <- glm(formula,family = 'binomial',data=total_var_cov)
            beta_exp1 <- coef(summary(model0))[2,1]
            se_exp1 <- coef(summary(model0))[2,2]
            Z_exp1 <- coef(summary(model0))[2,3]
            P_exp1 <- coef(summary(model0))[2,4]
            OR_exp1 <- exp(beta_exp1)
            lowCI_exp1 <- exp(beta_exp1-1.96*se_exp1)
            highCI_exp1 <- exp(beta_exp1+1.96*se_exp1)
            beta_exp2 <- coef(summary(model0))[3,1]
            se_exp2 <- coef(summary(model0))[3,2]
            Z_exp2 <- coef(summary(model0))[3,3]
            P_exp2 <- coef(summary(model0))[3,4]
            OR_exp2 <- exp(beta_exp2)
            lowCI_exp2 <- exp(beta_exp2-1.96*se_exp2)
            highCI_exp2 <- exp(beta_exp2+1.96*se_exp2)
            beta_inter <- coef(summary(model0))[17,1]
            se_inter <- coef(summary(model0))[17,2]
            Z_inter <- coef(summary(model0))[17,3]
            P_inter <- coef(summary(model0))[17,4]
            OR_inter <- exp(beta_inter)
            lowCI_inter <- exp(beta_inter-1.96*se_inter)
            highCI_inter <- exp(beta_inter+1.96*se_inter)
            outcome <- colnames(total_var_cov)[mm]
            exposure1 <- rownames(coef(summary(model0)))[2]
            exposure2 <- rownames(coef(summary(model0)))[3]
            inter <- rownames(coef(summary(model0)))[17]
            results <- data.frame(outcome,exposure1,beta_exp1,se_exp1,OR_exp1,lowCI_exp1,highCI_exp1,Z_exp1,P_exp1,
            exposure2,beta_exp2,se_exp2,OR_exp2,lowCI_exp2,highCI_exp2,Z_exp2,P_exp2,
            inter,beta_inter,se_inter,OR_inter,lowCI_inter,highCI_inter,Z_inter,P_inter)
            return(results)
        }) %>% do.call('rbind',.)
        return(logit2)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/interaction_logistic.csv',
        row.names=F,quote=F)

# metabolic diseases interaction analysis
meta_b <- c(5:14)

logit <- plyr::alply(15:21,1,function(mm){
    logit1 <- plyr::alply(c(5:13),1,function(kk){
        index <- kk+1
        logit2 <- plyr::alply(c(index:14),1,function(vv){
            formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+','`',
            colnames(total_var_cov)[vv],'`+','`',colnames(total_var_cov)[kk],'`*`',colnames(total_var_cov)[vv],'`+',cov_use)
            model0 <- glm(formula,family = 'binomial',data=total_var_cov)
            beta_exp1 <- coef(summary(model0))[2,1]
            se_exp1 <- coef(summary(model0))[2,2]
            Z_exp1 <- coef(summary(model0))[2,3]
            P_exp1 <- coef(summary(model0))[2,4]
            OR_exp1 <- exp(beta_exp1)
            lowCI_exp1 <- exp(beta_exp1-1.96*se_exp1)
            highCI_exp1 <- exp(beta_exp1+1.96*se_exp1)
            beta_exp2 <- coef(summary(model0))[3,1]
            se_exp2 <- coef(summary(model0))[3,2]
            Z_exp2 <- coef(summary(model0))[3,3]
            P_exp2 <- coef(summary(model0))[3,4]
            OR_exp2 <- exp(beta_exp2)
            lowCI_exp2 <- exp(beta_exp2-1.96*se_exp2)
            highCI_exp2 <- exp(beta_exp2+1.96*se_exp2)
            beta_inter <- coef(summary(model0))[17,1]
            se_inter <- coef(summary(model0))[17,2]
            Z_inter <- coef(summary(model0))[17,3]
            P_inter <- coef(summary(model0))[17,4]
            OR_inter <- exp(beta_inter)
            lowCI_inter <- exp(beta_inter-1.96*se_inter)
            highCI_inter <- exp(beta_inter+1.96*se_inter)
            outcome <- colnames(total_var_cov)[mm]
            exposure1 <- colnames(total_var_cov)[kk]
            exposure2 <- colnames(total_var_cov)[vv]
            inter <- rownames(coef(summary(model0)))[17]
            results <- data.frame(outcome,exposure1,beta_exp1,se_exp1,OR_exp1,lowCI_exp1,highCI_exp1,Z_exp1,P_exp1,
            exposure2,beta_exp2,se_exp2,OR_exp2,lowCI_exp2,highCI_exp2,Z_exp2,P_exp2,
            inter,beta_inter,se_inter,OR_inter,lowCI_inter,highCI_inter,Z_inter,P_inter)
            return(results)
        }) %>% do.call('rbind',.)
        return(logit2)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/interaction_logistic_metabolic.csv',
        row.names=F,quote=F)
# metabolic index interaction analysis
meta_index <- c(22:35)
logit <- plyr::alply(15:21,1,function(mm){
    logit1 <- plyr::alply(c(22:34),1,function(kk){
        cat(paste0(kk,':'))
        index <- kk+1
        logit2 <- plyr::alply(c(index:35),1,function(vv){
            cat(paste0(vv,':'))
            formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+','`',
            colnames(total_var_cov)[vv],'`+','`',colnames(total_var_cov)[kk],'`*`',colnames(total_var_cov)[vv],'`+',cov_use)
            model0 <- glm(formula,family = 'binomial',data=total_var_cov)
            beta_exp1 <- coef(summary(model0))[2,1]
            se_exp1 <- coef(summary(model0))[2,2]
            Z_exp1 <- coef(summary(model0))[2,3]
            P_exp1 <- coef(summary(model0))[2,4]
            OR_exp1 <- exp(beta_exp1)
            lowCI_exp1 <- exp(beta_exp1-1.96*se_exp1)
            highCI_exp1 <- exp(beta_exp1+1.96*se_exp1)
            beta_exp2 <- coef(summary(model0))[3,1]
            se_exp2 <- coef(summary(model0))[3,2]
            Z_exp2 <- coef(summary(model0))[3,3]
            P_exp2 <- coef(summary(model0))[3,4]
            OR_exp2 <- exp(beta_exp2)
            lowCI_exp2 <- exp(beta_exp2-1.96*se_exp2)
            highCI_exp2 <- exp(beta_exp2+1.96*se_exp2)
            if (nrow(coef(summary(model0)))==17) {
              beta_inter <- coef(summary(model0))[17,1]
              se_inter <- coef(summary(model0))[17,2]
              Z_inter <- coef(summary(model0))[17,3]
              P_inter <- coef(summary(model0))[17,4]
              inter <- rownames(coef(summary(model0)))[17]
            } else {
              beta_inter <- NA
              se_inter <- NA
              Z_inter <- NA
              P_inter <- NA
              inter <- NA
            }
            OR_inter <- exp(beta_inter)
            lowCI_inter <- exp(beta_inter-1.96*se_inter)
            highCI_inter <- exp(beta_inter+1.96*se_inter)
            outcome <- colnames(total_var_cov)[mm]
            exposure1 <- colnames(total_var_cov)[kk]
            exposure2 <- colnames(total_var_cov)[vv]
            results <- data.frame(outcome,exposure1,beta_exp1,se_exp1,OR_exp1,lowCI_exp1,highCI_exp1,Z_exp1,P_exp1,
            exposure2,beta_exp2,se_exp2,OR_exp2,lowCI_exp2,highCI_exp2,Z_exp2,P_exp2,
            inter,beta_inter,se_inter,OR_inter,lowCI_inter,highCI_inter,Z_inter,P_inter)
            return(results)
        }) %>% do.call('rbind',.)
        return(logit2)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/interaction_logistic_metabolic_index.csv',
        row.names=F,quote=F)


# metabolic diseases
## univariate logistic regression
logit <- plyr::alply(22:35,1,function(mm){
    logit1 <- plyr::alply(c(5:14,15:21),1,function(kk){
        formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+',cov_use)
        model0 <- glm(formula,family = 'binomial',data=total_var_cov)
        beta <- coef(summary(model0))[2,1]
        se <- coef(summary(model0))[2,2]
        Z <- coef(summary(model0))[2,3]
        P <- coef(summary(model0))[2,4]
        outcome <- colnames(total_var_cov)[mm]
        exposure <- rownames(coef(summary(model0)))[2]
        OR <- exp(beta)
        lowCI <- exp(beta-1.96*se)
        highCI <- exp(beta+1.96*se)
        results <- data.frame(outcome,exposure,beta,se,OR,lowCI,highCI,Z,P)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/metabolic_as_outcome.csv',
        row.names=F,quote=F)

logit <- plyr::alply(5:14,1,function(mm){
    logit1 <- plyr::alply(c(22:35,15:21),1,function(kk){
        formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+',cov_use)
        model0 <- lm(formula,data=total_var_cov)
        beta <- coef(summary(model0))[2,1]
        se <- coef(summary(model0))[2,2]
        Z <- coef(summary(model0))[2,3]
        P <- coef(summary(model0))[2,4]
        outcome <- colnames(total_var_cov)[mm]
        exposure <- rownames(coef(summary(model0)))[2]
        #OR <- exp(beta)
        lowCI <- beta-1.96*se
        highCI <- beta+1.96*se
        results <- data.frame(outcome,exposure,beta,se,lowCI,highCI,Z,P)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/metabolic_index_as_outcome.csv',
        row.names=F,quote=F)

## interaction analysis
logit <- plyr::alply(22:35,1,function(mm){
    logit1 <- plyr::alply(c(5:14),1,function(kk){
        logit2 <- plyr::alply(c(15:21),1,function(vv){
            formula <- paste0('`',colnames(total_var_cov)[mm],'`~','`',colnames(total_var_cov)[kk],'`+','`',
            colnames(total_var_cov)[vv],'`+','`',colnames(total_var_cov)[kk],'`*`',colnames(total_var_cov)[vv],'`+',cov_use)
            model0 <- glm(formula,family = 'binomial',data=total_var_cov)
            beta_exp1 <- coef(summary(model0))[2,1]
            se_exp1 <- coef(summary(model0))[2,2]
            Z_exp1 <- coef(summary(model0))[2,3]
            P_exp1 <- coef(summary(model0))[2,4]
            OR_exp1 <- exp(beta_exp1)
            lowCI_exp1 <- exp(beta_exp1-1.96*se_exp1)
            highCI_exp1 <- exp(beta_exp1+1.96*se_exp1)
            beta_exp2 <- coef(summary(model0))[3,1]
            se_exp2 <- coef(summary(model0))[3,2]
            Z_exp2 <- coef(summary(model0))[3,3]
            P_exp2 <- coef(summary(model0))[3,4]
            OR_exp2 <- exp(beta_exp2)
            lowCI_exp2 <- exp(beta_exp2-1.96*se_exp2)
            highCI_exp2 <- exp(beta_exp2+1.96*se_exp2)
            beta_inter <- coef(summary(model0))[17,1]
            se_inter <- coef(summary(model0))[17,2]
            Z_inter <- coef(summary(model0))[17,3]
            P_inter <- coef(summary(model0))[17,4]
            OR_inter <- exp(beta_inter)
            lowCI_inter <- exp(beta_inter-1.96*se_inter)
            highCI_inter <- exp(beta_inter+1.96*se_inter)
            outcome <- colnames(total_var_cov)[mm]
            exposure1 <- rownames(coef(summary(model0)))[2]
            exposure2 <- rownames(coef(summary(model0)))[3]
            inter <- rownames(coef(summary(model0)))[17]
            results <- data.frame(outcome,exposure1,beta_exp1,se_exp1,OR_exp1,lowCI_exp1,highCI_exp1,Z_exp1,P_exp1,
            exposure2,beta_exp2,se_exp2,OR_exp2,lowCI_exp2,highCI_exp2,Z_exp2,P_exp2,
            inter,beta_inter,se_inter,OR_inter,lowCI_inter,highCI_inter,Z_inter,P_inter)
            return(results)
        }) %>% do.call('rbind',.)
        return(logit2)
    }) %>% do.call('rbind',.)
    return(logit1)
}) %>% do.call('rbind',.)
write.csv(logit,file = '/public/home/biostat07/project/mr/interaction_logistic.csv',
        row.names=F,quote=F)
