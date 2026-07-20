ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("data.table", "tidyverse", "ggplot2", "stringr", "PredictABEL", "Epi", "ggrepel", "ivreg", "OneSampleMR", "TwoSampleMR" )
ipak(packages)



Datafile <- fread("DAR_2022_00374_CKB_CC4D_extra_data_convertr.csv" , header = TRUE)
Dosagefile <- fread("DAR-2022-00374.DAR-2022-00374_dosage.csv" , header = TRUE)
quest <- fread("data_baseline_questionnaires.csv", header=T)


gendat <- fread("data_gwas_genetics.csv", header=T)

SNPsEUR <- fread("241_EUR.txt" , header=TRUE)
SNPsBBJ <- fread("42_BBJ.txt" , header=TRUE)
liftedpositionsfile <- fread("Lifted_CC4D.txt", header = TRUE)
infotable <- fread("info.table" , header= T)
str(Datafile)
str(Dosagefile)

Datafile <- full_join(Datafile, quest, by = c("csid"="csid"))

Datafile <- full_join(Datafile, gendat, by = c("csid"="csid"))



Datafile <- Datafile[!(is.na(Datafile$mi) & is.na(Datafile$ihd) & is.na(Datafile$mce)),]

Datafile$anycase[Datafile$mi==1 | Datafile$ihd ==1 | Datafile$mce==1] <- 1
Datafile$anycase[Datafile$mi==0 & Datafile$ihd ==0 & Datafile$mce==0] <- 0
table(Datafile$anycase)


Datafile$included[Datafile$anycase ==1 | Datafile$is_in_gwas_population_subset ==1]  <- 1
table(Datafile$included)
Datafile <- Datafile[Datafile$included==1,]

colnames(Datafile)[165:173] <- c("RC_pc1" , "RC_pc2" , "RC_pc3", "RC_pc4" , "RC_pc5" , "RC_pc6" , "RC_pc7" , "RC_pc8" , "RC_pc9")


DosagefileNN <-as.data.frame(Dosagefile)

colnames(DosagefileNN)[2:258] <- gsub('_[A-Z]',"", colnames(Dosagefile)[2:258])

nameseursnp <-SNPsEUR$rsID[SNPsEUR$rsID %in% colnames(DosagefileNN)]
colnamesretain <- c("csid", nameseursnp)

DosagefileNN <- DosagefileNN[,colnamesretain]

Fulltable <- left_join(Datafile , DosagefileNN, by = c("csid"="csid"))


Fulltable <- Fulltable %>% mutate(eversmk = case_when(Fulltable$smoking_category>2 ~1, Fulltable$smoking_category<=2 ~0))

Fulltable$agesquared <- Fulltable$age_at_study_date_x100*Fulltable$age_at_study_date_x100




Eurweights <- SNPsEUR$Beta[SNPsEUR$rsID %in% colnames(DosagefileNN)]
names(Eurweights) <- SNPsEUR$rsID[SNPsEUR$rsID %in% colnames(DosagefileNN)]

colnames(SNPsBBJ)[10] <- "BetaC"

Jweights <- SNPsBBJ$BetaC[SNPsBBJ$rsID %in% colnames(DosagefileNN)]
names(Jweights) <- SNPsBBJ$rsID[SNPsBBJ$rsID %in% colnames(DosagefileNN)]

weights <- c(Eurweights, Jweights)

Fulltable$GRS <- riskScore(weights = Eurweights, data= Fulltable , cGenPreds=c(178:401), Type= "weighted")


quantv <- quantile(Fulltable$GRS, probs = seq(0,1,1/5))

Fulltable <- Fulltable %>% mutate(GRSQ = case_when(Fulltable$GRS<quantv[2] ~1, Fulltable$GRS<quantv[3] & Fulltable$GRS>=quantv[2] ~2, Fulltable$GRS<quantv[4] & Fulltable$GRS>=quantv[3] ~3, Fulltable$GRS<quantv[5] & Fulltable$GRS>=quantv[4] ~4,Fulltable$GRS>=quantv[5]~5))

FullmodelsMCE <- data.frame()
TFulltable <- as.data.frame(Fulltable)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmce <- as.data.frame(table(X$model$mce))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(mce)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmce[2,2], Controls =countsmce[1,2])
X$Region <- Regions[z]
FullmodelsMCE <- bind_rows(FullmodelsMCE, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$mce==1 & !is.na(countsGRSQ$mce),3], Controls=countsGRSQ[countsGRSQ$mce==0 & !is.na(countsGRSQ$mce),3]) 
floatres <- bind_rows(floatres,fl)
}


FullmodelsMCE$factor <-rownames(FullmodelsMCE)

GRSQresMCE <- FullmodelsMCE[FullmodelsMCE$factor %like% "GRSQ",]





metafunction <- function(b.est, se, cases=NULL, controls=NULL, regions=NULL, means = NA){
  #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
  b.F = sum(b.est / se^2) / sum(1 / se^2)
  se.F = 1 / sqrt(sum(1 / se^2))
  p.F = pchisq( (b.F / se.F)^2, df = 1, lower = F)
  sumcontrols= sum(controls)
  sumcases= sum(cases)
  nregions = length(regions[!is.na(regions)])
  means= sum(means*cases)/sum(cases)
  return(data.frame(beta = b.F , SE= se.F , OR = exp(b.F), CIL=exp(b.F-1.96*se.F)  , CIU = exp(b.F+1.96*se.F), P= p.F , Cases = sumcases , Controls = sumcontrols, Nregions = nregions, mean= means))
} 





GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares <- data.frame()

for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresMCE[GRSQresMCE$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares <- bind_rows(GRSQmares, SNPindres)
}


levs <- unique(floatres$level)


floatresmares <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres[floatres$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares <- bind_rows(floatresmares, SNPindres)
}
fwrite(floatresmares,"floatresmaresMCE.txt" , sep = "\t" , quote=FALSE, row.names=FALSE)


FullmodelsBMI <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
TFulltable <- TFulltable[!is.na(TFulltable$bmi_calc),]
Regions <- unique(TFulltable$region_code)
floatresBMI <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){

X <- lm(bmi_calc ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
fl <- float(X) 
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(bmi_calc))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsBMI <- bind_rows(FullmodelsBMI, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresBMI <- bind_rows(floatresBMI,fl)
}

FullmodelsBMI$factor <-rownames(FullmodelsBMI)


GRSQresBMI <- FullmodelsBMI[FullmodelsBMI$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresBMI <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresBMI[GRSQresBMI$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"],  )
SNPindres$level <- GRSQlev[z]
GRSQmaresBMI<- bind_rows(GRSQmaresBMI, SNPindres)
}


levs <- unique(floatresBMI$level)


floatresmaresBMI <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresBMI[floatresBMI$level==levs[z],]

SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"], means = SNPresdfu[,"mean"], cases= SNPresdfu[,"Counts"] )
SNPindres$level <- levs[z]
floatresmaresBMI <- bind_rows(floatresmaresBMI, SNPindres)
}



FullmodelsSBP <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatresSBP <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(sbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
fl <- float(X) 
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(sbp_mean))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsSBP <- bind_rows(FullmodelsSBP, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresSBP <- bind_rows(floatresSBP,fl)
}


FullmodelsSBP$factor <-rownames(FullmodelsSBP)


GRSQresSBP <- FullmodelsSBP[FullmodelsSBP$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresSBP <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresSBP[GRSQresSBP$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresSBP<- bind_rows(GRSQmaresSBP, SNPindres)
}


levs <- unique(floatresSBP$level)


floatresmaresSBP <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresSBP[floatresSBP$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]) , regions= SNPresdfu[,"region"], cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresSBP <- bind_rows(floatresmaresSBP, SNPindres)
}





FullmodelsDBP <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatresDBP <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(dbp_mean ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(dbp_mean))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsDBP <- bind_rows(FullmodelsDBP, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresDBP <- bind_rows(floatresDBP,fl)
}


FullmodelsDBP$factor <-rownames(FullmodelsDBP)



GRSQresDBP <- FullmodelsDBP[FullmodelsDBP$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresDBP <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresDBP[GRSQresDBP$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresDBP<- bind_rows(GRSQmaresDBP, SNPindres)
}


levs <- unique(floatresDBP$level)


floatresmaresDBP <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresDBP[floatresDBP$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"] , cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresDBP <- bind_rows(floatresmaresDBP, SNPindres)
}






GRSQmaresBMI$OR <- NULL

GRSQmaresBMI$Controls <- NULL
GRSQmaresBMI$CIL <- GRSQmaresBMI$beta-1.96*GRSQmaresBMI$SE
GRSQmaresBMI$CIU <- GRSQmaresBMI$beta+1.96*GRSQmaresBMI$SE

floatresmaresBMI$OR <- NULL

floatresmaresBMI$Controls <- NULL
floatresmaresBMI$CIL <- floatresmaresBMI$beta-1.96*floatresmaresBMI$SE
floatresmaresBMI$CIU <- floatresmaresBMI$beta+1.96*floatresmaresBMI$SE



GRSQmaresSBP$OR <- NULL

GRSQmaresSBP$Controls <- NULL
GRSQmaresSBP$CIL <- GRSQmaresSBP$beta-1.96*GRSQmaresSBP$SE
GRSQmaresSBP$CIU <- GRSQmaresSBP$beta+1.96*GRSQmaresSBP$SE

floatresmaresSBP$OR <- NULL

floatresmaresSBP$Controls <- NULL
floatresmaresSBP$CIL <- floatresmaresSBP$beta-1.96*floatresmaresSBP$SE
floatresmaresSBP$CIU <- floatresmaresSBP$beta+1.96*floatresmaresSBP$SE

GRSQmaresDBP$OR <- NULL

GRSQmaresDBP$Controls <- NULL
GRSQmaresDBP$CIL <- GRSQmaresDBP$beta-1.96*GRSQmaresDBP$SE
GRSQmaresDBP$CIU <- GRSQmaresDBP$beta+1.96*GRSQmaresDBP$SE

floatresmaresDBP$OR <- NULL

floatresmaresDBP$Controls <- NULL
floatresmaresDBP$CIL <- floatresmaresDBP$beta-1.96*floatresmaresDBP$SE
floatresmaresDBP$CIU <- floatresmaresDBP$beta+1.96*floatresmaresDBP$SE


GRSQmaresDBP$OR <- NULL

GRSQmaresDBP$Controls <- NULL
GRSQmaresDBP$CIL <- GRSQmaresDBP$beta-1.96*GRSQmaresDBP$SE
GRSQmaresDBP$CIU <- GRSQmaresDBP$beta+1.96*GRSQmaresDBP$SE

floatresmaresDBP$OR <- NULL

floatresmaresDBP$Controls <- NULL
floatresmaresDBP$CIL <- floatresmaresDBP$beta-1.96*floatresmaresDBP$SE
floatresmaresDBP$CIU <- floatresmaresDBP$beta+1.96*floatresmaresDBP$SE



write.csv(floatresmares,"floatresmares.csv", row.names=FALSE)
write.csv(GRSQmares,"GRSQmares.csv", row.names=FALSE)
write.csv(floatresmaresBMI, "floatresmares_BMI.csv", row.names=FALSE)
write.csv(floatresmaresSBP, "floatresmares_SBP.csv", row.names=FALSE)
write.csv(floatresmaresDBP, "floatresmares_DBP.csv", row.names=FALSE)
write.csv(GRSQmaresBMI, "GRSQmares_BMI.csv", row.names=FALSE)
write.csv(GRSQmaresSBP, "GRSQmares_SBP.csv", row.names=FALSE)
write.csv(GRSQmaresDBP, "GRSQmares_DBP.csv", row.names=FALSE)




FullmodelsMI <- data.frame()
TFulltable <- as.data.frame(Fulltable)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatresMI <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mi ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmi <- as.data.frame(table(X$model$mi))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(mi)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmi[2,2], Controls =countsmi[1,2])
X$Region <- Regions[z]
FullmodelsMI <- bind_rows(FullmodelsMI, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$mi==1 & !is.na(countsGRSQ$mi),3], Controls=countsGRSQ[countsGRSQ$mi==0 & !is.na(countsGRSQ$mi),3]) 
floatresMI <- bind_rows(floatresMI,fl)
}




FullmodelsMI$factor <-rownames(FullmodelsMI)

GRSQresMI <- FullmodelsMI[FullmodelsMI$factor %like% "GRSQ",]

GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresMI <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresMI[GRSQresMI$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresMI <- bind_rows(GRSQmaresMI, SNPindres)
}


levs <- unique(floatresMI$level)


floatresmaresMI <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatresMI[floatresMI$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmaresMI <- bind_rows(floatresmaresMI, SNPindres)
}



fwrite(GRSQmaresMI,"GRSQmaresMI.txt")
fwrite(floatresmaresMI,"floatresmaresMI.txt")





FullmodelsIHD <- data.frame()
TFulltable <- as.data.frame(Fulltable)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatresIHD <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(ihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsihd <- as.data.frame(table(X$model$ihd))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(ihd)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsihd[2,2], Controls =countsihd[1,2])
X$Region <- Regions[z]
FullmodelsIHD <- bind_rows(FullmodelsIHD, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$ihd==1 & !is.na(countsGRSQ$ihd),3], Controls=countsGRSQ[countsGRSQ$ihd==0 & !is.na(countsGRSQ$ihd),3]) 
floatresIHD <- bind_rows(floatresIHD,fl)
}




FullmodelsIHD$factor <-rownames(FullmodelsIHD)


GRSQresIHD <- FullmodelsIHD[FullmodelsIHD$factor %like% "GRSQ",]

GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresIHD <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresIHD[GRSQresIHD$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresIHD <- bind_rows(GRSQmaresIHD, SNPindres)
}


levs <- unique(floatresIHD$level)


floatresmaresIHD <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatresIHD[floatresIHD$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmaresIHD <- bind_rows(floatresmaresIHD, SNPindres)
}

fwrite(GRSQmaresIHD,"GRSQmaresIHD.txt")
fwrite(floatresmaresIHD,"floatresmaresIHD.txt")





phecodes1 <- fread("data_any_phecodes_v1_01.csv")
phecodes2 <- fread("data_any_phecodes_v1_02.csv")


phecodesparents <- fread("data_any_phecodes_v1_parents.csv", header=T)



meanGRS <- mean(Fulltable$GRS,na.rm=T)
sdGRS <- sd(Fulltable$GRS,na.rm=T)
Fulltable$GRS_std <- (Fulltable$GRS-meanGRS)/sdGRS

Allphecodes <- left_join(phecodes1, phecodes2, by = c("csid"="csid"))
Fulltable_phe <- left_join(Fulltable, Allphecodes, by=c("csid"="csid"))


table(colSums(Allphecodes[,2:length(Allphecodes)])>200)

table(colSums(phecodesparents[,2:length(phecodesparents)])>100)


colnames(Fulltable_phe)[154:164] <- paste0("PC",1:11)


phewasres<- data.frame()


Fulltable_phe <- as.data.frame(Fulltable_phe)
for (i in 408: length(colnames(Fulltable_phe))){
phecode <- colnames(Fulltable_phe)[i]
if(sum(Fulltable_phe[,i])>200){
X <- glm( Fulltable_phe[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11, data = Fulltable_phe, family = "binomial", subset= is_in_gwas_population_subset==1)
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$phecode <- phecode
phewasres<- bind_rows(phewasres,X)}
}

phewasres$factor <- rownames(phewasres)

phewasres_GRS <- phewasres[phewasres$factor %like% "GRS_std",]
colnames(phewasres_GRS)[4] <- "P"

phewasres_GRS$fdr <- p.adjust(phewasres_GRS$P, method="fdr")

phecodesmap <- fread("phecode_definitions1.2.csv", header= T )


phecodes_names_all <- data.frame(phecode=substring(colnames(Fulltable_phe)[408: length(colnames(Fulltable_phe))],5))  

phecodes_names_all$phecodenumeric <- sub("_",".",phecodes_names_all$phecode)
phecodes_names_all$phecodenumeric <- as.numeric(phecodes_names_all$phecodenumeric)

phecodesmapjoined <- left_join(phecodes_names_all, phecodesmaps , by = c("phecodenumeric"="phecode"))


phewasres_GRS <- phewasres_GRS %>% mutate(direction = case_when(phewasres_GRS$Estimate>0 ~"Positive", phewasres_GRS$Estimate<0~"Negative"))
phewasres_GRS$direction <- as.factor(phewasres_GRS$direction)

colnames(phewasres_GRS)[2] <-"SE"

phewasres_GRS <- phewasres_GRS[!phewasres_GRS$category== "NULL",]

phewasres_GRS$phecodenumeric <- sub("_",".", substring(phewasres_GRS$phecode,5))
phewasres_GRS$phecodenumeric <- as.numeric(phewasres_GRS$phecodenumeric)
phewasres_GRS <-  left_join(phewasres_GRS, phecodesmapjoined , by = c("phecodenumeric"="phecodenumeric"))

axisset <- phewasres_GRS %>% group_by(category) %>% summarize(count=n())

pdf("Phewas_MP_final.pdf", width=30, height=15)
ggplot(phewasres_GRS, aes(x=phecode.x, y=-log(P))) + geom_point(aes(col=category, shape=direction, size=abs(Estimate), fill=category))+scale_shape_manual(values=c(25,17)) + theme_classic() + theme(axis.text.x = element_blank(),  axis.ticks=element_blank()) + labs(x="phenotype",color="Category", size="Beta", y="log(p-value)") +geom_text_repel(data=. %>% mutate(label = ifelse(P < 0.05, as.character(phenotype), "")), aes(label=label), size=4.1, box.padding = unit(0.7, "lines")) + geom_hline(yintercept=-log(0.001), color="blue", size=1, alpha=0.5)+geom_hline(yintercept=-log(0.0004273504), color="red", size=1, alpha=0.5) + guides(fill="none") 
dev.off()

fwrite(phewasres, "phewas_GRS.txt")



olinkmeta <- fread("upd_rel18_1_data/olink_meta.csv", header=T)
olinkbaseline <- fread("data_baseline_olink.csv", header=T)  

olinkcardio <- fread("data_baseline_olink_cardiometabolic.csv", header=T)
olinkinf <- fread("data_baseline_olink_inflammation.csv", header=T)
olinkneuro <- fread("data_baseline_olink_neurology.csv",header=T)
olinkonco <- fread("data_baseline_olink_oncology.csv",header=T)


suppolink <- fread("DAR-2022-00374.olink_supp.csv", header=T)


duplicatedids <- c(colnames(olinkcardio), colnames(olinkinf), colnames(olinkneuro), colnames(olinkonco))[duplicated(c(colnames(olinkcardio), colnames(olinkinf), colnames(olinkneuro), colnames(olinkonco)))==T]




olinkall <- left_join(olinkcardio, olinkinf, by = c("csid"="csid"))
olinkall <- left_join(olinkall, olinkneuro, by = c("csid"="csid"))
olinkall <- left_join(olinkall, olinkonco, by = c("csid"="csid"))




suppolink <- suppolink[! suppolink$file=="for_import_olink_20190664_Clarke_NPX_2020-10-29.txt",]

suppolink$batch[suppolink$file=="for_import_olink_20180858_Chen_NPX_2021-12-28.txt"] <- 1
suppolink$batch[suppolink$file=="for_import_olink_Q-00533_Chen_NPX_2022-08-01.txt"] <- 2





coldetails <- fread("column_details.csv", header=T)


coldetails <- coldetails[coldetails$table %in% c("data_baseline_olink_cardiometabolic","data_baseline_olink_inflammation", "data_baseline_olink_neurology", "data_baseline_olink_oncology") ,]

for(i in 1: length(coldetails$notes)){
coldetails$panel_lot_nr[i] <- str_split(coldetails$notes[i], ",") [1][[1]][1]
}

coldetails$panel_lot_nr <- substring(coldetails$panel_lot_nr, 11)

coldetails <- coldetails[!coldetails$panel_lot_nr=="ion=1",]

coldetails$batch[coldetails$panel_lot_nr %in% c("B04411","B04412","B04413","B04414")] <- "1"
coldetails$batch[coldetails$panel_lot_nr %in% c("B20704","B20705","B20706","B20707")] <- "2"





Nucleus_olink <- suppolink %>% group_by(csid,plateid,batch, panel_lot_nr) %>% summarize(count=n())

protdetails <- coldetails %>% group_by(column_name, panel_lot_nr,batch) %>% summarize(count=n()) 


Bo4411ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B04411"] 
Bo4411ns[1] <- TRUE
olink_panelBo4411 <- olinkall[,..Bo4411ns]

Bo4412ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B04412"] 
Bo4412ns[1] <- TRUE
olink_panelBo4412 <- olinkall[,..Bo4412ns]

Bo4413ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B04413"] 
Bo4413ns[1] <- TRUE
olink_panelBo4413 <- olinkall[,..Bo4413ns]

Bo4414ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B04414"] 
Bo4414ns[1] <- TRUE
olink_panelBo4414 <- olinkall[,..Bo4414ns]


B20704ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B20704"] 
B20704ns[1] <- TRUE
olink_panelB20704 <- olinkall[,..B20704ns]

B20705ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B20705"] 
B20705ns[1] <- TRUE
olink_panelB20705 <- olinkall[,..B20705ns]

B20706ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B20706"] 
B20706ns[1] <- TRUE
olink_panelB20706 <- olinkall[,..B20706ns]

B20707ns <-colnames(olinkall) %in% protdetails$column_name[protdetails$panel_lot_nr=="B20707"] 
B20707ns[1] <- TRUE
olink_panelB20707 <- olinkall[,..B20707ns]

ascertainmentfile <- fread("data_baseline_ascertainments.csv", header=T)


Fulltable_prot <- left_join(Fulltable,ascertainmentfile, by = c("csid"="csid"))


colnames(Fulltable_prot)[154:164] <- paste0("PC",1:11)

 
Fulltable_prot_panel_Bo4411 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B04411",], by= c("csid"="csid"))
Fulltable_prot_panel_Bo4411 <- left_join(olink_panelBo4411, Fulltable_prot_panel_Bo4411 ,by=c("csid"="csid"))

GRSres_prot_pan1 <- data.frame()

Fulltable_prot_panel_Bo4411 <- as.data.frame(Fulltable_prot_panel_Bo4411)
for (i in 2:364){
X <- glm( Fulltable_prot_panel_Bo4411[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort+plateid , data = Fulltable_prot_panel_Bo4411, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4411)[i]
GRSres_prot_pan1 <- bind_rows(GRSres_prot_pan1,X) 
}
GRSres_prot_pan1$factor <- rownames(GRSres_prot_pan1)
Fulltable_prot_panel_Bo4412 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B04412",], by= c("csid"="csid"))
Fulltable_prot_panel_Bo4412 <- left_join(olink_panelBo4412, Fulltable_prot_panel_Bo4412 ,by=c("csid"="csid"))

GRSres_prot_pan2 <- data.frame()

Fulltable_prot_panel_Bo4412 <- as.data.frame(Fulltable_prot_panel_Bo4412)
Fulltable_prot_panel_Bo4412 <- Fulltable_prot_panel_Bo4412[!is.na(Fulltable_prot_panel_Bo4412$eversmk),]
for (i in 2:366){
X <- glm( Fulltable_prot_panel_Bo4412[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort+plateid , data = Fulltable_prot_panel_Bo4412, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4412)[i]
GRSres_prot_pan2 <- bind_rows(GRSres_prot_pan2,X) 
}
GRSres_prot_pan2$factor <- rownames(GRSres_prot_pan2)

Fulltable_prot_panel_Bo4413 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B04413",], by= c("csid"="csid"))
Fulltable_prot_panel_Bo4413 <- left_join(olink_panelBo4413, Fulltable_prot_panel_Bo4413 ,by=c("csid"="csid"))


GRSres_prot_pan3 <- data.frame()

Fulltable_prot_panel_Bo4413 <- as.data.frame(Fulltable_prot_panel_Bo4413)
for (i in 2:367){
X <- glm( Fulltable_prot_panel_Bo4413[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_Bo4413, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4413)[i]
GRSres_prot_pan3 <- bind_rows(GRSres_prot_pan3,X) 
}



GRSres_prot_pan3$factor <- rownames(GRSres_prot_pan3)

Fulltable_prot_panel_Bo4414 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B04414",], by= c("csid"="csid"))
Fulltable_prot_panel_Bo4414 <- left_join(olink_panelBo4414, Fulltable_prot_panel_Bo4414 ,by=c("csid"="csid"))

Fulltable_prot_panel_Bo4414$ol_casc4 <- NULL
Fulltable_prot_panel_Bo4414$ol_havcr2 <- NULL
Fulltable_prot_panel_Bo4414$ol_wisp2 <- NULL
Fulltable_prot_panel_Bo4414$ol_hars <- NULL
GRSres_prot_pan4 <- data.frame()

Fulltable_prot_panel_Bo4414 <- as.data.frame(Fulltable_prot_panel_Bo4414)
for (i in 2:365){
X <- glm( Fulltable_prot_panel_Bo4414[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_Bo4414, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4414)[i]
GRSres_prot_pan4 <- bind_rows(GRSres_prot_pan4,X) 
}
GRSres_prot_pan4$factor <- rownames(GRSres_prot_pan4)


Fulltable_prot_panel_B20704 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B20704",], by= c("csid"="csid"))
Fulltable_prot_panel_B20704 <- left_join(olink_panelB20704, Fulltable_prot_panel_B20704 ,by=c("csid"="csid"))


GRSres_prot_pan2_1 <- data.frame()

Fulltable_prot_panel_B20704 <- as.data.frame(Fulltable_prot_panel_B20704)
for (i in 2:366){
X <- glm( Fulltable_prot_panel_B20704[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20704, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20704)[i]
GRSres_prot_pan2_1 <- bind_rows(GRSres_prot_pan2_1,X) 
}
GRSres_prot_pan2_1$factor <- rownames(GRSres_prot_pan2_1)

Fulltable_prot_panel_B20705 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B20705",], by= c("csid"="csid"))
Fulltable_prot_panel_B20705 <- left_join(olink_panelB20705, Fulltable_prot_panel_B20705 ,by=c("csid"="csid"))



GRSres_prot_pan2_2 <- data.frame()

Fulltable_prot_panel_B20705 <- as.data.frame(Fulltable_prot_panel_B20705)
for (i in 2:366){
X <- glm( Fulltable_prot_panel_B20705[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20705, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20705)[i]
GRSres_prot_pan2_2 <- bind_rows(GRSres_prot_pan2_2,X) 
}
GRSres_prot_pan2_2$factor <- rownames(GRSres_prot_pan2_2)



Fulltable_prot_panel_B20706 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B20706",], by= c("csid"="csid"))
Fulltable_prot_panel_B20706 <- left_join(olink_panelB20706, Fulltable_prot_panel_B20706 ,by=c("csid"="csid"))



GRSres_prot_pan2_3 <- data.frame()

Fulltable_prot_panel_B20706 <- as.data.frame(Fulltable_prot_panel_B20706)
for (i in 2:364){
X <- glm( Fulltable_prot_panel_B20706[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20706, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20706)[i]
GRSres_prot_pan2_3 <- bind_rows(GRSres_prot_pan2_3,X) 
}
GRSres_prot_pan2_3$factor <- rownames(GRSres_prot_pan2_3)




Fulltable_prot_panel_B20707 <- left_join(Fulltable_prot, Nucleus_olink[Nucleus_olink$panel_lot_nr=="B20707",], by= c("csid"="csid"))
Fulltable_prot_panel_B20707 <- left_join(olink_panelB20707, Fulltable_prot_panel_B20707 ,by=c("csid"="csid"))



GRSres_prot_pan2_4 <- data.frame()

Fulltable_prot_panel_B20707 <- as.data.frame(Fulltable_prot_panel_B20707)
for (i in 2:363){
X <- glm( Fulltable_prot_panel_B20707[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20707, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20706)[i]
GRSres_prot_pan2_4 <- bind_rows(GRSres_prot_pan2_4,X) 
}
GRSres_prot_pan2_4$factor <- rownames(GRSres_prot_pan2_4)



GRSres_allpanels <- bind_rows(list(GRSres_prot_pan1,GRSres_prot_pan2,GRSres_prot_pan3,GRSres_prot_pan4,GRSres_prot_pan2_1,GRSres_prot_pan2_2,GRSres_prot_pan2_3,GRSres_prot_pan2_4)) 

GRSres_allpanels_GRS <- GRSres_allpanels[GRSres_allpanels$factor %like% "GRS_std",]

colnames(GRSres_allpanels_GRS)[4] <- "P"

GRSres_allpanels_GRS$fdr <- p.adjust(GRSres_allpanels_GRS$P, method="fdr") 
GRSres_allpanels_GRS$direction[GRSres_allpanels_GRS$Estimate>0] <-"Positive"
GRSres_allpanels_GRS$direction[GRSres_allpanels_GRS$Estimate<0] <-"Negative"

colnames(GRSres_allpanels_GRS)[2] <-"SE"

protein_positions <- fread("protein_pos.txt", header=T)
protein_positions <- protein_positions[-1,]
GRSres_allpanels_GRS <- left_join(GRSres_allpanels_GRS,protein_positions, by = c("prot"="Protein"))

colnames(GRSres_allpanels_GRS)[19] <- "max_end_position"
colnames(GRSres_allpanels_GRS)[18] <- "min_start_position"


GRSres_allpanels_GRS <- GRSres_allpanels_GRS[!is.na(GRSres_allpanels_GRS$chromosome_name),]

data_cum <- GRSres_allpanels_GRS %>% 
  group_by(chromosome_name) %>% 
  summarise(max_bp = max(min_start_position)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chromosome_name, bp_add)
  
GRSres_allpanels_GRS_plot <- GRSres_allpanels_GRS %>% 
  inner_join(data_cum, by = "chromosome_name") %>% 
  mutate(bp_cum = min_start_position + bp_add)

axis_set <- GRSres_allpanels_GRS_plot %>% 
  group_by(chromosome_name) %>% 
  summarize(center = mean(bp_cum))


pdf("Proteo_MP.pdf",width=20, height=10)
ggplot(GRSres_allpanels_GRS_plot, aes(x=bp_cum, y=-log(P))) + geom_point(aes(shape=direction, size=Estimate*10, color=as.factor(chromosome_name) , fill=as.factor(chromosome_name)), alpha=0.75)+ scale_shape_manual(values=c(25,17)) + theme_classic() + scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center)+ scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chromosome_name)))) + scale_fill_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chromosome_name))))+
labs(x="Chromsosome", size="Inverse SE", y="log(p-value)") +geom_text_repel(data=. %>% mutate(label = ifelse(fdr < 0.05, as.character(Protein_hcng), "")), aes(label=label), size=4.1, box.padding = unit(0.7, "lines")) + geom_hline(yintercept=-log(0.000017), color="red", size=1, alpha=0.5) + guides(fill="none",color="none")+theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 15, vjust = 0.5)
  )
dev.off()






biosamples <- fread("data_baseline_biosamples.csv", header=T, select= c("csid","ldl_mmoll", "ldl_invalid", "ldl_below_range", "ldl_note"))


Fulltable_bio <- left_join(Fulltable, biosamples, by = c("csid"="csid")) 



Fullmodelsldl <- data.frame()
TFulltable <- as.data.frame(Fulltable_bio[Fulltable_bio$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatresldl <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(ldl_mmoll ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(ldl_mmoll))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
Fullmodelsldl <- bind_rows(Fullmodelsldl, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresldl <- bind_rows(floatresldl,fl)
}



Fullmodelsldl$factor <-rownames(Fullmodelsldl)

GRSQresldl <- Fullmodelsldl[Fullmodelsldl$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresldl <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresldl[GRSQresldl$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresldl<- bind_rows(GRSQmaresldl, SNPindres)
}


levs <- unique(floatresldl$level)


floatresmaresldl <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresldl[floatresldl$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"] , cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresldl <- bind_rows(floatresmaresldl, SNPindres)
}




phewasres_regionadjusted<- data.frame()

Fulltable_phe <- as.data.frame(Fulltable_phe)
for (i in 408: length(colnames(Fulltable_phe))){
phecode <- colnames(Fulltable_phe)[i]
if(sum(Fulltable_phe[,i])>200){
X <- glm( Fulltable_phe[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + as.factor(region_code)+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11, data = Fulltable_phe, family = "binomial", subset= is_in_gwas_population_subset==1)
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$phecode <- phecode
phewasres_regionadjusted<- bind_rows(phewasres_regionadjusted,X)}
}

phewasres_regionadjusted$factor <- rownames(phewasres_regionadjusted)

phewasres_GRS_regadjusted <- phewasres_regionadjusted[phewasres_regionadjusted$factor %like% "GRS_std",]
colnames(phewasres_GRS_regadjusted)[4] <- "P"

phewasres_GRS_regadjusted$fdr <- p.adjust(phewasres_GRS_regadjusted$P, method="fdr")




phewasres_GRS_regadjusted <- phewasres_GRS_regadjusted %>% mutate(direction = case_when(phewasres_GRS_regadjusted$Estimate>0 ~"Positive", phewasres_GRS_regadjusted$Estimate<0~"Negative"))
phewasres_GRS_regadjusted$direction <- as.factor(phewasres_GRS_regadjusted$direction)

colnames(phewasres_GRS_regadjusted)[2] <-"SE"



phewasres_GRS_regadjusted$phecodenumeric <- sub("_",".", substring(phewasres_GRS_regadjusted$phecode,5))
phewasres_GRS_regadjusted$phecodenumeric <- as.numeric(phewasres_GRS_regadjusted$phecodenumeric)
phewasres_GRS_regadjusted <-  left_join(phewasres_GRS_regadjusted, phecodesmapjoined , by = c("phecodenumeric"="phecodenumeric"))

phewasres_GRS_regadjusted <- phewasres_GRS_regadjusted[!phewasres_GRS_regadjusted$category== "NULL",]

axisset <- phewasres_GRS_regadjusted %>% group_by(category) %>% summarize(count=n())


phewasres_GRS_regadjusted <- phewasres_GRS_regadjusted[order(phewasres_GRS_regadjusted$category),]


phewasres_GRS_regadjusted <- left_join(phewasres_GRS_regadjusted, axisset, by = c("category"="category"))

phewasres_GRS_regadjusted$nchar_cat <- nchar(phewasres_GRS_regadjusted$category)

phewasres_GRS_regadjusted <- phewasres_GRS_regadjusted[order(phewasres_GRS_regadjusted$count, phewasres_GRS_regadjusted$nchar_cat),]



phewasres_GRS_regadjusted$unique_position <- -log(c(1:length(phewasres_GRS_regadjusted$phecode.x)))





axis_set <- phewasres_GRS_regadjusted %>% 
  group_by(category) %>% 
  summarize(center = mean(unique_position))


camelCase = function(sv, upper=FALSE, capIsNew=FALSE, alreadyTrimmed=FALSE) {
  if (!is.character(sv)) stop("'sv' must be a string vector")
  if (!alreadyTrimmed) sv = gsub("[[:space:]]*$", "", gsub("^[[:space:]]*", "", sv))
  if (capIsNew) {
    sv = gsub("([A-Z])", " \\1", sv)
    sv = gsub("^[[:space:]]", "", sv)
    sv = tolower(sv)
  }
  apart = strsplit(sv, split="[[:space:][:punct:]]")
  apart = lapply(apart, tolower)
  capitalize = function(x){paste0(toupper(substring(x,1,1)), substring(x,2))}
  if (upper) {
    apart = lapply(apart, capitalize)
  } else {
    apart = lapply(apart, function(x) c(x[1], capitalize(x[-1])))
  }
  return(sapply(apart, paste, collapse=""))
}

phewasres_GRS_regadjusted$category <- camelCase(phewasres_GRS_regadjusted$category, upper=T, alreadyTrimmed=T, capIsNew=FALSE)



pdf("Phewas_MP_final_regadj_27Feb.pdf", width=30, height=15)
ggplot(phewasres_GRS_regadjusted, aes(x=unique_position, y=-log10(P))) + geom_point(aes(col=category, shape=direction, size=5*exp(Estimate), fill=category))+scale_shape_manual(values=c(25,17)) + theme_classic() + theme(axis.text.x = element_text(angle=90, vjust=0.5 ), text= element_text(size=25), legend.position = "right" , legend.box.background = element_rect(colour = "black", size = rel(2) , linewidth = 5)) + labs(x=NULL,color="Category", shape="Direction",size="OR", y="log(p-value)") +geom_text_repel(data=. %>% mutate(label = ifelse(P < 0.004273504, as.character(phenotype), "")), aes(label=label), size=7, box.padding = unit(0.7, "lines")) + 
geom_hline(yintercept=-log10(0.0004273504), color="red", size=1, alpha=0.5) + guides(fill="none", color="none", shape = guide_legend(override.aes = list(size = 4) ), size= guide_legend(override.aes=list(shape=c(17)))) +scale_x_continuous(label = axis_set$category, breaks = axis_set$center)  
dev.off()





 

GRSres_prot_pan1_regadj <- data.frame()

Fulltable_prot_panel_Bo4411 <- as.data.frame(Fulltable_prot_panel_Bo4411)
for (i in 2:364){
X <- glm( Fulltable_prot_panel_Bo4411[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code+PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort+plateid , data = Fulltable_prot_panel_Bo4411, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4411)[i]
GRSres_prot_pan1_regadj <- bind_rows(GRSres_prot_pan1_regadj,X) 
}
GRSres_prot_pan1_regadj$factor <- rownames(GRSres_prot_pan1_regadj)


GRSres_prot_pan2_regadj <- data.frame()

for (i in 2:366){
X <- glm( Fulltable_prot_panel_Bo4412[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort+plateid , data = Fulltable_prot_panel_Bo4412, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4412)[i]
GRSres_prot_pan2_regadj <- bind_rows(GRSres_prot_pan2_regadj,X) 
}
GRSres_prot_pan2_regadj$factor <- rownames(GRSres_prot_pan2_regadj)





GRSres_prot_pan3_regadj <- data.frame()

Fulltable_prot_panel_Bo4413 <- as.data.frame(Fulltable_prot_panel_Bo4413)
for (i in 2:367){
X <- glm( Fulltable_prot_panel_Bo4413[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_Bo4413, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4413)[i]
GRSres_prot_pan3_regadj <- bind_rows(GRSres_prot_pan3_regadj,X) 
}
GRSres_prot_pan3_regadj$factor <- rownames(GRSres_prot_pan3_regadj)

GRSres_prot_pan4_regadj <- data.frame()

Fulltable_prot_panel_Bo4414 <- as.data.frame(Fulltable_prot_panel_Bo4414)
for (i in 2:365){
X <- glm( Fulltable_prot_panel_Bo4414[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_Bo4414, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_Bo4414)[i]
GRSres_prot_pan4_regadj <- bind_rows(GRSres_prot_pan4_regadj,X) 
}
GRSres_prot_pan4_regadj$factor <- rownames(GRSres_prot_pan4_regadj)


GRSres_prot_pan2_1_regadj <- data.frame()

Fulltable_prot_panel_B20704 <- as.data.frame(Fulltable_prot_panel_B20704)
for (i in 2:366){
X <- glm( Fulltable_prot_panel_B20704[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20704, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20704)[i]
GRSres_prot_pan2_1_regadj <- bind_rows(GRSres_prot_pan2_1_regadj,X) 
}
GRSres_prot_pan2_1_regadj$factor <- rownames(GRSres_prot_pan2_1_regadj)



GRSres_prot_pan2_2_regadj <- data.frame()

Fulltable_prot_panel_B20705 <- as.data.frame(Fulltable_prot_panel_B20705)
for (i in 2:366){
X <- glm( Fulltable_prot_panel_B20705[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20705, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20705)[i]
GRSres_prot_pan2_2_regadj <- bind_rows(GRSres_prot_pan2_2_regadj,X) 
}
GRSres_prot_pan2_2_regadj$factor <- rownames(GRSres_prot_pan2_2_regadj)




GRSres_prot_pan2_3_regadj <- data.frame()

Fulltable_prot_panel_B20706 <- as.data.frame(Fulltable_prot_panel_B20706)
for (i in 2:364){
X <- glm( Fulltable_prot_panel_B20706[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + region_code + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20706, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20706)[i]
GRSres_prot_pan2_3_regadj <- bind_rows(GRSres_prot_pan2_3_regadj,X) 
}
GRSres_prot_pan2_3_regadj$factor <- rownames(GRSres_prot_pan2_3_regadj)





GRSres_prot_pan2_4_regadj <- data.frame()

Fulltable_prot_panel_B20707 <- as.data.frame(Fulltable_prot_panel_B20707)
for (i in 2:363){
X <- glm( Fulltable_prot_panel_B20707[,i]~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11+ olinkexp1536_chd_b1_subcohort + plateid , data = Fulltable_prot_panel_B20707, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Fulltable_prot_panel_B20706)[i]
GRSres_prot_pan2_4_regadj <- bind_rows(GRSres_prot_pan2_4_regadj,X) 
}
GRSres_prot_pan2_4_regadj$factor <- rownames(GRSres_prot_pan2_4_regadj)




GRSres_allpanels_regadj <- bind_rows(list(GRSres_prot_pan1_regadj,GRSres_prot_pan2_regadj,GRSres_prot_pan3_regadj,GRSres_prot_pan4_regadj,GRSres_prot_pan2_1_regadj,GRSres_prot_pan2_2_regadj,GRSres_prot_pan2_3_regadj,GRSres_prot_pan2_4_regadj)) 

GRSres_allpanels_GRS_regadj <- GRSres_allpanels_regadj[GRSres_allpanels_regadj$factor %like% "GRS_std",]

colnames(GRSres_allpanels_GRS_regadj)[4] <- "P"

GRSres_allpanels_GRS_regadj$fdr <- p.adjust(GRSres_allpanels_GRS_regadj$P, method="fdr") 
GRSres_allpanels_GRS_regadj$direction[GRSres_allpanels_GRS_regadj$Estimate>0] <-"Positive"
GRSres_allpanels_GRS_regadj$direction[GRSres_allpanels_GRS_regadj$Estimate<0] <-"Negative"

colnames(GRSres_allpanels_GRS_regadj)[2] <-"SE"

GRSres_allpanels_GRS_regadj <- left_join(GRSres_allpanels_GRS,protein_positions, by = c("prot"="Protein"))

colnames(GRSres_allpanels_GRS_regadj)[19] <- "max_end_position"
colnames(GRSres_allpanels_GRS_regadj)[18] <- "min_start_position"
colnames(GRSres_allpanels_GRS_regadj)[15] <- "chromosome_name"
colnames(GRSres_allpanels_GRS_regadj)[20] <- "Protein_hcng"


GRSres_allpanels_GRS_regadj <- GRSres_allpanels_GRS_regadj[!is.na(GRSres_allpanels_GRS_regadj$chromosome_name),]


data_cum <- GRSres_allpanels_GRS_regadj %>% 
  group_by(chromosome_name) %>% 
  summarise(max_bp = max(min_start_position)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chromosome_name, bp_add)
  
GRSres_allpanels_GRS_plot_regadj <- GRSres_allpanels_GRS_regadj %>% 
  inner_join(data_cum, by = "chromosome_name") %>% 
  mutate(bp_cum = min_start_position + bp_add)

axis_set <- GRSres_allpanels_GRS_plot %>% 
  group_by(chromosome_name) %>% 
  summarize(center = mean(bp_cum))


pdf("Proteo_MP_regadj.pdf",width=20, height=10)
ggplot(GRSres_allpanels_GRS_plot_regadj, aes(x=bp_cum, y=-log10(P))) + geom_point(aes(shape=direction, size=Estimate*10, color=as.factor(chromosome_name) , fill=as.factor(chromosome_name)), alpha=0.75)+ scale_shape_manual(values=c(25,17)) + theme_classic() + scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center)+ scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chromosome_name)))) + scale_fill_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chromosome_name))))+
labs(x="Chromsosome", size="Inverse SE", y="log(p-value)") +geom_text_repel(data=. %>% mutate(label = ifelse(fdr < 0.05, as.character(Protein_hcng), "")), aes(label=label), size=4.1, box.padding = unit(0.7, "lines")) + geom_hline(yintercept=-log10(0.000017), color="red", size=1, alpha=0.5) + guides(fill="none",color="none")+theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 15, vjust = 0.5)
  )
dev.off()




pQTLSdosage <- fread("DAR-2023-00152.DAR-2023-00152_dosage.csv", header=T)

pQTLsmap <- fread("PQTLs_CC4D.txt", header=T)

 
pQTLsmap[4,1] <- "rs145598072"

allelesdose <- colnames(pQTLSdosage)
colnames(pQTLSdosage)[2:7] <- gsub('_[A-Z]',"", colnames(pQTLSdosage)[2:7])


AllPQTLs <- fread("cojo_pQTL_info_table_annotated_ver4.csv", header=T)

AllPQTLs <- AllPQTLs[AllPQTLs$SNP %in% pQTLsmap$Markername]


Fulltable_pqtlanalysis <- left_join(Fulltable, pQTLSdosage,by = c("csid"="csid"))



FullmodelsMCE_pQTLs <- data.frame()
TFulltable <- as.data.frame(Fulltable_pqtlanalysis)

Regions <- unique(TFulltable$region_code)
floatres_pqtl <- data.frame()
for(x in 408:413) {
pqtlanalyzed <- colnames(TFulltable)[x]
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mce ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmce <- as.data.frame(table(X$model$mce))

X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmce[2,2], Controls =countsmce[1,2])
X$Region <- Regions[z]
X$pqtl <- pqtlanalyzed
FullmodelsMCE_pQTLs <- bind_rows(FullmodelsMCE_pQTLs, X)

}
}

head(FullmodelsMCE_pQTLs)
str(FullmodelsMCE_pQTLs)
str(floatres_pqtl)

FullmodelsMCE_pQTLs$factor <-rownames(FullmodelsMCE_pQTLs)


pqtlresMCE <- FullmodelsMCE_pQTLs[FullmodelsMCE_pQTLs$factor %like% "TFulltable",]


pqtllev <- colnames(pQTLSdosage)[2:7]
pqtlsmares <- data.frame()




for(z in 1:length(pqtllev)){

SNPresdfu <- pqtlresMCE[pqtlresMCE$pqtl %like% pqtllev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$pqtl <- pqtllev[z]
pqtlsmares <- bind_rows(pqtlsmares, SNPindres)
}

pqtlsmares <- left_join(pqtlsmares,pQTLsmap, by = c("pqtl"="rsID"))


Fulltable_SUSD2MR <- Fulltable_prot_panel_Bo4414[,c(1,323,366:454,518:528, 766:794)]
Fulltable_SUSD2MR <- left_join(Fulltable_SUSD2MR, pQTLSdosage[,c("csid","rs138546470")],by = c("csid"="csid"))


Fulltable_SUSD2MR <- Fulltable_SUSD2MR[!is.na(Fulltable_SUSD2MR$is_female),]

iv_mod_susd2_mce_SNP <- ivreg( mce~ ol_susd2 + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11, data = Fulltable_SUSD2MR)

iv_mod_susd2_mi_SNP <- ivreg( mi~ ol_susd2 + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11, data = Fulltable_SUSD2MR)

iv_mod_susd2_ihd_SNP <- ivreg( ihd~ ol_susd2 + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11, data = Fulltable_SUSD2MR)







iv_mod_mce_susd2_SNP <- ivreg( ol_susd2~ mce + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort, data = Fulltable_SUSD2MR)
iv_mod_mi_susd2_SNP <- ivreg( ol_susd2~ mi + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort, data = Fulltable_SUSD2MR)
iv_mod_ihd_susd2_SNP <- ivreg( ol_susd2~ ihd + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| rs138546470 + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort, data = Fulltable_SUSD2MR)




pqtl_susd2_2MR <- format_data(AllPQTLs, type="exposure", header= TRUE, phenotype_col = "Assay", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EA_FREQ", chr_col = "CHR", pos_col = "BP", effect_allele_col = "EA", other_allele_col = "AA", pval_col="P")


pqtlsmares <- left_join(pqtlsmares, infotable[,c("rsid","A1_freq")], by = c("Markername"="rsid"))
mce_susd2_2MR <- format_data(pqtlsmares, type="outcome", header= TRUE , snp_col = "Markername", beta_col = "beta", se_col = "SE", eaf_col = "A1_freq", chr_col = "Chr", pos_col = "Position", effect_allele_col = "A1", other_allele_col = "A2", pval_col="P" )

mce_susd2_2MR$outcome <- "mce"







mr_data_ad_MCE <- harmonise_data( exposure = pqtl_susd2_2MR, outcome = mce_susd2_2MR)

res_ad_MCE <- mr(mr_data_ad_MCE)


res_ad_single_MCE <- mr_singlesnp(mr_data_ad_MCE)


iv_mod_susd2_mce_GRS <- ivreg( mce~ ol_susd2 + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| GRS_std + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11, data = Fulltable_SUSD2MR)




iv_mod_mce_susd2_GRS <- ivreg( ol_susd2 ~ mce  + is_female+ age_at_study_date_x100  + eversmk+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| GRS_std + is_female+ age_at_study_date_x100 + eversmk +  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11, data = Fulltable_SUSD2MR)

iv_mod_mce_susd2_GRS <- ivreg( ol_susd2 ~ mce + is_female + age_at_study_date_x100 + eversmk + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort| GRS_std + is_female + age_at_study_date_x100 + eversmk + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+olinkexp1536_chd_b1_subcohort, data = Fulltable_SUSD2MR)



mr_data_ad$samplesize.exposure <- 3977
mr_data_ad$samplesize.outcome <- 80702
directionality <- directionality_test(mr_data_ad)


fp_multisnpmr <- forest_plot(res_ad_single[!res_ad_single$exposure=="GGT5",], exponentiate=TRUE, threshold=0.05, trans=log2)






mr_forest_plot <- function(singlesnp_results, exponentiate=FALSE)
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
		levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
		am <- grep("All", d$SNP, value=TRUE)
		d$up <- d$b + 1.96 * d$se
		d$lo <- d$b - 1.96 * d$se
		d$tot <- 0.01
		d$tot[d$SNP %in% am] <- 1
		d$SNP <- as.character(d$SNP)
		nom <- d$SNP[! d$SNP %in% am]
		nom <- nom[order(d$b)]
		d <- rbind(d, d[nrow(d),])
		d$SNP[nrow(d)-1] <- ""
		d$b[nrow(d)-1] <- NA
		d$up[nrow(d)-1] <- NA
		d$lo[nrow(d)-1] <- NA
		d$SNP <- ordered(d$SNP, levels=c(am, "", nom))

		xint <- 0
		if(exponentiate)
		{
			d$b <- exp(d$b)
			d$up <- exp(d$up)
			d$lo <- exp(d$lo)
			xint <- 1
		}

		ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
		ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
		ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
		ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
		ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
		ggplot2::scale_colour_manual(values=c("black", "red")) +
		ggplot2::scale_size_manual(values=c(0.3, 1)) +
		ggplot2::theme(
			legend.position="none", 
			axis.text.y=ggplot2::element_text(size=8), 
			axis.ticks.y=ggplot2::element_line(size=0),
			axis.title.x=ggplot2::element_text(size=8)) +
		ggplot2::labs(y="", x=paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}



fp_singlesnpmr_SUSD2 <- mr_forest_plot(res_ad_single[res_ad_single$exposure=="SUSD2",])

pdf("MRres_single_SUSD2.pdf")
fp_singlesnpmr_SUSD2
dev.off()




FullmodelsMI_pQTLs <- data.frame()
TFulltable <- as.data.frame(Fulltable_pqtlanalysis)

Regions <- unique(TFulltable$region_code)
for(x in 408:413) {
pqtlanalyzed <- colnames(TFulltable)[x]
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mi ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmi <- as.data.frame(table(X$model$mi))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmi[2,2], Controls =countsmi[1,2])
X$Region <- Regions[z]
X$pqtl <- pqtlanalyzed
FullmodelsMI_pQTLs <- bind_rows(FullmodelsMI_pQTLs, X)
}
}

FullmodelsMI_pQTLs$factor <-rownames(FullmodelsMI_pQTLs)

pqtlresMI <- FullmodelsMI_pQTLs[FullmodelsMI_pQTLs$factor %like% "TFulltable",]


pqtllev <- colnames(pQTLSdosage)[2:7]
pqtlsmares_MI <- data.frame()




for(z in 1:length(pqtllev)){

SNPresdfu <- pqtlresMI[pqtlresMI$pqtl %like% pqtllev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$pqtl <- pqtllev[z]
pqtlsmares_MI <- bind_rows(pqtlsmares_MI, SNPindres)
}

FullmodelsIHD_pQTLs <- data.frame()
TFulltable <- as.data.frame(Fulltable_pqtlanalysis)
Regions <- unique(TFulltable$region_code)
for(x in 408:413) {
pqtlanalyzed <- colnames(TFulltable)[x]
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(ihd ~ TFulltable[,x] + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsihd <- as.data.frame(table(X$model$ihd))

X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsihd[2,2], Controls =countsihd[1,2])
X$Region <- Regions[z]
X$pqtl <- pqtlanalyzed
FullmodelsIHD_pQTLs <- bind_rows(FullmodelsIHD_pQTLs, X)


FullmodelsIHD_pQTLs$factor <-rownames(FullmodelsIHD_pQTLs)

pqtlresIHD <- FullmodelsIHD_pQTLs[FullmodelsIHD_pQTLs$factor %like% "TFulltable",]


pqtllev <- colnames(pQTLSdosage)[2:7]

pqtlsmares_IHD <- data.frame()




for(z in 1:length(pqtllev)){

SNPresdfu <- pqtlresIHD[pqtlresIHD$pqtl %like% pqtllev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$pqtl <- pqtllev[z]
pqtlsmares_IHD <- bind_rows(pqtlsmares_IHD, SNPindres)
}


pqtlsmares_MI <- left_join(pqtlsmares_MI,pQTLsmap, by = c("pqtl"="rsID"))
pqtlsmares_IHD <- left_join(pqtlsmares_IHD,pQTLsmap, by = c("pqtl"="rsID"))

pqtl_susd2_2MR <- format_data(AllPQTLs, type="exposure", header= TRUE, phenotype_col = "Assay", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EA_FREQ", chr_col = "CHR", pos_col = "BP", effect_allele_col = "EA", other_allele_col = "AA", pval_col="P")

pqtlsmares_MI <- left_join(pqtlsmares_MI, infotable[,c("rsid","A1_freq")], by = c("Markername"="rsid"))
mi_susd2_2MR <- format_data(pqtlsmares_MI, type="outcome", header= TRUE , snp_col = "Markername", beta_col = "beta", se_col = "SE", eaf_col = "A1_freq", chr_col = "Chr", pos_col = "Position", effect_allele_col = "A1", other_allele_col = "A2", pval_col="P" )

mi_susd2_2MR$outcome <- "mi"

mr_data_ad_MI <- harmonise_data( exposure = pqtl_susd2_2MR, outcome = mi_susd2_2MR)

res_ad_MI <- mr(mr_data_ad_MI)


res_ad_single_MI <- mr_singlesnp(mr_data_ad_MI)


pqtl_susd2_2MR <- format_data(AllPQTLs, type="exposure", header= TRUE, phenotype_col = "Assay", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EA_FREQ", chr_col = "CHR", pos_col = "BP", effect_allele_col = "EA", other_allele_col = "AA", pval_col="P")

pqtlsmares_IHD <- left_join(pqtlsmares_IHD, infotable[,c("rsid","A1_freq")], by = c("Markername"="rsid"))
ihd_susd2_2MR <- format_data(pqtlsmares_IHD, type="outcome", header= TRUE , snp_col = "Markername", beta_col = "beta", se_col = "SE", eaf_col = "A1_freq", chr_col = "Chr", pos_col = "Position", effect_allele_col = "A1", other_allele_col = "A2", pval_col="P" )

ihd_susd2_2MR$outcome <- "ihd"

mr_data_ad_IHD <- harmonise_data( exposure = pqtl_susd2_2MR, outcome = ihd_susd2_2MR)

res_ad_IHD <- mr(mr_data_ad_IHD)


res_ad_single_IHD <- mr_singlesnp(mr_data_ad_IHD)



fp_multisnpmr_SUSD2_MI <- forest_plot(res_ad_MI[res_ad_MI$exposure=="SUSD2",])
fp_multisnpmr_SUSD2_IHD <- forest_plot(res_ad_IHD[res_ad_IHD$exposure=="SUSD2",])

pdf("MRres_Multi_SUSD2_MI.pdf")
fp_multisnpmr_SUSD2_MI
dev.off()


pdf("MRres_Multi_SUSD2_IHD.pdf")
fp_multisnpmr_SUSD2_IHD
dev.off()


res_ad_loo_MCE <- mr_leaveoneout(mr_data_ad_MCE)

res_ad_loo_MI <- mr_leaveoneout(mr_data_ad_MI)

res_ad_loo_IHD <- mr_leaveoneout(mr_data_ad_IHD)







FullmodelsBMI_std <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)





TFulltable <- TFulltable[!is.na(TFulltable$bmi_calc),]


TFulltable$BMI_std <- (TFulltable$bmi_calc-mean(TFulltable$bmi_calc))/sd(TFulltable$bmi_calc)

Regions <- unique(TFulltable$region_code)
floatresBMI_std <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){

X <- lm(BMI_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
fl <- float(X) 
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(BMI_std))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsBMI_std <- bind_rows(FullmodelsBMI_std, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresBMI_std <- bind_rows(floatresBMI_std,fl)
}


head(FullmodelsBMI_std)
str(FullmodelsBMI_std)
str(floatresBMI_std)





FullmodelsBMI_std$factor <-rownames(FullmodelsBMI_std)


GRSQresBMI_std <- FullmodelsBMI_std[FullmodelsBMI_std$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresBMI_std <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresBMI_std[GRSQresBMI_std$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"],  )
SNPindres$level <- GRSQlev[z]
GRSQmaresBMI_std<- bind_rows(GRSQmaresBMI_std, SNPindres)
}


levs <- unique(floatresBMI_std$level)


floatresmaresBMI_std <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresBMI_std[floatresBMI_std$level==levs[z],]

SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"], means = SNPresdfu[,"mean"], cases= SNPresdfu[,"Counts"] )
SNPindres$level <- levs[z]
floatresmaresBMI_std <- bind_rows(floatresmaresBMI_std, SNPindres)
}





FullmodelsSBP_std <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)

TFulltable$sbp_std <- (TFulltable$sbp_mean-mean(TFulltable$sbp_mean))/sd(TFulltable$sbp_mean)

Regions <- unique(TFulltable$region_code)
floatresSBP_std <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(sbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
fl <- float(X) 
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(sbp_std))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsSBP_std <- bind_rows(FullmodelsSBP_std, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresSBP_std <- bind_rows(floatresSBP_std,fl)
}


head(FullmodelsSBP_std)
str(FullmodelsSBP_std)
str(floatresSBP_std)





FullmodelsSBP_std$factor <-rownames(FullmodelsSBP_std)



GRSQresSBP_std <- FullmodelsSBP_std[FullmodelsSBP_std$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresSBP_std <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresSBP_std[GRSQresSBP_std$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresSBP_std<- bind_rows(GRSQmaresSBP_std, SNPindres)
}


levs <- unique(floatresSBP_std$level)


floatresmaresSBP_std <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresSBP_std[floatresSBP_std$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]) , regions= SNPresdfu[,"region"], cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresSBP_std <- bind_rows(floatresmaresSBP_std, SNPindres)
}


FullmodelsDBP_std <- data.frame()
TFulltable <- as.data.frame(Fulltable[Fulltable$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)

TFulltable$dbp_std <- (TFulltable$dbp_mean-mean(TFulltable$dbp_mean))/sd(TFulltable$dbp_mean)

Regions <- unique(TFulltable$region_code)
floatresDBP_std <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(dbp_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(dbp_std))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
FullmodelsDBP_std <- bind_rows(FullmodelsDBP_std, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresDBP_std <- bind_rows(floatresDBP_std,fl)
}


head(FullmodelsDBP_std)
str(FullmodelsDBP_std)
str(floatresDBP_std)




FullmodelsDBP_std$factor <-rownames(FullmodelsDBP_std)


GRSQresDBP_std <- FullmodelsDBP_std[FullmodelsDBP_std$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresDBP_std <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresDBP_std[GRSQresDBP_std$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresDBP_std<- bind_rows(GRSQmaresDBP_std, SNPindres)
}


levs <- unique(floatresDBP_std$level)


floatresmaresDBP_std <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresDBP_std[floatresDBP_std$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"] , cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresDBP_std <- bind_rows(floatresmaresDBP_std, SNPindres)
}





Fullmodelsldl_std <- data.frame()
TFulltable <- as.data.frame(Fulltable_bio[Fulltable_bio$is_in_gwas_population_subset==1,])
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)

TFulltable <- TFulltable[!is.na(TFulltable$ldl_mmoll),]
TFulltable$ldl_std <- (TFulltable$ldl_mmoll-mean(TFulltable$ldl_mmoll))/sd(TFulltable$ldl_mmoll)


Regions <- unique(TFulltable$region_code)
floatresldl_std <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- lm(ldl_std ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, subset= region_code == Regions[z])}
means <- X$model %>% group_by(GRSQ) %>% summarize(mean=mean(ldl_std))
countsGRSQ <- as.data.frame(table(X$model$GRSQ))
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$Region <- Regions[z]
Fullmodelsldl_std <- bind_rows(Fullmodelsldl_std, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], mean= means$mean, Counts= countsGRSQ[,2]) 
floatresldl_std <- bind_rows(floatresldl_std,fl)
}


head(Fullmodelsldl_std)
str(Fullmodelsldl_std)
str(floatresldl_std)




Fullmodelsldl_std$factor <-rownames(Fullmodelsldl_std)

GRSQresldl_std <- Fullmodelsldl_std[Fullmodelsldl_std$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmaresldl_std <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresldl_std[GRSQresldl_std$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], regions = SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmaresldl_std <- bind_rows(GRSQmaresldl_std, SNPindres)
}


levs <- unique(floatresldl_std$level)


floatresmaresldl_std <- data.frame()




for(z in 1:length(levs)){

SNPresdfu <- floatresldl_std[floatresldl_std$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), regions= SNPresdfu[,"region"] , cases= SNPresdfu[,"Counts"], means = SNPresdfu[,"mean"])
SNPindres$level <- levs[z]
floatresmaresldl_std <- bind_rows(floatresmaresldl_std, SNPindres)
}


GRSQprot_res <- data.frame()

Prot_analysis <- as.data.frame(Prot_analysis)
for (i in 4:1466){
X <- glm( Prot_analysis[,i]~ GRSQ + age_at_study_date_x100.x+ agesquared + is_female.x+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11, data = Prot_analysis, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Prot_analysis)[i]
GRSQprot_res <- bind_rows(GRSQprot_res,X) 
}


GRSQprot_res$factor <- rownames(GRSQprot_res)



fwrite(GRSQprot_res,"GRSQprotres.txt")

meanGRS <- mean(Prot_analysis$GRS,na.rm=T)
sdGRS <- sd(Prot_analysis$GRS,na.rm=T)
Prot_analysis$GRS_std <- (Prot_analysis$GRS-meanGRS)/sdGRS


GRSstdres <- data.frame()
for (i in 4:1466){
prot <- colnames(Prot_analysis)[i]
X <- glm( Prot_analysis[,i]~ GRS_std + age_at_study_date_x100.x+ agesquared + is_female.x+ eversmk + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9+PC10+PC11, data = Prot_analysis, family = "gaussian")
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)))
X$prot <- colnames(Prot_analysis)[i]
GRSstdres <- bind_rows(GRSstdres,X) 
}

GRSstdres$factor <- rownames(GRSstdres)

fwrite(GRSstdres,"GRSstdres.txt")







FullmodelsMCE_std <- data.frame()
TFulltable <- as.data.frame(Fulltable)

Regions <- unique(TFulltable$region_code)
floatres_std <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family="binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mce ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmce <- as.data.frame(table(X$model$mce))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmce[2,2], Controls =countsmce[1,2])
X$Region <- Regions[z]
FullmodelsMCE_std <- bind_rows(FullmodelsMCE_std, X)

}


head(FullmodelsMCE_std)
str(FullmodelsMCE_std)
str(floatres_std)

GRS_std_mares <- data.frame()

FullmodelsMCE_std$factor <-rownames(FullmodelsMCE_std)



SNPresdfu <- FullmodelsMCE_std[FullmodelsMCE_std$factor %like% "GRS_std",]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )

GRS_std_mares <- SNPindres

SNPresdfu <- FullmodelsMCE_std[FullmodelsMCE_std$factor %like% "Intercept",]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )

GRS_std_mares_INT <- SNPindres



FullmodelsMCE_Quint <- data.frame()
TFulltable <- as.data.frame(Fulltable)
Regions <- unique(TFulltable$region_code)
floatres_Quint <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family="binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mce ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsmce <- as.data.frame(table(X$model$mce))

X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsmce[2,2], Controls =countsmce[1,2])
X$Region <- Regions[z]
FullmodelsMCE_Quint <- bind_rows(FullmodelsMCE_Quint, X)

}


head(FullmodelsMCE_Quint)
str(FullmodelsMCE_Quint)

GRS_QUINT_mares <- data.frame()

FullmodelsMCE_Quint$factor <-rownames(FullmodelsMCE_Quint)



SNPresdfu <- FullmodelsMCE_Quint[FullmodelsMCE_Quint$factor %like% "GRSQ",]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRS_QUINT_mares <- SNPindres

SNPresdfu <- FullmodelsMCE_Quint[FullmodelsMCE_Quint$factor %like% "Intercept",]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )

GRS_QUINT_mares_INT <- SNPindres
CVD_Addendp <- fread("DAR-2023-00342.cardiovascularendpoints_convert.csv" , header=T)


Fulltable_additional <- left_join(Fulltable, CVD_Addendp, by = c("csid"="csid"))



Fullmodelsich <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_ich <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(ich ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsich <- as.data.frame(table(X$model$ich))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(ich)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsich[2,2], Controls =countsich[1,2])
X$Region <- Regions[z]
Fullmodelsich <- bind_rows(Fullmodelsich, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$ich==1 & !is.na(countsGRSQ$ich),3], Controls=countsGRSQ[countsGRSQ$ich==0 & !is.na(countsGRSQ$ich),3]) 
floatres_ich <- bind_rows(floatres_ich,fl)
}


head(Fullmodelsich)
str(Fullmodelsich)
str(floatres_ich)
Fullmodelsich$factor <-rownames(Fullmodelsich)

GRSQresich <- Fullmodelsich[Fullmodelsich$factor %like% "GRSQ",]



GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_ich <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresich[GRSQresich$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_ich <- bind_rows(GRSQmares_ich, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_ich <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_ich[floatres_ich$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_ich <- bind_rows(floatresmares_ich, SNPindres)
}



Fullmodelsistroke <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_istroke <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(istroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsistroke <- as.data.frame(table(X$model$istroke))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(istroke)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsistroke[2,2], Controls =countsistroke[1,2])
X$Region <- Regions[z]
Fullmodelsistroke <- bind_rows(Fullmodelsistroke, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$istroke==1 & !is.na(countsGRSQ$istroke),3], Controls=countsGRSQ[countsGRSQ$istroke==0 & !is.na(countsGRSQ$istroke),3]) 
floatres_istroke <- bind_rows(floatres_istroke,fl)
}


head(Fullmodelsistroke)
str(Fullmodelsistroke)
str(floatres_istroke)

Fullmodelsistroke$factor <-rownames(Fullmodelsistroke)


GRSQresistroke <- Fullmodelsistroke[Fullmodelsistroke$factor %like% "GRSQ",]


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_istroke <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresistroke[GRSQresistroke$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_istroke <- bind_rows(GRSQmares_istroke, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_istroke <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_istroke[floatres_istroke$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_istroke <- bind_rows(floatresmares_istroke, SNPindres)


GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_istroke <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresistroke[GRSQresistroke$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_istroke <- bind_rows(GRSQmares_istroke, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_istroke <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_istroke[floatres_istroke$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_istroke <- bind_rows(floatresmares_istroke, SNPindres)


Fullmodelsanystroke <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_anystroke <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(anystroke ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsanystroke <- as.data.frame(table(X$model$anystroke))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(anystroke)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsanystroke[2,2], Controls =countsanystroke[1,2])
X$Region <- Regions[z]
Fullmodelsanystroke <- bind_rows(Fullmodelsanystroke, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$anystroke==1 & !is.na(countsGRSQ$anystroke),3], Controls=countsGRSQ[countsGRSQ$anystroke==0 & !is.na(countsGRSQ$anystroke),3]) 
floatres_anystroke <- bind_rows(floatres_anystroke,fl)
}


head(Fullmodelsanystroke)
str(Fullmodelsanystroke)
str(floatres_anystroke)




Fullmodelsanystroke$factor <-rownames(Fullmodelsanystroke)



GRSQresanystroke <- Fullmodelsanystroke[Fullmodelsanystroke$factor %like% "GRSQ",]






GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_anystroke <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresanystroke[GRSQresanystroke$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_anystroke <- bind_rows(GRSQmares_anystroke, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_anystroke <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_anystroke[floatres_anystroke$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_anystroke <- bind_rows(floatresmares_anystroke, SNPindres)


Fullmodelshf <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_hf <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(hf ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countshf<- as.data.frame(table(X$model$hf))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(hf)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countshf[2,2], Controls =countshf[1,2])
X$Region <- Regions[z]
Fullmodelshf <- bind_rows(Fullmodelshf, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$hf==1 & !is.na(countsGRSQ$hf),3], Controls=countsGRSQ[countsGRSQ$hf==0 & !is.na(countsGRSQ$hf),3]) 
floatres_hf <- bind_rows(floatres_hf,fl)
}


head(Fullmodelshf)
str(Fullmodelshf)
str(floatres_hf)


Fullmodelshf$factor <-rownames(Fullmodelshf)
GRSQreshf <- Fullmodelshf[Fullmodelshf$factor %like% "GRSQ",]

GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_hf <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQreshf[GRSQreshf$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_hf <- bind_rows(GRSQmares_hf, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_hf <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_hf[floatres_hf$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_hf <- bind_rows(floatresmares_hf, SNPindres)
}






Fullmodelspulhd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_pulhd <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(pulmonaryhd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countspulhd <- as.data.frame(table(X$model$pulmonaryhd))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(pulmonaryhd)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countspulhd[2,2], Controls =countspulhd[1,2])
X$Region <- Regions[z]
Fullmodelspulhd <- bind_rows(Fullmodelspulhd, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$pulmonaryhd==1 & !is.na(countsGRSQ$pulmonaryhd),3], Controls=countsGRSQ[countsGRSQ$pulmonaryhd==0 & !is.na(countsGRSQ$pulmonaryhd),3]) 
floatres_pulhd <- bind_rows(floatres_pulhd,fl)
}


head(Fullmodelspulhd)
str(Fullmodelspulhd)
str(floatres_pulhd)
Fullmodelspulhd$factor <-rownames(Fullmodelspulhd)


GRSQrespulhd <- Fullmodelspulhd[Fullmodelspulhd$factor %like% "GRSQ",]

GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_pulhd <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQrespulhd[GRSQrespulhd$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_pulhd <- bind_rows(GRSQmares_pulhd, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_pulhd <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_pulhd[floatres_pulhd$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_pulhd <- bind_rows(floatresmares_pulhd, SNPindres)
} 




Fullmodelsvascmort<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_vascmort <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(vascularmortality ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsvascmort <- as.data.frame(table(X$model$vascularmortality))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(vascularmortality)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsvascmort[2,2], Controls =countsvascmort[1,2])
X$Region <- Regions[z]
Fullmodelsvascmort <- bind_rows(Fullmodelsvascmort, X)
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$vascularmortality==1 & !is.na(countsGRSQ$vascularmortality),3], Controls=countsGRSQ[countsGRSQ$vascularmortality==0 & !is.na(countsGRSQ$vascularmortality),3]) 
floatres_vascmort <- bind_rows(floatres_vascmort,fl)
}


head(Fullmodelsvascmort)
str(Fullmodelsvascmort)
str(floatres_vascmort)



Fullmodelsvascmort$factor <-rownames(Fullmodelsvascmort)



GRSQresvascmort <- Fullmodelsvascmort[Fullmodelsvascmort$factor %like% "GRSQ",]

GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_vascmort <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresvascmort[GRSQresvascmort$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_vascmort <- bind_rows(GRSQmares_vascmort, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_vascmort <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_vascmort[floatres_vascmort$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_vascmort <- bind_rows(floatresmares_vascmort, SNPindres)
} 




Fullmodelshtnhd<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_htnhd <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(hypertensivehd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countshtnhd<- as.data.frame(table(X$model$hypertensivehd))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(hypertensivehd)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countshtnhd[2,2], Controls =countshtnhd[1,2])
X$Region <- Regions[z]
Fullmodelshtnhd <- bind_rows(Fullmodelshtnhd, X)

if(Regions[z]==12){ 
dftemp <- data.frame(GRSQ = "2", hypertensivehd=1 , n =0)
countsGRSQ <- bind_rows(countsGRSQ, dftemp)
}
fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$hypertensivehd==1 & !is.na(countsGRSQ$hypertensivehd),3], Controls=countsGRSQ[countsGRSQ$hypertensivehd==0 & !is.na(countsGRSQ$hypertensivehd),3]) 
floatres_htnhd <- bind_rows(floatres_htnhd,fl)
}


head(Fullmodelshtnhd)
str(Fullmodelshtnhd)
str(floatres_htnhd)



Fullmodelshtnhd$factor <-rownames(Fullmodelshtnhd)

GRSQreshtnhd <- Fullmodelshtnhd[Fullmodelshtnhd$factor %like% "GRSQ",]






GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_htnhd <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQreshtnhd[GRSQreshtnhd$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_htnhd <- bind_rows(GRSQmares_htnhd, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_htnhd <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_htnhd[floatres_htnhd$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_htnhd <- bind_rows(floatresmares_htnhd, SNPindres)
} 


 fwrite(floatresmares_anystroke, "floatresmares_anystroke.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_hf, "floatresmares_hf.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_ich, "floatresmares_ich.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_pulhd, "floatresmares_pulhd.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_vascmort, "floatresmares_vascmort.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_istroke, "floatresmares_istroke.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 fwrite(floatresmares_htnhd, "floatresmares_htnhd.txt" , sep = "\t", quote=FALSE, row.names=FALSE)
 

Fullmodelsihdmort<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
TFulltable$GRSQ <- as.factor(TFulltable$GRSQ)
Regions <- unique(TFulltable$region_code)
floatres_ihdmort <- data.frame()
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(fatalihd ~ GRSQ + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countfihd<- as.data.frame(table(X$model$fatalihd))
countsGRSQ <- X$model %>% group_by(GRSQ) %>% count(fatalihd)
fl <- float(X) 
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countfihd[2,2], Controls =countfihd[1,2])
X$Region <- Regions[z]
Fullmodelsihdmort <- bind_rows(Fullmodelsihdmort, X)

fl <- data.frame(level = names(fl$coef) , coef= fl$coef, variance = fl$var, region= Regions[z], Cases = countsGRSQ[countsGRSQ$fatalihd==1 & !is.na(countsGRSQ$fatalihd),3], Controls=countsGRSQ[countsGRSQ$fatalihd==0 & !is.na(countsGRSQ$fatalihd),3]) 
floatres_ihdmort <- bind_rows(floatres_ihdmort,fl)
}


head(Fullmodelsihdmort)
str(Fullmodelsihdmort)
str(floatres_ihdmort)




Fullmodelsihdmort$factor <-rownames(Fullmodelsihdmort)

GRSQresihdmort <- Fullmodelsihdmort[Fullmodelsihdmort$factor %like% "GRSQ",]






GRSQlev <- c("GRSQ2", "GRSQ3", "GRSQ4", "GRSQ5")

GRSQmares_ihdmort <- data.frame()




for(z in 1:length(GRSQlev)){

SNPresdfu <- GRSQresihdmort[GRSQresihdmort$factor %like% GRSQlev[z],]
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
SNPindres$level <- GRSQlev[z]
GRSQmares_ihdmort <- bind_rows(GRSQmares_ihdmort, SNPindres)
}


levs <- unique(floatres$level)


floatresmares_ihdmort <- data.frame()


for(z in 1:length(levs)){

SNPresdfu <- floatres_ihdmort[floatres_ihdmort$level==levs[z],]
SNPindres <-  metafunction(SNPresdfu[,"coef"], sqrt(SNPresdfu[,3]), SNPresdfu[,5] , SNPresdfu[,6] ,SNPresdfu[,"region"] )
SNPindres$level <- levs[z]
floatresmares_ihdmort <- bind_rows(floatresmares_ihdmort, SNPindres)
} 


fwrite(floatresmares_ihdmort, "floatresmares_ihdmort.txt" , sep = "\t", quote=FALSE, row.names=FALSE)




Fullmodelsich_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable  , family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable , family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(ich ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
countsich <- as.data.frame(table(X$model$ich))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsich[2,2], Controls =countsich[1,2])
X$Region <- Regions[z]
Fullmodelsich_GRSstd <- bind_rows(Fullmodelsich_GRSstd, X)
}


head(Fullmodelsich_GRSstd)
str(Fullmodelsich_GRSstd)




Fullmodelsich_GRSstd$factor <-rownames(Fullmodelsich_GRSstd)

GRS_stdresich <- Fullmodelsich_GRSstd[Fullmodelsich_GRSstd$factor %like% "GRS_std",]







GRs_stdmares_ich <- data.frame()





SNPresdfu <- GRS_stdresich
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRs_stdmares_ich <- bind_rows(GRs_stdmares_ich, SNPindres)






Fullmodelsistroke_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(istroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsistroke <- as.data.frame(table(X$model$istroke))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsistroke[2,2], Controls =countsistroke[1,2])
X$Region <- Regions[z]
Fullmodelsistroke_GRSstd <- bind_rows(Fullmodelsistroke_GRSstd, X)
}


head(Fullmodelsistroke_GRSstd)
str(Fullmodelsistroke_GRSstd)




Fullmodelsistroke_GRSstd$factor <-rownames(Fullmodelsistroke_GRSstd)

GRSQresistroke_GRSstd <- Fullmodelsistroke_GRSstd[Fullmodelsistroke_GRSstd$factor %like% "GRS_std",]







GRSQmares_istroke_GRSstd <- data.frame()





SNPresdfu <- GRSQresistroke_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_istroke_GRSstd <- SNPindres






Fullmodelsanystroke_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable,  family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(anystroke ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
countsanystroke <- as.data.frame(table(X$model$anystroke))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsanystroke[2,2], Controls =countsanystroke[1,2])
X$Region <- Regions[z]
Fullmodelsanystroke_GRSstd <- bind_rows(Fullmodelsanystroke_GRSstd, X)
}


head(Fullmodelsanystroke_GRSstd)
str(Fullmodelsanystroke_GRSstd)




Fullmodelsanystroke_GRSstd$factor <-rownames(Fullmodelsanystroke_GRSstd)

GRSQresanystroke_GRSstd <- Fullmodelsanystroke_GRSstd[Fullmodelsanystroke_GRSstd$factor %like% "GRS_std",]







GRSQmares_anystroke_GRSstd <- data.frame()





SNPresdfu <- GRSQresanystroke_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_anystroke_GRSstd <- SNPindres







Fullmodelshf_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable,  family = "binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial",subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(hf ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countshf<- as.data.frame(table(X$model$hf))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countshf[2,2], Controls =countshf[1,2])
X$Region <- Regions[z]
Fullmodelshf_GRSstd <- bind_rows(Fullmodelshf_GRSstd, X)
}


head(Fullmodelshf_GRSstd)
str(Fullmodelshf_GRSstd)




Fullmodelshf_GRSstd$factor <-rownames(Fullmodelshf_GRSstd)

GRSQreshf_GRSstd <- Fullmodelshf_GRSstd[Fullmodelshf_GRSstd$factor %like% "GRS_std",]







GRSQmares_hf_GRSstd <- data.frame()





SNPresdfu <- GRSQreshf_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_hf_GRSstd <- SNPindres







Fullmodelspulhd_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(pulmonaryhd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countspulhd <- as.data.frame(table(X$model$pulmonaryhd))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countspulhd[2,2], Controls =countspulhd[1,2])
X$Region <- Regions[z]
Fullmodelspulhd_GRSstd <- bind_rows(Fullmodelspulhd_GRSstd, X)
}


head(Fullmodelspulhd_GRSstd)
str(Fullmodelspulhd_GRSstd)




Fullmodelspulhd_GRSstd$factor <-rownames(Fullmodelspulhd_GRSstd)

GRSQrespulhd_std <- Fullmodelspulhd_GRSstd[Fullmodelspulhd_GRSstd$factor %like% "GRS_std",]







GRSQmares_pulhd_GRSstd <- data.frame()





SNPresdfu <- GRSQrespulhd_std
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_pulhd_GRSstd <- SNPindres







Fullmodelsvascmort_GRSstd<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable,  family = "binomial",subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(vascularmortality ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countsvascmort <- as.data.frame(table(X$model$vascularmortality))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countsvascmort[2,2], Controls =countsvascmort[1,2])
X$Region <- Regions[z]
Fullmodelsvascmort_GRSstd <- bind_rows(Fullmodelsvascmort_GRSstd, X)
}


head(Fullmodelsvascmort_GRSstd)
str(Fullmodelsvascmort_GRSstd)




Fullmodelsvascmort_GRSstd$factor <-rownames(Fullmodelsvascmort_GRSstd)

GRSQresvascmort_GRSstd <- Fullmodelsvascmort_GRSstd[Fullmodelsvascmort_GRSstd$factor %like% "GRS_std",]







GRSQmares_vascmort_GRSstd <- data.frame()





SNPresdfu <- GRSQresvascmort_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_vascmort_GRSstd <- SNPindres







Fullmodelshtnhd_GRSstd <- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(hypertensivehd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
countshtnhd<- as.data.frame(table(X$model$hypertensivehd))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countshtnhd[2,2], Controls =countshtnhd[1,2])
X$Region <- Regions[z]
Fullmodelshtnhd_GRSstd <- bind_rows(Fullmodelshtnhd_GRSstd, X)

}


head(Fullmodelshtnhd_GRSstd)
str(Fullmodelshtnhd_GRSstd)




Fullmodelshtnhd_GRSstd$factor <-rownames(Fullmodelshtnhd_GRSstd)

GRSQreshtnhd_GRSstd <- Fullmodelshtnhd_GRSstd[Fullmodelshtnhd_GRSstd$factor %like% "GRS_std",]







GRSQmares_htnhd_GRSstd <- data.frame()





SNPresdfu <- GRSQreshtnhd_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_htnhd_GRSstd <- SNPindres


 
Fullmodelsihdmort_GRSstd<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(fatalihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
countfihd<- as.data.frame(table(X$model$fatalihd))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countfihd[2,2], Controls =countfihd[1,2])
X$Region <- Regions[z]
Fullmodelsihdmort_GRSstd <- bind_rows(Fullmodelsihdmort_GRSstd, X)

}


head(Fullmodelsihdmort_GRSstd)
str(Fullmodelsihdmort_GRSstd)




Fullmodelsihdmort_GRSstd$factor <-rownames(Fullmodelsihdmort_GRSstd)

GRSQresihdmort_GRSstd <- Fullmodelsihdmort_GRSstd[Fullmodelsihdmort_GRSstd$factor %like% "GRS_std",]







GRSQmares_ihdmort_GRSstd <- data.frame()





SNPresdfu <- GRSQresihdmort_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_ihdmort_GRSstd <- SNPindres
 

Fullmodelsmi_GRSstd<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(mi ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countmi<- as.data.frame(table(X$model$mi))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countmi[2,2], Controls =countmi[1,2])
X$Region <- Regions[z]
Fullmodelsmi_GRSstd <- bind_rows(Fullmodelsmi_GRSstd, X)

}


head(Fullmodelsmi_GRSstd)
str(Fullmodelsmi_GRSstd)




Fullmodelsmi_GRSstd$factor <-rownames(Fullmodelsmi_GRSstd)

GRSQresmi_GRSstd <- Fullmodelsmi_GRSstd[Fullmodelsmi_GRSstd$factor %like% "GRS_std",]







GRSQmares_mi_GRSstd <- data.frame()





SNPresdfu <- GRSQresmi_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_mi_GRSstd <- SNPindres

 

Fullmodelsihd_GRSstd<- data.frame()
TFulltable <- as.data.frame(Fulltable_additional)
Regions <- unique(TFulltable$region_code)
for(z in 1:length(Regions)){

if(Regions[z]==36 | Regions[z]==88){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8 + RC_pc9 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])

}

else if (Regions[z]==68 | Regions[z]==58){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6 + RC_pc7 + RC_pc8, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==26 | Regions[z]==78){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5, data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==52){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 + RC_pc5 + RC_pc6, data = TFulltable,  family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==46){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 + RC_pc4 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
else if (Regions[z]==12){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 + RC_pc3 , data = TFulltable, family = "binomial",  subset= region_code == Regions[z])}
else if (Regions[z]==16){
X <- glm(ihd ~ GRS_std + age_at_study_date_x100+ agesquared + is_female+ eversmk + RC_pc1 + RC_pc2 , data = TFulltable, family = "binomial", subset= region_code == Regions[z])}
countihd<- as.data.frame(table(X$model$ihd))
X <- cbind(as.data.frame(summary(X)$coefficients) , as.data.frame(confint.default(X)) , Cases = countihd[2,2], Controls =countihd[1,2])
X$Region <- Regions[z]
Fullmodelsihd_GRSstd <- bind_rows(Fullmodelsihd_GRSstd, X)

}


head(Fullmodelsihd_GRSstd)
str(Fullmodelsihd_GRSstd)




Fullmodelsihd_GRSstd$factor <-rownames(Fullmodelsihd_GRSstd)

GRSQresihd_GRSstd <- Fullmodelsihd_GRSstd[Fullmodelsihd_GRSstd$factor %like% "GRS_std",]







GRSQmares_ihd_GRSstd <- data.frame()





SNPresdfu <- GRSQresihd_GRSstd
SNPindres <-  metafunction(SNPresdfu[,1], SNPresdfu[,2], SNPresdfu[,"Cases"] , SNPresdfu[,"Controls"] , SNPresdfu[,"Region"] )
GRSQmares_ihd_GRSstd <- SNPindres
 
 
 
 
 
 
 
fwrite(GRs_stdmares_ich,"GRs_stdmares_ich.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_istroke_GRSstd,"GRSQmares_istroke_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_anystroke_GRSstd,"GRSQmares_anystroke_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_vascmort_GRSstd,"GRSQmares_vascmort_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_hf_GRSstd,"GRSQmares_hf_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_htnhd_GRSstd,"GRSQmares_htnhd_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_pulhd_GRSstd,"GRSQmares_pulhd_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_ihdmort_GRSstd,"GRSQmares_ihdmort_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRS_std_mares, "GRS_std_mares_MCE.txt" , quote=FALSE, row.names=FALSE, sep= "\t")
fwrite(GRSQmares_mi_GRSstd, "GRSQmares_mi_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t") 
fwrite(GRSQmares_ihd_GRSstd, "GRSQmares_ihd_GRSstd.txt" , quote=FALSE, row.names=FALSE, sep= "\t") 
 
