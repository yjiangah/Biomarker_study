###### Code for generating the age and sex adjusted levels of 21 proteins
#-----------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(plotly)

#---------------------------------------------------------------------------------------------
##### Input overall file
#---------------------------------------------------------------------------------------------
setwd('/PATH')

Merge_data = read.csv('./Merge_raw_protein_levels.csv', header = T)
Merge_Pheno_data=Merge_data[c(1:4)]
Merge_data_CN<-subset(Merge_data,Merge_data$Diagnosis=="CN")

### Examine the effects of age and gender on each protein level and generate the adjusted protein levels
library(robustbase)

Adjusted_Merge_data=data.frame(matrix(ncol = 21,nrow = length(Merge_data[,1])))
colnames(Adjusted_Merge_data)=names(Merge_data)[5:25]
Adjusted_Merge_data=cbind(Merge_Pheno_data,Adjusted_Merge_data)

Age_Sex_Effects_on_proteins=data.frame(matrix(ncol = 8,nrow = 21))
colnames(Age_Sex_Effects_on_proteins)=c("Age_Estimate","Age_SE","Age_T_value","Age_P_value","Gender_Estimate","Gender_SE","Gender_T_value","Gender_P_value")
rownames(Age_Sex_Effects_on_proteins)=names(Merge_data)[5:25]

for (k in 1:21){
  Protein_age_sex_test<-lmrob(Merge_data_CN[,k+4]~Age+Gender,data=Merge_data_CN,k.max=900000)
  
  Age_Sex_Effects_on_proteins[k,1:4]=as.vector(summary(Protein_age_sex_test)$coefficients[2,])
  Age_Sex_Effects_on_proteins[k,5:8]=as.vector(summary(Protein_age_sex_test)$coefficients[3,])
  
  Adjusted_Merge_data[,k+4]=Merge_data[,k+4]-Merge_data[,2]*as.numeric(Age_Sex_Effects_on_proteins[k,1])-Merge_data[,3]*as.numeric(Age_Sex_Effects_on_proteins[k,5])
}

write.csv(Age_Sex_Effects_on_proteins,'./Summary_AgeSex_Effects.csv',row.names = T)

rownames(Adjusted_Merge_data)=Merge_data$ID
write.csv(Adjusted_Merge_data,'./Merge_AgeSex_adjusted_protein_levels.csv',row.names = T)
