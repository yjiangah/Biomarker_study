######
###### Code for association tests between candidate protein levels and AD/MCI/related endophenotypes
#-----------------------------------------------------------------------------------------------------------------------
library(som)
library(robustbase)

##### Input overall file
setwd('/PATH')
Merge_data = read.csv('./Merged_data_for_asso_test.csv', header = T)

## Set NC as reference
Merge_data$Diagnosis<-relevel(factor(Merge_data$Diagnosis),ref = "NC")

## Normalization of protein levels
Merge_Pheno_data=Merge_data[c(1:23)]
Merge_Protein_data=Merge_data[c(24:44)]
Merge_normalized_Protein_data=normalize(Merge_Protein_data,byrow = F)
Merge_data=cbind(Merge_Pheno_data,Merge_normalized_Protein_data)

#-----------------------------------------------------------------------------------------------------------------------
### Test of protein levels ~ AD and MCI
Data_for_test=Merge_data

Phenotype_test_summary=data.frame(matrix(ncol=16,nrow=length(Data_for_test[1,c(24:44)])))
colnames(Phenotype_test_summary)<-c("AD_Estimate","AD_SE","AD_T_value","AD_P_value","MCI_Estimate","MCI_SE","MCI_T_value","MCI_P_value","Age_Estimate","Age_SE","Age_T_value","Age_P_value","Gender_Estimate","Gender_SE","Gender_T_value","Gender_P_value")
rownames(Phenotype_test_summary)=names(Data_for_test[c(24:44)])

for (i in 1:length(Data_for_test[1,c(24:44)]))
{
  tryCatch({
    Phenotype_test<-lmrob(as.vector(Data_for_test[,i+23])~Diagnosis+Age+Gender+BMI+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia,data=Data_for_test,k.max=900000)
    
    Phenotype_test_summary[i,1:4]=as.vector(summary(Phenotype_test)$coefficients[2,1:4])
    Phenotype_test_summary[i,5:8]=as.vector(summary(Phenotype_test)$coefficients[3,1:4])
    Phenotype_test_summary[i,9:12]=as.vector(summary(Phenotype_test)$coefficients[4,1:4])
    Phenotype_test_summary[i,13:16]=as.vector(summary(Phenotype_test)$coefficients[5,1:4])
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

Phenotype_test_summary$AD_FDR=p.adjust(Phenotype_test_summary$AD_P_value,method='fdr')
Phenotype_test_summary$MCI_FDR=p.adjust(Phenotype_test_summary$MCI_P_value,method='fdr')
Phenotype_test_summary$Age_FDR=p.adjust(Phenotype_test_summary$Age_P_value,method='fdr')
Phenotype_test_summary$Gender_FDR=p.adjust(Phenotype_test_summary$Gender_P_value,method='fdr')

write.csv(Phenotype_test_summary, './Summary_AD_MCI_Age_Gender_test.csv',row.names = T)

#-----------------------------------------------------------------------------------------------------------------------
### Test of protein levels ~ MoCA / Hippocampal / Gray matter / pTau181 / NfL / Ab42.40
setwd('/PATH')
Merge_data = read.csv('./Merged_endophenotype_data_for_asso_test.csv', header = T)
## Normalization of protein levels
Data_pheno=Merge_data[c(1,12:20)]
Data_protein=Merge_data[c(2:3,5:7,10:11,21,24:44)]
Data_normalized_protein=normalize(Data_protein,byrow = F)
Data_for_test=cbind(Data_pheno,Data_normalized_protein)

## Association tests between each endophenotypes and each protein levels 
Phenotype_test_summary=data.frame(matrix(ncol=24,nrow=length(Data_for_test[1,c(19:39)])))
colnames(Phenotype_test_summary)<-c("Gray_matter_Estimate","Gray_matter_SE","Gray_matter_T_value","Gray_matter_P_value","Hippocampu_Estimate","Hippocampu_SE","Hippocampu_T_value","Hippocampu_P_value","Ptau181_Estimate","Ptau181_SE","Ptau181_T_value","Ptau181_P_value","Nfl_Estimate","Nfl_SE","Nfl_T_value","Nfl_P_value","Ab42.40_Estimate","Ab42.40_SE","Ab42.40_T_value","Ab42.40_P_value","MoCA_Estimate","MoCA_SE","MoCA_T_value","MoCA_P_value")
rownames(Phenotype_test_summary)=names(Data_for_test[c(19:39)])

for (k in 1:6){
  Data_to_test=Data_for_test[c(k+10,2:10,19:39)]
  names(Data_to_test)[1]="Endophenotype"
  for (i in 1:21)
  {
    tryCatch({
      Phenotype_test<-lm(as.vector(Data_to_test[,i+10])~Endophenotype+Age+Gender+BMI+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia,data=Data_to_test)
      
      Phenotype_test_summary[i,(4*k-3):(4*k)]=as.vector(summary(Phenotype_test)$coefficients[2,1:4])
      
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
  }
}

Phenotype_test_summary$Gray_matter_FDR=p.adjust(Phenotype_test_summary$Gray_matter_P_value,method='fdr')
Phenotype_test_summary$Hippocampu_FDR=p.adjust(Phenotype_test_summary$Hippocampu_P_value,method='fdr')
Phenotype_test_summary$Ptau181_FDR=p.adjust(Phenotype_test_summary$Ptau181_P_value,method='fdr')
Phenotype_test_summary$Nfl_FDR=p.adjust(Phenotype_test_summary$Nfl_P_value,method='fdr')
Phenotype_test_summary$Ab42.40_FDR=p.adjust(Phenotype_test_summary$Ab42.40_P_value,method='fdr')
Phenotype_test_summary$MoCA_FDR=p.adjust(Phenotype_test_summary$MoCA_P_value,method='fdr')

write.csv(Phenotype_test_summary, './Summary_Endophenotype_test.csv',row.names = T)

