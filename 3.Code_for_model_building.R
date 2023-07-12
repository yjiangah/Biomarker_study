######
###### Code for model building using age and sex adjusted protein level data
#-----------------------------------------------------------------------------------------------------------------------
##### Input overall file
setwd('/PATH')

Merge_data = read.csv('./Merge_AgeSex_adjusted_protein_levels.csv', header = T)
Merge_data_CN<-subset(Merge_data,Merge_data$Diagnosis=="CN")
Merge_data_AD<-subset(Merge_data,Merge_data$Diagnosis=="AD")
Merge_data_CNAD_for_model=rbind(Merge_data_CN,Merge_data_AD)
Merge_data_CNAD_for_model$Diagnosis<-relevel(factor(Merge_data_CNAD_for_model$Diagnosis),ref = "CN")

#----------------------------------------------------------------------
##### Model building
#----------------------------------------------------------------------
library(pROC)

Data_for_model=Merge_data_CNAD_for_model[c(4,5:25)]
Data_for_model=na.omit(Data_for_model)

linmod<-glm(Diagnosis~.,data= Data_for_model,family = binomial)
pred=predict(linmod, type=c("response"))
Data_for_model$AD_risk_score=100*pred
graph<-plot.roc(Data_for_model$Diagnosis,Data_for_model$AD_risk_score,percent=TRUE,col="#008600")
plot(graph)
auc(graph)

Model_estimate=data.frame(matrix(ncol = 6,nrow = 22))
colnames(Model_estimate)=c("Protein_ID","Estimate","SE","Z_value","P_value")
Model_estimate[1,1]="Intercept"
Model_estimate[1,2:5]=as.vector(summary(linmod)$coefficients[1,1:4])
Model_estimate[2:22,1]=names(Data_for_model)[2:22]

for (q in (2:22)){
  Model_estimate[q,2:5]=as.vector(summary(linmod)$coefficients[q,1:4])
}

write.csv(Model_estimate,"./Summary_model_estimate.csv",row.names = F)
