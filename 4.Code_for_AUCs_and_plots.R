######
###### Code for evaluating the AUC and plots
#-----------------------------------------------------------------------------------------------------------------------
setwd('/PATH')

## Input data file
Data_to_calculate=read.csv('./Merge_AgeSex_adjusted_protein_levels.csv',header = T)
Model_estimate=read.csv('./Summary_model_estimate.csv',header = T)

#--------------------------
## Score calculation for each individual
#--------------------------
ID_score=data.frame(matrix(ncol=2,nrow = length(Data_to_calculate[,1])))
colnames(ID_score)=c("ID","AD_risk_score")
ID_score[,1]=Data_to_calculate[,1]


for (k in 1:length(ID_score[,1])){
  Individual_level=data.frame(t(Data_to_calculate[k,5:25]))
  colnames(Individual_level)=c("Protein_levels")

  if (sum(is.na(Individual_level))>0){
    ID_score[k,2]<-NA
  } else {
    ID_score[k,2]=(1/(1+exp(0-sum(Individual_level$Protein_levels[1:Number_proteins]*Model_estimate$Estimate[2:22])-Model_estimate$Estimate[1])))*100
  }
}

ID_score=cbind(Data_to_calculate,ID_score)
write.csv(ID_score,'./Merge_AD_risk_scores.csv',row.names = F)

#--------------------------
# AUC evaluation
#--------------------------
ID_score_CN<-subset(ID_score,ID_score$Diagnosis=="CN")
ID_score_MCI<-subset(ID_score,ID_score$Diagnosis=="MCI")
ID_score_AD<-subset(ID_score,ID_score$Diagnosis=="AD")

AUC_test_dataset=rbind(ID_score_AD,ID_score_CN)
AUC_test_dataset$Diagnosis<-relevel(factor(AUC_test_dataset$Diagnosis),ref = "CN")
AUC_test<-plot.roc(AUC_test_dataset$Diagnosis,AUC_test_dataset$AD_risk_score,percent=TRUE,col="#008600")
plot(AUC_test)
auc(AUC_test)

#--------------------------
## Calculation of scores and Radar plot for five biological processes
#--------------------------
setwd('/PATH')
# Input the reference file of the 5 proteins containing baseline levels (score = 100), bottom line levels (score = 0), and cut-off level (score = 60) for normal/abnormal range
Ref_5protein=read.csv('./Summary_5proteins_cutoff.csv',header = T)

Data_to_plot=ID_score
## Only include the adjusted levels of 5 proteins (i.e., NEFL,PARP1,CD33,LIFR,PPY)
Protein_position=na.omit(match(names(Ref_5protein)[2:6],names(Data_to_plot)))
Data_to_plot=Data_to_plot[,c(1,Protein_position)]

#-----
##### Calculate normalized protein scores
#-----
# Generate beta and intercept for normalization into a scale of 0 to 100
Ref_linear_asso_file=data.frame(matrix(ncol=10,nrow=2))
rownames(Ref_linear_asso_file)=c("Abnormal","Normal")
colnames(Ref_linear_asso_file)=c("Beta","Intercept","Beta","Intercept","Beta","Intercept","Beta","Intercept","Beta","Intercept")

for (q in 1:5){
  Dataset_Abnormal=Ref_5protein[1:2,c(1,q+1)]
  Dataset_Normal=Ref_5protein[2:3,c(1,q+1)]
  colnames(Dataset_Abnormal)=c("Scale","Protein")
  colnames(Dataset_Normal)=c("Scale","Protein")
  # Beta and intercept for the abnormal range
  Linear_asso_test_Abnormal<-lm(Scale~Protein,data = Dataset_Abnormal)
  Ref_linear_asso_file[1,(2*q-1)]=as.vector(t(summary(Linear_asso_test_Abnormal)$coefficients[2,1]))
  Ref_linear_asso_file[1,(2*q)]=as.vector(t(summary(Linear_asso_test_Abnormal)$coefficients[1,1]))
  # Beta and intercept for the Normal range
  Linear_asso_test_Normal<-lm(Scale~Protein,data = Dataset_Normal)
  Ref_linear_asso_file[2,(2*q-1)]=as.vector(t(summary(Linear_asso_test_Normal)$coefficients[2,1]))
  Ref_linear_asso_file[2,(2*q)]=as.vector(t(summary(Linear_asso_test_Normal)$coefficients[1,1]))
}

# Normalization of individual's 5 protein scores
Data_to_plot_normalized=Data_to_plot

for (x in 1:length(Data_to_plot_normalized[,1])){
  Individual_5protein_data=Data_to_plot_normalized[x,]
  for (y in 1:5){
    Individual_protein_data=Individual_5protein_data[1,(y+1)]
    # If the protein is decreased in abnormal situation
    if (Ref_linear_asso_file[1,(2*y-1)]>0){
      # If the protein level is NA
      if (is.na(Individual_protein_data)){
        Data_to_plot_normalized[x,(y+1)]<-NA
      }
      # If the protein level is below the lower limit
      else if (Individual_protein_data<=Ref_5protein[1,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Ref_5protein[1,1]
      }
      # If the protein level is above the upper limit
      else if (Individual_protein_data>=Ref_5protein[3,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Ref_5protein[3,1]
      }
      # If the protein level is within the abnormal range
      else if (Individual_protein_data<=Ref_5protein[2,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Data_to_plot_normalized[x,(y+1)]*Ref_linear_asso_file[1,(2*y-1)]+Ref_linear_asso_file[1,(2*y)]
      }
      # If the protein level is within the normal range
      else if (Individual_protein_data>Ref_5protein[2,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Data_to_plot_normalized[x,(y+1)]*Ref_linear_asso_file[2,(2*y-1)]+Ref_linear_asso_file[2,(2*y)]
      }
    }
    
    # If the protein is increased in abnormal situation
    if (Ref_linear_asso_file[1,(2*y-1)]<0){
      # If the protein level is NA
      if (is.na(Individual_protein_data)){
        Data_to_plot_normalized[x,(y+1)]<-NA
      }
      # If the protein level is below the lower limit
      else if (Individual_protein_data<=Ref_5protein[3,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Ref_5protein[3,1]
      }
      # If the protein level is above the upper limit
      else if (Individual_protein_data>=Ref_5protein[1,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Ref_5protein[1,1]
      }
      # If the protein level is within the abnormal range
      else if (Individual_protein_data>=Ref_5protein[2,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Data_to_plot_normalized[x,(y+1)]*Ref_linear_asso_file[1,(2*y-1)]+Ref_linear_asso_file[1,(2*y)]
      }
      # If the protein level is within the normal range
      else if (Individual_protein_data<Ref_5protein[2,(y+1)]){
        Data_to_plot_normalized[x,(y+1)]=Data_to_plot_normalized[x,(y+1)]*Ref_linear_asso_file[2,(2*y-1)]+Ref_linear_asso_file[2,(2*y)]
      }
    }
  }
}

write.csv(Data_to_plot_normalized,"./Merge_5biologicalprocesses_score.csv",row.names = F)

#--------------------------
##### Radar plot of mean levels of 5 protein scores of CN, MCI, AD
#--------------------------
setwd('/PATH')
Data_to_plot=read.csv("./Average_5biologicalprocesses_score.csv",header = T)
row.names(Data_to_plot)=Data_to_plot[,1]

library(fmsb)
# Define the variable ranges: maximum and minimum
max_min <- data.frame(
  Neurodegeneration = c(100, 0), Metabolic.activities = c(100, 0), Vascular.functions = c(100, 0),
  Innate.immunity = c(100, 0), Inflammation = c(100, 0)
)
rownames(max_min) <- c("Max", "Min")
# Bind the variable ranges to the data
df <- rbind(max_min, Data_to_plot)
df
### Define radar plot function
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 3, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 1.5,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
# Define colors and titles
colors <- c("#4B96E6", "#FD8008", "#FC6666")
titles <- c("CN", "MCI", "AD")

# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1,3))

# Create the radar chart
for(i in 1:3){
  create_beautiful_radarchart(
    data = df[c(1, 2, i+2), ], caxislabels = c(0, 25, 50, 75, 100),
    color = colors[i], title = titles[i]
  )
}
par(op)

