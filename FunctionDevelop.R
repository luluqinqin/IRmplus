install.packages("MplusAutomation")
install.packages("tidyverse")
library("MplusAutomation")
library("stringi")
library("tidyverse")

install.packages("devtools")
install.packages("roxygen2")
library("devtools")
library("roxygen2")


####try another example
example = read.csv("J:\\PSYC216\\Lecture6_SEM\\Example2.csv")
example = read.csv("C:\\Users\Lu\\Dropbox\\RFunctionDevelopment\\Example2.csv")
colnames(example) = NULL
write.table(example, "J:\\Research Studies\\RFunctionDevelopment\\try.example.dat")
#######################
#EXTRACT COVARIANCE PARAGRAPH!
txt = readLines("J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_nointeraction.out")
df <- data.frame(txt, paragraph = cumsum(txt == "")) 
# find  paragraphs with the search term
keep_paragraph <- df[grep("Covariances", df[, "txt"]), "paragraph"] 
# subset the data.frame
df <- df[df$paragraph %in% keep_paragraph,] 
df = df[grepl("GENDER",df[,1]),]
var2 = as.numeric(str_sub(last(df[,1]), -10))
############################################
readModels("J:\\Research Studies\\RFunctionDevelopment\\Model1_nointeraction_continuous.out", what="summaries")
readModels("J:\\Research Studies\\RFunctionDevelopment\\Model1_nointeraction_continuous.out", what="mode_indices")
readModels("J:\\Research Studies\\RFunctionDevelopment\\Model1_nointeraction_continuous.out", what="sampstat")
##for observed IV interact with the latent
#readModels("J:\\Research Studies\\RFunctionDevelopment\\Model1_nointeraction_continuous.out", what="sampstat")

############ TRY IT OUT!!!!
#### Continuous
DV = c("NER")
IV1 = c("DEGB")
IV2 = c("ICD")
INTER= c("INT")
model0_out = c("J:\\Research Studies\\RFunctionDevelopment\\model_lv_nointeraction.out")
model1_out = c("J:\\Research Studies\\RFunctionDevelopment\\model_lv_interaction.out")

m0_param = as.data.frame(readModels(model0_out , what="parameters")$parameters$unstandardized)
m0_betas = m0_param[grepl(paste0(DV,c("."), c("ON")), m0_param[, grepl("paramHeader", names(m0_param))]),]
m0_beta1 = m0_betas[grepl(paste0(IV1), m0_betas$param),]$est
m0_beta2 = m0_betas[grepl(paste0(IV2), m0_betas$param),]$est

m0_var = m0_param[grepl(paste0("Variances"), m0_param[, grepl("paramHeader", names(m0_param))]),]
m0_var1 = m0_var[grepl(paste0(c("^"),IV1,c("$")), m0_var$param),]$est
m0_var2 = m0_var[grepl(paste0(c("^"),IV2,c("$")), m0_var$param),]$est
m0_res = m0_var[grepl(paste0(c("^"),DV,c("$")), m0_var$param),]$est
m0_Rsq= ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2) / ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2*m0_res)

m1_param = as.data.frame(readModels(model1_out, what="parameters")$parameters$unstandardized)
m1_betas = m1_param[grepl(paste0(DV,c("."), c("ON")), m1_param[, grepl("paramHeader", names(m1_param))]),]
m1_beta1 = m1_betas[grepl(paste0(IV1), m1_betas$param),]$est
m1_beta2 = m1_betas[grepl(paste0(IV2), m1_betas$param),]$est
m1_beta3 = m1_betas[grepl(paste0(INTER), m1_betas$param),]$est

m1_var = m1_param[grepl(paste0("Variances"), m1_param[, grepl("paramHeader", names(m1_param))]),]
m1_var1 = m1_var[grepl(paste0(c("^"),IV1,c("$")), m1_var$param),]$est
m1_var2 = m1_var[grepl(paste0(c("^"),IV2,c("$")), m1_var$param),]$est
m1_res = m1_var[grepl(paste0(c("^"),DV,c("$")), m1_var$param),]$est
m1_cov = m1_param[grepl(paste0(c("WITH")), m1_param[, grepl("paramHeader", names(m1_param))]),]$est
m1_Rsq = ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2)) / ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2) + m1_res)

Interaction_Rsq = m1_Rsq - m0_Rsq
#######################


##############################################################################################
############################## STEP ONE - TWO WAY INTERACTION ################################
##############################################################################################

############################ 
#Latent Variable Interaction
############################
endogenous = c("NER")
exogenous = c("DEGB")
moderator = c("ICD")
interaction = c("INT")
M0 = c("J:\\Research Studies\\RFunctionDevelopment\\model_lv_nointeraction.out")
M1 = c("J:\\Research Studies\\RFunctionDevelopment\\model_lv_interaction.out")

LIR = function(M0, M1, endogenous, exogenous, moderator, interaction){
  m0_param = as.data.frame(readModels(M0, what="parameters")$parameters$unstandardized)
  m0_betas = m0_param[grepl(paste0(endogenous,c("."), c("ON")), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_beta1 = m0_betas[grepl(paste0(exogenous), m0_betas$param),]$est
  m0_beta2 = m0_betas[grepl(paste0(moderator), m0_betas$param),]$est
  
  m0_var = m0_param[grepl(paste0("Variances"), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_var1 = m0_var[grepl(paste0(c("^"),exogenous,c("$")), m0_var$param),]$est
  m0_var2 = m0_var[grepl(paste0(c("^"),moderator,c("$")), m0_var$param),]$est
  m0_res = m0_var[grepl(paste0(c("^"),endogenous,c("$")), m0_var$param),]$est
  m0_Rsq= ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2) / ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2 + m0_res)
  
  m1_param = as.data.frame(readModels(M1, what="parameters")$parameters$unstandardized)
  m1_betas = m1_param[grepl(paste0(endogenous,c("."), c("ON")), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_beta1 = m1_betas[grepl(paste0(exogenous), m1_betas$param),]$est
  m1_beta2 = m1_betas[grepl(paste0(moderator), m1_betas$param),]$est
  m1_beta3 = m1_betas[grepl(paste0(interaction), m1_betas$param),]$est
  
  m1_var = m1_param[grepl(paste0("Variances"), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_var1 = m1_var[grepl(paste0(c("^"),exogenous,c("$")), m1_var$param),]$est
  m1_var2 = m1_var[grepl(paste0(c("^"),moderator,c("$")), m1_var$param),]$est
  m1_res = m1_var[grepl(paste0(c("^"),endogenous,c("$")), m1_var$param),]$est
  m1_cov = m1_param[grepl(paste0(c("WITH")), m1_param[, grepl("paramHeader", names(m1_param))]),]$est
  m1_Rsq = ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2)) / ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2) + m1_res)
  
  Interaction_Rsq = m1_Rsq - m0_Rsq
  return(Interaction_Rsq)
  
}

LIR(M0 = "J:\\Research Studies\\RFunctionDevelopment\\model_lv_nointeraction.out", 
    M1 = "J:\\Research Studies\\RFunctionDevelopment\\model_lv_interaction.out", 
    endogenous="NER", exogenous = "DEGB", moderator = "ICD", interaction= "INT")  



#############################
#Latent & Observed Interaction
##############################
endogenous = c("NER")
exogenous = c("DEGB")
moderator = c("GENDER")
interaction = c("INT")
M0 = c("J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_nointeraction.out")
M1 = c("J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_interaction.out")

LOIR = function(M0, M1, endogenous, exogenous, moderator, interaction){
#M0 model
  m0_param = as.data.frame(readModels(M0, what="parameters")$parameters$unstandardized)
  m0_betas = m0_param[grepl(paste0(endogenous,c("."), c("ON")), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_beta1 = m0_betas[grepl(paste0(exogenous), m0_betas$param),]$est
  m0_beta2 = m0_betas[grepl(paste0(moderator), m0_betas$param),]$est
  m0_var = m0_param[grepl(paste0("Variances"), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_var1 = m0_var[grepl(paste0(c("^"),exogenous,c("$")), m0_var$param),]$est
  m0_res = m0_var[grepl(paste0(c("^"),endogenous,c("$")), m0_var$param),]$est
  #Get Observed Variable Variances
  txt1 = readLines(M0)
  df1 <- data.frame(txt1, paragraph = cumsum(txt1 == "")) 
  keep_paragraph1 <- df1[grep("Covariances", df1[, "txt1"]), "paragraph"] 
  df1 <- df1[df1$paragraph %in% keep_paragraph1,] 
  df1 = df1[grepl(moderator,df1[,1]),]
  m0_var2 = as.numeric(str_sub(last(df1[,1]), -10))
  m0_Rsq= ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2) / ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2 + m0_res)

#M1 Model
  m1_param = as.data.frame(readModels(M1, what="parameters")$parameters$unstandardized)
  m1_betas = m1_param[grepl(paste0(endogenous,c("."), c("ON")), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_beta1 = m1_betas[grepl(paste0(exogenous), m1_betas$param),]$est
  m1_beta2 = m1_betas[grepl(paste0(moderator), m1_betas$param),]$est
  m1_beta3 = m1_betas[grepl(paste0(interaction), m1_betas$param),]$est
  
  m1_var = m1_param[grepl(paste0("Variances"), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_var1 = m1_var[grepl(paste0(c("^"),exogenous,c("$")), m1_var$param),]$est
  m1_res = m1_var[grepl(paste0(c("^"),endogenous,c("$")), m1_var$param),]$est
  m1_cov = m1_param[grepl(paste0(c("WITH")), m1_param[, grepl("paramHeader", names(m1_param))]),]$est

  #Get Observed Variable Variances
  txt2 = readLines(M1)
  df2 <- data.frame(txt2, paragraph = cumsum(txt2 == "")) 
  keep_paragraph2 <- df2[grep("Covariances", df2[, "txt2"]), "paragraph"] 
  df2 <- df2[df2$paragraph %in% keep_paragraph2,]
  df2 = df2[grepl(moderator,df2[,1]),]
  m1_var2 = as.numeric(str_sub(last(df2[,1]), -10))
  
  m1_Rsq = ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2)) / ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2) + m1_res)
  
  Interaction_Rsq = m1_Rsq - m0_Rsq
  return(Interaction_Rsq)
  
}

LOIR(M0 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_nointeraction.out", 
    M1 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_interaction.out", 
    endogenous="NER", exogenous = "DEGB", moderator = "GENDER", interaction= "INT")  



##############################################################################################
############################## STEP TWO - THREE WAY INTERACTION ################################
##############################################################################################

#############################
#Latent & Observed Interaction
##############################
endogenous = "NER"
exogenous = "DEGB"
moderator = "GENDER"
interaction = "INT"
M0 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_nointeraction.out"
M1 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_interaction.out"

M0 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_3waynointeraction"
M1 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_3wayinteraction"


LOIR = function(M0, M1, endogenous, exogenous, moderator, interaction){
  #M0 model
  m0_param = as.data.frame(readModels(M0, what="parameters")$parameters$unstandardized)
  m0_betas = m0_param[grepl(paste0(endogenous,c("."), c("ON")), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_beta1 = m0_betas[grepl(paste0(exogenous), m0_betas$param),]$est
  m0_beta2 = m0_betas[grepl(paste0(moderator), m0_betas$param),]$est
  m0_var = m0_param[grepl(paste0("Variances"), m0_param[, grepl("paramHeader", names(m0_param))]),]
  m0_var1 = m0_var[grepl(paste0(c("^"),exogenous,c("$")), m0_var$param),]$est
  m0_res = m0_var[grepl(paste0(c("^"),endogenous,c("$")), m0_var$param),]$est
  #Observed Variable Variances
  txt1 = readLines(M0)
  df1 <- data.frame(txt1, paragraph = cumsum(txt1 == "")) 
  keep_paragraph1 <- df1[grep("Covariances", df1[, "txt1"]), "paragraph"] 
  df1 <- df1[df1$paragraph %in% keep_paragraph1,] 
  m0_var2 = as.numeric(str_sub(last(df1[,1]), -10))
  m0_Rsq= ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2) / ((m0_beta1)^2*m0_var1 + (m0_beta2)^2*m0_var2 + 2*m0_beta1*m0_beta2 + m0_res)
  
  #M1 Model
  m1_param = as.data.frame(readModels(M1, what="parameters")$parameters$unstandardized)
  m1_betas = m1_param[grepl(paste0(endogenous,c("."), c("ON")), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_beta1 = m1_betas[grepl(paste0(exogenous), m1_betas$param),]$est
  m1_beta2 = m1_betas[grepl(paste0(moderator), m1_betas$param),]$est
  m1_beta3 = m1_betas[grepl(paste0(interaction), m1_betas$param),]$est
  
  m1_var = m1_param[grepl(paste0("Variances"), m1_param[, grepl("paramHeader", names(m1_param))]),]
  m1_var1 = m1_var[grepl(paste0(c("^"),exogenous,c("$")), m1_var$param),]$est
  m1_res = m1_var[grepl(paste0(c("^"),endogenous,c("$")), m1_var$param),]$est
  m1_cov = m1_param[grepl(paste0(c("WITH")), m1_param[, grepl("paramHeader", names(m1_param))]),]$est
  
  #Observed Variable Variances
  txt2 = readLines(M1)
  df2 <- data.frame(txt2, paragraph = cumsum(txt2 == "")) 
  keep_paragraph2 <- df2[grep("Covariances", df2[, "txt2"]), "paragraph"] 
  df2 <- df2[df2$paragraph %in% keep_paragraph2,] 
  m1_var2 = as.numeric(str_sub(last(df2[,1]), -10))
  
  m1_Rsq = ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2)) / ((m1_beta1)^2*m1_var1 + (m1_beta2)^2*m1_var2 + 2*m1_beta1*m1_beta2 + m1_beta3^2*(m1_var1*m1_var2 + m1_cov^2) + m1_res)
  
  Interaction_Rsq = m1_Rsq - m0_Rsq
  return(Interaction_Rsq)
  
}

LOIR(M0 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_nointeraction.out", 
     M1 = "J:\\Research Studies\\RFunctionDevelopment\\model_lvandobserve_interaction.out", 
     endogenous="NER", exogenous = "DEGB", moderator = "GENDER", interaction= "INT")  

















