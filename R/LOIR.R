#' @title LOIR
#'
#' @description
#' R square of one two-way interactions between one latent factor and one observed covariate in structural equation model
#'
#' @param M0 the mplus output of the model without interaction;
#' @param M1 the mplus output of the model with interaction;
#' @param endogenous one endogenous latent factor
#' @param exogenous one exogenous latent factor
#' @param moderator one observed covariate
#' @param interaction one two-way interaction that are defined by the XWITH function in the mplus
#'
#' @author {Lu Qin, Howard University, \email{lu.qin@howard.edu}}
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@@uiowa.edu}}
#' @export


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
  m1_cov = m1_param[grepl(paste0(exogenous,c("."), c("WITH")), m1_param[, grepl("paramHeader", names(m1_param))]),]$est

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
