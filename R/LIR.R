#' @title LIR
#'
#' @description
#' R square of one two-way interactions between latent factors in structural equation model

#' @param M0 the mplus output of the model without interaction;
#' @param M1 the mplus output of the model with interaction;
#' @param endogenous one endogenous latent factor
#' @param exogenous one exogenous latent factor
#' @param moderator one latent moderator factor
#' @param interaction one two-way interaction that are defined by the XWITH function in the mplus
#'
#' @author {Lu Qin, Howard University, \email{lu.qin@howard.edu}}
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#' @export


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
