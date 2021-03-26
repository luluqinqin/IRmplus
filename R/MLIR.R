#' @title MLIR
#'
#' @description
#' R square of one or multiple two-way interactions between latent factors in structural equation model

#' @param Mplus_Output the mplus output of the SEM model with one or more two-way interactions
#' @param endogenous the endogenous latent factor
#' @param exogenous the exogenous latent factors that are regressed on the endogenous latent factor
#' @param Interaction_factor the latent factors that are used to define the two-way interactions
#' @param interaction_code interactions that are defined by the XWITH function in the mplus
#'
#' @author {Lu Qin, Howard University, \email{lu.qin@howard.edu}}
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@@uiowa.edu}}
#' @export


MLIR = function(Mplus_Output, endogenous, exogenous, Interaction_factor, interaction_code) {

  N = length(Interaction_factor)
  n = length(exogenous)

  int = matrix(0, (N-1),N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      int[i,j] = c(paste0(Interaction_factor[i],Interaction_factor[j]))
    }
  }
  interaction = list()
  for (i in 1:nrow(int)){
    interaction[[i]] = subset(int[i,],int[i,] != 0)
  }
  interaction = as.character(unlist(interaction))

  param = as.data.frame(readModels(Mplus_Output, what="parameters")$parameters$unstandardized)
  var = param[grepl(paste0("Variances"), param[, grepl("paramHeader", names(param))]),]
  var_end  = var[grepl(paste0(c("^"),endogenous,c("$")), var$param),]$est

  betas = param[grepl(paste0(endogenous,c("."), c("ON")), param[, grepl("paramHeader", names(param))]),]
  beta_exo = matrix(0,n,1)
  var_exo = matrix(0,n,1)
  for (i in 1:n){
    beta_exo[i] = betas[grepl(paste0(c("^"),exogenous[i],c("$")), betas$param),]$est
    var_exo[i] = var[grepl(paste0(c("^"),exogenous[i],c("$")), var$param),]$est
  }
  beta_int = matrix(0,(N*(N-1)/2),1)
  for (j in 1:(N*(N-1)/2)){
    beta_int[j] = betas[grepl(paste0(c("^"),interaction_code[j],c("$")), betas$param),]$est
  }
  beta = rbind(beta_exo, beta_int)
  rownames(beta) = c(exogenous,interaction)

  cov = param[grepl(paste0(c("WITH")), param[, grepl("paramHeader", names(param))]),]
  vcov_exo = matrix(0, n, n)
  for (p in 1:(n-2)){
    cov_exo = cov[grepl(paste0(c("^"),exogenous[p],c("$")), cov$param),]
    for (h in (1+p):n){
      cov_exo1 = cov_exo[grepl(paste0(exogenous[h],c(".WITH")), cov_exo$paramHeader),]$est
      vcov_exo[p,h] = cov_exo1
    }
  }
  vcov_exo[(n-1),n] = cov[grepl(paste0(c("^"),exogenous[n-1],c("$")), cov$param),]$est
  diag(vcov_exo) = var_exo
  vcov_exo[lower.tri(vcov_exo)] = vcov_exo[upper.tri(vcov_exo)]
  vcov_exo
  rownames(vcov_exo) = exogenous
  colnames(vcov_exo) = exogenous

  var_int = matrix(0, N, N)
  for (p in 1:(N-2)){
    for (h in (1+p):N){
      var_int[p,h] = var_exo[p] * var_exo[h] * vcov_exo[p,h]^2

    }
  }
  var_int[(N-1),N] = var_exo[(N-1)] * var_exo[N] * (vcov_exo[(N-1),N])^2
  D = list()
  for (i in 1:(nrow(var_int)-1)){
    D[[i]] = subset(var_int[i,],var_int[i,] != 0)
  }
  var_int = as.numeric(unlist(D))

  if (N >2) {

    A = matrix(0, length(rep(1:(N-1))), N)
    for (i in 1:(N-1)){
      for (j in (i+1):N) {
        A[i,j] = c(paste0(i,j))
      }
    }
    C = list()
    for (i in 1:nrow(A)){
      C[[i]] = subset(A[i,],A[i,] != 0)
    }

    IJ = as.numeric(unlist(C))
    I = as.numeric(substr(IJ,1,1))
    J = as.numeric(substr(IJ,2,2))
    H = as.numeric(substr(IJ,1,1))
    G = as.numeric(substr(IJ,2,2))

    vcov_int = matrix(0, length(IJ), length(IJ))
    for (k in 1:(length(IJ)-1)){
      for (p in (k+1):length(IJ)){
        ij = noquote(paste0(I[k],J[k]))
        hg = noquote(paste0(H[p],G[p]))
        vcov_int[k,p] = vcov_exo[I[k],H[p]]*vcov_exo[J[k],G[p]] + vcov_exo[I[k],G[p]]*vcov_exo[J[k],H[p]]
      }
    }
    diag(vcov_int) = var_int
    vcov_int[lower.tri(vcov_int)] = vcov_int[upper.tri(vcov_int)]
    vcov_int
    rownames(vcov_int) = interaction
    colnames(vcov_int) = interaction

    ve_1 = 0
    for (i in 1:n){
      temp = beta[exogenous[i],]^2 * vcov_exo[i,i]
      ve_1 = ve_1 + temp
    }
    ve_1

    ve_2 = 0
    for (i in 1: (n-1)){
      for (j in (i+1):n){
        temp = 2 * beta[exogenous[i],]*beta[exogenous[j],] * vcov_exo[i,j]
        ve_2 = ve_2 + temp
      }
    }
    ve_2

    ve_3 = 0
    for (i in 1:length(interaction)){
      temp = beta[interaction[i],]^2 * vcov_int[i,i]
      ve_3 = ve_3 + temp
    }
    ve_3

    ve_4 = 0
    for (i in 1: (length(interaction)-1)){
      for (j in (i+1):(length(interaction))){
        temp = 2 * beta[interaction[i],]*beta[interaction[j],] * vcov_int[i,j]
        ve_4 = ve_4 + temp
      }
    }
    ve_4

    R1_raw = round((ve_1 + ve_2 + ve_3 + ve_4) /(ve_1 + ve_2 + ve_3 + ve_4 + var_end),4)
    R1 = round(((ve_1 + ve_2 + ve_3 + ve_4) /(ve_1 + ve_2 + ve_3 + ve_4 + var_end))*100,2)

    R_int_raw = round((ve_3 + ve_4) / (ve_1 + ve_2 + ve_3 + ve_4 + var_end),4)
    R_int = round(((ve_3 + ve_4) / (ve_1 + ve_2 + ve_3 + ve_4 + var_end))*100,2)

  }

  else {

    ve_1 = 0
    for (i in 1:n){
      temp = beta[exogenous[i],]^2 * vcov_exo[i,i]
      ve_1 = ve_1 + temp
    }
    ve_1

    ve_2 = 0
    for (i in 1: (n-1)){
      for (j in (i+1):n){
        temp = 2 * beta[exogenous[i],]*beta[exogenous[j],] * vcov_exo[i,j]
        ve_2 = ve_2 + temp
      }
    }
    ve_2

    ve_3 = 0
    for (i in 1:length(interaction)){
      temp = beta[interaction[i],]^2 * var_int
      ve_3 = ve_3 + temp
    }
    ve_3

    R1_raw = round((ve_1 + ve_2 + ve_3) /(ve_1 + ve_2 + ve_3 + var_end),4)
    R1 = round(((ve_1 + ve_2 + ve_3) /(ve_1 + ve_2 + ve_3 + var_end))*100,2)

    R_int_raw = round((ve_3 ) / (ve_1 + ve_2 + ve_3 + var_end),4)
    R_int = round(((ve_3 ) / (ve_1 + ve_2 + ve_3 + var_end))*100,2)

  }


  return(c(paste0("R2 = ", R1_raw, "; Variance explained by the main and interaction effect is ", R1, "%"),
           paste0("R2 = ", R_int_raw, "; Variance explained by the multiple two-way interaction is ", R_int, "%")))

}






