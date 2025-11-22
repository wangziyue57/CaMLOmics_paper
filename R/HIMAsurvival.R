#High-dimensional Mediation Analysis in Survival Models
# load R packages
library(ncvreg)
library(ggm)

# SIS for alpha
sis_alpha <- function(X, M, COV, p){
  s_alpha <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s1 <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s2 <- summary(fit)$coef[2]           #coefficients for alpha
    s3 <- summary(fit)$coef[2,4]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_alpha)))
  alpha_sis <- t(dat)
  colnames(alpha_sis) = colnames(M)
  return(s_alpha=alpha_sis)
}

# SIS for beta (the regression Y~X+M)
sis_beta <- function(X, M, Y, COV, p){
  s_beta <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(Y=Y, M = M[, j], X = X)
    } else {
      MX <- data.frame(Y=Y, M = M[, j], X = X, COV = COV)
    }
    fit <- survival::coxph(Y ~., data = MX)
    s1 <- fit$var[1,1]                   #var for alpha
    s2 <- summary(fit)$coef[1]           #coefficients for alpha
    s3 <- summary(fit)$coef[1,5]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_beta)))
  beta_sis <- t(dat)
  colnames(beta_sis) = colnames(M)
  return(s_beta=beta_sis)
}

#main function
hmas <- function(X, Y, M, COV,
                 penalty = c("MCP", "SCAD", "lasso"),
                 path = c('MY', 'MX', "both"),
                 topN = NULL,
                 verbose = TRUE, 
                 ...) {
  penalty <- match.arg(penalty)
  
  n <- nrow(M)
  p <- ncol(M)
  
  if (is.null(topN)) {
    d <- ceiling(2*n/log(n))  #the top d mediators that associated with exposure
  } else {
    d <- topN  
  }
  
  if(verbose) message("Step 1: Prelimenary Screening...", "     (", Sys.time(), ")")
  
  if (path == 'MY'){
    beta_s <- sis_beta(X=X, M=M, Y=Y, COV=COV, p=ncol(M))
    SIS_beta <- beta_s[3,]
    SIS_beta_sort <- sort(SIS_beta)
    ID_SIS <- which(SIS_beta <= SIS_beta_sort[d])  # the index of top d significant mediators (Y~X+M)
  } else if (path == 'MX'){
    alpha_s <- sis_alpha(X, M, COV, p)
    SIS_alpha <- alpha_s[3,]
    SIS_alpha_sort <- sort(SIS_alpha)
    ID_SIS <- which(SIS_alpha <= SIS_alpha_sort[d])  # the index of top d significant mediators (M~X)
  } else{ 
    alpha_s <- sis_alpha(X, M, COV, p)
    SIS_alpha <- alpha_s[2,]
    
    beta_s <- sis_beta(X=X, M=M, Y=Y, COV=COV, p=ncol(M))
    SIS_beta <- beta_s[2,]
    
    SIS_ab <- SIS_alpha*SIS_beta
    ID_SIS  <- which(-abs(SIS_ab) <= sort(-abs(SIS_ab))[d]) # top d mediators with largest effect of alpha*beta
  }
  
  M_SIS <- M[, ID_SIS]
  
  if(verbose) cat("Top", length(ID_SIS), "mediators selected.", "\n")
  
  XM <- cbind(M_SIS, X)
  C_M <- colnames(XM)
  
  
  if(verbose) message("Step 2: Penalized Variable Selection (", penalty, ") ...", "  (", 
                      Sys.time(), ")")
  
  if (is.null(COV)) {
    fit <- ncvsurv(XM, Y,
                   penalty = penalty,
                   penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
    
    # cvfit <- cv.ncvsurv(XM, Y, 
    #                     penalty = penalty, 
    #                     penalty.factor = c(rep(1, ncol(M_SIS)), 0), 
    #                     nfolds = 10, ...)
  
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(XM, COV)
    
    fit <- ncvsurv(XM_COV, Y,
                   penalty = penalty,
                   penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))), ...)
    
    # cvfit <- cv.ncvsurv(XM_COV, Y, 
    #                     penalty = penalty, 
    #                     penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))), 
    #                     nfolds = 10, ...)
  }
  
  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) cat("lambda selected: ", lam, "\n")
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[1:length(ID_SIS)]
  ID_p_non <- which(est != 0)
  
  # lam <- cvfit$lambda.min
  # if(verbose) cat("lambda selected: ", lam, "\n")
  # Coefficients <- coef(cvfit)
  # est <- Coefficients[1:length(ID_SIS)]
  # ID_p_non <- which(est != 0)
  
  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", length(ID_p_non), "\n")
  
  if(length(ID_p_non) >= 1){
    beta_p <- est[ID_p_non]  # The non-zero MCP estimators of beta
    ID_p <- ID_SIS[ID_p_non]  # The index of the ID of non-zero beta in the Cox regression
    MCP_M <- names(ID_p_non)
    
    # the estimator for alpha
    alpha_est <- alpha_s[2, ID_p]   
    var_alpha <- alpha_s[1, ID_p]
    P_alpha <- alpha_s[3, ID_p]
    
    # # true estimation of alpha and beta from MCP penalized regression
    # beta_t <- beta_p
    # alpha_t <- alpha_est
    # ab_true <- alpha_t * beta_t
    
    # direct effect of X->Y given M and COV
    DE <- Coefficients[length(ID_SIS) + 1]
    # summary(cox_model)$coefficients[(length(ID_p)+1), 1]
    
    
    if(verbose) message("Step 3: The adjusted significance test ...", 
                        "     (", Sys.time(), ")")
    
    # refit the cox model for varaince estimation and inference (not the optimal solution!)
    if (is.null(COV)) {
      YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X, COV = COV)
    }
    cox_model <- survival::coxph(Y ~ ., data = YMX)
    
    # direct effect of X->Y given M and COV
    # DE <- Coefficients[length(ID_SIS) + 1]
    DE <- summary(cox_model)$coefficients[(length(ID_p)+1), 1]
    
    # the re-estimator of beta
    beta_est <- summary(cox_model)$coefficients[1: length(ID_p)]     
    var_beta <- (summary(cox_model)$coefficients[1:length(ID_p),3])^2
    P_beta <- summary(cox_model)$coefficients[1: length(ID_p), 5]  
    
    # the estimator of alpha*beta
    ab_est <- alpha_est * beta_est   
    
    # var(alpha*beta)
    var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2) + var_alpha * var_beta
    # var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2) 
    
    # confidence interval
    conf_low <- ab_est - 1.96 * sqrt(var_ab)
    conf_up <- ab_est + 1.96 * sqrt(var_ab)
    
    # option 1: sobel test for alpha*beta
    s.test <- abs(ab_est)/sqrt(var_ab)   #z-score for sobel test
    P_sobel <- 2 * (1-pnorm(s.test))     #p-value of sobel test
    P_fdr_sobel <- p.adjust(P_sobel, 'fdr', length(ID_p))
    P_fdr_sobel[P_fdr_sobel > 1] <- 1
    
    # option 2: joint test with mutiple-testing procedure
    PA <- cbind(t(P_alpha), t(P_beta))
    P_value <- apply(PA,1,max)  # the joint p-values for SIS variable
    
    ## the multiple-testing  procedure
    if(nrow(PA) > 1){
      N0 <- dim(PA)[1]*dim(PA)[2]
      
      input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
      nullprop <- HDMT::null_estimation(input_pvalues)
      fdrcut  <- HDMT::fdr_est(nullprop$alpha00,
                               nullprop$alpha01,
                               nullprop$alpha10, 
                               nullprop$alpha1,
                               nullprop$alpha2,
                               input_pvalues,
                               exact=0)
    }else{
      fdrcut <- P_value
    }
    
    results <- data.frame(alpha = alpha_est, beta = beta_est,
                          `alpha_est*beta_est` = ab_est, 
                          # beta_t = beta_t, `alpha_t*beta_t` = ab_true, 
                          conf_low=conf_low, conf_up=conf_up,
                          P_fdr_sobel=P_fdr_sobel,
                          P_joint=P_value, P_fdr_joint=fdrcut,
                          var_ab=var_ab, var_alpha=var_alpha, var_beta=var_beta,
                          DE, check.names = FALSE)
  }else{
    MCP_M <- NULL
    results <- data.frame(alpha = NA, beta = NA,
                          `alpha_est*beta_est` = NA, 
                          # beta_t = NA, `alpha_t*beta_t` = NA, 
                          conf_low=NA, conf_up=NA,
                          P_fdr_sobel=NA,
                          P_joint=NA, P_fdr_joint=NA,
                          var_ab=NA, var_alpha=NA, var_beta=NA,
                          DE=NA, check.names = FALSE)
  }
 
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  return(list(C_M, MCP_M, results))
}

