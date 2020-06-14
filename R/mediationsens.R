#' Simulation-Based Sensitivity Analysis Table for Causal Mediation Analysis
#' 
#' This function is used to evaluate the sensitivity of the estimated natural direrct and indirect effects to potential violations of the ignorability assumptions. It estimates the natural direct and indirect effects after adjusting for an unmeasured pretreatment confounder, U, with a specified degree of confounding. In a randomized experiment, the degree of confounding is evaluated via two sensitivity parameters, the coefficient of U in the mediator model and that in the outcome model, given the specified prior distribution of U. When the treatment is not randomized, an additional sensitivity parameter is introduced -- the coefficient of U in the treatment model. The mediator, outcome, and unmeasured pretreatment confounder could be either binary or continuous, while the treatment has to be binary. If the treatment is randomized, the sensitivity analysis method is applicable to various causal mediation analysis methods, including the ratio of mediator probability weights (RMPW) method, the regression-based method proposed by VanderWeele (2010), and the simulation-based method proposed by Imai et al. (2010a, 2010b). If the treatment is not randomized, the method is applicable to the RMPW method.
#'
#' @param y The name of the outcome variable (string).
#' @param t The name of the treatment variable (string). The treatment variable has to be binary.
#' @param t.rand TRUE if the treatment is randomized, and FALSE if not. The default is TRUE.
#' @param m The name of the mediator variable (string).
#' @param covariates A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. Transform each categorical variable into multiple binary indicators before specifying X.
#' @param scale.y The scale of the outcome variable (string). Can be "continuous" or "binary".
#' @param scale.m The scale of the mediator variable (string). Can be "continuous" or "binary".
#' @param scale.U The scale of the unobserved pretreatment confounder (string). Can be "continuous" or "binary".
#' @param b.t The value of the sensitivity parameter that represents the conditional association between the unobserved pretreatment confounder and the treatment. The default is NULL, because when the treatment is randomized, there is no need to specify this sensitivity parameter.
#' @param range.b.m The range of the sensitivity parameter that represents the conditional association between the unobserved pretreatment confounder and the mediator. The default is c(-2, 2).
#' @param range.b.y The range of the sensitivity parameter that represents the conditional association between the unobserved pretreatment confounder and the outcome. The default is c(-1, 1).
#' @param grid.b.m The horizontal dimension of the grid. The default is 5. Increase the number for more smooth curves. 
#' @param grid.b.y The vertical dimension of the grid. The default is 5. Increase the number for more smooth curves.
#' @param p.u The prior probability of the unobserved pretreatment confounder if it is binary. The default is 0.5.
#' @param sigma.u The standard deviation of the prior distribution of the unobserved pretreatment confounder if it is continuous. The default is 1.
#' @param iter The number of iterations in the stochastic EM algorithm for generating the unobserved pretreatment confounder. The default is 10. 
#' @param nsim The number of random draws of the unobserved pretreatment confounder in each cell. The default is 40. Increase the number for more smooth curves.
#' @param est The method to be used in the estimation of the causal mediation effects. "rmpw" uses the weighting-based method. It estimates the natural direct and indirect effects through mean contrasts of the outcome. When the outcome is coninuous, it controls for the pretreatment confounders in the outcome model for improving efficiency. When the outcome is binary, it cannot control for the pretreatment confounders. In addition to the scenario with a binary mediator and a continuous outcome, we use bootstrapping to estimate their standard errors and test their significance. "regression" uses the regression-based method developed by VanderWeele and Vansteelandt (2009). It is only applicable to continuous outcomes, because when the outcome is binary, the natural direct and indirect effects are identified on the odds ratio scale, rendering results difficult to report. Because the estimator of the natural direct effect or the natural indirect effect is a complex combination of regression coefficients and may not follow a normal distribution, we use bootstrapping to estimate their standard errors and test their significance. "simulation" uses the simulation-based method developed by Imai et al. (2010a, 2010b), and the mediation package is adopted. For a binary mediator or outcome, the link function is always specified to be logit. If t.rand = FALSE, i.e the treatment is not randomized, only rmpw is applicable.
#' @param B The number of bootstrapped samples in the causal mediation analysis if bootstrapping is used for estimating standard errors and testing significance.
#' @param data The data set for analysis.
#' @return A list containing
#' \item{results.old}{the causal effect estimates in the original analysis.}
#' \item{X.coef.plot}{the coefficients of the observed pretreatment covariates in the mediator and outcome models (which are used to calibrate the strength of the sensitivity parameters in the sensitivity plots).}
#' \item{range.b.m, range.b.y}{same as specified in the sens function.}
#' \item{b.m.all, b.y.all}{all the specific values of the sensitivity parameters for constructing the grids.}
#' \item{results.new}{tables for the causal effect estimates, standard errors, and p values in all the grids, after the adjustment of the unobserved pretreatment confounder.}
#' @author Xu Qin and Fan Yang
#' @references Qin, X., & Yang, F. (2020). Simulation-Based Sensitivity Analysis for Causal Mediation Studies.
#' @references VanderWeele, T. J. (2010). Bias formulas for sensitivity analysis for direct and indirect effects. Epidemiology, 21, 540. 
#' @references Imai, K., Keele, L., & Tingley, D. (2010a). A general approach to causal mediation analysis. Psychological Methods, 15, 309. 
#' @references Imai, K., Keele, L., & Yamamoto, T. (2010b). Identification, inference and sensitivity analysis for causal mediation effects. Statistical Science, 25, 51-71. 
#' @export
#' @importFrom stats as.formula binomial coef fitted glm lm pnorm dnorm rnorm rbinom predict model.matrix integrate qnorm sd sigma var
#' @importFrom distr AbscontDistribution r
#' @importFrom mediation mediate
#' @examples 
#' library(mediation)
#' data(jobs)
#' # Generate a binary outcome based on the median of the continuous outcome.
#' jobs$depress2_dich[jobs$depress2 >= median(jobs$depress2)] = 1
#' jobs$depress2_dich[jobs$depress2 < median(jobs$depress2)] = 0
#' 
#' sens.results = sens(y = "depress2", t = "treat", m = "job_seek", covariates = c("econ_hard", "depress1", "sex", "age"), scale.y = "continuous", scale.m = "continuous", scale.U = "binary", range.b.m = c(-2, 2), range.b.y = c(-1, 1), grid.b.y = 2, grid.b.m = 1, p.u = 0.5, iter = 10, nsim = 10, est = "rmpw", B = 2, data = jobs)


sens = function(y, t.rand = TRUE, t, m, covariates, scale.y, scale.m, scale.U, b.t = NULL, range.b.m = c(-2, 2), range.b.y = c(-1, 1), grid.b.y = 5, grid.b.m = 5, p.u = 0.5, sigma.u = 1, iter = 10, nsim = 40, est, B, data){
  X = covariates
  # Weighting-based method - RMPW
  rmpw = function(y, t, t.rand = TRUE, m, X, U = FALSE, scale.y, scale.m, b.m = NULL, b.t = NULL, B, data){
    if(scale.m == "binary" & scale.y == "continuous"){
      if(t.rand == TRUE){
        if(!U){
          l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), data = data, family = binomial)
          data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X])) %*% coef(l)))
          data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X])) %*% coef(l)))
          p = predict(l, data, type = "response")
          X0.matrix = as.matrix(cbind(1, 0, data[, X]))
          X1.matrix = as.matrix(cbind(1, 1, data[, X]))
          X.matrix = as.matrix(cbind(1, data[, c(t, X)]))
        }
        if(U){
          l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), offset = b.m * U, data = data, family = binomial)
          data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
          data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
          p = 1/(1 + exp(-as.matrix(cbind(1, data[, c(t, X)], U = data$U)) %*% c(coef(l), U = b.m)))
          X0.matrix = as.matrix(cbind(1, 0, data[, X]))
          X1.matrix = as.matrix(cbind(1, 1, data[, X]))
          X.matrix = as.matrix(cbind(1, data[, c(t, X)]))
        }
        h1 = X.matrix * as.vector((data[, m] - p))
        G11 = t(X.matrix) %*% diag(as.vector(p * (1 - p))) %*% X.matrix
        
        data$rmpw[data[, t] == 1 & data[, m] == 1] = (data$p0/data$p1)[data[, t] == 1 & data[, m] == 1]
        data$rmpw[data[, t] == 1 & data[, m] == 0] = ((1 - data$p0)/(1 - data$p1))[data[, t] == 1 & data[, m] == 0]
        data$rmpw[data[, t] == 0] = 1
        
        data_ctrl = data[data[, t] == 0, ]
        data_tr = data[data[, t] == 1, ]
        data_tr_dup = NULL
        for (j in 1:ncol(data_tr)) {
          data_tr_dup = cbind(data_tr_dup, rep(data_tr[, j], rep(2, nrow(data_tr))))
        }
        colnames(data_tr_dup) = colnames(data_ctrl)
        data_tr_dup = as.data.frame(data_tr_dup)
        data_tr_dup$rmpw[seq(2, length(data_tr_dup[, 1]), 2)] = 1
        d1 = c(rep(c(0, 1), nrow(data_tr)), rep(0, length(data_ctrl[, 1])))
        data_dup = cbind.data.frame(d1, rbind(data_tr_dup, data_ctrl))
        
        z = model.matrix(lm(as.formula(paste(m, "~", t, "+ d1 +", paste(X, collapse="+"))), data = data_dup))
        w = diag(data_dup$rmpw)
        if(!U){
          l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw)
          beta = coef(l.y)
        }
        if(U){
          l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw)
          beta = coef(l.y)
        }
        delta = beta[1]
        delta_counter = beta[1] + beta[2]
        delta_e = beta[1] + beta[2] + beta[3]
        lamda = beta[4:(length(beta))]
        
        if(!U){
          x.outcome = as.matrix(data[, X]) # x in the outcome model. It will be used in the moment function. 
        }
        if(U){
          x.outcome = as.matrix(data[, c(X, "U")]) # x in the outcome model. It will be used in the moment function. U is included, because its coefficient is not a sensitivity parameter and needs to be estimated.
        }
        
        n = nrow(data)
        outcome = as.matrix(data[,c(y, y, y)])
        outcome = matrix(0, n, 3)  
        outcome[,1] = data[, y]
        outcome[,3] = data[, y]
        outcome[,2] = data[, y]
        dimnames(outcome)[[2]] <- c( "delta","delta_*","delta_e")
        
        w = matrix(0, n, 3)
        w[,1] = (data[, t] == 0)
        w[,3] = (data[, t] == 1)
        w[,2] = (data[, t] == 1)*data$rmpw
        dimnames(w)[[2]] <- c( "delta","delta_*","delta_e")
        
        deltaEstimate = c(delta, delta_counter, delta_e)
        wSum = NULL
        hDelta = NULL #score function for delta, delta^*, delta^3  
        for (j in 1:3){
          delta = deltaEstimate[j]    
          hDelta = cbind(hDelta, (outcome[,j] - x.outcome %*% lamda - delta)*w[,j])
          wSum = cbind(wSum, sum(w[,j]))
        }  
        dimnames(wSum)[[2]] <- c( "delta","delta_*","delta_e")
        G22 = diag(as.numeric(wSum))
        dimnames(hDelta)[[2]] <- c( "delta","delta_*","delta_e")
        hX = (hDelta[ ,1] + hDelta[ ,2] + hDelta[ ,3]) * x.outcome
        hCombined = cbind(h1, hDelta, hX)
        B0 = t(hCombined)%*%hCombined
        
        #we calculate G21 matrix
        
        dp0_dbeta = as.vector(data$p0 * (1 - data$p0)) * X0.matrix
        dp1_dbeta = as.vector(data$p1 * (1 - data$p1)) * X1.matrix
        dw_dbeta = data[, m] * (dp0_dbeta * as.vector(data$p1) - as.vector(data$p0) * dp1_dbeta)/(as.vector(data$p1) ^ 2) + (1 - data[, m]) * (-dp0_dbeta * (1 - as.vector(data$p1)) + (1 - as.vector(data$p0)) * dp1_dbeta)/((1 - as.vector(data$p1)) ^ 2)
        dh.star_dbeta = -(data[, y] - as.numeric(x.outcome %*% lamda) - deltaEstimate[2]) * data[, t] * dw_dbeta
        nP = ncol(X.matrix)
        G21 = matrix(0, 3, nP)
        for(j in 1:ncol(G21)){
          G21[2, j] = sum(dh.star_dbeta[, j])
        }
        
        G23 = matrix(0, 3, length(lamda))
        G23[1, ] = apply((1 - data[, t]) * x.outcome, 2, sum)
        G23[2, ] = apply(data[, t] * w[, 2] * x.outcome, 2, sum)
        G23[3, ] = apply(data[, t] * x.outcome, 2, sum)
        G32 = t(G23)
        
        G31 = t(x.outcome) %*% diag(-data[, t] * (data[, y] -  as.numeric(x.outcome %*% lamda) - deltaEstimate[2])) %*% dw_dbeta
        
        G33 = t(x.outcome) %*% diag(1 + data[, t] * w[, 2]) %*% x.outcome
        
        A0 = matrix(0, (nP+3+ncol(x.outcome)), (nP+3+ncol(x.outcome)))
        A0[(nP+1):(nP+3), (nP+1):(nP+3)] = G22
        A0[1:(nP), 1:(nP)] = G11
        A0[(nP+1):(nP+3), 1:(nP)] = G21
        A0[(nP+3+1):(nP+3+ncol(x.outcome)), 1:(nP)] = G31
        A0[(nP+3+1):(nP+3+ncol(x.outcome)), (nP+1):(nP+3)] = G32
        A0[(nP+3+1):(nP+3+ncol(x.outcome)), (nP+3+1):(nP+3+ncol(x.outcome))] = G33
        A0[(nP+1):(nP+3), (nP+3+1):(nP+3+ncol(x.outcome))] = G23
        
        v_hw = solve(A0) %*% B0 %*% t(solve(A0))
        
        v_hw_ctrl_counterfacal_tr = v_hw[(nP+1):(nP+3), (nP+1):(nP+3)]  
        NDE = as.numeric(deltaEstimate[2] - deltaEstimate[1])
        NIE = as.numeric(deltaEstimate[3] - deltaEstimate[2])
        SE_NDE = sqrt(v_hw_ctrl_counterfacal_tr[2,2]+v_hw_ctrl_counterfacal_tr[1,1] - 2*v_hw_ctrl_counterfacal_tr[1,2])
        SE_NIE = sqrt(v_hw_ctrl_counterfacal_tr[3,3]+v_hw_ctrl_counterfacal_tr[2,2] - 2*v_hw_ctrl_counterfacal_tr[2,3])
        
        P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
        P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
      }
      
      if(t.rand == FALSE){
        est = function(data){
          if(!U){
            l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), data = data, family = binomial)
            data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X])) %*% coef(l)))
            data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X])) %*% coef(l)))
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X])) %*% coef(l.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X])) %*% coef(l.t)))
          }
          if(U){
            l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), offset = b.m * U, data = data, family = binomial)
            data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
            data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), offset = b.t * U, data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X], U = data[data[, t] == 1, "U"])) %*% c(coef(l.t), U = b.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X], U = data[data[, t] == 0, "U"])) %*% c(coef(l.t), U = b.t)))
          }
         
          data$rmpw[data[, t] == 1 & data[, m] == 1] = (data$p0/data$p1)[data[, t] == 1 & data[, m] == 1]
          data$rmpw[data[, t] == 1 & data[, m] == 0] = ((1 - data$p0)/(1 - data$p1))[data[, t] == 1 & data[, m] == 0]
          data$rmpw[data[, t] == 0] = 1
          data$iptw[data[, t] == 1] = (mean(data[, t])/data$pt)[data[, t] == 1]
          data$iptw[data[, t] == 0] = ((1 - mean(data[, t]))/data$pt)[data[, t] == 0]
          
          data_ctrl = data[data[, t] == 0, ]
          data_tr = data[data[, t] == 1, ]
          data_tr_dup = NULL
          for (j in 1:ncol(data_tr)) {
            data_tr_dup = cbind(data_tr_dup, rep(data_tr[, j], rep(2, nrow(data_tr))))
          }
          colnames(data_tr_dup) = colnames(data_ctrl)
          data_tr_dup = as.data.frame(data_tr_dup)
          data_tr_dup$rmpw[seq(2, length(data_tr_dup[, 1]), 2)] = 1
          d1 = c(rep(c(0, 1), nrow(data_tr)), rep(0, length(data_ctrl[, 1])))
          data_dup = cbind.data.frame(d1, rbind(data_tr_dup, data_ctrl))
          
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
          beta = coef(l.y)
          
          return(list(NDE = beta[2], NIE = beta[3]))
        }
        
        est.ori = est(data)
        NDE = est.ori$NDE
        NIE = est.ori$NIE
        
        NDE.boot = NULL
        NIE.boot = NULL
        for(i in 1:B){
          data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
          est.boot = est(data.boot)
          NDE.boot = c(NDE.boot, est.boot$NDE)
          NIE.boot = c(NIE.boot, est.boot$NIE)
        }
        
        SE_NDE = sd(NDE.boot)
        SE_NIE = sd(NIE.boot)
        
        P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
        P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
      }
    }
    
    if(scale.m == "binary" & scale.y == "binary"){
      est = function(data){
        if(!U){
          l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), data = data, family = binomial)
          data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X])) %*% coef(l)))
          data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X])) %*% coef(l)))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X])) %*% coef(l.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X])) %*% coef(l.t)))
          }
        }
        if(U){
          l = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), offset = b.m * U, data = data, family = binomial)
          data$p1 = 1/(1 + exp(-as.matrix(cbind(1, 1, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
          data$p0 = 1/(1 + exp(-as.matrix(cbind(1, 0, data[, X], U = data$U)) %*% c(coef(l), U = b.m)))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), offset = b.t * U, data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X], U = data[data[, t] == 1, "U"])) %*% c(coef(l.t), U = b.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X], U = data[data[, t] == 0, "U"])) %*% c(coef(l.t), U = b.t)))
          }
        }
        
        data$rmpw[data[, t] == 1 & data[, m] == 1] = (data$p0/data$p1)[data[, t] == 1 & data[, m] == 1]
        data$rmpw[data[, t] == 1 & data[, m] == 0] = ((1 - data$p0)/(1 - data$p1))[data[, t] == 1 & data[, m] == 0]
        data$rmpw[data[, t] == 0] = 1
        if(t.rand == FALSE){
          data$iptw[data[, t] == 1] = (mean(data[, t])/data$pt)[data[, t] == 1]
          data$iptw[data[, t] == 0] = ((1 - mean(data[, t]))/data$pt)[data[, t] == 0]
        }
        
        data_ctrl = data[data[, t] == 0, ]
        data_tr = data[data[, t] == 1, ]
        data_tr_dup = NULL
        for (j in 1:ncol(data_tr)) {
          data_tr_dup = cbind(data_tr_dup, rep(data_tr[, j], rep(2, nrow(data_tr))))
        }
        colnames(data_tr_dup) = colnames(data_ctrl)
        data_tr_dup = as.data.frame(data_tr_dup)
        data_tr_dup$rmpw[seq(2, length(data_tr_dup[, 1]), 2)] = 1
        d1 = c(rep(c(0, 1), nrow(data_tr)), rep(0, length(data_ctrl[, 1])))
        data_dup = cbind.data.frame(d1, rbind(data_tr_dup, data_ctrl))
        
        if(t.rand == TRUE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw)
        }
        if(t.rand == FALSE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
        }
        beta = coef(l.y)
        
        return(list(NDE = beta[2], NIE = beta[3]))
      }
      
      est.ori = est(data)
      NDE = est.ori$NDE
      NIE = est.ori$NIE
      
      NDE.boot = NULL
      NIE.boot = NULL
      for(i in 1:B){
        data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
        est.boot = est(data.boot)
        NDE.boot = c(NDE.boot, est.boot$NDE)
        NIE.boot = c(NIE.boot, est.boot$NIE)
      }
      
      SE_NDE = sd(NDE.boot)
      SE_NIE = sd(NIE.boot)
      
      P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
      P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
    }
    
    if(scale.m == "continuous" & scale.y == "continuous"){
      est = function(data){
        if(!U){
          l = lm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), data = data)
          data$p1 = dnorm(data[, m], mean = as.matrix(cbind(1, 1, data[, X])) %*% coef(l), sd = sigma(l))
          data$p0 = dnorm(data[, m], mean = as.matrix(cbind(1, 0, data[, X])) %*% coef(l), sd = sigma(l))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X])) %*% coef(l.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X])) %*% coef(l.t)))
          }
        }
        
        if(U){
          l = lm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), offset = b.m * U, data = data)
          data$p1 = dnorm(data[, m], mean = as.matrix(cbind(1, 1, data[, X], U = data$U)) %*% c(coef(l), U = b.m), sd = sigma(l))
          data$p0 = dnorm(data[, m], mean = as.matrix(cbind(1, 0, data[, X], U = data$U)) %*% c(coef(l), U = b.m), sd = sigma(l))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), offset = b.t * U, data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X], U = data[data[, t] == 1, "U"])) %*% c(coef(l.t), U = b.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X], U = data[data[, t] == 0, "U"])) %*% c(coef(l.t), U = b.t)))
          }
        }
        
        data$rmpw[data[, t] == 1] = (data$p0/data$p1)[data[, t] == 1]
        data$rmpw[data[, t] == 0] = 1
        if(t.rand == FALSE){
          data$iptw[data[, t] == 1] = (mean(data[, t])/data$pt)[data[, t] == 1]
          data$iptw[data[, t] == 0] = ((1 - mean(data[, t]))/data$pt)[data[, t] == 0]
        }
        
        data_ctrl = data[data[, t] == 0, ]
        data_tr = data[data[, t] == 1, ]
        data_tr_dup = NULL
        for (j in 1:ncol(data_tr)) {
          data_tr_dup = cbind(data_tr_dup, rep(data_tr[, j], rep(2, nrow(data_tr))))
        }
        colnames(data_tr_dup) = colnames(data_ctrl)
        data_tr_dup = as.data.frame(data_tr_dup)
        data_tr_dup$rmpw[seq(2, length(data_tr_dup[, 1]), 2)] = 1
        d1 = c(rep(c(0, 1), nrow(data_tr)), rep(0, length(data_ctrl[, 1])))
        data_dup = cbind.data.frame(d1, rbind(data_tr_dup, data_ctrl))
        
        if(t.rand == TRUE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw)
        }
        if(t.rand == FALSE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
        }
        beta = coef(l.y)
        
        return(list(NDE = beta[2], NIE = beta[3]))
      }
      
      est.ori = est(data)
      NDE = est.ori$NDE
      NIE = est.ori$NIE
      
      NDE.boot = NULL
      NIE.boot = NULL
      for(i in 1:B){
        data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
        est.boot = est(data.boot)
        NDE.boot = c(NDE.boot, est.boot$NDE)
        NIE.boot = c(NIE.boot, est.boot$NIE)
      }

      SE_NDE = sd(NDE.boot)
      SE_NIE = sd(NIE.boot)
      
      P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
      P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
    }
    
    if(scale.m == "continuous" & scale.y == "binary"){
      est = function(data){
        if(!U){
          l = lm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), data = data)
          data$p1 = dnorm(data[, m], mean = as.matrix(cbind(1, 1, data[, X])) %*% coef(l), sd = sigma(l))
          data$p0 = dnorm(data[, m], mean = as.matrix(cbind(1, 0, data[, X])) %*% coef(l), sd = sigma(l))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X])) %*% coef(l.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X])) %*% coef(l.t)))
          }
        }
        
        if(U){
          l = lm(as.formula(paste(m, "~", t, "+", paste(X, collapse="+"))), offset = b.m * U, data = data)
          data$p1 = dnorm(data[, m], mean = as.matrix(cbind(1, 1, data[, X], U = data$U)) %*% c(coef(l), U = b.m), sd = sigma(l))
          data$p0 = dnorm(data[, m], mean = as.matrix(cbind(1, 0, data[, X], U = data$U)) %*% c(coef(l), U = b.m), sd = sigma(l))
          if(t.rand == FALSE){
            l.t = glm(as.formula(paste(t, "~", paste(X, collapse="+"))), offset = b.t * U, data = data, family = binomial)
            data$pt[data[, t] == 1] = 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 1, X], U = data[data[, t] == 1, "U"])) %*% c(coef(l.t), U = b.t)))
            data$pt[data[, t] == 0] = 1 - 1/(1 + exp(-as.matrix(cbind(1, data[data[, t] == 0, X], U = data[data[, t] == 0, "U"])) %*% c(coef(l.t), U = b.t)))
          }
        }
        
        data$rmpw[data[, t] == 1] = (data$p0/data$p1)[data[, t] == 1]
        data$rmpw[data[, t] == 0] = 1
        if(t.rand == FALSE){
          data$iptw[data[, t] == 1] = (mean(data[, t])/data$pt)[data[, t] == 1]
          data$iptw[data[, t] == 0] = ((1 - mean(data[, t]))/data$pt)[data[, t] == 0]
        }
        
        data_ctrl = data[data[, t] == 0, ]
        data_tr = data[data[, t] == 1, ]
        data_tr_dup = NULL
        for (j in 1:ncol(data_tr)) {
          data_tr_dup = cbind(data_tr_dup, rep(data_tr[, j], rep(2, nrow(data_tr))))
        }
        colnames(data_tr_dup) = colnames(data_ctrl)
        data_tr_dup = as.data.frame(data_tr_dup)
        data_tr_dup$rmpw[seq(2, length(data_tr_dup[, 1]), 2)] = 1
        d1 = c(rep(c(0, 1), nrow(data_tr)), rep(0, length(data_ctrl[, 1])))
        data_dup = cbind.data.frame(d1, rbind(data_tr_dup, data_ctrl))
        
        if(t.rand == TRUE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw)
        }
        if(t.rand == FALSE){
          if(!U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(X, collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
          if(U)
            l.y = lm(as.formula(paste(y, "~", t, "+ d1 +", paste(c(X, "U"), collapse = "+"))), data = data_dup, weights = rmpw * data_dup$iptw)
        }
        beta = coef(l.y)
        
        return(list(NDE = beta[2], NIE = beta[3]))
      }
      
      est.ori = est(data)
      NDE = est.ori$NDE
      NIE = est.ori$NIE
      
      NDE.boot = NULL
      NIE.boot = NULL
      for(i in 1:B){
        data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
        est.boot = est(data.boot)
        NDE.boot = c(NDE.boot, est.boot$NDE)
        NIE.boot = c(NIE.boot, est.boot$NIE)
      }
      
      SE_NDE = sd(NDE.boot)
      SE_NIE = sd(NIE.boot)
      
      P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
      P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
    }
    
    results = c(NIE = as.numeric(NIE), SE_NIE = SE_NIE, P_NIE = P_NIE, NDE = as.numeric(NDE), SE_NDE = SE_NDE, P_NDE = P_NDE)
    
    return(results)
  }
  
  # Regression-based method - VanderWeele's approach
  reg_vand = function(y, t, m, X, U = FALSE, scale.y, scale.m, b.m = NULL, b.y = NULL, B, data){
    if(scale.y == "continuous" & scale.m == "binary"){
      est = function(data){
        if(!U){
          l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), family = binomial(link = "logit"), data = data)
          l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
          data1 = data0 = data
          data1[, t] = 1
          data0[, t] = 0
          p1 = predict(l.m, data1, type="response")
          p0 = predict(l.m, data0, type="response")
          exp1 = p1/(1 - p1)
          exp0 = p0/(1 - p0)
        }
        if(U){
          l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, family = binomial(link = "logit"), data = data)
          l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, data = data)
          modelmatrix1 = modelmatrix0 = model.matrix(l.m)
          modelmatrix1[, t] = 1
          modelmatrix0[, t] = 0
          exp1 = exp(modelmatrix1 %*% coef(l.m) + b.m * data$U)
          exp0 = exp(modelmatrix0 %*% coef(l.m) + b.m * data$U)
          p1 = exp1/(1 + exp1)
          p0 = exp0/(1 + exp0) # This is exactly equal to predict(l.m, data0, type="response"). I do not use it here is because, when I use it, b.m cannot be found when I run reg().
        }
        beta0 = coef(l.m)[1]
        beta1 = coef(l.m)[2]
        theta0 = coef(l.y)[1]
        theta1 = coef(l.y)[2]
        theta2 = coef(l.y)[3]
        theta3 = coef(l.y)[length(coef(l.y))]
        
        NDE = theta1 + theta3 * mean(p0)
        NIE = (theta2 + theta3) * (mean(p1) - mean(p0))
        
        return(list(NDE = NDE, NIE = NIE))
      }
      
      est.ori = est(data)
      NDE = est.ori$NDE
      NIE = est.ori$NIE
      
      NDE.boot = NULL
      NIE.boot = NULL
      for(i in 1:B){
        data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
        est.boot = est(data.boot)
        NDE.boot = c(NDE.boot, est.boot$NDE)
        NIE.boot = c(NIE.boot, est.boot$NIE)
      }
      
      SE_NDE = sd(NDE.boot)
      SE_NIE = sd(NIE.boot)
      
      P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
      P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
    }
    
    if(scale.y == "continuous" & scale.m == "continuous"){
      est = function(data){
        if(!U){
          l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), data = data)
          l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
          beta2.c = model.matrix(l.m)[, -c(1,2)] %*% coef(l.m)[-c(1, 2)]
        }
        if(U){
          l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, data = data)
          l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, data = data)
          beta2.c = model.matrix(l.m)[, -c(1,2)] %*% coef(l.m)[-c(1, 2)] + b.m * data$U
        }
        beta0 = coef(l.m)[1]
        beta1 = coef(l.m)[2]
        theta0 = coef(l.y)[1]
        theta1 = coef(l.y)[2]
        theta2 = coef(l.y)[3]
        theta3 = coef(l.y)[length(coef(l.y))]
        
        NDE = theta1 + theta3 * beta0 + theta3 * mean(beta2.c)
        NIE = (theta2 + theta3) * beta1
        
        return(list(NDE = NDE, NIE = NIE))
      }
      
      est.ori = est(data)
      NDE = est.ori$NDE
      NIE = est.ori$NIE
      
      NDE.boot = NULL
      NIE.boot = NULL
      for(i in 1:B){
        data.boot = data[sample(1:nrow(data), nrow(data), replace = T), ]
        est.boot = est(data.boot)
        NDE.boot = c(NDE.boot, est.boot$NDE)
        NIE.boot = c(NIE.boot, est.boot$NIE)
      }
      
      SE_NDE = sd(NDE.boot)
      SE_NIE = sd(NIE.boot)
      
      P_NDE = (1 - pnorm(abs(NDE/SE_NDE))) * 2
      P_NIE = (1 - pnorm(abs(NIE/SE_NIE))) * 2
    }
    
    results = c(NIE = as.numeric(NIE), SE_NIE = SE_NIE, P_NIE = P_NIE, NDE = as.numeric(NDE), SE_NDE = SE_NDE, P_NDE = P_NDE)
    
    return(results)
  }
  
  # Regression-based method - Imai's approach
  reg_imai = function(y, t, m, X, U = FALSE, scale.y, scale.m, b.m = NULL, b.y = NULL, data){
    if(scale.y == "continuous" & scale.m == "binary"){
      if(!U){
        l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), family = binomial(link = "logit"), data = data)
        l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
      }
      if(U){
        l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, family = binomial(link = "logit"), data = data)
        l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, data = data)
      }
      results = mediate(l.m, l.y, treat = t, mediator = m, robustSE = FALSE)
      NIE = results$d1
      SE_NIE = sd(results$d1.sim)
      P_NIE = results$d1.p
      NDE = results$z0
      SE_NDE = sd(results$z0.sim)
      P_NDE = results$z0.p
    }
    
    if(scale.y == "continuous" & scale.m == "continuous"){
      if(!U){
        l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), data = data)
        l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
      }
      if(U){
        l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, data = data)
        l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, data = data)
      }
      results = mediate(l.m, l.y, treat = t, mediator = m, robustSE = FALSE)
      NIE = results$d1
      SE_NIE = sd(results$d1.sim)
      P_NIE = results$d1.p
      NDE = results$z0
      SE_NDE = sd(results$z0.sim)
      P_NDE = results$z0.p
    }
    
    if(scale.y == "binary" & scale.m == "binary"){
      if(!U){
        l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), family = binomial(link = "logit"), data = data)
        l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), data = data)
      }
      if(U){
        l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, family = binomial(link = "logit"), data = data)
        l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), offset = b.y * U, data = data)
      }
      results = mediate(l.m, l.y, treat = t, mediator = m, robustSE = FALSE)
      NIE = results$d1
      SE_NIE = sd(results$d1.sim)
      P_NIE = results$d1.p
      NDE = results$z0
      SE_NDE = sd(results$z0.sim)
      P_NDE = results$z0.p
    }
    
    if(scale.y == "binary" & scale.m == "continuous"){
      if(!U){
        l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), data = data)
        l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), data = data)
      }
      if(U){
        l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), offset = b.m * U, data = data)
        l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), offset = b.y * U, data = data)
      }
      results = mediate(l.m, l.y, treat = t, mediator = m, robustSE = FALSE)
      NIE = results$d1
      SE_NIE = sd(results$d1.sim)
      P_NIE = results$d1.p
      NDE = results$z0
      SE_NDE = sd(results$z0.sim)
      P_NDE = results$z0.p
    }
    
    results = c(NIE = NIE, SE_NIE = SE_NIE, P_NIE = P_NIE, NDE = NDE, SE_NDE = SE_NDE, P_NDE = P_NDE)
    
    return(results)
  }
  
  # Generation of U based on the stochastic EM algorithm
  genU = function(y, t, t.rand = TRUE, m, X, scale.y, scale.m, scale.U, b.y, b.m, b.t, p.u, sigma.u, iter = 10, data) {
    # Generate U from its prior distribution first
    if(scale.U == "binary")
      data$U = rbinom(nrow(data), 1, p.u) 
    if(scale.U == "continuous")
      data$U = rnorm(nrow(data), 0, sigma.u) 
    coef.y.updated = NULL
    coef.m.updated = NULL
    if(t.rand == FALSE)
      coef.t.updated = NULL
    for(i in 1:iter) {
      if(t.rand == FALSE){
        l.t = glm(as.formula(paste(t, "~", paste(X, collapse = "+"))), offset = b.t * U, family = binomial(link = "logit"), data = data)
        coef.t = c(l.t$coef, U = b.t)
      }
      if(scale.y == "continuous"){
        l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, data = data)
        coef.y = c(l.y$coef, U = b.y)
        sd.y = sigma(l.y) # This is the sd of Y conditional on t, m, X, and U.
      }
      if(scale.y == "binary"){
        l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), offset = b.y * U, family = binomial(link = "logit"), data = data)
        coef.y = c(l.y$coef, U = b.y)
      }
      if(scale.m == "continuous"){
        l.m = lm(as.formula(paste(m, "~", t, "+", paste(X, collapse = "+"))), offset = b.m * U, data = data)
        coef.m = c(l.m$coef, U = b.m)
        sd.m = sigma(l.m) # This is the sd of M conditional on t, m, X, and U.
      }
      if(scale.m == "binary"){
        l.m = glm(as.formula(paste(m, "~", t, "+", paste(X, collapse = "+"))), offset = b.m * U, family = binomial(link = "logit"), data = data)
        coef.m = c(l.m$coef, U = b.m)
      }
      
      if(t.rand == FALSE){
        # Conditional probability of T
        pt1u1 = 1/(1 + exp(-(cbind(model.matrix(l.t), U = 1) %*% coef.t)))
        pt1u0 = 1/(1 + exp(-(cbind(model.matrix(l.t), U = 0) %*% coef.t)))
        ptu1 = pt1u1^data[, t] * (1 - pt1u1)^(1 - data[, t])
        ptu0 = pt1u0^data[, t] * (1 - pt1u0)^(1 - data[, t])
      }
      
      # For a binary U
      if(scale.U == "binary"){
        # Conditional probability of Y
        # If Y is continuous
        if(scale.y == "continuous"){
          mean.yu1 = cbind(model.matrix(l.y), U = 1) %*% coef.y
          sd.yu1 = sd.y
          mean.yu0 = cbind(model.matrix(l.y), U = 0) %*% coef.y
          sd.yu0 = sd.y
          pyu1 = dnorm(data[, y], mean = mean.yu1, sd = sd.yu1)
          pyu0 = dnorm(data[, y], mean = mean.yu0, sd = sd.yu0)
        }
        # If Y is binary
        if(scale.y == "binary"){
          py1u1 = 1/(1 + exp(-(cbind(model.matrix(l.y), U = 1) %*% coef.y)))
          py1u0 = 1/(1 + exp(-(cbind(model.matrix(l.y), U = 0) %*% coef.y)))
          pyu1 = py1u1^data[, y] * (1 - py1u1)^(1 - data[, y])
          pyu0 = py1u0^data[, y] * (1 - py1u0)^(1 - data[, y])
        }
        
        # Conditional probability of M
        # If M is continuous
        if(scale.m == "continuous"){
          mean.mu1 = cbind(model.matrix(l.m), U = 1) %*% coef.m
          sd.mu1 = sd.m
          mean.mu0 = cbind(model.matrix(l.m), U = 0) %*% coef.m
          sd.mu0 = sd.m
          pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
          pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
        }
        # If M is binary
        if(scale.m == "binary"){
          pm1u1 = 1/(1 + exp(-(cbind(model.matrix(l.m), U = 1) %*% coef.m)))
          pm1u0 = 1/(1 + exp(-(cbind(model.matrix(l.m), U = 0) %*% coef.m)))
          pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
          pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
        }
        
        # Conditional probability of U  
        if(t.rand == TRUE)
          p = pyu1 * pmu1 * p.u/(pyu1 * pmu1 * p.u + pyu0 * pmu0 * (1 - p.u))
        if(t.rand == FALSE)
          p = pyu1 * pmu1 * p.u * ptu1/(pyu1 * pmu1 * p.u * ptu1 + pyu0 * pmu0 * (1 - p.u) * ptu0)
        U = NULL
        for (k in 1:nrow(data)) U = c(U, rbinom(1, 1, p[k]))
        data$U = U
      }
      
      # For a continuous U
      if(scale.U == "continuous"){
        # If Y is continuous and M is binary
        if(scale.y == "continuous" & scale.m == "binary"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              integrand = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          
          # Generate random values of U
          U = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of U
            if(t.rand == TRUE){
              conditional.u = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              conditional.u = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            U = c(U, rdist(1))
          } 
          data$U = U
        }
        
        if(scale.y == "continuous" & scale.m == "continuous"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              integrand = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of U
          U = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of U
            if(t.rand == TRUE){
              conditional.u = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              conditional.u = function(u){
                dnorm(data[i, y], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            U = c(U, rdist(1))
          } 
          data$U = U       
        }
        
        if(scale.y == "binary" & scale.m == "binary"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              integrand = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of U
          U = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of U
            if(t.rand == TRUE){
              conditional.u = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              conditional.u = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^data[i, m] * (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - data[i, m]) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            U = c(U, rdist(1))
          } 
          data$U = U
        }
        
        if(scale.y == "binary" & scale.m == "continuous"){
          # Calculate the integral in the denominator
          integral = NULL
          for(i in 1:nrow(data)){
            if(t.rand == TRUE){
              integrand = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            if(t.rand == FALSE){
              integrand = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)
              }
            }
            integral = c(integral, integrate(integrand, lower = -10, upper = 10)$value)
          }
          # Generate random values of U
          U = NULL
          for (i in 1:nrow(data)){
            # Obtain the condition probability of U
            if(t.rand == TRUE){
              conditional.u = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            if(t.rand == FALSE){
              conditional.u = function(u){
                (1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^data[i, y] * (1 - 1/(1 + exp(-(c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y))))^(1 - data[i, y]) *
                  dnorm(data[i, m], mean = c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m, sd = sd.m) *
                  (1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^data[i, t] * (1 - 1/(1 + exp(-(c(model.matrix(l.t)[i, ] %*% l.t$coef) + u * b.t))))^(1 - data[i, t]) *
                  dnorm(u, mean = 0, sd = sigma.u)/integral[i]
              }
            }
            dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
            rdist = r(dist)# function to create random variates from conditional.u
            U = c(U, rdist(1))
          } 
          data$U = U
        }
      }
      
      # Check the convergence of the coefficients
      coef.y.updated = rbind(coef.y.updated, coef.y)
      coef.m.updated = rbind(coef.m.updated, coef.m)
    }
    
    return(list(U = U, coef.y.updated = coef.y.updated, coef.m.updated = coef.m.updated))
  }
  
  if(est == "rmpw")
    results.old = rmpw(y = y, t = t, t.rand = t.rand, m = m, X = X, scale.y = scale.y, scale.m = scale.m, B = B, data = data)
  if(est == "regression")
    results.old = reg_vand(y = y, t = t, m = m, X = X, scale.y = scale.y, scale.m = scale.m, B = B, data = data)
  if(est == "simulation")
    results.old = reg_imai(y = y, t = t, m = m, X = X, scale.y = scale.y, scale.m = scale.m, data = data)
  
  if(scale.y == "continuous" & scale.m == "binary"){
    l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), family = binomial(link = "logit"), data = data)
    l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
  }
  
  if(scale.y == "continuous" & scale.m == "continuous"){
    l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), data = data)
    l.y = lm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), data = data)
  }
  
  if(scale.y == "binary" & scale.m == "binary"){
    l.m = glm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), family = binomial(link = "logit"), data = data)
    l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), data = data)
  }
  
  if(scale.y == "binary" & scale.m == "continuous"){
    l.m = lm(as.formula(paste(m, "~", paste(c(t, X), collapse = "+"))), data = data)
    l.y = glm(as.formula(paste(y, "~", t, "*", m, "+", paste(X, collapse = "+"))), family = binomial(link = "logit"), data = data)
  }
  
  X.coef.plot = cbind(l.m$coef[3:(2 + length(X))], l.y$coef[4:(3 + length(X))])
  b.m.all = seq(range.b.m[1], range.b.m[2], length.out = grid.b.m)
  b.y.all = seq(range.b.y[1], range.b.y[2], length.out = grid.b.y)
  
  vals = matrix(NA, grid.b.m, grid.b.y)
  rownames(vals) = round(b.m.all, 2)
  colnames(vals) = round(b.y.all, 2)
  NDE.all = NIE.all = SE_NDE.all = SE_NIE.all = P_NDE.all = P_NIE.all = vals
  cell = 0
  for(i in 1:grid.b.m){
    for(j in 1:grid.b.y){
      cell = cell + 1
      results = NULL
      for(k in 1:nsim){
        if(t.rand == TRUE){
          data$U = genU(y = y, t = t, m = m, X = X, scale.y = scale.y, scale.m = scale.m, scale.U = scale.U, b.y = b.y.all[j], b.m = b.m.all[i], p.u = p.u, sigma.u = sigma.u, iter = iter, data = data)$U
          if(est == "rmpw")
            results = rbind(results, rmpw(y = y, t = t, m = m, X = X, U = TRUE, scale.y = scale.y, scale.m = scale.m, b.m = b.m.all[i], B = B, data = data))
          if(est == "regression")
            results = rbind(results, reg_vand(y = y, t = t, m = m, X = X, U = TRUE, scale.y = scale.y, scale.m = scale.m, b.y = b.y.all[j], b.m = b.m.all[i], B = B, data = data))
          if(est == "simulation")
            results = rbind(results, reg_imai(y = y, t = t, m = m, X = X, U = TRUE, scale.y = scale.y, scale.m = scale.m, b.y = b.y.all[j], b.m = b.m.all[i], data = data))
        }
        if(t.rand == FALSE){
          data$U = genU(y = y, t = t, t.rand == FALSE, m = m, X = X, scale.y = scale.y, scale.m = scale.m, scale.U = scale.U, b.y = b.y.all[j], b.m = b.m.all[i], b.t = b.t, p.u = p.u, sigma.u = sigma.u, iter = iter, data = data)$U
          if(est == "rmpw")
            results = rbind(results, rmpw(y = y, t = t, t.rand = FALSE, m = m, X = X, U = TRUE, scale.y = scale.y, scale.m = scale.m, b.m = b.m.all[i], b.t = b.t, B = B, data = data))
          if(est == "regression"|est == "simulation")
            stop("The estimation method is not applicable when the treatment is not randomized")
        }
      }
      NDE.all[i, j] = mean(results[, "NDE"])
      NIE.all[i, j] = mean(results[, "NIE"])
      SE_NDE.all[i, j] = sqrt(mean(results[, "SE_NDE"]^2) + (1 + 1/nsim) * var(results[, "NDE"]))
      SE_NIE.all[i, j] = sqrt(mean(results[, "SE_NIE"]^2) + (1 + 1/nsim) * var(results[, "NIE"]))
      P_NDE.all[i, j] = (1 - pnorm(abs(NDE.all[i, j]/SE_NDE.all[i, j]))) * 2
      P_NIE.all[i, j] = (1 - pnorm(abs(NIE.all[i, j]/SE_NIE.all[i, j]))) * 2
      message("Completed ", cell, " of ", grid.b.m * grid.b.y, " cells.\n")
    }
  }
  results.new = list(NDE = NDE.all, NIE = NIE.all, SE_NDE = SE_NDE.all, SE_NIE = SE_NIE.all, P_NDE = P_NDE.all, P_NIE = P_NIE.all)
  
  if(t.rand == TRUE)
    return(list(results.old = results.old, X.coef.plot = X.coef.plot, range.b.m = range.b.m, range.b.y = range.b.y, b.y.all = b.y.all, b.m.all = b.m.all, results.new = results.new))
  if(t.rand == FALSE)
    return(list(results.old = results.old, b.t = b.t, range.b.m = range.b.m, range.b.y = range.b.y, b.y.all = b.y.all, b.m.all = b.m.all, results.new = results.new))
}

#' Simulation-Based Sensitivity Analysis Plot for Causal Mediation Analysis
#' 
#' This function is used to visually represent the sensitivity analysis results when the treatment is randomized. Each black contour represents the combinations of sensitivity parameters that lead to the same indirect effect estimate as indicated by the number on the contour. The sensitivity parameters along the red dashed curves reduce the estimate to zero. In the region between the two blue dotted curves that bracket the red curve, the effect is changed to be insignificant at the significance level of 0.05. The larger the magnitudes of the sensitivity parameters are for removing the effects or changing their significance, the less sensitive the results are. Each dot corresponds to the conditional associations of each observed covariate with the outcome and the mediator, which are used to calibrate the strength of the sensitivity parameters. 
#' 
#' @param sens.results An output from the sens function.
#' @param effect The name of the effect whose sensitivity analysis results are to be plotted (string). "NIE" for natural indirect effect. "NDE" for natural direct effect. The default is "NIE".
#' @param est The method used in the estimation of the causal mediation effects. "rmpw" uses the weighting-based method. "regression" uses the regression-based method developed by VanderWeele and Vansteelandt (2009). "simulation" uses the simulation-based method developed by Imai et al. (2010a, 2010b).
#' @return Sensitivity analysis plots for the natural direct and in direct effects.
#' @author Xu Qin and Fan Yang
#' @references Qin, X., & Yang, F. (2020). Simulation-Based Sensitivity Analysis for Causal Mediation Studies.
#' @export
#' @importFrom graphics plot contour abline

sens.plot = function(sens.results, effect = "NIE", est){
  plot(sens.results$X.coef.plot[, 1], sens.results$X.coef.plot[, 2], pch = 19, xlim = sens.results$range.b.m, ylim = sens.results$range.b.y, xlab = "Partial Effect of U on logit of Pr(M=1)", ylab="Partial Effect of U on Y", main = paste("Sensitivity Analysis for", effect, "(", est, ")"))
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)
  contour(sens.results$b.m.all, sens.results$b.y.all, sens.results$results.new[[effect]], add = T, labcex = 1)
  # Add curves that represent the contour along which the treatment effect estimate is reduced to zero
  contour(sens.results$b.m.all, sens.results$b.y.all, sens.results$results.new[[effect]], levels = 0, add = T, col = "red", lwd = 2, lty = 2, labcex = 1)
  # Add curves for significance change
  contour(sens.results$b.m.all, sens.results$b.y.all, sens.results$results.new[[effect]]/sens.results$results.new[[paste0("SE_", effect)]], levels = c(-1, 1) * qnorm(0.05/2), labels = "Sig.Change", add = T, col = "blue", lwd = 2, lty = 3, labcex = 1)
}