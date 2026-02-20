# Dose-response function on log(Response) scale ----
DR_logRN_fun = function(Dose, Rmin, Rmax, beta, NEC){
  logyhat = Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)) 
  yhat = exp(logyhat)
  return(yhat)
}

DR_fun_log = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - exp(beta) * (Dose - NEC))
  return(log_yhat)
}

DR_fun_log_exp = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - exp(beta) * (Dose - NEC))
  yhat = exp(log_yhat)
  return(yhat)
}

