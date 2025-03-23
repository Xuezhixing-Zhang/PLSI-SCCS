#' Function for SCCS fit.
#' 
#' 
#' @param case A vector of case indices.
#' @param aevent A vector of event times.
#' @param adrug A vector or a matrix of ages at start of risk periods. One column per risk period. NA denotes no risk for this individual. Please make sure all individuals have the same length of risk period.
#' @param aedrug A vector or a matrix of ages at end of risk periods. One column per risk period. NA denotes no risk for this individual. Please make sure all individuals have the same length of risk period.
#' @param agebreaks A vector determines age intervals for each age group. This vector is rightmost and left open.
#' @param gamma Smoothing parameter for S1 curve.
#' @param M_1 Number of basis functions for generating B_1 spline. Default is 3. Can't be less than 3.
#' @param M_2 A vector of numbers of basis functions for B_2 splines.
#' @param same_expo True or False. If True, we assume all risk periods have the same curve shape.
#' @param warmstart A vector of lambda used as warm starts. Default is NA. 
#' @return A list of model estimates.
#' 
#' @keywords internal
#' @seealso \code{\link{SCCSFit}}

SCCSFit_nopar <- function(astart,
                          aend, 
                          case, 
                          aevent,
                          adrug, 
                          aedrug,
                          agebreaks,
                          M_1,
                          M_2,
                          gamma,
                          same_expo = FALSE,
                          ninitials = 1,
                          warmstart = NA){
  
  # Generate Observaiton list
  obs_list <- mapply(seq, from = astart, to = aend, SIMPLIFY = FALSE)
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  # Generate Age list
  age_labels <- paste(agebreaks[-length(agebreaks)], agebreaks[-1], sep = "-")
  age_labels[length(agebreaks)] <- paste0(agebreaks[length(agebreaks)],"+")
  age_list <- lapply(obs_list, function(x){
    return(findInterval(x, agebreaks, left.open = T, rightmost.closed = F))
  })
  
  # Optimization for all parameters
  if(same_expo){lambda_length <- M_2[1]}
  else{lambda_length <- sum(M_2)}
  p <- dim(cbind(adrug))[2]
  
  if(any(is.na(warmstart)) == T){
    
    para_0s <- matrix(runif(ninitials*lambda_length,-1,1), ncol = lambda_length)
    paras <- matrix(0, ncol = lambda_length, nrow = ninitials)
    values <-  NULL
    models <- list()
    
    for(i in 1:nrow(para_0s)){
      para_0 <- c(para_0s[i,])
      models[[i]] <- optim_all(warmstart = para_0,p,
                               astart,adrug,aedrug,aevent,
                               obs_list,age_list,
                               same_expo,
                               M_1,M_2,gamma)
      values[i] <- models[[i]]$obj
    }
    model <- models[[which.min(values)]]
  }
  else{
    model <- optim_all(warmstart = warmstart,p,
                       astart,adrug,aedrug,aevent,
                       obs_list,age_list,
                       same_expo,
                       M_1,M_2,gamma)
  }
  
  para_optimal <- model$para_optimal
  para_optimal_norm <- model$para_optimal_norm
  para_beta <- model$beta
  para_theta <- model$theta
  lambda_temp <- model$lambda_temp
  
  return(list(lambda = para_optimal, 
              lambda_normalized = para_optimal_norm,
              lambda_temp = lambda_temp,
              beta = para_beta, 
              theta = para_theta,
              `Boundary.knots` = model$`Boundary.knots`,
              knots = model$knots,
              knots_B2 = model$knots_B2,
              gamma = gamma,M_1 = M_1, M_2 = M_2,
              obj = model$obj,
              penalty = model$penalty))
  
}


#######################################################
#######################################################
#######################################################
#######################################################
#' An internal optimization function for lambda.
#' 
#' @keywords internal
#' @seealso \code{\link{SCCSFit_nopar}}
optim_all <- function(warmstart, p,astart,adrug,aedrug,aevent,obs_list,age_list,same_expo, M_1,M_2,gamma){
  
  ### Obtain optimal lambda
  if(any(is.na(warmstart)) == T){
    para_0 <- runif(lambda_length,-1,1)
  }
  else{para_0 <- warmstart}
  
  
  storage <- new.env()
  storage$lambda <- para_0
  storage$warm_para_beta_theta <- rep(0.5, M_1+length(agebreaks))
  # storage$warm_para_beta_theta <- c(rep(0.5,M_1),0.4,0.8,1.2)
  result <- optim(par=para_0,
                  fn=optim_lambda,
                  store = storage,
                  method = "BFGS",
                  gr=NULL,
                  p = p,
                  astart = astart,
                  adrug = adrug,
                  aedrug = aedrug,
                  aevent = aevent,
                  same_expo =same_expo,
                  M_1 = M_1,
                  M_2 = M_2,
                  gamma = gamma,
                  obs_list = obs_list,
                  age_list = age_list,
                  control=list(trace = 10,REPORT = 10,maxit=10000))
  value <- result$value
  paras <- result$par
  
  para_optimal <- paras
  if(same_expo){para_optimal <- rep(para_optimal,length(M_2))}
  para_optimal_norm <- para_optimal*sign(para_optimal[1]) 
  para_optimal_norm <- para_optimal_norm/sqrt(sum(para_optimal_norm^2))
  
  paras_beta_theta <- storage$warm_para_beta_theta
  para_beta <- paras_beta_theta[1:M_1]
  para_theta <- paras_beta_theta[-c(1:M_1)]
  
  knots_B2 <- storage$knots_B2
  Boundary.knots <- storage$Boundary.knots
  internal_knots <- storage$internal_knots
  
  ### Construct model results
  model <- list(para_optimal = para_optimal,
                para_optimal_norm = para_optimal_norm,
                lambda_temp = storage$lambda,
                beta = para_beta, 
                theta = para_theta,
                `Boundary.knots` = `Boundary.knots`,
                knots = internal_knots,
                knots_B2 = knots_B2,
                obj = -value,
                penalty = storage$penalty)
  
  return(model)
  
}


#' An internal optimization function for lambda.
#' 
#' @keywords internal
#' @seealso \code{\link{SCCSFit_nopar}}

optim_lambda <- function(lambda_coefs, store, p,astart,adrug,aedrug,aevent,obs_list,age_list,same_expo,M_1,M_2,gamma){
  
  time0 <- Sys.time()
  store$lambda <- lambda_coefs
  if(same_expo == T){
    lambda_normalized <- rep(lambda_coefs,length(M_2))
  }
  else{
    lambda_normalized <- lambda_coefs
  }
  lambda_normalized <- lambda_normalized*sign(lambda_normalized[1]) 
  lambda_normalized <- lambda_normalized/sqrt(sum(lambda_normalized^2))
  # print(lambda_normalized)
  ### Obtain internal knots for S2
  relative_event_time <- -sweep(adrug, 1, aevent,"-")
  knots_B2 <- list()
  u_values <- list()
  risk_times <- list()
  risk_length <- sapply(1:p, function(i) mean(aedrug[,i] - adrug[,i], na.rm = T))
  M_2_indices <- c(1,cumsum(M_2)+1)
  min_u_values <- NULL
  max_u_values <- NULL
  
  if(same_expo == FALSE){
    for(i in 1:p){
      relative_times <- c(relative_event_time[,i])
      relative_times <- relative_times[relative_times <= risk_length[i] & relative_times >= 0]
      relative_times <- relative_times[!is.na(relative_times)]
      knots_B2[[i]] <- quantile(relative_times, probs = seq(0, 1, length.out = M_2[i])[-c(1,M_2[i])])
      risk_times[[i]] <- seq(0,risk_length[i],by=1)
      u_values[[i]] <- ns(risk_times[[i]], 
                          Boundary.knots = range(0, risk_length[i]), 
                          knots = knots_B2[[i]],intercept = TRUE) %*% lambda_normalized[M_2_indices[i]:(M_2_indices[i+1]-1)]
      min_u_values[i] <- min(u_values[[i]])
      max_u_values[i] <- max(u_values[[i]])
    }
  }
  else{
    
    relative_times_all <- NULL
    for(i in 1:p){
      relative_times <- c(relative_event_time[,i])
      relative_times <- relative_times[relative_times <= risk_length[i] & relative_times >= 0]
      relative_times <- relative_times[!is.na(relative_times)]
      relative_times_all <- c(relative_times_all, relative_times)
    }
    knots <- quantile(relative_times_all, probs = seq(0, 1, length.out = M_2[i])[-c(1,M_2[i])])
    for(i in 1:p){
    knots_B2[[i]] <- quantile(relative_times_all, probs = seq(0, 1, length.out = M_2[i])[-c(1,M_2[i])])
    risk_times[[i]] <- seq(0,risk_length[i],by=1)
    u_values[[i]] <- ns(risk_times[[i]], 
                        Boundary.knots = range(0, risk_length[i]), 
                        knots = knots_B2[[i]],intercept = TRUE) %*% lambda_normalized[M_2_indices[i]:(M_2_indices[i+1]-1)]
    min_u_values[i] <- min(u_values[[i]])
    max_u_values[i] <- max(u_values[[i]])
    }
  }
 
  store$knots_B2 <- knots_B2
  # print(knots_B2)
  ### Obtain internal knots for S2
  u_value_list <- list()
  for(i in 1:p){
    u_value_list[[i]] <- mapply(function(obs, adrug) {
      rel_time <- obs - adrug
      rel_time <- ifelse(rel_time <0 | rel_time > risk_length[i], Inf, rel_time)
      u_value_vec <- u_values[[i]][match(rel_time, risk_times[[i]])]
      u_value_vec[is.na(u_value_vec)] <- 0
      return(u_value_vec)
    }, 
    obs = obs_list, 
    adrug = adrug[,i], 
    SIMPLIFY = FALSE)
  }
  # print(length(obs_list))
  # print(dim(adrug))
  sum_list <- Reduce(function(x, y) Map(`+`, x, y), u_value_list)
  combined_u <- unlist(sum_list)
  combined_u <- combined_u[combined_u != 0]
  unique_u <- unique(combined_u)
  
  ### Obtain internal knots and boundary knots for u values.
  Boundary.knots <- c(min(c(unique_u,
                            min_u_values,
                            sum(min_u_values*I(min_u_values<0)))), 
                      max(unique_u,
                          max_u_values,
                          sum(max_u_values*I(max_u_values>0))))
  
  internal_knots <- quantile(combined_u, probs = seq(0, 1, length.out = M_1-2)[-c(1,M_1-2)])
  store$Boundary.knots <- Boundary.knots
  store$internal_knots <- internal_knots
  ### Obtain B spline %*% beta for unique u values
  
  Bspline_obj <- bs(unique_u, 
                    knots = internal_knots, 
                    Boundary.knots = Boundary.knots,
                    degree = 3, intercept = T)
  ## Obtain B spline 
  all_knots <- c(Boundary.knots[1], internal_knots, Boundary.knots[2])
  nbasis <- length(internal_knots) + 4
  bspline_basis <- create.bspline.basis(rangeval = Boundary.knots,
                                        nbasis = nbasis,
                                        norder = 4,
                                        breaks = all_knots)
  # penalty_matrix_1 <- bsplinepen(bspline_basis, Lfdobj = 1)
  penalty_matrix_2 <- bsplinepen(bspline_basis, Lfdobj = 2)
  
  ### Obtain optimal beta and theta given lambda:
  
  para_beta_theta <- store$warm_para_beta_theta   #round(store$warm_para_beta_theta,2) 
  # para_beta_theta <- rep(0.5, length(store$warm_para_beta_theta))
  result <- optim(par= para_beta_theta,
                  fn=optim_beta,
                  method = "BFGS",
                  gr=NULL,
                  astart = astart,
                  aevent = aevent,
                  Bspline_obj = Bspline_obj,
                  unique_u = unique_u,
                  sum_list = sum_list,
                  age_list = age_list,
                  M_1= M_1,
                  gamma = gamma,
                  penalty_matrix_2 = penalty_matrix_2,
                  control=list(maxit=10000))
  
  value <- result$value
  para_beta_theta <- result$par
  store$warm_para_beta_theta <- para_beta_theta

  # print(store$warm_para_beta_theta)
  time2 <- Sys.time()
  
  beta_coefs <- para_beta_theta[c(1:M_1)]
  store$penalty <- gamma*(t(beta_coefs^2) %*% penalty_matrix_2 %*% beta_coefs^2)
  
  print(difftime(time2, time0, units = "secs"))
  return(value)
}

#' An internal optimization function for beta.
#' 
#' @keywords internal
#' @seealso \code{\link{SCCSFit_nopar}}
optim_beta <- function(para_beta_theta, astart,aevent,
                       Bspline_obj, unique_u,sum_list,age_list, penalty_matrix_2,M_1,gamma=gamma){
  
  beta_coefs <- para_beta_theta[1:M_1]
  theta_coefs <- para_beta_theta[-c(1:M_1)]
  
  unique_rr <- Bspline_obj %*% (beta_coefs^2)
  unique_u <- c(0, unique_u)
  unique_rr <- c(1, unique_rr)
  
  ### Map them into u_value_list
  rr_list <- lapply(sum_list, function(x){unique_rr[match(x, unique_u)]})
  
  ### Age risks
  theta_coefs_all <- c(0, theta_coefs)
  rr_age_list <- lapply(age_list, function(x){exp(theta_coefs_all[x+1])})
  
  ### Final relative risks
  RR_list <- Map(`*`, rr_list, rr_age_list)
  
  ### Obtain risks at event
  aevent_indices <- aevent - astart + 1
  event_log_risk <- log(mapply(function(RR, index) RR[index], RR_list, aevent_indices))
  
  ### Obtain risks over the whole observational period
  sum_log_risk <- sapply(RR_list, function(x)(log(sum(x))))
  log_lik <- sum(event_log_risk) - sum(sum_log_risk)
  
  obj <- log_lik - gamma*(t(beta_coefs^2) %*% penalty_matrix_2 %*% beta_coefs^2)
  
  return(-obj)
}


#######################################################
#######################################################
#######################################################
#######################################################
predict_SCCSFit <- function(astart,
                            aend, 
                            case, 
                            aevent,
                            adrug, 
                            aedrug,
                            agebreaks,
                            model){
  
  ##### Extract parameters
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  lambda_normalized <- model$lambda_normalized
  knots_B2 <- model$knots_B2
  M_2 <- model$M_2
  M_1 <- model$M_1
  Boundary.knots <- model$`Boundary.knots`
  internal_knots <- model$knots
  theta_coefs <- model$theta
  beta_coefs <- model$beta
  
  #####
  obs_list <- mapply(seq, from = astart, to = aend, SIMPLIFY = FALSE)
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  age_labels <- paste(agebreaks[-length(agebreaks)], agebreaks[-1], sep = "-")
  age_labels[length(agebreaks)] <- paste0(agebreaks[length(agebreaks)],"+")
  age_list <- lapply(obs_list, function(x){
    return(findInterval(x, agebreaks, left.open = T, rightmost.closed = F))
  })
  
  p <- dim(cbind(adrug))[2]
  
  ### Obtain internal knots for S2
  relative_event_time <- -sweep(adrug, 1, aevent,"-")
  u_values <- list()
  risk_times <- list()
  risk_length <- sapply(1:p, function(i) mean(aedrug[,i] - adrug[,i], na.rm = T))
  M_2_indices <- c(1,cumsum(M_2)+1)
  for(i in 1:p){
    relative_times <- c(relative_event_time[,i])
    relative_times <- relative_times[relative_times <= risk_length[i] & relative_times >= 0]
    relative_times <- relative_times[!is.na(relative_times)]
    risk_times[[i]] <- seq(0,risk_length[i],by=1)
    u_values[[i]] <- ns(risk_times[[i]], 
                        Boundary.knots = range(0, risk_length[i]), 
                        knots = knots_B2[[i]],intercept = TRUE) %*% lambda_normalized[M_2_indices[p]:(M_2_indices[p+1]-1)]
  }
  
  u_value_list <- list()
  for(i in 1:p){
    u_value_list[[i]] <- mapply(function(obs, adrug) {
      rel_time <- obs - adrug
      rel_time <- ifelse(rel_time <0 | rel_time > risk_length[i], Inf, rel_time)
      u_value_vec <- u_values[[i]][match(rel_time, risk_times[[i]])]
      u_value_vec[is.na(u_value_vec)] <- 0
      return(u_value_vec)
    }, 
    obs = obs_list, 
    adrug = adrug[,i], 
    SIMPLIFY = FALSE)
  }
  
  sum_list <- Reduce(function(x, y) Map(`+`, x, y), u_value_list)
  combined_u <- unlist(sum_list)
  combined_u <- combined_u[combined_u != 0]
  
  ### Obtain B spline %*% beta for unique u values
  unique_u <- unique(combined_u)
  Bspline_obj <- bs(unique_u, 
                    knots = internal_knots, 
                    Boundary.knots = Boundary.knots,
                    degree = 3, intercept = T)
  
  ## Obtain B spline 
  all_knots <- c(Boundary.knots[1], internal_knots, Boundary.knots[2])
  nbasis <- length(internal_knots) + 4
  bspline_basis <- create.bspline.basis(rangeval = Boundary.knots,
                                        nbasis = nbasis,
                                        norder = 4,
                                        breaks = all_knots)
  penalty_matrix_2 <- bsplinepen(bspline_basis, Lfdobj = 2)
  
  unique_rr <- Bspline_obj %*% (beta_coefs^2)
  unique_u <- c(0, unique_u)
  unique_rr <- c(1, unique_rr)
  unique_rr[I(unique_rr<0)] <- 0
  ### Map them into u_value_list
  rr_list <- lapply(sum_list, function(x){unique_rr[match(x, unique_u)]})
  
  ### Age risks
  theta_coefs_all <- c(0, theta_coefs)
  rr_age_list <- lapply(age_list, function(x){exp(theta_coefs_all[x+1])})
  
  ### Final relative risks
  RR_list <- Map(`*`, rr_list, rr_age_list)
  
  ### Obtain risks at event
  aevent_indices <- aevent - astart + 1
  event_risk <- mapply(function(RR, index) RR[index], RR_list, aevent_indices)
  
  integral_risk <- sapply(RR_list, function(x)(sum(x)))
  conditional_probs <- mapply(function(RR, sum_risk) RR/sum_risk, RR_list, integral_risk)
  conditional_event_probs <- event_risk/integral_risk
  
  ### Obtain CV score
  cv_score <- sum(1-conditional_event_probs) + sum(unlist(conditional_probs)) - sum(conditional_event_probs)
  return(cv_score)
}



SCCS_CV <- function(astart,
                    aend, 
                    case, 
                    aevent,
                    adrug, 
                    aedrug,
                    agebreaks,
                    M_1,
                    M_2,
                    gamma,
                    same_expo = FALSE,
                    nfolds = 5,
                    no_cores = 10){
  
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  # Divide event data into n folds.
  indices <- unique(case)
  fold_assignments <- caret::createFolds(indices, k = 5)
  
  # Fit the model
  cv_scores <- NULL
  
  # Parallelization
  cl <- makeCluster(no_cores) 
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("astart", "aend", "case", "aevent", "adrug", "aedrug",
                                "agebreaks", "M_1", "M_2", "gamma", "same_expo", "fold_assignments"),
                envir = environment())
  
  cv_scores <- foreach(i = 1:nfolds, .combine = c, .packages = c("dplyr")) %dopar% {
    
    library(splines)
    library(tidyverse)
    library(Matrix)
    library(SCCS)
    library(fda)
    library(caret)
    library(foreach)
    library(doParallel)
    source("/home/xz578/SCCS/SCCS_PLSI_Mspline.R")
    
    case_test <- fold_assignments[[i]]
    case_train <- setdiff(unique(case), case_test)
    
    astart_train <- astart[which(case %in% case_train)]
    aend_train <- aend[which(case %in% case_train)]
    aevent_train <- aevent[which(case %in% case_train)]
    adrug_train <- adrug[which(case %in% case_train),]
    aedrug_train <- aedrug[which(case %in% case_train),]
    
    astart_test <- astart[which(case %in% case_test)]
    aend_test <- aend[which(case %in% case_test)]
    aevent_test <- aevent[which(case %in% case_test)]
    adrug_test <- adrug[which(case %in% case_test),]
    aedrug_test <- aedrug[which(case %in% case_test),]
    
    model_fits <- SCCSFit_nopar(
      astart = astart_train,
      aend = aend_train, 
      case = case_train, 
      aevent = aevent_train,
      adrug = adrug_train, 
      aedrug = aedrug_train,
      agebreaks = agebreaks,
      M_1 = M_1,
      M_2 = M_2,
      gamma = gamma,
      same_expo = same_expo)
    
    prediction <- predict_SCCSFit(
      astart = astart_test,
      aend = aend_test, 
      case = case_test, 
      aevent = aevent_test,
      adrug = adrug_test, 
      aedrug = aedrug_test,
      agebreaks = agebreaks,
      model = model_fits)
    
    return(prediction)
  }
  
  stopCluster(cl)
  
  return(list(gamma = gamma, M_1 = M_1, M_2 = M_2, cv_score = mean(cv_scores, na.rm = T),
              all_scores = cv_scores))
  
}


#######################################################
#####################Inference#########################
#######################################################
#######################################################
## Inference for theta when model results have been obtained
model_summary_theta <- function(astart,
                                aend, 
                                case, 
                                aevent,
                                adrug, 
                                aedrug,
                                agebreaks,
                                model){
  
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  lambda_normalized <- model$lambda_normalized
  knots_B2 <- model$knots_B2
  M_2 <- model$M_2
  M_1 <- model$M_1
  Boundary.knots <- model$`Boundary.knots`
  internal_knots <- model$knots
  theta_coefs <- model$theta
  beta_coefs <- model$beta
  
  #####
  obs_list <- mapply(seq, from = astart, to = aend, SIMPLIFY = FALSE)
  adrug <- cbind(adrug)
  aedrug <- cbind(aedrug)
  
  age_labels <- paste(agebreaks[-length(agebreaks)], agebreaks[-1], sep = "-")
  age_labels[length(agebreaks)] <- paste0(agebreaks[length(agebreaks)],"+")
  age_list <- lapply(obs_list, function(x){
    return(findInterval(x, agebreaks, left.open = T, rightmost.closed = F))
  })
  
  p <- dim(cbind(adrug))[2]
  
  ### Obtain internal knots for S2
  relative_event_time <- -sweep(adrug, 1, aevent,"-")
  u_values <- list()
  risk_times <- list()
  risk_length <- sapply(1:p, function(i) mean(aedrug[,i] - adrug[,i], na.rm = T))
  M_2_indices <- c(1,cumsum(M_2)+1)
  for(i in 1:p){
    relative_times <- c(relative_event_time[,i])
    relative_times <- relative_times[relative_times <= risk_length[i] & relative_times >= 0]
    relative_times <- relative_times[!is.na(relative_times)]
    risk_times[[i]] <- seq(0,risk_length[i],by=1)
    u_values[[i]] <- ns(risk_times[[i]], 
                        Boundary.knots = range(0, risk_length[i]), 
                        knots = knots_B2[[i]],intercept = TRUE) %*% lambda_normalized[M_2_indices[p]:(M_2_indices[p+1]-1)]
  }
  
  u_value_list <- list()
  for(i in 1:p){
    u_value_list[[i]] <- mapply(function(obs, adrug) {
      rel_time <- obs - adrug
      rel_time <- ifelse(rel_time <0 | rel_time > risk_length[i], Inf, rel_time)
      u_value_vec <- u_values[[i]][match(rel_time, risk_times[[i]])]
      u_value_vec[is.na(u_value_vec)] <- 0
      return(u_value_vec)
    }, 
    obs = obs_list, 
    adrug = adrug[,i], 
    SIMPLIFY = FALSE)
  }
  
  sum_list <- Reduce(function(x, y) Map(`+`, x, y), u_value_list)
  combined_u <- unlist(sum_list)
  combined_u <- combined_u[combined_u != 0]
  
  ### Obtain B spline %*% beta for unique u values
  unique_u <- unique(combined_u)
  # return(unique_u)
  Bspline_obj <- bs(unique_u, 
                    knots = internal_knots, 
                    Boundary.knots = Boundary.knots,
                    degree = 3, intercept = T)
  # print(quantile(unique_u))
  ## Obtain B spline 
  all_knots <- c(Boundary.knots[1], internal_knots, Boundary.knots[2])
  nbasis <- length(internal_knots) + 4
  bspline_basis <- create.bspline.basis(rangeval = Boundary.knots,
                                        nbasis = nbasis,
                                        norder = 4,
                                        breaks = all_knots)
  penalty_matrix_2 <- bsplinepen(bspline_basis, Lfdobj = 2)
  
  unique_rr <- Bspline_obj %*% (beta_coefs^2)
  unique_u <- c(0, unique_u)
  unique_rr <- c(1, unique_rr)
  unique_rr[I(unique_rr<0)] <- 0
  ### Map them into u_value_list
  rr_list <- lapply(sum_list, function(x){unique_rr[match(x, unique_u)]})
  
  ### Age risks
  theta_coefs_all <- c(0, theta_coefs)
  rr_age_list <- lapply(age_list, function(x){exp(theta_coefs_all[x+1])})
  
  ### Final relative risks
  RR_list <- Map(`*`, rr_list, rr_age_list)
  
  ### Obtain risks at event
  aevent_indices <- aevent - astart + 1
  event_risk <- mapply(function(RR, index) RR[index], RR_list, aevent_indices)
  
  integral_risk <- sapply(RR_list, function(x)(sum(x)))
  conditional_probs <- mapply(function(RR, sum_risk) RR/sum_risk, RR_list, integral_risk)
  conditional_event_probs <- event_risk/integral_risk
  
  
  integral_risk_by_age <- mapply(function(RR, age_groups){tapply(RR, age_groups, sum)}, RR_list, age_list )
  
  if(!is.matrix(integral_risk_by_age)){
    
    integral_risk_by_age <- lapply(integral_risk_by_age,
                                   function(x){
                                     temp_vec <- rep(0, length(theta_coefs_all))
                                     names(temp_vec) <- 0:(length(theta_coefs_all)-1)
                                     temp_vec[as.numeric(names(x))+1] <- x
                                     return(temp_vec)
                                   })
    integral_risk_by_age <- matrix(unlist(integral_risk_by_age), 
                                   nrow = length(theta_coefs_all), byrow = F)
    
  }
  
  E_vec <- apply(integral_risk_by_age,2,function(x){x/sum(x)})
  E_vec <- E_vec[-1,]
  
  # V_vec <- apply(E_vec, 2, function(x) x-x^2)
  V_vec <- apply(E_vec, 2, function(x) diag(x)-x %*% t(x))
  # SEs <- sqrt(1/rowSums(V_vec))
  Var <- matrix(rowSums(V_vec), nrow = length(agebreaks))
  
  SEs <- sqrt(diag(ginv(Var)))
  
  # S_vec <- table(findInterval(aevent, agebreaks, left.open = T, rightmost.closed = F))
  # S_vec <- S_vec[-1]
  
  result <- data.frame(theta = theta_coefs, 
                       lci = theta_coefs - qnorm(0.975)*SEs,
                       uci = theta_coefs + qnorm(0.975)*SEs,
                       SE = SEs)
  return(result)
}
