
#genereare homo historical dataset
generate.data <- function(num_hist, n, mu, sigma, v) {
  m <- rnorm(num_hist, mu, v)  # historical data mean for num_hist trials
  data_list <- list()
  for (i in 1:num_hist) {
    data <- rnorm(n[i], m[i], sigma)
    data_list[[i]] <- data
  }
  mean_values <- sapply(data_list, mean)
  sd_values <- sapply(data_list, sd)
  se_values <- sd_values / sqrt(n)
  result <- data.frame(
    study = 1:num_hist,
    n = n,
    mean = mean_values,
    sd = sd_values,
    se = se_values
  )
  return(result)
}

#genereate heter historical dataset
generate.data2 <- function(num_hist, n, mu, sigma) {
  data_list <- list()
  for (i in 1:num_hist) {
    data <- rnorm(n[i], mu[i], sigma)
    data_list[[i]] <- data
  }
  mean_values <- sapply(data_list, mean)
  sd_values <- sapply(data_list, sd)
  se_values <- sd_values / sqrt(n)
  result <- data.frame(
    study = 1:num_hist,
    n = n,
    mean = mean_values,
    sd = sd_values,
    se = se_values
  )
  return(result)
}


#codes for ASAM_S(resample) and ASAM_D(ind) 
super_SAM<- function(hist, x_current, N, csd = 0.3, wtype = 'w1', method = 'resample') {
  hist$weight <- sapply(1:nrow(hist), function(i) {
    SAMprior::SAM_weight(
      if.prior = RBesT::mixnorm(c(1, hist$mean[i], hist$se[i]), sigma = hist$sd[i]),        
      theta.h = hist$mean[i], 
      delta = csd,
      method.w = 'LRT',
      data = x_current
    )
  })
  wmethod <- switch(wtype,
                    'w1' = hist$weight / sum(hist$weight),
                    'w2' = hist$weight * hist$n / sum(hist$weight * hist$n),
                    'w3' = hist$n / sum(hist$n),
                    'w4' = LaplacesDemon::rdirichlet(1, rep(1 / nrow(hist), nrow(hist)))
  )
  
  if (method == 'resample') {
    num_samples <- round(N * wmethod)
    num_samples[length(num_samples)] <- N - sum(num_samples[-length(num_samples)])
    result <- NULL
    for (i in 1:nrow(hist)) {
      sam <- SAMprior::SAM_prior(
        if.prior = RBesT::mixnorm(c(1, hist$mean[i], hist$se[i]), sigma = hist$sd[i]),
        nf.prior = RBesT::mixnorm(c(1, hist$mean[i], hist$sd[i]), sigma = hist$sd[i]),
        weight = hist$weight[i], sigma = hist$sd[i]
      )
      post <- RBesT::postmix(sam, x_current)
      if (num_samples[i] < 1) next
      post.mu <- RBesT::rmix(post, num_samples[i])
      result <- c(result, post.mu)
    }
    return(sample(result))
  }else{
    nw_if <- wmethod* hist$weight
    parameters <- list()
    # Append the parameters in the required format
    for (i in 1:nrow(hist)) {
      if (!is.na(nw_if[i])) parameters <- append(parameters, list(c(nw_if[i], hist$mean[i], hist$se[i]))) 
    }
    if (method == 'ind'){ 
      nw_nf <- wmethod*(1 - hist$weight)
      nw_nf[nrow(hist)]=1-sum(nw_if)-sum(nw_nf[-nrow(hist)])
      # Append the parameters in the required format
      for (i in 1:nrow(hist)) {
        if (!is.na(nw_nf[i]))  parameters <- append(parameters, list(c(nw_nf[i], hist$mean[i], hist$sd[i])))
      }
    }
    # Combine the priors using RBesT::mixnorm
    combined_prior <- do.call(RBesT::mixnorm, parameters)
    post <- RBesT::postmix(combined_prior, x_current)
    if (any(is.na(post[1, ]))) return(NULL)
    result <- RBesT::rmix(post, N)
    return(result)
  }
}



#Cumulative Self-adapting prior
## Sigma is Standard deviation in the current trial
ASA.prior <- function(hist, x_current, csd = 0.3, drop = 0.3) {
  mu_vec <- 0
  sig_vec <- 0
  
  for (i in 1:nrow(hist)) { 
    hist$weight[i] <- SAMprior::SAM_weight(
      if.prior = RBesT::mixnorm(c(1, hist$mean[i], hist$se[i]), sigma = hist$sd[i]),
      theta.h = hist$mean[i], 
      delta = csd,
      method.w = 'LRT',
      data = x_current
    )
  }
  
  # Remove rows where weight is below the drop threshold
  hist <- hist[hist$weight > drop, ]
  # If all rows are dropped, use the current
  if (nrow(hist) == 0) {
    ASA_prior <- RBesT::mixnorm(c(1, mean(x_current), sd(x_current)), sigma = sd(x_current))
  } else {
    # Compute mu_vec and sig_vec using the (filtered) hist
    for (i in 1:nrow(hist)) {
      mu_vec <- mu_vec + hist$n[i] * hist$weight[i] * hist$mean[i] / (hist$sd[i]^2)
      sig_vec <- sig_vec + hist$n[i] * hist$weight[i] / (hist$sd[i]^2)
    }
    ASA_prior <- RBesT::mixnorm(c(1, mu_vec / sig_vec, sqrt(1 / sig_vec)), sigma = sd(x_current))
  }
  
  return(ASA_prior)
}


#simulation codes:
NI_weight_combined <- function(method,Hist_r, Hist_p, pcsd, rcsd, ntrial, np, nr, ne, mur, sigmap, sigmar, sigmae, v, delta, lambda, cutoff) {
  options(warn = -1)
  pp.vec <- numeric(ntrial)
  
  # Generate MAP only once to save time
  if (method == 'MAP50'||method == 'MAPSAM') {
    sd_r <- sqrt(sum((Hist_r$n - 1) * Hist_r$sd^2) / sum(Hist_r$n - 1))
    sd_p <- sqrt(sum((Hist_p$n - 1) * Hist_p$sd^2) / sum(Hist_p$n - 1))
    set.seed(555) #fixed MAP prior
    map_r <- suppressWarnings(RBesT::gMAP(cbind(mean, se) ~ 1 | study, weights = n, data = Hist_r, family = gaussian, beta.prior = cbind(0, sd_r), tau.dist = "HalfNormal", tau.prior = cbind(0, sd_r / 2), iter = getOption("RBesT.MC.iter", 6000), warmup = getOption("RBesT.MC.warmup", 1000)))
    mix_r <- RBesT::automixfit(map_r)
    
    map_p <- suppressWarnings(RBesT::gMAP(cbind(mean, se) ~ 1 | study, weights = n, data = Hist_p, family = gaussian, beta.prior = cbind(0, sd_p), tau.dist = "HalfNormal", tau.prior = cbind(0, sd_p / 2), iter = getOption("RBesT.MC.iter", 6000), warmup = getOption("RBesT.MC.warmup", 1000)))
    mix_p <- RBesT::automixfit(map_p)
  }
  
  for (i in 1:ntrial) {
    set.seed(9000 + i)
    NI <- (1 - lambda) * delta
    AS <- delta
    mup <- mur - AS - v
    mue <- mur - NI + v / 2
    
    xp <- rnorm(np, mup, sigmap)
    xr <- rnorm(nr, mur, sigmar)
    xe <- rnorm(ne, mue, sigmae)
    
    post_e <- RBesT::mixnorm(c(1, mean(xe), sd(xe) / sqrt(ne)), sigma = sd(xe))
    post.mue <- RBesT::rmix(post_e, 5000)
    
    if (method == 'resample_w1') {
      post.mur <- super_SAM(Hist_r, xr, N = 5000, csd = rcsd, wtype = 'w1', method = 'resample')
      post.mup <- super_SAM(Hist_p, xp, N = 5000, csd = pcsd, wtype = 'w1', method = 'resample')
    } else if (method == 'w1') {
      post.mur <- super_SAM(Hist_r, xr, N = 5000, csd = rcsd, wtype = 'w1', method = 'ind')
      post.mup <- super_SAM(Hist_p, xp, N = 5000, csd = pcsd, wtype = 'w1', method = 'ind')
    } else if (method == 'w2') {
      post.mur <- super_SAM(Hist_r, xr, N = 5000, csd = rcsd, wtype = 'w2', method = 'ind')
      post.mup <- super_SAM(Hist_p, xp, N = 5000, csd = pcsd, wtype = 'w2', method = 'ind')
    } else if (method == 'w3') {
      post.mur <- super_SAM(Hist_r, xr, N = 5000, csd = rcsd, wtype = 'w3', method = 'ind')
      post.mup <- super_SAM(Hist_p, xp, N = 5000, csd = pcsd, wtype = 'w3', method = 'ind')
    } else if (method == 'w4') {
      post.mur <- super_SAM(Hist_r, xr, N = 5000, csd = rcsd, wtype = 'w4', method = 'ind')
      post.mup <- super_SAM(Hist_p, xp, N = 5000, csd = pcsd, wtype = 'w4', method = 'ind')
    } else if (method == 'NB') {
      post.mur <- RBesT::rmix(RBesT::mixnorm(c(1, mean(xr), sd(xr) / sqrt(nr)), sigma = sd(xr)), 5000)
      post.mup <- RBesT::rmix(RBesT::mixnorm(c(1, mean(xp), sd(xp) / sqrt(np)), sigma = sd(xp)), 5000)
    } else if (method == 'MAP50') {
      prior_r <- RBesT::robustify(mix_r, weight = 0.5, mean = mean(xr), sigma = sd(xr))
      prior_p <- RBesT::robustify(mix_p, weight = 0.5, mean = mean(xp), sigma = sd(xp))
      post_r <- RBesT::postmix(prior_r, xr)
      post_p <- RBesT::postmix(prior_p, xp)
      if (any(is.na(post_r[1, ])) || any(is.na(post_p[1, ]))) next
      post.mur <- RBesT::rmix(post_r, 5000)
      post.mup <- RBesT::rmix(post_p, 5000)
    } else if (method == 'MAPSAM') {
      rSAM= SAMprior::SAM_weight(if.prior=mix_r,theta.h =sum(mix_r[1,] * mix_r[2,]), delta = rcsd,method.w = 'LRT',data = xr)
      prior_r <-RBesT::robustify(mix_r, weight = 1-rSAM, mean = mean(xr), sigma = sd(xr))
      pSAM= SAMprior::SAM_weight(if.prior =mix_p,theta.h = sum(mix_p[1,] * mix_p[2,]), delta = pcsd,method.w = 'LRT',data = xp)
      prior_p <- RBesT::robustify(mix_p, weight = 1-pSAM, mean = mean(xp), sigma = sd(xp))
      
      post_r <- RBesT::postmix(prior_r, xr)
      post_p <- RBesT::postmix(prior_p, xp)
      #if (any(is.na(post_r[1, ])) || any(is.na(post_p[1, ]))) next
      post.mur <- RBesT::rmix(post_r, 5000)
      post.mup <- RBesT::rmix(post_p, 5000)
    }else if (method == 'CSA') {
      prior_r=ASA.prior(Hist_r,xr,rcsd, drop = 0.3)
      prior_p=ASA.prior(Hist_p,xp,pcsd, drop = 0.3)
      post_r <- RBesT::postmix(prior_r, xr)
      post_p <- RBesT::postmix(prior_p, xp)
      post.mur <- RBesT::rmix(post_r, 5000)
      post.mup <- RBesT::rmix(post_p, 5000)
    }else {
      stop("Unknown method")
    }
    if (is.null(post.mur) || is.null(post.mup)) {
      pp.vec[i]= NA
      print(c(i, method,mur,is.null(post.mur),is.null(post.mup)))
    }else{
      pp.vec[i] <- mean(post.mue - post.mur > -NI & post.mup - post.mur < -AS)
    }
    
  }
  #find calibrated cutoff
  if (is.null(cutoff)) {
    C.all <- seq(0, 1, by = 0.001)
    pp.vec.mat <- t(matrix(pp.vec, nrow = length(pp.vec), ncol = length(C.all)))
    # Handle NA values in pp.vec
    rej.mat <- pp.vec.mat > C.all
    # Calculate row means excluding NA values
    C.out <- C.all[which.min(abs(rowMeans(rej.mat, na.rm = TRUE) - 0.0025))]
    cutoff <- C.out
  }
  
  pwr = mean(pp.vec > cutoff, na.rm = TRUE)
  final <- data.frame(method = method, eplison = v, mur = mur, lambda = lambda, pwr = pwr, cutoff = cutoff)
  return(final)
}

