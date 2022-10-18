#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlab <- function(data, nu_0 = 0.0001) {
  ni <- length(unique(data$i))
  nu_1 <- 1000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  SSBCode <- nimble::nimbleCode({
    for (pi in 1:p) {
      delta[pi] ~ dbern(alpha)
      beta[pi] ~ dnorm(0, sd = omega * ((1 - delta[pi]) * nu_0 + delta[pi] * nu_1))
    }
    for (ii in 1:ni) {
      phi[ii] ~ dnorm(mu + inprod(beta[1:p], V[ii, 1:p]), sd = Gamma)
    }
    for (j in 1:N) {
      yhat[j] <- phi[i[j]] / (1 + exp(-t[j] / phi3))
      Y[j] ~ dnorm(yhat[j], sd = sigma)
    }
    sigma ~ dinvgamma(1 / 2, 1 / 2)
    Gamma ~ dinvgamma(1 / 2, 1 / 2)
    omega ~ dinvgamma(1 / 2, 1 / 2)
    alpha ~ dbeta(1, p)
    mu ~ dinvgamma(1 / 2, 1 / 2)

    sigma_prior ~ dinvgamma(1 / 2, 1 / 2)
    Gamma_prior ~ dinvgamma(1 / 2, 1 / 2)
    omega_prior ~ dinvgamma(1 / 2, 1 / 2)
    alpha_prior ~ dbeta(1, p)
    mu_prior ~ dinvgamma(1 / 2, 1 / 2)
    beta_prior ~ dnorm(0, sd = omega_prior * ((1 - delta_prior) * nu_0 + delta_prior * nu_1))
    delta_prior ~ dbern(alpha_prior)
  })

  SSBConsts <- list(
    N = N,
    p = p,
    phi3 = unique(data$phi3),
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    # sigma = 0.1,
    # delta = rep(1, p),
    # Gamma = 0.1,
    # beta = data %>% select(starts_with("beta")) %>% unique() %>% as.numeric(),
    # alpha = 0.2
    nu_0 = nu_0
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y
  )

  SSBInits <- list(
    beta = rep(1, p),
    delta = rep(1, p),
    alpha = 0.1,
    sigma = 0.1,
    Gamma = 1,
    omega = 0.1,
    mu = 15,
    phi = rep(0, ni),
    yhat = data$Y
  )

  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  # SSB$getNodeNames()

  mcmcConf <- nimble::configureMCMC(SSB)

  # mcmcConf$printSamplers()

  # Introducing thinning seems to cause a bug with message: matrix2mv halted because dimensions of variables do not match
  mcmc.out <- nimble::nimbleMCMC(
    code = SSBCode, constants = SSBConsts,
    data = SSBData, inits = SSBInits,
    nchains = 3, nburnin = 0, niter = 5000,
    summary = TRUE, WAIC = TRUE, monitors = c(
      "sigma",
      "beta",
      "delta",
      "Gamma",
      "alpha",
      "yhat",
      "omega",
      "mu",
      "sigma_prior",
      "Gamma_prior",
      "omega_prior",
      "alpha_prior",
      "mu_prior",
      "beta_prior",
      "delta_prior"
    ), samplesAsCodaMCMC = T
  )


  return(mcmc.out)

  # modelInfo <- list(
  #   code = SSBCode, name = "SSB", constants = SSBConsts,
  #   data = SSBData, inits = SSBInits
  # )
  # configure_nimble_slice <- function(model) {
  #   configureMCMC(model, onlySlice = TRUE)
  # }
  # res <- compareMCMCs::compareMCMCs(modelInfo,
  #   MCMCs = c("jags", "nimble", "nimble_slice"),
  #   nimbleMCMCdefs =
  #     list(nimble_slice = "configure_nimble_slice"),
  #   MCMCcontrol = list(
  #     inits = list(a = 1),
  #     niter = 2000,
  #     burnin = 100
  #   )
  # )
}

#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlab_fixed_sigmagamma <- function(data, nu_0 = 0.0001) {
  ni <- length(unique(data$i))
  nu_1 <- 1000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  SSBCode <- nimble::nimbleCode({
    for (pi in 1:p) {
      delta[pi] ~ dbern(alpha)
      beta[pi] ~ dnorm(0, sd = omega * ((1 - delta[pi]) * nu_0 + delta[pi] * nu_1))
    }
    for (ii in 1:ni) {
      phi[ii] ~ dnorm(mu + inprod(beta[1:p], V[ii, 1:p]), sd = Gamma)
    }
    for (j in 1:N) {
      yhat[j] <- phi[i[j]] / (1 + exp(-t[j] / phi3))
      Y[j] ~ dnorm(yhat[j], sd = sigma)
    }
    # sigma ~ dinvgamma(1 / 2, 1 / 2)
    # Gamma ~ dinvgamma(1 / 2, 1 / 2)
    omega ~ dinvgamma(1 / 2, 1 / 2)
    alpha ~ dbeta(1, p)
    mu ~ dinvgamma(1 / 2, 1 / 2)
  })

  SSBConsts <- list(
    N = N,
    p = p,
    phi3 = unique(data$phi3),
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    sigma = 0.1,
    # delta = rep(1, p),
    Gamma = 0.1,
    # beta = data %>% select(starts_with("beta")) %>% unique() %>% as.numeric(),
    # alpha = 0.2
    nu_0 = nu_0
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y
  )

  SSBInits <- list(
    beta = rep(0, p),
    delta = rep(1, p),
    alpha = 0.1,
    sigma = 0.1,
    Gamma = 1,
    omega = 0.1,
    mu = 15,
    phi = rep(0, ni),
    yhat = data$Y
  )

  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  # SSB$getNodeNames()

  mcmcConf <- nimble::configureMCMC(SSB)

  # mcmcConf$printSamplers()

  # Introducing thinning seems to cause a bug with message: matrix2mv halted because dimensions of variables do not match
  mcmc.out <- nimble::nimbleMCMC(
    code = SSBCode, constants = SSBConsts,
    data = SSBData, inits = SSBInits,
    nchains = 3, nburnin = 0, niter = 10000,
    summary = TRUE, WAIC = TRUE, monitors = c(
      # "sigma",
      "beta",
      "delta",
      # "Gamma",
      "alpha",
      "yhat",
      "omega",
      "mu"
    ), samplesAsCodaMCMC = T
  )


  return(mcmc.out)
}


#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlabDelayReparam <- function(data, nu_0 = 0.0001, niter = 5000) {
  nu_1 <- 1000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  SSBCode <- nimble::nimbleCode({


    # Extended model
    for (j in 1:N) {
      # yhat[j] <- Lmax*psi[1] / (1 + exp(-(t[j]-Tmax*psi[2]*phi[i[j]]) / (Tmax*psi[3])))
      yhat[j] <- Lmax * psi1 / (1 + exp(-(t[j] - Tmax * psi2 * phi[i[j]]) / (Tmax * psi3)))
      Y[j] ~ dnorm(yhat[j], sd = sqrt(sigma2))
    }

    for (ii in 1:ni) {
      phi[ii] ~ dnorm(1 + inprod(beta[1:p], V[ii, 1:p]), sd = sqrt(Gamma2))
    }

    # for (m in 1:3){
    #   # psi[m] ~ dnorm(eta[m], sd = sqrt(omega2[m]))
    #   psi[m] ~ dbeta(1, 1)
    # }

    psi1 ~ dbeta(1, 1)
    psi2 ~ dbeta(1, 1)
    psi3 ~ dbeta(1, 1)

    # Priors

    for (pi in 1:p) {
      beta[pi] ~ dnorm(0, sd = (1 - delta[pi]) * nu_0 + delta[pi] * nu_1)
    }

    Gamma2 ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)

    sigma2 ~ dinvgamma(nu_sigma / 2, nu_sigma * lambda_sigma / 2)

    for (pi in 1:p) {
      delta[pi] ~ dbern(alpha)
    }

    alpha ~ dbeta(a, b)

    # for (m in 1:3){
    #   eta[m] ~ dnorm(0, sd = rho_m)
    #   omega[m] ~ dinvgamma(nu_omega_m/2, nu_omega_m*lambda_omega_m/2)
    # }

    sigma2_prior ~ dinvgamma(nu_sigma / 2, nu_sigma * lambda_sigma / 2)
    Gamma2_prior ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)
    # omega_prior ~ dinvgamma(1 / 2, 1 / 2)
    alpha_prior ~ dbeta(a, b)
    # psi_prior ~ dbeta(1, 1)
    psi1_prior ~ dbeta(1, 1)
    psi2_prior ~ dbeta(1, 1)
    psi3_prior ~ dbeta(1, 1)
    # beta_prior ~ dnorm(0, sd = omega_prior * ((1 - delta_prior) * nu_0 + delta_prior * nu_1))
    # delta_prior ~ dbern(alpha_prior)
  })

  SSBConsts <- list(
    N = N,
    p = p,
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    nu_Gamma = 1,
    lambda_Gamma = 1 / 20,
    nu_sigma = 1,
    lambda_sigma = 1,
    a = 1,
    b = 1,
    rho_m = 0.5,
    nu_omega_m = 1,
    lambda_omega_m = 1 / 20,
    Lmax = data$Lmax %>% unique(),
    Tmax = data$Tmax %>% unique(),
    # psi = c(1., 0.25, .15),
    # psi3 = 0.15,
    # psi2 = 0.25,
    # psi1 = 1,
    # sigma2 = 0.1^2,
    # delta = rep(1, p),
    # Gamma2 = 0.1^2,
    # beta = data %>% select(starts_with("beta")) %>% unique() %>% as.numeric(),
    # alpha = 0.2
    nu_0 = nu_0
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y
  )

  SSBInits <- list(
    beta = rep(1, p),
    delta = rep(1, p),
    alpha = 0.1,
    sigma2 = 0.1,
    Gamma2 = 0.1,
    phi = rep(0, ni),
    # psi = rep(0.5, 3),
    psi1 = 0.5,
    psi2 = 0.5,
    psi3 = 0.5,
    yhat = data$Y,
    alpha_prior = 0.5,
    sigma2_prior = 0.1^2,
    psi1_prior = 0.5,
    psi2_prior = 0.5,
    psi3_prior = 0.5,
    Gamma2_prior = 0.1
  )

  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  # SSB$getNodeNames()

  mcmcConf <- nimble::configureMCMC(SSB)

  # mcmcConf$printSamplers()

  # Introducing thinning seems to cause a bug with message: matrix2mv halted because dimensions of variables do not match
  mcmc.out <- nimble::nimbleMCMC(
    code = SSBCode, constants = SSBConsts,
    data = SSBData, inits = SSBInits,
    nchains = 3, nburnin = 0, niter = niter,
    summary = TRUE, WAIC = TRUE, monitors = c(
      "sigma2",
      "delta",
      "Gamma2",
      "alpha",
      "yhat",
      # "psi",
      "psi1",
      "psi2",
      "psi3",
      # "omega",
      # "eta"
      "sigma2_prior",
      "Gamma2_prior",
      # "omega_prior",
      "alpha_prior",
      # "mu_prior",
      # "beta_prior",
      # "delta_prior"
      # "psi_prior",
      "psi1_prior",
      "psi2_prior",
      "psi3_prior",
      "beta"
    ), samplesAsCodaMCMC = T
  )


  return(mcmc.out)
}

#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlabDelay <- function(data, nu_0 = 0.0001, niter = 5000) {
  nu_1 <- 12000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  SSBCode <- nimble::nimbleCode({


    # Extended model
    for (j in 1:N) {
      # yhat[j] <- Lmax*psi[1] / (1 + exp(-(t[j]-Tmax*psi[2]*phi[i[j]]) / (Tmax*psi[3])))
      yhat[j] <- psi[1] / (1 + exp(-(t[j] - phi[i[j]]) / psi[2]))
      Y[j] ~ dnorm(yhat[j], sd = sqrt(sigma2))
    }

    for (ii in 1:ni) {
      phi[ii] ~ dnorm(mu + inprod(beta[1:p], V[ii, 1:p]), sd = sqrt(Gamma2))
    }

    for (m in 1:2) {
      psi[m] ~ T(dnorm(eta[m], sd = sqrt(omega2[m])), 0, )
      # psi[m] ~ dbeta(1, 1)
    }

    # Priors

    for (pi in 1:p) {
      beta[pi] ~ dnorm(0, sd = (1 - delta[pi]) * nu_0 + delta[pi] * nu_1)
    }


    mu ~ dnorm(0, sd = sqrt(sigma2_mu))

    Gamma2 ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)

    # sigma2 ~ dinvgamma(nu_sigma / 2, nu_sigma * lambda_sigma / 2)
    # sigma2 ~ T(dnorm(0, sd = 20), 0, )
    sigma2 ~ dunif(0, 200)


    # speedup: |delta| ~ dbinom(p, alpha)

    for (pi in 1:p) {
      delta[pi] ~ dbern(alpha)
    }

    alpha ~ dbeta(a, b)

    for (m in 1:2) {
      eta[m] ~ dnorm(0, sd = sqrt(rho_m2))
      #   omega2[m] ~ dinvgamma(nu_omega_m/2, nu_omega_m*lambda_omega_m/2)
    }

    # sigma2_prior ~ dinvgamma(nu_sigma / 2, nu_sigma * lambda_sigma / 2)
    # sigma2_prior ~ T(dnorm(0, sd = 20), 0, )
    sigma2_prior ~ dunif(0, 200)
    Gamma2_prior ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)
    # omega2_prior ~ dinvgamma(nu_omega_m/2, nu_omega_m*lambda_omega_m/2)
    alpha_prior ~ dbeta(a, b)
    eta_prior ~ dnorm(0, sd = sqrt(rho_m2))
    psi_prior ~ T(dnorm(eta_prior, sd = sqrt(omega2_prior)), 0, )
    mu_prior ~ dnorm(0, sd = sqrt(sigma2_mu))
    # psi1_prior ~ dbeta(1, 1)
    # psi2_prior ~ dbeta(1, 1)
    beta_prior ~ dnorm(0, sd = sqrt(omega2_prior) * ((1 - delta_prior) * nu_0 + delta_prior * nu_1))
    delta_prior ~ dbern(alpha_prior)
  })

  SSBConsts <- list(
    N = N,
    p = p,
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    nu_Gamma = 1,
    lambda_Gamma = 1,
    nu_sigma = 1.5, # Seem to have improved with respect to nu_sigma = 1
    lambda_sigma = 1,
    a = 1,
    b = p,
    rho_m2 = 1200,
    nu_omega_m = 1,
    lambda_omega_m = 1,
    sigma2_mu = 3000^2,
    # eta = rep(0, 2),
    omega2 = rep(20, 2),
    # eta_prior = 0,
    omega2_prior = 100^2,
    # psi = c(200, 300),
    # mu = 1200,
    # mu_prior = 1200,
    # psi_prior = c(200, 300),
    # psi3 = 0.15,
    # psi2 = 0.25,
    # psi1 = 1,
    # sigma2 = 30,
    # sigma2_prior = 30,
    # delta = rep(1, p),
    # Gamma2 = 200,
    # Gamma2_prior = 200,
    # beta = data %>% select(starts_with("beta")) %>% unique() %>% as.numeric(),
    # alpha = 0.2
    nu_0 = nu_0
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y
  )

  SSBInits <- list(
    # beta = rep(1, p),
    # beta = c(1400, rep(100, 10), rep(1, p - 10)),
    beta = c(rep(100, 10), rep(1, p - 10)),
    delta = rep(1, p),
    alpha = 0.1,
    sigma2 = 0.1,
    # sigma2 = 100,
    # Gamma2 = 0.1,
    Gamma2 = 5000,
    Gamma2 = 5,
    phi = rep(0, ni),
    psi = rep(0.5, 2),
    eta = rep(400, 2),
    # omega2 = rep(20, 2),
    # psi1 = 0.5,
    # psi2 = 0.5,
    yhat = data$Y,
    mu = 1000,
    # alpha_prior = 0.5,
    sigma2_prior = 0.1^2,
    # omega2_prior = 20,
    alpha_prior = 0.5,
    eta_prior = 400,
    psi_prior = 0.5,
    Gamma2_prior = 0.1,
    mu_prior = 100,
    beta_prior = 1,
    delta_prior = 1
  )

  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  # SSB$getNodeNames()

  mcmcConf <- nimble::configureMCMC(SSB)

  # mcmcConf$printSamplers()

  # Introducing thinning seems to cause a bug with message: matrix2mv halted because dimensions of variables do not match
  mcmc.out <- nimble::nimbleMCMC(
    code = SSBCode, constants = SSBConsts,
    data = SSBData, inits = SSBInits,
    nchains = 3, nburnin = 0, niter = niter,
    thin = max(1, round(niter / 5000)),
    summary = TRUE, WAIC = TRUE, monitors = c(
      "sigma2",
      "delta",
      "Gamma2",
      "alpha",
      "yhat",
      "psi",
      "phi",
      # "psi1",
      # "psi2",
      # "omega2",
      # "eta",
      "mu",
      # "eta_prior",
      "sigma2_prior",
      "Gamma2_prior",
      # "omega2_prior",
      "alpha_prior",
      "mu_prior",
      "beta_prior",
      "delta_prior",
      "psi_prior",
      # "psi1_prior",
      # "psi2_prior",
      "beta"
    ), samplesAsCodaMCMC = T
  )


  return(mcmc.out)
}


#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlabDelay_marginal_delta_slow <- function(data, nu_0 = 0.0001, niter = 5000) {
  nu_1 <- 12000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  # logsumexp <- nimbleFunction(
  #   run = function(logx = double(0), logy = double(0)) {
  #     returnType(double(0))
  #     logM <- max(logx, logy)
  #     logres <- logM + log(exp(logx - logM) + exp(logy - logM))
  #     return(logres)
  #   }
  # )

  # dmix2norm <- nimbleFunction(
  #   # Test with curve(dmix2norm(x, p = 0.3, mu1 = -1, mu2 = 1, sigma1 = 0.1, sigma2 = 0.02), from = -2, to = 2, n = 100)
  #   run = function(x = double(0), pp = double(0, default = 0.5),
  #                  mu1 = double(0, default = 0),
  #                  mu2 = double(0, default = 0),
  #                  sigma1 = double(0, default = 1),
  #                  sigma2 = double(0, default = 1),
  #                  log = integer(0, default = 0)) {
  #     returnType(double(0))
  #     print(pp)
  #     # print(log(pp))
  #     logProb1 <- log(pp) + dnorm(x = x, mean = mu1, sd = sigma1, log = TRUE)
  #     logProb2 <- log(1.0 - pp) + dnorm(x = x, mean = mu2, sd = sigma2, log = TRUE)
  #     logProb <- logsumexp(logProb1, logProb2)
  #     if (log) {
  #       return(logProb)
  #     } else {
  #       return(exp(logProb))
  #     }
  #   }
  # )

  dmix2norm_unsafe <- nimbleFunction(
    # Test with curve(dmix2norm(x, p = 0.3, mu1 = -1, mu2 = 1, sigma1 = 0.1, sigma2 = 0.02), from = -2, to = 2, n = 100)
    run = function(x = double(0), pp = double(0, default = 0.5),
                   mu1 = double(0, default = 0),
                   mu2 = double(0, default = 0),
                   sigma1 = double(0, default = 1),
                   sigma2 = double(0, default = 1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      # print(log(pp))
      Prob <- pp * dnorm(x = x, mean = mu1, sd = sigma1) + (1.0 - pp) * dnorm(x = x, mean = mu2, sd = sigma2)
      if (log) {
        return(log(Prob))
      } else {
        return(Prob)
      }
    }
  )

  rmix2norm <- nimbleFunction(
    # Test with hist(rmix2norm(n = 1000, p = 0.3, mu1 = -1, mu2 = 1, sigma1 = 0.1, sigma2 = 0.02))
    run = function(n = integer(0), pp = double(0, default = 0.5),
                   mu1 = double(0, default = 0),
                   mu2 = double(0, default = 0),
                   sigma1 = double(0, default = 1),
                   sigma2 = double(0, default = 1)) {
      returnType(double(0))
      if (n != 1) print("rmix2norm only allows n = 1; using n = 1.")

      u <- runif(1, 0, 1)
      # delta <- rbinom(n = n, size = 1, prob = p)

      # if(n == 1) {
      # res = rnorm(n = n, mean = )
      if (u < pp) {
        res <- rnorm(n = 1, mean = mu1, sd = sigma1)
      } else {
        res <- rnorm(n = 1, mean = mu2, sd = sigma2)
      }
      # }
      # else {
      #   res <- rnorm(n = n, mean = mu1 * delta + mu2 * (1-delta),
      #                sd = sigma1 * delta + sigma2 * (1-delta))
      # }
      return(res)
    }
  )

  # assign("logsumexp", dmix2norm, envir = .GlobalEnv)
  # assign("dmix2norm", dmix2norm, envir = .GlobalEnv)
  assign("dmix2norm_unsafe", dmix2norm_unsafe, envir = .GlobalEnv)
  # assign("rmix2norm", rmix2norm, envir = .GlobalEnv)
  assign("rmix2norm_unsafe", rmix2norm, envir = .GlobalEnv)


  SSBCode <- nimble::nimbleCode({

    # Extended model
    for (j in 1:N) {
      # yhat[j] <- Lmax*psi[1] / (1 + exp(-(t[j]-Tmax*psi[2]*phi[i[j]]) / (Tmax*psi[3])))
      yhat[j] <- psi[1] / (1 + exp(-(t[j] - phi[i[j]]) / psi[2]))
      Y[j] ~ dnorm(yhat[j], sd = sqrt(sigma2))
    }

    for (ii in 1:ni) {
      phi[ii] ~ dnorm(mu + inprod(beta[1:p], V[ii, 1:p]), sd = sqrt(Gamma2))
    }

    for (m in 1:2) {
      psi[m] ~ T(dnorm(eta[m], sd = sqrt(omega2[m])), 0, )
      # psi[m] ~ dbeta(1, 1)
    }

    # Priors

    for (pi in 1:p) {
      # beta[pi] ~ dmix2norm(0, sd = (1 - delta[pi]) * nu_0 + delta[pi] * nu_1)
      beta[pi] ~ dmix2norm_unsafe(pp = alpha, mu1 = 0, mu2 = 0, sigma1 = nu_1, sigma2 = nu_0)
    }


    mu ~ dnorm(0, sd = sqrt(sigma2_mu))

    Gamma2 ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)

    sigma2 ~ dunif(0, 200)

    alpha ~ dbeta(a, b)

    for (m in 1:2) {
      eta[m] ~ dnorm(0, sd = sqrt(rho_m2))
      #   omega2[m] ~ dinvgamma(nu_omega_m/2, nu_omega_m*lambda_omega_m/2)
    }

    # sigma2_prior ~ dinvgamma(nu_sigma / 2, nu_sigma * lambda_sigma / 2)
    # sigma2_prior ~ T(dnorm(0, sd = 20), 0, )
    sigma2_prior ~ dunif(0, 200)
    Gamma2_prior ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)
    # omega2_prior ~ dinvgamma(nu_omega_m/2, nu_omega_m*lambda_omega_m/2)
    alpha_prior ~ dbeta(a, b)
    eta_prior ~ dnorm(0, sd = sqrt(rho_m2))
    psi_prior ~ T(dnorm(eta_prior, sd = sqrt(omega2_prior)), 0, )
    mu_prior ~ dnorm(0, sd = sqrt(sigma2_mu))
    # psi1_prior ~ dbeta(1, 1)
    # psi2_prior ~ dbeta(1, 1)
    beta_prior ~ dnorm(0, sd = sqrt(omega2_prior) * ((1 - delta_prior) * nu_0 + delta_prior * nu_1))
    delta_prior ~ dbern(alpha_prior)
  })

  SSBConsts <- list(
    N = N,
    p = p,
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    nu_Gamma = 1,
    lambda_Gamma = 1,
    nu_sigma = 1, #1.5 Seem to have improved with respect to nu_sigma = 1
    lambda_sigma = 1,
    a = 1,
    b = p,
    rho_m2 = 1200,
    #nu_omega_m = 1,
    #lambda_omega_m = 1,
    sigma2_mu = 3000^2,
    # eta = rep(0, 2),
    omega2 = rep(20, 2),
    # eta_prior = 0,
    omega2_prior = 100^2,
    # psi = c(200, 300),
    # mu = 1200,
    # mu_prior = 1200,
    # psi_prior = c(200, 300),
    # psi3 = 0.15,
    # psi2 = 0.25,
    # psi1 = 1,
    # sigma2 = 30,
    # sigma2_prior = 30,
    # delta = rep(1, p),
    # Gamma2 = 200,
    # Gamma2_prior = 200,
    # beta = data %>% select(starts_with("beta")) %>% unique() %>% as.numeric(),
    # alpha = 0.2
    nu_0 = nu_0
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y
  )

  SSBInits <- list(
    # beta = rep(1, p),
    # beta = c(1400, rep(100, 10), rep(1, p - 10)),
    beta = c(rep(100, 10), rep(1, p - 10)),
    alpha = 0.1, #0.1
    sigma2 = 100, #0.1
    # sigma2 = 100,
    # Gamma2 = 0.1,
    # Gamma2 = 5000,
    Gamma2 = 500, #50
    eta = rep(400, 2),
    # omega2 = rep(20, 2),
    # psi1 = 0.5,
    # psi2 = 0.5,
    yhat = data$Y,
    mu = 1400, #1000
    # alpha_prior = 0.5,
    sigma2_prior = 0.1^2,
    # omega2_prior = 20,
    alpha_prior = 0.5,
    eta_prior = 400,
    psi_prior = 0.5,
    Gamma2_prior = 0.1,
    mu_prior = 100,
    beta_prior = 1,
    delta_prior = 1
  ) %>%
    (function(x){
      x$phi = rep(x$mu, ni)
      x$psi = x$eta
      return(x)
    })


  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  # SSB$getNodeNames()

  mcmcConf <- nimble::configureMCMC(SSB)

  # mcmcConf$printSamplers()

  # Introducing thinning seems to cause a bug with message: matrix2mv halted because dimensions of variables do not match
  mcmc.out <- nimble::nimbleMCMC(
    code = SSBCode, constants = SSBConsts,
    data = SSBData, inits = SSBInits,
    nchains = 3, nburnin = 0, niter = niter,
    thin = max(1, round(niter / 5000)),
    summary = TRUE, WAIC = TRUE, monitors = c(
      "sigma2",
      "Gamma2",
      "alpha",
      "yhat",
      "psi",
      "phi",
      # "psi1",
      # "psi2",
      # "omega2",
      # "eta",
      "mu",
      # "eta_prior",
      "sigma2_prior",
      "Gamma2_prior",
      # "omega2_prior",
      "alpha_prior",
      "mu_prior",
      "beta_prior",
      "psi_prior",
      # "psi1_prior",
      # "psi2_prior",
      "beta"
    ), samplesAsCodaMCMC = T
  )


  return(mcmc.out)
}


#' Title
#'
#' @param data
#' @param nu_0
#'
#' @return
#' @export
#'
#' @examples
OneVariableParamSpikeAndSlabDelay_marginal_delta_setup <- function(data, nu_0 = 0.0001, niter = 5000) {
  nu_1 <- 12000
  p <- data %>%
    names() %>%
    grepl("^V", x = .) %>%
    sum()
  ni <- length(unique(data$i))
  N <- nrow(data)

  dmix2norm_unsafe <- nimbleFunction(
    # Test with curve(dmix2norm(x, p = 0.3, mu1 = -1, mu2 = 1, sigma1 = 0.1, sigma2 = 0.02), from = -2, to = 2, n = 100)
    run = function(x = double(0), pp = double(0, default = 0.5),
                   mu1 = double(0, default = 0),
                   mu2 = double(0, default = 0),
                   sigma1 = double(0, default = 1),
                   sigma2 = double(0, default = 1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      # print(log(pp))
      Prob <- pp * dnorm(x = x, mean = mu1, sd = sigma1) + (1.0 - pp) * dnorm(x = x, mean = mu2, sd = sigma2)
      if (log) {
        return(log(Prob))
      } else {
        return(Prob)
      }
    }
  )

  rmix2norm <- nimbleFunction(
    # Test with hist(rmix2norm(n = 1000, p = 0.3, mu1 = -1, mu2 = 1, sigma1 = 0.1, sigma2 = 0.02))
    run = function(n = integer(0), pp = double(0, default = 0.5),
                   mu1 = double(0, default = 0),
                   mu2 = double(0, default = 0),
                   sigma1 = double(0, default = 1),
                   sigma2 = double(0, default = 1)) {
      returnType(double(0))
      if (n != 1) print("rmix2norm only allows n = 1; using n = 1.")

      u <- runif(1, 0, 1)
      if (u < pp) {
        res <- rnorm(n = 1, mean = mu1, sd = sigma1)
      } else {
        res <- rnorm(n = 1, mean = mu2, sd = sigma2)
      }
      return(res)
    }
  )

  assign("dmix2norm_unsafe", dmix2norm_unsafe, envir = .GlobalEnv)
  assign("rmix2norm_unsafe", rmix2norm, envir = .GlobalEnv)


  SSBCode <- nimble::nimbleCode({

    # Extended model
    for (j in 1:N) {
      yhat[j] <- psi[1] / (1 + exp(-(t[j] - phi[i[j]]) / psi[2]))
      Y[j] ~ dnorm(yhat[j], sd = sqrt(sigma2))
    }

    for (ii in 1:ni) {
      phi[ii] ~ dnorm(mu + inprod(beta[1:p], V[ii, 1:p]), sd = sqrt(Gamma2))
    }

    for (m in 1:2) {
      psi[m] ~ T(dnorm(eta[m], sd = sqrt(omega2[m])), 0, )
    }

    # Priors
    for (pi in 1:p) {
      beta[pi] ~ dmix2norm_unsafe(pp = alpha, mu1 = 0, mu2 = 0, sigma1 = nu_1, sigma2 = nu_0)
    }


    mu ~ dnorm(0, sd = sqrt(sigma2_mu))
    Gamma2 ~ dinvgamma(nu_Gamma / 2, nu_Gamma * lambda_Gamma / 2)
    sigma2 ~ dunif(0, 200)
    alpha ~ dbeta(a, b)

    for (m in 1:2) {
      eta[m] ~ dnorm(0, sd = sqrt(rho_m2))
    }
  })

  SSBConsts <- list(
    N = N,
    p = p,
    ni = ni,
    i = data$i,
    nu_1 = nu_1,
    nu_Gamma = 1,
    lambda_Gamma = 1,
    nu_sigma = 1, # 1.5 Seem to have improved with respect to nu_sigma = 1
    lambda_sigma = 1,
    a = 1,
    b = p,
    rho_m2 = 1200,
    #nu_omega_m = 1,
    #lambda_omega_m = 1,
    sigma2_mu = 3000^2,
    omega2 = rep(20, 2)
    #omega2_prior = 100^2
  )

  SSBData <- list(
    t = data$times,
    V = data %>%
      names() %>%
      grepl("^V", x = .) %>%
      (function(mask) data[, mask]) %>%
      unique(),
    Y = data$Y,
    nu_0 = nu_0
  )

  SSBInits <- list(
    beta = c(rep(100, 10), rep(1, p - 10)),
    alpha = 0.1, #0.1
    sigma2 = 100, #0.1
    Gamma2 = 500, #50
    #phi = rep(0, ni),
    #psi = rep(0.5, 2),
    eta = rep(400, 2),
    yhat = data$Y,
    mu = 1400 #1000
  ) %>%
    (function(x){
      x$phi = rep(x$mu, ni)
      x$psi = x$eta
      return(x)
    })

  SSB <- nimble::nimbleModel(
    code = SSBCode, name = "SSB", constants = SSBConsts,
    data = SSBData, inits = SSBInits
  )

  mcmcConf <- configureMCMC(model = SSB, monitors = c(
    "sigma2",
    "Gamma2",
    "alpha",
    "yhat",
    "psi",
    "phi",
    "mu",
    "beta"
  ))

  Rmcmc <- buildMCMC(SSB)

  Cmcmc <- compileNimble(Rmcmc, SSB)

  return(Cmcmc)
}

#' Title
#'
#' @param Cmcmc
#' @param nu_0
#' @param niter
#'
#' @return
#' @export
#'
#' @examples
TimeOneVariableParamSpikeAndSlabDelay_marginal_delta <- function(Cmcmc, nu_0 = 0.0001, niter = 5000) {
  Cmcmc$SSB$nu_0 <- nu_0
  Cmcmc$Rmcmc$run(
    niter = niter,
    time = TRUE,
    thin = max(1, round(niter / 1000)), nburnin = floor(niter / 2)
  )
  return(Cmcmc$Rmcmc$getTimes() %>% sum())
}

#' Get a sample from the posterior distribution given the data in Cmcmc
#'
#' The function should not output more than 1000 iterations, to save on memory storage. The first half of the iterations are discarded for warmup, and the chain is thinned to output around 1000 iterations.
#'
#' @param Cmcmc
#' @param nu_0
#' @param niter
#'
#' @return
#' @export
#'
#' @examples
FitOneVariableParamSpikeAndSlabDelay_marginal_delta <- function(Cmcmc, nu_0 = 0.0001, niter = 5000) {
  Cmcmc$SSB$nu_0 <- nu_0
  runMCMC(mcmc = Cmcmc$Rmcmc,
          niter = niter,
          samplesAsCodaMCMC = T,
          thin = max(1, round(niter / (2*1000))), nburnin = floor(niter / 2)
  ) %>%
    as_tibble()
}
