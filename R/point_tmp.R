# Mean für log und log-shift aus aggregierten Kovariaten mittels der Schätzung von Ed

est_area_1_g <- function(data_smp, area_size, area_size_names, domains, gewicht_o_grenze, gewicht_u_grenze,
                         x_mean_d, x_sd_d, x_cov_d, formel, sel){


  # optpar bestimmen ------
  transformation <- choose_trafo(sel)

  smp_domains = unlist(data_smp[domains])
  X_smp       = model.matrix(formel, data_smp)
  Y_smp       = as.matrix(data_smp[paste(formel[2])])

  par_bound=c(-1,2)

  if(sel != 2 && sel !=5 && sel != 12 && sel != 10 && sel != 13 && sel !=6) # if no estimation, no optmimization
  {
    par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, sel = sel)
    optpar    = stats::optimize(generic_opt, y = Y_smp , dat = X_smp, form = formel, idD = smp_domains,
                                meth = method, interval = par_bound, maximum = F, sel = sel)$minimum
  }
  if(sel == 12)
  {
    optpar <- optim(par = c(1,1), fn = generic_opt, dat = X_smp, form=formel, idD=smp_domains, meth=method, y=Y_smp)$par
  }
  if(sel == 10 | sel == 13 ) #johnson
  {
    m <- select_john(z = 1, y = resid(lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains), method="REML")))
    if(sel == 10)
    {
      par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, m = m, sel = sel)
      optpar    = sbplx(x0 = c(1,1), fn = generic_opt, upper = c(Inf,Inf), lower = c(-Inf,0.001),
                        control = list(xtol_rel = 1e-16), nl.info = FALSE, m = m, dat = X_smp, form = formel,
                        idD = smp_domains, meth = method, y = Y_smp)$par
    }
    else
    {
      optpar    = optim(par = c(1,1,1,1), fn = generic_opt, upper = c(Inf,Inf,Inf,Inf), lower = c(-Inf, 0.01, -Inf, 0.00001),
                        method = "L-BFGS-B", m = m, dat = X_smp, form = formel, idD = smp_domains,
                        meth = method, y = Y_smp)$par
    }
  }
  if(sel == 2 | sel == 5 | sel == 6)
  {
    optpar = NULL
  }

  # Spezfikationen
  smp_domains = unlist(data_smp[domains][,1])
  #n_smp       = as.vector(table(smp_domains))
  n_smp <- NULL
  for(i in 1:(length(area_size))){
    n_smp[i] <- nrow(data_smp[data_smp[domains] == area_size_names[i],])
  }

  X_smp       = model.matrix(formel,data_smp)

  Y_smp       = as.matrix(data_smp[paste(formel[2])])
  tmp         = transformation(y=Y_smp, inv=F)
  Y_smp       = tmp$y
  par_m       = tmp$m

  # Schätzen des Modells
  nlm_smp     = nlme::lme(Y_smp ~ -1 + X_smp,random = ~1|as.factor(smp_domains),method="REML", control=(msMaxIter=500))

  beta        = nlme::fixed.effects(nlm_smp)
  u_est_d     = nlme::random.effects(nlm_smp)[,1]
  sigmae2est  = nlm_smp$sigma^2
  sigmau2est  = as.numeric(nlme::VarCorr(nlm_smp)[1,1])
  gamma_est_d = sigmau2est / (sigmau2est + sigmae2est / n_smp)

  l <- 1
  u_est_d_tmp <- rep(NA, length(area_size))

  for(i in 1:length(area_size)){
    if(n_smp[i] == 0){
      u_est_d_tmp[i] <- 0
    } else{
      u_est_d_tmp[i] <- u_est_d[l]
      l <- l+1
    }
  }

  u_est_d <- u_est_d_tmp

  # Schätzen von Ed
  input_est_Ed        = data.frame(X_smp %*% beta, data_smp[domains])
  names(input_est_Ed) = c("x", "idD")

  x_mean_input        = x_mean_d %*% beta
  #x_sd_input_t        = sqrt(x_sd_d^2 %*% beta^2) # nur Varianzen
  x_sd_input          = sqrt(x_cov_d %*% as.numeric(beta %*% t(beta)))


  Res <- Ed_estimation_gruppenweise(data_smp = input_est_Ed, domains = domains,
                                    gewichtete_den_grenze = gewicht_o_grenze, gewichtete_den_grenze_u = gewicht_u_grenze,
                                    area_size = area_size, area_size_names = area_size_names, x_mean_d = x_mean_input,
                                    x_sd_d = x_sd_input,
                                    plot_TF = F)

  # Bestimmen Area means
  alpha_est_d         <- (sigmau2est * (1 - gamma_est_d) + sigmae2est)/2
  tau_est_area_st1    <- 1/area_size * (Res$E_d_density_est_kor_pop * exp(u_est_d + alpha_est_d))

  out <- list(tau_est_area_st1, sigmae2est, sigmau2est, nlm_smp)
  names(out) <- c("tau_est_area_st1", "sigmae2est_b", "sigmau2est_b", "model_trafo")
  return(out)
}
est_area_1_g_ls <- function(data_smp, area_size, area_size_names, domains, gewicht_o_grenze, gewicht_u_grenze,
                            x_mean_d, x_sd_d, x_cov_d, formel, sel, method = 6){

  # optpar bestimmen ------
  transformation <- choose_trafo(sel)

  smp_domains = unlist(data_smp[domains])
  X_smp       = model.matrix(formel, data_smp)
  Y_smp       = as.matrix(data_smp[paste(formel[2])])

  par_bound=c(-1,2)

  if(sel != 2 && sel !=5 && sel != 12 && sel != 10 && sel != 13 && sel !=6) # if no estimation, no optmimization
  {
    par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, sel = sel)
    optpar    = try(stats::optimize(generic_opt, y = Y_smp , dat = X_smp, form = formel, idD = smp_domains,
                                    meth = method, interval = par_bound, maximum = F, sel = sel)$minimum)
  }
  if(sel == 12)
  {
    optpar <- optim(par = c(1,1), fn = generic_opt, dat = X_smp, form=formel, idD=smp_domains, meth=method, y=Y_smp)$par
  }
  if(sel == 10 | sel == 13 ) #johnson
  {
    m <- select_john(z = 1, y = resid(lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains), method="REML")))
    if(sel == 10)
    {
      par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, m = m, sel = sel)
      optpar    = sbplx(x0 = c(1,1), fn = generic_opt, upper = c(Inf,Inf), lower = c(-Inf,0.001),
                        control = list(xtol_rel = 1e-16), nl.info = FALSE, m = m, dat = X_smp, form = formel,
                        idD = smp_domains, meth = method, y = Y_smp)$par
    }
    else
    {
      optpar    = optim(par = c(1,1,1,1), fn = generic_opt, upper = c(Inf,Inf,Inf,Inf), lower = c(-Inf, 0.01, -Inf, 0.00001),
                        method = "L-BFGS-B", m = m, dat = X_smp, form = formel, idD = smp_domains,
                        meth = method, y = Y_smp)$par
    }
  }
  if(sel == 2 | sel == 5 | sel == 6)
  {
    optpar = NULL
  }

  if(class(optpar) == "try-error"){
    print("ein Bootstrap ist nicht konvergiert")
    tau_est_area_st1 <- rep(NA, length(area_size))
    optpar <- NA
    sigmae2est_b <- NA
    sigmau2est_b <- NA
    out <- list(tau_est_area_st1, optpar, sigmae2est_b, sigmau2est_b)
    names(out) <- c("tau_est_area_st1", "optpar_b", "sigmae2est_b", "sigmau2est_b")
    return(out)
  }else{

    # Spezfikationen
    smp_domains = unlist(data_smp[domains][,1])
    #n_smp       = as.vector(table(smp_domains))
    n_smp <- NULL
    for(i in 1:(length(area_size))){
      n_smp[i] <- nrow(data_smp[data_smp[domains] == area_size_names[i],])
    }

    X_smp       = model.matrix(formel,data_smp)

    Y_smp       = as.matrix(data_smp[paste(formel[2])])
    tmp         = transformation(l = optpar, y=Y_smp, inv=F)
    Y_smp       = tmp$y
    par_m       = tmp$m

    # Schätzen des Modells
    nlm_smp     = nlme::lme(Y_smp ~ -1 + X_smp,random = ~1|as.factor(smp_domains),method="REML", control=(msMaxIter=5000))

    beta        = nlme::fixed.effects(nlm_smp)
    u_est_d     = nlme::random.effects(nlm_smp)[,1]
    sigmae2est  = nlm_smp$sigma^2
    sigmau2est  = as.numeric(nlme::VarCorr(nlm_smp)[1,1])
    gamma_est_d = sigmau2est / (sigmau2est + sigmae2est / n_smp)

    l <- 1
    u_est_d_tmp <- rep(NA, length(area_size))

    for(i in 1:length(area_size)){
      if(n_smp[i] == 0){
        u_est_d_tmp[i] <- 0
      } else{
        u_est_d_tmp[i] <- u_est_d[l]
        l <- l+1
      }
    }

    u_est_d <- u_est_d_tmp

    # Schätzen von Ed
    input_est_Ed        = data.frame(X_smp %*% beta, data_smp[domains])
    names(input_est_Ed) = c("x", "idD")

    x_mean_input        = x_mean_d %*% beta
    #x_sd_input_t        = sqrt(x_sd_d^2 %*% beta^2) # nur Varianzen
    x_sd_input          = sqrt(x_cov_d %*% as.numeric(beta %*% t(beta)))

    Res <- Ed_estimation_gruppenweise(data_smp = input_est_Ed, domains = domains, area_size_names = area_size_names,
                                      gewichtete_den_grenze = gewicht_o_grenze, gewichtete_den_grenze_u = gewicht_u_grenze,
                                      area_size = area_size, x_mean_d = x_mean_input, x_sd_d = x_sd_input,
                                      plot_TF = F)

    # Bestimmen Area means
    alpha_est_d         <- (sigmau2est * (1 - gamma_est_d) + sigmae2est)/2
    tau_est_area_st1    <- 1/area_size * (Res$E_d_density_est_kor_pop * exp(u_est_d + alpha_est_d)) - optpar

    out <- list(tau_est_area_st1, optpar, sigmae2est, sigmau2est, nlm_smp)
    names(out) <- c("tau_est_area_st1", "optpar_b", "sigmae2est_b", "sigmau2est_b", "model_trafo")
    return(out)}
}
# hier wird das SAE package angewandt
est_area_1_bhf <- function(data_smp, area_size, area_size_names, domains, x_mean_d, formel, sel, method = 6){

  # optpar bestimmen ------
  transformation <- choose_trafo(sel = 2)

  smp_domains = unlist(data_smp[domains])
  X_smp       = model.matrix(formel, data_smp)
  Y_smp       = as.matrix(data_smp[paste(formel[2])])

  if(sel == 2 | sel == 5 | sel == 6)
  {
    optpar = NULL
  }

  # Spezfikationen
  smp_domains = unlist(data_smp[domains][,1])

  X_smp       = model.matrix(formel,data_smp)
  n_smp <- NULL
  for(i in 1:(length(area_size))){
    n_smp[i] <- nrow(data_smp[data_smp[domains] == area_size_names[i],])
  }

  Y_smp       = as.matrix(data_smp[paste(formel[2])])
  tmp         = transformation(y=Y_smp, inv=F)
  Y_smp       = tmp$y
  par_m       = tmp$m

  # Schätzen des Modells
  nlm_smp     = nlme::lme(Y_smp ~ -1 + X_smp,random = ~1|as.factor(smp_domains),method="REML")
  beta        = nlme::fixed.effects(nlm_smp)

  data_smp$domains_ <- data_smp[,domains]
  tau_est_area_bhf <-
    sae::eblupBHF(formula = formel, dom = domains_,
                  popnsize = cbind(1:nrow(x_mean_d), area_size),
                  meanxpop = cbind(1:nrow(x_mean_d), x_mean_d[,2:ncol(x_mean_d)]),
                  data = data_smp, selectdom =1:length(area_size))

  return(tau_est_area_bhf$eblup$eblup)
}

est_area_naive <- function(data_smp, domains, optpar, formel, x_mean_d, transformation, sel, id, method = 6){

  # optpar bestimmen ------
  transformation <- choose_trafo(sel)

  smp_domains = unlist(data_smp[domains])
  X_smp       = model.matrix(formel, data_smp)
  Y_smp       = as.matrix(data_smp[paste(formel[2])])

  par_bound=c(-1,2)

  if(sel != 2 && sel !=5 && sel != 12 && sel != 10 && sel != 13 && sel !=6) # if no estimation, no optmimization
  {
    par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, sel = sel)
    optpar    = stats::optimize(generic_opt, y = Y_smp , dat = X_smp, form = formel, idD = smp_domains,
                                meth = method, interval = par_bound, maximum = F, sel = sel)$minimum
  }
  if(sel == 12)
  {
    optpar <- optim(par = c(1,1), fn = generic_opt, dat = X_smp, form=formel, idD=smp_domains, meth=method, y=Y_smp)$par
  }
  if(sel == 10 | sel == 13 ) #johnson
  {
    m <- select_john(z = 1, y = resid(lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains), method="REML")))
    if(sel == 10)
    {
      par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, m = m, sel = sel)
      optpar    = sbplx(x0 = c(1,1), fn = generic_opt, upper = c(Inf,Inf), lower = c(-Inf,0.001),
                        control = list(xtol_rel = 1e-16), nl.info = FALSE, m = m, dat = X_smp, form = formel,
                        idD = smp_domains, meth = method, y = Y_smp)$par
    }
    else
    {
      optpar    = optim(par = c(1,1,1,1), fn = generic_opt, upper = c(Inf,Inf,Inf,Inf), lower = c(-Inf, 0.01, -Inf, 0.00001),
                        method = "L-BFGS-B", m = m, dat = X_smp, form = formel, idD = smp_domains,
                        meth = method, y = Y_smp)$par
    }
  }
  if(sel == 2 | sel == 5 | sel == 6)
  {
    optpar = NULL
  }

  # Spefizikationen-----
  area_size <- as.numeric(table(data_pop[domains]))
  area_size_names <- names(table(data_pop[domains]))

  smp_domains = unlist(data_smp[domains][,1])

  n_smp       = c()
  for(i in 1:length(area_size)){
    n_smp[i]  = length(which(data_smp[domains] == area_size_names[i]))
  }

  X_smp       = model.matrix(formel,data_smp)

  Y_smp       = as.matrix(data_smp[paste(formel[2])])
  tmp         = transformation(y = Y_smp, l = optpar, inv=F)
  Y_smp       = tmp$y
  par_m       = tmp$m

  #X_r <- as.data.frame(matrix(NA, nrow = nrow(data_pop) - nrow(data_smp) , ncol = ncol(data_pop)))
  #names(X_r) <- names(data_pop)

  # l <- 0
  # for(i in 1 : length(area_size)){
  #   tmp_pop <- data_pop[data_pop[, domains] == area_size_names[i],]
  #   tmp_smp <- data_smp[data_smp[, domains] == area_size_names[i],]
  #   if(n_smp[i] != 0){
  #     X_r[(l+1) : (l+(area_size[i]-n_smp[i])),] <- tmp_pop[- which(tmp_pop$id %in% tmp_smp$id),]
  #   }else{
  #     X_r[(l+1) : (l+(area_size[i]-n_smp[i])),] <- tmp_pop
  #   }
  #   l <- l+(area_size[i] - n_smp[i])
  # }

  # Schätzen des Modells
  nlm_smp     = nlme::lme(Y_smp ~ -1 + X_smp,random = ~1|as.factor(smp_domains),method="REML")

  beta        = nlme::fixed.effects(nlm_smp)
  u_est_d     = nlme::random.effects(nlm_smp)
  sigmae2est  = nlm_smp$sigma^2
  sigmau2est  = as.numeric(nlme::VarCorr(nlm_smp)[1,1])
  gamma_est_d = sigmau2est / (sigmau2est + sigmae2est / n_smp)

  #m_di <- NULL
  w_est_ohnebc_dri_trafo <- NULL

  for(j in 1: length(area_size)){
    if(n_smp[j] != 0){
      w_est_ohnebc_dri_trafo[j] <- x_mean_d[j,] %*% beta + u_est_d[row.names(u_est_d) == j, ]
    }else{
      w_est_ohnebc_dri_trafo[j] <- x_mean_d[j,] %*% beta
    }
  }

  tau_est_ohnebc_d <- NULL

  w_est_ohnebc_dr <- transformation(optpar, w_est_ohnebc_dri_trafo, inv = TRUE, m = 0)$y

  for(j in 1: length(area_size)){
    w_ds_ <- data_smp[data_smp[, domains] == (1:length(area_size))[j],]
    w_ds  <- as.matrix(w_ds_[paste(formel[2])])
    w_est_ohnebc_dr_total <- w_est_ohnebc_dr[j] * (area_size[j] - length(w_ds))

    tau_est_ohnebc_d[j] <- 1 / (area_size[j]) * (sum(w_ds) + w_est_ohnebc_dr_total)
  }

  return(tau_est_ohnebc_d)
}


est_area_naive_bc <- function(data_smp, domains, optpar, formel, x_mean_d, transformation, sel, id, method = 6){

  # optpar bestimmen ------
  transformation <- choose_trafo(sel)

  smp_domains = unlist(data_smp[domains])
  X_smp       = model.matrix(formel, data_smp)
  Y_smp       = as.matrix(data_smp[paste(formel[2])])

  par_bound=c(-1,2)

  if(sel != 2 && sel !=5 && sel != 12 && sel != 10 && sel != 13 && sel !=6) # if no estimation, no optmimization
  {
    par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, sel = sel)
    optpar    = stats::optimize(generic_opt, y = Y_smp , dat = X_smp, form = formel, idD = smp_domains,
                                meth = method, interval = par_bound, maximum = F, sel = sel)$minimum
  }
  if(sel == 12)
  {
    optpar <- optim(par = c(1,1), fn = generic_opt, dat = X_smp, form=formel, idD=smp_domains, meth=method, y=Y_smp)$par
  }
  if(sel == 10 | sel == 13 ) #johnson
  {
    m <- select_john(z = 1, y = resid(lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains), method="REML")))
    if(sel == 10)
    {
      par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, m = m, sel = sel)
      optpar    = sbplx(x0 = c(1,1), fn = generic_opt, upper = c(Inf,Inf), lower = c(-Inf,0.001),
                        control = list(xtol_rel = 1e-16), nl.info = FALSE, m = m, dat = X_smp, form = formel,
                        idD = smp_domains, meth = method, y = Y_smp)$par
    }
    else
    {
      optpar    = optim(par = c(1,1,1,1), fn = generic_opt, upper = c(Inf,Inf,Inf,Inf), lower = c(-Inf, 0.01, -Inf, 0.00001),
                        method = "L-BFGS-B", m = m, dat = X_smp, form = formel, idD = smp_domains,
                        meth = method, y = Y_smp)$par
    }
  }
  if(sel == 2 | sel == 5 | sel == 6)
  {
    optpar = NULL
  }

  # Spefizikationen-----
  area_size <- as.numeric(table(data_pop[domains]))
  area_size_names <- names(table(data_pop[domains]))

  smp_domains = unlist(data_smp[domains][,1])

  n_smp       = c()
  for(i in 1:length(area_size)){
    n_smp[i]  = length(which(data_smp[domains] == area_size_names[i]))
  }

  X_smp       = model.matrix(formel,data_smp)

  Y_smp       = as.matrix(data_smp[paste(formel[2])])
  tmp         = transformation(y = Y_smp, l = optpar, inv=F)
  Y_smp       = tmp$y
  par_m       = tmp$m

  # Schätzen des Modells
  nlm_smp     = nlme::lme(Y_smp ~ -1 + X_smp,random = ~1|as.factor(smp_domains),method="REML")

  beta        = nlme::fixed.effects(nlm_smp)
  u_est_d     = nlme::random.effects(nlm_smp)
  sigmae2est  = nlm_smp$sigma^2
  sigmau2est  = as.numeric(nlme::VarCorr(nlm_smp)[1,1])
  gamma_est_d = sigmau2est / (sigmau2est + sigmae2est / n_smp)
  alpha_est_d = (sigmau2est * (1 - gamma_est_d) + sigmae2est)/2

  #m_di <- NULL
  w_est_ohnebc_dri_trafo <- NULL

  for(j in 1: length(area_size)){
    if(n_smp[j] != 0){
      w_est_ohnebc_dri_trafo[j] <- x_mean_d[j,] %*% beta + u_est_d[row.names(u_est_d) == j, ] + alpha_est_d[j]
    }else{
      w_est_ohnebc_dri_trafo[j] <- x_mean_d[j,] %*% beta + alpha_est_d[j]
    }
  }

  tau_est_ohnebc_d <- NULL

  w_est_ohnebc_dr <- transformation(optpar, w_est_ohnebc_dri_trafo, inv = TRUE, m = 0)$y

  for(j in 1: length(area_size)){
    w_ds_ <- data_smp[data_smp[, domains] == (1:length(area_size))[j],]
    w_ds  <- as.matrix(w_ds_[paste(formel[2])])
    w_est_ohnebc_dr_total <- w_est_ohnebc_dr[j] * (area_size[j] - length(w_ds))

    tau_est_ohnebc_d[j] <- 1 / (area_size[j]) * (sum(w_ds) + w_est_ohnebc_dr_total)
  }

  return(tau_est_ohnebc_d)
}
