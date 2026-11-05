# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'foreach', 'metaplus')
lapply(packagenames, library, character.only = TRUE)

outdir <- "metaplus"
dir.create(outdir)
f <- list.files(outdir, full.names = TRUE)
file.remove(f)

# Factorial design
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = 2L,
  vi = c(.001, .01, .1),
  tau2_ratio = c(2, 4, 6),
  p_ratio = c(1, 1.5)
)
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
summarydata$rownum <- 1:nrow(summarydata)

# Prepare seeds
RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
set.seed(6994)
seeds <- vector("list", nrow(summarydata))
seeds[[1]] <- .Random.seed
for (i in 2:length(seeds)) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
summarydata$seed <- seeds
saveRDS(summarydata, file = file.path(outdir, "summarydata_metaplus.RData"))

sim_data <- function(k, classes, p_ratio, vi, r2) {

  # Class proportions
  ps <- c(p_ratio, rep(1, classes - 1L))
  ps <- ps / sum(ps)

  # Common mean
  mu <- rep(0, classes)

  # Class-specific heterogeneity variances
  # Class 2 has greater heterogeneity than class 1
  tau2 <- if (classes == 1L) {
    1
  } else {
    c(1, r2)
  }

  # True class memberships
  C <- sample(
    seq_len(classes),
    size = k,
    replace = TRUE,
    prob = ps
  )

  # Study-specific sampling variances
  vi_indiv <- runif(
    k,
    min = 0,
    max = vi * min(tau2)
  )

  # Generate observed effect sizes
  y <- rnorm(
    k,
    mean = mu[C],
    sd = sqrt(vi_indiv + tau2[C])
  )

  data.frame(
    yi = y,
    vi = vi_indiv,
    C = C,
    tau2 = tau2[C]
  )
}

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
summarydata <- summarydata[sample.int(nrow(summarydata), length(cl), replace = FALSE), ]
write(as.character(Sys.time()), file = file.path(outdir, "time.log"), append = TRUE)
foreach(rownum = 1:nrow(summarydata), .packages = c("metamixture", "OpenMx", "tidySEM", "metaplus"), .errorhandling = "remove", .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]
  assign(".Random.seed", summarydata$seed[[rownum]], envir = .GlobalEnv)
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    sim_data(k = k, classes = classes, p_ratio = p_ratio, vi = vi, r2 = tau2_ratio)
  })
  true_k <- summarydata$classes[rownum]
  true_pi <- summarydata$p_ratio[rownum] / (summarydata$p_ratio[rownum] + true_k - 1L)

  conv_mp <- FALSE
  time_mp <- m2ll_mp <- NA_real_
  res_mu_mp <- res_tau_mp <- res_pi_mp <- NA_real_
  ci_mu_mp <- ci_tau_mp <- ci_pi_mp <- c(NA_real_, NA_real_)

  # METAPLUS
  mp_try <- try({

    tstart <- Sys.time()
    fit_mp <- metaplus::metaplus(
      yi = df$yi,
      sei = sqrt(df$vi),
      random = "mixture",
      cores = 1,
      notrials = 20,
      justfit = FALSE
    )
    time_mp <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    m2ll_mp <- -2 * as.numeric(logLik(fit_mp))

    # Extract estimates and confidence intervals
    res_mu_mp <- fit_mp$results["muhat", "Est."]
    res_tau_mp <- fit_mp$results["tau2", "Est."]
    res_pi_mp <- (1-fit_mp$results["Outlier prob.", "Est."])

    ci_mu_mp <- fit_mp$results["muhat", c(2:3)]
    ci_tau_mp <- fit_mp$results["tau2", c(2:3)]
    ci_pi_mp <- 1 - rev(as.numeric(fit_mp$results["Outlier prob.", 2:3]))

    conv_mp <- is.finite(m2ll_mp)

  }, silent = TRUE)

  # OPENMX
  conv_mx <- FALSE
  time_mx <- m2ll_mx <- ent <- NA_real_
  res_mu_mx <- res_tau_mx <- res_pi_mx <- NA_real_
  ci_mu_mx <- c(NA_real_, NA_real_)
  ci_tau_mx_wald <- c(NA_real_, NA_real_)
  ci_tau_mx_profile <- c(NA_real_, NA_real_)
  ci_pi_mx <- c(NA_real_, NA_real_)

  mx_try <- try({
    tstart <- Sys.time()
    res_mx <- metamixture::remix(
      yi = df$yi,
      vi = df$vi,
      classes = 2L,
      mu = "equal",
      tau2 = "varying",
      run = TRUE
    )
    time_mx <- as.numeric(difftime(Sys.time(), tstart, units = "secs"))

    conv_mx <- tidySEM::mx_convergence(res_mx, check = c("convergence", "fit"))$checks_pass

    # Retain likelihood even for nonconverged fits
    m2ll_mx <- as.numeric(res_mx@output$Minus2LogLikelihood)

    if (conv_mx) {
      cp <- tidySEM:::extract_postprob(res_mx)
      ent <- metamixture:::entropy(cp)
      res_mu_mx <- res_mx$M$values[1, 1]
      ci_mu_mx <- res_mu_mx + c(-1.96, 1.96) * res_mx@output$standardErrors["mu", 1]

      # Identify smaller-variance class
      tau_mx <- diag(res_mx$S$values)
      small_tau_class <- which.min(tau_mx)
      res_tau_mx <- tau_mx[small_tau_class]

      res_pi_mx <- as.numeric(OpenMx::mxEvalByName("weights", res_mx))[small_tau_class]

      # Wald CI for smaller heterogeneity variance
      tau_name <- paste0("tau2_", small_tau_class)

      ci_tau_mx_wald <- res_tau_mx + c(-1.96, 1.96) * res_mx@output$standardErrors[tau_name, 1]
      ci_tau_mx_profile <- as.numeric(res_mx@output$confidenceIntervals[tau_name, c("lbound", "ubound")])
      ci_pi_mx <- res_pi_mx + c(-1.96, 1.96) * OpenMx::mxSE("weights", res_mx, silent = TRUE)[1, small_tau_class]
    }

  }, silent = TRUE)

  if (inherits(mx_try, "try-error")) {
    conv_mx <- FALSE
  }

  out <- data.frame(
    rownum = summarydata$rownum[rownum],

    # remix
    conv_mx = conv_mx,
    m2ll_mx = m2ll_mx,
    time_mx = time_mx,
    res_mu_mx = res_mu_mx,
    res_tau_mx = res_tau_mx,
    res_pi_mx = res_pi_mx,
    ent_mx = ent,
    ci_mu_mx1 = ci_mu_mx[1],
    ci_mu_mx2 = ci_mu_mx[2],
    ci_tau_mx_wald1 = ci_tau_mx_wald[1],
    ci_tau_mx_wald2 = ci_tau_mx_wald[2],
    ci_tau_mx_profile1 = ci_tau_mx_profile[1],
    ci_tau_mx_profile2 = ci_tau_mx_profile[2],
    ci_pi_mx1 = ci_pi_mx[1],
    ci_pi_mx2 = ci_pi_mx[2],

    # metaplus
    conv_mp = conv_mp,
    m2ll_mp = m2ll_mp,
    time_mp = time_mp,
    res_mu_mp = res_mu_mp,
    res_tau_mp = res_tau_mp,
    res_pi_mp = res_pi_mp,

    ci_mu_mp1 = ci_mu_mp[1],
    ci_mu_mp2 = ci_mu_mp[2],
    ci_tau_mp1 = ci_tau_mp[1],
    ci_tau_mp2 = ci_tau_mp[2],
    ci_pi_mp1 = ci_pi_mp[1],
    ci_pi_mp2 = ci_pi_mp[2]
  )
  stopifnot(nrow(out) == 1L)
  saveRDS(out, file = file.path(outdir, paste0("out_", rownum, ".rds")), compress = FALSE)
}
write(as.character(Sys.time()), file = file.path(outdir, "time.log"), append = TRUE)
#Close cluster
stopCluster(cl)
rm(cl)

# Read files --------------------------------------------------------------
library(data.table)

res <- as.data.table(readRDS(list.files(path = "metaplus", pattern = "summarydata_metaplus", full.names = TRUE)))
#res[, seed := NULL]

f <- list.files(outdir, pattern = "^out_[0-9]+\\.rds$", full.names = TRUE)
tab <- lapply(f, readRDS)

tab <- as.data.table(do.call(rbind, tab))

setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
setorderv(res, cols = c("rownum"), order=1L, na.last=FALSE)

merged <- merge(res, tab, by = "rownum", all.x = TRUE, sort = FALSE)

saveRDS(merged, file.path(outdir, paste0("sim_results_", Sys.Date(), ".RData")))

# Analysis

library(data.table)
dat <- readRDS(list.files("metaplus/", pattern = "sim_results", full.names = T))
dat <- as.data.table(dat)
# Remove duplicated column names
dat <- dat[, !duplicated(names(dat)), with = FALSE]

# True values for the current two-class simulation
dat[, `:=`(
  true_mu = 0,
  true_tau = 1,
  true_pi = p_ratio / (p_ratio + 1)
)]

# Bias for both estimation methods
for (method in c("mx", "mp")) {
  for (par in c("mu", "tau", "pi")) {

    dat[, paste0("bias_", par, "_", method) :=
          get(paste0("res_", par, "_", method)) -
          get(paste0("true_", par))]
  }
}

# CI coverage and width
ci_methods <- c(
  "mu_mx",
  "tau_mx_wald",
  "tau_mx_profile",
  "pi_mx",
  "mu_mp",
  "tau_mp",
  "pi_mp"
)

for (method in ci_methods) {

  par <- sub("_.*", "", method)
  truth <- dat[[paste0("true_", par)]]

  lb <- dat[[paste0("ci_", method, "1")]]
  ub <- dat[[paste0("ci_", method, "2")]]

  valid <- is.finite(lb) &
    is.finite(ub) &
    lb <= ub

  coverage <- rep(NA, nrow(dat))
  coverage[valid] <-
    lb[valid] <= truth[valid] &
    truth[valid] <= ub[valid]

  width <- rep(NA_real_, nrow(dat))
  width[valid] <- ub[valid] - lb[valid]

  dat[, paste0("cov_", method) := coverage]
  dat[, paste0("width_", method) := width]
}

prop.table(table(dat$cov_tau_mx_wald))
prop.table(table(dat$cov_tau_mx_profile))
dat[, (grep("wald", names(dat), value = TRUE)) := NULL]
setnames(dat, grep("profile", names(dat), value = TRUE), gsub("_profile", "", grep("profile", names(dat), value = TRUE)))
# 1. Optimization performance ----------------------------

optimization <- dat[, .(
  n = .N,

  # Number of converged models (1-4)
  mean_conv_mx = mean(conv_mx, na.rm = TRUE),
  mean_conv_mp = mean(conv_mp, na.rm = TRUE),

  # Number of paired likelihood comparisons
  n_ll = sum(is.finite(m2ll_mx) & is.finite(m2ll_mp)),

  # Which method obtained the lower minus2LL?
  better_mx = mean(m2ll_mx < m2ll_mp - 1e-4, na.rm = TRUE),
  better_mp = mean(m2ll_mp < m2ll_mx - 1e-4, na.rm = TRUE),
  equal_ll = mean(abs(m2ll_mx - m2ll_mp) <= 1e-4, na.rm = TRUE),

  # Difference: positive means MX has better fit
  #mean_ll_diff = mean(m2ll_mp - m2ll_mx, na.rm = TRUE),
  #median_ll_diff = median(m2ll_mp - m2ll_mx, na.rm = TRUE),

  # Runtime
  #mean_time_mx = mean(time_mx, na.rm = TRUE),
  #mean_time_mp = mean(time_mp, na.rm = TRUE),
  median_time_mx = median(time_mx, na.rm = TRUE),
  median_time_mp = median(time_mp, na.rm = TRUE)
)]

outdir <- "metaplus/results"
dir.create(outdir)

write.csv(optimization, file.path(outdir, "optimization.csv"), row.names = FALSE)

# Bias

# Reshape bias estimates into long format
bias_long <- melt(
  dat,
  measure.vars = patterns(
    mx = "^bias_(mu|tau|pi)_mx$",
    mp = "^bias_(mu|tau|pi)_mp$"
  ),
  variable.name = "parameter",
  value.name = c("mx", "mp")
)

bias_long[, parameter := c("mu", "tau2", "pi")[parameter]]

bias_table <- bias_long[, .(
  n_mx = sum(!is.na(mx)),
  n_mp = sum(!is.na(mp)),

  mean_bias_mx = mean(mx, na.rm = TRUE),
  mean_bias_mp = mean(mp, na.rm = TRUE),

  mean_abs_bias_mx = mean(abs(mx), na.rm = TRUE),
  mean_abs_bias_mp = mean(abs(mp), na.rm = TRUE),

  rmse_mx = sqrt(mean(mx^2, na.rm = TRUE)),
  rmse_mp = sqrt(mean(mp^2, na.rm = TRUE)),

  # Among paired, nonmissing estimates
  mx_more_accurate = mean(abs(mx) < abs(mp), na.rm = TRUE),
  mp_more_accurate = mean(abs(mp) < abs(mx), na.rm = TRUE)
), by = parameter]

write.csv(bias_table, file.path(outdir, "bias.csv"), row.names = FALSE)

# Coverage
coverage_long <- melt(
  dat,
  measure.vars = patterns(
    mx = "^cov_(mu|tau|pi)_mx$",
    mp = "^cov_(mu|tau|pi)_mp$"
  ),
  variable.name = "parameter",
  value.name = c("mx", "mp")
)

coverage_long[, parameter := c("mu", "tau2", "pi")[parameter]]

coverage_table <- coverage_long[, .(
  n_mx = sum(!is.na(mx)),
  n_mp = sum(!is.na(mp)),

  coverage_mx = if (all(is.na(mx))) NA_real_
  else mean(mx, na.rm = TRUE),

  coverage_mp = if (all(is.na(mp))) NA_real_
  else mean(mp, na.rm = TRUE)

  # Deviation from nominal 95% coverage
  #deviation_mx = if (all(is.na(mx))) NA_real_
  #else mean(mx, na.rm = TRUE) - 0.95,

  #deviation_mp = if (all(is.na(mp))) NA_real_
  #else mean(mp, na.rm = TRUE) - 0.95
), by = parameter]


write.csv(coverage_table, file.path(outdir, "coverage.csv"), row.names = FALSE)

summary_table <- data.table(
  method = c("remix", "metaplus"),

  better_ll = c(
    optimization$better_mx,
    optimization$better_mp
  ),

  runtime = c(
    median(dat$time_mx, na.rm = TRUE),
    median(dat$time_mp, na.rm = TRUE)
  ),

  MAE_mu = c(
    mean(abs(dat$bias_mu_mx), na.rm = TRUE),
    mean(abs(dat$bias_mu_mp), na.rm = TRUE)
  ),

  MAE_tau2 = c(
    mean(abs(dat$bias_tau_mx), na.rm = TRUE),
    mean(abs(dat$bias_tau_mp), na.rm = TRUE)
  ),

  MAE_pi = c(
    mean(abs(dat$bias_pi_mx), na.rm = TRUE),
    mean(abs(dat$bias_pi_mp), na.rm = TRUE)
  ),


  cov_mu = c(
    mean(dat$cov_mu_mx, na.rm = TRUE),
    mean(dat$cov_mu_mp, na.rm = TRUE)
  ),

  cov_tau2 = c(
    mean(dat$cov_tau_mx, na.rm = TRUE),
    NA_real_
  ),

  cov_pi = c(
    mean(dat$cov_pi_mx, na.rm = TRUE),
    NA_real_
  ),
  t(bias_table[, .SD, .SDcols = grep("more_accurate", names(bias_table))])
)
names(summary_table)[ncol(summary_table)-c(2:0)] <- paste0("Acc. ", c("mu", "tau2", "pi"))

write.csv(summary_table, file.path(outdir, "summary_table.csv"), row.names = FALSE)
