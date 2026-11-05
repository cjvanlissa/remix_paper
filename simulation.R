# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'OpenMx')

lapply(packagenames, library, character.only = TRUE)

# Factorial design
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = c(2L:3L),
  vi = c(.001, .01, .1),
  es = c(2L:4L),
  p_ratio = c(1, 1.5, 2))
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = 1L,
  vi = c(.001, .01, .1),
  es = 0,
  p_ratio = 1)
summarydata_oneclass <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
summarydata <- rbind(summarydata, summarydata_oneclass)
summarydata$rownum <- 1:nrow(summarydata)

# Prepare seeds
RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
set.seed(824)
seeds <- vector("list", nrow(summarydata))
seeds[[1]] <- .Random.seed
for (i in 2:length(seeds)) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
summarydata$seed <- seeds

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
# summarydata <- summarydata[sample.int(nrow(summarydata), length(cl), replace = FALSE), ]
f <- list.files("paper/output", full.names = TRUE)
file.remove(f)
saveRDS(summarydata, file = "paper/summarydata.RData")
write(as.character(Sys.time()), file = "paper/output/time.log", append = TRUE)
max_classes <- max(summarydata$classes)+1L
foreach(rownum = 1:nrow(summarydata), .packages = packagenames, .errorhandling = "remove", .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]

  .Random.seed <- summarydata$seed[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    metamixture::sim_remix_fixed_tau(k = k, classes = classes, p_ratio = p_ratio, vi_range = c(.001, vi), es = es)
  })
  true_k <- length(unique(df$C))
  true_pi <- sum(df$C==1)/nrow(df)

  # OPENMX
  mx_try <- try({
    tstart <- Sys.time()
    res_mx <- lapply(
      seq_len(max_classes),
      metamixture:::remix,
      yi = df$yi,
      vi = df$vi,
      mu = "varying",
      tau2 = "equal",
      run = TRUE
    )
    time <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    conv <- sapply(
      res_mx,
      function(x) metamixture:::mx_convergence(x)$checks_pass
    )

    valid <- seq_len(sum(cumprod(conv)))

    if (length(valid) > 0L) {
      aic <- which.min(sapply(res_mx[valid], function(x) {
        tryCatch(
          x@output$Minus2LogLikelihood +
            2 * length(x@output$estimate),
          error = function(e) Inf
        )
      }))

      bic <- which.min(sapply(res_mx[valid], function(x) {
        tryCatch(
          x@output$Minus2LogLikelihood +
            log(x@data$numObs) * length(x@output$estimate),
          error = function(e) Inf
        )
      }))

      icl <- which.min(
        sapply(res_mx[valid], metamixture:::icl_mx)
      )

      if (length(valid) > 1L) {
        pmcs <- metamixture:::pmc_mx2(res_mx[valid])
        pmcs <- pmcs[1L:(length(valid)-1L), , drop = FALSE]

        pmc <- sum(cumprod(
          rowSums(sign(pmcs[, c("lb", "ub"), drop = FALSE])) == -2
        )) + 1L

        cv <- metamixture:::cv_enumeration(res_mx[valid])$best_class
      } else {
        pmc <- 1L
        cv <- 1L
      }

    } else {
      cv <- pmc <- aic <- bic <- icl <- NA
    }

    # Parameter estimates and coverage for the true model
    if (true_k %in% valid) {
      m2ll <- res_mx[[true_k]]@output$Minus2LogLikelihood
      cp <- metamixture:::class_prob_individual_remix(model = res_mx[[true_k]])
      ent <- metamixture:::entropy(cp)
      mus <- res_mx[[true_k]]@matrices$M$values[1, ]
      min_mu <- which.min(abs(mus))
      res_mu <- mus[min_mu]
      ci_mu <- res_mu +
        c(-1.96, 1.96) * res_mx[[true_k]]@output$standardErrors[paste0("mu", min_mu), 1]

      # Extract mixing proportion
      if(true_k > 1L){
        res_pi <- res_mx[[true_k]]$weights$result[1, min_mu]
        ci_pi <- res_pi +
            c(-1.96, 1.96) *
            OpenMx::mxSE("weights", res_mx[[true_k]], silent = TRUE)[1, min_mu]
      } else {
        res_pi <- ci_pi <- NA
      }
      res_tau <- res_mx[[true_k]]@matrices$S$values[1, 1]
      ci_tau_wald <- res_mx[[true_k]]$S$values[1,1]+ c(-1.96, 1.96)*res_mx[[true_k]]$output$standardErrors["tau2", 1]
      ci_tau_profile <- res_mx[[true_k]]$output$confidenceIntervals["tau2", c(1,3)]
    } else {
      m2ll <- NA
      ci_mu <- ci_tau_profile <- ci_tau_wald <- ci_pi <- NA
      res_mu <- res_tau <- res_pi <- NA
    }

  }, silent = TRUE)


  if (inherits(mx_try, "try-error")) {
    time <- m2ll <- NA
    cv <- aic <- bic <- icl <- pmc <- NA_real_
    ci_mu <- ci_tau_profile <- ci_tau_wald <- ci_pi <- NA
    res_mu <- res_tau <- res_pi <- NA
    valid <- NA
  }


  out <- data.frame(
    rownum = summarydata$rownum[rownum],
    true_pi = true_pi,

    valid = length(valid),
    m2ll = m2ll,
    ent = ent,
    aic = aic,
    bic = bic,
    icl = icl,
    pmc = pmc,
    cv = cv,
    time = time,
    res_mu = res_mu,
    res_tau = res_tau,
    res_pi = res_pi,
    ci_mu = ci_mu,
    ci_tau_wald = ci_tau_wald,
    ci_tau_profile = ci_tau_profile,
    ci_pi = ci_pi
  )
  saveRDS(out, file = paste0("paper/output/out_", rownum, ".rds"), compress = FALSE)
}
write(as.character(Sys.time()), file = "paper/output/time.log", append = TRUE)
#Close cluster
stopCluster(cl)
rm(cl)
## stop("End of simulation")

#save.image(paste0("tab", Sys.Date(), ".RData"))

# Read files --------------------------------------------------------------
library(data.table)

res <- as.data.table(readRDS(list.files(path = "paper", pattern = "summarydata", full.names = TRUE)))


f <- list.files("paper/output", pattern = "\\.rds", full.names = TRUE)
tab <- lapply(f, readRDS)
tab <- lapply(tab, function(x){
  c(unlist(x[1, c("rownum", "true_pi", "valid", "m2ll", "ent", "aic", "bic",
    "icl", "pmc", "cv", "time", "res_mu", "res_tau", "res_pi"), drop = T]),
    unlist(x[, c("ci_mu",
                 "ci_tau_wald", "ci_tau_profile", "ci_pi")]))
})
tab <- as.data.table(do.call(rbind, tab))
setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
setorderv(res, cols = c("rownum"), order=1L, na.last=FALSE)
if(!all(res$rownum == tab$rownum)){
  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$rownum)]
  stop()
}
merged <- cbind(res, tab)

saveRDS(merged, paste0("paper/sim_results_", Sys.Date(), ".RData"))

