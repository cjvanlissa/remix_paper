# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'OpenMx')

lapply(packagenames, library, character.only = TRUE)

outdir <- "paper/ci_profile"
dir.create(outdir)
f <- list.files(outdir)
file.remove(f)
# Read summarydata
summarydata <- readRDS(file = "paper/summarydata.RData")

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
summarydata <- summarydata[sample.int(nrow(summarydata), 10*length(cl), replace = FALSE), ]

max_classes <- max(summarydata$classes)+1L
ci_profile <- foreach(rownum = 1:nrow(summarydata), .packages = packagenames, .combine = rbind, .errorhandling = "remove") %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]

  .Random.seed <- summarydata$seed[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    metamixture::sim_remix_fixed_tau(k = k, classes = classes, p_ratio = p_ratio, vi_range = c(0, vi), es = es)
  })
  true_k <- length(unique(df$C))
  true_pi <- sum(df$C==1)/nrow(df)

  # OPENMX
  res_mx <- metamixture:::remix(
      yi = df$yi,
      vi = df$vi,
      classes = true_k,
      mu = "varying",
      tau2 = "equal",
      run = TRUE
    )
  c(res_mx$S$values[1,1]+ c(-1.96, 1.96)*res_mx$output$standardErrors["tau2", 1], res_mx$output$confidenceIntervals["tau2", c(1,3)])
}
#Close cluster
stopCluster(cl)
rm(cl)

colMeans(t(apply(ci_profile, 1, function(x){c(x[2]-x[1], x[4]-x[3])})), na.rm = TRUE)

colMeans(t(apply(ci_profile, 1, function(x){c(x[2] > 1 && x[1] < 1, x[4] > 1 && x[3] < 1)})), na.rm = TRUE)

## stop("End of simulation")
saveRDS(ci_profile, "outdir/ci_profile.RData")

