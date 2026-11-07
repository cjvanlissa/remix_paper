library(data.table)
library(ggplot2)

outdir <- "simulation"
ignore <- c("reps", "rownum", "seed", "true_pi")
factors <- c('k', 'classes', 'vi', 'es', 'p_ratio')
enum <- c('aic', 'bic', 'icl', 'pmc', 'cv')
pars <- c('mu', 'tau', 'pi')
true_tau = 1

outfile <- file.path(outdir, "output.txt")

f <- list.files(outdir, pattern = "^sim_results_", full.names = TRUE)
f <- f[which.max(sapply(f, function(x) file.info(x)$mtime))]

dat <- as.data.table(readRDS(f))
setnames(dat, make.unique(names(dat)))
if ('rownum.1' %in% names(dat)) dat[, rownum.1 := NULL]
dat[, (ignore) := NULL]
#dat[, (grep("_(mc|cpp)$", names(dat), value = TRUE)) := NULL]
#data.table::setnames(dat, grep("_mx$", names(dat), value = TRUE), gsub("_mx$", "", grep("_mx$", names(dat), value = TRUE)))

# True parameter values
dat[, `:=`(
  true_mu = 0,
  true_tau = 1,
  true_pi = p_ratio / (p_ratio + classes - 1)
)]

# Bias
for (par in c("mu", "tau", "pi")) {
  dat[, paste0("bias_", par) :=
        get(paste0("res_", par)) - get(paste0("true_", par))]
}

# Confidence interval coverage and width
ci_methods <- c("mu", "tau_wald", "tau_profile", "pi")

for (method in ci_methods) {

  par <- sub("_.*", "", method)
  truth <- dat[[paste0("true_", par)]]

  lb <- dat[[paste0("ci_", method, "1")]]
  ub <- dat[[paste0("ci_", method, "2")]]

  dat[, paste0("cov_", method) :=
        lb <= truth & truth <= ub]

  dat[, paste0("width_", method) := ub - lb]
}

# Coverage for tau
prop.table(table(dat$cov_tau_wald))
prop.table(table(dat$cov_tau_profile))
dat[, (grep("wald", names(dat), value = TRUE)) := NULL]
setnames(dat, grep("profile", names(dat), value = TRUE), gsub("_profile", "", grep("profile", names(dat), value = TRUE)))

# Check convergence
conv_tab <- prop.table(table(dat$valid))

dat[, paste0("bias_abs_", pars) := lapply(.SD, abs), .SDcols = paste0("bias_", pars)]
dat[, paste0("acc_", enum) := lapply(.SD, function(x) x == classes), .SDcols = enum]
dat[, paste0("over_", enum) := lapply(.SD, function(x) x > classes), .SDcols = enum]
dat[, paste0("under_", enum) := lapply(.SD, function(x) x < classes), .SDcols = enum]

dat[, (factors) := lapply(.SD, function(v) ordered(v, levels = unique(v)[order(as.numeric(unique(v)))])), .SDcols = factors]


cols <- setdiff(names(dat), c(ignore, factors))

sum_fun <- function(x, vnam){
  x_obs <- x[!is.na(x)]

  data.table(
    outcome = vnam,
    n = 100 * length(x_obs) / length(x),
    mean = mean(x_obs),
    median = median(x_obs),
    sd = sd(x_obs),
    p10 = quantile(x_obs, 0.10),
    p90 = quantile(x_obs, 0.90)
  )
}
summary_stats <- rbindlist(lapply(cols, function(vnam) {
  x <- dat[[vnam]]
  sum_fun(x, vnam)
}))


sum_by_factor <- rbindlist(lapply(factors, function(fac) {

  rbindlist(lapply(cols, function(vnam) {
    result <- dat[, {
      x <- as.numeric(get(vnam))
      sum_fun(x, vnam)
    }, by = fac]
    setnames(result, fac, "level")
    result[, level := as.character(level)]
    result[, "outcome" := vnam]
    result[, "factor" := fac]
    result
  }), fill = TRUE)
}), fill = TRUE)


fwrite(summary_stats, file.path(outdir, 'performance_overall.csv'))
fwrite(sum_by_factor, file.path(outdir, 'performance_by_condition.csv'))

# Three plots: bias (signed + absolute), coverage, and enumeration
# Prepare plotting data

plot_theme <- function(data, filter, ribbons = FALSE, lines_by = NULL, panel = TRUE) {
  p <- ggplot(
    data,
    aes(
      x = level,
      y = mean,
      group = 1
    )
  )
  if(ribbons){
    p <- p + geom_ribbon(aes(ymin = p10, ymax = p90), alpha = .4)
  }
  if(!is.null(lines_by)){
    p <- p +
      geom_line(aes(linetype = .data[[lines_by]], group = .data[[lines_by]]), na.rm = TRUE) +
      geom_point(aes(shape = .data[[lines_by]], group = .data[[lines_by]]), na.rm = TRUE, size = 2)
  } else {
    p <- p +
      geom_line(na.rm = TRUE) +
      geom_point(na.rm = TRUE, size = 2)
  }
  if(panel){
    p <- p +
      facet_grid(
        outcome ~ factor,
        scales = "free",
        space = "free_x"
      ) +
      ggh4x::facetted_pos_scales(
        y = list(

          # Proportions
          outcome %in% c(
            "cov_mu", "cov_tau", "cov_pi",
            "acc_aic", "acc_bic", "acc_icl", "acc_pmc", "acc_cv",
            "over_aic", "over_bic", "over_icl", "over_pmc", "over_cv",
            "under_aic", "under_bic", "under_icl", "under_pmc", "under_cv"
          ) ~ scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.25)
          ),

          # Valid
          outcome == "valid" ~ scale_y_continuous(
            limits = c(0, 4),
            breaks = 0:4
          ),

          # Bias: symmetric around zero
          outcome %in% c(
            "bias_mu", "bias_tau", "bias_pi"
          ) ~ scale_y_continuous(
            limits = function(lims) {
              m <- max(abs(lims), na.rm = TRUE)
              c(-m, m)
            },
            breaks = scales::breaks_pretty(n = 5)
          ),

          # Absolute bias
          outcome %in% c(
            "bias_abs_mu", "bias_abs_tau", "bias_abs_pi"
          ) ~ scale_y_continuous(
            limits = c(0, NA),
            breaks = scales::breaks_pretty(n = 5)
          ),

          # Number of classes
          outcome %in% c(
            "aic", "bic", "icl", "pmc", "cv"
          ) ~ scale_y_continuous(
            limits = c(1, 4),
            breaks = 1:4
          ),

          # Minus 2 log-likelihood
          outcome == "m2ll" ~ scale_y_continuous(
            breaks = scales::breaks_pretty(n = 5)
          ),

          # Runtime
          outcome == "time" ~ scale_y_continuous(
            limits = c(0, NA),
            breaks = scales::breaks_pretty(n = 5)
          )
        )
      )
  }
  p <- p +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom"
    )

  # ggsave(file.path(outdir, paste0("plot_", make.names(selection), ".pdf")),
  #                  p, width = 15, height = 8)
  return(p)
}


df_plot <- copy(sum_by_factor)
df_plot <- df_plot[!level == "0", ]
df_plot$level <- ordered(df_plot$level, levels = levels(df_plot$level)[order(as.numeric(levels(df_plot$level)), decreasing = F)])
df_plot$level <- droplevels(df_plot$level)
# 1. Bias and absolute bias
#ordered_levels <- unique(unlist(factor_levels, use.names = FALSE))
df_bias <- copy(df_plot[grep("bias_(abs_)?(mu|ta|pi)", outcome), ])
df_bias[, type := ordered(grepl("abs", outcome), levels = c("FALSE", "TRUE"), labels = c("bias", "absolute bias"))]
df_bias[, outcome := ordered(gsub("bias_(abs_)?", "", outcome), levels = c("mu", "tau", "pi"))]
my_labs <- c(
  "k" = "K",
  "classes" = "C",
  "vi" = "v[i]",
  "es" = "Delta[mu]",
  "p_ratio" = "p[plain(ratio)]",
  "mu" = "mu",
  "tau" = "tau^2",
  "pi" = "pi"
)
p <- plot_theme(df_bias, ribbons = F, lines_by = "type") + geom_hline(yintercept = 0, linetype = 2)
p <- p + facet_grid(
  outcome ~ factor,
  scales = "free",
  space = "free_x",
  labeller = as_labeller(my_labs, default = label_parsed)
)
ggsave(file.path(outdir, paste0("plot_bias.pdf")), p, width = 29, height = 10, units = "cm")



# 2. Coverage
# Coverage variables and simulation factors
df_cov <- copy(df_plot[grep("^cov_", outcome), ])
df_cov[, outcome := ordered(outcome, levels = paste0("cov_", c("mu", "tau", "pi")), labels = c("mu", "tau", "pi"))]
names(my_labs) <- gsub("^bias", "cov", names(my_labs))

p <- plot_theme(df_cov, ribbons = F) + geom_hline(yintercept = .95, linetype = 2)
p <- p + facet_grid(
  outcome ~ factor,
  scales = "free",
  space = "free_x",
  labeller = as_labeller(my_labs, default = label_parsed)
)
ggsave(file.path(outdir, paste0("plot_cov.pdf")), p, width = 29, height = 10, units = "cm")
breaks <- unique(df_cov$outcome)
labs <- parse(text = unname(my_labs[as.character(breaks)]))

p <- plot_theme(df_cov, ribbons = F, lines_by = "outcome", panel = F) +
  geom_hline(yintercept = .95, linetype = 2) +
  scale_y_continuous(limits = c(0,1))+
  scale_linetype_discrete(
    name = "Parameter",
    breaks = breaks,
    labels = labs
  ) +
  scale_shape_discrete(
    name = "Parameter",
    breaks = breaks,
    labels = labs
  ) +
  guides(
    linetype = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  )+
  facet_grid(
  ~ factor,
  scales = "free",
  space = "free_x",
  labeller = as_labeller(my_labs, default = label_parsed))

ggsave(file.path(outdir, paste0("plot_cov2.pdf")), p, width = 29, height = 10, units = "cm")




# 3. Class enumeration

df_enum <- copy(df_plot[grep("^(acc|under|over)", outcome), ])
df_enum[, algorithm := ordered(gsub(".+?_", "", outcome), levels = c("aic", "bic", "icl", "pmc", "cv"), labels = c("AIC", "BIC", "ICL", "PMC", "CV"))]
df_enum[, outcome := gsub("_.*$", "", outcome)]
my_labs <- levels(df_enum$algorithm)
names(my_labs) <- my_labs
my_labs <- c(my_labs, c("acc" = "Accurate", "over" = "Overfit", "under" = "Underfit",
                        "k" = "K",
                        "classes" = "C",
                        "vi" = "v[i]",
                        "es" = "Delta[mu]",
                        "p_ratio" = "p[plain(ratio)]"))

p <- plot_theme(df_enum, ribbons = F, lines_by = "outcome", panel = F) +
  geom_hline(yintercept = .25, linetype = 2) +
  scale_linetype_discrete(
    name = "Outcome",
    labels = function(x) parse(text = unname(my_labs[x]))
  ) +
  scale_shape_discrete(
    name = "Outcome",
    labels = function(x) parse(text = unname(my_labs[x]))
  )+
  facet_grid(
    algorithm ~ factor,
    scales = "free",
    space = "free_x",
    labeller = as_labeller(my_labs, default = label_parsed))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = .25))

ggsave(file.path(outdir, paste0("plot_enum.pdf")), p, width = 29, height = 17, units = "cm")

# Descriptive partial eta-squared. Treat simulation factors as categorical.
# Sequential sums of squares depend on term order if the design is unbalanced.

# Select outcomes, excluding identifiers and factors
outcomes <- setdiff(names(dat), c("m2ll", ignore, factors))


partial_eta <- function(fit) {
  a <- anova(fit)
  ss <- a[["Sum Sq"]]
  names(ss) <- rownames(a)

  residual_ss <- ss["Residuals"]
  terms <- setdiff(names(ss), "Residuals")

  data.table(
    factor = terms,
    partial_eta2 = as.numeric(
      ss[terms] / (ss[terms] + residual_ss)
    )
  )
}

effects <- list()

for (outcome in outcomes) {

  d <- copy(dat[
    !is.na(get(outcome)),
    c(factors, outcome),
    with = FALSE
  ])

  # Convert logical outcomes to numeric
  if (is.logical(d[[outcome]])) {
    set(d, j = outcome, value = as.numeric(d[[outcome]]))
  }

  # Skip nonnumeric outcomes
  if (!is.numeric(d[[outcome]])) next

  # Remove nonfinite observations
  d <- d[is.finite(get(outcome))]

  #if (nrow(d) < 10L ||
  #    uniqueN(d[[outcome]]) < 2L) next

  # Treat simulation conditions as categorical
  d[, (factors) := lapply(.SD, factor), .SDcols = factors]

  # Skip factors with fewer than two observed levels
  valid_fac <- factors[
    vapply(d[, ..factors], uniqueN, integer(1)) > 1L
  ]

  if (!length(valid_fac)) next

  # Fit additive model
  fit <- try(
    lm(
      reformulate(valid_fac, response = outcome),
      data = d
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) next

  z <- partial_eta(fit)

  z[, `:=`(
    outcome = outcome,
    n = nobs(fit)
  )]

  effects[[length(effects) + 1L]] <- z
}

# Combine results
effect_sizes <- rbindlist(effects, fill = TRUE)

setcolorder(
  effect_sizes,
  c("outcome", "factor", "partial_eta2", "n")
)

effect_sizes <- dcast(
  effect_sizes,
  outcome ~ factor,
  value.var = "partial_eta2"
)

fwrite(
  effect_sizes,
  file.path(outdir, "factor_effect_sizes.csv")
)
