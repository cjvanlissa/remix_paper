  Handling between-study differences in effect sizes, or heterogeneity, is a key challenge in meta-analysis.
  The common assumption that heterogeneity is normally distributed implies an absence of systematic sources of   heterogeneity, or known and measured sources of heterogeneity that can be controlled using moderator analysis.
  These assumptions fall short when effect sizes originate from multiple unknown subpopulations - for example, due to different research traditions, unmeasured moderators, publication bias, and outliers.
  It has been argued that finite mixture models - also known as latent class analysis (LCA) - can capture these unobserved subgroups,
  but there is currently no generalized framework for estimating such models.
  We introduce open source software for estimating random-effects mixture, or "remix" meta-analyses.
  Implemented in `OpenMx`, these structural equation models can be fully customized and further extended.
  We conducted a simulation study to evaluate the `remix` algorithm's convergence, bias and confidence interval coverage,
  and the performance of several class enumeration metrics: the AIC, BIC, and ICL information criteria, predictive model comparison (PMC), and cross-validation.
  Results show that the method is fast and has low bias, but confidence interval coverage is below nominal levels,
  and class enumeration remains challenging - which are known limitations of LCA.
  Nevertheless, remix meta-analysis is the appropriate method when effect sizes originate from a mixture of unobserved subpopulations; a scenario which is arguably commonplace in meta-analysis.
