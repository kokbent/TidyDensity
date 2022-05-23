# TidyDensity (development version)

## Breaking Changes
None

## New Features
1. Fix #181 - Add functions `color_blind()` `td_scale_fill_colorblind()` and 
`td_scale_color_colorblind()`
2. Fix #187 - Add functions `ci_lo()` and `ci_hi()`
3. Fix #189 - Add function `tidy_bootstrap()`

## Minor Fixes and Improvements
1. Fix #176 - Update `_autoplot` functions to include cumulative mean MCMC chart
by taking advantage of the `.num_sims` parameter of `tidy_` distribution
functions.
2. Fix #184 - Update `tidy_empirical()` to add a parameter of `.distribution_type`
3. Fix #183 - `tidy_empirical()` is now again plotted by `_autoplot` functions.
4. Fix #188 - Add the `.num_sims` parameter to `tidy_empirical()`

# TidyDensity 1.1.0

## Breaking Changes
None

## New Features
1. Fix #119 - Add function `tidy_four_autoplot()` - This will auto plot the density,
qq, quantile and probability plots to a single graph.
2. Fix #125 - Add function `util_weibull_param_estimate()`
3. Fix #126 - Add function `util_uniform_param_estimate()`
4. Fix #127 - Add function `util_cauchy_param_estimate()`
5. Fix #130 - Add function `tidy_t()` - Also add to plotting functions.
6. Fix #151 - Add function `tidy_mixture_density()`
7. Fix #150 - Add function `util_geometric_stats_tbl()`
8. Fix #149 - Add function `util_hypergeometric_stats_tbl()`
9. Fix #148 - Add function `util_logistic_stats_tbl()`
10. Fix #147 - Add function `util_lognormal_stats_tbl()`
11. Fix #146 - Add function `util_negative_binomial_stats_tbl()`
12. Fix #145 - Add function `util_normal_stats_tbl()`
13. Fix #144 - Add function `util_pareto_stats_tbl()`
14. Fix #143 - Add function `util_poisson_stats_tbl()`
15. Fix #142 - Add function `util_uniform_stats_tbl()`
16. Fix #141 - Add function `util_cauchy_stats_tbl()`
17. Fix #140 - Add function `util_t_stats_tbl()`
18. Fix #139 - Add function `util_f_stats_tbl()`
19. Fix #138 - Add function `util_chisquare_stats_tbl()`
20. Fix #137 - Add function `util_weibull_stats_tbl()`
21. Fix #136 - Add function `util_gamma_stats_tbl()`
22. Fix #135 - Add function `util_exponential_stats_tbl()`
23. Fix #134 - Add function `util_binomial_stats_tbl()`
24. Fix #133 - Add function `util_beta_stats_tbl()`

## Minor Fixes and Improvements
1. Fix #110 - Bug fix, correct the `p` calculation in `tidy_poisson()` that will
now produce the correct probability chart from the auto plot functions.
2. Fix #112 - Bug fix, correct the `p` calculation in `tidy_hypergeometric()` that
will no produce the correct probability chart from the auto plot functions.
3. Fix #115 - Fix spelling in Quantile chart.
4. Fix #117 - Fix probability plot x axis label.
4. Fix #118 - Fix fill color on combined auto plot
5. Fix #122 - The `tidy_distribution_summary_tbl()` function did not take the 
output of `tidy_multi_single_dist()` 
6. Fix #166 - Change in all plotting functions `ggplot2::xlim(0, max_dy)` to 
`ggplot2::ylim(0, max_dy)`
7. Fix #169 - Fix the computation of the `q` column 
8. Fix #170 - Fix the graphing of the quantile chart due to #169

# TidyDensity 1.0.1

## Breaking Changes
1. Fix #91 - Bug fix, change `tidy_gamma()` parameter of `.rate` to `.scale
Fix `tidy_autoplot_` functions to incorporate this change. Fix `util_gamma_param_estimate()`
to say `scale` instead of `rate` in the returned estimated parameters.

## New Features
None

## Minor Fixes and Improvements
1. Fix #90 - Make sure when `.geom_smooth` is set to TRUE that `ggplot2::xlim(0, max_dy)`
is set.
2. Fix #100 - `tidy_multi_single_dist()` failed on distribution with single parameter
like `tidy_poisson()`
3. Fix #96 - Enhance all `tidy_` distribution functions to add an attribute of 
either discrete or continuous that helps in the autoplot process.
4. Fix #97 - Enhance `tidy_autoplot()` to use histogram or lines for density plot
depending on if the distribution is discrete or continuous.
5. Fix #99 - Enhance `tidy_multi_dist_autoplot()` to use histogram or lines for
density plot depending on if the distribution is discrete or continuous.

# TidyDensity 1.0.0

## Breaking Changes
None

## New Features
1. Fix #27 - Add function `tidy_binomial()`
2. Fix #32 - Add function `tidy_geometric()`
3. Fix #33 - Add function `tidy_negative_binomial()`
4. Fix #34 - Add function `tidy_zero_truncated_poisson()`
5. Fix #35 - Add function `tidy_zero_truncated_geometric()`
6. Fix #36 - Add function `tidy_zero_truncated_binomial()`
7. Fix #37 - Add function `tidy_zero_truncated_negative_binomial()`
8. Fix #41 - Add function `tidy_pareto1()`
9. Fix #42 - Add function `tidy_pareto()`
10. Fix #43 - Add function `tidy_inverse_pareto()`
11. Fix #58 - Add function `tidy_random_walk()`
12. Fix #60 - Add function `tidy_random_walk_autoplot()`
13. Fix #47 - Add function `tidy_generalized_pareto()`
14. Fix #44 - Add function `tidy_paralogistic()`
15. Fix #38 - Add function `tidy_inverse_exponential()`
16. Fix #45 - Add function `tidy_inverse_gamma()`
17. Fix #46 - Add function `tidy_inverse_weibull()`
18. Fix #48 - Add function `tidy_burr()`
19. Fix #49 - Add function `tidy_inverse_burr()`
20. Fix #50 - Add function `tidy_inverse_normal()`
21. Fix #51 - Add function `tidy_generalized_beta()`
22. Fix #26 - Add function `tidy_multi_single_dist()`
23. Fix #62 - Add function `tidy_multi_dist_autoplot()`
24. Fix #66 - Add function `tidy_combine_distributions()`
25. Fix #69 - Add functions `tidy_kurtosis_vec()`, `tidy_skewness_vec()`, and
`tidy_range_statistic()`
26. Fix #75 - Add function `util_beta_param_estimate()`
27. Fix #76 - Add function `util_binomial_param_estimate()`
28. Fix #77 - Add function `util_exponential_param_estimate()`
29. Fix #78 - Add function `util_gamma_param_estimate()`
30. Fix #79 - Add function `util_geometric_param_estimate()`
31. Fix #80 - Add function `util_hypergeometric_param_estimate()`
32. Fix #81 - Add function `util_lognormal_param_estimate()`
33. Fix #89 - Add function `tidy_scale_zero_one_vec()`
34. Fix #87 - Add function `tidy_combined_autoplot()`
35. Fix #82 - Add function `util_logistic_param_estimate()`
36. Fix #83 - Add function `util_negative_binomial_param_estimate()`
37. Fix #84 - Add function `util_normal_param_estimate()`
38. Fix #85 - Add function `util_pareto_param_estimate()`
39. Fix #86 - Add function `util_poisson_param_estimate()`

## Fixes and Minor Improvements
1. Fix #30 - Move `crayon`, `rstudioapi`, and `cli` from Suggests to Imports due to `pillar`
no longer importing.
2. Fix #52 - Add parameter `.geom_rug` to `tidy_autoplot()` function
3. Fix #54 - Add parameter `.geom_point` to `tidy_autoplot()` function
4. Fix #53 - Add parameter `.geom_smooth` to `tidy_autoplot()` function
5. Fix #55 - Add parameter `.geom_jitter` to `tidy_autoplot()` function
6. Fix #57 - Fix `tidy_autoplot()` for when the distribution is `tidy_empirical()`
the legend argument would fail.
7. Fix #56 - Add attributes of .n and .num_sims (1L for now) to `tidy_empirical()`
8. Fix #61 - Update `_pkgdown.yml` file to update site.
9. Fix #67 - Add `param_grid`, `param_grid_txt`, and `dist_with_params` to the
attributes of all `tidy_` distribution functions.
10. Fix #70 - Add `...` as a grouping parameter to `tidy_distribution_summary_tbl()`
11. Fix #88 - Make the column `dist_type` a factor for `tidy_combine_distributions()`

# TidyDensity 0.0.1

## Breaking Changes
None

## New Features
1. Fix #1 - Add function `tidy_normal()`
2. Fix #4 - Add function `tidy_gamma()`
3. Fix #5 - Add function `tidy_beta()`
4. Fix #6 - Add function `tidy_poisson()`
5. Fix #2 - Add function `tidy_autoplot()`
6. Fix #11 - Add function `tidy_distribution_summary_tbl()`
7. Fix #10 - Add function `tidy_empirical()`
8. Fix #13 - Add function `tidy_uniform()`
9. Fix #14 - Add function `tidy_exponential()`
10. Fix #15 - Add function `tidy_logistic()`
11. Fix #16 - Add function `tidy_lognormal()`
12. Fix #17 - Add function `tidy_weibull()`
13. Fix #18 - Add function `tidy_chisquare()`
14. Fix #19 - Add function `tidy_cauchy()`
15. Fix #20 - Add function `tidy_hypergeometric()`
16. Fix #21 - Add function `tidy_f()`

## Minor Fixes and Improvements
None

# TidyDensity 0.0.0.9000

## Breaking Changes
None

## New Features
* Added a `NEWS.md` file to track changes to the package.

## Fixes and Minor Improvements
None
