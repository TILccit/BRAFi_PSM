#’ Run propensity score matching & plot survival for pairwise multiple cohort comparisons
#’
#’ @description
#’ Takes a data frame with a cohort‐defining column and runs
#’ nearest‐neighbor propensity score matching (caliper = 0.2, standardized)
#’ for each pair you supply.  It then returns:
#’ 1) balance tables, 2) love plots, and 3) Kaplan–Meier curves annotated
#’    with the hazard ratio (HR) & 95% CI, optionally truncating at a
#’    specified follow‐up time.
#’
#’ @param dataDF A `data.frame` containing at least:
#’   * your cohort column (e.g. `Cohort` or `Tumor`)  
#’   * `OS_months`, `Dead`            — for overall survival  
#’   * `PFS_months`, `Progressed`     — for progression‐free survival  
#’   * plus whatever covariates your `matching_formula` uses  
#’
#’ @param comparisons A `list` of character vectors (length 2).  Each element
#’   names the two cohorts to compare, e.g.
#’   `list(c("Melanoma","Renal cancer"), c("Lung cancer","Melanoma"))`.  
#’
#’ @param covariate_to_split A single variable name (string) or simple RHS
#’   formula (e.g. `"Tumor"`, `"Treatment"`) used to define the strata in
#’   the KM curves and Cox model.  
#’
#’ @param endpoint Character, either `"OS"` (overall survival) or `"PFS"`
#’   (progression‐free survival).  Default: `"OS"`.  
#’
#’ @param truncate_month Optional numeric.  If non‐NULL, all times > this
#’   value are set to `truncate_month` and status to 0 (administrative
#’   censoring), and both curves & HRs are computed on the truncated data.
#’   Default: `NULL` (no truncation).  
#’
#’ @return A list with three elements:
#’  * `balance_tables`: list of before/after matching tables  
#’  * `love_plots`    : list of `ggplot` love plots  
#’  * `surv_plots`    : list of `ggsurvfit` survival plots (w/ HR & 95% CI)  
#’
#’ @examples
#’ \dontrun{
#’ matching_formula <- Cohort ~ Age + Sex + Stage
#’ comps <- list(c("Melanoma","Renal cancer"), c("Lung cancer","Melanoma"))
#’ # no truncation, OS
#’ out1 <- run_psm_and_survival(myData, comps, "Cohort", "OS")
#’ # truncate at 24 months, PFS
#’ out2 <- run_psm_and_survival(myData, comps, "Cohort", "PFS", truncate_month = 24)
#’ print(out2$surv_plots[[1]])
#’ }
#’ @export
run_psm_and_survival <- function(dataDF,
                                 comparisons,
                                 covariate_to_split,
                                 endpoint = c("OS", "PFS"),
                                 truncate_month = NULL) {
  
  endpoint <- match.arg(endpoint)
  stopifnot(require(MatchIt), require(cobalt),
            require(survival), require(ggsurvfit),
            require(ggplot2), require(dplyr))
  
  # 1) PSM & love plots
  balance_tbls <- vector("list", length(comparisons))
  love_plots   <- vector("list", length(comparisons))
  matched_dfs  <- vector("list", length(comparisons))
  
  for (i in seq_along(comparisons)) {
    pair   <- comparisons[[i]]
    df_sub <- subset(dataDF,
                     dataDF[[covariate_to_split]] %in% pair
    )
    
    fit_i <- matchit(
      matching_formula,
      data        = df_sub,
      method      = "nearest",
      distance    = "logit",
      ratio       = 1,
      caliper     = 0.2,
      std.caliper = TRUE
    )
    
    balance_tbls[[i]] <- bal.tab(fit_i, un = TRUE)
    
    love_plots[[i]] <- love.plot(x = fit_i, abs = TRUE) +
      ggtitle(sprintf("SMD Before & After Matching: %s vs %s",
                      pair[1], pair[2])) +
      geom_vline(xintercept = 0.2,
                 linetype   = "dashed",
                 color      = "black",
                 size       = 1)
    
    matched_dfs[[i]] <- match.data(fit_i)
  }
  
  # 2) Survival plots with Cox‐HR & optional truncation
  surv_plots <- vector("list", length(matched_dfs))
  
  for (i in seq_along(matched_dfs)) {
    md   <- matched_dfs[[i]]
    pair <- comparisons[[i]]
    
    # pick the right columns + title
    if (endpoint == "OS") {
      timevar    <- "OS_months";   statusvar <- "Dead"
      main_title <- "Propensity-Matched Overall Survival"
    } else {
      timevar    <- "PFS_months"; statusvar <- "Progressed"
      main_title <- "Propensity-Matched PFS"
    }
    
    # administrative censoring if requested
    if (!is.null(truncate_month)) {
      orig_t <- md[[timevar]]
      md[[timevar]]   <- pmin(orig_t, truncate_month)
      md[[statusvar]] <- ifelse(orig_t > truncate_month, 0, md[[statusvar]])
    }
    
    # build the Surv‐formula
    surv_formula <- as.formula(
      paste0("Surv(", timevar, ", ", statusvar, ") ~ ",
             covariate_to_split)
    )
    
    # Kaplan–Meier fit
    fit_km <- survfit(surv_formula, data = md)
    names(fit_km$strata) <- sub(
      paste0(covariate_to_split, "="), "",
      names(fit_km$strata)
    )
    
    # Cox model for HR & CI
    fit_cox <- coxph(surv_formula, data = md)
    s_cox   <- summary(fit_cox)
    hr      <- s_cox$coefficients[, "exp(coef)"]
    ci_lo   <- s_cox$conf.int[, "lower .95"]
    ci_hi   <- s_cox$conf.int[, "upper .95"]
    hr_label <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)",
                        hr[1], ci_lo[1], ci_hi[1])
    
    # determine the x‐axis max
    x_max <- if (is.null(truncate_month)) {
      max(md[[timevar]], na.rm = TRUE)
    } else truncate_month
    
    # build the plot
    p <- ggsurvfit(fit_km, size = 1.5) +
      add_censor_mark() +
      add_risktable(
        size = 5,
        theme = theme_risktable_default(
          axis.text.y.size = 15,
          plot.title.size  = 15
        ),
        risktable_stats = "{n.risk} ({cum.event})"
      ) +
      labs(
        x        = "Months",
        y        = "Survival Probability",
        title    = main_title,
        subtitle = sprintf("%s vs %s", pair[1], pair[2])
      ) +
      scale_x_continuous(
        breaks = seq(0, x_max, by = 12),
        limits = c(0, x_max)
      ) +
      scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1)) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5, size = 18),
        plot.subtitle   = element_text(hjust = 0.5, size = 15),
        axis.title.x    = element_text(size = 20),
        axis.title.y    = element_text(size = 20),
        axis.text.x     = element_text(size = 15),
        axis.text.y     = element_text(size = 15),
        legend.position = "bottom",
        legend.direction= "horizontal",
        legend.text     = element_text(size = 18),
        legend.key.size = unit(10, "bigpts"),
        legend.title    = element_blank()
      ) +
      annotate(
        "text",
        x     = x_max * 0.1,
        y     = 0.5,
        label = hr_label,
        size  = 6,
        hjust = 0
      )
    
    surv_plots[[i]] <- p
  }
  
  list(
    balance_tables = balance_tbls,
    love_plots      = love_plots,
    surv_plots      = surv_plots,
    matched_dfs = matched_dfs
  )
}
