packages <- c("mr.raps")
invisible(lapply(packages, library, character.only = TRUE))


global_uv <- function(api_use_ex, api_use_out, exposure_id, outcome_id,
                      exposure_useful_cols = "", outcome_useful_cols = "",
                      snplist = NULL, action = 2, allow_palindromes = 1,
                      clump_r2 = 0.001, kb = 10000, pval = 5e-8,
                      nboot = 10000, over_dispersion = FALSE, run_mr_presso = FALSE,
                      paper_results = NULL,
                      method_list = full_method_list)
{
    # Import data
    dat <- data_import(api_use_ex, api_use_out, exposure_id, outcome_id, 
                       exposure_useful_cols, outcome_useful_cols, 
                       snplist, action = action, allow_palindromes = allow_palindromes, 
                       clump_r2 = clump_r2, kb = kb, pval = pval)

    # Save it into a temp file
    # dat_temp_file <- tempfile()
    # write.csv(dat, file = dat_temp_file, row.names = FALSE)
    # print(paste("Data saved in", dat_temp_file))

    # MR Methods
    res <- run_all_mr(dat = dat, nboot = nboot, 
                      over_dispersion = over_dispersion, 
                      run_mr_presso = run_mr_presso,
                      method_list = method_list)

    res <- plotting_and_sensitivity_analysis(res = res, dat = dat,
                                             api_use_ex = api_use_ex,
                                             param_mr = param_mr,
                                             method_list = method_list)
    print(compare_to_paper(res, paper_results, api_use_ex))

    return(list(data = dat, result = res))
}

global_mv <- function(run_all_uvmr, run_mvmr,
                      api_use, exposure_id_list, outcome_id, 
                      exposure_useful_cols_list = NULL, outcome_useful_cols = NULL, 
                      snplist = NULL, action = 2, allow_palindromes = 1, 
                      clump_r2 = 0.001, clump_kb = 10000, pval = 5e-8,
                      nboot = 10000, over_dispersion = FALSE, run_mr_presso = FALSE,
                      intercept = FALSE, instrument_specific = TRUE, plots = FALSE,
                      pop = "EUR", plink_bin = genetics.binaRies::get_plink_binary(),
                      bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR"))
{
    dat <- data_import_mv(api_use, exposure_id_list, outcome_id, 
                          exposure_useful_cols_list, outcome_useful_cols,
                          snplist, action, allow_palindromes,
                          clump_r2, clump_kb, pval,
                          pop, plink_bin, bfile)
    uv_dat = dat$uv_dat
    mv_dat = dat$mv_dat

    # Perform standard MR analysis
    set.seed(global_seed)
    param_mr <- TwoSampleMR::default_parameters()
    param_mr$nboot <- nboot
    param_mr$over.dispersion <- over_dispersion

    if (run_all_uvmr) {
        method_list <- c("mr_wald_ratio", "mr_egger_regression", "mr_simple_median", 
                        "mr_weighted_median", "mr_ivw", "mr_raps")

        res <- lapply(uv_dat, function(x) {
            run_all_mr(dat = x, nboot = nboot, over_dispersion = over_dispersion, 
                    run_mr_presso = run_mr_presso, param_mr = param_mr, method_list = method_list)
        })

        # Plotting and sensitivity analysis
        for(i in seq_along(res)) {
            res[[i]] <- plotting_and_sensitivity_analysis(res[[i]], uv_dat[[i]], 
                                                            api_use_ex = api_use,
                                                            param_mr = param_mr, 
                                                            method_list = method_list)
            print(compare_to_paper(res[[i]]))
        }
    }
    
    if (run_mvmr) {
        res_list <- lapply(uv_dat, function(x) {
            x <- filter(x, mr_keep == TRUE)
            res_ivw <- TwoSampleMR::mr_ivw(b_exp = x$beta.exposure, b_out = x$beta.outcome, 
                                        se_exp = x$se.exposure, se_out = x$se.outcome) %>% 
                generate_odds_ratios()
            res_ivw$id.exposure <- x$id.exposure[1]
            return(res_ivw)
        })
        res_bind <- bind_rows(res_list)
        res_bind$mv <- FALSE

        # ----------
        # Computation of MVMR on mv_dat
        res_mv <- TwoSampleMR::mv_multiple(mv_dat, intercept, instrument_specific, pval, plots)$result %>% 
                generate_odds_ratios()
        res_mv$mv <- TRUE

        res <- export_mvmr_plot(res_mv, res_bind, api_use)
    }
    return(res)
}