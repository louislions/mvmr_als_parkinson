packages = c('MRPRESSO', 'gsmr2')
invisible(lapply(packages, library, character.only = TRUE))



compute_gsmr <- function(dat) {
    # Get LD Matrix
    ld_matrix <- ld_matrix(dat$SNP,
                            with_alleles = FALSE,
                            pop = "EUR",
                            plink_bin = genetics.binaRies::get_plink_binary(),
                            bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR"))

    gsmr_res <- gsmr(bzx = dat$beta.exposure, bzx_se = dat$se.exposure, bzx_pval = dat$pval.exposure, 
                    bzy = dat$beta.outcome, bzy_se = dat$se.outcome, bzy_pval = dat$pval.outcome,
                    ldrho = ld_matrix, snpid = colnames(ld_matrix), n_ref = dat$samplesize.outcome[1], heidi_outlier_flag=T,
                    gwas_thresh=5e-8, single_snp_heidi_thresh=0.01, multi_snps_heidi_thresh=0.01, 
                    nsnps_thresh=10, ld_r2_thresh=0.05, ld_fdr_thresh=0.05, gsmr2_beta=1)
                    
    return(data.frame(b = gsmr_res$bxy, se = gsmr_res$bxy_se, pval = gsmr_res$bxy_pval))
}

compute_mr_presso <- function(dat, SignifThreshold = 0.05) {
    presso_res <- mr_presso(data = dat,
                            BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", SdExposure = "se.exposure",
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                            NbDistribution = nrow(dat)/SignifThreshold,  
                            SignifThreshold = SignifThreshold, seed = global_seed)
    presso_res <- presso_res[[1]][2,]
    presso_res$method <- c('MR PRESSO')

    presso_res <- presso_res %>%
                    select(-Exposure, -'MR Analysis', -'T-stat') %>%
                    rename(b = 'Causal Estimate',
                        se = Sd,
                        pval = 'P-value')
    return(presso_res)
}

run_all_mr <- function(dat, nboot = 10000, over_dispersion = FALSE, 
                       run_mr_presso = FALSE, 
                       param_mr = TwoSampleMR::default_parameters(),
                       method_list = full_method_list)
{   
    # MR Egger analysis
    res_egger_int <- NULL
    if ("mr_egger_regression" %in% method_list) {
        res_egger <- mr_egger_regression(b_exp = dat$beta.exposure,b_out = dat$beta.outcome, 
                                         se_exp = dat$se.exposure, se_out = dat$se.outcome,
                                         parameters = param_mr)
        res0 <- c(id.exposure = dat$id.exposure[1],
                    id.outcome = dat$id.outcome[1],
                    outcome = dat$outcome[1],
                    exposure = dat$exposure[1])

        res_egger_int <- rbind(
                            c(res0,
                            method = "MR Egger",
                            nsnp = res_egger$nsnp,
                            b = res_egger$b,
                            se = res_egger$se,
                            pval = res_egger$pval),
                            c(res0,
                            method = "MR Egger (intercept)",
                            nsnp = res_egger$nsnp,
                            b = res_egger$b_i,
                            se = res_egger$se_i,
                            pval = res_egger$pval_i))

        method_list <- method_list[! method_list %in% c("mr_egger_regression")]
    }
    res <- rbind(TwoSampleMR::mr(dat = dat, parameters = param_mr, method_list = method_list),
                 res_egger_int)

    # MR Presso analysis
    if (run_mr_presso) {
        presso_res <- compute_mr_presso(dat)
        diff_col <- setdiff(colnames(res), colnames(presso_res))
        presso_res[diff_col] <- res[1, diff_col]
        res <- merge(res, presso_res, all = TRUE)
    }

    # Add GSMR method
    # Essayer de modif manuellement ld_matrix pr la rendre compatible avec gsmr avec duplicates
    #res <- add_result_to_df(compute_gsmr(dat), 'GSMR', res)

    res[, c("nsnp", "b", "se", "pval")] <- lapply(res[, c("nsnp", "b", "se", "pval")], as.numeric)
    return(res)
}

mvmr_robust <- function(dat, pval_thresh = 5e-8){
    bx <- dat$exposure_beta
    by <- dat$outcome_beta
    seby <- dat$outcome_se

    robmod = robust::lmRob(by ~ bx - 1, weights = seby^-2)
    coefficients = summary(robmod)$coef[, 1]
    se = summary(robmod)$coef[, 2] / min(summary(robmod)$sigma, 1)
    # nsnps <- (dat$exposure_pval < pval_thresh) %>%
    #     apply(MARGIN=2, FUN=sum) %>% as.vector()
    nsnps <- nrow(dat$exposure_pval)

    res_robust <- data.frame(id.exposure = colnames(dat$exposure_beta), 
                        id.outcome = dat$outname$id.outcome, outcome=dat$outname$outcome, 
                        nsnp = nsnps, b = coefficients, se = se,
                        pval = 2 * stats::pnorm(abs(coefficients)/se, lower.tail = FALSE), 
                        stringsAsFactors = FALSE)
    return(merge(dat$expname, res_robust))
}

mvmr_egger <- function(dat) {
    res_egger <- mr_mvegger(object = mr_mvinput(bx = dat$exposure_beta, 
                                                bxse = dat$exposure_se,
                                                by = dat$outcome_beta,
                                                byse = dat$outcome_se),
                            orientate = 1,
                            correl = FALSE,
                            distribution = "normal",
                            alpha = 0.05)

    res_egger <- data.frame(id.exposure = colnames(dat$exposure_beta), 
            id.outcome = dat$outname$id.outcome, outcome=dat$outname$outcome, 
            nsnp = res_egger$SNPs, b = res_egger$Estimate, se = res_egger$StdError.Est,
            pval = res_egger$Pvalue.Est, Q = res_egger$Heter.Stat[1],
            Q_df = res_egger$SNPs - length(res_egger$Estimate) - 1, 
            Q_pval = res_egger$Heter.Stat[2],
            stringsAsFactors = FALSE) %>%
    merge(dat$expname, .)
    
    return(res_egger)
}