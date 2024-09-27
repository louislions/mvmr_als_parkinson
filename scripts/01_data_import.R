packages = c()
invisible(lapply(packages, library, character.only = TRUE))



import_format_local_data <- function(file, useful_cols, snplist = NULL, type = "exposure") {
    df <- fread(file) %>% data.frame()
    # Check if the data contains the sample size for each SNP
    if (length(useful_cols) > 11) {
        df <- format_data(df, 
                        type = type,
                        snps = snplist,
                        header = TRUE,
                        snp_col = useful_cols[1],
                        chr_col = useful_cols[2],
                        pos_col = useful_cols[3],
                        effect_allele_col = useful_cols[4],
                        other_allele_col = useful_cols[5],
                        eaf_col = useful_cols[6],
                        beta_col = useful_cols[7],
                        se_col = useful_cols[8],
                        pval_col = useful_cols[9],
                        id_col = useful_cols[10],
                        phenotype_col = useful_cols[11],
                        samplesize_col = useful_cols[12],
                        min_pval = 1e-200)}
    else {
        df <- format_data(df, 
                        type = type,
                        snps = snplist,
                        header = TRUE,
                        snp_col = useful_cols[1],
                        chr_col = useful_cols[2],
                        pos_col = useful_cols[3],
                        effect_allele_col = useful_cols[4],
                        other_allele_col = useful_cols[5],
                        eaf_col = useful_cols[6],
                        beta_col = useful_cols[7],
                        se_col = useful_cols[8],
                        pval_col = useful_cols[9],
                        id_col = useful_cols[10],
                        phenotype_col = useful_cols[11],
                        min_pval = 1e-200)}

    if (!is.numeric(df$'pos')) {
        df$'pos' <- comma_char_to_num(df$'pos')
    }
    return(df)
}

clump_snp_data <- function(df, p = 5e-8, kb = 10000, r2 = 0.001) {
    #panel <- ld_reflookup(rsid = df$SNP, pop = 'EUR', opengwas_jwt = get_opengwas_jwt())
    # 1. Select rsids that are present in the remote LD reference panel
    # 2. Perform LD clumping
    df %>%
    # filter(SNP %in% panel) %>%
        filter(pval.exposure < p) %>%
        rename(rsid = SNP) %>%
        rename(pval = pval.exposure) %>%
        ld_clump(clump_kb = kb,
                clump_r2 = r2,
                clump_p = p,
                pop = "EUR",
                plink_bin = genetics.binaRies::get_plink_binary(),
                bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR")) %>%
        rename(SNP = rsid) %>%
        rename(pval.exposure = pval)
}

data_import <- function(api_use_ex, api_use_out, exposure_id, outcome_id, 
                        exposure_useful_cols = "", outcome_useful_cols = "", 
                        snplist = NULL, action = 2, allow_palindromes = 1, 
                        clump_r2 = 0.001, kb = 10000, pval = 5e-8) 
{
    # Extract the instruments for exposure
    if (api_use_ex) { # Use API of TwoSampleMR & ieugwasr for extraction
        if (is.null(snplist)) {
            expd <- TwoSampleMR::extract_instruments(exposure_id,
                                                     p1 = pval,
                                                     clump = TRUE,
                                                     p2 = pval,
                                                     r2 = clump_r2,
                                                     kb = kb)
        }
        else {
            expd <- TwoSampleMR::extract_instruments(exposure_id,
                                                     p1 = pval,
                                                     clump = FALSE) %>%
                    subset(SNP %in% snplist)
        }
    }
    else {
        expd <- import_format_local_data(exposure_id, exposure_useful_cols, snplist = snplist, type = "exposure")
        if (is.null(snplist)) {expd <- clump_snp_data(expd, p = pval)}
    }
    
    # Extract those SNP effects for outcome
    if (api_use_out) {
        outd <- TwoSampleMR::extract_outcome_data(expd$SNP, 
                                                outcome_id, 
                                                proxies = TRUE,        #Look for LD tags?
                                                rsq = 0.8,             #Minimum LD rsq value
                                                align_alleles = 1,       
                                                palindromes = allow_palindromes,       #Allow palindromic SNPs
                                                maf_threshold = 0.3,   #MAF threshold to try to infer palindromic SNPs
                                                splitsize = 10000,
                                                proxy_splitsize = 500)
    }
    else {
        outd <- import_format_local_data(file = outcome_id, 
                                         useful_cols = outcome_useful_cols, 
                                         snplist = expd$SNP,
                                         type = "outcome")
    }

    dat <- TwoSampleMR::harmonise_data(expd,
                                       outd,
                                       action = action)
    if (!is.null(snplist)) {dat$mr_keep <- TRUE}

    return(distinct(.data = dat, SNP, .keep_all = TRUE)) # Remove duplicates SNPs
}

mv_extract_data_locally <- function(exposure_id_list, outcome_id, 
                                    exposure_useful_cols_list = NULL, outcome_useful_cols = NULL,
                                    snplist = NULL, action = 2,
                                    clump_r2 = 0.001, clump_kb = 10000, pval = 5e-8,
                                    pop = "EUR", plink_bin = genetics.binaRies::get_plink_binary(),
                                    bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR"),
                                    second_clumping = TRUE) 
{
    ## Import exposure data of all exposures
    # List containing all available ivs for each exposure
    l_full <- list()
    # Subset of l_full only containing relevant ivs
    l_inst <- list()
    for(i in seq_along(exposure_id_list))
    {   
        # First import data as outcome-formated data
        l_full[[i]] <- import_format_local_data(file = exposure_id_list[i],
                                                useful_cols = exposure_useful_cols_list[[i]],
                                                type = "outcome")
        # Give name to exposures if not given
        if(l_full[[i]]$outcome[1] == "outcome") l_full[[i]]$outcome <- paste0("exposure", i)
        # Only keep snps under the significance threshold
        l_inst[[i]] <- subset(l_full[[i]], pval.outcome < pval)
        # Remove duplicated SNPs and convert to exposure format
        l_inst[[i]] <- subset(l_inst[[i]], !duplicated(SNP)) %>% convert_outcome_to_exposure()
        
        # If no predefined SNP list to extract given, clump the SNP according to given parameters
        if (is.null(snplist[[i]])) {
            l_inst[[i]] <- clump_data(l_inst[[i]], clump_p1=pval, clump_r2=clump_r2, 
                                    clump_kb=clump_kb, bfile=bfile, plink_bin=plink_bin, 
                                    pop=pop)
        }
        # Or extract predefined SNP list if given
        else {
            l_inst[[i]] <- subset(l_inst[[i]], SNP %in% snplist[[i]])
        }
        # Print number of SNP extracted
        message("Identified ", nrow(l_inst[[i]]), " hits for trait ", l_inst[[i]]$exposure[1])
    }
    # Combine data of all exposures
    exposure_dat <- dplyr::bind_rows(l_inst)
    # Get the unique exposure IDs of each exposure
    id_exposure <- unique(exposure_dat$id.exposure)

    ## Import outcome data on all the selected SNPs
    outd <- import_format_local_data(file = outcome_id, 
                                        useful_cols = outcome_useful_cols, 
                                        snplist = unique(exposure_dat$SNP),
                                        type = "outcome")

    ## Harmonize the data for each univariate analysis on each exposure
    uv_dat <- list()
    for (i in seq_along(exposure_id_list)) {  
        dat <- l_inst[[i]]
        uv_dat[[i]] <- TwoSampleMR::harmonise_data(dat, outd, action = action)
        
        # Force keep all the SNPs in MR analysis when predefined list of SNPs given
        if (!is.null(snplist[[i]])) {uv_dat[[i]]$mr_keep <- TRUE}
        # Remove duplicates SNPs
        uv_dat[[i]] <- distinct(.data = uv_dat[[i]], SNP, .keep_all = TRUE)
    }

    ## Processing and Harmonization of multivariate data
    # Preprocess
    exposure_dat$id.exposure <- 1
    exposure_dat <- exposure_dat[order(exposure_dat$pval.exposure, decreasing=FALSE), ]
    exposure_dat <- subset(exposure_dat, !duplicated(SNP))
    # Clump a second time the combined exposure data as data from a single exposure
    if (second_clumping) {
        exposure_dat <- clump_data(exposure_dat, clump_p1=pval, clump_r2=clump_r2, 
                            clump_kb=clump_kb, bfile=bfile, plink_bin=plink_bin, pop=pop)
    }
    message("Identified ", length(unique(exposure_dat$SNP)), " variants to include")

    # Extract the relevant and combined SNPs in all exposures and combine them
    d1 <- lapply(l_full, function(x) {
        subset(x, SNP %in% exposure_dat$SNP)
        }) %>% dplyr::bind_rows()

    # Make sure that we extracted data for each exposure
    stopifnot(length(unique(d1$id.outcome)) == length(unique(id_exposure)))
    d1 <- subset(d1, mr_keep.outcome)

    # Extract data from all but the 1st exposure
    d2 <- subset(d1, id.outcome != id_exposure[1])
    # Extract data from the first exposure and put it in exposure format
    d1 <- subset(d1, id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()

    # Harmonise against the first exposure id
    d <- harmonise_data(d1, d2, action=action)

    # Only keep SNPs that are present in each secondary exposure
    tab <- table(d$SNP)
    keepsnps <- names(tab)[tab == length(id_exposure)-1]
    d <- subset(d, SNP %in% keepsnps)

    # Re-align all exposure data
    dh1 <- subset(d, id.outcome == id.outcome[1],
                select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, mr_keep.exposure))
    dh2 <- subset(d, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome, mr_keep.outcome))
    names(dh2) <- gsub("outcome", "exposure", names(dh2))

    expd <- rbind(dh1, dh2) %>%
        rename(mr_keep = mr_keep.exposure)

    mv_dat <- TwoSampleMR::mv_harmonise_data(expd, outd, harmonise_strictness = action)

    return(list(uv_dat = uv_dat, mv_dat = mv_dat))
}

mv_extract_data_api <- function(exposure_id_list, outcome_id,
                                snplist = NULL, action = 2, allow_palindromes = 1,
                                clump_r2 = 0.001, clump_kb = 10000, pval = 5e-8,
                                pop = "EUR", plink_bin = genetics.binaRies::get_plink_binary(),
                                bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR"))
{
    if (is.null(snplist)) {
        expd <- TwoSampleMR::mv_extract_exposures(exposure_id_list,
                                                    clump_r2 = clump_r2,
                                                    clump_kb = clump_kb,
                                                    harmonise_strictness = action,
                                                    pval_threshold = pval,
                                                    pop = pop,
                                                    plink_bin = plink_bin,
                                                    bfile = bfile)

        outd <- extract_outcome_data(unique(expd$SNP), 
                                        outcome_id, 
                                        proxies = TRUE,                     #Look for LD tags?
                                        rsq = 0.8,                          #Minimum LD rsq value
                                        align_alleles = 1,       
                                        palindromes = allow_palindromes,    #Allow palindromic SNPs
                                        maf_threshold = 0.3,                #MAF threshold to try to infer palindromic SNPs
                                        splitsize = 10000,
                                        proxy_splitsize = 500)
        
        uv_dat <- list()
        for (i in seq_along(exposure_id_list)) {  
            expd_uv <- TwoSampleMR::extract_instruments(exposure_id_list[[i]],
                                                    p1 = pval,
                                                    clump = TRUE,
                                                    p2 = pval,
                                                    r2 = clump_r2)
            print(nrow(expd_uv))
            ####
            uv_dat[[i]] <- TwoSampleMR::harmonise_data(expd_uv, outd, action = action)
            if (!is.null(snplist[[i]])) {uv_dat[[i]]$mr_keep <- TRUE}
        }
        # Harmonise the exposure and outcome data
        mv_dat <- TwoSampleMR::mv_harmonise_data(expd, outd, harmonise_strictness = action)
    }
    #TODO
    else stop("API use is not yet implemented with a predetermined list of SNPs")

    return(list(mv_dat = mv_dat, uv_dat = uv_dat))
}

data_import_mv <- function(api_use, exposure_id_list, outcome_id, 
                           exposure_useful_cols_list = NULL, outcome_useful_cols = NULL,
                           snplist = NULL, action = 2, allow_palindromes = 1, 
                           clump_r2 = 0.001, clump_kb = 10000, pval = 5e-8,
                           pop = "EUR", plink_bin = genetics.binaRies::get_plink_binary(),
                           bfile = paste0(normalizePath(paste0(data_dir, "/LD_ref")), "/EUR"),
                           second_clumping = TRUE)
{
    # Get the number of exposures
    n <- length(exposure_id_list)
    if (!(n > 1)) warning("At least two exposures are needed for multivariate analysis")

    if (api_use) {
        dat <- mv_extract_data_api(exposure_id_list, outcome_id,
                                   snplist, action, allow_palindromes,
                                   clump_r2, clump_kb, pval,
                                   pop, plink_bin, bfile)
    }
    else {
        dat <- mv_extract_data_locally(exposure_id_list, outcome_id, 
                                       exposure_useful_cols_list, outcome_useful_cols,
                                       snplist, action,
                                       clump_r2, clump_kb, pval,
                                       pop, plink_bin, bfile, second_clumping)
    }
    return(dat)
}