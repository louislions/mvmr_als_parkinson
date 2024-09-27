packages = c("stringr")
invisible(lapply(packages, library, character.only = TRUE))



print_col_names <- function(df) {
    cat(paste(shQuote(colnames(data.frame(df)), type = "cmd"), collapse = ", "))
}

save_plot <- function(p, path, filename) {
    # Create the directory if it doesn't exist
    dir.create(path, showWarnings = FALSE)
    filepath <- paste0(path, filename)
    ggsave(filename = filepath, plot = p, width = 10, height = 8, dpi = 300)
}

add_result_to_df <- function(res, method_name, df) {
    res[["method"]] <- method_name
    missing_columns <- setdiff(names(df), names(res))
    for (col in missing_columns) {
        res[[col]] <- df[[col]][nrow(df)]
    }
    rbind(df, res)
}

remove_ids <- function(df) {
    # Delete columns id.exposure, id.outcome, outcome and exposure
    df %>%
        subset(select = -c(id.exposure, id.outcome, outcome, exposure))
}

merge_on_method <- function(df, res) {
    merge(df, res, by = "method", suffixes = c("_true", "_replicate"), all = TRUE)
}

add_comparaison_to_true_results <- function(res, outcome_results) {
    # Results from the paper's supplementary (doi.org/10.1038/s41588-021-00973-1)
    res0 <- fread(outcome_results)
    res <- merge_on_method(res0, res) %>%
            remove_ids()
    
    # Select the columns with suffix _true and _replicate
    true_cols <- grep("_true$", colnames(res), value = TRUE)
    replicate_cols <- grep("_replicate$", colnames(res), value = TRUE)

    # Compute the ratio for each column
    res[, paste0(colnames(res)[grep("_true$", colnames(res))], "_ratio%")] <- round(abs((res[, ..true_cols] - res[, ..replicate_cols]) / res[, ..true_cols]) * 100, digits = 2)

    return(res)
}

compare_to_paper <- function(res, paper_results = NULL, api_use) {

    if (!is.null(paper_results)) {
        plot_paper_comparaison(res, paper_results, api_use)
        #res <- add_comparaison_to_true_results(res, paste0(data_dir, "/paper-results/", paper_results)) 
    }
    else {
        res <- res[,!grepl("exposure$|outcome$",names(res))] %>%
            subset(method != "Wald ratio") %>%
            subset(select = -c(b, se, Q, Q_df, Q_pval, 
                                intercept, intercept_se, intercept_pval, SNP)) %>%
            kable()
    }
}

print_mvmr_res <- function(res) {
    res %>%
        select(id.exposure, nsnp, or, or_lci95, or_uci95, pval, mv) %>%
        kable() %>%
        print()
}

export_results <- function(res, dir) {
    filepath <- paste0(dir, "/mr_results", device = ".csv")
    write.csv(generate_odds_ratios(res), filepath, row.names = FALSE)
}

# conversion function
S4_to_dataframe <- function(s4obj) {
  nms <- slotNames(s4obj)

  lst <- lapply(nms, function(nm) slot(s4obj, nm))
  as.data.frame(setNames(lst, nms))
}

add_pheno_id_in_file <- function(file_path, df = NULL, id = NULL, pheno = NULL) {
    # Read the file if df is not provided
    if (is.null(df)) df <- fread(file_path) 
    # Add id
    if (!is.null(id)) {
        df <- mutate(df, id = rep(id, n()))
    }
    # Add phenotype name
    if (!is.null(pheno)) {
        pheno <- paste0(pheno, " || id:", df$id[1])
        df <- mutate(df, Phenotype = rep(pheno, n()))
    }

    file_path_no_gz <- gsub(".gz$", "", file_path)
    # Write the modified dataframe to a TSV file
    moiR::write.tsv(table = df, file = file_path_no_gz)
    # Overwrite the original TSV file with the modified dataframe
    R.utils::gzip(filename = file_path_no_gz, overwrite = TRUE)
}

get_row_as_vector <- function(df, row) {
    df[row, ] %>% as.vector() %>% as.character()
}

download_file <- function(id, destfile) {
    tophits(id, clump = 0, pval = 1) %>%
        write.csv(destfile, row.names = FALSE)
}

comma_char_to_num <- function(x) as.numeric(gsub(",", "", x))

export_snp_list_from_df <- function(df) {
    exp_name <- df$exposure[1] %>%
        gsub("[|]", "", .) %>%
        gsub(" ", "_", .)
    file_path <- paste0(data_dir, '/iv-lists/', exp_name, '.csv')
    write.csv(x = unique(subset(df, method == "Wald ratio")$SNP), 
              file = file_path,
              row.names = FALSE)
}

compute_f_stat <- function(dat, res) {
    # Check if one or more sample size is missing
    if (!('samplesize.exposure' %in% colnames(dat)) | (sum(is.na(dat$samplesize.exposure)) > 0)) { 
        f_stat <- dat$beta.exposure^2 / dat$se.exposure^2 # Approximation without sample size
    }
    else {
        n_iv <- dim(dat)[1]
        pve <- dat$beta.exposure^2 / (dat$beta.exposure^2 + dat$se.exposure^2 * dat$samplesize.exposure)
        f_stat <- (pve * (dat$samplesize.exposure - 1 - n_iv)) / ((1 - pve))
    }
    # p0 <- 0.085 # overall prevalence
    # p1 <- 0.095 # proportion of cases
    # z0 <- dnorm(qnorm(p0))
    # correction_factor <- p0^2 * (1 - p0)^2 / (p1 * (1 - p1) * z0^2)
    # f_stat <- f_stat * correction_factor

    res$F_min <- min(f_stat)
    res$F_max <- max(f_stat)

    return(res)
}

#unweightedIsq
compute_I_squared <- function(dat, res) {
    res$I <- TwoSampleMR::Isq(dat$beta.exposure %>% abs(), dat$se.exposure)
    return(res)
}

export_genetic_variant_associations <- function(dat, api_use) {
    if (api_use) ext <- "_API"
    else ext <- "_local"

    # Create the output directory if it doesn't exist
    dir <- paste0(output_dir,
                    gsub("[ /]", "-", dat$outcome[1]), 
                    "/", 
                    gsub("[ /]", "-", dat$exposure[1]),
                    ext)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    filepath <- paste0(dir, "/genetic_variant_associations", device = ".xlsx")

    dat %>% 
        select(SNP, effect_allele.exposure, eaf.exposure, 
               beta.exposure, se.exposure, pval.exposure, 
               beta.outcome, se.outcome, pval.outcome) %>% 
        writexl::write_xlsx(filepath)
} 
