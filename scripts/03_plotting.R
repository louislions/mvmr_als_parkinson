packages = c('geni.plots')
invisible(lapply(packages, library, character.only = TRUE))


### Plotting exposure gwas ----------------------------------------------------------
plot_manhattan_exp <- function(dat, ex_dat) {
    n_ex <- nrow(ex_dat)
    manhattan_dat <- data.frame(chr = ex_dat$chr,
                                pos = ex_dat$position,
                                pvalue = ex_dat$p,
                                highlight = to_vec(for(i in 1:n_ex) if(ex_dat$rsid[i] %in% dat$SNP) 1 else 0),
                                highlight_shape = to_vec(for(i in 1:n_ex) if(ex_dat$rsid[i] %in% dat$SNP) 4 else 0),
                                label = to_vec(for(i in 1:n_ex) if(ex_dat$rsid[i] %in% dat$SNP) ex_dat$rsid[i] else "")) %>%
                                filter(pvalue != 0) # filter out null p-values making the plot fail

    p <- fig_manhattan(data = manhattan_dat,
                        block_thresh = 5e-8,
                        trunc = min(manhattan_dat$pvalue),
                        label_size = 2.5,
                        label_box = TRUE)

    dir <- paste0(output_dir, "/Manhattan_plots/")
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    set_panel_size(p = p, 
                    width=unit(max(nrow(dat) / 2, 20),"cm"), 
                    height=unit(-log(min(manhattan_dat$pvalue)) / 10, "cm"),
                    file = paste0(dir, gsub(" ", "-", dat$exposure[1]), device = ".png"))
    return
}

### Plotting outcome gwas ----------------------------------------------------------
plot_manhattan_snp <- function(snp, dat, out_dat) {
    n_ex <- nrow(out_dat)
    manhattan_dat <- data.frame(chr = out_dat$chr, pos = out_dat$position, pvalue = out_dat$p, 
                                highlight = to_vec(for(i in 1:n_ex) if((out_dat$rsid[i] %in% snp) | (out_dat$rsid[i] %in% dat$SNP)) 1 else 0),
                                highlight_shape = to_vec(for(i in 1:n_ex) if((out_dat$rsid[i] %in% snp) | (out_dat$rsid[i] %in% dat$SNP)) 4 else 0),
                                label = to_vec(for(i in 1:n_ex) if(out_dat$rsid[i] %in% snp) 
                                                                {
                                                                    if (out_dat$rsid[i] %in% dat$SNP) 
                                                                    {
                                                                        paste(out_dat$rsid[i], " (both)", sep = "") # in snp & dat
                                                                    } 
                                                                    else 
                                                                    {
                                                                        paste(out_dat$rsid[i], " (original)", sep = "") # in snp
                                                                    }
                                                                }   
                                                                else
                                                                {
                                                                    if (out_dat$rsid[i] %in% dat$SNP) 
                                                                    {
                                                                        paste(out_dat$rsid[i], " (new)", sep = "") # in dat
                                                                    } 
                                                                    else "" # neither
                                                                })) %>%
                            filter(pvalue != 0) # filter out null p-values making the plot fail

    p <- fig_manhattan(data = manhattan_dat,
                        block_thresh = 5e-8,
                        trunc = min(manhattan_dat$pvalue),
                        label_size = 2.5,
                        label_box = TRUE)

    w <- max(nrow(dat) / 2, 20)
    h <- -log(min(manhattan_dat$pvalue)) / 10

    dir <- paste0(output_dir, "/Manhattan_plots/")
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    set_panel_size(p = p, width=unit(w,"cm"), 
                    height=unit(h, "cm"),
                    file = paste0(dir, gsub(" ", "-", dat$exposure[1]), device = ".png"))
    return
}


### Other plots ---------------------------------------------------------------------
export_mr_scatter_plots <- function(res, dat, dir) {
    l <- mr_scatter_plot(res, dat)
    # Export each plot as a PNG image
    for (i in seq_along(l)) {
        set_panel_size(p = l[[i]], 
                    width=unit(10,"cm"), height=unit(10, "cm"),
                    file = paste0(dir, "/scatter_plot_", i, device = ".png"))
    }
}

export_mr_results_plot <- function(res, dir) {
    res <- generate_odds_ratios(res)
    plot <- ggplot(res, aes(x = or, y = method)) +
        geom_point(shape = 15, size = 4, color = 'red') +  # Fixed size for points
        geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.1, color = 'red') +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
        labs(x = paste("Odds ratio (95% CI) for ", 
                        res$exposure[1], 
                        " on ",
                        res$outcome[1]), 
                y = "",
                title = "Methods") +
        #theme_bw()
        theme_minimal() +
        theme(plot.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(), # Remove grid
                panel.grid.minor = element_blank(), # Remove minor grid
                panel.border = element_blank(), # Remove border
                axis.line.x = element_line(color = "black"), # Keep only x-axis
                axis.text.y = element_text(size = 10, hjust = 0), # Display y-axis text for traits
                axis.ticks.y = element_blank(), # Remove y-axis ticks
                legend.position = "none", # Remove legend
                strip.text = element_text(face = "bold", size = 14, hjust = 0), # Shift titles left
                plot.title = element_text(face = "bold", size = 16, hjust = -0.7)) # Shift plot title left

    # Prepare the table with Odds Ratios [95% CI] and P-values
    data_table <- res %>%
        mutate('Odds Ratio [95% CI]' = sprintf("%.2f [%.2f, %.2f]", or, or_lci95, or_uci95),
                'P-value' = format(pval, digits = 4),
                'SNPs' = paste0(nsnp, " SNPs")) %>%
        select(method, 'Odds Ratio [95% CI]', 'P-value', 'SNPs') %>%
        # Convert to long format for table presentation
        pivot_longer(cols = c('Odds Ratio [95% CI]', 'P-value', 'SNPs'), 
                    names_to = "Metric", values_to = "Value")

    # Generate the table plot
    table_plot <- ggplot(data_table, aes(x = Metric, y = method, label = Value)) +
        geom_text(size = 4, hjust = 0) +
        theme_void() +
        theme(plot.background = element_rect(fill = "white"),
                strip.text = element_blank(),
                panel.spacing.y = unit(1, "lines"),
                plot.margin = margin(0, 0, 0, 0, "cm"),
                legend.position = "none") # Remove legend

    # Add column headers manually by creating a small table
    headers <- data.frame(
    Metric = c("Odds Ratio [95% CI]", "P-value", "SNPs"),
    y = 1,  # Dummy y value
    label = c("OR [95% CI]", "P-value", "SNPs")
    )

    header_plot <- ggplot(headers, aes(x = Metric, y = y, label = label)) +
        geom_text(size = 5, hjust = 0, fontface = "bold") +
        theme_void() +
        theme(plot.background = element_rect(fill = "white"),
                plot.margin = margin(0, 0, 0, 0, "cm"))

    # Combine the headers, table, and plot
    final_plot <- plot_grid(
        plot, 
        plot_grid(header_plot, table_plot, nrow = 2, rel_heights = c(0.1, 1), align = 'v'), 
        nrow = 1,
        rel_widths = c(1.5, 1))

    filename <- paste0(dir, "/mr_results_plot.png")

    # Save or display the plot
    output_width <- 12  # Width in inches
    output_height <- 6  # Height in inches
    ggsave(filename, final_plot, width = output_width, height = output_height, dpi = 300)
}

export_funnel_plot <- function(dat, dir) {
    res_single <- mr_singlesnp(dat)
    p <- mr_funnel_plot(res_single)

    h <- 
    set_panel_size(p = p[[1]], 
                    width=unit(10,"cm"), height=unit(10, "cm"),
                    file = paste0(dir , "/funnel_plot", device = ".png"))
}

export_forest_plot <- function(dat, dir) {
    n <- nrow(dat)
    if (n > 150) {
        warning("Too many SNPs to plot")
    }
    else {
        res_single <- mr_singlesnp(dat)
        p <- mr_forest_plot(res_single)
        
        set_panel_size(p = p[[1]], 
                        width=unit(12,"cm"), height=unit(n / 3, "cm"),
                        file = paste0(dir, "/forest_plot", device = ".png")) %>%
                        suppressWarnings()
    }
}

export_loo_plot <- function(dat, dir) {
    n <- nrow(dat)
    if (n > 150) {
        warning("Too many SNPs to plot")
    }
    else {
        res_loo <- mr_leaveoneout(dat)
        p <- mr_leaveoneout_plot(res_loo)
        
        set_panel_size(p = p[[1]], 
                        width=unit(12,"cm"), height=unit(n / 3, "cm"),
                        file = paste0(dir, "/loo_plot", device = ".png")) %>%
                        suppressWarnings()
    }
}

export_volcano_plot <- function(dat, dir) {
    res_single <- mr_singlesnp(dat)
    p <- ggplot(data=res_single, aes(x=b, y=-log10(p))) + geom_point() + theme_bw()

    set_panel_size(p = p,
                    width=unit(10,"cm"), height=unit(10, "cm"),
                    file = paste0(dir, "/volcano_plot", device = ".png"))
}

export_radial_plot <- function(dir, r_input, egger) {
    ivw <- RadialMR::ivw_radial(r_input = r_input, alpha = 0.05, weights = 3, summary = FALSE)
    p <- RadialMR::plot_radial(r_object = c(ivw, egger), radial_scale = FALSE, 
                        show_outliers = TRUE, scale_match = TRUE)

    set_panel_size(p = p, 
                    width=unit(10,"cm"), height=unit(10, "cm"),
                    file = paste0(dir, "/radial_plot", device = ".png"))
}


### Compute sensitivity analysis ---------------------------------------------------
plotting_and_sensitivity_analysis <- function(res, dat, api_use_ex,
                                              param_mr = TwoSampleMR::default_parameters(),
                                              method_list = full_method_list)
{
    if (api_use_ex) {
            ext <- "_API"
    } else {
            ext <- "_local"
    }

    # Create the output directory if it doesn't exist
    dir <- paste0(output_dir,
                    gsub("[ /]", "-", dat$outcome[1]), 
                    "/", 
                    gsub("[ /]", "-", dat$exposure[1]),
                    ext)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    # Add Heterogeneity analysis
    het <- mr_heterogeneity(dat, parameters = param_mr,
                            method_list = intersect(subset(mr_method_list(), heterogeneity_test)$obj,
                                                           method_list))[, c("method", "Q", "Q_df", "Q_pval")]

    res <- merge(res, het, by = 'method', all.x = TRUE) 

    r_input <- RadialMR::tsmr_to_rmr_format(dat)
    egger_radial <- RadialMR::egger_radial(r_input = r_input, alpha = 0.05, weights = 3, summary = FALSE)
    res_radial <- c(method = "MR Egger Radial", id.exposure = res$id.exposure[1], id.outcome = res$id.outcome[1],
                    outcome = res$outcome[1], exposure = res$exposure[1], nsnp = res$nsnp[1],
                    b = egger_radial$coef[2,1], se = egger_radial$coef[2,2], pval = egger_radial$coef[2,4], 
                    Q = egger_radial$qstatistic, Q_df = egger_radial$df, Q_pval = pchisq(egger_radial$qstatistic, df = egger_radial$df, lower.tail = FALSE))
    res <- rbind(res, as.data.frame(t(res_radial)))
    res[, c("nsnp", "b", "se", "pval", "Q", "Q_df", "Q_pval")] <- lapply(res[, c("nsnp", "b", "se", "pval", "Q", "Q_df", "Q_pval")], as.numeric)
    res <- generate_odds_ratios(res)

    res <- compute_f_stat(dat, res)
    res <- compute_I_squared(dat, res)

    # Plotting
    export_mr_scatter_plots(res, dat, dir)
    export_mr_results_plot(res, dir)
    export_funnel_plot(dat, dir)
    export_forest_plot(dat, dir)
    export_loo_plot(dat, dir)
    export_volcano_plot(dat, dir)
    export_radial_plot(dir, r_input, egger_radial)

    # Export MR results as a CSV file
    export_results(res, dir)

    #res <- rename(res, method = Method)
    return(res)
}

### Multivariable MR plot ---------------------------------------------------------
export_mvmr_plot <- function(res_mv, res_bind, api_use,
                             model_names = c('TRUE' = "Multivariable", 'FALSE' = "Univariable"),
                             colors = c("blue", "red"))
{
    res <- bind_rows(res_mv, res_bind) %>%
            group_by(id.exposure) %>%
            fill(everything(), .direction = "down")

    res$mv <- factor(res$mv, levels = c(TRUE, FALSE))
    res$exposure <- factor(res$exposure, levels = res_mv$exposure)

    # Create the plot
    plot <- ggplot(res, aes(x = or, y = exposure, color = mv)) +
    geom_point(shape = 15, size = 4) +  # Fixed size for points
    geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_manual(values = colors) +
    facet_wrap(~ mv, scales = "free_y", ncol = 1, strip.position = "right", 
                labeller = labeller(mv = model_names)) +
    labs(x = paste("Odds Ratio (with 95% CI) on ", res$outcome[1], " for each trait"),
        y = "",
        title = "Traits") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), # Remove grid
            panel.grid.minor = element_blank(), # Remove minor grid
            panel.border = element_blank(), # Remove border
            axis.line.x = element_line(color = "black"), # Keep only x-axis
            axis.text.y = element_text(size = 10, hjust = 0), # Display y-axis text for traits
            axis.ticks.y = element_blank(), # Remove y-axis ticks
            legend.position = "none", # Remove legend
            strip.text = element_text(face = "bold", size = 14, hjust = 0), # Shift titles left
            plot.title = element_text(face = "bold", size = 16, hjust = -0.3)) # Shift plot title left


    # Prepare the table with Odds Ratios [95% CI] and P-values
    data_table <- res %>%
        mutate('Odds Ratio [95% CI]' = sprintf("%.2f [%.2f, %.2f]", or, or_lci95, or_uci95),
                'P-value' = format(pval, scientific = TRUE, digits = 3),
                'SNPs' = paste0(nsnp, " SNPs")) %>%
        select(exposure, mv, 'Odds Ratio [95% CI]', 'P-value', 'SNPs') %>%
        # Convert to long format for table presentation
        pivot_longer(cols = c('Odds Ratio [95% CI]', 'P-value', 'SNPs'), 
                    names_to = "Metric", values_to = "Value")

    # Generate the table plot
    table_plot <- ggplot(data_table, aes(x = Metric, y = exposure, label = Value, color = mv)) +
    geom_text(size = 4, hjust = 0) +
    facet_wrap(~ mv, scales = "free_y", ncol = 1) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white"),
            strip.text = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            plot.margin = margin(0, 0, 0, 0, "cm"),
            legend.position = "none") # Remove legend

    # Add column headers manually by creating a small table
    headers <- data.frame(
    Metric = c("Odds Ratio [95% CI]", "P-value", "SNPs"),
    y = 1,  # Dummy y value
    label = c("OR [95% CI]", "P-value", "SNPs")
    )

    header_plot <- ggplot(headers, aes(x = Metric, y = y, label = label)) +
    geom_text(size = 5, hjust = 0, fontface = "bold") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white"),
            plot.margin = margin(0, 0, 0, 0, "cm"))

    # Combine the headers, table, and plot
    final_plot <- plot_grid(
    plot, 
    plot_grid(header_plot, table_plot, nrow = 2, rel_heights = c(0.1, 1), align = 'v'), 
    nrow = 1,
    rel_widths = c(1.5, 1))

    # Export each plot as a PNG image
    if (api_use) {
            ext <- "API"
    } else {
            ext <- "local"
    }

    dir <- paste0(output_dir, 
                  gsub("[ /]", "-", res$outcome[1]),
                  "/MVMR_",
                  gsub(" ", "_", sapply(strsplit(na.omit(levels(res$exposure)), "\\|"), 
                       function(x) x[1]) %>% paste(collapse = "")),
                  ext)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    filename <- paste0(dir,"/",model_names[[1]],"_vs_",model_names[[2]],"_plot.png")

    # Save or display the plot
    output_width <- 12  # Width in inches
    output_height <- 6  # Height in inches
    ggsave(filename, final_plot, width = output_width, height = output_height, dpi = 300)

    return(res)
}

plot_paper_comparaison <- function(res, paper_results, api_use) 
{
    paper_df <- fread(paste0(data_dir, "/paper-results/", paper_results)) 
    
    if (!("or" %in% names(paper_df))) {
        paper_df <- generate_odds_ratios(paper_df) 
    }
    res0 <- subset(res, method != "Wald ratio")

    if (api_use) {
            ext <- "_API"
    } else {
            ext <- "_local"
    }

    # Create the output directory if it doesn't exist
    dir <- paste0(output_dir,
                    gsub("[ /]", "-", res$outcome[1]), 
                    "/", 
                    gsub("[ /]", "-", res$exposure[1]),
                    ext)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    res0$mv <- FALSE
    paper_df$mv <- TRUE

    res <- subset(paper_df, method %in% res0$method) %>%
            bind_rows(res0) %>%
            group_by(method) %>%
            fill(everything(), .direction = "up")

    res$mv <- factor(res$mv, levels = c(TRUE, FALSE))
    res$method <- factor(res$method, levels = paper_df$method)

    res <- res[!is.na(res$method), ]

    model_names = c('TRUE' = "Paper results", 'FALSE' = "Reproduction")

    # Create the plot
    plot <- ggplot(res, aes(x = or, y = method, color = mv)) +
            geom_point(shape = 15, size = 4) +  # Fixed size for points
            geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.1) +
            scale_x_log10() +
            geom_vline(xintercept = 1, linetype = "dashed") +
            scale_color_manual(values = c("blue", "red")) +
            facet_wrap(~ mv, scales = "free_y", ncol = 1, strip.position = "right", 
                    labeller = labeller(mv = model_names)) +
            labs(x = paste("Odds Ratio (with 95% CI) on ", res$outcome[1], " for each trait"),
            y = "",
            title = "Methods") +
            theme_minimal() +
            theme(plot.background = element_rect(fill = "white"),
                    panel.grid.major = element_blank(), # Remove grid
                    panel.grid.minor = element_blank(), # Remove minor grid
                    panel.border = element_blank(), # Remove border
                    axis.line.x = element_line(color = "black"), # Keep only x-axis
                    axis.text.y = element_text(size = 10, hjust = 0), # Display y-axis text for traits
                    axis.ticks.y = element_blank(), # Remove y-axis ticks
                    legend.position = "none", # Remove legend
                    strip.text = element_text(face = "bold", size = 14, hjust = 0), # Shift titles left
                    plot.title = element_text(face = "bold", size = 16, hjust = -0.3)) # Shift plot title left


    # Prepare the table with Odds Ratios [95% CI] and P-values
    data_table <- res %>%
        mutate('Odds Ratio [95% CI]' = sprintf("%.2f [%.2f, %.2f]", or, or_lci95, or_uci95),
                'P-value' = format(pval, scientific = TRUE, digits = 3),
                'SNPs' = paste0(nsnp, " SNPs")) %>%
        select(method, mv, 'Odds Ratio [95% CI]', 'P-value', 'SNPs') %>%
        # Convert to long format for table presentation
        pivot_longer(cols = c('Odds Ratio [95% CI]', 'P-value', 'SNPs'), 
                    names_to = "Metric", values_to = "Value")

    # Generate the table plot
    table_plot <- ggplot(data_table, aes(x = Metric, y = method, label = Value, color = mv)) +
    geom_text(size = 4, hjust = 0) +
    facet_wrap(~ mv, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("blue", "red")) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white"),
            strip.text = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            plot.margin = margin(0, 0, 0, 0, "cm"),
            legend.position = "none") # Remove legend

    # Add column headers manually by creating a small table
    headers <- data.frame(
    Metric = c("Odds Ratio [95% CI]", "P-value", "SNPs"),
    y = 1,  # Dummy y value
    label = c("OR [95% CI]", "P-value", "SNPs")
    )

    header_plot <- ggplot(headers, aes(x = Metric, y = y, label = label)) +
    geom_text(size = 5, hjust = 0, fontface = "bold") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white"),
            plot.margin = margin(0, 0, 0, 0, "cm"))

    # Combine the headers, table, and plot
    final_plot <- plot_grid(
    plot, 
    plot_grid(header_plot, table_plot, nrow = 2, rel_heights = c(0.1, 1), align = 'v'), 
    nrow = 1,
    rel_widths = c(1.5, 1))

    filename <- paste0(dir, "/mr_results_plot_comparaison.png")

    # Save or display the plot
    output_width <- 12  # Width in inches
    output_height <- 6  # Height in inches
    ggsave(filename, final_plot, width = output_width, height = output_height, dpi = 300)
}


plot_rucker_model <- function(res, alpha = 0.05, line_length = 50) 
{
    Q_stat <- subset(res, method == "Inverse variance weighted")$Q
    Q_rucker <- subset(res, method == "MR Egger Radial")$Q
    Q_R <- Q_rucker / Q_stat
    df <- subset(res, method == "Inverse variance weighted")$Q_df

    chi_1 <- qchisq(1 - alpha, df = 1)
    chi_L_1 <- qchisq(1 - alpha, df = df)
    chi_L_2 <- qchisq(1 - alpha, df = df - 1)
    
    p <- ggplot() +
        # Line for Q' = Q
        geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1) +
        geom_segment(aes(x = chi_L_1, y = chi_L_1 - chi_1, xend = chi_L_1 + line_length, yend = chi_L_1 - chi_1 + line_length), linetype = "solid", color = "black", size = 1) +
        geom_segment(aes(x = chi_L_1, y = 0, xend = chi_L_1, yend = chi_L_1), linetype = "solid", color = "black", size = 1) +
        geom_segment(aes(x = chi_L_2 + chi_1, y = chi_L_2, xend = chi_L_2 + chi_1 + line_length, yend = chi_L_2), linetype = "solid", color = "black", size = 1) +
        
        geom_point(aes(x = Q_stat, y = Q_rucker), size = 2, color = "red") +
        geom_segment(aes(x = 0, y = 0, xend = Q_stat, yend = Q_rucker), linetype = "dashed", color = "red", size = .5) +

        annotate("text", x = chi_L_1, y = -5, label = expression(chi^2[1-alpha, L-1])) +
        
        # Axis labels
        labs(x = expression(Q), y = expression(Q*"'")) +
        theme_bw()

    # Display the plot
    suppressWarnings(print(p))

    return(data.frame(Q_stat = Q_stat, Q_rucker = Q_rucker, Q_R = Q_R, chi_L_1 = chi_L_1, chi_L_2 = chi_L_2, chi_1 = chi_1))
}