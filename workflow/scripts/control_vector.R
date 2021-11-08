require(tidyverse)
require(ggplot2)
require(patchwork)
require(viridisLite)

import_data <- function(files) {
    df <- NULL

    for (f in files) {
        if (file.info(f)$size > 0) {
            name <- gsub("\\.bed", "", basename(f))
            print(cat("Loading", f, name, "\n"))

            my_bed <- read.delim(f, sep="\t", header=FALSE)
            colnames(my_bed) <- c("chr", "start", "end", "beta", "depth", "context")
            my_bed$sample <- rep(name, nrow(my_bed))

            df <- rbind(df, my_bed)
        }
    }

    return(df)
}

create_plot <- function(lam_files, puc_files, out_files) {
    lam <- import_data(lam_files)
    puc <- import_data(puc_files)
    print(lam)
    print(puc)

    n_samples_l <- length(unique(lam$sample))
    n_samples_p <- length(unique(puc$sample))
    print(n_samples_l)
    print(n_samples_p)

    # Unmethylated control
    if (!is.null(lam)) {
        cat("making lambda plots\n")
        topleft <- ggplot(lam, aes(x=sample, y=depth)) +
            geom_boxplot(color='#357BA2FF') +
            theme_bw() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=25),
                plot.title = element_text(size=25, hjust=0.5),
                plot.subtitle = element_text(size=15, hjust=0.5),
            ) +
            scale_color_manual() +
            ggtitle("Unmethylated", subtitle=paste("N =", n_samples_l, "Samples")) +
            ylab("Coverage") +
            xlab("")

        bottomleft <- ggplot(lam, aes(x=sample, y=beta)) +
            geom_boxplot(color='#357BA2FF') +
            theme_bw() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=25),
                plot.title = element_text(size=25, hjust=0.5),
                plot.subtitle = element_text(size=15, hjust=0.5),
            ) +
            scale_color_manual() +
            ylim(c(0, 1)) +
            ylab("Beta") +
            xlab("")
    }

    # Methylated control
    if (!is.null(puc)) {
        cat("making puck plots\n")
        topright <- ggplot(puc, aes(x=sample, y=depth)) +
            geom_boxplot(color='#357BA2FF') +
            theme_bw() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=25),
                plot.title = element_text(size=25, hjust=0.5),
                plot.subtitle = element_text(size=15, hjust=0.5),
            ) +
            scale_color_manual() +
            ggtitle("Methylated", subtitle=paste("N =", n_samples_p, "Samples")) +
            ylab("") +
            xlab("")

        bottomright <- ggplot(puc, aes(x=sample, y=beta)) +
            geom_boxplot(color='#357BA2FF') +
            theme_bw() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=25),
                plot.title = element_text(size=25, hjust=0.5),
                plot.subtitle = element_text(size=15, hjust=0.5),
            ) +
            scale_color_manual() +
            ylim(c(0, 1)) +
            ylab("") +
            xlab("")
    }

    # Create plot
    if (exists("topleft") & exists("topright")) {
        cat("trying to make patchwork plot\n")
        layout <- '
        AB
        CD
        '

        pdf(file=out_files, width=7, height=10)

        topleft + topright + bottomleft + bottomright +
            patchwork::plot_annotation(tag_levels = 'A', title = 'Control Vectors') +
            patchwork::plot_layout(design = layout) & theme(plot.tag = element_text(face = 'bold')) & theme(plot.title = element_text(size=25, hjust=0.5))

        dev.off()
    } else {
        cat("trying to make unfilled plot\n")
        pdf(file=out_files, width=7, height=10)

        plot.new()
        text(x=0.5, y=0.5, "NO PLOT CREATED.\nLIKELY REASON: NO DEPTH FOR CONTROL VECTOR(S)\nIN SUPPLED BED FILES")

        dev.off()
    }
}

create_plot(snakemake@input[["lambda_files"]], snakemake@input[["puc19_files"]], snakemake@output[["pdf"]])
