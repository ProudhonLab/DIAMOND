#!/usr/local/bin/Rscript
# ****************************************************
# Title      : plot_read_or_uid_length.R
# Author     : Marc Michel
# Project    : https://github.com/michel-m
# Date       : November 2018
# Description: n/a
# ****************************************************
if (!requireNamespace("tools", quietly = TRUE))
  install.packages("tools")
library(tools)
library(optparse)
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(plotly))

option_list = list(
  make_option(c("-t", "--title"), type = "character", default = "Untitled",
              help = "Plot title; defaults to '%default'", metavar = "'a title'"),
  make_option(c("-o", "--out"), type = "character", default = "~/read_or_uid_length_plot.png",
              help = "Output file path; suffix will decide file type (html, png, jpg...); defaults to %default",
              metavar = "path"),
  make_option(c("-c", "--control"), type = "character", default = NULL, help = "control sample file path")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, positional_arguments = 1)

set.seed(1)
pdf(NULL)

input_files <- opt$args
plot_title <- opt$options$title
output_file_path <- opt$options$out

df <- read.csv2(input_files[1], header = F, sep = ',', stringsAsFactors = F)
tbl <- table(df)
df_final <- as.data.frame(tbl)
names(df_final) <- c("Position", "Frequency")
df_final$Type <- "tumor"
df_final$Position <- levels(droplevels(df_final$Position))
df_final$Position <- as.numeric(df_final$Position)
legend_position <- "none"

if (!is.null(opt$options$control)) {
  df2 <- read.csv2(opt$options$control, header = F, sep = ',', stringsAsFactors = F)
  tbl2 <- table(df2)
  df_final2 <- as.data.frame(tbl2)
  names(df_final2) <- c("Position", "Frequency")
  df_final2$Type <- "control"
  df_final2$Position <- levels(droplevels(df_final2$Position))
  df_final2$Position <- as.numeric(df_final2$Position)
  legend_position <- "right"

  df_final <- bind_rows(df_final, df_final2)
}

p <- ggplot(df_final, aes(x = Position, y = Frequency, fill = Type,
                          label = ifelse(Frequency > (max(df_final$Frequency) * .33),
                                         prettyNum(Frequency, big.mark = ',', scientific = F), ""),
                          text = sprintf("Frequency: %s<br>Position in read: %s",
                                         prettyNum(Frequency, big.mark = ',', scientific = F), Position))) +
  geom_bar(stat = "identity", position = "dodge", color = "black", cex = .05) +
  scale_x_continuous(name = "Read length (n)",
                     breaks = c(df_final[df_final$Frequency == max(df_final$Frequency), ]$Position,
                                seq(0, max(df_final$Position), by = 10))) +
  scale_y_continuous(name = "Frequency", labels = scales::comma) +
  geom_text(nudge_y = max(df_final$Frequency) / 40) +
  theme_classic() +
  theme(legend.position = legend_position) +
  ggtitle(plot_title)

if (!is.null(opt$options$control)) {
  p <- p + scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                             name = "Sample type",
                             breaks = c("control", "tumor"),
                             labels = c("Control", "Tumor"))
}

output_file_ext <- file_ext(output_file_path)
if (output_file_ext == "html") {
  htmlwidgets::saveWidget(ggplotly(p, tooltip = c('text')), file = output_file_path)
} else {
  ggsave(filename = output_file_path, plot = p, scale = 1.1, height = 5.73,
         width = 9.11, units = "in", dpi = 320, device = output_file_ext)
}
