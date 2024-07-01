#!/usr/local/bin/Rscript
# ****************************************************
# Title      : adapter_stats.R
# Author     : Marc Michel
# Project    : https://github.com/michel-m
# Date       : May 2018
# Description: n/a
# ****************************************************
if (!requireNamespace("tools", quietly = TRUE))
  install.packages("tools")
library(tools)
library(optparse)
library(ggplot2)
library(readr)
suppressMessages(library(dplyr))
suppressMessages(library(plotly))
suppressMessages(library(plyr))

option_list = list(
  make_option(c("-t", "--title"), type = "character", default = "Untitled",
              help = "Plot title; defaults to '%default'", metavar = "'a title'"),
  make_option(c("-o", "--out"), type = "character", default = "~/read_or_uid_length_plot.png",
              help = "Output file path; suffix will decide file type (html, png, jpg...); defaults to %default",
              metavar = "path")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, positional_arguments = 1)

set.seed(1)
pdf(NULL)

input_files <- opt$args
plot_title <- opt$options$title
output_image_path <- opt$options$out

if (length(input_files) == 0) {
  stop("The path of a file containing adapter sequences must be provided.", call. = F)
}

input_file_path <- input_files[1]
output_file_path <- sub('\\.UMIs.length_filtered.csv$', '', input_file_path)

df <- suppressMessages(read_csv2(input_file_path, col_names = F))
df <- count(df, 'X1')
df <- arrange(df, -freq)
names(df) <- c('umi', 'frequency')
write.csv(df, file = paste0(output_file_path, '.UMIs.length_filtered.frequency.csv'), row.names = F, quote = F)

tbl <- table(df$frequency)
df_final <- as.data.frame(tbl)
names(df_final) <- c("Frequency", "Count")
df_final$Frequency <- levels(droplevels(df_final$Frequency))
df_final$Frequency <- as.numeric(df_final$Frequency)

p <- ggplot(df_final, aes(x = Frequency, y = Count)) +
  geom_bar(stat = "identity", cex = .05, fill = "#F8766D") +
  scale_x_continuous(name = "UMIs frequency") +
  scale_y_continuous(name = "Count", labels = scales::comma) +
  theme_classic() +
  ggtitle(plot_title)

output_image_ext <- file_ext(output_image_path)
if (output_image_ext == "html") {
  htmlwidgets::saveWidget(ggplotly(p, tooltip = c('text')), file = output_image_path)
} else {
  ggsave(filename = output_image_path, plot = p, scale = 1.1, height = 5.73,
         width = 9.11, units = "in", dpi = 320, device = output_image_ext)
}
