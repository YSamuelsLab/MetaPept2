#!/usr/bin/env -S Rscript --vanilla
# -S variable splits the string

############################## Cumulative de novo QC ############################## 

library(argparse)
library(sys)

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(reshape2)
require(gridExtra)
library(ggExtra)
library(extrafont)

# install.packages("extrafont")
# font_import()
# loadfonts(device="all")

### usage:
### ./denovo_QC.R --name MAPP --input ./input --descr ./input/sample_description.csv --max_len 15 --output ./output 

if(!interactive()) pdf(NULL)  # To avoid the creating of Rscript.pdf file
denovo_filename <- 'de novo peptides.csv'

##### Argument Parser
parser <- ArgumentParser(description='QC for de novo peptides')
parser$add_argument('--name', help='the project name')
parser$add_argument('--input', default='./input', help='input directory with PEAKS de novo peptide files')
parser$add_argument('--descr', default='./input/sample_description.csv', help='Sample description file')
parser$add_argument('--min_len', default=8, help='min length of peptides')
parser$add_argument('--max_len', default=0, help='max length of peptides to plot separately')
parser$add_argument('--output', default='./output', help='output directory')
argv <- parser$parse_args()
#####

project <- argv$name
input_dir <- argv$input
descr_file <- argv$descr
min_len <- argv$min_len
max_len <- argv$max_len
output_dir <- argv$output

### DEBUG
#setwd('/data/storwis/dmitryma/Projects2/MAPP')
#project <- 'MAPP'
###

path <- file.path(input_dir, paste0('*/', denovo_filename, collapse =''))
output_dir <- sub('/$', '', output_dir)
denovo_files <- Sys.glob(path)

sample_desc <- read.csv(file = descr_file, header = T, sep = '\t')
alc <- data.frame()
for (f in denovo_files) {
  denovo_df <- read.csv(file = f, header = T, sep = ",", check.names = F)
  alc <- rbind(alc, denovo_df %>% select('Source File', 'ALC (%)', 'length'))
}
data2plot <- alc %>% inner_join(sample_desc, by=c('Source File'='Source_File'))
data2plot$Sample_Name <- as.factor(data2plot$Sample_Name)
data2plot$Name <- paste(data2plot$Sample_Name, data2plot$Sample_Replica, sep="_")
if (max_len == 0){
  max_len <- max(data2plot$length)
}

# unique(alc$`Source File`)
# colnames(denovo_df)

#################### ALC ####################
p <- ggplot(data2plot %>% filter(`ALC (%)` >= 50), aes(x=`Name`, y=`ALC (%)`)) + 
  labs(x="Sample_Replica", y="ALC (%)", title=paste(project), subtitle="ALC distribution") + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5, size=28, face = "bold"),
        plot.subtitle=element_text(size=26, hjust=0.5),
        text = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text=element_text(size=24),
        axis.text.x = element_text(size=20, angle=45, hjust=1))
p

file <- paste(output_dir, '/', project, '.ALC_dist.pdf', sep='')
ggsave(file, width = 15, height = 15, dpi = 600, device = "pdf")

#################### Length ####################
agg_tbl <- data2plot %>% select(Sample_Name, Sample_Replica, length) %>% 
  filter(length >= min_len) %>%
  group_by(Sample_Name, Sample_Replica, length) %>%
  summarise(Count = n(), .groups = 'drop')

agg_tbl_detailed <- agg_tbl %>% filter(length <= max_len)
agg_tbl_condenced <- agg_tbl %>% filter(length > max_len) %>% 
  group_by(Sample_Name, Sample_Replica) %>%
  summarise(Count = sum(length), .groups = 'drop')
agg_tbl_condenced$length = paste(max_len+1, '+', sep='')
data2plot_agg <- rbind(agg_tbl_detailed, agg_tbl_condenced)

data2plot_agg_total <- data2plot_agg %>% 
  group_by(Sample_Name, Sample_Replica) %>%
  summarise(Total = sum(Count), .groups = 'drop')
data2plot_agg_norm <- data2plot_agg %>% inner_join(data2plot_agg_total, by=c('Sample_Name', 'Sample_Replica'))
data2plot_agg_norm$Proportion <- data2plot_agg_norm$Count/data2plot_agg_norm$Total

lengths <- append(as.character(seq(min_len, max_len, by=1)), paste(as.character(max_len+1), '+', sep=''))
data2plot_agg_norm$length <- factor(data2plot_agg_norm$length, levels=lengths)

agg_tbl_sum <- data2plot_agg_norm %>% group_by(Sample_Name, length) %>%
  summarize(Mean = mean(Proportion), Sem = sd(Proportion) / sqrt(n()), .groups = 'drop')

# saving the plot
p <- ggplot(agg_tbl_sum, aes(x=as.factor(length), y=Mean, fill = Sample_Name)) + 
  labs(x="Peptide length", y="", fill='Sample', title=paste(project), subtitle="peptide length distribution") + 
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_errorbar(
    aes(ymin = Mean - Sem, ymax = Mean + Sem),
    width = 0.4, size=0.5,
    position = position_dodge(width = 0.9)
  ) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"),
        plot.subtitle=element_text(size=18, hjust=0.5),
        text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=18),
        axis.text.x = element_text(size=18, angle=0, hjust=1))
p

file <- paste(output_dir, '/', project, '.pept_length_dist.pdf', sep='')
ggsave(file, width = round(15 * (max_len+1)/16, digits = 0), height = 15, units = 'cm', dpi = 600, device = "pdf")

# saving the plot
p <- ggplot(data2plot_agg_norm, aes(x=as.factor(length), y=Proportion, fill = Sample_Name)) + 
  labs(x="Peptide length", y="Proportion", fill='Sample', title=paste(project), subtitle="peptide length distribution") + 
  geom_boxplot(lwd=0.3) +
  theme(plot.title = element_text(hjust = 0.5, size=22, face = "bold"),
        plot.subtitle=element_text(size=24, hjust=0.5),
        text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=18),
        axis.text.x = element_text(size=18, angle=0, hjust=1))
p

file <- paste(output_dir, '/', project, '.pept_length_dist.box.pdf', sep='')
ggsave(file, width = round(15 * (max_len+1)/16, digits = 0), height = 15, units = 'cm', dpi = 600, device = "pdf")
#unlink(file)

