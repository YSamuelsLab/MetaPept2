#!/usr/bin/env -S Rscript --vanilla
# -S variable splits the string

############################## Cumulative de novo seq ALC distribution ############################## 

library(argparse)

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
### ./Q_vs_ALC.WdeplMice.R --name Microbiom --proxy ProxyTrp --comb combined_scan_integration.csv --prism prism_unique_scan_integration.csv --output output 

if(!interactive()) pdf(NULL)  # To avoid the creating of Rscript.pdf file

##### Argument Parser
parser <- ArgumentParser(description='QC for Metapept')
parser$add_argument('--name', help='the name of the project')
parser$add_argument('--comb', help='CSV file with PRISM/MQ combined hits')
parser$add_argument('--prism', help='CSV file with PRISM unique hits')
parser$add_argument('--proxy', help='proxy amino acid (for example: ProxyPhe, ProxyTrp et.)')
parser$add_argument('--output', help='output directory')
argv <- parser$parse_args()
#####

sample <- argv$name
proxy <- argv$proxy
comb_file <- argv$comb
prism_file <- argv$prism
output_dir <- argv$output

### DEBUG
#sample <- 'microbiom'
#proxy <- 'ProxyTrp'
#comb_file <- '/data/storwis/dmitryma/Projects2/ProxyW_mice/MetaPept/output/combined_scan_integration.csv'
#prism_file <- '/data/storwis/dmitryma/Projects2/ProxyW_mice/MetaPept/output/prism_unique_scan_integration.csv'
#output_dir <- '/data/storwis/dmitryma/Projects2/ProxyW_mice/Rational_filtering/'
###

output_dir <- sub('/$', '', output_dir)

Combined_df <- read.csv(file = comb_file, header = T, sep = "\t")
PRISM_df <- read.csv(file = prism_file, header = T, sep = "\t")

Combined_df <- Combined_df %>% select(Sequence, Best_Q, Best_ALC, Categories, Length)
PRISM_df <- PRISM_df %>% select(Sequence, Best_Q, Best_ALC, Categories, Length)

# Create a copy of the filtered dataframe
Combined_df_dup2 <- Combined_df
PRISM_df_dup2 <- PRISM_df

# Initialize the "Category" column with empty strings
PRISM_df$Category <- ""
PRISM_df_dup2$Category <- ""

Combined_df$Category <- ""
Combined_df_dup2$Category <- ""

# Step 1: Set "Canonical" for rows with "CDS" in Categories column
PRISM_df$Category[PRISM_df$Categories == "CDS"] <- "Canonical"
PRISM_df_dup2$Category[PRISM_df_dup2$Categories == "CDS"] <- "Canonical"

Combined_df$Category[Combined_df$Categories == "CDS"] <- "Canonical"
Combined_df_dup2$Category[Combined_df_dup2$Categories == "CDS"] <- "Canonical"

# Step 2: Add rows for condition 2 (Cryptic)
condition2 <- c(
  "OffFrame", "UTR5", "Frameshift", "OffFrame,Extra",
  "Frameshift,Extra", "UTR5,Extra", "Extra,Frameshift",
  "Extra,UTR5", "Extra,OffFrame"
)
PRISM_df$Category[PRISM_df$Categories %in% condition2] <- "Cryptic"
Combined_df$Category[Combined_df$Categories %in% condition2] <- "Cryptic"

# Step 3: Add rows for condition 3 (Proxy-associated)
condition3 <- c(
  "OffFrame,Extra", "Frameshift,Extra", "UTR5,Extra",
  "Extra,Frameshift", "Extra", "Extra,UTR5", "Extra,OffFrame"
)
PRISM_df_dup2$Category[PRISM_df_dup2$Categories %in% condition3] <- proxy
Combined_df_dup2$Category[Combined_df_dup2$Categories %in% condition3] <- proxy

# Combine the original and duplicated dataframes
results_table_PRISM <- bind_rows(PRISM_df, PRISM_df_dup2)
results_table_Combined <- bind_rows(Combined_df, Combined_df_dup2)

# Filter out empty values in the "Category" column
results_table_PRISM <- results_table_PRISM %>% filter(Category != "")
results_table_Combined <- results_table_Combined %>% filter(Category != "")

# Print the unique values in the "Category" column
paste(unique(results_table_PRISM$Category))

linecolors <- c("#714C02", "#01587A", "#a60c0c")
fillcolors <- c("#9D6C06", "#077DAA", "#d60d0d")

plot_data1 <- results_table_Combined %>% 
  select(Sequence, Best_Q, Best_ALC, Category, Length) %>% 
  filter(Length > 7) %>% distinct() 
plot_data2 <- results_table_PRISM %>% 
  select(Sequence, Best_Q, Best_ALC, Category, Length) %>% 
  filter(Length > 7) %>% distinct()


plot1 <- ggplot(plot_data1, aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  # the position_jitter is removed due to Best_Q=0 values are misrepresented!
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors) +
  scale_fill_manual(values=fillcolors) +
  ggtitle("\nCombined data") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.title=element_blank(), legend.position = c(1,1), legend.justification = c(1, 1),  legend.box.margin = margin(2, 2, 2, 2), legend.spacing.y = unit(-1, "mm"), legend.spacing.x = unit(0, "mm"), legend.key.height = unit(3, 'mm'), legend.key.width = unit(3, 'mm')) + 
  xlim(0, 0.2) + ylim(50, 100)

plot1 <- ggMarginal(plot1, margins = 'y', type="density", fill = "gray")

plot2 <- ggplot(plot_data2, aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors) +
  scale_fill_manual(values=fillcolors) +
  ggtitle("\nPRISM unique data") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.title=element_blank(), legend.position = c(1,1), legend.justification = c(1, 1), legend.box.margin = margin(2, 2, 2, 2), legend.spacing.y = unit(-1, "mm"), legend.spacing.x = unit(0, "mm"), legend.key.height= unit(3, 'mm'), legend.key.width = unit(3, 'mm')) + 
  xlim(0, 0.2) + ylim(50, 100)

plot2 <- ggMarginal(plot2, margins = 'y', type="density", fill = "gray")

### Proxy-associated
linecolors_phe <- c("#a60c0c")
fillcolors_phe <- c("#d60d0d")

plot3 <- ggplot(plot_data1 %>% filter(Category %in% c(proxy)), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +  # !!! JITTER KILLS VALUES CLOSE TO 0 !!!
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors_phe) +
  scale_fill_manual(values=fillcolors_phe) +
  ggtitle(paste("\nCombined data\n(", proxy, ")")) +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot3 <- ggMarginal(plot3, margins = 'y', type="density", fill = "#f56262")
 
plot4 <- ggplot(plot_data2 %>% filter(Category %in% c(proxy)), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors_phe) +
  scale_fill_manual(values=fillcolors_phe) +
  ggtitle(paste("\nPRISM unique data\n(", proxy, ")")) +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot4 <- ggMarginal(plot4, margins = 'y', type="density", fill = "#f56262")

### Cryptic
linecolors_phe <- c("#01587A")
fillcolors_phe <- c("#077DAA")

plot5 <- ggplot(plot_data1 %>% filter(Category %in% c("Cryptic")), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors_phe) +
  scale_fill_manual(values=fillcolors_phe) +
  ggtitle("\nCombined data\n(Cryptic)") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot5 <- ggMarginal(plot5, margins = 'y', type="density", fill = "#6cb5f5")

plot6 <- ggplot(plot_data2 %>% filter(Category %in% c("Cryptic")), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors_phe) +
  scale_fill_manual(values=fillcolors_phe) +
  ggtitle("\nPRISM unique data\n(Cryptic)") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot6 <- ggMarginal(plot6, margins = 'y', type="density", fill = "#6cb5f5")

### Canonical
linecolors_can <- c("#714C02")
fillcolors_can <- c("#9D6C06")

plot7 <- ggplot(plot_data1 %>% filter(Category %in% c("Canonical")), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size =1) +
  geom_point(shape = 21, alpha = 0.5, size =1) +
  scale_color_manual(values=linecolors_can) +
  scale_fill_manual(values=fillcolors_can) +
  ggtitle("\nCombined data\n(Canonical)") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot7 <- ggMarginal(plot7, margins = 'y', type="density", fill = "#f5b942")

plot8 <- ggplot(plot_data2 %>% filter(Category %in% c("Canonical")), aes(Best_Q, Best_ALC, colour = Category, fill = Category)) +
  #geom_point(position=position_jitter(h=0.0001, w=0.0001), shape = 21, alpha = 0.5, size = 1) +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  scale_color_manual(values=linecolors_can) +
  scale_fill_manual(values=fillcolors_can) +
  ggtitle("\nPRISM unique data\n(Canonical)") +
  theme_bw() + theme(plot.title = element_text(family="Arial", size=9, margin=margin(b=0, unit="pt")), text=element_text(family="Arial", size=9), legend.position = "none") +
  xlim(0, 0.2) + ylim(50, 100)

plot8 <- ggMarginal(plot8, margins = 'y', type="density", fill = "#f5b942")

margin = theme(plot.margin = unit(c(0,0,0,2), "mm"))
g <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2, nrow=4)
ggsave(paste(output_dir, '/',  sample, '.combined_VS_PRISM_unique.png', sep=''), width = 16, height = 32, dpi = 600, units = "cm", plot=g)
ggsave(paste(output_dir, '/',  sample, '.combined_VS_PRISM_unique.pdf', sep=''), width = 16, height = 32, dpi = 600, units = "cm", plot=g)

#####################################################################################################

