library(ggpubr)
library(readxl)
library(xlsx)
library(tidyverse)
        
setwd("/Users/irene/Desktop/DPANN/Excel_files") 

N2O_read <- read_excel("N2O_electrode_test.xlsx", 
                           sheet="normalized")

N2O_values<-gather(
  N2O_read,
  key = "strain",
  value = "N2O",
  na.rm = TRUE)

N2O_bar<-ggbarplot(N2O_values, x = "strain", y = "N2O", title = "N2O concentrations for each strain",
          ylab= "Normalized N2O concentrations", xlab = "Strain", color = "strain", fill = "strain",
          add = c("mean_se", "jitter"), palette = "jco",col="black",
          position = position_dodge(.8)) +
  scale_color_manual(values = c(V = "black", Z = "black"))+
  geom_hline(yintercept=1) + theme_bw()
N2O_bar

#t.test(N2O_read$`delta nosZ`, N2O_read$putnos2, paired=FALSE)

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="N2O_barplot.svg", N2O_bar)
