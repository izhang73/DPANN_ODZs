library(cowplot)
library(phyloseq)
library(MetBrewer)
library(ggplot2)
library(vegan)
library(dplyr)
library(xlsx)
library(paletteer)
library(Polychrome)

setwd("/Users/irene/Desktop/DPANN/Excel_files/")

#read in all MAG tables from CoverM
MAG_abund_table_raw <- read.xlsx2("CoverM_GTDB_joined.xlsx", 
                                  sheetIndex=1,row.names = 1)

MAG_abund_table_raw <- mutate_all(MAG_abund_table_raw, function(x) as.numeric(as.character(x)))
MAG_matrix <- as.matrix(MAG_abund_table_raw)

tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=2, row.names=1)
tax_mat <- as.matrix(tax_mat)
meta <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=3, row.names=1)

OTU = otu_table(MAG_matrix, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(meta)

#make phyloseq object
MAG_phyloseq = phyloseq(OTU, TAX, samples)

#separate only archaea
Archaea <- subset_taxa(MAG_phyloseq, Domain=="Archaea")

#subset phyloseq object by Location
Arabian <- subset_samples(MAG_phyloseq, Location=="Arabian Sea")
ETNP <- subset_samples(MAG_phyloseq, Location=="Eastern Tropical North Pacific Ocean")
ETSP <- subset_samples(MAG_phyloseq, Location=="Eastern Tropical South Pacific Ocean")

#for datasets with particle fractions, remove all particle and filtered fractions
ETNP_bulk <- subset_samples(ETNP, type=="untreated seawater")
ETSP_bulk <- subset_samples(ETSP, type=="untreated seawater")

#for datasets with particle fractions, subset all particle and filtered fractions
ETNP_particle <- subset_samples(ETNP, type=="particle")
ETSP_particle <- subset_samples(ETSP, type=="particle")

#subset only Archaea from different regions
Archaea_AS <- subset_taxa(Arabian, Domain=="Archaea")
Archaea_ETNP <- subset_taxa(ETNP, Domain=="Archaea")
Archaea_ETNP_bulk <- subset_taxa(ETNP_bulk, Domain=="Archaea")
Archaea_ETSP <- subset_taxa(ETSP, Domain=="Archaea")
Archaea_ETSP_bulk <- subset_taxa(ETSP_bulk, Domain=="Archaea")

#subset only Archaea from particles
Archaea_ETNP_particle <- subset_taxa(ETNP_particle, Domain=="Archaea")
Archaea_ETSP_particle <- subset_taxa(ETSP_particle, Domain=="Archaea")

desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")

desired_ETNP_order <- list("AMAL5","AMAL17","Glass_SRR1509790","Glass_SRR1509797","AMAL18","AMAL1","AMAL13",
                           "Fuchsman_SRR4465037","Tsementzi_SRR3718412","Fuchsman_SRR4465034","Glass_SRR1509798",
                           "Glass_SRR1509792","Fuchsman_SRR4465025","AMAL14","Fuchsman_SRR4465024","Fuchsman_SRR4465031",
                           "Glass_SRR1509793","Fuchsman_SRR4465027","AMAL2","Fuchsman_SRR4465026","Fuchsman_SRR4465030",
                           "Fuchsman_SRR4465032","Tsementzi_SRR3718413","Glass_SRR1509794","Glass_SRR1509799",
                           "Fuchsman_SRR4465029","Fuchsman_SRR4465033","Fuchsman_SRR4465028","Fuchsman_SRR4465036","AMAL7",
                           "AMAL3","AMAL15","AMAL8","AMAL20","Fuchsman_SRR4465035","Glass_SRR1509796","Glass_SRR1509800")

desired_ETNP_bulk_order <- list("AMAL5","AMAL17","Glass_SRR1509790","Glass_SRR1509797","AMAL18","AMAL1","AMAL13",
                                "Fuchsman_SRR4465037","Tsementzi_SRR3718412","Fuchsman_SRR4465034","Glass_SRR1509798",
                                "Glass_SRR1509792","Fuchsman_SRR4465025","AMAL14","Fuchsman_SRR4465024",
                                "Glass_SRR1509793","Fuchsman_SRR4465027","AMAL2","Fuchsman_SRR4465026",
                                "Tsementzi_SRR3718413","Glass_SRR1509794","Glass_SRR1509799",
                                "Fuchsman_SRR4465029","Fuchsman_SRR4465028","Fuchsman_SRR4465036","AMAL7",
                                "AMAL3","AMAL15","AMAL8","AMAL20","Fuchsman_SRR4465035","Glass_SRR1509796","Glass_SRR1509800")

desired_ETSP_order <- list("Stewart_SRR304684","Stewart_SRR304671","Stewart_SRR064444","Stewart_SRR304672","Stewart_SRR304674",
                           "Stewart_SRR304656","Stewart_SRR070081","Ganesh_SRR960580","Ganesh_SRR961671","Stewart_SRR064448",
                           "Stewart_SRR304673","Stewart_SRR304680","Ganesh_SRR961673","Ganesh_SRR961675","Stewart_SRR070084",
                           "Stewart_SRR064450","Stewart_SRR070082","Ganesh_SRR961676","Ganesh_SRR961677","Stewart_SRR304668",
                           "Stewart_SRR304683","Ganesh_SRR961679","Ganesh_SRR961680")

desired_ETSP_bulk_order <- list("Stewart_SRR304684","Stewart_SRR304671","Stewart_SRR064444","Stewart_SRR304672","Stewart_SRR304674",
                                "Stewart_SRR304656","Stewart_SRR070081","Ganesh_SRR960580","Stewart_SRR064448",
                                "Stewart_SRR304673","Stewart_SRR304680","Ganesh_SRR961673","Stewart_SRR070084",
                                "Stewart_SRR064450","Stewart_SRR070082","Ganesh_SRR961676","Stewart_SRR304668",
                                "Stewart_SRR304683","Ganesh_SRR961679")

relabel_ETSP_bulk <- c('Stewart_2008_15m','Stewart_2008_35m','Stewart_2008_50m','Stewart_2008_50m','Stewart_2008_50m',
                       'Stewart_2008_65m','Stewart_2008_70m','Ganesh_2010_70m','Stewart_2008_110m','Stewart_2008_110m','Stewart_2008_110m',
                       'Ganesh_2010_110m','Stewart_2008_150m','Stewart_2008_200m','Stewart_2008_200m','Ganesh_2010_200m','Stewart_2008_500m',
                       'Stewart_2008_800m','Ganesh_2010_1000m','Ganesh_2010_1000m')
relabel_ETNP_bulk <- c('2016_10m','2018_16m','Glass_2013_30m','Glass_2013_30m','2018_45m','2016_53m','2018_60m','Fuchsman_2012_60m','Tsementzi_2013_68m','Fuchsman_2012_70m',
                       'Glass_2013_80m','Glass_2013_85m','Fuchsman_2012_90m','2018_95m','Fuchsman_2012_100m','Glass_2013_100m','Fuchsman_2012_110m','2016_120m','Fuchsman_2012_120m',
                       'Tsementzi_2013_120m','Glass_2013_125m','Glass_2013_125m','Fuchsman_2012_140m','Fuchsman_2012_160m','Fuchsman_2012_180m','2016_185m',
                       '2016_200m','2018_200m','2016_215m','2018_250m','Fuchsman_2012_300m','Glass_2013_300m','Glass_2013_300m')
relabel_AS <- c('2007_130m','2007_150m','2007_200m','2007_400m')

Archaea_AS_plot <- plot_bar(Archaea_AS, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Lakota", 9)) + 
  geom_bar( stat="identity", position="stack") + xlab("Arabian Sea") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_AS)
Archaea_AS_plot$data$Sample <- factor(Archaea_AS_plot$data$Sample, levels = desired_AS_order)
Archaea_AS_plot <- (Archaea_AS_plot) + theme(legend.position = "none")+theme(axis.text.x=element_text(angle=45,hjust=1))
print(Archaea_AS_plot)

Archaea_ETNP_plot <- plot_bar(Archaea_ETNP, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Lakota", 12)) + 
  geom_bar( stat="identity", position="stack")
Archaea_ETNP_plot$data$Sample <- factor(Archaea_ETNP_plot$data$Sample, levels = desired_ETNP_order)
print(Archaea_ETNP_plot)

Archaea_ETNP_bulk_plot <- plot_bar(Archaea_ETNP_bulk, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Lakota", 9)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical North Pacific") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_ETNP_bulk)
Archaea_ETNP_bulk_plot$data$Sample <- factor(Archaea_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
Archaea_ETNP_bulk_plot <- (Archaea_ETNP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1)) +theme(legend.position = "none")
#Archaea_ETNP_bulk_plot <- (Archaea_ETNP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1)) +theme(legend.position = "right")
print(Archaea_ETNP_bulk_plot)

Archaea_ETSP_plot <- plot_bar(Archaea_ETSP, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Lakota", 9)) + 
  geom_bar( stat="identity", position="stack")
Archaea_ETSP_plot$data$Sample <- factor(Archaea_ETSP_plot$data$Sample, levels = desired_ETSP_order)
print(Archaea_ETSP_plot)

Archaea_ETSP_bulk_plot <- plot_bar(Archaea_ETSP_bulk, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Lakota", 9)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical South Pacific") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_ETSP_bulk)
Archaea_ETSP_bulk_plot$data$Sample <- factor(Archaea_ETSP_bulk_plot$data$Sample, levels = desired_ETSP_bulk_order)
Archaea_ETSP_bulk_plot <- (Archaea_ETSP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1))
print(Archaea_ETSP_bulk_plot)

DPANN_plots <- plot_grid(Archaea_AS_plot,Archaea_ETNP_bulk_plot, Archaea_ETSP_bulk_plot, 
          labels = "A",
          align = 'h',
          ncol = 3,
          rel_widths = c(0.5,2,2.75))

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="CoverM_DPANN_ETNP.svg", plot=DPANN_plots)
ggsave(filename="CoverM_DPANN_ETNP.png", plot=DPANN_plots)


#transform into proportion diagrams - basically scale it all up to 100%
scaled_Archaea_ETSPbulk = transform_sample_counts(Archaea_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_Archaea_ETNPbulk = transform_sample_counts(Archaea_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_Archaea_AS = transform_sample_counts(Archaea_AS, function(x) (x / sum(x))*100 )

Archaea_ETSP_bulk_plot_scaled <- plot_bar(scaled_Archaea_ETSPbulk, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Nizami", 5)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical South Pacific") + ylab("Relative Abundance Archaeal Community (%)")+
  scale_x_discrete(labels = relabel_ETSP_bulk)
Archaea_ETSP_bulk_plot_scaled$data$Sample <- factor(Archaea_ETSP_bulk_plot_scaled$data$Sample, levels = desired_ETSP_bulk_order)
Archaea_ETSP_bulk_plot_scaled <- (Archaea_ETSP_bulk_plot_scaled) +theme(axis.text.x=element_text(angle=45,hjust=1))
print(Archaea_ETSP_bulk_plot_scaled)

Archaea_ETNP_bulk_plot_scaled <- plot_bar(scaled_Archaea_ETNPbulk, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Nizami", 5)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical North Pacific") + ylab("Relative Abundance Archaeal Community (%)")+
  scale_x_discrete(labels = relabel_ETNP_bulk)
Archaea_ETNP_bulk_plot_scaled$data$Sample <- factor(Archaea_ETNP_bulk_plot_scaled$data$Sample, levels = desired_ETNP_bulk_order)
#Archaea_ETNP_bulk_plot_scaled <- (Archaea_ETNP_bulk_plot_scaled) +theme(axis.text.x=element_text(angle=45,hjust=1))
Archaea_ETNP_bulk_plot_scaled <- (Archaea_ETNP_bulk_plot_scaled) +theme(axis.text.x=element_text(angle=45,hjust=1))+theme(legend.position = "none")
print(Archaea_ETNP_bulk_plot_scaled)

Archaea_AS_plot_scaled <- plot_bar(scaled_Archaea_AS, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Nizami", 5)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Arabian Sea") + ylab("Relative Abundance Archaeal Community (%)")+
  scale_x_discrete(labels = relabel_AS)
Archaea_AS_plot_scaled$data$Sample <- factor(Archaea_AS_plot_scaled$data$Sample, levels = desired_AS_order)
Archaea_AS_plot_scaled <- (Archaea_AS_plot_scaled) +theme(axis.text.x=element_text(angle=45,hjust=1))+theme(legend.position = "none")
print(Archaea_AS_plot_scaled)

scaled_DPANN_plots <- plot_grid(Archaea_AS_plot_scaled,Archaea_ETNP_bulk_plot_scaled, Archaea_ETSP_bulk_plot_scaled, 
                         labels = "A",
                         align = 'h',
                         ncol = 3,
                         rel_widths = c(0.5,2.1,2.4))

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="Scaled_CoverM_DPANN_ETNP.svg", plot=scaled_DPANN_plots)
ggsave(filename="Scaled_CoverM_DPANN_ETNP.png", plot=scaled_DPANN_plots)

scaled_Archaea_ETSPpart = transform_sample_counts(Archaea_ETSP_particle, function(x) (x / sum(x))*100 )
scaled_Archaea_ETNPpart = transform_sample_counts(Archaea_ETNP_particle, function(x) (x / sum(x))*100 )

Archaea_ETNP_particle_plot <- plot_bar(scaled_Archaea_ETNPpart, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Nizami", 5)) + 
  geom_bar( stat="identity", position="stack")
desired_ETNP_particle_order <- list("Fuchsman_SRR4465031","Fuchsman_SRR4465030","Fuchsman_SRR4465033")
Archaea_ETNP_particle_plot$data$Sample <- factor(Archaea_ETNP_particle_plot$data$Sample, levels = desired_ETNP_particle_order)
Archaea_ETNP_particle_plot <- (Archaea_ETNP_particle_plot) + theme(legend.position = "none")
print(Archaea_ETNP_particle_plot)

Archaea_ETSP_particle_plot <- plot_bar(scaled_Archaea_ETSPpart, fill = "Superphylum") +  scale_fill_manual(values=met.brewer("Nizami", 5)) + 
  geom_bar( stat="identity", position="stack")
desired_ETSP_particle_order <- list("Ganesh_SRR961671","Ganesh_SRR961675","Ganesh_SRR961677","Ganesh_SRR961680")
Archaea_ETSP_particle_plot$data$Sample <- factor(Archaea_ETSP_particle_plot$data$Sample, levels = desired_ETSP_particle_order)
print(Archaea_ETSP_particle_plot)

plot_grid(Archaea_ETNP_particle_plot,Archaea_ETSP_particle_plot, 
          labels = "AUTO",
          align = 'h',
          ncol = 2,
          rel_widths = c(2,3.5))


#subset DPANN out
#subset only Archaea from different regions
DPANN_AS <- subset_taxa(Arabian, Superphylum=="DPANN")
DPANN_ETNP_bulk <- subset_taxa(ETNP_bulk, Superphylum=="DPANN")
DPANN_ETSP_bulk <- subset_taxa(ETSP_bulk, Superphylum=="DPANN")

#subset only Archaea from particles
DPANN_ETNP_particle <- subset_taxa(ETNP_particle, Superphylum=="DPANN")
DPANN_ETSP_particle <- subset_taxa(ETSP_particle, Superphylum=="DPANN")

DPANN_AS_plot <- plot_bar(DPANN_AS, fill = "Phylum") +  scale_fill_manual(values=met.brewer("Juarez", 6)) + 
  geom_bar( stat="identity", position="stack") + xlab("Arabian Sea") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_AS)
DPANN_AS_plot$data$Sample <- factor(DPANN_AS_plot$data$Sample, levels = desired_AS_order)
DPANN_AS_plot <- (DPANN_AS_plot) + theme(legend.position = "none")+theme(axis.text.x=element_text(angle=45,hjust=1))
print(DPANN_AS_plot)

DPANN_ETNP_bulk_plot <- plot_bar(DPANN_ETNP_bulk, fill = "Phylum") +  scale_fill_manual(values=met.brewer("Juarez", 6)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical North Pacific") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_ETNP_bulk)
DPANN_ETNP_bulk_plot$data$Sample <- factor(DPANN_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
DPANN_ETNP_bulk_plot <- (DPANN_ETNP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1)) +theme(legend.position = "none")
#DPANN_ETNP_bulk_plot <- (DPANN_ETNP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1)) +theme(legend.position = "right")
print(DPANN_ETNP_bulk_plot)

DPANN_ETSP_bulk_plot <- plot_bar(DPANN_ETSP_bulk, fill = "Phylum") +  scale_fill_manual(values=met.brewer("Juarez", 6)) + 
  geom_bar( stat="identity", position="stack")+
  xlab("Eastern Tropical South Pacific") + ylab("Relative Abundance (%)")+
  scale_x_discrete(labels = relabel_ETSP_bulk)
DPANN_ETSP_bulk_plot$data$Sample <- factor(DPANN_ETSP_bulk_plot$data$Sample, levels = desired_ETSP_bulk_order)
DPANN_ETSP_bulk_plot <- (DPANN_ETSP_bulk_plot) +theme(axis.text.x=element_text(angle=45,hjust=1))
print(DPANN_ETSP_bulk_plot)

Phylum_DPANN_plots <- plot_grid(DPANN_AS_plot,DPANN_ETNP_bulk_plot, DPANN_ETSP_bulk_plot, 
                         labels = "A",
                         align = 'h',
                         ncol = 3,
                         rel_widths = c(0.5,2.3,2.3))
Phylum_DPANN_plots

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="Phylum_CoverM_DPANN_ETNP.svg", plot=Phylum_DPANN_plots)
#ggsave(filename="Phylum_CoverM_DPANN_ETNP.png", plot=Phylum_DPANN_plots)

#run co-occurrence analysis with cutoff of 0.3
jg = make_network(MAG_phyloseq, "taxa", "jaccard", 0.3)
cooccurence_plot <- plot_network(jg, MAG_phyloseq, "taxa", color = "Phylum",line_weight = 0.5, label = NULL,)
cooccurence_plot <- cooccurence_plot + scale_color_paletteer_d("pals::glasbey")
cooccurence_plot
