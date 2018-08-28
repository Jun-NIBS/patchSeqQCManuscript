library(homologene)
library(biomaRt)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(cowplot)

# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# human_annot<-getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mart) %>% 
#   mutate(gene= hgnc_symbol)
# 
# chen_expr = read.csv('data-raw/chen/tpmMatrix.genes', row.names = 1, sep = '\t')
# 
# 
# chen_expr = chen_expr + 1

chen_expr = read.csv('data-raw/chen/GSE77564_single_neuron_exp.txt', row.names = 1)
chen_expr = chen_expr + 1
rownames(chen_expr) %<>% make.names(., unique = T)

chen_gene_names = rownames(chen_expr)

chen_ephys = read_excel('data-raw/chen/13238_2016_247_MOESM1_ESM.xlsx', na = 'ND', col_names = T, col_types = 'numeric', )[1:20, -1]
chen_ephys_names = read_excel('data-raw/chen/13238_2016_247_MOESM1_ESM.xlsx', na = 'ND', col_names = T, col_types = 'text')[1:20, 1]
chen_ephys = cbind(chen_ephys_names, chen_ephys) %>% as.data.frame

chen_ephys %<>% rename(sample_id = Sample_ID)

chen_df = cbind(chen_ephys, chen_expr %>% t()) %>% as.data.frame()
chen_df$cell_id = chen_df$sample_id


chen_meta_stats_df = calculateDatasetQCMeasures(chen_df, chen_gene_names, 
                           protein_genes = human_annot %>% 
                             filter(gene_biotype == 'protein_coding') %>% pull(gene) %>% make.names,
                           mito_genes = human_annot %>% filter(gene_biotype == 'Mt_tRNA') %>% pull(gene) %>% make.names() %>%
                             intersect(., chen_gene_names), min_expr_val = 1) 

chen_meta_stats_df %<>% rename(sample_id = cell_id)


chen_plot_df = chen_df[chen_df$AP == 1, ]

chen_plot_df[chen_plot_df$AP == 1, 'major_type'] = 'MatureNeuron'
chen_plot_df[chen_plot_df$AP == 1, 'contam_type'] = 'MatureNeuron'


source('~/patchSeqQC/R/plotMarkerHeatmaps.R')

ann_colors = list(Markers = c(Ndnf_on = "indianred1", Ndnf = "indianred1", Pyramidal = "turquoise", Pyramidal_on = "turquoise",
                              Microglia = "gray50", Sncg_on = "purple", Astrocyte = "green4", Pvalb = "red",
                              Inhibitory = "red",
                              Endothelial = "brown",
                              Oligodendrocyte = "sandybrown",MatureNeuron = "red",
                              OPC = 'darkorange4'),
                  CellTypes = c(Ndnf = "indianred1", Pyramidal = "turquoise", Pvalb = "red", Sncg = "purple", Astrocyte = "green4", 
                                MatureNeuron = "red"))


chen_heatmap = plotMarkerHeatmap(markerlist = plot_marker_list, # named list of lists defining which markers to show
                                  expr_matrix = chen_plot_df, # data frame that combines gene expr with metadata, each row is a single-cell
                                  show_legend = T, # show color bar ?
                                  show_cell_labels = T, # show sample names in heatmap (defined by rownames of expr_matrix)
                                 max_log_expr = 8
                                 
)

# chen_heatmap = plot_grid(chen_heatmap$gtable, bardy_astro_expr_plot, nrow = 1, rel_widths = c(4, 1.1))

ggsave(filename = 'plots/chen_heatmap.pdf', plot = chen_heatmap$gtable, units = "in", width = 8, height = 5)



