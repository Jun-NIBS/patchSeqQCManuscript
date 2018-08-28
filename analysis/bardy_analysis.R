library(homologene)
library(biomaRt)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(cowplot)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
human_annot<-getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mart) %>% 
  mutate(gene= hgnc_symbol)

bardy_expr = read.csv('data-raw/bardy/tpmMatrix.genes', row.names = 1)
bardy_expr = bardy_expr + 1

bardy_gene_names = bardy_expr[, 1:2] %>% as.data.frame %>% tibble::rownames_to_column(var = 'ensembl_gene_id')

bardy_gene_names = left_join(bardy_gene_names, human_annot %>% distinct(ensembl_gene_id, .keep_all = T), by = 'ensembl_gene_id')
bardy_gene_names$gene_name_norm = apply(bardy_gene_names, 1, function(r){
  if (!is.na(r[4]) && (r[4] != '')){
    return(r[4])
  } else{
    return(r[1])
  }
})
rownames(bardy_expr) = bardy_gene_names$gene_name_norm %>% make.names(., unique = T)

bardy_meta = read.csv('data-raw/bardy/metadata_56cells.csv')

bardy_df = cbind(bardy_meta, 
                 bardy_expr  %>% t()) #
rownames(bardy_df) = bardy_df$cell.number
colnames(bardy_df) = make.names(colnames(bardy_df), unique = T)

bardy_star_qc_df = read.csv('data-raw/bardy/star_qc_data_frame.csv')

bardy_star_qc_df = bardy_star_qc_df
bardy_star_qc_df$sample_id = bardy_star_qc_df$Sample
bardy_star_qc_df$major_type = 'Pyramidal'
bardy_star_qc_df$contam_type = 'Pyramidal'

bardy_star_qc_df %<>% rename(read_count = total_reads, ercc_pct = spike_in_pct)

bardy_ercc_count_mat = read.csv('data-raw/bardy/ercc_count_mat.genes', row.names = 1)


bardy_df = bardy_expr %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
bardy_df$sample_id = bardy_star_qc_df$Sample

bardy_df = merge(bardy_df, bardy_star_qc_df , by = 'sample_id')
bardy_df$cell_id = bardy_df$sample_id

bardy_meta_stats_df = calculateDatasetQCMeasures(bardy_df, 
                                                 dataset_genes = bardy_gene_names$gene %>% make.names, 
                                                                 protein_genes = human_annot %>% 
                                                                   filter(gene_biotype == 'protein_coding') %>% pull(gene) %>% make.names,
                                                 mito_genes = human_annot %>% filter(gene_biotype == 'Mt_tRNA') %>% pull(gene) %>%
                                                                intersect(., bardy_gene_names$gene), min_expr_val = 1) 

bardy_meta_stats_df %<>% rename(sample_id = cell_id)

bardy_ercc_counts = bardy_ercc_count_mat %>% colSums %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
colnames(bardy_ercc_counts)[2] = 'ercc_sum'

bardy_star_qc_df = merge(bardy_star_qc_df, bardy_ercc_counts, by = 'sample_id') %>% mutate(has_ercc = ercc_sum > 10000, ercc_pct = 100 * (ercc_sum / read_count),
                                                                                           exon_pct = 100 * (exon_sums / (read_count - spike_in_sums)))

bardy_star_qc_df = merge(bardy_star_qc_df, bardy_meta_stats_df)

intersecting_cols = intersect(colnames(bardy_star_qc_df),colnames(bardy_df))

use_cols = setdiff(colnames(bardy_df), intersecting_cols)

bardy_df = merge(cbind('sample_id' = bardy_df$sample_id, bardy_df[, use_cols]), bardy_star_qc_df , by = 'sample_id')

humanMarkers = fullMarkerList

humanMarkers$Astrocyte = c('GJA1', 'SOX9', 'EDNRB', 'MLC1', 
                           'DIO2', 'FGFR3', 'MLC1',
                           'SLC4A4', 'GFAP', 'AQP4', 'ALDH1L1'
                           # 'SLC1A3'
                           )

astro_sorted_genes = bardy_df[bardy_df$FTC_simple == 0, humanMarkers$Astrocyte] %>% colMeans %>% sort(decreasing = T) %>% names %>% unlist
humanMarkers$Astrocyte = astro_sorted_genes

humanMarkers$Oligodendrocyte = c('GPR37', 'SOX10', 'OPALIN', 'MOG', 'MOBP', 'MBP', 'PLP1')
humanMarkers$Microglia = c('RUNX1', 'CD83', 'TMEM119', 'CX3CR1', 'TNF', 'CCL3', 'CCL2', 'IL1A', 'TLR2')

bardy_df %<>% arrange(FTC_simple)

bardy_df$major_type = 'Pyramidal'
bardy_df$contam_type = 'Pyramidal'
rownames(bardy_df) = bardy_df$sample_id

bardy_df$contam_sum = bardy_df[intersect(humanMarkers$Astrocyte, bardy_gene_names %>% pull(gene))] %>% log2() %>% rowSums()

plot_cell_types_bardy = c('Astrocyte', 'Microglia')
plot_marker_list = c(humanMarkers[plot_cell_types_bardy])

bardy_plot_df = bardy_df[bardy_df$FTC_simple %in% c(0, 5) & !bardy_df$sample_id %in% c('C130', 'C132'), ]

bardy_plot_df[bardy_plot_df$FTC_simple == 0, 'major_type'] = 'Astrocyte'
bardy_plot_df[bardy_plot_df$FTC_simple == 0, 'contam_type'] = 'Astrocyte'

bardy_plot_df[bardy_plot_df$FTC_simple %in% c( 4, 5), 'major_type'] = 'MatureNeuron'
bardy_plot_df[bardy_plot_df$FTC_simple %in% c(4, 5), 'contam_type'] = 'MatureNeuron'

# bardy_plot_df = bardy_plot_df[bardy_plot_df$num_genes_exon_tpm_5 > 2000, ]

source('~/patchSeqQC/R/plotMarkerHeatmaps.R')
bardy_heatmap = plotMarkerHeatmap(markerlist = plot_marker_list, # named list of lists defining which markers to show
                  expr_matrix = bardy_plot_df, # data frame that combines gene expr with metadata, each row is a single-cell
                  show_legend = T, # show color bar ?
                  show_cell_labels = T # show sample names in heatmap (defined by rownames of expr_matrix)
)

ann_colors = list(Markers = c(Ndnf_on = "indianred1", Ndnf = "indianred1", Pyramidal = "turquoise", Pyramidal_on = "turquoise",
                              Microglia = "gray50", Sncg_on = "purple", Astrocyte = "green4", Pvalb = "red",
                              Inhibitory = "red",
                              Endothelial = "brown",
                              Oligodendrocyte = "sandybrown",MatureNeuron = "red",
                              OPC = 'darkorange4'),
                  CellTypes = c(Ndnf = "indianred1", Pyramidal = "turquoise", Pvalb = "red", Sncg = "purple", Astrocyte = "green4", 
                                MatureNeuron = "red"))



bardy_spikes = read_xlsx('~/ephys_analysis/data/bardy_patchseq/Bardy_2016_MolPsychiatry_spikeins.xlsx', skip = 1)

bardy_plot_df = merge(bardy_plot_df, bardy_spikes %>% dplyr::rename(sample_id = 'X__1', has_ercc = "Stock Dilution\r\nof ERCC RNA Spike-In Mix Used"    ) )



# df = read_excel('~/Genes vs cell size.xlsx')
# cap_d = df[df$X__1 == 'Capacitance_Meanfit[pF]', ]
# cap_d = cap_d[, -1] %>% as.data.frame()
# rownames(cap_d) = 'cap'
# cap_d = cap_d %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
# 
# bardy_plot_df = left_join(bardy_plot_df, cap_d)
bardy_plot_df$Astro_sum  = bardy_plot_df[intersect(humanMarkers$Astrocyte, bardy_gene_names %>% pull(gene))] %>% log2() %>% rowSums()

bardy_astro_expr_plot = ggplot(bardy_plot_df , 
       aes(x = contam_type, y = Astro_sum,  group = contam_type)) + geom_violin() + geom_jitter(width = .25) + 
  ylab('Astrocyte marker expr \n(log2 TPM+1)') + xlab('Firing phenotype class') + xlab('') + 
  theme(axis.text.x=element_text(angle=30, hjust=1)) 
  
bardy_heatmap = plot_grid(bardy_heatmap$gtable, bardy_astro_expr_plot, nrow = 1, rel_widths = c(4, 1.1))

ggsave(filename = 'plots/bardy_heatmap.pdf', plot = bardy_heatmap, units = "in", width = 10, height = 5)
# 
# 
# bardy_plot_df$Astro_sum  = bardy_plot_df[intersect(humanMarkers$Astrocyte, bardy_gene_names %>% pull(gene))] %>% log2() %>% rowSums()
# 
# ggplot(bardy_df %>% 
#          dplyr::select(one_of('Astro_sum', 'num_genes_exon_tpm_5_resids', 'FTC_simple', 'num_genes_exon_tpm_5')) %>% 
#          filter(FTC_simple > 3), 
#        aes(x = Astro_sum, num_genes_exon_tpm_5_resids, color = FTC_simple)) + geom_point() + geom_smooth(method = 'lm')
# 
# bardy_plot_df$FTC_simple = factor(bardy_plot_df$FTC_simple, levels = c(0, 4, 5))
# 
# bardy_plot_df = merge(bardy_plot_df, bardy_spikes %>% dplyr::rename(sample_id = 'X__1', has_ercc = "Stock Dilution\r\nof ERCC RNA Spike-In Mix Used"    ) )
# 
# 
# df = read_excel('~/Genes vs cell size.xlsx')
# cap_d = df[df$X__1 == 'Capacitance_Meanfit[pF]', ]
# cap_d = cap_d[, -1] %>% as.data.frame()
# rownames(cap_d) = 'cap'
# cap_d = cap_d %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
# 
# bardy_plot_df = left_join(bardy_plot_df, cap_d)
# 
# ggplot(bardy_plot_df , 
#        aes(x = contam_type, y = Astro_sum,  group = contam_type)) + geom_violin() + geom_jitter(width = .1) + 
#   ylab('Astrocyte marker expr \n(log2 TPM+1)') + xlab('Firing phenotype class')
# 
# 
# 
# p0 =  ggplot(bardy_plot_df[bardy_plot_df$has_ercc & bardy_plot_df$FTC_simple %in% c(5), ], 
#              aes(x = read_count, y = num_genes)) +geom_point() + geom_smooth(method = "lm", se = F) +
#   ylab('Detected genes') + xlab('Library size (sequenced reads)')
# p1 = ggplot(bardy_plot_df[bardy_plot_df$has_ercc & bardy_plot_df$FTC_simple %in% c(5), ],
#             aes(x = Astro_sum, y = num_genes)) +geom_point() + geom_smooth(method = "lm", se = F) + 
#   ylab('Detected genes') + xlab('Astrocyte marker sum')
# 
# p2 = ggplot(bardy_plot_df[bardy_plot_df$has_ercc & bardy_plot_df$FTC_simple %in% c(5), ], 
#             aes(x = ercc_pct, y = num_genes)) +geom_point() + geom_smooth(method = "lm", se = F) + 
#   ylab('Detected genes') + xlab('Spike-in mRNA ratio (%)')
# 
# plot_grid(p0, p2, p1, ncol = 3)
# 
# 
# 
# ### model capacitance and num genes after regressing out confounds
# 
# mod = lm(scale(num_genes) ~ scale(log10(read_count)) + scale(log10(ercc_pct)) + scale(contam_sum) , data = 
#             bardy_plot_df[bardy_plot_df$FTC_simple %in% c(1, 2, 3, 4, 5) & bardy_plot_df$has_ercc == T, ], na.action = na.omit)
# summary(mod)
# 
# 
# afss <- anova(mod)$"Sum Sq"
# rsq = summary(mod)$r.squared
# 
# var_exp_df_bardy = data.frame(variable = c('Library size', 'Spike-in ratio',  'Astrocyte contam', 'Capacitance', 'Residual'))
# var_exp_df_bardy = cbind(var_exp_df_bardy, percent_exp =(afss/sum(afss))/rsq*100, var_exp = (afss/sum(afss)) * 100)
# var_exp_df_bardy$dataset = 'Bardy'
# 
# 
# var_exp_comb = rbind(var_exp_df_bardy)
# var_exp_comb$dataset = factor(var_exp_comb$dataset, levels = c('Bardy'))
# # var_exp_comb$dataset = factor(var_exp_comb$dataset, levels = c('Tasic Ndnf', 'Cadwell'))
# 
# # drop residual term
# var_exp_comb = var_exp_comb[var_exp_comb$variable != 'Residual', ]
# 
# var_exp_comb$variable = factor(var_exp_comb$variable, levels = c('Library size', 'Spike-in ratio', 'Astrocyte contam', 'Capacitance'))
# 
# # make a figure that summarizes variance explained
# var_exp_fig = ggplot(var_exp_comb, aes(x = variable, y = var_exp)) + geom_bar(stat = 'identity', position = 'dodge') + 
#   ylab('Var. Explained (norm. %)')+ theme(legend.position="top") + xlab('') 
# # scale_fill_manual("", values = c("Tasic Ndnf" = "grey80", "Cadwell" = "grey30"))
