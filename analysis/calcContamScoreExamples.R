
cadwell_ndnf_contam = patch_seq_datasets$cadwell$contam_scores %>% as.data.frame() %>% mutate(dataset = 'Cadwell') %>% filter(major_type %in% c('eNGC'))

## calculate contam scores for zeisel dataset
zeisel_df = zeiselExprDataDf %>% select(one_of(c('cell_id' )))
rownames(zeisel_df) = zeisel_df$cell_id

zeisel_df = cbind(zeisel_df, zeisel_contam_all_sub)
zeisel_contam_df = normalizeContam(zeisel_df, zeisel_med_exprs, setdiff(compare_cell_types_inh, 'OPC'), replace_zeros = T) %>% 
  t() %>% 
  as.data.frame %>% 
  mutate(dataset = 'Zeisel')

zeisel_contam_df$cell_id = zeisel_df$cell_id
zeisel_contam_df$contam_type = zeisel_df$contam_type

## calculate contam scores for tasic dataset
tasic_df = aibsExprDataDf %>% select(one_of(c('sample_title' )))
rownames(tasic_df) = tasic_df$sample_title

tasic_df = cbind(tasic_df, aibs_contam_all_sub)
tasic_contam_df = normalizeContam(tasic_df, aibs_med_exprs, compare_cell_types_inh, replace_zeros = T) %>% 
  t() %>% 
  as.data.frame %>% 
  mutate(dataset = 'Tasic')

tasic_contam_df$sample_title = tasic_df$sample_title
tasic_contam_df$contam_type = tasic_df$contam_type

ndnf_pyramidal_cell_df = bind_rows(cadwell_ndnf_contam, 
                                   tasic_contam_df %>% filter(contam_type == 'Ndnf_on'), 
                                   zeisel_contam_df %>% filter(contam_type == 'Ndnf_on'))

ndnf_pyramidal_cell_contam_scores_plot = ndnf_pyramidal_cell_df %>% 
  ggplot(aes(x = Pyramidal)) + 
  geom_vline(xintercept = .5, linetype='dashed') + 
  geom_histogram(bins = 15) + 
  ggtitle('Ndnf cells') + 
  xlab('Contamination score,\npyramidal (ratio)') + ylab('Cell count') + 
  scale_x_continuous(breaks=seq(0, 1.1, .25), limits = c(-.1, 1.1)) + 
  facet_wrap(~dataset, ncol = 1, scales = "free_y") 

foldy_sncg_contam = patch_seq_datasets$foldy$contam_scores %>% as.data.frame() %>% mutate(dataset = 'Földy') %>% filter(major_type %in% c('RS-INT'))

sncg_microglia_cell_df = bind_rows(foldy_sncg_contam, 
                                   tasic_contam_df %>% filter(contam_type == 'Sncg_on'), 
                                   zeisel_contam_df %>% filter(contam_type == 'Sncg_on'))

sncg_microglia_cell_contam_scores_plot = sncg_microglia_cell_df %>% ggplot(aes(x = Microglia)) + 
  geom_vline(xintercept = .5, linetype='dashed') + 
  geom_histogram(bins = 15) + 
  xlab('Contamination score,\nmicroglial (ratio)') + ylab('Cell count') + 
  ggtitle('Sncg cells') + 
  scale_x_continuous(breaks=seq(0, 1.1, .25), limits = c(-.1, 1.1)) + 
  facet_wrap(~dataset, ncol = 1, scales = "free_y") 

contam_scores_example_plots = plot_grid(ndnf_pyramidal_cell_contam_scores_plot, sncg_microglia_cell_contam_scores_plot, nrow = 1, align = 'v')
ggsave(filename = 'plots/contam_scores_example_plots.pdf', plot = contam_scores_example_plots, units = "in", height = 5.5, width = 9)


contam_threshold = .5
lapply(patch_seq_datasets, function(dataset) {
  num_contam_cells = dataset$contam_scores %>% 
    filter(Astrocyte > contam_threshold | 
             Pyramidal > contam_threshold | 
             Microglia > contam_threshold | 
             # OPC > contam_threshold |
             Endothelial > contam_threshold | 
             Oligodendrocyte > contam_threshold |
             Inhibitory > contam_threshold
    ) %>% nrow()
  total_cells = dataset$contam_scores %>% nrow()
  print(num_contam_cells)
  print(total_cells)
  # print(num_contam_cells / total_cells)
})

