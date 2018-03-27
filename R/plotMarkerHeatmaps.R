library(pheatmap)
library(RColorBrewer)


# generate heatmap showing top n marker expression per cell type per patch seq dataset
plotMarkerHeatmap = function(markerlist, expr_matrix, show_legend = T, show_cell_labels = F){
  
  
  # markerGeneList = c(ndnf_markers %>% head(10), l_affy$Pyramidal %>% head(10), l_affy$Microglia %>% head(10) )
  trimmed_marker_list = lapply(markerlist, function(l) l %>% getValidGenes(., colnames(expr_matrix)) %>% make.names() %>% head(10)) 
  markers_count = lapply(trimmed_marker_list, function(l) length(l)) %>% unlist
  
  all_trimmed_markers = unlist(trimmed_marker_list) %>% as.character()
  # all_trimmed_markers = getValidGenes(all_trimmed_markers)
  
  print(all_trimmed_markers)
  
  expr_mat = expr_matrix[, all_trimmed_markers]

  
  gaps_row = cumsum(markers_count)
  
  annotation_row = data.frame(Markers = factor(rep(names(markers_count),
                                                   markers_count), levels = names(markers_count)))
  rownames(annotation_row) = colnames(expr_mat[, all_trimmed_markers])
  
  # annotation_col = data.frame(CellTypes = factor(scde_ob$joined_df$major_type))
  # rownames(annotation_col) = rownames(expr_mat)
  
  MAXLOGEXPR = 12
  breaksList = seq(0, MAXLOGEXPR, by = 1)
  
  expr_mat[expr_mat > 2^MAXLOGEXPR] = 2^MAXLOGEXPR
  
  expr_mat = t(expr_mat)

  plot_heatmap = pheatmap(expr_mat %>% log2,              
                             cluster_rows=F,
                             cluster_cols=F, gaps_row = gaps_row, annotation_colors = ann_colors,
                             annotation_row = annotation_row, annotation_names_row = F, annotation_legend = F, show_colnames = show_cell_labels,
                             breaks = breaksList, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                             legend = show_legend
  )
  
}

# sets colors for cell types
ann_colors = list(Markers = c(Ndnf_on = "indianred1", Pyramidal = "turquoise", Pyramidal_on = "turquoise", 
                              Microglia = "gray50", Sncg_on = "purple", Astrocyte = "green4", Pvalb = "red", 
                              Inhibitory = "red",
                              Endothelial = "brown", 
                              Oligodendrocyte = "sandybrown",
                              OPC = 'darkorange4'))


contam_plot_types = factor(c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'Pyramidal'), levels = c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'OPC', 'Pyramidal'))

getContamMatrixFromPatchSeq = function(patch_seq_dataset, cell_type_name, contam_show_types = contam_plot_types){
  
  
  dissoc_corr = patch_seq_dataset$joined_df[patch_seq_dataset$joined_df$major_type == cell_type_name, ]$dissoc_corr
  # on_marker_sum = sumExpression(cadwell_expr, ndnf_markers)
  sorted_cell_ids = order(-dissoc_corr)
  contam_mat = patch_seq_dataset$contam_scores %>% filter(str_detect(cell_id, make.names(cell_type_name))) %>% select(contam_show_types %>% unlist %>% as.character) %>% repWZero() %>% t()
  colnames(contam_mat) = patch_seq_dataset$contam_scores %>% filter(str_detect(cell_id, make.names(cell_type_name))) %>% select(cell_id) %>% unlist

  return(contam_mat[, sorted_cell_ids])
}

plotContamHeatmap = function(contam_matrix, show_cell_labels = F){
  
  MAXCONTAM = 1
  breaksList = seq(0, MAXCONTAM, by = .1)
  
  contam_matrix[contam_matrix > MAXCONTAM] = MAXCONTAM
  contam_heatmap = pheatmap(contam_matrix, cluster_rows = F, cluster_cols = F, show_colnames = show_cell_labels,
                            breaks = breaksList, color = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(length(breaksList)))
  return(contam_heatmap)
}

