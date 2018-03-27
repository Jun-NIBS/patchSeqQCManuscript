library(ogbox)

### Load Tasic 2016 data

print('Loading Tasic 2016 data')

# define clusters of cell subtypes
aibsColorDf <- read.csv("data-raw/tasic2016/aibs_color_legend.csv")
source('data-raw/tasic2016/cellColors.R')
library(viridis)
cell.colors = cellColors()

aibsColors = aibsColorDf %>% 
  separate(legend_label, into = c('tg', 'rest'), sep= c("\\s")) %>% 
  select(tg, color) %>% 
  group_by(tg) %>%
  filter(row_number() == 1)
aibsColors = dplyr::rename(aibsColors, sample.name = tg, colors = color)
aibsColors$sample.name = aibsColors$sample.name %>% make.names
aibsColors$colors = as.character(aibsColors$colors)
aibsColors[aibsColors$sample.name == "Rorb",]$colors = "#74add1"
aibsColors[aibsColors$sample.name == "Gad2",]$colors = "#fddcff"
aibsColors[aibsColors$sample.name == "Chat",]$colors = "#F800FF"
aibsColors[aibsColors$sample.name == "Ndnf",]$colors = "#CD6090"
aibsColors[aibsColors$sample.name == "Nkx2",]$colors = "#FF6A6A"
aibsColors[aibsColors$sample.name == "Chrna2",]$colors = "#E54E03"
#aibsColors[aibsColors$sample.name == "Pvalb",]$colors = "#74add1"
aibsColors[aibsColors$sample.name == "Pvalb",]$colors = cell.colors['GabaPV']
aibsColors[aibsColors$sample.name == "Sst",]$colors = cell.colors['GabaSSTReln']
aibsColors[aibsColors$sample.name == "Htr3a",]$colors = cell.colors['GabaVIPReln']
aibsColors[aibsColors$sample.name == "Vip",]$colors = "purple"
aibsColors[aibsColors$sample.name == "Rbp4",]$colors = "turquoise"


aibsColors$colors <- factor(aibsColors$colors)
rownames(aibsColors) = make.names(aibsColors$sample.name)

aibsColorSubDf <- read.csv("data-raw/tasic2016/cluster_metadata.csv")
# vignette label / # geo label
levels(aibsColorSubDf$Tasic_et_al_2016_label) <- gsub("L5 Hsd11b1","L5a Hsd11b1", levels(aibsColorSubDf$Tasic_et_al_2016_label))
levels(aibsColorSubDf$Tasic_et_al_2016_label) <- gsub("Oligo 9630013A20Rik ","Oligo 96_Rik", levels(aibsColorSubDf$Tasic_et_al_2016_label))
levels(aibsColorSubDf$Tasic_et_al_2016_label) <- gsub("Astro Aqp4","Astro Gja1", levels(aibsColorSubDf$Tasic_et_al_2016_label))
levels(aibsColorSubDf$Tasic_et_al_2016_label) <- gsub("L5b Chrna6","L5 Chrna6", levels(aibsColorSubDf$Tasic_et_al_2016_label))
expr_cell_colors = aibsColorSubDf %>% select(Tasic_et_al_2016_label, vignette_color)
colnames(expr_cell_colors) = c('primary_type', 'subtype_colors')


# load expression data
aibsExprMeta <- read.csv("data-raw/tasic2016/GSE71585_Clustering_Results.csv")
# aibsExprData <- read.csv("~/ephys_analysis/aibs_cre_analysis/GSE71585_RefSeq_TPM.csv")

aibs_soft_file = "data-raw/tasic2016/GSE71585_family.soft"
aibs_soft = softParser(aibs_soft_file)
aibs_soft = aibs_soft %>% dplyr::rename(sample_title = '!Sample_title', geo_id ='!Sample_geo_accession')

aibsExprMeta = merge(aibsExprMeta, aibs_soft %>% select(sample_title, geo_id))

aibsExprData <- read.csv("data-raw/tasic2016/tpmMatrix.genes", sep = '\t')

aibsExprData = merge(annot, aibsExprData, by.x = 'ensembl_gene_id', by.y = 'genes')

rownames(aibsExprData) = make.names(aibsExprData$mgi_symbol, unique = T)
aibs_gene_names = rownames(aibsExprData)

gsm_names = colnames(aibsExprData)[str_detect(colnames(aibsExprData), 'GSM')]
aibsExprData = aibsExprData[, gsm_names]

aibsExprDataTrans = aibsExprData %>% t() %>% as.data.frame

aibsExprDataTrans = aibsExprDataTrans + 1

aibsExprDataTrans$geo_id = rownames(aibsExprDataTrans)

aibsExprDataDf = merge(aibsExprMeta, aibsExprDataTrans, by = 'geo_id')
rownames(aibsExprDataDf) = aibsExprData$geo_id

# rename mouse_line as sample.name
aibsExprDataDf = dplyr::rename(aibsExprDataDf, sample.name = mouse_line)
colnames(aibsExprDataDf) = make.names(colnames(aibsExprDataDf))


aibsExprDataDf = aibsExprDataDf %>% filter(pass_qc_checks == 'Y', broad_type != 'Unclassified') 
aibsExprDataDf = merge(aibsExprDataDf, aibsColors, by = "sample.name")
aibsExprDataDf = left_join(aibsExprDataDf, expr_cell_colors, by = "primary_type")


# load raw counts from tasic
tasic_counts = read.csv('data-raw/tasic2016/GSE71585_RefSeq_counts.csv')
tasic_ercc_counts = read.csv('data-raw/tasic2016/GSE71585_ERCC_and_tdTomato_counts.csv')

read_count = colSums(tasic_counts[-1])
ercc_counts = colSums(tasic_ercc_counts[-1])
num_genes = tasic_counts[-1] %>% apply(., 2, function(x) sum(x > 0))
td_tomato_pct = t(tasic_ercc_counts[tasic_ercc_counts$gene == 'tdTomato', -1]) / read_count 
colnames(td_tomato_pct) = 'td_tomato_pct'

# create a data frame
tasic_df = data.frame(sample_title = names(read_count), num_genes, read_count, ercc_counts, ercc_pct = ercc_counts/read_count, td_tomato_pct = td_tomato_pct)

aibsExprDataDf = left_join(aibsExprDataDf , tasic_df, by = 'sample_title')

detach("package:dplyr")
library("dplyr")

saveRDS(aibsExprDataDf, file = 'data/tasic_expression.rda')
saveRDS(aibs_gene_names, file = 'data/tasic_gene_names.rda')
saveRDS(expr_cell_colors, file = 'data/tasic_expr_cell_colors.rda')
saveRDS(aibsColors, file = 'data/aibs_colors.rda')

