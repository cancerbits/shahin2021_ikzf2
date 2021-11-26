## ---- SMART_seq_split_celltype
#split data based on cell type
samples_lst <- split(as.character(smart_seq_coldata$sample_name),smart_seq_coldata$cell_type_short)
meta_data_lst <- split(smart_seq_coldata,smart_seq_coldata$cell_type_short)
counts_lst <- lapply(samples_lst,function(x){
  smart_seq_counts$counts[,x]
  })

## ---- run_deseq_SMART_seq
#compare_all_biological_replicates_using_pid_as_a_reference
#deseq object
dds_lst <- lapply(names(samples_lst),function(cell_type){
  #subset samples belong to the cell type
  counts <- counts_lst[[cell_type]]
  meta_data <- meta_data_lst[[cell_type]]
  #remove levels not used
  meta_data$bio_rep <- droplevels(meta_data$bio_rep)
  #DESeq object from matrix
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta_data,
                                design = ~ bio_rep)
  #filter low expressing genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  }) %>%
  setNames(.,names(samples_lst))

dds_lst <- lapply(dds_lst, function(cell_type){
  DESeq(cell_type)
  })
## ---- compare_bio_rep
#differential expression analysis analysis results 
res_lst <- lapply(dds_lst, function(cell_type){
  normal_donor <- levels(cell_type$bio_rep) %>% str_subset("nd")
  lapply(normal_donor, function(nd){
    results(cell_type,
            contrast= c("bio_rep","pid",nd))%>%
      as.data.frame() %>% 
      rownames_to_column("gene_symbol") %>% 
      mutate(comparison = glue("pid_vs_{nd}"))
    }) %>%
    bind_rows()
  })%>% 
  setNames(.,names(dds_lst))

#Common, differentially expressed genes
deg_lst <- lapply(res_lst, function(cell_type){
  cell_type %>%
    filter(padj<0.05 & abs(log2FoldChange)>=1) %>%
    group_by(gene_symbol) %>%
    filter(n()>=2,#the gene is differentially expressed in at least two out of the three comparisons
           abs(sum(sign(log2FoldChange)))>=2#the gene maintain the same direction of regulation in at least two out of the three comparisons
           ) %>%
    ungroup() %$%
    unique(gene_symbol)
  })
#transformed expression values
vsd_lst <- lapply(dds_lst, function(cell_type){
  assay(vst(cell_type, blind=FALSE))
  })

## ---- heatmap_of_common_DEGs
#heatmap of common DEGs
ph_lst <- map(names(vsd_lst), function(cell_type){
  degs <- deg_lst[[cell_type]]
  vsd <- vsd_lst[[cell_type]]
  deg_expr <- vsd[degs,]
  ph_col_annot <- meta_data_lst[[cell_type]][, c("genotype","bio_rep")] %>% mutate(across(.fns = as.character))
  ph_annot_color <- list(genotype = c( WT="#525252", `I325V-Hom` = "#fc4e2a"),
                         bio_rep = c(C1 = "#d9d9d9", C5= "#969696",C6= "#737373", P = "#fd8d3c"))
  ph <- pheatmap::pheatmap(deg_expr,
                           scale = "row",
                           color = colorRampPalette(rev(brewer.pal(7, "PiYG")))(50),
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           annotation_col = ph_col_annot,
                           annotation_colors= ph_annot_color,
                           cutree_cols = 2,
                           cutree_rows = 2,
                           cellwidth = 18,
                           main = glue("{toupper(cell_type)}-cell\n DEGs:{length(degs)}"))
  }) %>% 
  setNames(.,names(vsd_lst))

pdf("heatmap_t.pdf")
ph_lst$t
dev.off()
pdf("heatmap_b.pdf")
ph_lst$b
dev.off()
## ---- enricher_analysis
#cut and label branches of hierarchical clustering
b_branches <- cutree(ph_lst$b$tree_row,2)%>% split(names(.),.)
names(b_branches) <- c("up","down")

t_branches <- cutree(ph_lst$t$tree_row,2)%>% split(names(.),.)
names(t_branches) <- c("down","up")

enrich_wrap <- function(gene_set, padj){
  addResp <- enrichrAddList(gene_set)
  #selected databases
  libs <- c("WikiPathways_2019_Human",
            "WikiPathways_2019_Mouse",
            "TRRUST_Transcription_Factors_2019",
            "TRANSFAC_and_JASPAR_PWMs",
            "Reactome_2016",
            "OMIM_Expanded",
            "OMIM_Disease",
            "Panther_2016",
            "NCI-Nature_2016",
            "MGI_Mammalian_Phenotype_2017",
            "MGI_Mammalian_Phenotype_Level_4_2019",
            "KEGG_2019_Human",
            "KEGG_2019_Mouse",
            "KEA_2015",
            "Jensen_TISSUES",
            "Jensen_DISEASES",
            "Jensen_COMPARTMENTS",
            "Human_Phenotype_Ontology",
            "GO_Biological_Process_2018",
            "GO_Cellular_Component_2018",
            "GO_Molecular_Function_2018",
            "ChEA_2016",
            "DisGeNET",
            "BioPlex_2017",
            "BioPlanet_2019",
            "BioCarta_2016"
  )
  res_enrich <- enrichrEnrich(addResp$userListId,
                              libNames = libs)
  # Check if any requests failed
  if (length(res_enrich$failure) > 0) {
    cat(length(res_enrich$failure), " requests failed")
    allLibs <- getEnrichrLibNames()
    # Libraries that failed
    setdiff(allLibs, names(res$success))
  }
  # Filter the data sets
  y <- lapply(res_enrich$success, function(x){
    x %>%
      dplyr::select(Term,Combined.Score,Adjusted.p_value,nOverlapping.Genes,Overlapping.Genes) %>%
      mutate(nOverlapping.Genes_percent = round(as.numeric(nOverlapping.Genes)/length(gene_set),2),
             Combined.Score = as.double(Combined.Score),
             Adjusted.p_value = as.double(Adjusted.p_value),
             nOverlapping.Genes = as.numeric(nOverlapping.Genes)) %>%
      mutate_if(is.logical,as.character)
  }) %>%
    bind_rows(.id = "db") %>%
    group_by(db) %>%
    arrange(desc(Combined.Score)) %>%
    ungroup()
  y
}
#run enricher for both cell types
enrichr_res_b <- lapply(b_branches, enrich_wrap, padj = Inf)
enrichr_res_t <- lapply(t_branches, enrich_wrap, padj = Inf)

## ---- save_enricher_output
walk(names(enrichr_res_b), function(b_regulation){
  write_csv(enrichr_res_b[[b_regulation]], glue("enricher_b_{b_regulation}.csv"))
  })
walk(names(enrichr_res_t), function(t_regulation){
  write_csv(enrichr_res_t[[t_regulation]], glue("enricher_t_{t_regulation}.csv"))
  })

db_of_interest <- c("BioPlanet_2019",
                    "WikiPathways_2019_Human",
                    "KEGG_2019_Human",
                    "GO_Biological_Process_2018",
                    "GO_Molecular_Function_2018")

bind_rows(enrichr_res_b,.id = "regulation") %>%
  filter(db %in% db_of_interest,
         Adjusted.p_value<0.05) %>%
  write_csv(path = glue("enricher_b_sig_db.csv"))

bind_rows(enrichr_res_t,.id = "regulation") %>%
  filter(db %in% db_of_interest,
         Adjusted.p_value<0.05) %>%
  write_csv(path = glue("enricher_t_sig_db.csv"))

## ---- top20_significant_enricher_terms
#subset of GO_Biological_Process
db_of_interest <- c("GO_Biological_Process_2018")

#plot enricher results for b cells
b_enricher_db_lst <- lapply(enrichr_res_b, function(gene_set){
  lapply(db_of_interest, function(database){
    gene_set %>% 
      filter(str_detect(db, database)) %>% 
      as.data.frame() %>% 
      filter(Adjusted.p_value<0.05) %>% #significant enrichments
      arrange(desc(Combined.Score)) %>%
      head(20) %>% # top 20 terms
      mutate(Term = str_replace(Term, "..GO:.+|WP[0-9]+$",""),
             sig_level_short = case_when(
               Adjusted.p_value < 0.0005 ~ "***",
               Adjusted.p_value < 0.005 ~ "**",
               Adjusted.p_value < 0.05 ~ "*",
               TRUE ~ " "
             ),
             sig_level = case_when(
               Adjusted.p_value < 0.0005 ~ "<0.0005",
               Adjusted.p_value < 0.005 ~ "<0.005",
               Adjusted.p_value < 0.05 ~ "<0.05",
               TRUE ~ ">=0.05"
             )) %>% 
      dplyr::select(db,Term, Combined.Score,nOverlapping.Genes,sig_level,sig_level_short)
  })  %>% 
    bind_rows() 
})

b_enricher_plot_df <- bind_rows(b_enricher_db_lst,.id = "regulation") %>%
  arrange(db,ifelse(regulation=="up",-Combined.Score,Combined.Score)) %>%
  mutate(Term = factor(Term, levels = unique(Term)),
         Combined.Score = ifelse(regulation == "down", Combined.Score * -1, Combined.Score),
  ) %>%
  group_by(Term) %>% 
  filter(Combined.Score == max(Combined.Score)) %>% 
  ungroup() %>% 
  split(.,.$db)

pdf("b_enricher_temp.pdf", height = 20,width = 15)
b_enricher_plot <- lapply(names(b_enricher_plot_df), function(database){
  df <- b_enricher_plot_df[[database]]
  db_color <- c("<0.0005" = "#c51b8a",
                "<0.005" = "#fa9fb5",
                "<0.05" = "#fde0dd",
                ">=0.05" = "#bdbdbd")
  df %>% 
    ggplot()+
    geom_col(aes(Term, Combined.Score, fill = sig_level),
             color = "black",
             size = 0.5,
             width = 0.7)+
    geom_text(data = df %>% filter(Combined.Score > 0) %>% mutate(Combined.Score = Combined.Score+0.03),
              aes(Term, Combined.Score, label = sig_level_short),
              size = 5,
              nudge_x = -0.2,
              hjust = 0)+
    geom_text(data = df %>% filter(Combined.Score < 0) %>% mutate(Combined.Score = Combined.Score-0.03),
              aes(Term, Combined.Score, label = sig_level_short),
              size = 5,
              nudge_x = -0.2,
              hjust = 1)+
    coord_flip()+
    scale_fill_manual(values = db_color)+
    facet_grid(db+regulation~., space="free", scales="free")+
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 10),
      axis.title.y = element_text(size = 30),
      axis.title.x = element_text(size = 10),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 20))
  
}) %>% setNames(., names(b_enricher_plot_df))
patchwork::wrap_plots(b_enricher_plot)
dev.off()
# save the table used to create b_enricher plot
write_csv(b_enricher_plot_df$GO_Biological_Process_2018, "b_enricher_GO_Biological_Process_2018.csv")

## ---- DESeq_output_both_cell_types
#label genes based on the cell type/s they were detected in
DEG_regulation <- rbind(stack(b_branches),
                        stack(t_branches)) %>%
  group_by(values) %>%
  summarize(regulation = paste(ind,collapse = "_")) %>% #concatenate cell types for genes detected in both
  ungroup() %>% 
  dplyr::rename(gene_symbol=values)

DEG_celltype <- stack(deg_lst) %>% 
  as.data.frame() %>% 
  group_by(values) %>%
  summarize(DEG_group = paste(ind,collapse = "_")) %>% #concatenate cell types for genes detected in both
  ungroup()%>%
  dplyr::rename(gene_symbol=values)
#tidy/long output
deseq_out_tidy <- lapply(names(res_lst),function(cell_type){
  genes <- deg_lst[[cell_type]]
  res_lst[[cell_type]] %>%
    mutate(cell_type = cell_type) %>%
    dplyr::select(gene_symbol, log2FoldChange, padj, comparison, cell_type)
})%>%
  bind_rows()

deseq_out_tidy %>%
  left_join(DEG_celltype) %>%
  left_join(DEG_regulation) %>%
  write_csv("DEseq_output_long.csv")
#untidy output
comp_id_lst <- deseq_out_tidy %>%
  unite("comparison_id", c(comparison, cell_type )) %>%
  split(.,.$comparison_id)

deseq_out_untidy <- lapply(names(comp_id_lst), function(comparison){
  comp_id_lst[[comparison]] %>%
    as.data.frame() %>%
    rename_with(~ paste0(.,"_",comparison), where(is.numeric))%>%
    dplyr::select(-comparison_id)%>% 
    pivot_longer(-gene_symbol, names_to = "cols", values_to = "values")
  })%>%
  bind_rows()%>% 
  pivot_wider( names_from = "cols", values_from = "values")%>%
  left_join(DEG_celltype) %>%
  left_join(DEG_regulation) %>%
  dplyr::select(gene_symbol, DEG_group, regulation, everything())%>%
  arrange(DEG_group, regulation)

deseq_out_untidy %>%
  write_csv("DEseq_output_wide.csv")
