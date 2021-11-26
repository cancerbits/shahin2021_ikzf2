## ---- read_metadata
clindata <- read_csv(file.path(metadata_dir,"metadata_ikzf2.csv"))
## ---- find_input_path
sample_path <- dir_ls(path=data_dir,
                      recurs=TRUE,
                      type ="directory",
                      glob = "*/COUNT/*") %>%
  as.character() %>% 
  setNames(., basename(.) %>% str_replace("_","-"))
## ---- read_10x
sample_expr <- Read10X(sample_path)
seurat_obj <- CreateSeuratObject(counts = sample_expr,project = project_name)

## ---- read_smart_seq
#bam files
smart_seq_bams <- list.files(glue("{out_dir}/results/pipeline"),
                             recursive=TRUE,
                             pattern = ".Aligned.out_sorted.bam$",full.names = TRUE) %>% 
  setNames(.,basename %>%str_replace(".Aligned.out_sorted.bam$",""))
ref_gtf <- glue("{genome_dir}/Human_Ensembl98_GRCh38.p13/Human_Ensembl98_GRCh38.p13.gtf")
#count reads
smart_seq_counts <- Rsubread::featureCounts(files = smart_seq_bams,
                                            annot.ext = ref_gtf,
                                            isGTFAnnotationFile = TRUE,
                                            GTF.attrType = "gene_name")
#metadata
smart_seq_coldata  <- read_tsv(file.path(metadata_dir,"ikzf2_rna_stats_summary.tsv")) %>% 
  as.data.frame
rownames(smart_seq_coldata) <- smart_seq_coldata$sample_name
smart_seq_coldata <- smart_seq_coldata %>% 
  mutate(across(where(is.character), as.factor))
#rename samples
colnames(smart_seq_counts$counts) <- str_replace(colnames(smart_seq_counts$counts), "_Aligned.out_sorted.bam", "")
#the same order of samples in count and metadata tables
smart_seq_counts$counts <- smart_seq_counts$counts[,rownames(smart_seq_coldata)]

