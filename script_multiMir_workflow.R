params <- list(
     circ = "your_value", # Change this as needed
     l = 10000, #changed to 10000
     QUANTILE1 = "FALSE",
     thr1 = 0.95,
     score_miRNA = 140,
     energy_miRNA = -20,
     QUANTILE2 = "FALSE",
     thr2 = 0.95,
     dGduplex_miRNA = -15,
     dGopen_miRNA = -15,
     QUANTILE3 = "FALSE",
     thr3 = 0.9,
     voteFrac_RBP = 0.15,
     orgdb = "org.Rn.eg.db",
     meshdb = "MeSH.Rn.eg.db",
     symbol2eg = "org.Rn.egSYMBOL2EG",
     eg2uniprot = "org.Rn.egUNIPROT",
     org = "rnorvegicus"
 )
library("knitr")

library(BiocManager)

library(data.table)
library(devtools)
library(dplyr)
library(DT)
library(ggplot2)
library(gridExtra)
library(gtable)
library(igraph)
library(lattice)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(rmdformats)
library(tibble)
library(viridis)
library(grid)
library(easyGgplot2)
library(multiMiR)


#adjust as appropriate
setwd("")

#########input files##########

# circRNA length
file_circ_length <- "backsplice_circRNA_length_1.txt"

# backsplice gene names
file_gene_names <- "backsplice_gene_names.txt"

# file of parameters
file_parameters <- "params_R.txt"

# miRanda output
file_miRanda <- "output_miRanda_per_R.txt"

# probe name file
file_probe_name <- "probe_names.csv"




parameters <- read.table(file_parameters, header=F)
parameters <- as.data.table(parameters)


## miRNA prediction

# miRanda
miRanda_pred <- read.table(file_miRanda, header=T, sep="\t")
miRanda_pred <- as.data.table(miRanda_pred)
miRanda_pred$circ_group <- sub("::.*", "", miRanda_pred$circ_id)


modify_string <- function(string) {
  pattern <- "([^:]+):(\\d+)-(\\d+)"
  match <- regexpr(pattern, string)
  
  if (match[1] != -1) {  # Check if there is a match
    components <- regmatches(string, match)
    
    prefix <- sub(":.*", "", components)  # Extract the prefix
    first_int <- as.numeric(sub(".*:(\\d+)-.*", "\\1", components)) - 1  # Subtract 1 from the first integer
    second_int <- sub(".*-(\\d+)", "\\1", components)  # Extract the second integer
    
    # Reconstruct the new string
    new_string <- paste0(prefix, ":", first_int, "-", second_int)
    return(new_string)
  } else {
    return(string)  # Return the original string if no match
  }
}

##read probe names
probe_names <- read_delim(file_probe_name, 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

probe_names <- probe_names[,-1]

circ_length <- read.table(file_circ_length, header=T, sep="\t")
circ_length <- as.data.table(circ_length)


gene_names <- read.table(file_gene_names, header=T, sep="\t")
gene_names <- as.data.table(gene_names)

gene_names <- merge(gene_names, probe_names, by.x = "gene_names", by.y = "probeID")

gene_names <- gene_names %>%
  mutate(new_circ_group = sapply(circ_id, modify_string))

miRanda_pred <-left_join(miRanda_pred, gene_names, by = c("circ_group" = "new_circ_group"), suffix = c("", ".y"))

miRanda_pred$GeneSymbol <- paste0("circ",miRanda_pred$GeneSymbol)

output_dir <- getwd()

# remove circRNAs with length > "l" and < 60
circ_l <- circ_length[length < params$l & length > 60, circ_id ]

setkey(gene_names, circ_id)
setkey(circ_length, circ_id)
circ_info <- merge(gene_names, circ_length)


miRanda_pred_l <- miRanda_pred[circ_id %in% circ_l,]

# filter on "score" ed "energy"
if (params$QUANTILE1 == "TRUE") {
  score_miRNA <- quantile(miRanda_pred_l$score, params$thr1)
  energy_miRNA <- quantile(miRanda_pred_l$energy, 1-params$thr1)
} else {
  score_miRNA <- params$score_miRNA
  energy_miRNA <- params$energy_miRNA
}

miRanda_pred_l_sel <- miRanda_pred_l[(score > score_miRNA & energy < energy_miRNA),]

miRanda_pred_l_sel$circ_miRNA_end_id <- paste0(miRanda_pred_l_sel$circ_id, "__", miRanda_pred_l_sel$miRNA_id, "__", miRanda_pred_l_sel$circ_end)
miRNA_pred <- miRanda_pred_l_sel

miRNA_pred <- data.table("circ_id" = miRNA_pred$circ_id, "miRNA_id" = miRNA_pred$miRNA_id, "putative_circ_start" = miRNA_pred$putative_circ_start, "circ_end" = miRNA_pred$circ_end, "circ_id" = miRNA_pred$circ_id, "miRanda_score" = miRNA_pred$score, "miRanda_energy" = miRNA_pred$energy, "gene_symbol" = miRNA_pred$GeneSymbol)

miRNA_pred$gene_circ_names <- paste0("circ", miRNA_pred$gene_names, "_", miRNA_pred$circ_id)

bar1a <- miRNA_pred[, .N, by = circ_id]

bar1a$gene_names <- gsub("circ", "", bar1a$gene_circ_names)
bar1a$gene_names <- gsub("_.*$", "", bar1a$gene_names)
bar1a$circ_id <- gsub("^.*_", "", bar1a$gene_circ_names)

bar1a2 <- miRNA_pred[, .(length(unique(miRNA_id))), by = circ_id]
setkey(bar1a, gene_circ_names)
setkey(bar1a2, gene_circ_names)
bar1a <- bar1a[bar1a2, nomatch = 0]
colnames(bar1a)[5] <- "N_unique"

setkey(bar1a, circ_id)
setkey(circ_length, circ_id)
table1a <- merge(bar1a, circ_length)
table1a <- table1a[, c(1,4,3,5,6)]

datatable(data = table1a,
          rownames = F,
          colnames = c("CircRNA", "Gene name", "Number of miRNA binding sites", "Number of miRNAs", "CircRNA length"),
          style = "bootstrap",
          class = "compact display",
          fillContainer = F,
          autoHideNavigation = T,
          filter = "top",
          escape = FALSE,
          extensions = c('ColReorder', 'Buttons', 'KeyTable'),
          options = list(colReorder = TRUE,
                         searching = TRUE,
                         pageLength = 10,
                         dom = 'Bfrtip',
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                         keys = TRUE
          )
)

bar1b <- miRNA_pred[, .N, by = .(circ_id, miRNA_id, miRanda_score, miRanda_energy, gene_symbol)]

setkey(bar1b, circ_id)
setkey(circ_length, circ_id)
bar1c <- bar1b[circ_length, nomatch = 0]

colnames(bar1c)[3] <- "num_siti"

bar1c <- data.table(bar1c, "density" = bar1c$num_siti/bar1c$length*100)

setkey(bar1b, circ_id)
setkey(gene_names, circ_id)
bar1b <- bar1b[gene_names, nomatch = 0]

bar1b$gene_circ_names <- paste0("circ", bar1b$gene_names, "_", bar1b$circ_id)

bar1b$circ_group <- sub("::.*", "", bar1b$circ_id)

if (dim(bar1b)[1] > 0) {
  num_miRNA_per_circRNA <- as.data.table(dcast(bar1b, miRNA_id ~ gene_symbol, value.var = "N", fill = 0))
}

setkey(bar1c, circ_id)
setkey(gene_names, circ_id)

####miRNA pred counts###########
miRNA_pred_binding_sites_filtered <- miRNA_pred


miRNA_pred_binding_sites_filtered$circ_group <- sub("::.*", "", miRNA_pred_binding_sites_filtered$circ_id)
unique_mirna_pred <- as.data.frame(unique(miRNA_pred_binding_sites_filtered$miRNA_id))
colnames(unique_mirna_pred) <- "miRNA_id"
unique_mirna_pred$num_of_circ_ids_per_mirna <- 0

for(i in 1:length(unique_mirna_pred[,1])) {
  mirna_of_interest <- miRNA_pred_binding_sites_filtered[i,]$miRNA_id
  num_of_circ_ids <- length(unique(miRNA_pred_binding_sites_filtered[miRNA_id == mirna_of_interest]$circ_group))
  unique_mirna_pred[i,]$num_of_circ_ids_per_mirna <- num_of_circ_ids
}


merged_df <- merge(bar1c, unique_mirna_pred, by = "miRNA_id", all = TRUE)

bar1c <- merged_df

if (parameters[2] == "hsa" | parameters[2] == "mmu" | parameters[2] == "rno") {
  
  ## multiMir
  
  if (dim(bar1b)[1] > 0) {
    
    multimir_results_all <- data.table()
    
    for (circ in circ_length$circ_id) {
      miRNA_id <- bar1c[circ_id == circ,]$miRNA_id
      
      if (length(miRNA_id) > 0) {
        print(length(miRNA_id))
      }
    }
    
    for (circ in circ_length$circ_id) {
      miRNA_id <- bar1c[circ_id == circ,]$miRNA_id
      print(length(miRNA_id))
      if (length(miRNA_id) > 0) {
        
        miRanda_info <- bar1b[circ_id == circ,]
        filtered_df <- miRanda_info  %>%
          group_by(miRNA_id) %>%
          slice(which.max(miRanda_score)) %>%
          ungroup()
        
        filtered_df <- filtered_df %>%
          rename(mature_mirna_id = miRNA_id)
        
        specie <- case_when(
          parameters[2] == "hsa" ~ "hsa",
          parameters[2] == "mmu" ~ "mmu",
          parameters[2] == "rno" ~ "rno"
        )
        
        multimir_results <- get_multimir(org = specie, mirna = miRNA_id, table = "validated", predicted.site = "miranda")
        multimir_results <- as.data.table(multimir_results@data)
        
        if(multimir_results[, .N] > 0){
          merged_df <- merge(multimir_results, filtered_df, by = c("mature_mirna_id"))
          multimir_results <- merged_df
        }
        
        if (multimir_results[, .N] != 0) {
          multimir_results <- data.table("circ_id" = circ, "miRNA_id" = multimir_results$mature_mirna_id, "target_gene" = multimir_results$target_symbol)
          multimir_results_all <- rbind(multimir_results_all, multimir_results)
          
        }
        
      }
      
    }
    
    if (multimir_results_all[, .N] > 0) {
      
      multimir_results_all <- unique(multimir_results_all)
      multimir_results_all <- multimir_results_all[target_gene != "",]
      
      setkey(multimir_results_all, circ_id)
      setkey(gene_names, circ_id)
      multimir_results_all <- multimir_results_all[gene_names, nomatch = 0]
      
      multimir_results_all <- multimir_results_all[, c(1,4,2,3)]
      
      # Save table
      write.table(multimir_results_all, paste0(output_dir, "All_validated_TGs.csv"), sep = ",", row.names = F, quote = F, dec = ",")
      
    }
    
  }
  
}

merged_df <- merge(multimir_results_all, bar1b[, c("circ_id", "gene_symbol")], by = "circ_id", all.x = TRUE)

bar1c_small <- bar1c[,c("circ_id", "gene_symbol")]
bar1c_small <- unique(bar1c_small)

merged_df <- multimir_results_all %>%
  left_join(bar1c_small, by = "circ_id")

filtered_multimir_results_all <- multimir_results_all#[miranda_score > 155,]
target_genes <- unique(filtered_multimir_results_all$target_gene)

multimir_results_all$circ_group <- sub("::.*", "", multimir_results_all$circ_id)
bar1c$circ_group <- sub("::.*", "", bar1c$circ_id)

eg <- bitr(target_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=params$orgdb)



sessionInfo()



# Our session info
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-redhat-linux-gnu
# Running under: Red Hat Enterprise Linux 9.4 (Plow)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C               LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8     LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
# [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: America/Winnipeg
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.6.5        viridisLite_0.4.2    tibble_3.2.1         RColorBrewer_1.1-3   lattice_0.22-6       gtable_0.3.5         gridExtra_2.3        ggplot2_3.5.1       
# [9] BiocManager_1.30.23  knitr_1.49           GO.db_3.19.1         AnnotationDbi_1.66.0 Biobase_2.64.0       Biostrings_2.72.1    GenomeInfoDb_1.40.1  XVector_0.44.0      
# [17] IRanges_2.38.0       S4Vectors_0.42.0     BiocGenerics_0.50.0  stringr_1.5.1        dplyr_1.1.4          readr_2.1.5          multiMiR_1.26.0      data.table_1.16.4   
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.44.1         xfun_0.50               tzdb_0.4.0              vctrs_0.6.5             tools_4.4.3             bitops_1.0-9            generics_0.1.3         
# [8] fansi_1.0.6             RSQLite_2.3.11          blob_1.2.4              pkgconfig_2.0.3         lifecycle_1.0.4         GenomeInfoDbData_1.2.12 compiler_4.4.3         
# [15] munsell_0.5.1           htmltools_0.5.8.1       RCurl_1.98-1.17         pillar_1.9.0            crayon_1.5.3            cachem_1.1.0            tidyselect_1.2.1       
# [22] digest_0.6.37           stringi_1.8.4           purrr_1.0.2             fastmap_1.2.0           colorspace_2.1-0        cli_3.6.3               magrittr_2.0.3         
# [29] XML_3.99-0.18           utf8_1.2.4              withr_3.0.0             UCSC.utils_1.0.0        scales_1.3.0            bit64_4.6.0-1           httr_1.4.7             
# [36] bit_4.5.0.1             png_0.1-8               hms_1.1.3               evaluate_1.0.3          memoise_2.0.1           rlang_1.1.4             glue_1.7.0             
# [43] DBI_1.2.3               rstudioapi_0.17.1       jsonlite_1.8.8          R6_2.5.1                zlibbioc_1.50.0     
