## FILTRADO DE 200 GENES - SGA

#Set the working directory
setwd("C:/Users/IBFG/Desktop/DOCTORADO/SGA")

# Load the installed packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(writexl)
# 
# data <- read.csv("top100_down100_data.csv", sep = ",")
# data[is.na(data)] <- 0
# 
# genes <- data$X
# indices <- match(genes, data_2$gene2)
# genes_ids <- data_2$strain[indices]
# write.table(
#   genes_ids, file = "top100_down100_IDs.txt", row.names = FALSE,
#   col.names = FALSE, quote = FALSE)


#Load the data

data <- read.csv("data_heatmap_completo_replica123_sergio_log2fold.csv", sep = ";")
data[is.na(data)] <- 0

#Select 72h YES HU 4601 for every replicate
data_A <- data %>%
  filter(
    str_detect(label, "A_Batch1") &
      str_detect(label, "72hYESHU4601"))
data_B <- data %>%
  filter(
    str_detect(label, "B_Batch1") &
      str_detect(label, "72hYESHU4601"))
data_C <- data %>%
  filter(str_detect(label, "C_Batch1") &
           str_detect(label, "72hYESHU4601"))

#Select TOP 200 for every replicate
genes_top_200_A <- data_A %>%
  arrange(desc(log2fold_value_normalized2)) %>%
  slice_head(n = 200)

genes_top_200_B <- data_B %>%
  arrange(desc(log2fold_value_normalized2)) %>%
  slice_head(n = 200)

genes_top_200_C <- data_C %>%
  arrange(desc(log2fold_value_normalized2)) %>%
  slice_head(n = 200)

#EXTRACT TOP500 IDs AND NAMES
genes.id_A <- genes_top_200_A$strain
genes.id_A_no_rep <- unique(genes.id_A)
genes_A <- genes_top_200_A$gene2
genes_A_no_rep <- unique(genes_A)

genes.id_B <- genes_top_200_B$strain
genes.id_B_no_rep <- unique(genes.id_B)
genes_B <- genes_top_200_B$gene2
genes_B_no_rep <- unique(genes_B)

genes.id_C <- genes_top_200_C$strain
genes.id_C_no_rep <- unique(genes.id_C)
genes_C <- genes_top_200_C$gene2
genes_C_no_rep <- unique(genes_C)

#GENERATE INTERSECTIONS IN THE TOP700 OF 2 REPLICATES
interAB <- intersect(genes.id_A_no_rep, genes.id_B_no_rep)
interAC <- intersect(genes.id_A_no_rep, genes.id_C_no_rep)
interBC <- intersect(genes.id_B_no_rep, genes.id_C_no_rep)

genesAB <- intersect(genes_A_no_rep, genes_B_no_rep)
genesAC <- intersect(genes_A_no_rep, genes_C_no_rep)
genesBC <- intersect(genes_B_no_rep, genes_C_no_rep)

genes_finales <- unique(c(interAB, interAC, interBC))
names_finales <- unique(c(genesAB, genesAC, genesBC))
length(genes_finales)

write.table(
  genes_finales, file = "genes_4601_50_UP.txt", row.names = FALSE,
  col.names = FALSE,quote = FALSE)

write.table(
  names_finales, file = "nombres_4601_50_UP.txt", row.names = FALSE,
  col.names = FALSE, quote = FALSE)

#CALCULATE THE AVERAGE OF THE COLONY AREA AND SORT THE LIST
results_A <- data_A %>%
  filter(strain %in% genes_finales) %>%
  dplyr::select(strain, gene2, log2fold_value_normalized2)

results_B <- data_B %>%
  filter(strain %in% genes_finales) %>%
  dplyr::select(strain, gene2, log2fold_value_normalized2)

results_C <- data_C %>%
  filter(strain %in% genes_finales) %>%
  dplyr::select(strain, gene2, log2fold_value_normalized2)

df_total <- bind_rows(results_A, results_B, results_C)
results_total <- df_total %>%
  group_by(strain) %>%
  summarise(area_total = (sum(log2fold_value_normalized2, na.rm = TRUE))/3,
            .groups = "drop") %>%
  arrange(area_total)

results_50 <- results_total[1:50,]
genes_50 <- as.character(results_50$strain)
names_50 <- data_A %>%
  semi_join(results_50, by = "strain") %>%
  pull(gene2) %>%
  as.character()
names_200 <- unique(names_50)

write.table(
  genes_50, file = "genes_4601_50_UP.txt", row.names = FALSE,
  col.names = FALSE,quote = FALSE)

write.table(
  names_50, file = "nombres_4601_50_UP.txt", row.names = FALSE,
  col.names = FALSE, quote = FALSE)

match("ctr1", names_200) #Posición de X en la lista

####FILTRAR GENES DE ACTINA, MICROTÚBULOS Y CITOESQUELETO
####CONSEGUIR TÉRMINOS DE PROCESOS BIOLÓGICOS DE S. POMBE
library(biomaRt)
ensembl_fungi <- useMart("fungi_mart", host = "https://fungi.ensembl.org")
pombe_mart <- useDataset("spombe_eg_gene", mart = ensembl_fungi)
anotaciones <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                    "go_id", "name_1006"),
                     filters = "ensembl_gene_id",
                     values = genes_200,
                     mart = pombe_mart)
###ACTINA
genes_actina <- anotaciones[grep("actin", anotaciones$name_1006, ignore.case = TRUE), ]
genes_actina <- unique(genes_actina[genes_actina$go_id != "", ])
genes_actina <- genes_actina %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarise(
    all_go_ids = paste(unique(go_id), collapse = "; "),
    all_descriptions = paste(unique(name_1006), collapse = "; "),
    .groups = "drop"
  )

###CITOESQUELETO
genes_citoes <- anotaciones[grep("cytoskeleton", anotaciones$name_1006, ignore.case = TRUE), ]
genes_citoes <- unique(genes_citoes[genes_citoes$go_id != "", ])
genes_citoes <- genes_citoes %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarise(
    all_go_ids = paste(unique(go_id), collapse = "; "),
    all_descriptions = paste(unique(name_1006), collapse = "; "),
    .groups = "drop"
  )

###MICROTÚBULOS
genes_microt <- anotaciones[grep("microtubule", anotaciones$name_1006, ignore.case = TRUE), ]
genes_microt <- unique(genes_microt[genes_microt$go_id != "", ])
genes_microt <- genes_microt %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarise(
    all_go_ids = paste(unique(go_id), collapse = "; "),
    all_descriptions = paste(unique(name_1006), collapse = "; "),
    .groups = "drop"
  )

##CONCENTRACIÓN VERTICAL DE LOS DISTINTOS DATAFRAMES
genes_total <- rbind(genes_actina, genes_citoes, genes_microt)
##ELIMINAR DUPLICACIONES DE GENES
genes_total <- genes_total[!duplicated(genes_total$ensembl_gene_id), ]

#AÑADIR A GENES_TOTAL UNA NUEVA COLUMNA CON SU POSICIÓN EN EL TOP200 Y ORDENAR
genes_total$position <- match(genes_total$ensembl_gene_id, results_200$strain)
genes_total <- genes_total %>%
  arrange(position)

write.csv(genes_total, file = "genes_actin_cytosk_microt_new.csv", row.names = FALSE)

#####HEATMAP TOP200
library(pheatmap)
