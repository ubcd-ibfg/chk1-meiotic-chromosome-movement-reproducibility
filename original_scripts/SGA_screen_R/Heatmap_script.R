# =========================
# LIBRERÍAS
# =========================
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(tibble)

data2 <- read.csv("data_heatmap_completo_replica123_sergio_log2fold.csv", sep = ";")

df <- data2 %>%
  semi_join(results_50, by = "strain") ## DataFrame solo con datos TOP200

# =========================
# 1. FILTRAR SOLO 72h Y 4601
# =========================
df_72h <- df %>%
  dplyr::filter(stringr::str_detect(label, "72h")) %>%
  dplyr::filter(stringr::str_detect(label, "4601"))

# =========================
# 2. LIMPIAR ETIQUETAS
# =========================
df_72h <- df_72h %>%
  dplyr::mutate(
    clean_label = stringr::str_replace(label, "^[0-9]+_", "")
  )

df_72h <- df_72h %>%
  mutate(
    gene2 = if_else(
      gene2 %in% c("0", "0/NC"),
      strain,
      gene2
    )
  ) %>%
  mutate(
    gene2 = sub("/.*$", "", gene2)
  )


df_72h_heatmap <- df_72h %>%
  group_by(clean_label, strain, gene2) %>%
  summarise(
    log2fold_value_normalized2 = if (all(is.na(log2fold_value_normalized2))) {
      NA_real_
    } else {
      mean(log2fold_value_normalized2, na.rm = TRUE)
    },
    .groups = "drop"
  )



# =========================
# 3. LONG → WIDE
# =========================
heatmap_df <- df_72h_heatmap %>%
  dplyr::select(gene2, clean_label, log2fold_value_normalized2) %>%
  tidyr::pivot_wider(
    names_from  = clean_label,
    values_from = log2fold_value_normalized2,
    values_fn   = mean
  )

# =========================
# 4. CONVERTIR A MATRIZ
# =========================
heatmap_matrix_subset <- heatmap_df %>%
  tibble::column_to_rownames("gene2") %>%
  as.matrix()

# asegurar que todo sea numérico
heatmap_matrix_subset <- apply(
  heatmap_matrix_subset,
  2,
  as.numeric
)
rownames(heatmap_matrix_subset) <- heatmap_df$gene2
colnames(heatmap_matrix_subset) <- c("72hYES4601_rep1", "72hYES4601_rep2",
                                     "72hYES4601_rep3", "72hYESHU4601_rep1",
                                     "72hYESHU4601_rep2", "72hYESHU4601_rep3")

custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)
color_breaks <- seq(-2, 2, length.out=101)
library(RColorBrewer)

custom_palette <- c(
  "black",                       # para -10
  colorRampPalette(
    brewer.pal(11, "RdBu")
  )(99)
)
# =========================
# 5. HEATMAP
# =========================
pheatmap::pheatmap(
  heatmap_matrix_subset,
  color = custom_palette,
  breaks = color_breaks,
  main = "TOP 50 genes UP",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = colnames(heatmap_matrix_subset),
  labels_row = rownames(heatmap_matrix_subset),
  show_rownames = TRUE,
  fontsize_col = 7,
  fontsize_row = 5,
  angle_col = 45,
  cellwidth = 30,
  cellheight = 5,
  na_col = "black"
)





###MINIPLOT: GENES ACTINA/CITOESQUELETO/MICROTÚBULOS
actin_csv <- read.csv("C:/Users/IBFG/Desktop/DOCTORADO/SGA/PAPER/genes_actin_cytosk_microt_new.csv", sep = ";")
data_actin <- inner_join(data2, actin_csv, by = c("strain" = "ensembl_gene_id"))

# =========================
# 1. FILTRAR SOLO 72h Y 4601
# =========================
df_72h <- data_actin %>%
  dplyr::filter(stringr::str_detect(label, "72h")) %>%
  dplyr::filter(stringr::str_detect(label, "4601"))

# =========================
# 2. LIMPIAR ETIQUETAS
# =========================
df_72h <- df_72h %>%
  dplyr::mutate(
    clean_label = stringr::str_replace(label, "^[0-9]+_", "")
  )

df_72h <- df_72h %>%
  mutate(
    gene2 = if_else(
      gene2 %in% c("0", "0/NC"),
      strain,
      gene2
    )
  ) %>%
  mutate(
    gene2 = sub("/.*$", "", gene2)
  )


df_72h_heatmap <- df_72h %>%
  group_by(clean_label, strain, gene2) %>%
  summarise(
    log2fold_value_normalized2 = if (all(is.na(log2fold_value_normalized2))) {
      NA_real_
    } else {
      mean(log2fold_value_normalized2, na.rm = TRUE)
    },
    .groups = "drop"
  )



# =========================
# 3. LONG → WIDE
# =========================
heatmap_df <- df_72h_heatmap %>%
  dplyr::select(gene2, clean_label, log2fold_value_normalized2) %>%
  tidyr::pivot_wider(
    names_from  = clean_label,
    values_from = log2fold_value_normalized2,
    values_fn   = mean
  )

# =========================
# 4. CONVERTIR A MATRIZ
# =========================
heatmap_matrix_subset <- heatmap_df %>%
  tibble::column_to_rownames("gene2") %>%
  as.matrix()

# asegurar que todo sea numérico
heatmap_matrix_subset <- apply(
  heatmap_matrix_subset,
  2,
  as.numeric
)
rownames(heatmap_matrix_subset) <- heatmap_df$gene2
colnames(heatmap_matrix_subset) <- c("72hYES4601_rep1", "72hYES4601_rep2",
                                     "72hYES4601_rep3", "72hYESHU4601_rep1",
                                     "72hYESHU4601_rep2", "72hYESHU4601_rep3")

custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)
color_breaks <- seq(-2, 2, length.out=101)

# =========================
# 5. HEATMAP
# =========================
pheatmap::pheatmap(
  heatmap_matrix_subset,
  color = custom_palette,
  breaks = color_breaks,
  main = "Actin genes",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = colnames(heatmap_matrix_subset),
  labels_row = rownames(heatmap_matrix_subset),
  show_rownames = TRUE,
  fontsize_col = 7,
  fontsize_row = 6,
  angle_col = 45,
  cellwidth = 30,
  cellheight = 7,
  na_col = "black"
)
