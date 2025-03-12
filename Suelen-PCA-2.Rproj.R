rm(list = ls())     # clear objects  
graphics.off() 
#######################################
###### Suelen PCA 2 ############
#######################################


# Packages ----------------------------------------------------------------
inst <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("tidyverse","cluster", "factoextra","NbClust","tidyr", 
              "ggplot2", "ggpubr", "broom", "AICcmodavg", "ggcorrplot", 
              "fpc","cluster", "readxl", "magrittr","hrbrthemes",
              "multipanelfigure","klaR","psych","MASS","ggord","devtools",
              "reshape2","RColorBrewer","SensoMineR","FactoMineR","stats",
              "dplyr","writexl","gtools","ggbiplot","ggrepel","pheatmap", 
              "ggcorrplot", "CCA")
inst(packages)
theme_set(theme_minimal())


# -------------------------- #
#    Importing Dataset       #
# -------------------------- #
(df <- read_excel("Suelen-PCA-2.xlsx"))


# -------------------------- #
#    Phenolics               #
# -------------------------- #
pc <- df %>%
  select(!c(Samples, `Beta-glucans`, `Total phenolics`, `Antioxidant activity`,
            AGS, Caco2, `MCF-7`, `NCI-H460`, VERO,
            Asp, Glu, Ser, Gly, `His*`, Tau, Arg, `Thr*`, Ala, Pro, Tyr,
            `Val*`, `Met*`, Cys, `Ile*`, `Leu*`, `Phe*`, `Lys*`, Hyp, `Trp*`, `Total aa`))

# Convert to matrix and normalize data
matrix_data <- as.matrix(pc)
rownames(matrix_data) <- df$Samples  # Set sample names as row labels

# Min-Max function to normalization (0 to 3)
scale_minmax <- function(x, min_val = 0, max_val = 1) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max_val - min_val) + min_val)
}

# Normalization of data
matrix_scaled <- apply(matrix_data, 2, scale_minmax)

pheatmap(matrix_data,  
         clustering_distance_rows = "euclidean",    
         clustering_distance_cols = "euclidean",  
         clustering_method = "complete",
         main = "Heatmap of Phenolic Compounds (mg/g)",  
         color = colorRampPalette(c("#ECEFF1", "#78909C"))(50),
         display_numbers = TRUE,   
         angle_col = 45,
         fontsize_number = 8,  # Ajusta el tamaño de los números
         number_color = "#000000", # Color principal del número
         number_format = "%.2f") # Formato para mostrar solo dos decimales




# Perform PCA
pc_df <- as.data.frame(pc)
rownames(pc_df) <- df$Samples
pca_result <- PCA(pc_df, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(pca_result, 
                repel = TRUE, # Evita superposición de etiquetas
                col.ind = "#263239", # Color de individuos (muestras)
                col.var = "#D35400",  # Color de variables (compuestos fenólicos)
                labelsize = 5) + 
  labs(title = "Biplot PCA of Phenolic Compounds")

fviz_pca_var(pca_result, col.var = "contrib")


# Extract principal coordinates and merge with sample names
pca_df <- as.data.frame(pca_result$ind$coord)
pca_df$Sample <- df$Samples

# Plot PCA
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(aes(color = Sample), size = 4) +
  geom_text(vjust = 1.5, hjust = 1) +
  theme_minimal() +
  labs(title = "PCA of Phenolic Compounds", x = "Dimension 1", y = "Dimension 2")


# Convert data to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = -Samples, names_to = "Compound", values_to = "Value")

# Boxplot by compound and sample
ggplot(df_long, aes(x = Samples, y = Value, fill = Samples)) +
  geom_boxplot() +
  facet_wrap(~Compound, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Phenolic Compounds by Sample", 
       y = "Concentration", x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



# -------------------------- #
#    aminoacids              #
# -------------------------- #
aa <- df %>%
  select(!c(Samples, `Beta-glucans`, `Total phenolics`, `Antioxidant activity`,
            AGS, Caco2, `MCF-7`, `NCI-H460`, VERO,
            `Protocatechuic acid`, `Gallic acid`, `Chlorogenic acid`,     
            `4-Hydroxybenzoic acid`, `(-)-Epicatechin`, `Vanillic acid`, 
            `Sinapic acid`, `(+)-Catechin`, `Ellagic acid`, `Ferulic aacid`,         
            `Quinic acid`, `Syringic acid`, `Total aa`, Tau, Hyp)) %>%
  slice(-c(1,2))

# Convert to matrix and normalize data
matrix_data_aa <- as.matrix(aa)
rownames(matrix_data_aa) <- df$Samples[-c(1,2)]  # Set sample names as row labels

# Min-Max function to normalization (0 to 3)
scale_minmax <- function(x, min_val = 0, max_val = 1) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max_val - min_val) + min_val)
}

# Normalization of data
matrix_scaled_aa <- apply(matrix_data_aa, 2, scale_minmax)

# Generate heatmap
pheatmap(t(matrix_data_aa), 
         clustering_distance_rows = "euclidean",    
         clustering_distance_cols = "euclidean",  
         clustering_method = "complete",
         main = "Hierarchical Heatmap of Amino Acids Content (%)",   
         color = colorRampPalette(c("#ECEFF1", "#78909C"))(50),
         display_numbers = TRUE,   
         angle_col = 45,
         fontsize_number = 8,
         number_color = "#000000", 
         number_format = "%.2f") 


# Perform PCA
aa_df <- as.data.frame(aa)
rownames(aa_df) <- df$Samples[-c(1,2)]
pca_result_aa <- PCA(aa_df, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(pca_result_aa, 
                repel = TRUE, # Evita superposición de etiquetas
                col.ind = "#263239", # Color de individuos (muestras)
                col.var = "#D35400",  # Color de variables (compuestos fenólicos)
                labelsize = 5) + 
  labs(title = "Biplot PCA of Amino acids")

fviz_pca_var(pca_result_aa, col.var = "contrib")


# Extract principal coordinates and merge with sample names
pca_df_aa <- as.data.frame(pca_result_aa$ind$coord)
pca_df_aa$Sample <- df$Samples

# Plot PCA
ggplot(pca_df_aa, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(aes(color = Sample), size = 4) +
  geom_text(vjust = 1.5, hjust = 1) +
  theme_minimal() +
  labs(title = "PCA of Amino acids", x = "Dimension 1", y = "Dimension 2")



# -------------------------- #
#    betaglucans             #
# -------------------------- #
bg <- df %>%
  select(!c(Samples,
            `Protocatechuic acid`, `Gallic acid`, `Chlorogenic acid`,     
            `4-Hydroxybenzoic acid`, `(-)-Epicatechin`, `Vanillic acid`, 
            `Sinapic acid`, `(+)-Catechin`, `Ellagic acid`, `Ferulic aacid`,         
            `Quinic acid`, `Syringic acid`,
            Asp, Glu, Ser, Gly, `His*`, Tau, Arg, `Thr*`, Ala, Pro, Tyr,
            `Val*`, `Met*`, Cys, `Ile*`, `Leu*`, `Phe*`, `Lys*`, Hyp, `Trp*`, `Total aa`))

# Convert to matrix and normalize data
matrix_data_bg <- as.matrix(bg)
rownames(matrix_data_bg) <- df$Samples  # Set sample names as row labels

# Min-Max function to normalization (0 to 3)
scale_minmax <- function(x, min_val = 0, max_val = 1) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max_val - min_val) + min_val)
}

# Normalization of data
matrix_scaled_bg <- apply(matrix_data_bg, 2, scale_minmax)

pheatmap(matrix_scaled_bg, 
         clustering_distance_rows = "euclidean",    
         clustering_distance_cols = "euclidean",  
         clustering_method = "complete",
         main = "Heatmap of Global Responses",    
         color = colorRampPalette(c("#ECEFF1", "#78909C"))(50),
         display_numbers = TRUE,   
         angle_col = 45,
         fontsize_number = 8,
         number_color = "#000000", 
         number_format = "%.2f") 



  # Perform PCA
bg_df <- as.data.frame(bg)
rownames(bg_df) <- df$Samples
pca_result_bg <- PCA(bg_df, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(pca_result_bg, 
                repel = TRUE, # Evita superposición de etiquetas
                col.ind = "#263239", # Color de individuos (muestras)
                col.var = "#D35400",  # Color de variables (compuestos fenólicos)
                labelsize = 5) + 
  labs(title = "Biplot PCA of Global Responses")

fviz_pca_var(pca_result_bg, col.var = "contrib")


# Extract principal coordinates and merge with sample names
pca_df_bg <- as.data.frame(pca_result_bg$ind$coord)
pca_df_bg$Sample <- df$Samples

# Plot PCA
ggplot(pca_df_bg, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(aes(color = Sample), size = 4) +
  geom_text(vjust = 1.5, hjust = 1) +
  theme_minimal() +
  labs(title = "PCA of activities", x = "Dimension 1", y = "Dimension 2")



# -------------------------- #
#    Correlations            #
# -------------------------- #
data_cor <- df[, c("Beta-glucans", "Total phenolics", 
                   "Antioxidant activity", "AGS", "Caco2", "MCF-7", 
                   "NCI-H460", "VERO")]

# Calcular la matriz de correlación
cor_matrix <- cor(data_cor, method = "pearson")

# Visualizar con heatmap
ggcorrplot(cor_matrix, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           colors = c("#D35400", "#ECEFF1", "#78909C")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1)) +
  labs(x = NULL, y = NULL)


# -------------------------- #
#    Cell IS                 #
# -------------------------- #
# Calcular el Índice de Selectividad (IS) para cada muestra
df$IS_AGS <- df$VERO / df$AGS
df$IS_Caco2 <- df$VERO / df$Caco2
df$IS_MCF7 <- df$VERO / df$`MCF-7`
df$IS_NCIH460 <- df$VERO / df$`NCI-H460`

# Ver los resultados
df[, c("Samples", "IS_AGS", "IS_Caco2", "IS_MCF7", "IS_NCIH460")]

# Reorganizar el dataframe para formato largo (long format)
df_long <- df %>%
  select(Samples, IS_AGS, IS_Caco2, IS_MCF7, IS_NCIH460) %>%
  pivot_longer(cols = -Samples, names_to = "Cell_Line", values_to = "IS")

# Crear el heatmap
ggplot(df_long, aes(x =Cell_Line, y =Samples, fill = IS)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(IS, 2)), size = 3, fontface = "bold") +
  scale_fill_gradientn(colors = c("#78909C", "#ECEFF1", "#D35400"), 
                       values = scales::rescale(c(0, 1, 3, max(df_long$IS, na.rm = TRUE))), 
                       name = "Selectivity\nIndex") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Selectivity Index Heatmap", 
       x = "Cell Line", 
       y = "Samples")


# -------------------------- #
#    Correlations TPC        #
# -------------------------- #
data_cor_tpc <- df[, c("Protocatechuic acid", "Gallic acid", "Chlorogenic acid",     
                   "4-Hydroxybenzoic acid", "(-)-Epicatechin", "Vanillic acid", 
                   "Sinapic acid", "(+)-Catechin", "Ellagic acid", "Ferulic aacid",         
                   "Quinic acid", "Syringic acid", "Antioxidant activity")]

# Calcular la matriz de correlación
cor_matrix_tpc <- cor(data_cor_tpc, method = "pearson")

# Visualizar con heatmap
ggcorrplot(cor_matrix_tpc, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           colors = c("#D35400", "#ECEFF1", "#78909C")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1)) +
  labs(x = NULL, y = NULL)
