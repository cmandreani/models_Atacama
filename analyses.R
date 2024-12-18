# Multivariate analyses
# 18 Dec 2024

#Load global packages (needed in almost every section)
library(ggplot2)
library(vegan)

#### 1) BACKGROUND OF THE TLT ####

#       1.1 Environmental metadata ####

#Load data
metadata <- read.table("MX_environm_metadata.csv", header = TRUE, dec = ".",sep = ";")

#Format data
rownames(metadata) <- metadata[,1]
metadata_real_S1.S6 <- as.data.frame(scale(metadata[,-1])) #all sites
metadata_real_S1.S5 <- as.data.frame(scale(metadata[-nrow(metadata),-1])) #exclude S6 for Fig. S1

#Extract names of sites (n=6) and environmental variables (n=17)
env_vars <- colnames(metadata_real_S1.S6)
sites_6 <- rownames(metadata_real_S1.S6)
sites_5 <- rownames(metadata_real_S1.S5)

#Define colors of sites for visualization
site_colorBlind_6 <- c("#FFAD65", "#FF6E00", "#00FCF6", "#00BCB8", "#006666", "#4B0082")
site_colorBlind_5 <- site_colorBlind_6[-length(site_colorBlind_6)]
colors_6 <- setNames(site_colorBlind_6, sites_6)
colors_5 <- setNames(site_colorBlind_5, sites_5)

#Compute PCAs
pca_scaled_6_sites <- prcomp(metadata_real_S1.S6, scale = T)
pca_scaled_5_sites <- prcomp(metadata_real_S1.S5, scale = T)

#Explore projections
pca.summ_6_sites <- summary(pca_scaled_6_sites)
pca.summ_5_sites <- summary(pca_scaled_5_sites)

library(factoextra)
 
#Optional
fviz_eig(pca_scaled_6_sites, addlabels = TRUE, ylim = c(0, 60), geom = "bar",
         barfill = "lightgrey", barcolor = "lightgrey", linecolor = "black")
fviz_eig(pca_scaled_5_sites, addlabels = TRUE, ylim = c(0, 60), geom = "bar",
         barfill = "lightgrey", barcolor = "lightgrey", linecolor = "black")

#Extract PC1 and PC2
dimensions_6 <- as.data.frame(pca_scaled_6_sites$x)[1:2]
dimensions_5 <- as.data.frame(pca_scaled_5_sites$x)[1:2]

#Visualize PCA all sites (Fig. 1A)
fviz_pca_biplot(pca_scaled_6_sites, 
                repel = TRUE, 
                axes = c(1, 2),
                label = c("ind", "var"),
                col.ind = "white", 
                col.var = "black",
                geom.ind = c("text", "point"), 
                geom.var = c("text", "point"),
                pointsize = 1.4) +
  geom_text(aes(label = rownames(dimensions_6), color = sites_6,
                size = 4.3, fontface = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "PC1 (55.9%)", y = "PC2 (24.7%)", title = "Soil metadata") +
  theme_minimal() +
  scale_color_manual(values = colors_6) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        axis.text.x = element_text(margin = margin(t = 10)), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size = 11),
        axis.title.y = element_text(margin = margin(r = 10), size = 11),
        text = element_text(size = 10))

#Visualize PCA excluding S6 (Fig. S1)
fviz_pca_biplot(pca_scaled_5_sites, repel = T, axes = c(1, 2),
                label = c("ind", "var"),
                col.ind = "white", 
                col.var = "black",
                geom.ind = c("text", "point"), 
                geom.var = c("text", "point"),
                pointsize = 1.4) +
  geom_text(aes(label = sites_5, color = sites_5, 
                size = 4.3, fontface = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "PC1 (45.9%)", y = "PC2 (22.8%)", title = "Soil metadata") +
  theme_minimal() +
  scale_color_manual(values = colors_5) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        axis.text.x = element_text(margin = margin(t = 10)), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size = 11),
        axis.title.y = element_text(margin = margin(r = 10), size = 11),
        text = element_text(size = 10))

#Format data for boxplots
metadata_long_6 <- stack(metadata_real_S1.S6)
colnames(metadata_long_6)[colnames(metadata_long_6) == "ind"] <- "environmental_conditions"
colnames(metadata_long_6)[colnames(metadata_long_6) == "values"] <- "scaled_measurements"
metadata_long_6$sites <- rep(sites_6, length(env_vars))

#Visualize boxplots (Fig. 4B)
ggplot(metadata_long_6, aes(x=environmental_conditions, y=scaled_measurements)) +
  geom_boxplot(outliers = T, outlier.size = 5, outlier.alpha = 1, 
               outlier.shape = 1, outlier.color = "black", color = "black") +
  geom_point(aes(colour = sites), size = 3, alpha=0.8) +
  scale_color_manual(values = colors_6) +
  labs(x = "Soil physicochemical parameters measured", y = "Scaled measurements",
       title = "") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11))

#       1.2 Alpha Diversity ####

library(phyloseq)

#Load data
abund <- read.table("raw_abundances_TAXA_metagenomes.csv", sep = ";", dec = ".",
                    header = TRUE)   

#Optional: Extract dictionary of OTUs' taxonomy
phyla_refs <- as.data.frame(abund$ref_phylum)
rownames(phyla_refs) <- abund$ref_ID
phyla_refs <- abund[-nrow(phyla_refs),] #remove 'unassigned'

#Format data
rownames(abund) <- abund$ref_ID
abund <- abund[-nrow(abund),-1] #remove 'unassigned' and ref_ID (stored in rownames)
abund_t <- as.data.frame(t(abund[,-ncol(abund)])) #transpose df removing phyla_refs

#Rarefy counts to the least sequenced sample
phylo_temp <- phyloseq(otu_table(t(abund_t), taxa_are_rows = TRUE))
phylo_rrfy <- rarefy_even_depth(phylo_temp, rngseed = 711)
rrfy_df = as.data.frame(as(otu_table(phylo_rrfy), "matrix"))
colnames(rrfy_df) <- sites_6

#Compute of richness, shannon, evenness, and fisher
abund_rarefied <- rarefy(abund_t, min(rowSums(abund_t)))
shannon <- diversity(t(rrfy_df)) 
even_pielou <- shannon/log(specnumber(t(rrfy_df))) #av : 0.5887846
fisher <- fisher.alpha(t(rrfy_df)) # alpha parameter of Fisher's log-series

#       1.3 PCoA upon taxonomic and functional abundances ####

#Load data
temp <- list.files(pattern="*relative.csv")
rel_abunds <- lapply(setNames(temp, make.names(gsub("*_relative.csv$", "", temp))), 
                     function(f) {
                       read.table(f, header = TRUE, sep = ";", dec = ".", 
                                  stringsAsFactors = T)
                       })

#Format data
sites <- sites_6
rel_abunds2 <- list()

for (i in seq_along(rel_abunds)) {
  df <- rel_abunds[[i]]
  rownames(df) <- df$annot
  df2 <- df[-1:ncol(df)]
  rel_abunds2[[i]] <- df2
  names(rel_abunds2)[i] <- names(rel_abunds)[i]
}

#Compute HELLINGER-TRANSFORMED ABUNDANCES 
fx_sqrt <- function(x) {
  y <- as.data.frame(sqrt(x))
  return(y)
  }
  
hllg_abunds <- lapply(rel_abunds2, fx_sqrt)

#Compute PCoAs
pcoa_coords_hllg <- list()
pcoa_eig_hllg <- list()

for (i in seq_along(hllg_abunds)) {
  df <- hllg_abunds[[i]]
  dist_abund <- vegdist(t(df), method = "bray")
  pcoa_abund <- cmdscale(dist_abund, eig = T, add = T)
  colnames(pcoa_abund$points) <- c("PCo1", "PCo2")
  pcoa_coords_hllg[[i]] <- pcoa_abund$points
  names(pcoa_coords_hllg)[i] <- names(rel_abunds2)[i]
  pcoa_eig_hllg[[i]] <- pcoa_abund$eig
  names(pcoa_eig_hllg)[i] <- names(rel_abunds2)[i]
  print(paste(names(rel_abunds2)[i], ":", round(100 * pcoa_abund$eig / sum(pcoa_abund$eig),2)[1:2]))
}

#Visualize PCoAs (Figs. 1B, 1C, S2, and S3)
pcoa_plots_hllg <- list()

for (i in seq_along(pcoa_coords_hllg)) {
  coords <- pcoa_coords_hllg[[i]]
  pcoa_plot <- ggplot(coords, mapping =aes(x=PCo1, y=PCo2, color = sites)) +
    #geom_point() +
    geom_text(aes(label = rownames(coords)), hjust = 0, vjust = 0, size = 4.3, fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = paste0(names(pcoa_coords_hllg[i])),
         x = paste0("PCo1 (",  round(100 * pcoa_eig_hllg[[i]] / sum(pcoa_eig_hllg[[i]]),1)[1], "%)"),
         y = paste0("PCo2 (",  round(100 * pcoa_eig_hllg[[i]] / sum(pcoa_eig_hllg[[i]]),1)[2], "%)")) +
    scale_color_manual(values = site_colorBlind_6) +
    theme_minimal()+
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white"),
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
          axis.text.x = element_text(margin = margin(t = 10)), 
          axis.text.y = element_text(margin = margin(r = 10)), 
          axis.title.x = element_text(margin = margin(t = 10), size= 11),
          axis.title.y = element_text(margin = margin(r = 10), size= 11),
          text = element_text(size = 10))
  
  plot(pcoa_plot)
  pcoa_plots_hllg[[i]] <- pcoa_plot
  names(pcoa_plots_hllg)[i] <- names(pcoa_coords_hllg)[i]
}

#### 2) METABOLISM OF THE TLT ####

#       2.1 Extraction of scopes for metagenomic data (MetaG-GEMs) ####

#Data upload (output of Menetools)
temp_metaG = list.files(pattern=".txt")

#Extraction of lists
fx_readFile <- function(out_txt) {
  myfile <- readLines(out_txt)
  list1_start <- grep('"produced_seeds": \\[', myfile)
  list2_start <- grep('"scope": \\[', myfile)
  produced_seeds <- trimws(gsub("[\n\",]", "", (myfile[(list1_start + 1):(list2_start - 2)])))
  scope <- trimws(gsub("[\n\",]", "", (myfile[(list2_start + 1):(length(myfile) - 2)])))
  mylists <- list(produced_seeds, scope)
  names(mylists) <- c("produced_seeds", "scope")
  return(mylists)
}

myfiles_metaG = lapply(setNames(temp_metaG, make.names(gsub("*.txt$", "", temp_metaG))), fx_readFile)

#Data mngmnt and download
metaG_data <- list()

for (i in seq_along(myfiles_metaG)) {
  filename_produced_seeds <- paste("producedSeeds_", names(myfiles_metaG)[i], ".csv", sep = "")
  filename_scope <- paste(names(myfiles_metaG)[i], "__metaG.csv", sep = "") 
  df <- as.data.frame(myfiles_metaG[[i]]$scope)
  colnames(df)[colnames(df) == "myfiles_metaG[[i]]$scope"] <- "compound"
  parts <- unlist(strsplit(names(myfiles_metaG)[i], "_"))
  site <- parts[[1]] #Note that file names follow the logic site_seed_dataset, adapt to your needs
  seed <- parts[[2]]
  df$site <- site
  df$seed <- seed
  df$source_list <- "metaG"
  metaG_data[[i]] <- df
  names(metaG_data)[i] <- names(myfiles_metaG)[i]
  write.table(metaG_data[[i]], file = filename_scope, row.names = FALSE, sep = ";")
}

#       2.2 Extraction of scopes for genomic data (MAG-GEMs) ####

#File renaming (output of Mentage2Metabo)
main_dir <- getwd()
subdirectories_p <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE) 
subdirectories_c <- list.dirs(subdirectories_p, full.names = TRUE, recursive = FALSE) #use if you need to extract iscope.tsv

#Move scope files of all simulations into current directory
fx_rename.move <- function(subdir) { 
  folder_name <- basename(subdir)
  files <- list.files(subdir, full.names = TRUE)
  log_file <- files[grep("m2m_metacom", files)]   
  if (length(log_file) > 0) {
    new_log_file <- file.path(main_dir, paste(folder_name, "log.txt", sep = "_"))
    file.rename(log_file, new_log_file)
  }
}

lapply(subdirectories_p, fx_rename.move)

#Data upload
temp_log = list.files(pattern="*_log.txt")
my_logs = lapply(setNames(temp_log, make.names(gsub("*_log.txt$", "", temp_log))), readLines)

#Extraction of lists
my_logs_df <- lapply(my_logs, function(log_file) {
  data_frame <- data.frame(content = log_file, line_number = seq_along(log_file))
  names(data_frame) <- names(log_file)
  return(data_frame)
})

fx_extract.cpds <- function(data_frame) {
  metabolite_rows <- data_frame[grepl("M_[a-zA-Z0-9_-]+_[ce]$", data_frame[,1]), ]
  return(metabolite_rows)
}

metabolite_rows_list <- lapply(my_logs_df, fx_extract.cpds)

fx_split.chunks <- function(df) {
  line_numbers <- df[,2]
  interrupt_indices <- which(diff(as.numeric(line_numbers)) > 1)
  chunks <- split(df, cumsum(seq_along(line_numbers) %in% c(1, interrupt_indices + 1)))
  chunk_names <- c("core", "iscope", "coop", "targeted") #according to the order in the logs
  named_chunks <- setNames(chunks, chunk_names)
  if (length(named_chunks) != 4) {
    warning("Number of chunks is not equal to 4!")
  } 
  return(named_chunks)
}

split_named_chunks_list <- lapply(metabolite_rows_list, fx_split.chunks)

#Data mngmnt and download of scopes
MAGs_data <- list()

for (i in seq_along(split_named_chunks_list)) {
  df <- split_named_chunks_list[[i]]
  coop_df <- df$coop
  iscope_df <- df$iscope
  parts <- unlist(strsplit(names(split_named_chunks_list)[i], "_"))
  prefix <- parts[[1]] #Note that file names follow the logic site_seed_dataset, adapt to your needs
  key <- parts[[2]]
  coop_df$site <- prefix
  coop_df$seed <- key
  coop_df$source_list <- "coop"
  iscope_df$site <- prefix
  iscope_df$seed <- key
  iscope_df$source_list <- "iscope"
  scope_df <- rbind(coop_df, iscope_df)
  colnames(scope_df)[1] <- "compound"
  scope_df[2] <- NULL
  MAGs_data[[i]] <- scope_df
  names(MAGs_data)[i] <- names(split_named_chunks_list)[i]
}

#Download of scope lists

for (i in names(MAGs_data)) {
  filename <- paste0(i, "__MAGs.csv")
  write.table(MAGs_data[[i]], file = filename, sep = ";", row.names = FALSE)
}

# /!\ DONT RUN :  download splitted lists

for (unique_name in names(split_named_chunks_list)) {
  named_chunks <- split_named_chunks_list[[unique_name]]
  chunk_names <- c("core", "iscope", "coop", "targeted") #according to the order in the logs
  
  for (i in seq_along(chunk_names)) {
    chunk_name <- chunk_names[i]
    filename <- paste(paste(unique_name, chunk_name, sep = "_"), ".csv", sep = "")
    data_frame <- named_chunks[[i]]
    write.csv(data_frame[,1], file = filename, row.names = FALSE)
  }
}

#       2.3 Conversion of lists into MXs ####

#Data mngmnt
metaG.cpds_list <- list() 

for (i in seq_along(metaG_data)) {
  simulation <- metaG_data[[i]]
  metaG.cpds_list[[i]] <- simulation$compound
  cpds_metaG <- unique(unlist(metaG.cpds_list))
  names(metaG.cpds_list)[i] <- names(metaG_data)[i]
}

MAGs.cpds_list <- list() 

for (i in seq_along(MAGs_data)) {
  simulation <- MAGs_data[[i]]
  MAGs.cpds_list[[i]] <- simulation$compound
  cpds_MAGs <- unique(unlist(MAGs.cpds_list))
  names(MAGs.cpds_list)[i] <- names(MAGs_data)[i]
}

all.cpds <- unique(c(cpds_metaG, cpds_MAGs))

#  Construction of summary MXs
MX_metaG <- data.frame(compound = all.cpds)

for (i in seq_along(metaG_data)) {
  simulation <- metaG_data[[i]]
  simulation_name <- names(metaG_data)[i]
  presence_metaG <- as.integer(all.cpds %in% simulation$compound)
  MX_metaG[[simulation_name]] <- presence_metaG
}

colnames(MX_metaG)[-1:ncol(MX_metaG)] <- paste(colnames(MX_metaG[-1:ncol(MX_metaG)]), 
                                               "_metagenomes", sep = "")

MX_MAGs <- data.frame(compound = all.cpds)

for (i in seq_along(scope_lists)) {
  simulation <- scope_lists[[i]]
  simulation_name <- names(scope_lists)[i]
  presence_MAGs <- as.integer(all.cpds %in% simulation$compound)
  MX_MAGs[[simulation_name]] <- presence_MAGs
}

colnames(MX_MAGs)[-1:ncol(MX_MAGs)] <- paste(colnames(MX_MAGs[-1:ncol(MX_MAGs)]), 
                                             "_MAGs", sep = "")

MX_combined <- cbind(MX_metaG, MX_MAGs[-1:ncol(MX_MAGs)])
colnames(MX_combined[1]) <- "compound"

#Download binary matrices
write.csv2(MX_metaG.run9, file = "MX_metaG_run9.csv", row.names = FALSE)
write.csv2(MX_MAGs, file = "MX_MAGs_run9.csv", row.names = FALSE)
write.csv2(MX_combined, file = "MX_allData_run9.csv", row.names = FALSE)

#       2.4 Draw boxplots (Figs. 3A-C) ####

#Data mngmnt
col_elements <- strsplit(colnames(MX_combined)[-1], "_")
col_metadata <- do.call(rbind, col_elements)

scopes_stacked <- data.frame(
  Sites = col_metadata[, 1],       
  Dataset = col_metadata[, 3],   
  Seeds = col_metadata[, 2],  
  Scopes = colSums(MX_combined[, -1]))

scopes_stacked$Sites <- as.factor(scopes_stacked$Sites)
scopes_stacked$Dataset <- as.factor(scopes_stacked$Dataset)
scopes_stacked$Seeds <- as.factor(scopes_stacked$Seeds)
scopes_stacked$Scopes <- as.numeric(scopes_stacked$Scopes)


# /!\ Run only if not in environment already
#sites_6 <- unique(scopes_df$Site)
#site_colorBlind_6 <- c("#FFAD65", "#FF6E00", "#00FCF6", "#00BCB8", "#006666", "#4B0082")
#colors_6 <- setNames(site_colorBlind_6, sites_6)

seeds_5 <- unique(scopes_stacked$Seeds)
shapes_seeds <- c(21,22,23,24,25)
shapes_5 <- setNames(shapes_seeds, seeds_5)

#Visualize boxplots by dataset (Fig. 3A))=()
ggplot(scopes_stacked, aes(x=Dataset, y=Scopes, alpha = 0.4)) +
  geom_point(aes(colour = Sites, shape = Seeds), 
             size = 5, alpha = 0.7, stroke = 0.9) + 
  scale_color_manual(values = colors_6) +
  #scale_fill_manual(values = colors_6) + 
  scale_shape_manual(values = shapes_5) +
  scale_y_continuous(limits = c(500, 1500)) +
  geom_boxplot(aes(group = Dataset), outliers = T, 
               outlier.size = 2, outlier.shape = 20, 
               outlier.color = "black", outlier.fill = "white",
               alpha = 0.4) +
  labs(x = "Dataset", y = "Scope sizes (number of producible metabolites predicted)",
       title = "") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(margin = margin(t = 10)), #check this line
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

#Visualize boxplots by condition for MetaG-GEMs (Fig. 3B)
metagenomes <- scopes_stacked[scopes_stacked$Dataset == "Metagenomes", ]

ggplot(metagenomes, aes(x=Seeds, y=Scopes)) +
  geom_point(aes(colour = Sites, fill = Sites, shape = Seeds), size = 5, alpha=0.7) +
  scale_color_manual(values = colors_6) +
  scale_fill_manual(values = colors_6) +
  scale_shape_manual(values = shapes_5) +
  scale_y_continuous(limits = c(1000, 1500), breaks = seq(1000, 1500, by = 250)) +
  geom_boxplot(aes(group = Seeds), outliers = T, 
               outlier.size = 2, outlier.shape = 20, 
               outlier.color = "black", outlier.fill = "white",
               alpha = 0.5) +
  labs(x = "Condition", y = "Scope sizes",
       title = "MetaG-GEMs") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        #legend.text = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) 

#Visualize boxplots by condition for MAG-GEMs (Fig. 3C)
MAGs <- scopes_stacked[scopes_stacked$Dataset == "MAGs", ]

ggplot(MAGs, aes(x=Seeds, y=Scopes)) +
  geom_point(aes(colour = Sites, fill = Sites, shape = Seeds), size = 5, alpha=0.7) +
  scale_color_manual(values = colors_6) +
  scale_fill_manual(values = colors_6) +
  scale_shape_manual(values = shapes_5) +
  scale_y_continuous(limits = c(500, 1000), breaks = seq(500, 1000, by = 250)) +
  geom_boxplot(aes(group = Seeds), outliers = T, 
               outlier.size = 2, outlier.shape = 20, 
               outlier.color = "black", outlier.fill = "white",
               alpha = 0.5) +
  labs(x = "Condition", y = "Scope sizes",
       title = "MAG-GEMs") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        #legend.text = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) 

#       2.5 Draw flowcharts (Fig. S4) ####

library("diagram")

names <- paste(c("basal medium", "simple sugars", "complex sugars",
                 "non-sulfured \n amino acids", "all \n amino acids", #Check that the order follows data structure
                 "", "", "", "", "",
                 "", "", "", "", "",
                 "", "", "", "", ""))

M <- matrix(nrow = 20, ncol = 20, byrow = TRUE, data = 0)
M[2,1] <- M[3,1] <- M[4,1] <- M[5,4] <- "flow"

core <- 469 #core scope across simulations (928 in MetaG ∩ 472 in MAGs)

curves <- matrix(nrow = 20, ncol = 20, byrow = TRUE, data = 0)
curves[2,1] <- 0
curves[3,1] <- 0.3
curves[4,1] <- -0.3
curves[5,4] <- 0.3

flechas <- matrix(nrow = 20, ncol = 20, byrow = TRUE, data = 0)
flechas[2,1] <- 0.6
flechas[3,1] <- 0.7
flechas[4,1] <- 0.7
flechas[5,4] <- 0.65

positions <- cbind(c(0.4,0.4,0.1,0.7,0.9,
                    0.4,0.4,0.1,0.7,0.9,
                    0.4,0.4,0.1,0.7,0.9,
                    0.4,0.4,0.1,0.7,0.9),
                  c(0.7,0.38,0.4,0.4,0.15,
                    0.7,0.38,0.4,0.4,0.15,
                    0.7,0.38,0.4,0.4,0.15,
                    0.7,0.38,0.4,0.4,0.15))

#Visualize flowcharts for each site

for (site in unique(scopes_stacked$Sites)) {
  site_data <- subset(scopes_stacked, Sites == site)
  metaG_data <- site_data[site_data$Dataset == "Metagenomes", ]
  MAGs_data <- site_data[site_data$Dataset == "MAGs", ]
  box_size <- c((metaG_data$Scopes - rep(core, 5)) / 9000,          # Scopes for Metagenomes by condition
                rep(metaG_data$Scopes[1] - core, 5) / 9000,         # Basal medium Metagenomes
                (MAGs_data$Scopes - rep(core, 5)) / 9000,           # Scopes for MAGs by condition
                rep(MAGs_data$Scopes[1] - core, 5) / 9000)          # Basal medium MAGs
  plotmat(M,                        
          pos = positions, 
          curve = curves, 
          name = "", 
          shadow.size = 1, shadow.col = "grey",
          cex.txt = 0, txt.col = "white", 
          txt.xadj = 0, txt.yadj = 0,
          lwd = 2, box.lwd = 0.1, 
          box.type = "circle", box.cex = 0, 
          box.size = box_size, 
          box.col = c("#023047","#023047","#023047","#023047","#023047",
                      "#64b5f6","#64b5f6","#64b5f6","#64b5f6","#64b5f6",
                      "#023047","#023047","#023047","#023047","#023047",
                      "#d9ed92","#d9ed92","#d9ed92","#d9ed92","#d9ed92"),
          box.lcol = "white",
          arr.col = "black", arr.lcol = "black", 
          arr.type = "triangle", 
          arr.length = 0.3, arr.width = 0.3,
          arr.pos = flechas, 
          endhead = TRUE, 
          cex.main = 1,
          main = paste("Flowchart for", site))
}

#       2.6 Draw PCoA (Figs. 3D,E) ####

#Calculate the (dis)similarity index for each dataset
dist_metagenomes <- vegdist(t(MX_metaG[,-1:ncol(MX_metaG)]), method = "jaccard")
dist_MAGs <- vegdist(t(MX_MAGs[,-1:ncol(MX_MAGs)]), method = "jaccard")

#Compute PCoAs
pcoa_metagenomes <- cmdscale(dist_metagenomes, eig = T, add = T)
pcoa_MAGs <- cmdscale(dist_MAGs, eig = T, add = T)

#Data mngmnt
colnames(pcoa_metagenomes$points) <- c("PCo1", "PCo2")
colnames(pcoa_MAGs$points) <- c("PCo1", "PCo2")

#Extract coordenates of first two dimensions
pcoa_coords_metagenomes <- as.data.frame(pcoa_metagenomes$points)
pcoa_coords_MAGs <- as.data.frame(pcoa_MAGs$points)

#Add source info for legends: colors
pcoa_coords_metagenomes$sites <- sort(rep(sites_6,5))
pcoa_coords_MAGs$sites <- sort(rep(sites_6,5))

#Add source info for legends: shapes
pcoa_coords_metagenomes$seeds <- rep(seeds_5,6)
pcoa_coords_MAGs$seeds <- rep(seeds_5_2,6)

#Extract explained var
pcoa_eig_metagenomes <- pcoa_metagenomes$eig
pcoa_eig_MAGs <- pcoa_MAGs$eig

#Cut groups for unsupervised clustering
k_metagenomes <- kmeans(pcoa_coords_metagenomes[, c("PCo1", "PCo2")], centers = 7)
k_MAGs <- kmeans(pcoa_coords_MAGs[, c("PCo1", "PCo2")], centers = 7)

#Assign unsupervised clusters
pcoa_coords_metagenomes$cluster <- as.factor(k_metagenomes$cluster)
pcoa_coords_MAGs$cluster <- as.factor(k_MAGs$cluster)

# Visualize PCoA of MetaG-GEMs (Fig. 3D)
ggplot(pcoa_coords_metagenomes, mapping = aes(x=PCo1, y=PCo2, color = sites, shape = seeds)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray") + 
  stat_ellipse(aes(x = PCo1, y = PCo2, group = cluster), type = "t", level = 0.95, linetype = "dashed", color = "blue") + 
  geom_point(aes(colour = sites, shape = seeds, alpha = 0.7), size = 5, stroke = 0.9) +  
  labs(title = "MetaG-GEMs",
       x = paste0("PCo1 (",  round(100 * pcoa_eig_metagenomes / sum(pcoa_eig_metagenomes),1)[1], "%)"),
       y = paste0("PCo2 (",  round(100 * pcoa_eig_metagenomes / sum(pcoa_eig_metagenomes),1)[2], "%)")) +
  scale_color_manual(values = colors_6) +
  scale_shape_manual(values = shapes_5) +
  theme_minimal()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        axis.text.x = element_text(margin = margin(t = 10)), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5))

# Visualize PCoA of MAG-GEMs(Fig. 3E)
ggplot(pcoa_coords_MAGs, mapping = aes(x=PCo1, y=PCo2, color = sites, shape = seeds)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
  stat_ellipse(aes(x = PCo1, y = PCo2, group = cluster), type = "t", level = 0.95, linetype = "dashed", color = "blue") +
  geom_point(aes(colour = sites, shape = seeds, alpha = 0.7), size = 5, stroke = 0.9) + 
  labs(title = "MAG-GEMs",
       x = paste0("PCo1 (",  round(100 * pcoa_eig_MAGs / sum(pcoa_eig_MAGs),1)[1], "%)"),
       y = paste0("PCo2 (",  round(100 * pcoa_eig_MAGs / sum(pcoa_eig_MAGs),1)[2], "%)")) +
  scale_color_manual(values = colors_6) +
  scale_shape_manual(values = shapes_5) +
  theme_minimal()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
        axis.text.x = element_text(margin = margin(t = 10)), 
        axis.text.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10), size= 11),
        axis.title.y = element_text(margin = margin(r = 10), size= 11),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5))

#       2.7 Draw heatmap (Fig. S5) ####

#Create column with binary vector of occurrence profiles across simulations
fx_create_pattern <- function(row) {
  paste(row, collapse = "")
}

MX_combined$Pattern_ID <- apply(data[, -1], 1, fx_create_pattern)

#Assign dummy name (rename later as you wish) and store in empty df
pattern_mapping <- data.frame(
  Pattern_ID = unique(MX_combined$Pattern_ID),
  Group_name = paste("X", seq_along(unique(MX_combined$Pattern_ID)), sep = "_"))

#Fill df with dereplicated Metabolite groups (n=269) and download
MX_by_metGroup <- merge(MX_combined, pattern_mapping, by = "Group_name")
write.csv2(MX_by_metGroup, "MX_metGroups_269x60simulations.csv")

#Modify if in environment already
site_colorBlind <-  c("#FFAD65", "#FF6E00", "#00FCF6", "#00BCB8", "#006666", "#4B0082","white")
sites_6 <- c("S1", "S2", "S3", "S4", "S5", "S6", "0") #We add empty value and "white" for absence
colors_7 <- structure(site_colorBlind, names= sites)

#Data mngmnt
rownames(MX_by_metGroup) <- MX_by_metGroup$Group_name
MX_by_metGroup <- MX_by_metGroup[-1:ncol(MX_by_metGroup)]

library(dplyr)

#Create copy of binary matrix
MX_num <- MX_by_metGroup

#Replace 1s for site of source
fx_replace_1s <- function(df, column) {
  prefix <- gsub("^(S[0-9]+).*", "\\1", names(df)[column])
  df[[column]][df[[column]] == 1 ] <- prefix
  return(df)
}

for (col in 1:ncol(MX_by_metGroup)) {
  MX_by_metGroup <- fx_replace_1s(MX_by_metGroup, col)
}

library(ComplexHeatmap)

#Draw heatmap 
Heatmap(as.matrix(t(MX_by_metGroup)), 
        col = colors_7, 
        name = "Sites",
        clustering_distance_columns = function(x) {
        vegdist(MX_num, method = "jaccard")},
        clustering_method_columns = "ward.D2",
        cluster_columns = TRUE,
        show_column_names = F,
        show_row_names = F,
        row_dend_side = "left", 
        #row_names_side = "left", 
        column_dend_side = "bottom", 
        #column_names_side = "top", 
        column_dend_height = unit(15, "mm"),
        row_dend_width = unit(10, "mm"),
        column_split = 3,
        column_gap = unit(3, "mm"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        #top_annotation = HeatmapAnnotation(ontology = onto_combi,
        #                                   col = list(ontology = category_colors)))

#Explore tree
ht_drawn <- draw(ht_jac)
column_order(ht_drawn)[1] #16 NON.SULFURED
column_order(ht_drawn)[2] #42 SULFURED
column_order(ht_drawn)[3] #74 METAG DISTINTIONS     
column_order(ht_drawn)[4] #42 METAG GENE RESERVOIR  
column_order(ht_drawn)[5] #95 METAG GENE RESERVOIR  

#Download lists of compounds associated to all and non-sulfured amino acids
write.table(column_order(ht_drawn)[1], "list_NITROGENATED_patterns_16.txt", sep = ";")
write.table(column_order(ht_drawn)[2], "list_SULFURED_patterns_42.txt", sep = ";")

#### 3) Selection of key metabolites with an elastic net regression ####

#       3.1 Creation of the elastic net ####

library(dplyr)

#If not in environment already, load data (see subsection 1.1)
#metadata <- read.table("MX_environm_metadata.csv", header = TRUE, dec = ".",sep = ";")
#rownames(metadata) <- metadata[,1]
#metadata_real_S1.S6 <- as.data.frame(scale(metadata[,-1]))
#MX_by_metGroup <- read.table("MX_metGroups_269x60simulations.csv", header = TRUE, dec = ".",sep = ";")
#rownames(MX_by_metGroup) <- MX_by_metGroup$Group_name
#MX_by_metGroup <- MX_by_metGroup[-1]

#Add site to environmental metadata
metadata_real_S1.S6$site <- rownames(metadata_real_S1.S6)

#Transpose metabolic data and add site of source
t_MX_by_metGroup <- as.data.frame(t(MX_by_metGroup))
t_MX_by_metGroup$site <- gsub("^(S[1-6]).*", "\\1", rownames(t_MX_by_metGroup))
t_MX_by_metGroup$site <- as.factor(MX_pttrns_temp$site)

#Optional: Add seed to metabolic data
#t_MX_combined$seed <- gsub(".*_(.*)_.*", "\\1", rownames(MX_pttrns_temp))
#t_MX_combined$seed <- as.factor(MX_pttrns_temp$seed)

#Concatenate environmental and metabolic data
MX_elasticNet <- left_join(t_MX_by_metGroup, metadata_real_S1.S6, by = "site")
rownames(MX_elasticNet) <- rownames(t_MX_by_metGroup) #simulations  

#Extract input data for feature selection with Elastic Net (alpha=0.85)
data_269 <- MX_elasticNet[,1:269]  #predictor variables
data_17 <- MX_elasticNet[,271:287] #response variables
response_vars <- colnames(data_17)
data_list <- lapply(response_vars, function(var) {
  data_17[[var]]
})

names(data_list) <- response_vars

library(glmnet)

#Create function for elastic net
fx_Lasso <- function(x, y, environm_var, alpha_value) {

  #Cross validation
  cv_model_mse <- cv.glmnet(as.matrix(x), y, alpha = alpha_value, type.measure = "mse", 
  family = "gaussian", intercept = FALSE)
  cv_model_mse_0.5 <- cv.glmnet(as.matrix(x), y, alpha = 0.5, type.measure = "mse", 
  family = "gaussian", intercept = FALSE)
  cv_model_mse_0.75 <- cv.glmnet(as.matrix(x), y, alpha = 0.75, type.measure = "mse", 
  family = "gaussian", intercept = FALSE)
  cv_model_mse_1 <- cv.glmnet(as.matrix(x), y, alpha = 1, type.measure = "mse", 
  family = "gaussian", intercept = FALSE)
  
  #plot(cv_model_mse)
  plot(log(cv_model_mse$lambda), cv_model_mse$cvm , pch = 19, col = "blue",
       xlab = "log(Lambda)", ylab = cv_model_mse$name)
  points(log(cv_model_mse_1$lambda), cv_model_mse_1$cvm , pch = 19, col = "red")
  points(log(cv_model_mse_0.75$lambda), cv_model_mse_0.75$cvm , pch = 19, col = "darkgrey")
  points(log(cv_model_mse_0.5$lambda), cv_model_mse_0.5$cvm, pch = 19, col = "lightgrey")
  
  legend("topleft", legend = c("alpha= 1.0", paste0("alpha= ", alpha_value), "alpha= 0.75", "alpha= 0.5"),
         pch = 19, col = c("red", "blue","darkgrey","lightgrey"))
  title(main = paste(environm_var))
 
  #/!\ you gotta choose between min and 1se 
  model_lambda_mse <- glmnet(as.matrix(x), y, alpha = alpha_value, lambda = cv_model_mse$lambda.min,
                             intercept = FALSE) #asumes mse

  #Extract coefficients                        
  coefs_mse <- as.matrix(coef(model_lambda_mse, model_lambda_mse$lambda))
  nonzero_coefs_mse <- coefs_mse[coefs_mse!=0,]
  print(paste("Elastic net model (by minimum square error) for", environm_var, ": average absolute lambda of", 
              length(nonzero_coefs_mse), "non-zero values =", abs(unique(ave(nonzero_coefs_mse)))))

  #Data mngmnt and download
  df_coefs <- data.frame(met_Group = names(nonzero_coefs_mse), 
                         coefficients = nonzero_coefs_mse,
                         absolute_values = abs(nonzero_coefs_mse))
  df_coefs_ord <- df_coefs[order(df_coefs$absolute_values, decreasing = TRUE),]
  write.table(df_coefs_ord, paste0("ElasticNet_", environm_var, ".tsv"), sep = "\t", dec = ".", row.names = FALSE)
  return(df_coefs_ord)
}

#Compute elastic net
alpha_value = 0.85 #Customize this value according to your needs
#/!\ Remember that alpha=1 is a lasso regression and alpha=0 is a ridge rergession

lassos <- lapply(names(data_list), function(var) {
  fx_Lasso(data_269, data_list[[var]], var, alpha_value)
})

names(lassos) <- names(data_list)

#Extract nonzero env variables
nonzero_coefs <- data.frame(met_Group = character(0),
                         coefficients = character(0),
                         absolute_values = character(0))

for (i in seq_along(lassos)) {
  met_Group <- lassos[[i]]
  met_Group$environm_var <- names(lassos)[i]
  nonzero_coefs <- rbind(nonzero_coefs, met_Group)
}

#Download data
write.table(nonzero_coefs, "nonzero_metGroups_elasticNet_alpha0.85.csv", sep = ";", dec = ".", row.names = FALSE)

#       3.2 Summary (Fig. 4A) ####

#Build dataframe for heatmap
nonzero_data <- unstack(nonzero_coefs, coefficients ~ environm_var)
names_data <- unstack(nonzero_coefs, met_Group ~ environm_var)
metGroups_coeffs <- list()

for (i in seq_along(nonzero_data)) {
  coeff <- nonzero_data[[i]]
  metGroup <- names_data[[i]]
  names(coeff) <- metGroup
  metGroups_coeffs[[i]] <- coeff
  names(metGroups_coeffs)[i] <- names(nonzero_data)[i] 
}

nonzero_metGroups <- unique(nonzero_coefs$met_Group)
nonzero_df <- data.frame(met_Group = nonzero_metGroups)

for (i in seq_along(metGroups_coeffs)) {
  coeff <- metGroups_coeffs[[i]]
  environm_var <- names(nonzero_data)[i]
  coeff_vector <- ifelse(nonzero_metGroups %in% names(coeff), coeff[match(nonzero_metGroups, names(coeff))], 0)
  nonzero_df[[environm_var]] <- coeff_vector
}
        
#Filter nonzero coefficients to define key features
nonzero_df[nonzero_df  >= -0.3 & nonzero_df <= 0.3] <- 0
key_metGroups <- nonzero_df[!apply(nonzero_df[-1] == 0, 1, all), ] 
write.table(nonzero_coefs, "key_metGroups_elasticNet_alpha0.85_forHeatmap.csv", sep = ";", dec = ".", row.names = FALSE)

palette <- colorRampPalette(c("blue", "white"))(18) 
palette <- c(palette, colorRampPalette(c("white", "red"))(15))
palette <- unique(palette)

library(ComplexHeatmap)

Heatmap(as.matrix(key_metGroups), 
        col = palette, 
        name = "coef.",
        clustering_distance_rows = function(x) {
          vegdist(abs(key_metGroups), method = "bray")},
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = function(x) {
        vegdist(t(abs(key_metGroups)), method = "bray")},
        clustering_method_columns = "ward.D2",
        cluster_columns = T,
        cluster_rows = T,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_names_rot=TRUE,
        #row_names_side = "left",
        #row_dend_side = "right", 
        row_dend_width = unit(20, "mm"),
        column_names_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 9),
        column_names_centered = TRUE)
