library(tidyverse)
library(CorLevelPlot)
library(WGCNA)
library(gridExtra)
library(CorLevelPlot) 
library(Hmisc)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)
library(RColorBrewer)
library(phyloseq)
#setting
print_flag <- TRUE
pdf_width <- 8
pdf_height <- 6
Depth_names <- c("0_5cm","5_10cm","10_20cm","20_40cm",
                 "40_60cm","60_80cm","80_100cm","100_120cm",
                 "120_140cm","140_160cm","160_180cm","180_200cm")
Site_names <- c("AZ","SG","LC","NB")
# load data --------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1,
                  header = TRUE,
                  na.strings = "NA",
                  check.names = FALSE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                header = TRUE,
                check.names = FALSE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  row.names = 1,
                  header = TRUE)
# clean data
env.idna <- env %>%
  filter(dna_typen == 0, 
         sample_id %in% colnames(asv)) %>%
  data.frame()

asv.idna <- asv[, colnames(asv) %in% env.idna$sample_id] %>%
  filter(rowSums(.) > 0)

tax.idna <- tax[rownames(tax) %in% rownames(asv.idna), ]
# 
# tax.idna_filtered <- tax.idna %>% 
#   filter(Phylum != "p_NA")
# asv.idna_filtered <- asv.idna[rownames(asv.idna) %in% rownames(tax.idna_filtered), ]
# 
# tax.idna <- tax.idna_filtered
# asv.idna <- asv.idna_filtered
# Averaging
asv_tmp <- as.data.frame(t(asv.idna))
env_tmp <- env.idna[match(rownames(asv_tmp), env.idna$sample_id), ]
asv_tmp$site_dna_depth <- env_tmp$site_dna_depth
#average asv
# asv.idna.av <- asv_tmp %>% 
#   group_by(site_dna_depth) %>% 
#   summarise(across(everything(), mean, na.rm = TRUE)) %>% #MD什么鬼，原来能运行，更新后就不行啦。shit
#   column_to_rownames(var = "site_dna_depth")

asv.idna.av <- asv_tmp %>% 
  group_by(site_dna_depth) %>% 
  summarise_all(list(mean = ~mean(., na.rm = TRUE))) %>% 
  column_to_rownames(var = "site_dna_depth")


env.idna.av <- env_tmp %>%
  group_by(site_dna_depth) %>%
  #summarise(across(c(depth, pH:N, C, CN:Sio, NH4:Mnp), mean, na.rm = TRUE)) %>%
  summarise_at(vars(c(depth, pH:N, C, CN:Sio, NH4:Mnp)), mean, na.rm = TRUE) %>%
  column_to_rownames(var = "site_dna_depth")


# get site and depth order
Site <- distinct(env_tmp, 
                 site_dna_depth, 
                 site_num) %>% 
  column_to_rownames(var = "site_dna_depth")

Depth <- distinct(env_tmp, 
                  site_dna_depth, 
                  depth_order) %>% 
  column_to_rownames(var = "site_dna_depth")

#get relative abundance
asv.idna.av.rel <- 100 * asv.idna.av / rowSums(asv.idna.av)
#script from website
#https://github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R
#https://ramellose.github.io/networktutorials/wgcna.html
# WGCNA data prepare -------------------------------------------------------------------
asv <- asv.idna.av.rel
env <- env.idna.av
tax <- tax.idna
asv <- asv[colMeans(asv) > 0.02]# Remove ASVs mean abundance < 0.02%
asv[asv < 0.05] <- 0# remove low abundance ASV (less than 0.05%)
asv <- asv / rowSums(asv)# rescale the relative abundance
#rename rownames (remove "iDNA") and save relative abundance to CSV file
rownames(asv) <- sub("iDNA", "_", rownames(asv)) #replace iDNA with _
write.csv(asv, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/relative_abundance_filtered.csv")
rownames(env) <- sub("iDNA", "_", rownames(env))
rownames(Site) <- sub("iDNA", "_", rownames(Site))
rownames(Depth) <- sub("iDNA", "_", rownames(Depth))
# WGCNA analysis ----------------------------------------------------------
# 1) Quality check ----
gsg = goodSamplesGenes(asv, verbose = 3)
gsg$allOK # Note: if it is false, then activate the exclude outlayers step below
# Show the clustering quality - hierarchical clustering with traits
sampleTree = hclust(dist(asv), 
                    method = "average")
traitColors = numbers2colors(env, 
                             signed = FALSE)
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/simple_dendrogram_trait_heatmap.pdf",
    width = pdf_width, 
    height = pdf_height)

plotDendroAndColors(sampleTree,
                    traitColors, 
                    groupLabels = names(env),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# Show the clustering quality - PCA method 
pca <- prcomp(asv)
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, 
                         digits = 2)
pca.dat <- as.data.frame(pca.dat)
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/PCA_ASV.pdf",
    width = pdf_width, 
    height = pdf_height)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()
# Exclude outlayer ASVs(Currently there is no outliers)
if (!gsg$allOK) {
  asv <- asv[,gsg$goodGenes == TRUE]
}
# 2) Normalization --------------------------------------------------------
asv <- sqrt(asv)
# 3) Network Construction -------------------------------------------------
# Choose a set of soft-threshold powers
power <- c(1:10, 
           seq(12, 50, 
               by = 2))
# Call the network topology analysis function
#sft = scale-free topology criteria. 
sft <- pickSoftThreshold(asv,
                         powerVector = power,
                         networkType = "unsigned",
                         verbose = 5)
sft.data <- sft$fitIndices
# visualization to pick power
p1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

p2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/R2_connectivity.pdf",
    width = pdf_width, 
    height = pdf_height)
grid.arrange(p1, p2, nrow = 2)
dev.off()
# Network construction and module detection
soft_power <- 7  # picked from the graph above (high R² and low connectivity)
# co-occurrence network base on WGCNA (changed by Xiuling)-------------------------------------
temp_cor <- cor
cor <- WGCNA::cor
bwnet <- blockwiseModules(asv,
                          maxBlockSize = 14000,
                          TOMType = "unsigned",
                          networkType="unsigned",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# 4) Module Eigengenes ----------------------------------------------------------
module_eigengenes <- bwnet$MEs#特征基因（18个模块）
head(module_eigengenes) # preview
table(bwnet$colors)     # get number of genes for each module
#grey module = others gene
# Plot the dendrogram and the module colors before and after merging underneath
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/dendrogram_modules.pdf",
    width = pdf_width,
    height = pdf_height * 2)#double height to add the legend

plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, 
                          bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
# plot legend
moduleColors <- unique(bwnet$colors)
par(mar = c(1, 1, 1, 1)) # 设置较小的边距
plot(1, 
     type = "n", 
     xlim = c(1, length(moduleColors)), 
     ylim=c(0, 2), 
     xlab="", 
     ylab="", 
     xaxt="n", 
     yaxt="n", 
     bty="n")  # 创建一个空图用于绘制图例
legend("topright",
       legend = moduleColors,  
       fill = moduleColors,
       cex = 0.8, 
       bty = "n", 
       x.intersp = 0.5,  # 调整图例内元素间的水平间距
       y.intersp = 1)  # 调整图例内元素间的垂直间距
dev.off()

# module colors
module_colors <- bwnet$colors
table(module_colors)# 18 modules
# prepare adjacency_matrix for network
adjacency_matrix <- adjacency(as.matrix(asv), 
                              power = soft_power,
                              type="unsigned")
heatmap(adjacency_matrix, 
        labRow = FALSE, 
        labCol = FALSE)

TOM <- TOMsimilarity(adjacency_matrix, # TOM：Topological Overlap Matrix
                     TOMType = "unsigned")
heatmap(as.matrix(TOM), 
        labRow=FALSE, 
        labCol=FALSE)

net <- graph_from_adjacency_matrix(TOM, 
                                   mode = "undirected", 
                                   weighted = TRUE)


tax <- tax.idna %>% 
  filter(rownames(.) %in% colnames(asv)) %>% 
  as.data.frame()

V(net)$Phylum <- tax$Phylum
V(net)$name <- rownames(tax)

net <- delete_edges(net, E(net)[E(net)$weight < 0.1])
V(net)$color <- module_colors # module_colors <- bwnet$colors
net <- simplify(net)  # removes self-loops
net <- delete.vertices(net,degree(net) == 0)

# Network parameters:
net_nodes <- vcount(net)
print(net_nodes)

net_edges <- ecount(net)
print(net_edges)

net_diameter <- diameter(net, 
                         directed = FALSE)
print(net_diameter)

clustering_coefficient <- transitivity(net, 
                                       type = "local")
print(clustering_coefficient)

average_clustering <- mean(clustering_coefficient, 
                           na.rm = TRUE)
print(average_clustering)

avg_path_length <- average.path.length(net)
print(avg_path_length)

graph_density <- graph.density(net)
print(graph_density)

node_degrees <- degree(net)
print(node_degrees)

avg_degree <- mean(degree(net))
print(avg_degree)

closeness_centrality <- closeness(net)
print(closeness_centrality)

betweenness_centrality <- betweenness(net)
print(betweenness_centrality)


network_features <- data.frame(
  Nodes = net_nodes,
  Edges = net_edges,
  Diameter = net_diameter,
  Average_Clustering = average_clustering,
  Avg_Path_Length = avg_path_length,
  Graph_Density = graph_density,
  Average_Degree = avg_degree
  # 添加其他网络特征
)


write.csv(network_features, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/Figure for papers😊/table/TableS4.csv")

#查看顶点属性
vertex_attributes <- vertex_attr_names(net)
vertex_attributes

# Node topological features ----
node.topology <- data.frame(
  name = V(net)$name,
  Phylum = V(net)$Phylum,
  node.degree = degree(net),
  closeness_centrality = closeness(net),
  betweenness_centrality = betweenness(net)
)
head(node.topology)
unique_phyla <- unique(node.topology$Phylum)
num_colors <- length(unique_phyla)
num_colors
colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_colors)
# Create a color dictionary with class-color mappings
color_dict <- setNames(colors, unique_phyla)

write.csv(node.topology, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/Figure for papers😊/table/TableS5.csv")
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/Keystone_Species.pdf",
    width = pdf_width,
    height = pdf_height)

ggplot(node.topology, 
       aes(x = node.degree, 
           y = betweenness_centrality)) +
  geom_point(aes(color = as.factor(Phylum)), size = 3) +
  scale_color_manual(values = colors) +
  labs(x = "Degree", y = "Betweenness") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(linetype = "dotted", color = "gray"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) 

ggplot(node.topology, aes(x = node.degree, 
                          y = betweenness.centrality)) +
  geom_point(aes(color = Phylum), size = 3) +
  scale_color_manual(values = colors) +
  geom_text(aes(label = rownames(node.topology)), 
            size = 2, 
            hjust = 0, 
            vjust = 0) +
  labs(x = "Degree", 
       y = "Betweenness") +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(linetype = "dotted", color = "gray"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

ggplot(node.topology, 
       aes(x = node.degree, 
           y = closeness_centrality)) +
  geom_point(aes(color = as.factor(Phylum)), size = 3) +
  scale_color_manual(values = colors) +
  labs(x = "Degree", y = "Closeness") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(linetype = "dotted", color = "gray"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) 

dev.off()
# test plot
par(mar=c(0,0,0,0))
plot(net, 
     layout=layout.fruchterman.reingold(net), 
     edge.arrow.size = 0.2)
# plot with ggraph
color_palette <- c(
  "yellow" = "yellow",
  "tan" = "tan",
  "pink" = "pink",
  "greenyellow" = "greenyellow",
  "cyan" = "cyan",
  "salmon" = "salmon",
  "midnightblue" = "midnightblue",
  "magenta" = "magenta",
  "blue" = "blue",
  "purple" = "purple",
  "lightcyan" = "lightcyan",
  "green" = "green",
  "brown" = "brown",
  "turquoise" = "turquoise",
  "red" = "red",
  "grey60" = "grey60",
  "black" = "black",
  "grey" = "grey"
)
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/co-ocurrence network.pdf",
    width = pdf_width,
    height = pdf_height)

set.seed(104)#good

ggraph(net, layout = 'fr') + 
  geom_edge_link(aes(width = E(net)$weight),  # 使用边的权重作为宽度的基础
                 color = 'lightgray', 
                 alpha = 0.9) +  
  geom_node_point(aes(color = V(net)$color), 
                  size = 3) +
  scale_edge_width(range = c(0.1, 2)) +  # 定义边宽度的范围
  scale_color_manual(values = color_palette) +
  theme_void() +
  ggtitle("Co-occurrence Network") +
  guides(color = guide_legend(override.aes = list(size = 4)))  # 调整图例色块大小

dev.off()
#
site_annotations <- c("AZ" = "purple",
                      "AZ" = "lightcyan",
                      "AZ" = "green",
                      "AZ" = "brown",
                      
                      "SG" = "turquoise",
                      "SG" = "red",
                      "SG" = "grey60",
                      "SG" = "black",
                      
                      "LC" = "yellow",
                      "LC" = "tan",
                      "LC" = "pink",
                      "LC" = "greenyellow",
                      "LC" = "cyan",
                      
                      "NB" = "salmon",
                      "NB" = "midnightblue",
                      "NB" = "magenta",
                      "NB" = "blue",
                      
                      "GR" = "grey")
cor <- temp_cor
# 6) Relate modules to traits (environmental variables) -------------------------
# report correlate using Pearson correlation and p-value
if (print_flag) {
  module.trait.corr <- cor(module_eigengenes, env, 
                           use = 'p')
  module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nrow(asv))
  write.csv(module.trait.corr, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/Module_trait_r.csv")
  write.csv(module.trait.corr.pvals, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/Module_trait_p-value.csv")
}

## rearrange the modules by sorting them alphabetically and remove MEgry
# Get the column names that start with "ME"

heatmap.data <- merge(module_eigengenes, env, by = 'row.names') %>%
  column_to_rownames(var = 'Row.names')

me_columns <- names(heatmap.data)[grepl("^ME", names(heatmap.data))]
sorted_me_columns <- sort(me_columns)

heatmap.data <- heatmap.data %>%
  select(sort(grep("^ME", names(.))), everything()) %>%
  select(-MEgrey)
#colnames(heatmap.data)  # show names to choose the x,y names below
nMod <- length(module_eigengenes)-1
nEnd <- length(colnames(heatmap.data))
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/environment_modules_heatmap.pdf",
    width = pdf_width+10, 
    height = pdf_height)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[(nMod+1):nEnd],
             y = names(heatmap.data)[1:nMod],
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()
# 7) Relate modules to sites and depths -----------------------------------------
# Sites
# TODO: Add description
Site_factors <- Site
Site_factors$site_num <- factor(Site$site_num, 
                                levels = c(1, 2, 3, 4))

Site.data <- binarizeCategoricalColumns(Site_factors$site_num,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)
Site.data$data.1.vs.all <- ifelse(Site_factors$site_num == 1, 1, 0)
Site.data <- Site.data %>%
  select(data.1.vs.all, everything())

rownames(Site.data) <- rownames(Site)
colnames(Site.data) <- Site_names

heatmap.site.data <- merge(module_eigengenes, 
                           Site.data, 
                           by = 'row.names')
heatmap.site.data <- heatmap.site.data %>% 
  column_to_rownames(var = 'Row.names')

# sort the modules similar to the environmental factors above
other_columns <- setdiff(colnames(heatmap.site.data), 
                         me_columns)
reordered_columns <- c(sorted_me_columns, 
                       other_columns)

# Reorder the data frame based on the new column order
heatmap.site.data <- heatmap.site.data[, reordered_columns]
heatmap.site.data <- heatmap.site.data[, !colnames(heatmap.site.data) %in% "MEgrey"]

colnames(heatmap.site.data) # to show the x,y index below
nEnd <- length(colnames(heatmap.site.data))

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/sites_modules_heatmap.pdf",
    width = pdf_width-3, 
    height = pdf_height)

CorLevelPlot(heatmap.site.data,
             x = names(heatmap.site.data)[(nMod+1):nEnd],
             y = names(heatmap.site.data)[1:nMod],
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()
# here we got site_annotations
site_annotations <- c("AZ" = "purple",
                      "AZ" = "lightcyan",
                      "AZ" = "green",
                      "AZ" = "brown",
                      
                      "SG" = "turquoise",
                      "SG" = "red",
                      "SG" = "grey60",
                      "SG" = "black",
                      
                      "LC" = "yellow",
                      "LC" = "tan",
                      "LC" = "pink",
                      "LC" = "greenyellow",
                      "LC" = "cyan",
                      
                      "NB" = "salmon",
                      "NB" = "midnightblue",
                      "NB" = "magenta",
                      "NB" = "blue",
                      
                      "GR" = "grey")
# Depth
# TODO: Add description
Depth_factors <- Depth
Depth_factors$depth_order <- factor(Depth$depth_order,
                                    levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

Depth.data <- binarizeCategoricalColumns(Depth_factors$depth_order,
                                         includePairwise = FALSE,
                                         includeLevelVsAll = TRUE,
                                         minCount = 1)
Depth.data$data.1.vs.all <- ifelse(Depth_factors$depth_order == 1, 1, 0)
Depth.data <- Depth.data %>%
  select(data.1.vs.all, everything())

rownames(Depth.data) <- rownames(Depth)
colnames(Depth.data) <- Depth_names

heatmap.Depth.data <- merge(module_eigengenes, Depth.data, by = 'row.names')
heatmap.Depth.data <- heatmap.Depth.data %>% 
  column_to_rownames(var = 'Row.names')

colnames(heatmap.Depth.data) # to show the x,y index below
nEnd <- length(colnames(heatmap.Depth.data))

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_figure/depth_modules_heatmap.pdf",
    width = pdf_width+8, 
    height = pdf_height)

CorLevelPlot(heatmap.Depth.data,
             x = names(heatmap.Depth.data)[(nMod+2):nEnd],               # depth
             y = names(heatmap.Depth.data)[1:nMod],                # Modules
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()
##---- Outputs ----
# 1) table of ASV: asv, module, tax, rel. ab.
# ASV - Modules
module.asv.mapping <- as.data.frame(bwnet$colors)
write.csv(module.asv.mapping, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/ASV_module.csv")
# ASV - Modules - taxa
tax_order <- match(rownames(module.asv.mapping),rownames(tax))
tax.idna.order <- tax[tax_order,]
data <- cbind(module.asv.mapping, tax.idna.order)
# order the table
sorted_data <- data[order(data[, 1], data[, 2], data[, 3],
                          data[, 4], data[, 5], data[, 6], data[, 7]), ]
write.csv(sorted_data, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/ASV_module_taxa.csv")
# ASV - rel. ab. - Modules
asv_t <- as.data.frame(t(asv))
asv_order <- match(rownames(module.asv.mapping),rownames(asv_t))
asv.mod <- asv_t[asv_order,]
combined_df <- cbind(module.asv.mapping, tax.idna.order,asv.mod)
write.csv(combined_df, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/ASV_module_relative_abu_taxa.csv")

# 2) Modules overview
# Modules
table_output <- table(bwnet$colors)
# top ASV in each module
hub_output <- chooseTopHubInEachModule(asv, module.asv.mapping$'bwnet$colors')
# Create a character vector combining the outputs
output_lines <- c("Modules:",
                  paste(rownames(table_output), 
                        table_output, sep = ": "),
                  "",
                  "Top ASV in each module:",
                  paste(rownames(as.data.frame(hub_output)),hub_output, sep = ": "))
# Save the combined output to a text file
writeLines(output_lines, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/Modules_overview.txt")

# 3) Top hub in each module with taxa
ASV_rows <- match(hub_output,rownames(tax))
TopHub_ASV <- tax[ASV_rows, ]
TopHub_ASV$Module <- rownames(as.data.frame(hub_output))
write.csv(TopHub_ASV, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/scripts/01 🌸Chile Bacteria/Rahma_network/mid_data/ASV_top_hub.csv")
# Example of how to show a certain module
module.asv.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()