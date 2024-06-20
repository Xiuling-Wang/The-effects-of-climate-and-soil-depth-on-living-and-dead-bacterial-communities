library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rdacca.hp)
library(RColorBrewer)
# load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names =1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                row.names =1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  row.names = 1,
                  header = TRUE)
# 01. top ASV for all sample-------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv),
         depth < 200) %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]

filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,] 
env2 <- env1 %>% 
  mutate(group_col = paste(site_dna,layer_label,sep = "_"))
  
group_info <- data.frame(col_name = rownames(env2),
                         group = env2$group_col)
asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", value = "value", -id) %>%  # Convert to long format
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env2$group_col)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]

total_reads <- colSums(asv2.8)
relative_abundance <- 100*asv2.8 / total_reads

rankRA <- relative_abundance %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance))
  
top_50_asv <- rankRA %>%
  slice_head(n = 50) %>% 
  select(-mean_abundance)

filt5 <- rownames(tax) %in% rownames(top_50_asv)
tax_1 <- tax[filt5,]

# 为 asv-table 和 tax-table 添加行名
asv_table <- top_50_asv %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax_1 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)
# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)
unique(long_data$sample_id)
tmp123 <- c("AZ__iDNA_0_40cm","AZ__iDNA_40_120cm","AZ__iDNA_120_200cm",
            "AZ_eDNA_0_40cm","AZ_eDNA_40_120cm","AZ_eDNA_120_200cm",
            "SG__iDNA_0_40cm","SG__iDNA_40_120cm","SG__iDNA_120_200cm",
            "SG_eDNA_0_40cm","SG_eDNA_40_120cm","SG_eDNA_120_200cm",
            "LC__iDNA_0_40cm","LC__iDNA_40_120cm","LC__iDNA_120_200cm",
            "LC_eDNA_0_40cm","LC_eDNA_40_120cm","LC_eDNA_120_200cm",
            "NB__iDNA_0_40cm","NB__iDNA_40_120cm","NB__iDNA_120_200cm",
            "NB_eDNA_0_40cm","NB_eDNA_40_120cm","NB_eDNA_120_200cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,
                            levels = tmp123,
                            labels = c("Ai0-40cm","Ai40-120cm","Ai120-200cm",
                                       "Ae0-40cm","Ae40-120cm","Ae120-200cm",
                                       "Si0-40cm","Si40-120cm","Si120-200cm",
                                       "Se0-40cm","Se40-120cm","Se120-200cm",
                                       "Li0-40cm","Li40-120cm","Li120-200cm",
                                       "Le0-40cm","Le40-120cm","Le120_200cm",
                                       "Ni0-40cm","Ni40-120cm","Ni120-200cm",
                                       "Ne0-40cm","Ne40-120cm","Ne120-200cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]

# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
# samples-to-highlight <- c("AZiDNA0-5cm","AZiDNA5-10cm","AZiDNA10-20cm","AZiDNA20-40cm",
#                             "LCiDNA0-5cm","LCiDNA10-20cm","LCiDNA20-40cm","LCiDNA5-10cm")
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/04 Species structure barplot/top 50 ASV for allV2.pdf",
    width = 14,
    height = 10)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), breaks = c(0.1,1,4,8,16,24)) +
  labs(x = "Sample", y = "", title = "top 50 ASV") +
  guides(size = guide_legend(title = "Relative Abundance"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 02. top ASV for 0-40 sample-------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv),
         depth< 40) %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]

filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,] 
env2 <- env1 %>% 
  mutate(group_col = site_dna_depth)

group_info <- data.frame(col_name = rownames(env2),
                         group = env2$group_col)

asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", 
         value = "value", 
         -id) %>%  # Convert to long format
  left_join(group_info, 
            by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", 
         value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env2$group_col)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]

total_reads <- colSums(asv2.8)
relative_abundance <- 100*asv2.8 / total_reads

rankRA <- relative_abundance %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance))

top_50_asv <- rankRA %>%
  slice_head(n = 50) %>% 
  select(-mean_abundance)

filt5 <- rownames(tax) %in% rownames(top_50_asv)
tax_1 <- tax[filt5,]

# 为 asv-table 和 tax-table 添加行名
asv_table <- top_50_asv %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax_1 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)

unique(long_data$sample_id)

tmp123 <- c("AZiDNA0_5cm", "AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm",
            "AZeDNA0_5cm", "AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm",
            "SGiDNA0_5cm", "SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm",
            "SGeDNA0_5cm", "SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm",
            "LCiDNA0_5cm", "LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm",
            "LCeDNA0_5cm", "LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm",
            "NBiDNA0_5cm", "NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm",
            "NBeDNA0_5cm", "NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZiDNA0_5cm", "AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm",
                                       "AZeDNA0_5cm", "AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm",
                                       "SGiDNA0_5cm", "SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm",
                                       "SGeDNA0_5cm", "SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm",
                                       "LCiDNA0_5cm", "LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm",
                                       "LCeDNA0_5cm", "LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm",
                                       "NBiDNA0_5cm", "NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm",
                                       "NBeDNA0_5cm", "NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
# samples-to-highlight <- c("AZiDNA0-5cm","AZiDNA5-10cm","AZiDNA10-20cm","AZiDNA20-40cm",
#                             "LCiDNA0-5cm","LCiDNA10-20cm","LCiDNA20-40cm","LCiDNA5-10cm")
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/04 Species structure barplot/top50ASV_0-40V2.pdf",
    width = 14,
    height = 10)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 8), breaks = c(0.1, 1, 2, 4, 6)) +
  labs(x = "Sample", y = "", title = "top 50 ASV") +
  guides(size = guide_legend(title = "Relative Abundance"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 03. top ASV for 40_120 sample-------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv), depth > 40 & depth < 120) %>%
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]

filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,] 
env2 <- env1 %>% 
  mutate(group_col = site_dna_depth)

group_info <- data.frame(col_name = rownames(env2),group = env2$group_col)

asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", value = "value", -id) %>%  # Convert to long format
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env2$group_col)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]

total_reads <- colSums(asv2.8)
relative_abundance <- 100*asv2.8 / total_reads

rankRA <- relative_abundance %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance))

top_50_asv <- rankRA %>%
  slice_head(n = 50) %>% 
  select(-mean_abundance)

filt5 <- rownames(tax) %in% rownames(top_50_asv)
tax_1 <- tax[filt5,]

# 为 asv-table 和 tax-table 添加行名
asv_table <- top_50_asv %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax_1 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)

unique(long_data$sample_id)

tmp123 <- c("AZiDNA40_60cm", "AZiDNA60_80cm","AZiDNA80_100cm","AZiDNA100_120cm",
            "AZeDNA40_60cm", "AZeDNA60_80cm","AZeDNA80_100cm","AZeDNA100_120cm",
            "SGiDNA40_60cm", "SGiDNA60_80cm","SGiDNA80_100cm","SGiDNA100_120cm",
            "SGeDNA40_60cm", "SGeDNA60_80cm","SGeDNA80_100cm","SGeDNA100_120cm",
            "LCiDNA40_60cm", "LCiDNA60_80cm","LCiDNA80_100cm","LCiDNA100_120cm",
            "LCeDNA40_60cm", "LCeDNA60_80cm","LCeDNA80_100cm","LCeDNA100_120cm",
            "NBiDNA40_60cm", "NBiDNA60_80cm","NBiDNA80_100cm","NBiDNA100_120cm",
            "NBeDNA40_60cm", "NBeDNA60_80cm","NBeDNA80_100cm","NBeDNA100_120cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZiDNA40_60cm", "AZiDNA60_80cm","AZiDNA80_100cm","AZiDNA100_120cm",
                                       "AZeDNA40_60cm", "AZeDNA60_80cm","AZeDNA80_100cm","AZeDNA100_120cm",
                                       "SGiDNA40_60cm", "SGiDNA60_80cm","SGiDNA80_100cm","SGiDNA100_120cm",
                                       "SGeDNA40_60cm", "SGeDNA60_80cm","SGeDNA80_100cm","SGeDNA100_120cm",
                                       "LCiDNA40_60cm", "LCiDNA60_80cm","LCiDNA80_100cm","LCiDNA100_120cm",
                                       "LCeDNA40_60cm", "LCeDNA60_80cm","LCeDNA80_100cm","LCeDNA100_120cm",
                                       "NBiDNA40_60cm", "NBiDNA60_80cm","NBiDNA80_100cm","NBiDNA100_120cm",
                                       "NBeDNA40_60cm", "NBeDNA60_80cm","NBeDNA80_100cm","NBeDNA100_120cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
# samples-to-highlight <- c("AZiDNA0-5cm","AZiDNA5-10cm","AZiDNA10-20cm","AZiDNA20-40cm",
#                             "LCiDNA0-5cm","LCiDNA10-20cm","LCiDNA20-40cm","LCiDNA5-10cm")
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/04 Species structure barplot/top50ASV_40-120V2.pdf",
    width = 16,
    height = 10)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), breaks = c(0.1, 1, 4, 8,16,30)) +
  labs(x = "Sample", y = "", title = "top 50 ASV") +
  guides(size = guide_legend(title = "Relative Abundance"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 04. top ASV for 120_200 sample-------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv), depth > 120 & depth < 200) %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]

filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,] 
env2 <- env1 %>% 
  mutate(group_col = site_dna_depth)

group_info <- data.frame(col_name = rownames(env2),group = env2$group_col)

asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", value = "value", -id) %>%  # Convert to long format
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env2$group_col)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]

total_reads <- colSums(asv2.8)
relative_abundance <- 100*asv2.8 / total_reads

rankRA <- relative_abundance %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance))

top_50_asv <- rankRA %>%
  slice_head(n = 50) %>% 
  select(-mean_abundance)

filt5 <- rownames(tax) %in% rownames(top_50_asv)
tax_1 <- tax[filt5,]

# 为 asv-table 和 tax-table 添加行名
asv_table <- top_50_asv %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax_1 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)

unique(long_data$sample_id)

tmp123 <- c("AZiDNA120_140cm", "AZiDNA140_160cm","AZiDNA160_180cm","AZiDNA180_200cm",
            "AZeDNA120_140cm", "AZeDNA140_160cm","AZeDNA160_180cm","AZeDNA180_200cm",
            "SGiDNA120_140cm", "SGiDNA140_160cm","SGiDNA160_180cm","SGiDNA180_200cm",
            "SGeDNA120_140cm", "SGeDNA140_160cm","SGeDNA160_180cm","SGeDNA180_200cm",
            "LCiDNA120_140cm", "LCiDNA140_160cm","LCiDNA160_180cm","LCiDNA180_200cm",
            "LCeDNA120_140cm", "LCeDNA140_160cm","LCeDNA160_180cm","LCeDNA180_200cm",
            "NBiDNA120_140cm", "NBiDNA140_160cm","NBiDNA160_180cm","NBiDNA180_200cm",
            "NBeDNA120_140cm", "NBeDNA140_160cm","NBeDNA160_180cm","NBeDNA180_200cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZiDNA120_140cm", "AZiDNA140_160cm","AZiDNA160_180cm","AZiDNA180_200cm",
                                       "AZeDNA120_140cm", "AZeDNA140_160cm","AZeDNA160_180cm","AZeDNA180_200cm",
                                       "SGiDNA120_140cm", "SGiDNA140_160cm","SGiDNA160_180cm","SGiDNA180_200cm",
                                       "SGeDNA120_140cm", "SGeDNA140_160cm","SGeDNA160_180cm","SGeDNA180_200cm",
                                       "LCiDNA120_140cm", "LCiDNA140_160cm","LCiDNA160_180cm","LCiDNA180_200cm",
                                       "LCeDNA120_140cm", "LCeDNA140_160cm","LCeDNA160_180cm","LCeDNA180_200cm",
                                       "NBiDNA120_140cm", "NBiDNA140_160cm","NBiDNA160_180cm","NBiDNA180_200cm",
                                       "NBeDNA120_140cm", "NBeDNA140_160cm","NBeDNA160_180cm","NBeDNA180_200cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
# samples-to-highlight <- c("AZiDNA0-5cm","AZiDNA5-10cm","AZiDNA10-20cm","AZiDNA20-40cm",
#                             "LCiDNA0-5cm","LCiDNA10-20cm","LCiDNA20-40cm","LCiDNA5-10cm")
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/04 Species structure barplot/top50ASV_120-200V2.pdf",
    width = 16,
    height = 10)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), breaks = c(0.1, 1, 4, 8,15)) +
  labs(x = "Sample", y = "", title = "top 50 ASV") +
  guides(size = guide_legend(title = "Relative Abundance"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()