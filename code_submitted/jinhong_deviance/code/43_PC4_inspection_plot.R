library(ggplot2)
library(ggpubr)
library(data.table)
library(org.Hs.eg.db)
library(PMA)
library(MASS)
library(highmean)
library(GO.db)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)

jinhong_folder <- './jinhong_deviance/'

cell_type <- 'T4'
result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_2.rds')
result_list <- readRDS(file = result_file_name)


pc1_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                  latent_fator_index = 1)
pc2_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 2)
pc3_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 3)
pc4_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 4)
pc5_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 5)
pc6_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 6)
pc7_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 7)
pc8_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 8)
pc9_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 9)
pc10_gene_name <- summarize_pc_name(result_list$debiased_pc_result, method = 'majority voting',
                                   latent_fator_index = 10)


result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_1.rds')
result_list <- readRDS(file = result_file_name)
effect_gene_name <- summarize_feature_name(result_list$simple_lasso_test_result, method = 'majority voting')
print(effect_gene_name)


raw_gene <- data.table(read.csv(paste0(jinhong_folder, '/data/resid_', cell_type,'.csv')))
raw_label <- data.table(read.csv(paste0(jinhong_folder, '/data/case_ctrl_', cell_type,'.csv')))
raw_merged <- merge(raw_label, raw_gene, by = 'X')
rownames(raw_merged) <- raw_merged$X
raw_merged <- raw_merged[, -1]
colnames(raw_merged)[1] <- 'case'
# raw_merged$case <- sample(raw_merged$case)

group1_data <- raw_merged[case == 1,]
group1_data <- group1_data[, -1]
group0_data <- raw_merged[case == 0,]
group0_data <- group0_data[, -1]
group0_data <- data.frame(group0_data)
group1_data <- data.frame(group1_data)

pooled_data <- rbind(group1_data, group0_data)
pooled_data <- as.data.frame(apply(pooled_data,2,mynorm <- function(x){return (x/sd(x))}))

# group0_data <- as.data.frame(apply(group0_data,2,mynorm <- function(x){return (x/sd(x))}))
# group1_data <- as.data.frame(apply(group1_data,2,mynorm <- function(x){return (x/sd(x))}))
# pooled_data <- group1_data


group1_data <- pooled_data[1:nrow(group1_data),]
group0_data <- pooled_data[(nrow(group1_data)+ 1):nrow(pooled_data),]


library(ggplot2)
library(reshape2)
library(ggpubr)
# library(viridis)
# library(hrbrthemes)

####overall plot
gene_order <- unique(c(pc1_gene_name,pc2_gene_name,pc3_gene_name,pc4_gene_name,pc5_gene_name,pc6_gene_name,pc7_gene_name, pc8_gene_name,pc9_gene_name,pc10_gene_name,colnames(pooled_data)))[1:200]
gene_order <- rev(gene_order)
plot_gene_data <- data.frame(pooled_data[,gene_order])


plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)
color_limit <- c(-1,1)
names(heatmap_data)[3] <- 'cor'
overall_plot <- ggplot(heatmap_data,
                       aes(x=Var1,
                           y=Var2,
                           fill=cor))+geom_tile()+
  # scale_fill_viridis(discrete=FALSE) +
  scico::scale_fill_scico(palette = "vik", limits=color_limit,breaks = c(-1, -0.5,0,0.5,1), labels = c('-1', '-0.5', '0','0.5','1'))+
  theme(plot.margin = margin(2.5, 1, 2.5, 1, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('gene 1')+
  ylab('gene 2')+
  scale_x_discrete(breaks = NULL)+
  scale_y_discrete(breaks = NULL)
  # geom_rect(xmin = 182, xmax = 200, ymin = 182, ymax = 200,color = "#6c4dee", fill = NA, size = 1 , linejoin = 'round')
####

####PC4 zoom-in plot
gene_order <- intersect(pc4_gene_name, colnames(pooled_data))
# gene_order <- intersect(unique(c(pc_gene_name, colnames(pooled_data))), colnames(pooled_data))[1:21]

gene_order <- rev(gene_order)
plot_gene_data <- data.frame(pooled_data[,gene_order])

plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)
names(heatmap_data)[3] <- 'cor'
zoom_in_pc4 <- ggplot(heatmap_data,
                       aes(x=Var1,
                           y=Var2,
                           fill=cor))+geom_tile()+
  # scale_fill_viridis(discrete=FALSE) +
  scico::scale_fill_scico(palette = "vik", limits=color_limit)+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab(NULL)+
  ylab(NULL)

gene_order <- gene_order[!grepl("^MT\\.", gene_order)]
plot_gene_data <- data.frame(pooled_data[,gene_order])

plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)
names(heatmap_data)[3] <- 'cor'
zoom_in_pc4_without_MT <- ggplot(heatmap_data,
                       aes(x=Var1,
                           y=Var2,
                           fill=cor))+geom_tile()+
  # scale_fill_viridis(discrete=FALSE) +
  scico::scale_fill_scico(palette = "vik", limits=color_limit)+
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab(NULL)+
  ylab(NULL)

ggarrange(zoom_in_pc4, zoom_in_pc4_without_MT, 
          nrow = 1, common.legend = TRUE, legend = 'bottom', labels = c('A', 'B'))
PC4_inspection_plot_file <- paste0(jinhong_folder, 
                            '/result/PC4_inspection.pdf')
ggsave(file = PC4_inspection_plot_file, width = 10, height = 5)
