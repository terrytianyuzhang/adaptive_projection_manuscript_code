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
library(ggplot2)
library(reshape2)
library(ggpubr)

jinhong_folder <- './jinhong_deviance/'
source('./code_paper/supporting_function.R')

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

####zoom in pc plot
gene_order <- intersect(pc1_gene_name, colnames(pooled_data))
# gene_order <- intersect(unique(c(pc_gene_name, colnames(pooled_data))), colnames(pooled_data))[1:21]

gene_order <- rev(gene_order)
plot_gene_data <- data.frame(pooled_data[,gene_order])

plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)

zoom_in_plot <- ggplot(heatmap_data,
                       aes(x=Var1,
                           y=Var2,
                           fill=value))+geom_tile()+
  # scale_fill_viridis(discrete=FALSE) +
  scico::scale_fill_scico(palette = "vik", limits=color_limit,breaks = c(-1, -0.5,0,0.5,1), labels = c('-1', '-0.5', '0','0.5','1'))+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')+
  xlab(NULL)+
  ylab(NULL)

###bar plot for PC1
gene_order <- intersect(pc1_gene_name, colnames(pooled_data))
# gene_order <- intersect(unique(c(pc_gene_name, colnames(pooled_data))), colnames(pooled_data))[1:21]

gene_order <- rev(gene_order)
pc_bar_plot <- data.frame(Gene = gene_order,
                          PC1_load = -extract_pc(result_list$debiased_pc_result)[[3]][gene_order,1])
pc_bar_plot <- pc_bar_plot[order(pc_bar_plot$PC1_load, decreasing = T),]
pc_bar_plot$Gene <- factor(pc_bar_plot$Gene, levels = pc_bar_plot$Gene)
PC1_load_bar <- ggplot(data = pc_bar_plot)+
  geom_bar(aes(x = Gene, y = PC1_load,  fill = PC1_load),color = "#FFFFFF",stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+
  scale_fill_gradient(name = ' ', low = "#D1ADEF", high = "#3F0020", na.value = NA)+
  xlab('')+
  ylab('PC1 Loading')+
  ylim(c(0,0.5))

####lasso plot
effect_gene_name <- c( 'ARRDC3',"RPS27","STAT1","B2M", "LTB")
gene_order <- intersect(effect_gene_name, colnames(pooled_data))
# gene_order <- intersect(poisson_names, colnames(pooled_data))
# gene_order <- gene_order[grepl("^I", gene_order) | grepl("^E", gene_order) ]

# gene_order <- rev(gene_order)
plot_gene_data <- data.frame(pooled_data[,gene_order])

plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)
# heatmap(plot_gene_data_cor, symm = T)

lasso_plot <- ggplot(heatmap_data,
                     aes(x=Var1,
                         y=Var2,
                         fill=value))+geom_tile()+
  # scale_fill_viridis(discrete=FALSE) +
  scico::scale_fill_scico(palette = "vik", limits=color_limit)+
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none')+
  xlab(NULL)+
  ylab(NULL)

###HISTOGRAM FOR PC1 SCORE
result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_2.rds')
result_list <- readRDS(file = result_file_name)

PC1_index <- 1
PC_benchmark_1 <- extract_pc(result_list$debiased_pc_result)[[2]][,PC1_index] ##make sure all the pc estimates are of roughly the same direction

pc_plot_df <- data.frame()
for(split_index in 1:length(result_list$debiased_pc_result$split_data)){
  group1_data_subset <- group1_data[result_list$debiased_pc_result$split_data[[split_index]]$sample_1_cross_index, ]
  group0_data_subset <- group0_data[result_list$debiased_pc_result$split_data[[split_index]]$sample_2_cross_index, ]
  
  PC_vector <- extract_pc(result_list$debiased_pc_result)[[split_index]][,PC1_index]
  PC_vector <- as.numeric(sign(crossprod(PC_vector, PC_benchmark_1))) * PC_vector
  group1_score_1 <- as.matrix(group1_data_subset) %*% PC_vector
  group0_score_1 <- as.matrix(group0_data_subset) %*% PC_vector
  
  pc_plot_df <- rbind(pc_plot_df,
                      data.frame(Group = as.factor(c(rep(1, nrow(group1_data_subset)), rep(0, nrow(group0_data_subset)))),
                                 split_index = as.factor(c(rep(split_index, nrow(group1_data_subset)), rep(split_index, nrow(group0_data_subset)))),
                                 pc_score = c(group1_score_1, group0_score_1)))
}
pc_plot_df <- data.table(pc_plot_df)
pc_plot_df[Group == 1, Group := 'Case']
pc_plot_df[Group == 0, Group := 'Control']
pc_score_plot <- ggplot(data = pc_plot_df, aes(x = pc_score)) + 
  geom_histogram(aes(y = ..density.., group = Group, fill = Group), position = 'identity', alpha = 0.5, color = 'black')+
  # geom_density(aes(group = group_index), alpha = 0.2, color = '#3333CC')+
  theme_minimal()+
  xlab(paste0('PC', PC1_index,' Score'))+
  ylab('Density')+
  scale_fill_manual(name = 'Group', values = c("#00203F", "#ADEFD1", "grey"))+
  theme(legend.position = 'bottom')+
  xlim(c(-3*sd(pc_plot_df$pc_score), 3*sd(pc_plot_df$pc_score)))

####histogram for lasso score
result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_1.rds')
result_list <- readRDS(file = result_file_name)
effect_gene_name <- summarize_feature_name(result_list$simple_lasso_test_result, method = 'majority voting')
print(effect_gene_name)


lasso_plot_df <- data.frame()
for(split_index in 1:length(result_list$debiased_pc_result$split_data)){
  lasso_vector <- -extract_lasso_coef(result_list$simple_lasso_test_result)[[split_index]]
  
  group1_data_subset <- group1_data[-result_list$simple_lasso_test_result$split_data[[split_index]]$sample_1_cross_index, ]
  group0_data_subset <- group0_data[-result_list$simple_lasso_test_result$split_data[[split_index]]$sample_2_cross_index, ]
  group1_score_lasso <- as.matrix(group1_data_subset) %*% lasso_vector
  group0_score_lasso <- as.matrix(group0_data_subset) %*% lasso_vector
  
  lasso_plot_df <- rbind(lasso_plot_df,
                         data.frame(group_index = as.factor(c(rep(1, nrow(group1_data_subset)), rep(0, nrow(group0_data_subset)))),
                                    score_lasso = c(group1_score_lasso, group0_score_lasso)))
}
lasso_plot_df <- data.table(lasso_plot_df)
lasso_plot_df[group_index == 1, group_index := 'Case']
lasso_plot_df[group_index == 0, group_index := 'Control']
lasso_score_plot <- ggplot(data = lasso_plot_df, aes(x = score_lasso)) + 
  geom_histogram(aes(y = ..density.., group = group_index, fill = group_index), position = 'identity', alpha = 0.5, color = 'black')+
  # geom_density(aes(group = group_index), alpha = 0.2, color = '#3333CC')+
  theme_minimal()+
  xlab('Lasso Score')+
  ylab('Density')+
  scale_fill_manual(name = 'Group', values = c("#00203F", "#ADEFD1", "grey"))+
  theme(legend.position = 'bottom')+
  xlim(c(-3*sd(lasso_score_plot$score_lasso), 3*sd(lasso_score_plot$score_lasso)))


right_panel <- ggarrange(PC1_load_bar, pc_score_plot, lasso_score_plot, labels = c('C', 'D', 'E'), nrow = 3, ncol = 1)
left_panel <- ggarrange(zoom_in_plot, lasso_plot, labels = c('A', 'B'),
          nrow = 2, ncol = 1, common.legend = FALSE)
ggarrange(left_panel, right_panel, nrow = 1, ncol = 2)
interpretation_plot_file <- paste0(jinhong_folder, 
                            '/result/interpretation_plot.pdf')
ggsave(file = interpretation_plot_file, width = 10, height = 10)
  # geom_rect(xmin = 182, xmax = 200, ymin = 182, ymax = 200,color = "#6c4dee", fill = NA, size = 1 , linejoin = 'round')
####
