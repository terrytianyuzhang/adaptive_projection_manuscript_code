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

group1_data <- pooled_data[1:nrow(group1_data),]
group0_data <- pooled_data[(nrow(group1_data)+ 1):nrow(pooled_data),]


result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_2.rds')
result_list <- readRDS(file = result_file_name)


###extract the estimated PC vector
PC1_index <- 1
PC2_index <- 4
PC_benchmark_1 <- extract_pc(result_list$debiased_pc_result)[[2]][,PC1_index]
PC_benchmark_2 <- extract_pc(result_list$debiased_pc_result)[[2]][,PC2_index]

pc_plot_df <- data.frame()
for(split_index in 1:length(result_list$debiased_pc_result$split_data)){
  group1_data_subset <- group1_data[result_list$debiased_pc_result$split_data[[split_index]]$sample_1_cross_index, ]
  group0_data_subset <- group0_data[result_list$debiased_pc_result$split_data[[split_index]]$sample_2_cross_index, ]
  
  PC_vector <- extract_pc(result_list$debiased_pc_result)[[split_index]][,PC1_index]
  PC_vector <- as.numeric(sign(crossprod(PC_vector, PC_benchmark_1))) * PC_vector
  group1_score_1 <- as.matrix(group1_data_subset) %*% PC_vector
  group0_score_1 <- as.matrix(group0_data_subset) %*% PC_vector
  
  PC_vector <- extract_pc(result_list$debiased_pc_result)[[split_index]][,PC2_index]
  PC_vector <- as.numeric(sign(crossprod(PC_vector, PC_benchmark_2))) * PC_vector
  group1_score_2 <- as.matrix(group1_data_subset) %*% PC_vector
  group0_score_2 <- as.matrix(group0_data_subset) %*% PC_vector

  pc_plot_df <- rbind(pc_plot_df,
                      data.frame(Group = as.factor(c(rep(1, nrow(group1_data_subset)), rep(0, nrow(group0_data_subset)))),
                                 split_index = as.factor(c(rep(split_index, nrow(group1_data_subset)), rep(split_index, nrow(group0_data_subset)))),
                           PC1 = c(group1_score_1, group0_score_1),
                           PC4 = c(group1_score_2, group0_score_2)))
}
pc_plot_df <- data.table(pc_plot_df)
pc_plot_df[Group == 1, Group := 'Case']
pc_plot_df[Group == 0, Group := 'Control']
pc_plot <- ggplot(data = pc_plot_df) + 
  geom_point(aes(x = PC1, y = PC4, group = Group, color = Group, shape = Group),size = 3, alpha = 0.8)+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))
  

########

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

gene_order <- unique(c(pc1_gene_name,pc2_gene_name,pc3_gene_name,pc4_gene_name,pc5_gene_name,pc6_gene_name,pc7_gene_name, pc8_gene_name,pc9_gene_name,pc10_gene_name,colnames(pooled_data)))[1:200]
gene_order <- rev(gene_order)
plot_gene_data <- data.frame(pooled_data[,gene_order])


plot_gene_data_cor <- cor(plot_gene_data[sapply(plot_gene_data,is.numeric)])
heatmap_data <- melt(plot_gene_data_cor)
color_limit <- c(min(heatmap_data$value),1)
color_limit <- c(-1,1)
names(heatmap_data)[3] <- 'cor'
overall_plot <- ggplot(heatmap_data,
                       aes(ordered(Var1, levels = rev(sort(unique(Var1), decreasing = F))), Var2, fill = cor)
                       # aes(x=Var1,
                       #     y=Var2,
                       #     fill=cor)
                       )+
                geom_tile()+
  
                # scale_fill_viridis(discrete=FALSE) +
                scico::scale_fill_scico(name = 'gene-gene correlation ', palette = "vik", limits=color_limit,breaks = c(-1, -0.5,0,0.5,1), labels = c('-1', '-0.5', '0','0.5','1'))+
                theme(plot.margin = margin(1, 1, 1, 1, "cm"),
                      axis.text.x = element_text(angle = 90, hjust = 1),
                      legend.position = 'bottom')+
                xlab('')+
                ylab('')+
                scale_x_discrete(breaks = NULL)+
                scale_y_discrete(breaks = NULL)
                
overall_plot
ggarrange(overall_plot, pc_plot, nrow = 1, labels = c('A', 'B'))

motivation_plot_file <- paste0(jinhong_folder, 
                            '/result/motivation_', cell_type,'.pdf')
ggsave(file = motivation_plot_file, width = 9, height = 5)



result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_1.rds')
result_list <- readRDS(file = result_file_name)
effect_gene_name <- summarize_feature_name(result_list$simple_lasso_test_result, method = 'majority voting')
print(effect_gene_name)


lasso_plot_df <- data.frame()
for(split_index in 1:length(result_list$debiased_pc_result$split_data)){
  lasso_vector <- extract_lasso_coef(result_list$simple_lasso_test_result)[[split_index]]
  
  group1_data_subset <- group1_data[result_list$simple_lasso_test_result$split_data[[split_index]]$sample_1_cross_index, ]
  group0_data_subset <- group0_data[result_list$simple_lasso_test_result$split_data[[split_index]]$sample_2_cross_index, ]
  group1_score_lasso <- as.matrix(group1_data_subset) %*% lasso_vector
  group0_score_lasso <- as.matrix(group0_data_subset) %*% lasso_vector
  
  lasso_plot_df <- rbind(lasso_plot_df,
                         data.frame(group_index = as.factor(c(rep(1, nrow(group1_data_subset)), rep(0, nrow(group0_data_subset)))),
                              score_lasso = c(group1_score_lasso, group0_score_lasso)))
}
lasso_plot_df$PC1 <- pc_plot_df$PC1
lasso_plot_df <- data.table(lasso_plot_df)
lasso_plot_df[group_index == 1, group_index := 'Case']
lasso_plot_df[group_index == 0, group_index := 'Control']
ggplot(data = lasso_plot_df) + 
  geom_point(aes(x = PC1, y = score_lasso, group = group_index, color = group_index, shape = group_index))+
  theme_minimal()
