rm(list = ls())
library(ggplot2)
library(ggpubr)
library(data.table)
library(VennDiagram)
library(gridExtra)
library(UpSetR)

jinhong_folder <- './jinhong_deviance/'

source('./code_paper/v2high_dimensional_mean_testing_function.R')
source('./code_paper/supporting_function.R')

cell_types <- c('T4', 'T8', 'cM', 'B', 'NK')
all_structured_result <- data.table()
upset_pre_name <- c("All" , "PC1 " , 'glm', "Lasso")

  
  upset_plot_list <- list()
  upset_plot_index <- 0
  for(cell_type in cell_types){
    
        upset_plot_index <- upset_plot_index + 1
        
        result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_2.rds')
        result_list <- readRDS(file = result_file_name)
        
        all_gene_names <- rownames(result_list$debiased_pc_result$split_data[[1]]$nuisance_collection$estimate_leading_pc)
        gene_in_pc1 <- summarize_pc_name(result_list$debiased_pc_result,latent_fator_index = 1,method = 'majority voting')
        
        gene_in_pc2 <- summarize_pc_name(result_list$debiased_pc_result,latent_fator_index = 2,method = 'majority voting')
        gene_in_pc3 <- summarize_pc_name(result_list$debiased_pc_result,latent_fator_index = 3,method = 'majority voting')
        gene_in_pc4 <- summarize_pc_name(result_list$debiased_pc_result,latent_fator_index = 4,method = 'majority voting')
        
        ###
        result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_1.rds')
        result_list_lasso <- readRDS(file = result_file_name)
        gene_in_lasso <- summarize_feature_name(result_list_lasso$simple_lasso_test_result, method = 'majority voting')
        
        nb_glm_file <-  paste0(jinhong_folder,'data/pairwise_test/test_nb_', cell_type,'.csv')
        nb_glm <- data.table(read.csv(nb_glm_file))
        nb_names <- nb_glm[p_values < 0.05,]$X
        nb_names <- intersect(nb_names, all_gene_names)
        
        nb_names_adjusted <- nb_glm[p_values < 0.05/nrow(nb_glm),]$X
        nb_names_adjusted <- intersect(nb_names_adjusted, all_gene_names)
        
        # names_for_upset <- list(all_gene_names, gene_in_pc1, gene_in_pc2, poisson_names, gene_in_lasso)
        names_for_upset <- list(all_gene_names, gene_in_pc1, gene_in_pc2, nb_names, gene_in_lasso, nb_names_adjusted)
        
        
        structured_result <- data.table(cell_type = cell_type,
                                        simple_p = 2*pnorm(-abs(result_list$simple_pc_result$test_statistics[1])),
                                        simple_p_2 = 2*pnorm(-abs(result_list$simple_pc_result$test_statistics[2])),
                                        simple_p_3 = 2*pnorm(-abs(result_list$simple_pc_result$test_statistics[3])),
                                        simple_p_4 = 2*pnorm(-abs(result_list$simple_pc_result$test_statistics[4])),
                                        simple_p_5 = 2*pnorm(-abs(result_list$simple_pc_result$test_statistics[5])),
                                        debias_p = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[1])),
                                        debias_p_2 = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[2])),
                                        debias_p_3 = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[3])),
                                        debias_p_4 = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[4])),
                                        debias_p_5 = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[5])),
                                        # debias_p_2 = 2*pnorm(-abs(result_list$debiased_pc_result$test_statistics[2])),
                                        lasso_p = 2*pnorm(-abs(result_list_lasso$simple_lasso_test_result$test_statistics[1])))
                                        # chen2010_p = result_list$chen2010_result$pval,
                                        # t_p = min(unlist(result_list$t_test_result$p.value)) * nrow(result_list$t_test_result))
        # t_p = min(unlist(result_list$t_test_result$p.value)) )
        all_structured_result <- rbind(all_structured_result, structured_result)

        listInput <- list(
          Lasso = as.character(gene_in_lasso),
          PC1 = as.character(gene_in_pc1),
          PC2 = as.character(gene_in_pc2),
          PC3 = as.character(gene_in_pc3),
          PC4 = as.character(gene_in_pc4),
          GLM_NB = as.character(nb_names),
          GLM_NB_ADJ = as.character(nb_names_adjusted)
          )
      
        my_format_p_value <- function(numerical_p){
          candidate <- paste0('=',sprintf(numerical_p, fmt = '%#.3f'))
          if(candidate == '=0.000'){
            candidate <- '<0.001'
          }
          return(candidate)
        }
        
        names(listInput)[1] <- paste0('Lasso(p', my_format_p_value(structured_result$lasso_p),')')
        names(listInput)[2] <- paste0('PC1(p', my_format_p_value(structured_result$debias_p) ,')')
        names(listInput)[3] <- paste0('PC2(p', my_format_p_value(structured_result$debias_p_2),')')
        names(listInput)[4] <- paste0('PC3(p', my_format_p_value(structured_result$debias_p_3),')')
        names(listInput)[5] <- paste0('PC4(p', my_format_p_value(structured_result$debias_p_4),')')
        
        
        pdf(paste0(jinhong_folder,
                   '/result/upset_', cell_type, '_with_negative_binomial.pdf'), width = 5, height = 3, onefile = FALSE)
        upset_plot <- upset(fromList(listInput), nsets = 8, text.scale = 0.8, keep.order = T, sets = rev(names(listInput)))
        print(upset_plot)
        grid.text(paste0(cell_type, ' cell (', length(all_gene_names), ' genes)'), hjust = 1.9, vjust = -16, gp = gpar(fontsize = 8))
        
        dev.off()
        
        # print(length(extract_pc(result_list$debiased_pc_result)[[3]][,1]))
    
  }

