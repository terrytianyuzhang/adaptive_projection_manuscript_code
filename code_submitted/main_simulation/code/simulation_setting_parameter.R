if(SIMULATION_BATCH %in% c('InfGlobal', 'InfGlobal_hard', 'InfGlobal_hard_no_truncation')){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1
  non_zero_mean_2_part2 <- 0
}

if(SIMULATION_BATCH %in% c('InfProject', 'InfProject_hard', 'InfProject_hard_no_truncation')){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1
  non_zero_mean_2_part2 <- 5
}

if(SIMULATION_BATCH %in% c('InfAlternative', 'InfAlternative_hard', 'InfAlternative_hard_no_truncation')){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1.2
  non_zero_mean_2_part2 <- 0.9
}

if(SIMULATION_BATCH %in% c('InfAlternative_hard_weak')){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1.2
  non_zero_mean_2_part2 <- 0.5 ###changed this from InfAlternative_hard
}


if(SIMULATION_BATCH %in% c('NormalGlobal', 'NormalGlobal_hard')){
  num_features <- c(100)
  num_repeat <- 1000
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- non_zero_mean_1_part2 <- non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}

if(SIMULATION_BATCH %in% c('NormalProject', 'NormalProject_hard')){
  num_features <- c(100)
  num_repeat <- 1000
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- non_zero_mean_1_part2 <- 5
  non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}

if(SIMULATION_BATCH %in% c('NormalAlternative', 'NormalAlternative_hard')){
  num_features <- c(100)
  num_repeat <- 1000
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- 1.1
  non_zero_mean_1_part2 <- 0.9
  non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}

if(SIMULATION_BATCH == 'NormalGlobal_Chen2010'){
  num_features <- c(100)
  num_repeat <- 1e3
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- non_zero_mean_1_part2 <- non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}

if(SIMULATION_BATCH == 'NormalProject_Chen2010'){
  num_features <- c(100)
  num_repeat <- 1e3
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- non_zero_mean_1_part2 <- 5
  non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}

if(SIMULATION_BATCH == 'NormalAlternative_Chen2010'){
  num_features <- c(100)
  num_repeat <- 1e3
  non_zero_number_features <- c(10)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  pc_weight <- 3
  
  non_zero_mean_1_part1 <- 1.1
  non_zero_mean_1_part2 <- 0.9
  non_zero_mean_2_part1 <- non_zero_mean_2_part2 <- 1
}


if(SIMULATION_BATCH == 'InfGlobal_Chen2010'){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1
  non_zero_mean_2_part2 <- 0
}

if(SIMULATION_BATCH == 'InfProject_Chen2010'){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1
  non_zero_mean_2_part2 <- 5
}

if(SIMULATION_BATCH == 'InfAlternative_Chen2010'){
  mask_probability <- 0.5
  
  num_features <- c(1000)
  num_repeat <- 1e3
  non_zero_number_features <- c(20)
  candidate_sample_size_1 <- c(100, 300, 500)
  noise_variance <- 1
  
  non_zero_mean_1 <- 1
  non_zero_mean_2_part1 <- 1.2
  non_zero_mean_2_part2 <- 0.9
}

# 
# if(SIMULATION_BATCH == '1'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   # candidate_sample_size_2 is the same as candidate_sample_size_1
#   # candidate_feature_number is candidate_sample_size_1/2, candidate_sample_size_1, and candidate_sample_size_1*2
#   # number of non-zero feature for pc is feature_number/2, feature_number^(1/2), feature_number^(1/3)
#   num_repeat <- 20
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 4.96
#   pc1_non_zero <- 1
#   pca_method <- 'sparse_pca'
# }else if(SIMULATION_BATCH == '2'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   # candidate_sample_size_2 is the same as candidate_sample_size_1
#   # candidate_feature_number is candidate_sample_size_1/2, candidate_sample_size_1, and candidate_sample_size_1*2
#   # number of non-zero feature for pc is feature_number/2, feature_number^(1/2), feature_number^(1/3)
#   num_repeat <- 30
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 4.96
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
# }else if(SIMULATION_BATCH == '3'){
#   candidate_sample_size_1 <- c(50, 200, 400)
#   # candidate_sample_size_2 is the same as candidate_sample_size_1
#   # candidate_feature_number is candidate_sample_size_1/2, candidate_sample_size_1, and candidate_sample_size_1*2
#   # number of non-zero feature for pc is feature_number/2, feature_number^(1/2), feature_number^(1/3)
#   num_repeat <- 200
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 5
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
# }else if(SIMULATION_BATCH == '4'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   # candidate_sample_size_2 is the same as candidate_sample_size_1
#   # candidate_feature_number is candidate_sample_size_1, candidate_sample_size_1*3
#   # number of non-zero feature for pc is feature_number^(1/2), feature_number^(1/3)
#   # number of pc non-zero feature is 0.6 of  non-zero feature for mean diff 
#   num_repeat <- 30
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 4.96
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
#   mean_method <- 'lasso'
# }else if(SIMULATION_BATCH == '5'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   # candidate_sample_size_2 is the same as candidate_sample_size_1
#   # candidate_feature_number is candidate_sample_size_1, candidate_sample_size_1*3
#   # number of non-zero feature for pc is feature_number^(1/2), feature_number^(1/3)
#   # number of pc non-zero feature is 0.6 of  non-zero feature for mean diff 
#   num_repeat <- 200
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 5
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
#   mean_method <- 'lasso'
# }else if(SIMULATION_BATCH == '6'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 4.96
#   spike_signal <- 4.8
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
#   mean_method <- 'lasso'
# }else if(SIMULATION_BATCH == '7'){
#   candidate_sample_size_1 <- c(100, 200, 300, 400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 200
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 4.95
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
#   mean_method <- 'lasso'
# }else if(SIMULATION_BATCH == '8'){
#   candidate_sample_size_1 <- c(50, 100, 200, 400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 400
#   noise_variance <- 0.1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 5
#   pc1_non_zero <- 1
#   pca_method <- 'dense_pca'
#   mean_method <- 'lasso'
# }else if(SIMULATION_BATCH == '9'){
#   candidate_sample_size_1 <- c(100, 200, 300, 400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 1
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'naive'
#   num_features <- c(1000)
# }else if(SIMULATION_BATCH == '10'){
#   candidate_sample_size_1 <- c(100, 200, 300, 400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'naive'
#   num_features <- c(1000)
# }else if(SIMULATION_BATCH == '11'){
#   candidate_sample_size_1 <- c(100, 400, 800)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 1
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'naive'
#   num_features <- c(1000)
# }else if(SIMULATION_BATCH == '12'){
#   candidate_sample_size_1 <- c(100, 400, 800)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'naive'
#   num_features <- c(1000)
# }else if(SIMULATION_BATCH == '13'){
#   candidate_sample_size_1 <- c(100, 400, 800)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 1
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   num_features <- c(1000)
# }else if(SIMULATION_BATCH == '14'){
#   candidate_sample_size_1 <- c(100, 400, 1600, 6400)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   mean_2_non_zero <- 1
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   num_features <- c(100)
# }else if(SIMULATION_BATCH == '15'){
#   candidate_sample_size_1 <- seq(150, 700, by = 50)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   sparse_pca_method <- 'sparse_pca'
#   num_features <- c(1000)
#   non_zero_number_features <- c(15)
# }else if(SIMULATION_BATCH == '16'){
#   candidate_sample_size_1 <- seq(150, 700, by = 50)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 1
#   mean_1_non_zero <- 0.1
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'naive'
#   sparse_pca_method <- 'sparse_pca'
#   num_features <- c(1000)
#   non_zero_number_features <- c(15)
# }else if(SIMULATION_BATCH == '17'){
#   candidate_sample_size_1 <- seq(50, 1000, by = 50)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   sparse_pca_method <- 'sparse_pca'
#   num_features <- c(100)
#   non_zero_number_features <- c(5)
# }else if(SIMULATION_BATCH == '18'){
#   candidate_sample_size_1 <- seq(150, 700, by = 50)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   sparse_pca_method <- 'sparse_pca'
#   num_features <- c(1000)
#   non_zero_number_features <- c(5)
# }else if(SIMULATION_BATCH == '19'){
#   ##this is similar to setting 15
#   
#   candidate_sample_size_1 <- seq(150, 700, by = 50)
#   #this is basically SIMULATION_BATCH = 4, except for an additionally spike signal
#   num_repeat <- 30
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   sparse_pca_method <- 'hard'
#   num_features <- c(1000)
#   non_zero_number_features <- c(15)
# }else if(SIMULATION_BATCH == '20'){
#   ##this is similar to setting 17
#   
#   candidate_sample_size_1 <- seq(50, 1000, by = 50)
#   
#   num_repeat <- 100
#   noise_variance <- 1
#   mean_1_non_zero <- 5
#   pc1_non_zero <- 1
#   portion_of_pc <- 5 
#   sample_mean_method <- 'hard'
#   sparse_pca_method <- 'hard'
#   num_features <- c(100)
#   non_zero_number_features <- c(5)
# }
