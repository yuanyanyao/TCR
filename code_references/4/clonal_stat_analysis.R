library(ggplot2)
library(magrittr)
library(dplyr)
library(stats)

tcr1_cluster_list <- readRDS("/Users/fahadparyani/Downloads/tcr1_cluster_list.rds")
cl_patient<-do.call(rbind,sapply(tcr1_cluster_list, function(x) {table(as.character(x$patients1))}))
tab_cond_patient <- table(tcr1_cluster_list[1]$YQK_4_15$patients1, tcr1_cluster_list[1]$YQK_4_15$Condition)

variance <- apply(cl_patient, 1, var)

mean_pt <- apply(cl_patient,1,function(x){mean(x)})

mean_var_df <- data.frame('mean' = mean_pt,
                          'variance' = variance)

ggplot(mean_var_df, aes(x = mean, y = variance)) +
  geom_point()

condition_rand <- mean_pt > 1.5
condition_rand <- (mean_pt > 1 & variance > 2)

cl_patient_red <- cl_patient[condition_rand,]


#this just reorganizes the data so that it outputs
# the frequency of clonal type 
clonal_type_shuffle <- apply(cl_patient_red,2,function(col){
  
  clonal_type <- names(col)
  # print(clonal_type)
  seq <- mapply(function(cl_type, rep){
    rep(cl_type, rep)
  }, cl_type = clonal_type, rep = col) %>% 
    unlist() %>% 
    as.vector() %>% 
    sample() #shuffles the data
  
  return(seq)
})

store_sn_p_value <- list()
store_ctx_p_value <- list()
for(i in 1:100){
  print(i)
  hold_ct <- lapply(clonal_type_shuffle, function(x){sample(x,1000)}) %>% unlist() %>% table()
  
  condition_rand <- names(hold_ct[log(1+hold_ct) > 3])
  
  cl_patient_analysis <- cl_patient[condition_rand,]
  
  # the frequency of clonal type 
  temp <- apply(cl_patient_analysis,2,function(col){
    
    clonal_type <- names(col)
    
    seq <- mapply(function(cl_type, rep){
      rep(cl_type, rep)
    }, cl_type = clonal_type, rep = col) %>% 
      unlist() %>% 
      as.vector() %>% 
      sample() #shuffles the data
    
    return(seq)
  })
  
  obs_sample <- lapply(temp, function(x){
    table(x)
  })
  
  all_comp <- lapply(obs_sample, function(x){
    ks.test(hold_ct,x)
  })
  
  all_p_value <- lapply(all_comp, function(patient){
    patient$p.value 
  }) %>% unlist()
  
  all_d_value <- lapply(all_comp, function(patient){
    patient$statistic
  }) %>% unlist()
  
  
  df_res <- apply(tab_cond_patient,1,function(x){
    colnames(tab_cond_patient)[which(x != 0)]
  }) %>% as.data.frame()
  colnames(df_res)[1] <- "Condition"
  
  # df_res$significant <- sig_res
  df_res$p_value <- all_p_value
  df_res$D_value <- all_d_value
  
  hold_p_sn <- t.test(df_res$D_value[df_res$Condition == "sn_ctrl"],
                      df_res$D_value[df_res$Condition == "sn_pd"]
  )
  store_sn_p_value[[i]] <- hold_p_sn$p.value
  
  hold_p_ctx <- t.test(df_res$D_value[df_res$Condition == "ctx_ctrl"],
                       df_res$D_value[df_res$Condition == "ctx_pd"]
  )
  store_ctx_p_value[[i]] <- hold_p_ctx$p.value
  
}

df_hist <- data.frame("p_value" = c(unlist(store_sn_p_value),
                                    unlist(store_ctx_p_value)),
                      "Region" = c(rep("SN", 100),
                                   rep("Ctx", 100)))



ggplot(df_hist, aes(p_value, fill = Region)) + 
  geom_histogram(alpha = 0.9, aes(y = ..density..), position = 'identity')
