runDiffLME_alt <- function(peak_mat,
                           gene_mat,
                           p2gs_to_test,
                           meta_data,
                           covariates,
                           random_effect,
                           cell_type,
                           cell_type_1,
                           cell_type_2,
                           meta_cells,
                           cores,
                           out_path,
                           id){
  
  if(all.equal(colnames(peak_mat),colnames(gene_mat))){
    print("Passed cell names check between peak_mat and gene_mat")
  }else{
    stop("Failed cell names check, make sure peak_mat and gene_mat have the same column order of cell names")
  }
  
  if(all.equal(rownames(meta_data),colnames(gene_mat))){
    print("Passed cell names check between meta_data and gene_mat")
  }else{
    stop("Failed cell names check, make sure the rownames of meta_data match the colnames of gene_mat")
  }
  
  if(all.equal(rownames(meta_data),colnames(peak_mat))){
    print("Passed cell names check between meta_data and peak_mat")
  }else{
    stop("Failed cell names check, make sure the rownames of meta_data match the colnames of peak_mat")
  }
  
  if(!is.null(covariates)){
    
    for(i in covariates){
      if(i %in% colnames(meta_data)){
        print("Passed covariate-meta data check")
      }else{
        stop("Failed covariate-meta data check, make sure covariate(s) are present as columns in meta_data")
      }
    }
  }
  
  if(cell_type %in% colnames(meta_data)){
    print("Passed cell type-meta data check")
  }else{
    stop("Failed cell type-meta data check, make sure celltype is present as a column in meta_data")
  }
  
  if(random_effect %in% colnames(meta_data)){
    print("Passed random effect-meta data check")
  }else{
    stop("Failed random effect-meta data check, make sure random_effect is present as a column in meta_data")
  }
  
  
  # Using code adapted from Huang et al.: 
  # (https://github.com/QinqinHuang/CAS_eQTL/blob/519ac9d3c68631e931cf93fd7616d1dbe398afc2/Response_eQTLs/1_Interaction_test_topeSNPs_permutation.R#L1-L372)
  
  # Interaction model formulas
  if(is.null(covariates)){
    interaction_formula_h0 <- as.formula(paste0("gene ~ peak + ",cell_type,paste0(" + (1|",random_effect,")")))
    print(paste0("Interaction model h0: ",paste0("gene ~ peak + ",cell_type,paste0(" + (1|",random_effect,")"))))
    interaction_formula_h1 <- as.formula(paste0("gene ~ peak + ",cell_type,paste0(" + (1|",random_effect,")"," + ","peak:",cell_type)))
    print(paste0("Interaction model h1: ",paste0("gene ~ peak + ",cell_type,paste0(" + (1|",random_effect,")"," + ","peak:",cell_type))))
  }else{
    interaction_formula_h0 <- as.formula(paste0("gene ~ peak + ",paste0(c(cell_type,covariates),collapse = " + "),paste0(" + (1|",random_effect,")")))
    print(paste0("Interaction model h0: ",paste0("gene ~ peak + ",paste0(c(cell_type,covariates),collapse = " + "),paste0(" + (1|",random_effect,")"))))
    interaction_formula_h1 <- as.formula(paste0("gene ~ peak + ",paste0(c(cell_type,covariates),collapse = " + "),paste0(" + (1|",random_effect,")"," + ","peak:",cell_type)))
    print(paste0("Interaction model h1: ",paste0("gene ~ peak + ",paste0(c(cell_type,covariates),collapse = " + "),paste0(" + (1|",random_effect,")"," + ","peak:",cell_type))))
  }
  
  # Cell type 1 model formulas
  if(is.null(covariates)){
    cell_type_1_formula_h0 <- as.formula(paste0("gene ~ ",paste0("(1|",random_effect,")")))
    print(paste0(cell_type_1," model h0: ",paste0("gene ~ ",paste0("(1|",random_effect,")"))))
    cell_type_1_formula_h1 <- as.formula(paste0("gene ~ ",paste0("(1|",random_effect,")"," + peak")))
    print(paste0(cell_type_1," model h1: ",paste0("gene ~ ",paste0("(1|",random_effect,")"," + peak"))))
  }else{
    cell_type_1_formula_h0 <- as.formula(paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")")))
    print(paste0(cell_type_1," model h0: ",paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"))))
    cell_type_1_formula_h1 <- as.formula(paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"," + peak")))
    print(paste0(cell_type_1," model h1: ",paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"," + peak"))))
  }
  
  # Cell type 2 model formulas
  if(is.null(covariates)){
    cell_type_2_formula_h0 <- as.formula(paste0("gene ~ ",paste0("(1|",random_effect,")")))
    print(paste0(cell_type_2," model h0: ",paste0("gene ~ ",paste0("(1|",random_effect,")"))))
    cell_type_2_formula_h1 <- as.formula(paste0("gene ~ ",paste0("(1|",random_effect,")"," + peak")))
    print(paste0(cell_type_2," model h1: ",paste0("gene ~ ",paste0("(1|",random_effect,")"," + peak"))))
  }else{
    cell_type_2_formula_h0 <- as.formula(paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")")))
    print(paste0(cell_type_2," model h0: ",paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"))))
    cell_type_2_formula_h1 <- as.formula(paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"," + peak")))
    print(paste0(cell_type_2," model h1: ",paste0("gene ~ ",paste0(covariates,collapse = " + "),paste0(" + (1|",random_effect,")"," + peak"))))
  }
  
  # OLS model formula (full data, cell type 1, and cell type 2)
  OLS_formula_h1 <- as.formula("gene ~ peak")
  print("OLS model h1: gene ~ peak")

  # Go through each peak-gene pair
  print("Fitting P2G models in parallel...")
  library(doParallel)
  library(foreach)
  #library(doMC)
  registerDoParallel(cores = cores)
  #registerDoMC(cores = cores)
  results <- foreach(ii = 1:nrow(p2gs_to_test), .combine = rbind, .inorder = TRUE) %dopar% {
    
    # Make dataframe for each model
    if(is.null(covariates)){
      df <- data.frame(gene = gene_mat[grep(paste0("^",p2gs_to_test[ii,2],"$"),rownames(gene_mat)),],
                       peak = peak_mat[grep(paste0("^",p2gs_to_test[ii,1],"$"),rownames(peak_mat)),])
      df <- cbind(df,meta_data[,colnames(meta_data) %in% c(cell_type,random_effect)])
      
      df1 <- df[df$cell_type == cell_type_1,]
      df1$gene <- base::scale(df1$gene)
      if(meta_cells){
        df1$peak <- base::scale(df1$peak)
      }
      df2 <- df[df$cell_type == cell_type_2,]
      df2$gene <- base::scale(df2$gene)
      if(meta_cells){
        df2$peak <- base::scale(df2$peak)
      }
      
      df[[cell_type]] <- ifelse(df[[cell_type]] == cell_type_1,1,0)
      df$gene <- base::scale(df$gene)
      if(meta_cells){
        df$peak <- base::scale(df$peak)
      }
    }else{
      df <- data.frame(gene = gene_mat[grep(paste0("^",p2gs_to_test[ii,2],"$"),rownames(gene_mat)),],
                       peak = peak_mat[grep(paste0("^",p2gs_to_test[ii,1],"$"),rownames(peak_mat)),])
      df <- cbind(df,meta_data[,colnames(meta_data) %in% c(cell_type,random_effect,covariates)])
      
      df1 <- df[df$cell_type == cell_type_1,]
      df1$gene <- base::scale(df1$gene)
      if(meta_cells){
        df1$peak <- base::scale(df1$peak)
      }
      for( j in covariates ){
        df1[[j]] <- base::scale(df1[[j]])
      }
      df2 <- df[df$cell_type == cell_type_2,]
      df2$gene <- base::scale(df2$gene)
      if(meta_cells){
        df2$peak <- base::scale(df2$peak)
      }
      for( j in covariates ){
        df2[[j]] <- base::scale(df2[[j]])
      }
      
      df[[cell_type]] <- ifelse(df[[cell_type]] == cell_type_1,1,0)
      df$gene <- base::scale(df$gene)
      if(meta_cells){
        df$peak <- base::scale(df$peak)
      }
      for( j in covariates ){
        df[[j]] <- base::scale(df[[j]])
      }
    }
    # print(paste0("Mean value scaled gene expression (full data): ",mean(df$gene)))
    # print(paste0("SD scaled gene expression (full data): ",sd(df$gene)))
    # print(paste0("Mean value scaled peak accessibility (full data): ",mean(df$peak)))
    # print(paste0("SD scaled peak accessibility (full data): ",sd(df$peak)))
    # 
    # print(paste0("Mean value scaled gene expression (cell type 1): ",mean(df1$gene)))
    # print(paste0("SD scaled gene expression (cell type 1): ",sd(df1$gene)))
    # print(paste0("Mean value scaled peak accessibility (cell type 1): ",mean(df1$peak)))
    # print(paste0("SD scaled peak accessibility (cell type 1): ",sd(df1$peak)))
    # 
    # print(paste0("Mean value scaled gene expression (cell type 2): ",mean(df2$gene)))
    # print(paste0("SD scaled gene expression (cell type 2): ",sd(df2$gene)))
    # print(paste0("Mean value scaled peak accessibility (cell type 2): ",mean(df2$peak)))
    # print(paste0("SD scaled peak accessibility (cell type 2): ",sd(df2$peak)))
    
    # Perform interaction model first ########################################## 
    
    # Return results 
    return_results_interaction <- p2gs_to_test[ii,]
    
    if(any(is.na(df))){
      # Statistics for peak peak term (Satterthwaite's method)
      return_results_interaction$peak_effect_size <- NA
      return_results_interaction$peak_se <- NA
      return_results_interaction$peak_df <- NA
      return_results_interaction$peak_tvalue <- NA
      return_results_interaction$peak_pval <- NA
      
      # ANOVA
      return_results_interaction$peak_ANOVA_Chisq <- NA
      return_results_interaction$peak_ANOVA_pval <- NA
      
      # Statistics for peak peak by condition interaction term (Satterthwaite's method)
      return_results_interaction$peak_cell_type_interaction_effect_size <- NA
      return_results_interaction$peak_cell_type_interaction_se <- NA
      return_results_interaction$peak_cell_type_interaction_df <- NA
      return_results_interaction$peak_cell_type_interaction_tvalue <- NA
      return_results_interaction$peak_cell_type_interaction_pval <- NA
      
      # ANOVA
      return_results_interaction$peak_cell_type_interaction_ANOVA_Chisq <- NA
      return_results_interaction$peak_cell_type_interaction_ANOVA_pval <- NA
      
      # Singular fit
      return_results_interaction$peak_cell_type_interaction_singular_fit <- NA
      
      # Random effect variance and standard deviation
      return_results_interaction$random_effect_variance <- NA
      return_results_interaction$random_effect_sd <- NA
      
    }else{
      # Test peak accessibility in linear mixed effects model
      interaction_model_H0 <- lmerTest::lmer(interaction_formula_h0, data = df, REML = FALSE)
      interaction_model_H1 <- lmerTest::lmer(interaction_formula_h1, data = df, REML = FALSE)
      
      # Get random effect variance
      vc <- as.data.frame(lme4::VarCorr(interaction_model_H1))
      
      # Compare two models using likelihood ratio test (anova)
      anova <- anova(interaction_model_H0, interaction_model_H1)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(interaction_model_H1,ddf = "Satterthwaite")
      
      # Statistics for peak peak term (Satterthwaite's method)
      if( length(grep("peak:cell_type",rownames(modelSummary$coefficients)))>0 ){
        
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_interaction$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_interaction$peak_se <- modelSummary$coefficients["peak",2]
        return_results_interaction$peak_df <- modelSummary$coefficients["peak",3]
        return_results_interaction$peak_tvalue <- modelSummary$coefficients["peak",4]
        return_results_interaction$peak_pval <- modelSummary$coefficients["peak",5]
        
        # ANOVA
        return_results_interaction$peak_ANOVA_Chisq <- anova$Chisq[2]
        return_results_interaction$peak_ANOVA_pval <- anova$`Pr(>Chisq)`[2]
        
        # Statistics for peak peak by condition interaction term (Satterthwaite's method)
        return_results_interaction$peak_cell_type_interaction_effect_size <- modelSummary$coefficients["peak:cell_type",1]
        return_results_interaction$peak_cell_type_interaction_se <- modelSummary$coefficients["peak:cell_type",2]
        return_results_interaction$peak_cell_type_interaction_df <- modelSummary$coefficients["peak:cell_type",3]
        return_results_interaction$peak_cell_type_interaction_tvalue <- modelSummary$coefficients["peak:cell_type",4]
        return_results_interaction$peak_cell_type_interaction_pval <- modelSummary$coefficients["peak:cell_type",5]
        
        # ANOVA
        return_results_interaction$peak_cell_type_interaction_ANOVA_Chisq <- anova$Chisq[2]
        return_results_interaction$peak_cell_type_interaction_ANOVA_pval <- anova$`Pr(>Chisq)`[2]
        
        # Singular fit
        return_results_interaction$peak_cell_type_interaction_singular_fit <- ifelse(lme4::isSingular(interaction_model_H1), TRUE, FALSE)
        
        # Random effect variance and standard deviation
        return_results_interaction$random_effect_variance <- vc[1,4]
        return_results_interaction$random_effect_sd <- vc[1,5]
        
      }else{
        
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_interaction$peak_effect_size <- NA
        return_results_interaction$peak_se <- NA
        return_results_interaction$peak_df <- NA
        return_results_interaction$peak_tvalue <- NA
        return_results_interaction$peak_pval <- NA
        
        # ANOVA
        return_results_interaction$peak_ANOVA_Chisq <- NA
        return_results_interaction$peak_ANOVA_pval <- NA
        
        # Statistics for peak peak by condition interaction term (Satterthwaite's method)
        return_results_interaction$peak_cell_type_interaction_effect_size <- NA
        return_results_interaction$peak_cell_type_interaction_se <- NA
        return_results_interaction$peak_cell_type_interaction_df <- NA
        return_results_interaction$peak_cell_type_interaction_tvalue <- NA
        return_results_interaction$peak_cell_type_interaction_pval <- NA
        
        # ANOVA
        return_results_interaction$peak_cell_type_interaction_ANOVA_Chisq <- NA
        return_results_interaction$peak_cell_type_interaction_ANOVA_pval <- NA
        
        # Singular fit
        return_results_interaction$peak_cell_type_interaction_singular_fit <- NA
        
        # Random effect variance and standard deviation
        return_results_interaction$random_effect_variance <- NA
        return_results_interaction$random_effect_sd <- NA
      }
    }
    
    # Perform cell_type_1 model ################################################
    
    # Return results 
    return_results_cell_type_1 <- p2gs_to_test[ii,]
    
    if( any(is.na(df1)) ){
      
      return_results_cell_type_1$peak_effect_size <- NA
      return_results_cell_type_1$peak_se <- NA
      return_results_cell_type_1$peak_df <- NA
      return_results_cell_type_1$peak_tvalue <- NA
      return_results_cell_type_1$peak_pval <- NA
      
      # ANOVA
      return_results_cell_type_1$peak_ANOVA_Chisq <- NA
      return_results_cell_type_1$peak_ANOVA_pval <- NA
      
      # Singular fit
      return_results_cell_type_1$peak_singular_fit <- NA
      
      # Random effect variance and standard deviation
      return_results_cell_type_1$random_effect_variance <- NA
      return_results_cell_type_1$random_effect_sd <- NA
      
    }else{
      # Test peak accessibility in linear mixed effects model
      cell_type_1_model_H0 <- lmerTest::lmer(cell_type_1_formula_h0, data = df1, REML = FALSE)
      cell_type_1_model_H1 <- lmerTest::lmer(cell_type_1_formula_h1, data = df1, REML = FALSE)
      
      # Get random effect variance
      vc <- as.data.frame(lme4::VarCorr(cell_type_1_model_H1))
      
      # Compare two models using likelihood ratio test (anova)
      anova <- anova(cell_type_1_model_H0, cell_type_1_model_H1)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(cell_type_1_model_H1,ddf = "Satterthwaite")
      
      if( length(grep("peak",rownames(modelSummary$coefficients)))>0 ){       
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_cell_type_1$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_cell_type_1$peak_se <- modelSummary$coefficients["peak",2]
        return_results_cell_type_1$peak_df <- modelSummary$coefficients["peak",3]
        return_results_cell_type_1$peak_tvalue <- modelSummary$coefficients["peak",4]
        return_results_cell_type_1$peak_pval <- modelSummary$coefficients["peak",5]
        
        # ANOVA
        return_results_cell_type_1$peak_ANOVA_Chisq <- anova$Chisq[2]
        return_results_cell_type_1$peak_ANOVA_pval <- anova$`Pr(>Chisq)`[2]
        
        # Singular fit
        return_results_cell_type_1$peak_singular_fit <- ifelse(lme4::isSingular(cell_type_1_model_H1), TRUE, FALSE)
        
        # Random effect variance and standard deviation
        return_results_cell_type_1$random_effect_variance <- vc[1,4]
        return_results_cell_type_1$random_effect_sd <- vc[1,5]
        
      }else{
        return_results_cell_type_1$peak_effect_size <- NA
        return_results_cell_type_1$peak_se <- NA
        return_results_cell_type_1$peak_df <- NA
        return_results_cell_type_1$peak_tvalue <- NA
        return_results_cell_type_1$peak_pval <- NA
        
        # ANOVA
        return_results_cell_type_1$peak_ANOVA_Chisq <- NA
        return_results_cell_type_1$peak_ANOVA_pval <- NA
        
        # Singular fit
        return_results_cell_type_1$peak_singular_fit <- NA
        
        # Random effect variance and standard deviation
        return_results_cell_type_1$random_effect_variance <- NA
        return_results_cell_type_1$random_effect_sd <- NA
        
      }
    }
    
    # Perform cell_type_2 model ################################################
    
    # Return results 
    return_results_cell_type_2 <- p2gs_to_test[ii,]
    
    if( any(is.na(df2)) ){
      
      return_results_cell_type_2$peak_effect_size <- NA
      return_results_cell_type_2$peak_se <- NA
      return_results_cell_type_2$peak_df <- NA
      return_results_cell_type_2$peak_tvalue <- NA
      return_results_cell_type_2$peak_pval <- NA
      
      # ANOVA
      return_results_cell_type_2$peak_ANOVA_Chisq <- NA
      return_results_cell_type_2$peak_ANOVA_pval <- NA
      
      # Singular fit
      return_results_cell_type_2$peak_singular_fit <- NA
      
      # Random effect variance and standard deviation
      return_results_cell_type_2$random_effect_variance <- NA
      return_results_cell_type_2$random_effect_sd <- NA
      
    }else{
      # Test peak accessibility in linear mixed effects model
      cell_type_2_model_H0 <- lmerTest::lmer(cell_type_2_formula_h0, data = df2, REML = FALSE)
      cell_type_2_model_H1 <- lmerTest::lmer(cell_type_2_formula_h1, data = df2, REML = FALSE)
      
      # Get random effect variance
      vc <- as.data.frame(lme4::VarCorr(cell_type_2_model_H1))
      
      # Compare two models using likelihood ratio test (anova)
      anova <- anova(cell_type_2_model_H0, cell_type_2_model_H1)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(cell_type_2_model_H1,ddf = "Satterthwaite")
      
      if( length(grep("peak",rownames(modelSummary$coefficients)))>0 ){       
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_cell_type_2$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_cell_type_2$peak_se <- modelSummary$coefficients["peak",2]
        return_results_cell_type_2$peak_df <- modelSummary$coefficients["peak",3]
        return_results_cell_type_2$peak_tvalue <- modelSummary$coefficients["peak",4]
        return_results_cell_type_2$peak_pval <- modelSummary$coefficients["peak",5]
        
        # ANOVA
        return_results_cell_type_2$peak_ANOVA_Chisq <- anova$Chisq[2]
        return_results_cell_type_2$peak_ANOVA_pval <- anova$`Pr(>Chisq)`[2]
        
        # Singular fit
        return_results_cell_type_2$peak_singular_fit <- ifelse(lme4::isSingular(cell_type_2_model_H1), TRUE, FALSE)
        
        # Random effect variance and standard deviation
        return_results_cell_type_2$random_effect_variance <- vc[1,4]
        return_results_cell_type_2$random_effect_sd <- vc[1,5]
        
      }else{
        return_results_cell_type_2$peak_effect_size <- NA
        return_results_cell_type_2$peak_se <- NA
        return_results_cell_type_2$peak_df <- NA
        return_results_cell_type_2$peak_tvalue <- NA
        return_results_cell_type_2$peak_pval <- NA
        
        # ANOVA
        return_results_cell_type_2$peak_ANOVA_Chisq <- NA
        return_results_cell_type_2$peak_ANOVA_pval <- NA
        
        # Singular fit
        return_results_cell_type_2$peak_singular_fit <- NA
        
        # Random effect variance and standard deviation
        return_results_cell_type_2$random_effect_variance <- NA
        return_results_cell_type_2$random_effect_sd <- NA
      }
    }
    
    # Perform OLS with full data ###############################################
    
    # Return results 
    return_results_OLS <- p2gs_to_test[ii,]
    
    if( any(is.na(df)) ){
      
      return_results_OLS$peak_effect_size <- NA
      return_results_OLS$peak_se <- NA
      return_results_OLS$peak_tvalue <- NA
      return_results_OLS$peak_pval <- NA
      
    }else{
      
      # Test peak accessibility in linear model
      OLS_model_H1 <- stats::lm(OLS_formula_h1,data=df)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(OLS_model_H1)
      
      if( length(grep("peak",rownames(modelSummary$coefficients)))>0 ){       
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_OLS$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_OLS$peak_se <- modelSummary$coefficients["peak",2]
        return_results_OLS$peak_tvalue <- modelSummary$coefficients["peak",3]
        return_results_OLS$peak_pval <- modelSummary$coefficients["peak",4]
        
      }else{
        return_results_OLS$peak_effect_size <- NA
        return_results_OLS$peak_se <- NA
        return_results_OLS$peak_tvalue <- NA
        return_results_OLS$peak_pval <- NA
      }
    }
    
    # Perform OLS cell type 1 ###############################################
    
    # Return results 
    return_results_OLS_cell_type_1 <- p2gs_to_test[ii,]
    
    if( any(is.na(df1)) ){
      
      return_results_OLS_cell_type_1$peak_effect_size <- NA
      return_results_OLS_cell_type_1$peak_se <- NA
      return_results_OLS_cell_type_1$peak_tvalue <- NA
      return_results_OLS_cell_type_1$peak_pval <- NA
      
    }else{
      
      # Test peak accessibility in linear model
      OLS_model_H1 <- stats::lm(OLS_formula_h1,data=df1)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(OLS_model_H1)
      
      if( length(grep("peak",rownames(modelSummary$coefficients)))>0 ){       
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_OLS_cell_type_1$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_OLS_cell_type_1$peak_se <- modelSummary$coefficients["peak",2]
        return_results_OLS_cell_type_1$peak_tvalue <- modelSummary$coefficients["peak",3]
        return_results_OLS_cell_type_1$peak_pval <- modelSummary$coefficients["peak",4]
      }else{
        return_results_OLS_cell_type_1$peak_effect_size <- NA
        return_results_OLS_cell_type_1$peak_se <- NA
        return_results_OLS_cell_type_1$peak_tvalue <- NA
        return_results_OLS_cell_type_1$peak_pval <- NA
      }
    }
    
    # Perform OLS cell type 2 ###############################################
    
    # Return results 
    return_results_OLS_cell_type_2 <- p2gs_to_test[ii,]
    
    if( any(is.na(df2)) ){
      
      return_results_OLS_cell_type_2$peak_effect_size <- NA
      return_results_OLS_cell_type_2$peak_se <- NA
      return_results_OLS_cell_type_2$peak_tvalue <- NA
      return_results_OLS_cell_type_2$peak_pval <- NA
      
    }else{
      
      # Test peak accessibility in linear model
      OLS_model_H1 <- stats::lm(OLS_formula_h1,data=df2)
      
      # Evaluate peak accessibility significance using Satterthwaite's method
      modelSummary <- summary(OLS_model_H1)
      
      if( length(grep("peak",rownames(modelSummary$coefficients)))>0 ){       
        # Statistics for peak peak term (Satterthwaite's method)
        return_results_OLS_cell_type_2$peak_effect_size <- modelSummary$coefficients["peak",1]
        return_results_OLS_cell_type_2$peak_se <- modelSummary$coefficients["peak",2]
        return_results_OLS_cell_type_2$peak_tvalue <- modelSummary$coefficients["peak",3]
        return_results_OLS_cell_type_2$peak_pval <- modelSummary$coefficients["peak",4]
      }else{
        return_results_OLS_cell_type_2$peak_effect_size <- NA
        return_results_OLS_cell_type_2$peak_se <- NA
        return_results_OLS_cell_type_2$peak_tvalue <- NA
        return_results_OLS_cell_type_2$peak_pval <- NA
      }
    }
    
    # Concatenate statistics from each model
    colnames(return_results_interaction) <- paste0("LME_interaction_model_",colnames(return_results_interaction))
    colnames(return_results_cell_type_1) <- paste0("LME_cell_type_1_model_",colnames(return_results_cell_type_1))
    colnames(return_results_cell_type_2) <- paste0("LME_cell_type_2_model_",colnames(return_results_cell_type_2))
    colnames(return_results_OLS) <- paste0("OLS_model_",colnames(return_results_OLS))
    colnames(return_results_OLS_cell_type_1) <- paste0("OLS_cell_type_1_model_",colnames(return_results_OLS_cell_type_1))
    colnames(return_results_OLS_cell_type_2) <- paste0("OLS_cell_type_2_model_",colnames(return_results_OLS_cell_type_2))
    
    res <- cbind(return_results_interaction,
                 return_results_cell_type_1,
                 return_results_cell_type_2,
                 return_results_OLS,
                 return_results_OLS_cell_type_1,
                 return_results_OLS_cell_type_2)
    
    # Write to out_path
    data.table::fwrite(res,
                       file = paste0(out_path,"/interactionLMM_univariateLMM_and_OLS_results-",id,".csv"),
                       append = TRUE,
                       quote = "auto",
                       sep = ",",
                       row.names = FALSE,
                       col.names = ifelse(ii == 1,TRUE,FALSE))
    
    print(paste0("Completed P2G: ",ii))
    
    return(NULL)
    
  }
  
  return(paste0(out_path,"/interactionLMM_univariateLMM_and_OLS_results-",id,".csv"))
  
}