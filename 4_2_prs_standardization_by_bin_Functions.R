####################
# Create functions
###################

# Function to set admixture intervals
create_intervals <- function(x, use.decile=T){
  if(use.decile){
    x$AFR_intervals <- ceiling(x$AFR_ances_prop*10)*10
    x$AFR_intervals[x$AFR_ances_prop<0.50] <- 50
    x$AFR_intervals <- paste0("< ",x$AFR_intervals, "%")
    x$AFR_intervals[x$AFR_ances_prop==0] <- "0%"
    x$AFR_intervals <- factor(x$AFR_intervals,
                              levels = c("0%", "< 50%", "< 60%",
                                         "< 70%","< 80%", "< 90%", "< 100%"))
    levels(x$AFR_intervals) <- c("0%", "< 50%", "50-60%",
                                 "60-70%","70-80%", "80-90%", "> 90%")
  } else{
    x$AFR_intervals <- NA
    x$AFR_intervals[x$AFR_ances_prop>=0.70] <- "70-100%"
    x$AFR_intervals[x$AFR_ances_prop<0.70] <- "< 70%"
    x$AFR_intervals[x$Race=="W"] <- "EA"
    x$AFR_intervals <- factor(x$AFR_intervals,
                              levels = c("EA", "< 70%", "70-100%"))
  }
  return(x)
}

# Function to plot PRS by admixture intervals
plot_intervals <- function(x, prs_value="grs.wt", save_plot = FALSE){
  p1 <- x %>% 
    rename("prs_value"={{prs_value}}) %>%
    mutate(AFR_ances_intervals = factor(AFR_ances_intervals, levels=levels(AFR_ances_intervals))) %>%
    ggplot(., aes(x=prs_value, color=AFR_ances_intervals)) +
    geom_density(linewidth=1.2) + 
    viridis::scale_color_viridis(name = "% African Ancestry", discrete = TRUE) +
    labs(y = "Density",  x = paste0("PRS ")) +
    theme_classic() 
  if(save_plot){
    plot_filename <- paste0(results_dir, "/img/REGARDS_GRS-", prs_value, ".png")
    ggsave(plot_filename, p1, height = 5, width = 9)
  } else{
    return(p1)
  }
}

# Function to standardize PRS by bin
std_intervals <- function(x, prs_value="grs.wt"){
  
  # Create model to see how much PRS changes as function of % Afr in AAs
  st1_lm <- lm(paste0(prs_value, " ~ AFR_ances_prop"), data=x[x$Race=="B",])
  
  # Calculate mean based on model, and sd by bin
  x <- x %>% 
    rename("prs_col"={{prs_value}}) %>%
    add_predictions(st1_lm, var = "mean", type = NULL) %>%
    group_by(Race, AFR_ances_intervals) %>%
    mutate(sd = sd(prs_col)) %>%  
    ungroup
  
  # Use mean and sd to standardize PRS 
  x <- x %>%
    mutate(prs_std = (prs_col - mean)/sd) %>%
    rename_with(~paste0(prs_value, "_stdxbin"), prs_std) %>%
    rename_with(~prs_value, prs_col) %>%
    select(!c(mean, sd))
  return(x)
}

# Function to calculate AUC 
calc_auc <- function(test_model, x, return_AIC=FALSE){
  
  # Ensure formula is correct type
  model1 <- as.formula(test_model)
  
  # Try to run binomial model, return NA if not enough cases/top10
  assoc1 <- tryCatch( {
    glm(model1, data=x, family='binomial')
    },
    error=function(e) {
      message('Not enough observations')
      return(NA)
    }
    )
  if(any(is.na(assoc1))){return(NA)}
  
  
  # calculate ROC and AUC
  outcome <- str_split(test_model, pattern="\\s+")[[1]][1]
  predicted1 <- predict.glm(assoc1, x, type="response")
  if(length(unique(x[[outcome]]))<2){
    # skip if not enough classes
    roc1 <- auc1 <- NA
  } else{
    # calculate ROC
    roc1 <- roc(x[[outcome]], predicted1)
    # calculate AUC
    auc1 <- auc(roc1)
    # calculate AIC
    if(return_AIC){
      aic1 <- AIC(assoc1)
      stopifnot(all.equal(AIC(assoc1),
                          AIC(logLik(assoc1))))
    }
  }
  
  # create object to return AUC and AIC (if requested)
  if(return_AIC){
    model_results <- paste0(auc1, "_", aic1)
  } else{
    model_results <- auc1
  }
  
  # return AUC
  return(model_results)
}


# Function to calculate OR between quantiles
calc_or <- function(x, outcome=pheno_name, prs_variable=prs_name, 
                    prs_quantiles=seq(0, 1, 0.1), do_plot=T){
  
  if(length(unique(x[[outcome]]))>2){return(NA)}
  
  # Compare quantiles of risk of PRS
  quantile_lower <- prs_quantiles[-length(prs_quantiles)]
  quantile_upper <- prs_quantiles[-1]
  df_quant <- data.frame(prs_variable=prs_variable, 
                       quantiles=paste(quantile_lower, quantile_upper, sep="-"), 
                       cutoff1=quantile(unlist( x[[each_prs]] ), probs=quantile_lower),
                       cutoff2=quantile(unlist( x[[each_prs]] ), probs=quantile_upper),
                       or=NA, or_ci1=NA, or_ci2=NA, stringsAsFactors=F)
  df_quant$or[1] <- 1 # baseline
  
  for (q in 2:nrow(df_quant)){
    # get cutoff interval to test
    ref_quantile_cutoff <- df_quant$cutoff2[2]
    test_quantile_cutoffs <- c(df_quant$cutoff1[q],df_quant$cutoff2[q])
    
    # mutate data to flag individuals in reference and testing PRS cutoffs
    x <- x %>% 
      mutate(prs_ref = .data[[each_prs]] <= ref_quantile_cutoff) %>%
      mutate(prs_test = .data[[each_prs]] > test_quantile_cutoffs[1] &
               .data[[each_prs]] <= test_quantile_cutoffs[2]) %>%
      mutate(prs_ref = ifelse(prs_ref, 1, 0),
             prs_test = ifelse(prs_test, 1, 0)) %>% 
      mutate(across(all_of(c("prs_ref", "prs_test")), as.factor)) 
    
    #Recode variables from 0/1 to T/F
    ls_outcome <- ifelse(x[[outcome]] == 1, T, F) 
    ls_var1 <- ifelse(x[["prs_ref"]] == 1, T, F) 
    ls_var2 <- ifelse(x[["prs_test"]]==1, T, F) 

    # Create table with var1/var2 vs. CHD/not-CHD
    or_df <- matrix(c(sum(ls_var1 & !ls_outcome), sum(ls_var2 & !ls_outcome), # no disease
                      sum(ls_var1 & ls_outcome), sum(ls_var2 & ls_outcome)),  # with disease
                    nrow=2, ncol=2, byrow=TRUE)
    dimnames(or_df) <- list('Predictor'=c("PRS_ref", "PRS_test"), 
                            'Outcome'=c("neg", "pos"))
    
    # Estimate OR for quantile
    if(any(or_df==0)){next}
    or_results <- oddsratio(or_df)
    
    # Format data
    or_results <- cbind(or_results$measure, or_results$p.value)
    or_results <- round(or_results, digits=4)
    or_results <- as.data.frame(or_results)
    
    # Add to results
    df_quant$or[q] <- or_results["PRS_test", "estimate"]
    df_quant$or_ci1[q] <- or_results["PRS_test", "lower"]
    df_quant$or_ci2[q] <- or_results["PRS_test", "upper"]
  }
  
  # Plot results
  if(do_plot){
    df_quant %>%
      rownames_to_column(var="quantile") %>%
      ggplot(., aes(x=quantile, y=or)) +
      geom_linerange(aes(ymin=or_ci1, ymax=or_ci2), lwd = 1/2, 
                     position = position_dodge(width = 1/2)) +
      geom_pointrange(aes(ymin=or_ci1, ymax=or_ci2), shape = 21, lwd = 1/2, 
                      position = position_dodge(width = 1/2), fill = "white") + 
      labs(x="Quantiles", y="Odds Ratio") +
      theme_classic()
  }
  
  # Format OR results
  or_list <- paste(df_quant$or, collapse="_")
  or_ci_list <- paste(paste(df_quant$or_ci1, df_quant$or_ci2, sep="-"), collapse = "_")
  or_result <- c(OR=or_list, OR_ci=or_ci_list)
  
  # Output OR results
  return(or_result)
}

# Function to estimate full and partial R2
calc_r2 <- function(x, model_full){
  
  # Run linear regression model with PRS
  assoc_full <- lm(as.formula(model_full), data = x)
  result_full <- summary(assoc_full)
  
  # Run linear regression model WITHOUT PRS
  model2 <- gsub(paste0("[+] ", prs_name, ".*$"), "", model_full)
  assoc_base <- lm(as.formula((model2)), data = x)
  result_base <- summary(assoc_base)
  
  # Calculate adjusted partial R squared
  rsquared <- rsq.partial(assoc_full, assoc_base, adj = T)
  
  # Add results to objects
  r2_full <- round(result_full$adj.r.squared*100, digits=4) 
  r2_prs <- round(rsquared$partial.rsq*100, digits=4)
  
  # join results
  all_prs <- list(r2_full=r2_full, r2_prs=r2_prs)
  
  # return list
  return(all_prs)
}


# Function to estimate AUC, AIC and OR in a given dataset
calc_model_fit <- function(x, model1, 
                           calculate_ci=FALSE, n_bootstrap=500){
  
  # Set variables for fit estimation
  model1_outcome <- gsub(" \\~.*$", "", model1) 
  prs_predictor <- gsub("^.*\\+ ", "", model1)
  
  # Set all performance statistics empty
  auc_test_ci <- aic_test_ci <- NA
  or_q10 <- or_q25 <- or_top10 <- or_top25 <- NA
  r2_full <- r2_prs <- NA
  
  # If categorical data, calc AUC and OR
  is_categorical <- length(unique(x[[model1_outcome]])) <= 2
  if(is_categorical){
    
    # Calculate AUC and AIC by bin
    auc_aic_test <- calc_auc(model1, x=x, return_AIC=TRUE) %>%
      str_split(., "_") %>% unlist %>% as.numeric
    auc_test <- auc_aic_test[1]
    aic_test <- auc_aic_test[2]
    
    # Estimate confidence intervals for AUC and AIC
    if(calculate_ci){
      if(!is.na(auc_test)){
        # bootstrap
        foo1 <- boot(x, function(df,indices) calc_auc(model1, x=df[indices,],return_AIC=TRUE), R=n_bootstrap)
        # get 95%CI intervals
        foo1_aucs <- as.numeric(unlist(lapply(str_split(foo1$t, "_"), getFirst)))
        auc_test_ci <- paste(round(quantile(foo1_aucs,c(0.025,0.975), na.rm = T), digits=4), collapse="-")
        foo1_aics <- as.numeric(unlist(lapply(str_split(foo1$t, "_"), getLast)))
        aic_test_ci <- paste(round(quantile(foo1_aics,c(0.025,0.975), na.rm = T), digits=4), collapse="-")
      }
    }
    
    # Calculate OR for first vs. every subsequent quantile, 10% and 25%
    if(!grepl("top", prs_predictor)){
      or_q10 <- calc_or(x, outcome = model1_outcome, prs_variable=prs_predictor, 
                        prs_quantiles=seq(0, 1, 0.1))
      or_q25 <- calc_or(x, outcome = model1_outcome, prs_variable=prs_predictor, 
                        prs_quantiles=seq(0, 1, 0.25))
      # Calculate OR for first vs. all the rest, >90% and > 75%
      or_top10 <- calc_or(x, outcome = model1_outcome, prs_variable=prs_predictor, 
                          prs_quantiles=c(0, 0.9, 1))
      or_top10 <- gsub("^[1]_|^NA\\-NA[_]","", or_top10)
      or_top25 <- calc_or(x, outcome = model1_outcome, prs_variable=prs_predictor, 
                          prs_quantiles=c(0, 0.75, 1))
      or_top25 <- gsub("^[1]_|^NA\\-NA[_]","", or_top25)
    }
  }
  
  # If continuous outcome, calculate R2
  is_numeric <- is.numeric(x[[model1_outcome]])
  if(is_numeric){
    
    # Calculate R2
    all_r2 <- calc_r2(x=x, model1)
    print(each_pop)
    print(each_bin)
    print(all_r2)
    
  }
  
  # Create data frame to save results
  bin_auc_results <- data.frame(pop=each_pop, bin=each_bin,
                                prs_in_model=prs_predictor,
                                # AUCs and AICs
                                auc_value=auc_test, auc_ci=auc_test_ci,
                                aic_value=aic_test, aic_ci=aic_test_ci,
                                # ORs at 10%
                                or_top10 = or_top10[1], or_top10_ci = or_top10[2],
                                or_q10 = or_q10[1], or_q10_ci = or_q10[2],
                                # ORs at 25%
                                or_top25 = or_top25[1], or_top25_ci = or_top25[2],
                                or_q25 = or_q25[1], or_q25_ci = or_q25[2],
                                # R2
                                r2_full = r2_full,
                                r2_prs = r2_prs)


  # Format result
  bin_auc_results <- bin_auc_results %>%
    mutate(bin = factor(bin, levels=bin_list))
  ci_cols <- grep("ci", names(bin_auc_results), value = T)
  ci_cols <- grep("q", ci_cols, value = T, invert = T)
  for(each_col in ci_cols){
    bin_auc_results <- bin_auc_results %>%
      separate(col=each_col, into=paste0(each_col, c(1,2)), sep = "-|_")
  }
  row.names(bin_auc_results) <- 1:nrow(bin_auc_results)

  # Output results
  return(bin_auc_results)
}

## End
