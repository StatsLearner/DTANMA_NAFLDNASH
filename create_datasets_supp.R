# #######################################################################################
# - This is a R script to be sourced for creating a dataset.
# - For meta-analysis, this is to pick up one result per study.
# - For network meta-analysis, the study number needs to be re-numbered to be used in stan.
# - Created by Tetsuro Oda
# ########################################################################################

create_dataset_separate_MA <-
  function(.data, .testname, .fibrosisgrade) {
    
    .out <-
      
    .data %>%
      # Pick one outcome
      filter(Test == .testname,
             TargetBoundaryGrade == .fibrosisgrade) %>%
      
      # Pick up one result per study
      mutate(sens_spe = Sens + Spe) %>%
      group_by(Study_Num) %>%
      mutate(sens_spe_max = max(sens_spe)) %>%
      filter(sens_spe == sens_spe_max) %>%
      select(-sens_spe,-sens_spe_max) %>%
      ungroup()
    
    return(.out)
    
  }

create_dataset_NMA <-
  function(.data,.fibrosisgrade) {
    
    .tmp1 <-
      
      .data %>%
      # Pick one grade but include all types of diagnosis
      filter(TargetBoundaryGrade == .fibrosisgrade) 
      
    # Do re-numbering for each study for NMA
    .tmp2 <-
      .tmp1 %>%
      filter(TargetBoundaryGrade == .fibrosisgrade) %>%
      distinct(Study_Author,Study_Num) %>%
      arrange(Study_Num) %>%
      rowid_to_column() %>%
      select(-Study_Num) %>%
      rename(Study_Num = rowid)
    
    # Merge the study number
    .out <- 
      .tmp1 %>%
      select(-Study_Num) %>%
      left_join(.tmp2, by = "Study_Author") 
    
    
    return(.out)
    
  }

