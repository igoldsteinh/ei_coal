# Utility functions
library(dplyr)
library(lubridate)
library(tidybayes)
set.seed(1234)



# make fixed posterior samples --------------------------------------------

make_fixed_posterior_samples <- function(posterior_gq) {
  posterior_gq_samples <- posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]", negate = T))
  
  return(posterior_gq_samples)
}


# make time varying posterior quantiles -----------------------------------
make_timevarying_posterior_quantiles <- function(posterior_gq) {
  timevarying_posterior_quantiles <-
    posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
  
  return(timevarying_posterior_quantiles)
  
}

# make prior posterior plot -----------------------------------------------

make_fixed_param_plot <- function(posterior_gq_samples, prior_gq_samples, sim_data = TRUE){
  
  if (sim_data == TRUE) {
    priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type, true_value), prior_gq_samples)
    
    
    param_plot <- priors_and_posteriors %>%
      ggplot(aes(value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
      facet_wrap(. ~ name, scales = "free_x") +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Fixed parameter prior and posteriors")
    
  } else {
    priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type), prior_gq_samples)
    
    
    param_plot <- priors_and_posteriors %>%
      ggplot(aes(value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      facet_wrap(. ~ name, scales = "free_x") +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Fixed parameter prior and posteriors")
    
  }
  
  return(param_plot)
  
}

# rt_metrics --------------------------------------------------------------
# used for calculating frequentist characteristics of inference for rt
rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}
