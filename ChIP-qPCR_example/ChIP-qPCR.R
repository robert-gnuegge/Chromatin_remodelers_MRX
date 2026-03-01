###############################################################################
# Functions to analyze ChIP-qPCR data

# necessary libraries
require("propagate")


#--------------------------------------------------------------------------
Calc_percent_input <- function(Cq_Input, Cq_IP, Input_dil, E){
  percent_input <- E^(Cq_Input - log(1 / Input_dil, E) - Cq_IP) * 100
  return(percent_input)
}

#--------------------------------------------------------------------------
Calc_fold_over_ctrl <- function(Cq_Input_target, Cq_IP_target, Cq_Input_ctrl, Cq_IP_ctrl, Input_dil, E_target, E_ctrl){
  fold_over_ctrl <- E_target^(Cq_Input_target - log(1 / Input_dil, E_target) - Cq_IP_target) / E_ctrl^(Cq_Input_ctrl - log(1 / Input_dil, E_ctrl) - Cq_IP_ctrl)
  return(fold_over_ctrl)
}

#--------------------------------------------------------------------------
Calc_percent_input_err_propagated <- function(Cq_Input_mean_and_sd, Cq_IP_mean_and_sd, Input_dil, E){
  
  # create evaluation expression
  eval_expr <- as.expression(substitute(expr = E^(Cq_Input - log(1 / Input_dil, E) - Cq_IP) * 100,
                                        env = list(Input_dil = Input_dil, E = E)))
  # calc moments with error propagation
  propagate_res <- propagate(expr = eval_expr,
                             data = cbind(Cq_Input = as.numeric(Cq_Input_mean_and_sd),
                                          Cq_IP = as.numeric(Cq_IP_mean_and_sd)),
                             do.sim = FALSE)
  
  # construct result df
  df <- data.frame(mean = propagate_res$prop["Mean.1"], sd = propagate_res$prop["sd.1"], 
                   conf.2.5 = propagate_res$prop["2.5%"], conf.97.5 = propagate_res$prop["97.5%"],
                   row.names = NULL)

  return(df)

}

#--------------------------------------------------------------------------
Calc_fold_over_ctrl_err_propagated <- function(Cq_Input_target_mean_and_sd, Cq_IP_target_mean_and_sd, Cq_Input_ctrl_mean_and_sd, Cq_IP_ctrl_mean_and_sd, Input_dil, E_target, E_ctrl){
  
  # create evaluation expression
  eval_expr <- as.expression(substitute(expr = E_target^(Cq_Input_target - log(1 / Input_dil, E_target) - Cq_IP_target) / E_ctrl^(Cq_Input_ctrl - log(1 / Input_dil, E_ctrl) - Cq_IP_ctrl),
                                        env = list(Input_dil = Input_dil, E_target = E_target, E_ctrl = E_ctrl)))
  
  # calc moments with error propagation
  propagate_res <- propagate(expr = eval_expr,
                             data = cbind(Cq_Input_target = as.numeric(Cq_Input_target_mean_and_sd),
                                          Cq_IP_target = as.numeric(Cq_IP_target_mean_and_sd),
                                          Cq_Input_ctrl = as.numeric(Cq_Input_ctrl_mean_and_sd),
                                          Cq_IP_ctrl = as.numeric(Cq_IP_ctrl_mean_and_sd)),
                             do.sim = FALSE)

  # construct result df
  df <- data.frame(mean = propagate_res$prop["Mean.1"], sd = propagate_res$prop["sd.1"],
                   conf.2.5 = propagate_res$prop["2.5%"], conf.97.5 = propagate_res$prop["97.5%"],
                   row.names = NULL)
  
  return(df)
  
}