
## set working directory
setwd("~/case_study_binary_endpoints/")
## source file of functions
source("RAR_unweighted_functions_binary.r")
## load packages
library(asd)
library(multcomp)
library(multxpert)
library(cubature)
library(data.table)
library(doParallel)
library(ggplot2)

n.cluster = 39 ## number of clusters for parallel computing
one.sided.alpha = 0.025 ## one-sided significance level
n.itt.in = 10^5 ## number of simulation iterations
start.time = Sys.time()

#####################################################################################
## generate output file
NRAR.table.output.com = TRAR.table.output.com = NULL

for (scen.ind in 1){
   
     ## N-M
     n.look.in = 90
     
     ## different underlying reponse rates
     if (scen.ind==1) resp.vec.in = c(0.151, 0.282, 0.4)
     if (scen.ind==2) resp.vec.in = rep(0.1, 3)
     if (scen.ind==3) resp.vec.in = rep(0.2, 3)
     if (scen.ind==4) resp.vec.in = rep(0.3, 3)
     if (scen.ind==5) resp.vec.in = rep(0.4, 3)
     if (scen.ind==6) resp.vec.in = rep(0.5, 3)
     if (scen.ind==7) resp.vec.in = rep(0.6, 3)
     if (scen.ind==8) resp.vec.in = rep(0.7, 3)
     if (scen.ind==9) resp.vec.in = rep(0.8, 3)
     if (scen.ind==10) resp.vec.in = rep(0.9, 3)
     
     ## number of active treatment groups 
    n.endpoints = length(resp.vec.in)-1
    
    closure.name.vec = rep(NA, 2^n.endpoints-1)
    col.name.ind = 1
    for (i in 1:n.endpoints){
      combn.temp = combn(2:(n.endpoints+1), i)
      col.name.length = dim(combn.temp)[2]
      closure.name.vec[col.name.ind:(col.name.ind+col.name.length-1)] = 
        as.vector(apply(combn.temp, 2, function(x){paste0(x, collapse = "_")}))
      col.name.ind = col.name.ind+col.name.length
    }
    
    ## initial output
    n.NRAR.ind = 2
    n.TRAR.ind = 1
    NRAR.table.output = matrix(NA, nrow = n.NRAR.ind, ncol = 11+10*n.endpoints+length(closure.name.vec))
    TRAR.table.output = matrix(NA, nrow = n.TRAR.ind, ncol = 16+14*n.endpoints+length(closure.name.vec))
    
    ## burn-in period size M
    m.total.in = 90
    ## change randomization probabilities for every future enrolled subject. 
    h.in = 1
    ## total sample size
    n.RAR.in = m.total.in+h.in*n.look.in
    
    ## implement RABR
    for (NRAR.ind in 1:n.NRAR.ind){

      if (NRAR.ind==1) block.in = c(4, 4, 4)
      if (NRAR.ind==2) block.in = c(7, 7, 1)

      n.active.arm = length(resp.vec.in)-1
      n.initial.asn.m.vec = m.total.in*c(block.in[1]/sum(block.in),
                          rep((1-block.in[1]/sum(block.in))/n.active.arm, n.active.arm))
      n.target.asn.vec = n.initial.asn.m.vec + block.in/sum(block.in)*h.in*n.look.in
      
        NRAR.input.in = list("resp.vec" = resp.vec.in,
                             "n.RAR" = n.RAR.in,
                             "m.total" = m.total.in,
                             "h" = h.in,
                            "block" = block.in,
                            "n.itt" = n.itt.in,
                            "one.sided.alpha" = one.sided.alpha,
                            "n.target.asn.vec" = n.target.asn.vec
        )

        sim.NRAR.fit = sim.NRAR.func(NRAR.input.in)
        NRAR.table.output[NRAR.ind, ] = as.matrix(sim.NRAR.fit)

    }
    
    ## implement DBCD
    for (TRAR.ind in 1){
      # traditional RAR
      if (TRAR.ind==1){eta.RAR.in = 9}

      TRAR.input.in = list("resp.vec" = resp.vec.in,
                           "n.RAR" = n.RAR.in,
                           "m.total" = m.total.in,
                           "h" = h.in,
                           "eta.RAR" = eta.RAR.in,
                           "n.itt" = n.itt.in,
                           "phi.function" = "phi",
                           "one.sided.alpha" = one.sided.alpha,
                           "n.target.asn.vec" = n.target.asn.vec
      )
      sim.TRAR.fit = sim.TRAR.func(TRAR.input.in)
      TRAR.table.output[TRAR.ind, ] = as.matrix(sim.TRAR.fit)
    }

    colnames(NRAR.table.output) = colnames(sim.NRAR.fit)
    colnames(TRAR.table.output) = colnames(sim.TRAR.fit)
    
    NRAR.table.output.com = rbind(NRAR.table.output.com, NRAR.table.output)
    TRAR.table.output.com = rbind(TRAR.table.output.com, TRAR.table.output)
}

## write results to file
NRAR.table.output.com = data.frame(NRAR.table.output.com)
TRAR.table.output.com = data.frame(TRAR.table.output.com)

NRAR.table.output.com$ASN_target_pbo= n.target.asn.vec[1]
TRAR.table.output.com$ASN_target_pbo= n.target.asn.vec[1]

for (i in 1:n.active.arm){
  
  eval(parse(text = paste0("NRAR.table.output.com$ASN_target_", i, 
                           " = n.target.asn.vec[i+1]")))
  eval(parse(text = paste0("TRAR.table.output.com$ASN_target_", i,
                           " = n.target.asn.vec[i+1]")))
}


NRAR.table.output.re = NRAR.table.output.com[, c(
  "resp", "block", "m_RAR", "h_RAR","n_looks_RAR", 
  paste0("reject_adj_", 1:n.active.arm), 
  as.vector(paste0("sd_reject_adj_", 1:n.active.arm)),
  "reject_adj",  "reject_adj_al_2",
  as.vector(paste0("reject_adj_s_", 1:n.active.arm)), 
  "ASN_target_pbo", 
  paste0("ASN_target_", 1:n.active.arm), 
  paste0("ASN_s_", 1:n.active.arm), 
  paste0("ASN_s_sig_", 1:n.active.arm), 
  "ASN_pbo", paste0("ASN_", 1:n.active.arm), 
  paste0("reject_unadj_", 1:n.active.arm),
  paste0("bias_", 1:(n.active.arm+1)),
  paste0("selected_dose_", 1:n.active.arm),
  "n_total_RAR",
  "ASN_pbo_sd", "ASN_1_sd", "ASN_2_sd"
)]



TRAR.table.output.re = TRAR.table.output.com[, c(
  "resp", paste0("reject_adj_", 1:n.active.arm),
  as.vector(paste0("sd_reject_adj_", 1:n.active.arm)),
  "reject_adj","reject_adj_al_2",
  as.vector(paste0("reject_adj_s_", 1:n.active.arm)),
  "ASN_target_pbo",
  paste0("ASN_target_", 1:n.active.arm),
  "eta_RAR",
  paste0("ASN_s_", 1:n.active.arm),
  paste0("ASN_s_sig_", 1:n.active.arm),
  "ASN_pbo", paste0("ASN_", 1:n.active.arm),
  paste0("reject_", 1:n.active.arm),
  paste0("bias_", 1:(n.active.arm+1)),
  paste0("selected_dose_", 1:n.active.arm),
  "mean_lim_prob_pbo",
  paste0("mean_lim_prob_", 1:(n.active.arm)),
  "sd_lim_prob_pbo",
  paste0("sd_lim_prob_", 1:(n.active.arm)),
  paste0("prop_", 1:(n.active.arm+1)),
  paste0("sd_prop_", 1:(n.active.arm+1)),
  "func",
  "n_total_RAR","m_RAR", "h_RAR","n_looks_RAR",
  "ASN_pbo_sd", "ASN_1_sd", "ASN_2_sd"
)]

write.csv(NRAR.table.output.re, paste0("results/NRAR_bin", 
                                    "_endpoints_", n.endpoints, "_sample_size_", n.RAR.in, 
                                    "_nitt_", n.itt.in, ".csv"))
write.csv(TRAR.table.output.re, paste0("results/TRAR_bin",
                                    "_endpoints_", n.endpoints, "_sample_size_", n.RAR.in,
                                    "_nitt_", n.itt.in, ".csv"))



print(Sys.time()-start.time)



































