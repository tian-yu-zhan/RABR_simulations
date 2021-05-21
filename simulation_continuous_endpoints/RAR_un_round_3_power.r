
## set working directory
setwd("~/simulation_continuous_endpoints/")
## source file of functions
source("RAR_unweighted_functions_round_3.r")
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

for (scen.ind in c(1:3)){
   
     ## N-M
     n.look.in = 60

     ## different response mean
     if (scen.ind==1) resp.vec.in = c(0.43, 0.48, 0.63, 1.2)
     if (scen.ind==2) resp.vec.in = c(0.43, 0.68, 0.93, 1.2)
     if (scen.ind==3) resp.vec.in = c(0.43, 1, 1.15, 1.2)
     
     ## number of active treatment groups.  
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
    n.TRAR.ind = 3
    NRAR.table.output = matrix(NA, nrow = n.NRAR.ind, ncol = 11+10*n.endpoints+length(closure.name.vec))
    TRAR.table.output = matrix(NA, nrow = n.TRAR.ind, ncol = 17+14*n.endpoints+length(closure.name.vec))
    
    ## M for burn-in period
    m.total.in = 60
    ## change randomization probabilities for every future enrolled subject. 
    h.in = 1
    ## total number of subjects
    n.RAR.in = m.total.in+h.in*n.look.in
    
    ## implement RABR
    for (NRAR.ind in c(1:2)){

      if (NRAR.ind==1) {block.in = c(5, 5, 5, 5); RAR.func.in = 2}
      if (NRAR.ind==2) {block.in = c(9, 9, 1, 1); RAR.func.in = 2}
 
        NRAR.input.in = list("resp.vec" = resp.vec.in,
                             "n.RAR" = n.RAR.in,
                             "m.total" = m.total.in,
                             "h" = h.in,
                            "block" = block.in,
                            "n.itt" = n.itt.in,
                            "one.sided.alpha" = one.sided.alpha,
                            "RAR.func" = RAR.func.in
        )

        sim.NRAR.fit = sim.NRAR.func(NRAR.input.in)
        NRAR.table.output[NRAR.ind, ] = as.matrix(sim.NRAR.fit)

    }
    
    ## implement DBCD
    for (TRAR.ind in 1:3){
      r0.RAR.in = 5/20 ## equal randomization in burn-in period
      
      ## parameter \lambda for DBCD
      if (TRAR.ind==1){eta.RAR.in = -2}
      if (TRAR.ind==2){eta.RAR.in = 0}
      if (TRAR.ind==3){eta.RAR.in = 2}

      TRAR.input.in = list("resp.vec" = resp.vec.in,
                           "n.RAR" = n.RAR.in,
                           "m.total" = m.total.in,
                           "h" = h.in,
                           "eta.RAR" = eta.RAR.in,
                           "r0.RAR" = r0.RAR.in,
                           "n.itt" = n.itt.in,
                           "phi.function" = "phi",
                           "n.target.asn.vec" = 
      as.numeric(sim.NRAR.fit[, as.vector(c("ASN_target_pbo",paste0("ASN_target_", 1:n.endpoints)))]),
                           "one.sided.alpha" = one.sided.alpha
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

NRAR.table.output.re = NRAR.table.output.com[, c(
  "resp", "block", "m_RAR", "h_RAR","n_looks_RAR", 
  paste0("reject_adj_", 1:n.endpoints), 
  as.vector(paste0("sd_reject_adj_", 1:n.endpoints)),
  "reject_adj",  "reject_adj_al_2",
  as.vector(paste0("reject_adj_s_", 1:n.endpoints)), 
  "ASN_target_pbo", 
  paste0("ASN_target_", 1:n.endpoints), 
  paste0("ASN_s_", 1:n.endpoints), 
  paste0("ASN_s_sig_", 1:n.endpoints), 
  "ASN_pbo", paste0("ASN_", 1:n.endpoints), 
  paste0("reject_unadj_", 1:n.endpoints),
  paste0("bias_", 1:(n.endpoints+1)),
  paste0("selected_dose_", 1:n.endpoints),
  "n_total_RAR"
)]



TRAR.table.output.re = TRAR.table.output.com[, c(
  "resp", paste0("reject_adj_", 1:n.endpoints), 
  as.vector(paste0("sd_reject_adj_", 1:n.endpoints)), 
  "reject_adj","reject_adj_al_2",
  as.vector(paste0("reject_adj_s_", 1:n.endpoints)), 
  "ASN_target_pbo", 
  paste0("ASN_target_", 1:n.endpoints), 
  "eta_RAR",
  paste0("ASN_s_", 1:n.endpoints), 
  paste0("ASN_s_sig_", 1:n.endpoints), 
  "ASN_pbo", paste0("ASN_", 1:n.endpoints),
  "r0_RAR", 
  paste0("reject_", 1:n.endpoints),
  paste0("bias_", 1:(n.endpoints+1)),
  paste0("selected_dose_", 1:n.endpoints),
  "mean_lim_prob_pbo",
  paste0("mean_lim_prob_", 1:(n.endpoints)),
  "sd_lim_prob_pbo",
  paste0("sd_lim_prob_", 1:(n.endpoints)),
  paste0("prop_", 1:(n.endpoints+1)),
  paste0("sd_prop_", 1:(n.endpoints+1)),
  "func",
  "n_total_RAR","m_RAR", "h_RAR","n_looks_RAR"
)]

write.csv(NRAR.table.output.re, paste0("results/NRAR", 
                                    "_endpoints_", n.endpoints, "_sample_size_", n.RAR.in, 
                                    "_nitt_", n.itt.in, ".csv"))
write.csv(TRAR.table.output.re, paste0("results/TRAR", 
                                    "_endpoints_", n.endpoints, "_sample_size_", n.RAR.in, 
                                    "_nitt_", n.itt.in, ".csv"))



print(Sys.time()-start.time)



































