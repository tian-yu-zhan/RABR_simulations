
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

for (scen.ind in c(1:4)){

    ## four simulation scenarios with different response mean, M, and N-M.
    if (scen.ind==1) {resp.vec.in = rep(0, 4); m.total.in = 60; n.look.in = 60}
    if (scen.ind==2) {resp.vec.in = rep(1, 4); m.total.in = 60; n.look.in = 60}
    if (scen.ind==3) {resp.vec.in = rep(0, 4); m.total.in = 20; n.look.in = 20}
    if (scen.ind==4) {resp.vec.in = rep(1, 4); m.total.in = 20; n.look.in = 20}
     
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
    n.NRAR.ind = 5
    n.TRAR.ind = 3
    NRAR.table.output = matrix(NA, nrow = n.NRAR.ind, ncol = 11+10*n.endpoints+length(closure.name.vec))
    TRAR.table.output = matrix(NA, nrow = n.TRAR.ind, ncol = 17+14*n.endpoints+length(closure.name.vec))
    
    ## change randomization probabilities for every future enrolled subject. 
    h.in = 1

    ## total number of subjects
    n.RAR.in = m.total.in+h.in*n.look.in
    
    ## proposed RABR
    for (NRAR.ind in c(1:5)){

      ## different vector R. 
      if (NRAR.ind==1) {block.in = c(8, 4, 4, 4); RAR.func.in = 2}
      if (NRAR.ind==2) {block.in = c(8, 5, 4, 3); RAR.func.in = 2}
      if (NRAR.ind==3) {block.in = c(8, 7, 4, 1); RAR.func.in = 2}
      if (NRAR.ind==4) {block.in = c(8, 5, 5, 2); RAR.func.in = 2}
      if (NRAR.ind==5) {block.in = c(9, 9, 1, 1); RAR.func.in = 2}
 
        NRAR.input.in = list("resp.vec" = resp.vec.in,
                             "n.RAR" = n.RAR.in,
                             "m.total" = m.total.in,
                             "h" = h.in,
                            "block" = block.in,
                            "n.itt" = n.itt.in,
                            "one.sided.alpha" = one.sided.alpha,
                            "RAR.func" = RAR.func.in
        )

        ## implement RABR
        sim.NRAR.fit = sim.NRAR.func(NRAR.input.in)
        NRAR.table.output[NRAR.ind, ] = as.matrix(sim.NRAR.fit)

    }

    colnames(NRAR.table.output) = colnames(sim.NRAR.fit)

    NRAR.table.output.com = rbind(NRAR.table.output.com, NRAR.table.output)
}

## write results to file
NRAR.table.output.com = data.frame(NRAR.table.output.com)

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


write.csv(NRAR.table.output.re, paste0("results/NRAR_error", 
                                    "_endpoints_", n.endpoints, "_sample_size_", n.RAR.in, 
                                    "_nitt_", n.itt.in, ".csv"))

print(Sys.time()-start.time)



































