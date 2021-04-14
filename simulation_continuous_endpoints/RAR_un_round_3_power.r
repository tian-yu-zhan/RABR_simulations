
setwd("~/simulation_continuous_endpoints/")
source("RAR_unweighted_functions_round_3.r")

library(asd)
library(multcomp)
library(multxpert)
library(cubature)
library(data.table)
library(doParallel)
library(ggplot2)

n.cluster = 4
one.sided.alpha = 0.025
n.itt.in = 10^5
start.time = Sys.time()

#####################################################################################
## function for parallel computing

## new RAR method
sim.NRAR.func = function(NRAR.input){
  resp.vec = NRAR.input$resp.vec
  m.total = NRAR.input$m.total
  n.arm = length(resp.vec)
  n.active.arm = n.arm -1
  h = NRAR.input$h  
  n.total = NRAR.input$n.RAR
  n.look = floor((n.total-m.total)/h)
  h.last = h + n.total -m.total-h*n.look
  block = NRAR.input$block
  n.itt = NRAR.input$n.itt
  one.sided.alpha = NRAR.input$one.sided.alpha
  RAR.func = NRAR.input$RAR.func
  
  n.initial.asn.m.vec = m.total.in/length(resp.vec.in)
  n.target.asn.vec = n.initial.asn.m.vec + block/sum(block)*h*n.look
  
  closure.name.vec = rep(NA, 2^n.active.arm-1)
  col.name.ind = 1
  for (i in 1:n.active.arm){
    combn.temp = combn(2:(n.active.arm+1), i)
    col.name.length = dim(combn.temp)[2]
    closure.name.vec[col.name.ind:(col.name.ind+col.name.length-1)] = 
      as.vector(apply(combn.temp, 2, function(x){paste0(x, collapse = "_")}))
    col.name.ind = col.name.ind+col.name.length
  }
  
  p.mat.RAR = matrix(0, nrow = n.itt, ncol = 4+8*n.active.arm+
                       length(closure.name.vec))
  colnames(p.mat.RAR) = c(
    as.vector(paste0("reject_unadj_", 1:n.active.arm)),
    as.vector(paste0("reject_adj_", 1:n.active.arm)),
    as.vector(paste0("reject_adj_s_", 1:n.active.arm)),
    "reject_adj", "reject_adj_al_2", 
    as.vector(paste0("reject_component_", closure.name.vec)),
    as.vector(paste0("bias_", 1:n.arm)),
    "ASN_pbo",
    as.vector(paste0("ASN_", 1:n.active.arm)),
    as.vector(paste0("ASN_s_", 1:n.active.arm)),
    as.vector(paste0("ASN_s_sig_", 1:n.active.arm)),
    as.vector(paste0("selected_dose_", 1:n.active.arm))
   )
  p.mat.RAR = data.frame(p.mat.RAR)

  cl <- makeCluster(n.cluster)
  registerDoParallel(cl)
  pred = foreach(itt = 1:n.itt) %dopar% {
    
    source("RAR_unweighted_functions_round_3.r")
    library(multcomp)
    if ((itt%%100)==0) print(paste("itt:", itt, "/", n.itt))
    set.seed(itt)
    
    #######################################################################################
    ## RAR with dunnett test
    rand.initial.m.vec = rep(1/n.arm, n.arm)

    n.initial.m.vec = round(rand.initial.m.vec*m.total)
    
    current.sample.list = lapply(1:n.arm, function(x){rnorm(n.initial.m.vec[x], resp.vec[x], 1)})
    
    rand.size.mat = rand.ratio.mat = matrix(NA, nrow = n.arm, ncol = n.look)
    
    h.current = h
    
    for (i in 1:n.look){
      
      if (RAR.func==1){
        stand.mean.vec = sapply(2:n.arm, function(x){(mean(current.sample.list[[x]]) -
                                                        mean(current.sample.list[[1]]))*
            sqrt(length(current.sample.list[[x]]))})
      } else if (RAR.func==2){
        stand.mean.vec = sapply(2:n.arm, function(x){(mean(current.sample.list[[x]]))*
            sqrt(length(current.sample.list[[x]]))/sd(current.sample.list[[x]])})
      } else if (RAR.func==3){
        stand.mean.vec = sapply(2:n.arm, function(x){(mean(current.sample.list[[x]]))})
      } else if (RAR.func==4){
        stand.mean.vec = sapply(2:n.arm, function(x){(mean(current.sample.list[[x]]) -
                                                        mean(current.sample.list[[1]]))/
            sqrt(1/length(current.sample.list[[x]])+1/length(current.sample.list[[1]]))})
      }
      
      ## break tie
      if (!(length(unique(stand.mean.vec))==n.active.arm)){
        stand.mean.vec = stand.mean.vec + (n.active.arm:1)*0.00001
      }
      
      RAR.ratio.unadj = c(block[1], block[-1][match(stand.mean.vec,sort(stand.mean.vec, decreasing = TRUE))])
      new.sample.n = as.vector(rmultinom(1, h.current, RAR.ratio.unadj/sum(RAR.ratio.unadj)))
      
      new.sample.list = lapply(1:n.arm, function(x){rnorm(new.sample.n[x], resp.vec[x], 1)})
      
      current.sample.list = mapply(c, current.sample.list, new.sample.list, SIMPLIFY=FALSE)
      
      rand.ratio.mat[, i] = sapply(current.sample.list, length)/
        sum(sapply(current.sample.list, length))
      rand.size.mat[, i] = sapply(current.sample.list, length)
    }
    
    dunnett.fit = RAR.dunnett.func(current.sample.list, n.arm)
    RAR.dec.adj = as.numeric(dunnett.fit<=one.sided.alpha)
    
    RAR.selected.arm = which.min(dunnett.fit)
    ## ASN 
    ASN_pbo = length(unlist(current.sample.list[[1]]))
    ASN_trt_vec = sapply(current.sample.list, length)[-1]
    
    ASN_trt_sel_vec = rep(NA, length(ASN_trt_vec))
    if (RAR.dec.adj[RAR.selected.arm]==1){
      ASN_trt_sel_vec[RAR.selected.arm] = ASN_trt_vec[RAR.selected.arm]
    }
    
    RAR.p.unadj = RAR.p.unadj.func(current.sample.list, n.active.arm, one.sided.alpha)
    RAR.dec.unadj = as.numeric(RAR.p.unadj<=one.sided.alpha)
    
    RAR.dec.adj.select = rep(0, n.active.arm)
    RAR.dec.adj.select[RAR.selected.arm] = RAR.dec.adj[RAR.selected.arm]

    ASN_s_vec = ASN_trt_vec[order(dunnett.fit, decreasing = FALSE)]
    ASN_s_sig_vec = ASN_s_vec
    ASN_s_sig_vec[RAR.dec.adj==0] = NA
 
    rand.ratio.mat.order = rand.ratio.mat[c(1, 1+order(dunnett.fit, decreasing = FALSE)), ]
    rand.size.mat.order = rand.size.mat[c(1, 1+order(dunnett.fit, decreasing = FALSE)), ]
    
    clustlist = list("reject_unadj_vec" = RAR.dec.unadj,
                     "reject_adj_vec" = RAR.dec.adj,
                     "reject_adj_s_vec" = RAR.dec.adj.select, 
                     # "reject_closure_component" = closure.func.fit$closure.dec.vec,
                     "ASN_pbo" = ASN_pbo,
                     "ASN_trt_vec" = ASN_trt_sel_vec,
                     "ASN_s_vec" = ASN_s_vec,
                     "ASN_s_sig_vec" = ASN_s_sig_vec,
                     "RAR.selected.arm" = RAR.selected.arm,
                     "rand_ratio_mat" = rand.ratio.mat.order, 
                     "rand_size_mat" = rand.size.mat.order, 
                     "bias" = (sapply(1:n.arm, function(x){
                       mean(current.sample.list[[x]])})-resp.vec)[
                         c(1, order(dunnett.fit, decreasing = FALSE)+1)
                         ]
    )
    return(clustlist)
    
  }
  
  stopCluster(cl)
  
  rand_ratio_array = rand_size_array = array(NA, dim = c(n.arm, n.look, n.itt))
  for (itt in 1:n.itt){
    
    pred.temp = pred[[itt]]
    rand_ratio_array[,,itt]  = pred.temp$rand_ratio_mat
    rand_size_array[,,itt]  = pred.temp$rand_size_mat
    
    p.mat.RAR[itt, as.vector(paste0("reject_unadj_", 1:n.active.arm))] = 
      pred.temp$reject_unadj_vec
    p.mat.RAR[itt, as.vector(paste0("reject_adj_", 1:n.active.arm))] = 
      pred.temp$reject_adj_vec
    p.mat.RAR[itt, as.vector(paste0("reject_adj_s_", 1:n.active.arm))] = 
      pred.temp$reject_adj_s_vec
    
    p.mat.RAR[itt, "reject_adj"] = as.numeric(sum(pred.temp$reject_adj_vec)>0)
    p.mat.RAR[itt, "reject_adj_al_2"] = as.numeric(sum(pred.temp$reject_adj_vec)>1)
       
    p.mat.RAR[itt, as.vector(paste0("bias_", 1:n.arm))] = pred.temp$bias
    
    p.mat.RAR[itt, "ASN_pbo"] = pred.temp$ASN_pbo
    p.mat.RAR[itt, as.vector(paste0("ASN_", 1:n.active.arm))] = pred.temp$ASN_trt_vec
    
    p.mat.RAR[itt, as.vector(paste0("ASN_s_", 1:n.active.arm))] = pred.temp$ASN_s_vec
    p.mat.RAR[itt, as.vector(paste0("ASN_s_sig_", 1:n.active.arm))] = pred.temp$ASN_s_sig_vec
 
    
    p.mat.RAR[itt, as.vector(paste0("selected_dose_", 1:n.active.arm))] = as.vector(rep(0, n.active.arm))
    p.mat.RAR[itt, as.vector(paste0("selected_dose_", pred.temp$RAR.selected.arm))] = 1
    
  }
  
  ## write randomization ratio summary file
  rand_ratio_mean_table = apply(rand_ratio_array, 1:2, mean)
  rand_ratio_sd_table = apply(rand_ratio_array, 1:2, sd)
  
  rand_size_mean_table = apply(rand_size_array, 1:2, mean)
  rand_size_sd_table = apply(rand_size_array, 1:2, sd)
  
  write.csv(rand_ratio_mean_table,
            paste0("results/plot_data/NRAR_mean_rand_ratio_", paste(resp.vec,collapse = ", "),
                   ".csv"))
  
  write.csv(rand_ratio_sd_table,
            paste0("results/plot_data/NRAR_sd_rand_ratio_", paste(resp.vec,collapse = ", "),
                   ".csv"))
  
  write.csv(rand_size_mean_table,
            paste0("results/plot_data/NRAR_mean_rand_size_", paste(resp.vec,collapse = ", "),
                   ".csv"))
  
  write.csv(rand_size_sd_table,
            paste0("results/plot_data/NRAR_sd_rand_size_", paste(resp.vec,collapse = ", "),
                   ".csv"))
  
  rand_ratio_plot_data_in = data.frame("ratio" = as.vector(t(rand_ratio_mean_table)),
                                    "sd" = as.vector(t(rand_ratio_sd_table)),
                                    "group" = rep(c("Placebo", paste0("S_", 1:n.active.arm)), 
                                                  each = n.look),
                                    "look" = rep(1:n.look, n.arm),
                                    "target" = rep(n.target.asn.vec/sum(n.target.asn.vec),
                                                   each = n.look),
                                    "text" = paste0(round(as.vector(t(rand_ratio_mean_table)), 3),
                                                    " (", 
                                                    round(as.vector(t(rand_ratio_sd_table)), 3), ")"),
                                    "vadj" = rep(c(-0.4, 1.2), each = n.look)
                                    )
  
  plot.ratio.func(rand_ratio_plot_data_in, 
                  paste0("NRAR_", paste(resp.vec,collapse = ", ")))
  
 
  
  ##########################################################
  ## put summary to table
  summary.table = apply(p.mat.RAR, 2, function(x){mean(x,na.rm = TRUE)})
  summary.table = data.frame(t(summary.table))
  
  summary.table[, as.vector(paste0("sd_reject_adj_", 1:n.active.arm))] =
    apply(p.mat.RAR[, c(as.vector(paste0("reject_adj_", 1:n.active.arm)))], 2, sd)
  
  summary.table$n_total_RAR = n.total
  summary.table$m_RAR = m.total
  summary.table$h_RAR = h
  summary.table$n_looks_RAR = n.look
  summary.table$block = paste(NRAR.input$block,collapse = ", ")
  summary.table$resp = paste(NRAR.input$resp.vec,collapse = ", ")
  
  summary.table[, as.vector(c("ASN_target_pbo",paste0("ASN_target_", 1:n.active.arm)))] = 
    n.target.asn.vec
  
  return(summary.table)
}

sim.TRAR.func = function(TRAR.input){
  resp.vec = TRAR.input$resp.vec
  m.total = TRAR.input$m.total
  n.arm = length(resp.vec)
  n.active.arm = n.arm -1
  h = TRAR.input$h  
  n.total = TRAR.input$n.RAR
  n.look = floor((n.total-m.total)/h)
  h.last = h + n.total -m.total-h*n.look
  r0.RAR = TRAR.input$r0.RAR
  n.itt = TRAR.input$n.itt
  eta.RAR= TRAR.input$eta.RAR
  phi.function = TRAR.input$phi.function
  one.sided.alpha = TRAR.input$one.sided.alpha
  n.target.asn.vec = TRAR.input$n.target.asn.vec
  
  closure.name.vec = rep(NA, 2^n.active.arm-1)
  col.name.ind = 1
  for (i in 1:n.active.arm){
    combn.temp = combn(2:(n.active.arm+1), i)
    col.name.length = dim(combn.temp)[2]
    closure.name.vec[col.name.ind:(col.name.ind+col.name.length-1)] = 
      as.vector(apply(combn.temp, 2, function(x){paste0(x, collapse = "_")}))
    col.name.ind = col.name.ind+col.name.length
  }
  
  p.mat.RAR = matrix(0, nrow = n.itt, ncol = 6+n.active.arm*10+length(closure.name.vec))
  colnames(p.mat.RAR) = c(
    as.vector(paste0("reject_", 1:n.active.arm)),
    as.vector(paste0("reject_adj_", 1:n.active.arm)),
    as.vector(paste0("reject_adj_s_", 1:n.active.arm)), 
    "reject_adj", "reject_adj_al_2", 
    as.vector(paste0("reject_component_", closure.name.vec)),
    as.vector(paste0("bias_", 1:n.arm)), 
    as.vector(paste0("prop_", 1:n.arm)), 
    "ASN_pbo",
    as.vector(paste0("ASN_", 1:n.active.arm)),
    as.vector(paste0("ASN_s_", 1:n.active.arm)),
    as.vector(paste0("ASN_s_sig_", 1:n.active.arm)),
    "mean_lim_prob_pbo",
    as.vector(paste0("mean_lim_prob_", 1:n.active.arm)),
    as.vector(paste0("selected_dose_", 1:n.active.arm)))
  p.mat.RAR = data.frame(p.mat.RAR)
  

  cl <- makeCluster(n.cluster)
  registerDoParallel(cl)
  pred = foreach(itt = 1:n.itt) %dopar% {
    
    source("RAR_unweighted_functions_round_3.r")
    library(multcomp)
    if ((itt%%100)==0) print(paste("itt:", itt, "/", n.itt))
    set.seed(itt)
    rand.ratio.mat = rand.size.mat = matrix(NA, nrow = n.arm, ncol = n.look)
 
    rand.initial.m.vec = c(r0.RAR, rep((1-r0.RAR)/n.active.arm, n.active.arm))
    
    n.initial.m.vec = round(rand.initial.m.vec*m.total)
    
    current.sample.list = lapply(1:n.arm, 
              function(x){rnorm(n.initial.m.vec[x], resp.vec[x], 1)})

    h.current = h
    
    for (i in 1:n.look){
      
      tempx.vec = sapply(current.sample.list, length)

      tempy.vec = sapply(current.sample.list, 
                  function(t){sqrt(pnorm((mean(t)-eta.RAR)/sd(t)))/sd(t)})
      
      tempx = tempx.vec / sum(tempx.vec)
      tempy = tempy.vec / sum(tempy.vec)
      
      stand.mean.vec = tempy*(tempy/tempx)^2
      
      rand.factor.vec = stand.mean.vec
      rand.ratio.vec = rand.factor.vec/sum(rand.factor.vec)
 
      if (i==n.look) h.current = h.last
      new.sample.n = as.vector(rmultinom(1, size = h.current, prob = rand.ratio.vec))
 
      new.sample.list = lapply(1:n.arm, function(x){rnorm(new.sample.n[x], resp.vec[x], 1)})
      
      current.sample.list = mapply(c, current.sample.list, new.sample.list, SIMPLIFY=FALSE)
      
      rand.ratio.mat[, i] = sapply(current.sample.list, length)/
        sum(sapply(current.sample.list, length))
      rand.size.mat[, i] = sapply(current.sample.list, length)
    }

    dunnett.fit = RAR.dunnett.func(current.sample.list, n.arm)
    RAR.dec.adj = as.numeric(dunnett.fit<=one.sided.alpha)
    
    print(dunnett.fit)
    RAR.selected.arm = which.min(dunnett.fit)
    ## ASN 
    ASN_pbo = length(unlist(current.sample.list[[1]]))
    ASN_trt_vec = sapply(current.sample.list, length)[-1]
  
    ASN_trt_sel_vec = rep(NA, length(ASN_trt_vec))
    if (RAR.dec.adj[RAR.selected.arm]==1){
      ASN_trt_sel_vec[RAR.selected.arm] = ASN_trt_vec[RAR.selected.arm]
    }
    
    RAR.p.unadj = RAR.p.unadj.func(current.sample.list, n.active.arm, one.sided.alpha)
    RAR.dec.unadj = as.numeric(RAR.p.unadj<=one.sided.alpha)
    
    RAR.dec.adj.select = rep(0, n.active.arm)
    RAR.dec.adj.select[RAR.selected.arm] = RAR.dec.adj[RAR.selected.arm]
  
    ASN_s_vec = ASN_trt_vec[order(dunnett.fit, decreasing = FALSE)]
    ASN_s_sig_vec = ASN_s_vec
    ASN_s_sig_vec[RAR.dec.adj==0] = NA

    rand.ratio.mat.order = rand.ratio.mat[c(1, 1+order(dunnett.fit, decreasing = FALSE)), ]
    rand.size.mat.order = rand.size.mat[c(1, 1+order(dunnett.fit, decreasing = FALSE)), ]
    
    lim_prob = rand.ratio.vec
    
    newlist = list("ASN_pbo" = ASN_pbo,
                   "ASN_trt_vec" = ASN_trt_sel_vec,
                   "ASN_s_vec" = ASN_s_vec,
                   "ASN_s_sig_vec" = ASN_s_sig_vec,
                   "bias" = (sapply(1:n.arm, function(x){
                     mean(current.sample.list[[x]])})-resp.vec)[
                       c(1, order(dunnett.fit, decreasing = FALSE)+1)
                       ],
                     "prop" = sapply(current.sample.list, length)/n.total,
                   "RAR.dec.unadj" = RAR.dec.unadj,
                   "reject_adj_vec" = RAR.dec.adj,
                   "reject_adj_s_vec" = RAR.dec.adj.select, 
                   "RAR.selected.arm" = RAR.selected.arm,
                   "rand_ratio_mat" = rand.ratio.mat.order, 
                   "rand_size_mat" = rand.size.mat.order, 
                   "lim_prob" = lim_prob
    )
    return(newlist)
  }
  
  stopCluster(cl)
  
  rand_ratio_array = rand_size_array = array(NA, dim = c(n.arm, n.look, n.itt))
  for (itt in 1:n.itt){
    
    pred.temp = pred[[itt]]
    rand_ratio_array[,,itt]  = pred.temp$rand_ratio_mat
    rand_size_array[,,itt]  = pred.temp$rand_size_mat
    
    p.mat.RAR[itt, as.vector(paste0("reject_", 1:n.active.arm))] = pred.temp$RAR.dec.unadj
    p.mat.RAR[itt, as.vector(paste0("reject_adj_", 1:n.active.arm))] = 
      pred.temp$reject_adj_vec
    p.mat.RAR[itt, "reject_adj"] = as.numeric(sum(pred.temp$reject_adj_vec)>0)
    p.mat.RAR[itt, "reject_adj_al_2"] = as.numeric(sum(pred.temp$reject_adj_vec)>1)
    p.mat.RAR[itt, as.vector(paste0("reject_adj_s_", 1:n.active.arm))] = 
      pred.temp$reject_adj_s_vec

    p.mat.RAR[itt, as.vector(paste0("bias_", 1:n.arm))] = pred.temp$bias
    p.mat.RAR[itt, as.vector(paste0("prop_", 1:n.arm))] = pred.temp$prop
    
    p.mat.RAR[itt, "ASN_pbo"] = pred.temp$ASN_pbo
    p.mat.RAR[itt, as.vector(paste0("ASN_", 1:n.active.arm))] = pred.temp$ASN_trt_vec

    p.mat.RAR[itt, as.vector(paste0("ASN_s_", 1:n.active.arm))] = pred.temp$ASN_s_vec
    p.mat.RAR[itt, as.vector(paste0("ASN_s_sig_", 1:n.active.arm))] = pred.temp$ASN_s_sig_vec
    
    p.mat.RAR[itt, c("mean_lim_prob_pbo",as.vector(paste0("mean_lim_prob_", 1:n.active.arm)))] = pred.temp$lim_prob
    
    p.mat.RAR[itt, as.vector(paste0("selected_dose_", 1:n.active.arm))] = as.vector(rep(0, n.active.arm))
    p.mat.RAR[itt, as.vector(paste0("selected_dose_", pred.temp$RAR.selected.arm))] = 1
    
  }
  
  ## write randomization ratio summary file
  rand_ratio_mean_table = apply(rand_ratio_array, 1:2, mean)
  rand_ratio_sd_table = apply(rand_ratio_array, 1:2, sd)
  
  rand_size_mean_table = apply(rand_size_array, 1:2, mean)
  rand_size_sd_table = apply(rand_size_array, 1:2, sd)
  
  write.csv(rand_ratio_mean_table,
            paste0("results/plot_data/TRAR_mean_rand_ratio_", paste(resp.vec,collapse = ", "),
                   "_eta_", eta.RAR, 
                   ".csv"))
  
  write.csv(rand_ratio_sd_table,
            paste0("results/plot_data/TRAR_sd_rand_ratio_", paste(resp.vec,collapse = ", "),
                   "_eta_", eta.RAR,
                   ".csv"))
  
  write.csv(rand_size_mean_table,
            paste0("results/plot_data/TRAR_mean_rand_size_", paste(resp.vec,collapse = ", "),
                   "_eta_", eta.RAR, 
                   ".csv"))
  
  write.csv(rand_size_sd_table,
            paste0("results/plot_data/TRAR_sd_rand_size_", paste(resp.vec,collapse = ", "),
                   "_eta_", eta.RAR,
                   ".csv"))
  
  rand_ratio_plot_data_in = data.frame("ratio" = as.vector(t(rand_ratio_mean_table)),
                                       "sd" = as.vector(t(rand_ratio_sd_table)),
                                       "group" = rep(c("Placebo", paste0("S_", 1:n.active.arm)), 
                                                     each = n.look),
                                       "look" = rep(1:n.look, n.arm),
                                       "target" = rep(n.target.asn.vec/sum(n.target.asn.vec),
                                                      each = n.look),
                                       "text" = paste0(round(as.vector(t(rand_ratio_mean_table)), 3),
                                                       " (", 
                                          round(as.vector(t(rand_ratio_sd_table)), 3), ")"),
                                       "vadj" = rep(c(-0.4, 1.2), each = n.look)
                                       )
  
  plot.ratio.func(rand_ratio_plot_data_in, 
                  paste0("TRAR_", paste(resp.vec,collapse = ", "), "_eta_", eta.RAR))
  
  ##########################################################
  ## put summary to table
  summary.table = apply(p.mat.RAR, 2, function(x){mean(x,na.rm = TRUE)})
  summary.table = data.frame(t(summary.table))
  
  summary.table[,c("sd_lim_prob_pbo",as.vector(paste0("sd_lim_prob_", 1:n.active.arm)))] = 
    apply(p.mat.RAR[,c("mean_lim_prob_pbo",as.vector(paste0("mean_lim_prob_", 1:n.active.arm)))], 2, sd)
  
  summary.table[, as.vector(paste0("sd_prop_", 1:n.arm))] =
    apply(p.mat.RAR[, c(as.vector(paste0("prop_", 1:n.arm)))], 2, sd)
  
  summary.table[, as.vector(paste0("sd_reject_adj_", 1:n.active.arm))] =
    apply(p.mat.RAR[, c(as.vector(paste0("reject_adj_", 1:n.active.arm)))], 2, sd)
  
  summary.table$r0_RAR = r0.RAR
  summary.table$eta_RAR = eta.RAR
  
  summary.table$n_total_RAR = n.total
  summary.table$m_RAR = m.total
  summary.table$h_RAR = h
  summary.table$n_looks_RAR = n.look
  summary.table$func = phi.function
  summary.table$resp = paste(TRAR.input$resp.vec,collapse = ", ")
  
  summary.table[, as.vector(c("ASN_target_pbo",paste0("ASN_target_", 1:n.active.arm)))] = 
    n.target.asn.vec
  
  return(summary.table)
}

#####################################################################################

NRAR.table.output.com = TRAR.table.output.com = NULL

for (scen.ind in c(1:3)){
   
     n.look.in = 60

     if (scen.ind==1) resp.vec.in = c(0.43, 0.48, 0.63, 1.2)
     if (scen.ind==2) resp.vec.in = c(0.43, 0.68, 0.93, 1.2)
     if (scen.ind==3) resp.vec.in = c(0.43, 1, 1.15, 1.2)
     
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
    
    n.NRAR.ind = 2
    n.TRAR.ind = 3
    NRAR.table.output = matrix(NA, nrow = n.NRAR.ind, ncol = 11+10*n.endpoints+length(closure.name.vec))
    TRAR.table.output = matrix(NA, nrow = n.TRAR.ind, ncol = 17+14*n.endpoints+length(closure.name.vec))
    
    m.total.in = 60
    h.in = 1

    n.RAR.in = m.total.in+h.in*n.look.in
    
    for (NRAR.ind in c(1:2)){
        # new RAR
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
    
    for (TRAR.ind in 1:3){
      r0.RAR.in = 5/20
      
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



































