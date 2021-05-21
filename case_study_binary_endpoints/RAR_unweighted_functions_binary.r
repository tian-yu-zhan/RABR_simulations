
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
  n.target.asn.vec = NRAR.input$n.target.asn.vec
  
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
    
    source("RAR_unweighted_functions_binary.r")
    library(multcomp)
    if ((itt%%100)==0) print(paste("itt:", itt, "/", n.itt))
    set.seed(itt+sum(resp.vec)*n.itt)
    
    rand.initial.m.vec = rep(1/n.arm, n.arm)
    
    n.initial.m.vec = round(rand.initial.m.vec*m.total)
    
    ## sample from binomial distribution 
    current.sample.list = lapply(1:n.arm, function(x){rbinom(size = 1, prob = resp.vec[x], n = n.initial.m.vec[x])})
    
    rand.ratio.mat = matrix(NA, nrow = n.arm, ncol = n.look)
    
    h.current = h
    
    for (i in 1:n.look){
      
      stand.mean.vec.temp = pmin(0.99,pmax(0.01, sapply(current.sample.list, mean)))
      stand.n.vec.temp = sapply(current.sample.list, length)
      
      stand.mean.vec = ((stand.mean.vec.temp[2:n.arm])/
                          (sqrt(stand.mean.vec.temp[2:n.arm]*(1-stand.mean.vec.temp[2:n.arm]))))*
        sqrt(stand.n.vec.temp[2:n.arm])
      
      ## break tie
      if (!(length(unique(stand.mean.vec))==n.active.arm)){
        stand.mean.vec = stand.mean.vec + (n.active.arm:1)*0.00001
      }
      
      RAR.ratio.unadj = c(block[1], block[-1][match(stand.mean.vec,sort(stand.mean.vec, decreasing = TRUE))])
      new.sample.n = as.vector(rmultinom(1, h.current, RAR.ratio.unadj/sum(RAR.ratio.unadj)))
      
      ## binary
      new.sample.list = lapply(1:n.arm, function(x){rbinom(n = new.sample.n[x], prob = resp.vec[x], size = 1)})
      
      current.sample.list = mapply(c, current.sample.list, new.sample.list, SIMPLIFY=FALSE)
      
      rand.ratio.mat[, i] = sapply(current.sample.list, length)/
        sum(sapply(current.sample.list, length))
    }
    
    ## bonferroni for binary
    bonferroni.fit = RAR.bonferroni.func(current.sample.list, n.arm)
    RAR.dec.adj = as.numeric(bonferroni.fit$p.adj<=one.sided.alpha)
    
    ##########################################################################
    RAR.selected.arm = which.min(bonferroni.fit$p.unadj)
    ## ASN 
    ASN_pbo = length(unlist(current.sample.list[[1]]))
    ASN_trt_vec = sapply(current.sample.list, length)[-1]
    
    ASN_trt_sel_vec = rep(NA, length(ASN_trt_vec))
    if (RAR.dec.adj[RAR.selected.arm]==1){
      ASN_trt_sel_vec[RAR.selected.arm] = ASN_trt_vec[RAR.selected.arm]
    }
    
    RAR.dec.unadj = as.numeric(bonferroni.fit$p.unadj<=one.sided.alpha)
    
    RAR.dec.adj.select = rep(0, n.active.arm)
    RAR.dec.adj.select[RAR.selected.arm] = RAR.dec.adj[RAR.selected.arm]
    
    ASN_s_vec = ASN_trt_vec[order(bonferroni.fit$p.unadj, decreasing = FALSE)]
    ASN_s_sig_vec = ASN_s_vec
    ASN_s_sig_vec[RAR.dec.adj==0] = NA
    
    rand.ratio.mat.order = rand.ratio.mat[c(1, 1+order(bonferroni.fit$p.unadj, decreasing = FALSE)), ]
    
    clustlist = list("reject_unadj_vec" = RAR.dec.unadj,
                     "reject_adj_vec" = RAR.dec.adj,
                     "reject_adj_s_vec" = RAR.dec.adj.select, 
                     "ASN_pbo" = ASN_pbo,
                     "ASN_trt_vec" = ASN_trt_sel_vec,
                     "ASN_s_vec" = ASN_s_vec,
                     "ASN_s_sig_vec" = ASN_s_sig_vec,
                     "RAR.selected.arm" = RAR.selected.arm,
                     "rand_ratio_mat" = rand.ratio.mat.order, 
                     "bias" = (sapply(1:n.arm, function(x){
                       mean(current.sample.list[[x]])})-resp.vec)[
                         c(1, order(bonferroni.fit$p.unadj, decreasing = FALSE)+1)
                         ]
    )
    return(clustlist)
    
  }
  
  stopCluster(cl)
  
  rand_ratio_array = array(NA, dim = c(n.arm, n.look, n.itt))
  for (itt in 1:n.itt){
    
    pred.temp = pred[[itt]]
    rand_ratio_array[,,itt]  = pred.temp$rand_ratio_mat
    
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
  
  write.csv(rand_ratio_mean_table,
            paste0("results/plot_data/NRAR_bin_mean_rand_ratio_", paste(resp.vec,collapse = ", "),
                   ".csv"))
  
  write.csv(rand_ratio_sd_table,
            paste0("results/plot_data/NRAR_bin_sd_rand_ratio_", paste(resp.vec,collapse = ", "),
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
                                                       round(as.vector(t(rand_ratio_sd_table)), 3), ")")
  )
  
  plot.ratio.func(rand_ratio_plot_data_in, 
                  paste0("NRAR_", paste(resp.vec,collapse = ", ")))
  
  
  
  ##########################################################
  ## put summary to table
  summary.table = apply(p.mat.RAR, 2, function(x){mean(x,na.rm = TRUE)})
  summary.table = data.frame(t(summary.table))
  
  summary.table[, as.vector(paste0("sd_reject_adj_", 1:n.active.arm))] =
    apply(p.mat.RAR[, c(as.vector(paste0("reject_adj_", 1:n.active.arm)))], 2, sd)
  
  summary.table$ASN_pbo_sd = sd(p.mat.RAR$ASN_pbo)
  summary.table$ASN_1_sd = sd(p.mat.RAR$ASN_s_1)
  summary.table$ASN_2_sd = sd(p.mat.RAR$ASN_s_2)
  
  summary.table$n_total_RAR = n.total
  summary.table$m_RAR = m.total
  summary.table$h_RAR = h
  summary.table$n_looks_RAR = n.look
  summary.table$block = paste(NRAR.input$block,collapse = ", ")
  summary.table$resp = paste(NRAR.input$resp.vec,collapse = ", ")
  
  return(summary.table)
}

## DBCD method
sim.TRAR.func = function(TRAR.input){
  resp.vec = TRAR.input$resp.vec
  m.total = TRAR.input$m.total
  n.arm = length(resp.vec)
  n.active.arm = n.arm -1
  h = TRAR.input$h  
  n.total = TRAR.input$n.RAR
  n.look = floor((n.total-m.total)/h)
  h.last = h + n.total -m.total-h*n.look
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
    as.vector(paste0("selected_dose_", 1:n.active.arm))
  )
  p.mat.RAR = data.frame(p.mat.RAR)
  
  cl <- makeCluster(n.cluster)
  registerDoParallel(cl)
  pred = foreach(itt = 1:n.itt) %dopar% {
    
    source("RAR_unweighted_functions_binary.r")
    library(multcomp)
    if ((itt%%100)==0) print(paste("itt:", itt, "/", n.itt))
    set.seed(itt+sum(resp.vec)*n.itt)
    rand.ratio.mat = matrix(NA, nrow = n.arm, ncol = n.look)
    
    rand.initial.m.vec = rep(1/n.arm, n.arm)
    
    n.initial.m.vec = round(rand.initial.m.vec*m.total)
    
    ## sample from binomial distribution 
    current.sample.list = lapply(1:n.arm, function(x){rbinom(size = 1, prob = resp.vec[x], n = n.initial.m.vec[x])})
    
    h.current = h
    
    for (i in 1:n.look){
      
      tempx.vec = sapply(current.sample.list, length)
      
      tempy.vec.temp = pmin(0.99, pmax(0.01, sapply(current.sample.list, mean)))
      tempy.vec = sqrt(tempy.vec.temp*(1-tempy.vec.temp))
      
      tempx = tempx.vec / sum(tempx.vec)
      tempy = tempy.vec / sum(tempy.vec)
      
      stand.mean.vec = tempy*(tempy/tempx)^2
      
      
      rand.factor.vec = stand.mean.vec
      
      rand.ratio.vec = rand.factor.vec/sum(rand.factor.vec)
      
      if (i==n.look) h.current = h.last
      new.sample.n = as.vector(rmultinom(1, size = h.current, prob = rand.ratio.vec))
      
      ## binary
      new.sample.list = lapply(1:n.arm, function(x){rbinom(n = new.sample.n[x], prob = resp.vec[x], size = 1)})
      
      current.sample.list = mapply(c, current.sample.list, new.sample.list, SIMPLIFY=FALSE)
      
      rand.ratio.mat[, i] = sapply(current.sample.list, length)/
        sum(sapply(current.sample.list, length))
    }
    
    ## bonferroni for binary
    bonferroni.fit = RAR.bonferroni.func(current.sample.list, n.arm)
    RAR.dec.adj = as.numeric(bonferroni.fit$p.adj<=one.sided.alpha)
    
    RAR.selected.arm = which.min(bonferroni.fit$p.unadj)
    ## ASN 
    ASN_pbo = length(unlist(current.sample.list[[1]]))
    ASN_trt_vec = sapply(current.sample.list, length)[-1]
    
    ASN_trt_sel_vec = rep(NA, length(ASN_trt_vec))
    if (RAR.dec.adj[RAR.selected.arm]==1){
      ASN_trt_sel_vec[RAR.selected.arm] = ASN_trt_vec[RAR.selected.arm]
    }
    
    
    RAR.dec.unadj = as.numeric(bonferroni.fit$p.unadj<=one.sided.alpha)
    
    
    RAR.dec.adj.select = rep(0, n.active.arm)
    RAR.dec.adj.select[RAR.selected.arm] = RAR.dec.adj[RAR.selected.arm]
    
    ASN_s_vec = ASN_trt_vec[order(bonferroni.fit$p.unadj, decreasing = FALSE)]
    ASN_s_sig_vec = ASN_s_vec
    ASN_s_sig_vec[RAR.dec.adj==0] = NA
    
    # rand.ratio.mat = cbind(n.initial.m.vec/m.total, rand.ratio.mat)
    rand.ratio.mat.order = rand.ratio.mat[c(1, 1+order(bonferroni.fit$p.unadj, decreasing = FALSE)), ]
    
    lim_prob = rand.ratio.vec
    
    newlist = list("ASN_pbo" = ASN_pbo,
                   "ASN_trt_vec" = ASN_trt_sel_vec,
                   "ASN_s_vec" = ASN_s_vec,
                   "ASN_s_sig_vec" = ASN_s_sig_vec,
                   "bias" = (sapply(1:n.arm, function(x){
                     mean(current.sample.list[[x]])})-resp.vec)[
                       c(1, order(bonferroni.fit$p.unadj, decreasing = FALSE)+1)
                       ],
                   "prop" = sapply(current.sample.list, length)/n.total,
                   "RAR.dec.unadj" = RAR.dec.unadj,
                   "reject_adj_vec" = RAR.dec.adj,
                   "reject_adj_s_vec" = RAR.dec.adj.select, 
                   "RAR.selected.arm" = RAR.selected.arm,
                   "rand_ratio_mat" = rand.ratio.mat.order, 
                   "lim_prob" = lim_prob
    )
    return(newlist)
  }
  
  stopCluster(cl)
  
  rand_ratio_array = array(NA, dim = c(n.arm, n.look, n.itt))
  for (itt in 1:n.itt){
    
    pred.temp = pred[[itt]]
    rand_ratio_array[,,itt]  = pred.temp$rand_ratio_mat
    
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
  
  write.csv(rand_ratio_mean_table,
            paste0("results/plot_data/TRAR_bin_mean_rand_ratio_", paste(resp.vec,collapse = ", "),
                   "_eta_", eta.RAR, 
                   ".csv"))
  
  write.csv(rand_ratio_sd_table,
            paste0("results/plot_data/TRAR_bin_sd_rand_ratio_", paste(resp.vec,collapse = ", "),
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
                                                       round(as.vector(t(rand_ratio_sd_table)), 3), ")")
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
  
  summary.table$ASN_pbo_sd = sd(p.mat.RAR$ASN_pbo)
  summary.table$ASN_1_sd = sd(p.mat.RAR$ASN_s_1)
  summary.table$ASN_2_sd = sd(p.mat.RAR$ASN_s_2)
  
  summary.table$eta_RAR = eta.RAR
  
  summary.table$n_total_RAR = n.total
  summary.table$m_RAR = m.total
  summary.table$h_RAR = h
  summary.table$n_looks_RAR = n.look
  summary.table$func = phi.function
  summary.table$resp = paste(TRAR.input$resp.vec,collapse = ", ")
  
  return(summary.table)
}

plot.ratio.func = function(rand_ratio_plot_data, plot.name){
  png(paste0("results/plot_data/bin_", plot.name, ".png"), 
      width = 2400, height = 1600)
  ggplot.fit = ggplot(rand_ratio_plot_data, aes(x = look, y = ratio, color = group)) +
    geom_point(size = 8, aes(color=group))+
    geom_text(aes(label=ifelse(look == max(rand_ratio_plot_data$look),as.character(text),'')),
              hjust= -0.05 ,vjust=0,size = 15) + 
    # geom_point()+
    # geom_line()+
    # ylim(0, 0.6)+
    # xlim(1, max(rand_ratio_plot_data$look)+10)+
    scale_x_continuous(limits = c(1, max(rand_ratio_plot_data$look)+10),
                       breaks = c(1, 20, 40, 60))+
    scale_y_continuous(breaks = c(0, round(unique(rand_ratio_plot_data$target), 2), 0.6), 
                       limits = c(0, 0.6))+
    geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=1, size = 2)+
    geom_line(aes(x = look, y = target, color = group), size = 2 ,linetype="dashed") + 
    # geom_hline(aes(yintercept = target, color = group)) + 
    #scale_x_continuous(breaks = c(1, 100, 200, 300, 400, sim.RAR.fit$summary.table$n_looks_RAR))+
    # labs(title = "test") +
    ylab ("Proportion") + xlab("RAR look") +
    theme_bw()+
    theme(plot.background = element_rect(fill = "transparent"),
          plot.margin = unit(c(2,0,1,1),units="lines"),
          text = element_text(size=45),
          axis.text.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=45,angle=0,hjust=1,vjust=0,face="plain"),
          axis.title.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=45,angle=90,hjust=.5,vjust=.5,face="plain"),
          legend.text = element_text(colour="black", size = 37, face = "plain"),
          legend.title = element_text(colour="black", size = 42, face = "plain"),
          legend.key.size = unit(2,"line"),
          legend.position="right", plot.title = element_text(hjust = 0.5))
  
  print(ggplot.fit)
  dev.off()
}

lm.func = function(comb.final.inn,stratum.vec.inn){
  
  n.stratum.inn = length(stratum.vec.inn)
  
  comb.final.inn$group=as.factor(comb.final.inn$group)
  comb.final.inn$stratum.ind=as.factor(comb.final.inn$stratum.ind)
  
  lm.fit = lm(resp~group+stratum.ind+group:stratum.ind, data=comb.final.inn)
  
  stratum.count = 
    aggregate(comb.final.inn$resp, list(comb.final.inn$group, comb.final.inn$stratum.ind), 
              length)
  
  # stratum.count.mat = xtabs( x ~ Group.1 + Group.2, stratum.count)
  # stratum.count.vec = apply(stratum.count.mat, 2, function(x){x[1]*x[2]/(x[1]+x[2])})
  # stratum.weight = stratum.count.vec/sum(stratum.count.vec)
  
  stratum.weight = sqrt(1/stratum.vec.inn)
  stratum.weight = stratum.weight / sum(stratum.weight)
  
  lm.contrast = matrix(c(0, 1,rep(0, n.stratum.inn-1), stratum.weight[-1]), 1)
  lm.inter = glht(lm.fit, linfct=lm.contrast)
  lm.sum = summary(lm.inter)
  
  lm.pvalue = pt(lm.sum$test$tstat, lm.sum$df, lower.tail = FALSE)
  #print(lm.fit)
  #print(lm.sum)
  return(lm.pvalue)
}

RAR.bonferroni.func = function(current.sample.func, n.arm.func){
  p.unadj.func.vec = rep(NA, n.arm.func-1)
  for (i in 2:n.arm.func){
    n.vec.func = c(length(current.sample.func[[1]]), length(current.sample.func[[i]]))
    x.vec.func = c(sum(current.sample.func[[1]]), sum(current.sample.func[[i]]))
    
    prop.test.fit = prop.test(x = x.vec.func, n = n.vec.func, alternative = "less",
                              correct = FALSE)
    p.unadj.func.vec[i-1] = prop.test.fit$p.value
  }
  p.adj.func.vec = p.adjust(p.unadj.func.vec, method = "bonferroni")
  
  newlist = list("p.unadj" = p.unadj.func.vec,
                 "p.adj" = p.adj.func.vec)
  
  return(newlist)
}
























