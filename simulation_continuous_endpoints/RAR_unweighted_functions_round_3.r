
plot.ratio.func = function(rand_ratio_plot_data, plot.name){
  png(paste0("results/plot_data/", plot.name, ".png"), 
      width = 2400, height = 1600)
  ggplot.fit = ggplot(rand_ratio_plot_data, aes(x = look, y = ratio, color = group)) +
    geom_point(size = 8, aes(color=group))+
    geom_text(aes(label=ifelse(look == max(rand_ratio_plot_data$look),as.character(text),''),
                  vjust= vadj),
              hjust = -0.3,size = 15) + 
    # geom_point()+
    # geom_line()+
    # ylim(0, 0.6)+
    # xlim(1, max(rand_ratio_plot_data$look)+10)+
    scale_x_continuous(limits = c(1, max(rand_ratio_plot_data$look)+10),
                       breaks = c(1, 20, 40, 60))+
    scale_y_continuous(breaks = c(0, round(unique(rand_ratio_plot_data$target), 2), 0.6), 
                       limits = c(0, 0.6))+
    # geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=1, size = 2)+
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

RAR.p.unadj.func = function(current.sample.list, n.active.arm, one.sided.alpha){
  pbo.vec = current.sample.list[[1]]
  dec.vec.out = rep(NA, n.active.arm)
  for (i in 1:n.active.arm){
    trt.vec = current.sample.list[[1+i]]  
    
    ## known var
    # z.stats = (mean(trt.vec) - mean(pbo.vec))/sqrt(1/length(trt.vec)+1/length(pbo.vec))
    # p.out = 1-pnorm(z.stats)
    
    ## unknown var, t-test
    p.out = t.test(trt.vec, pbo.vec, alternative = "greater")$p.value
    
    dec.vec.out[i] = p.out
  }
  return(dec.vec.out)
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

RAR.dunnett.func = function(current.sample.func, n.arm.func){
  # set.seed(1)
  
  ## Dunnett test 
  data.dunnett.stage.1 = data.frame("x" = unlist(current.sample.func),
                                    "group" = as.character(
                                      unlist(sapply(1:n.arm.func, function(x){rep(x, length(unlist(current.sample.func[[x]])))}))     
                                    ))
  
  dunnett.stage.1.anova = aov(x ~ group, data.dunnett.stage.1)
  dunnett.stage.1.fit = glht(dunnett.stage.1.anova, linfct = mcp(group = "Dunnett"),
                             alternative = "greater")
  dunnett.stage.1.summary = summary(dunnett.stage.1.fit,
                                    test=(adjusted(type = "free")))
  
  p.value.dunnett.vec = dunnett.stage.1.summary$test$pvalues
  p.value.dunnett.vec = p.value.dunnett.vec + runif(n.arm.func-1, 0, 1)*10^(-6)
  return(p.value.dunnett.vec)
}


























