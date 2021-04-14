
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
























