
setwd("~/case_study_binary_endpoints/results/")
library(xtable)
library(ggplot2)

#####################################
## type I error
data.error = read.csv("NRAR_error_bin_endpoints_2_sample_size_180_nitt_1e+05.csv")
data.error = data.frame(data.error, stringsAsFactors = FALSE)

data.error.plot.in = data.frame("reject_unadj" = c(data.error$reject_unadj_1, data.error$reject_unadj_2,
                                                   data.error$reject_adj),
                                "group" = c(rep("pairwise low exposure", 9), rep("pairwise high exposure", 9),
                                            rep("Bonferroni at least one", 9)),
                                "rate" = rep(0.1*(1:9), 3),
                                stringsAsFactors = FALSE

)

png("error_bin.png",
    width = 1500, height = 800)
ggplot.fit = ggplot(data.error.plot.in, aes(x = rate, y = reject_unadj)) +
  geom_point(size = 8)+
  geom_line(size = 2, aes(linetype = group))+
  scale_linetype_manual(values = c("solid", "longdash", "dotted"))+
  # geom_text(aes(label=ifelse(look == max(rand_ratio_plot_data_com$look),as.character(text),'')),
  #           hjust= -0.05 ,vjust=0,size = 16) +
  scale_x_continuous(limits = c(0.1, 0.9),
                     breaks = 0.1*(1:9))+
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025),
                     limits = c(0, 0.025))+
  # geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=1, size = 2)+
  # geom_line(aes(x = look, y = target, color = Group), size = 2 ,linetype="dashed") +
  # geom_hline(aes(yintercept = target, color = group)) +
  #scale_x_continuous(breaks = c(1, 100, 200, 300, 400, sim.RAR.fit$summary.table$n_looks_RAR))+
  # labs(title = "test") +
  ylab ("Type I error rate") + xlab("Response rate") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,0,1,1),units="lines"),
        text = element_text(size=40),
        axis.text.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=30,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=35,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=35,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text = element_text(colour="black", size = 30, face = "plain"),
        legend.title = element_text(colour="black", size = 30, face = "plain"),
        legend.key.size = unit(6,"line"),
        legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 10))+
  guides(col=guide_legend(nrow=1,byrow=TRUE))

print(ggplot.fit)
dev.off()

######################################################
## power table
power.table.1.temp = read.csv("NRAR_bin_endpoints_2_sample_size_180_nitt_1e+05.csv")
power.table.2 = read.csv("TRAR_bin_endpoints_2_sample_size_180_nitt_1e+05.csv")

power.table.1.temp = data.frame(power.table.1.temp, stringsAsFactors = FALSE)
power.table.1 = rbind(power.table.1.temp[2,], power.table.1.temp[1,])

power.table.2 = data.frame(power.table.2, stringsAsFactors = FALSE)

table.power.out = data.frame("method" = c("RABR", "Fixed randomization", "DBCD"),
                             "reject_1" = paste0(sprintf('%.2f',
                                        c(power.table.1$reject_adj_s_1, power.table.2$reject_adj_s_1)*100), "%"),
                             "reject_2" = paste0(sprintf('%.2f',
                                          c(power.table.1$reject_adj_s_2, power.table.2$reject_adj_s_2)*100), "%"),
                             "reject" = paste0(sprintf('%.2f',
                                          c(power.table.1$reject_adj, power.table.2$reject_adj)*100), "%"),
                             "ASN_pbo" = paste0(sprintf('%.2f',
                                        c(power.table.1$ASN_pbo, power.table.2$ASN_pbo)), " (",
                                        sprintf('%.2f',
                                                c(power.table.1$ASN_pbo_sd, power.table.2$ASN_pbo_sd)), 
                                        ")"),
                             "ASN_1" = paste0(sprintf('%.2f',
                                        c(power.table.1$ASN_s_1, power.table.2$ASN_s_1)),
                                        " (",
                                        sprintf('%.2f',
                                                c(power.table.1$ASN_1_sd, power.table.2$ASN_1_sd)), 
                                        ")"
                                        ),
                             "ASN_2" = paste0(sprintf('%.2f',
                                      c(power.table.1$ASN_s_2, power.table.2$ASN_s_2)),
                                      " (",
                                      sprintf('%.2f',
                                              c(power.table.1$ASN_2_sd, power.table.2$ASN_2_sd)), 
                                      ")"
                                      )

                             )

print(xtable(table.power.out), include.rownames = FALSE)















