
setwd("~/simulation_continuous_endpoints/results/")
library(xtable)
library(ggplot2)

######################################################################################
## type I error

data.NRAR.in = read.csv("NRAR_error_endpoints_3_sample_size_40_nitt_1e+05.csv")
data.NRAR.in = data.frame(data.NRAR.in,stringsAsFactors = FALSE)
data.NRAR.in$resp_num = c(0, rep("", 4), 1, rep("", 4), 0, rep(NA, 4), 1, rep(NA, 4))
data.NRAR.in$n = c(120, rep("", 9), 40, rep("", 9))

table.NRAR.in = data.NRAR.in[, c("n", "resp_num", "block",
                                 "reject_unadj_1", "reject_unadj_2", "reject_unadj_3",
                                 "reject_adj_1", "reject_adj_2", "reject_adj_3",
                                 "reject_adj"
                                 )]

table.NRAR.in$reject_unadj_1 = paste0(sprintf('%.2f',table.NRAR.in$reject_unadj_1*100), "%")
table.NRAR.in$reject_unadj_2 = paste0(sprintf('%.2f',table.NRAR.in$reject_unadj_2*100), "%")
table.NRAR.in$reject_unadj_3 = paste0(sprintf('%.2f',table.NRAR.in$reject_unadj_3*100), "%")
table.NRAR.in$reject_adj_1 = paste0(sprintf('%.2f',table.NRAR.in$reject_adj_1*100), "%")
table.NRAR.in$reject_adj_2 = paste0(sprintf('%.2f',table.NRAR.in$reject_adj_2*100), "%")
table.NRAR.in$reject_adj_3 = paste0(sprintf('%.2f',table.NRAR.in$reject_adj_3*100), "%")
table.NRAR.in$reject_adj = paste0(sprintf('%.2f',table.NRAR.in$reject_adj*100), "%")

print(xtable(table.NRAR.in,digits = 4), include.rownames = FALSE)

#######################################################################################
## power
power.com.table.final = NULL
power.NRAR.com = read.csv("NRAR_endpoints_3_sample_size_120_nitt_1e+05.csv")
power.TRAR.com = read.csv("TRAR_endpoints_3_sample_size_120_nitt_1e+05.csv")
for (i in 1:3){

  power.NRAR.temp.reverse = data.frame(power.NRAR.com[(i-1)*2+(1:2), ])
  power.TRAR.temp = data.frame(power.TRAR.com[(i-1)*3+(1:3), ])

  power.NRAR.temp = rbind(power.NRAR.temp.reverse[2, ], power.NRAR.temp.reverse[1, ])

  power.NRAR.table = data.frame("mu" = paste0("$mu_",c("A", "B", "C")[i], "$"),
                                "method" = c("RABR", "Fixed"),
                                "lambda" = 0,
                                "reject_1" =
                                  paste0(sprintf('%.2f',power.NRAR.temp$reject_adj_s_1*100), "%"),
                                "reject_2" =
                                paste0(sprintf('%.2f',power.NRAR.temp$reject_adj_s_2*100), "%"),
                                "reject_3" =
                                paste0(sprintf('%.2f',power.NRAR.temp$reject_adj_s_3*100), "%"),
                                "overall" =
                                  paste0(sprintf('%.2f',power.NRAR.temp$reject_adj*100), "%"),

                                # "ASN_pbo" = sprintf('%.2f',power.NRAR.temp$ASN_pbo),
                                # "ASN_1" = sprintf('%.2f',power.NRAR.temp$ASN_s_1),
                                # "ASN_2" = sprintf('%.2f',power.NRAR.temp$ASN_s_2),
                                # "ASN_3" = sprintf('%.2f',power.NRAR.temp$ASN_s_3),
                                stringsAsFactors = FALSE
                                )

  power.TRAR.table = data.frame("mu" = "",
                                "method" = c("DBCD", "", ""),
                                "lambda" = power.TRAR.temp$eta_RAR,
                                "reject_1" =
                                  paste0(sprintf('%.2f',power.TRAR.temp$reject_adj_s_1*100), "%"),
                                "reject_2" =
                                  paste0(sprintf('%.2f',power.TRAR.temp$reject_adj_s_2*100), "%"),
                                "reject_3" =
                                  paste0(sprintf('%.2f',power.TRAR.temp$reject_adj_s_3*100), "%"),
                                "overall" =
                                  paste0(sprintf('%.2f',power.TRAR.temp$reject_adj*100), "%"),
                                # "ASN_pbo" = sprintf('%.2f',power.TRAR.temp$ASN_pbo),
                                # "ASN_1" = sprintf('%.2f',power.TRAR.temp$ASN_s_1),
                                # "ASN_2" = sprintf('%.2f',power.TRAR.temp$ASN_s_2),
                                # "ASN_3" = sprintf('%.2f',power.TRAR.temp$ASN_s_3),
                                stringsAsFactors = FALSE
  )


  power.com.table = rbind(power.NRAR.table, power.TRAR.table)
  power.com.table$lambda[1:2] = "-"
  power.com.table$mu[2] = ""

  power.com.table.final = rbind(power.com.table.final, power.com.table)
}

# power.com.table.final[2, ] = paste0("\textbf{", as.character(power.com.table.final[2, ]),
#                                     "}")

print(xtable(power.com.table.final,digits = 4), include.rownames = FALSE)

#####################################################################################
## ASN table
power.com.table.final = NULL
power.NRAR.com = read.csv("NRAR_endpoints_3_sample_size_120_nitt_1e+05.csv")
power.TRAR.com = read.csv("TRAR_endpoints_3_sample_size_120_nitt_1e+05.csv")
for (i in 1:3){

  power.NRAR.temp.reverse = data.frame(power.NRAR.com[(i-1)*2+(1:2), ])
  power.TRAR.temp = data.frame(power.TRAR.com[(i-1)*3+(1:3), ])

  power.NRAR.temp = rbind(power.NRAR.temp.reverse[2, ], power.NRAR.temp.reverse[1, ])

  power.NRAR.table = data.frame("mu" = paste0("$mu_",c("A", "B", "C")[i], "$"),
                                "method" = c("RABR", "Fixed"),
                                "lambda" = 0,
                                "ASN_pbo" = sprintf('%.2f',power.NRAR.temp$ASN_pbo),
                                "ASN_1" = sprintf('%.2f',power.NRAR.temp$ASN_s_1),
                                "ASN_2" = sprintf('%.2f',power.NRAR.temp$ASN_s_2),
                                "ASN_3" = sprintf('%.2f',power.NRAR.temp$ASN_s_3),
                                stringsAsFactors = FALSE
  )

  power.TRAR.table = data.frame("mu" = "",
                                "method" = c("DBCD", "", ""),
                                "lambda" = power.TRAR.temp$eta_RAR,
                                "ASN_pbo" = sprintf('%.2f',power.TRAR.temp$ASN_pbo),
                                "ASN_1" = sprintf('%.2f',power.TRAR.temp$ASN_s_1),
                                "ASN_2" = sprintf('%.2f',power.TRAR.temp$ASN_s_2),
                                "ASN_3" = sprintf('%.2f',power.TRAR.temp$ASN_s_3),
                                stringsAsFactors = FALSE
  )


  power.com.table = rbind(power.NRAR.table, power.TRAR.table)
  power.com.table$lambda[1:2] = "-"
  power.com.table$mu[2] = ""

  power.com.table.final = rbind(power.com.table.final, power.com.table)
}

# power.com.table.final[2, ] = paste0("\textbf{", as.character(power.com.table.final[2, ]),
#                                     "}")

print(xtable(power.com.table.final,digits = 4), include.rownames = FALSE)

#####################################################################################
## plot for sample proportions. 

rand_ratio_plot_data_com = NULL
for (i in 1:6){
  response_name_vec = c("0.43, 0.48, 0.63, 1.2",
                        "0.43, 0.68, 0.93, 1.2",
                        "0.43, 1, 1.15, 1.2"
                        )

  response_name = response_name_vec[(i+1)%/%2]

  if ((i%%2)==1){
    method_name = "NRAR"
    method_new_name = "RABR"
    rand_ratio_mean_table = read.csv(
      paste0("plot_data/", method_name,"_mean_rand_ratio_", response_name,".csv"))[,-1]
    rand_ratio_sd_table = read.csv(
      paste0("plot_data/", method_name,"_sd_rand_ratio_", response_name,".csv"))[,-1]
  }else {
    method_name = "TRAR"
    method_new_name = "DBCD"
    rand_ratio_mean_table = read.csv(
      paste0("plot_data/", method_name,"_mean_rand_ratio_", response_name,"_eta_2.csv"))[,-1]
    rand_ratio_sd_table = read.csv(
      paste0("plot_data/", method_name,"_sd_rand_ratio_", response_name,"_eta_2.csv"))[,-1]
  }

  n.active.arm = dim(rand_ratio_mean_table)[1]-1
  n.look = dim(rand_ratio_mean_table)[2]
  n.arm = dim(rand_ratio_mean_table)[2]

  test.file = read.csv("NRAR_endpoints_3_sample_size_120_nitt_1e+05.csv")[2,]
  n.target.asn.vec = as.vector(c(test.file$ASN_target_pbo, test.file$ASN_target_1,
                                 test.file$ASN_target_2, test.file$ASN_target_3))

  rand_ratio_plot_data_in = data.frame("ratio" = as.vector(t(rand_ratio_mean_table)),
                                       "sd" = as.vector(t(rand_ratio_sd_table)),
                                       "Group" = rep(c("Placebo", paste0("S_", 1:n.active.arm)),
                                                     each = n.look),
                                       "look" = rep(1:n.look, n.arm),
                                       "target" = rep(n.target.asn.vec/sum(n.target.asn.vec),
                                                      each = n.look),
                                       "text" = paste0(round(as.vector(t(rand_ratio_mean_table)), 3)),
                                       "i" = paste0(method_new_name, " mean vector = (",
                                                    response_name, ")"),
                                       "vadj" = rep(c(-0.3, 0.6), each = n.look)
  )
  # if (method_new_name=="DBCD"){
  #   rand_ratio_plot_data_in$vadj = 0
  # }


  rand_ratio_plot_data_com = rbind(rand_ratio_plot_data_com, rand_ratio_plot_data_in)
}

png("plot_data/power_com_props.png",
    width = 3000, height = 3000)
ggplot.fit = ggplot(rand_ratio_plot_data_com, aes(x = look, y = ratio)) +
  geom_point(size = 7, aes(shape=Group))+
  scale_shape_manual(values = c(16, 17, 8, 3))+
  # geom_line(size = 5, aes(linetype=Group))+
  # scale_shape_manual(values = c(16, 17, 8, 3))+
  geom_text(aes(vjust= vadj, label=ifelse(look == max(rand_ratio_plot_data_com$look),as.character(text),'')),
            hjust= -0.3, size = 16) +
  # geom_text(label = text,
  #           hjust= -0.05 ,vjust=0,size = 15) +
  # geom_point()+
  # geom_line()+
  # ylim(0, 0.6)+
  scale_x_continuous(limits = c(1, max(rand_ratio_plot_data_com$look)+20),
                     breaks = c(1, 20, 40, 60))+
  scale_y_continuous(breaks = c(0.1, 0.25,
                                round(unique(rand_ratio_plot_data_com$target), 2), 0.4),
                     limits = c(0.1, 0.4))+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  # geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=1, size = 2)+
  geom_segment(aes(x=1,xend=60,y= unique(rand_ratio_plot_data_com$target)[1],
                   yend= unique(rand_ratio_plot_data_com$target)[1]),
               size = 2, linetype="dashed")+
  geom_segment(aes(x=1,xend=60,y= unique(rand_ratio_plot_data_com$target)[2],
                   yend= unique(rand_ratio_plot_data_com$target)[2]),
               size = 2, linetype="dashed")+
  # geom_hline(yintercept = unique(rand_ratio_plot_data_com$target), size = 2, linetype="dashed") +
  # geom_line(aes(x = look, y = target, color = Group), size = 2 ,linetype="dashed") +
  # geom_hline(aes(yintercept = target, color = group)) +
  #scale_x_continuous(breaks = c(1, 100, 200, 300, 400, sim.RAR.fit$summary.table$n_looks_RAR))+
  # labs(title = "test") +
  ylab ("Proportion of sample size") + xlab("Response adaptive randomization checkpoint") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,0,1,1),units="lines"),
        text = element_text(size=70),
        axis.text.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=.5,face="plain",
                                   margin=margin(20,0,0,0)),
        axis.text.y = element_text(colour="black",size=55,angle=0,hjust=1,vjust=0,face="plain",
                                   margin=margin(0,20,0,0)),
        axis.title.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=0,face="plain",
                                    margin=margin(50,0,0,0)),
        axis.title.y = element_text(colour="black",size=60,angle=90,hjust=.5,vjust=.5,face="plain",
                                    margin=margin(0,50,0,0)),
        # legend.position = "none"
        legend.text = element_text(colour="black", size = 60, face = "plain"),
        legend.title = element_text(colour="black", size = 55, face = "plain"),
        legend.key.size = unit(8,"line"),
        legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 60)
        )+
  theme(plot.margin=unit(c(0,0,2,2),"cm"))+
  # guides(size = guide_legend(20))+
  facet_wrap( ~i, scales = "fixed", ncol = 2, nrow = 3)

print(ggplot.fit)
dev.off()









