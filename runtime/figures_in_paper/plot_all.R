library(ggplot2)
library(gridExtra)
library(ggpubr)

load("all_results.Rdata")

fdr1.mthd = fdr1.mthd[, c(1,6,7,3,5,2)]
etp1.mthd = etp1.mthd[, c(1,6,7,3,5,2)]
fdr2.mthd = fdr2.mthd[, c(1,6,7,3,5,2)]
etp2.mthd = etp2.mthd[, c(1,6,7,3,5,2)]


df.setting21.fdr = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(2,4,by=0.5),6), values=c(fdr1.mthd))
fig_sim_21_fdr = ggplot(data = df.setting21.fdr, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab(expression(mu)) + ylab('FDR') + ylim(0,0.1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(A) FDR Comparison (Setting 1)")
#fig_sim_21_fdr

df.setting21.etp = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(2,4,by=0.5),6), values=c(etp1.mthd))
fig_sim_21_etp = ggplot(data = df.setting21.etp, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab(expression(mu)) + ylab('Power') + ylim(0,1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(B) Power Comparison (Setting 1)")
#fig_sim_21_etp

df.setting22.fdr = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(0.5,0.9,by=0.1),6), values=c(fdr2.mthd))
fig_sim_22_fdr = ggplot(data = df.setting22.fdr, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab(expression(pi[l])) + ylab('FDR') + ylim(0,0.1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(C) FDR Comparison (Setting 2)")
#fig_sim_22_fdr

df.setting22.etp = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(0.5,0.9,by=0.1),6), values=c(etp2.mthd))
fig_sim_22_etp = ggplot(data = df.setting22.etp, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab(expression(pi[l])) + ylab('Power') + ylim(0,1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(D) Power Comparison (Setting 2)")
#fig_sim_22_etp

gc()
fig3 = ggpubr::ggarrange(fig_sim_21_fdr, fig_sim_21_etp, fig_sim_22_fdr, fig_sim_22_etp, ncol=2, nrow=2, common.legend = T, legend="bottom")
ggsave(fig3, file="fig3.eps", device="eps", width=7, height=7)



prior.strengths <- c('None', 'Moderate', 'Strong')
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

df.param = as.data.frame(cbind(rep(prior.strengths, each=3),rep(sig.densities,3)))
df.param = as.data.frame(cbind(rep(df.param$V1, each=length(sig.strengths)),rep(df.param$V2, each=length(sig.strengths)),rep(sig.strengths,9)))
np = dim(df.param)[1]

list.fig.sim1 = list()
for(i in c(5,6,8,9)) {
plot_label = switch(EXPR=paste(i), '5'='(A) ', '6'='(B) ', '8'='(B) ', '9'='(B) ')
prior.strength = df.param[length(sig.strengths)*i, 1]
sig.density = df.param[length(sig.strengths)*i, 2]

df.setting1i.fdr = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(2,6,by=1),6), values=c(fdr.mthd[(length(sig.strengths)*(i-1)+2):(length(sig.strengths)*i), ]))
df.setting1i.etp = data.frame(procs=factor((0:29)%/%5+1), inds=rep(seq(2,6,by=1),6), values=c(etp.mthd[(length(sig.strengths)*(i-1)+2):(length(sig.strengths)*i), ]))
fig_sim_1i_fdr = ggplot(data = df.setting1i.fdr, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab('desired FDR level') + ylab('num of rejections') + xlab(expression(k[t])) + ylab('FDR') + ylim(0,0.1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle=paste0(plot_label, prior.strength, ', ', sig.density))
fig_sim_1i_etp = ggplot(data = df.setting1i.etp, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line(linewidth=0.7) + geom_point(size=1.5) + scale_shape_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,2,6,4,5,3), labels=c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab(expression(k[t])) + ylab('Power') + ylim(0,1) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle=paste0(plot_label, prior.strength, ', ', sig.density))
#print(fig_sim_1i_fdr)
#print(fig_sim_1i_etp)

list.fig.sim1i = list(fig_sim_1i_fdr, fig_sim_1i_etp)
list.fig.sim1 = append(list.fig.sim1, list.fig.sim1i)
}

gc()
fig1 = ggpubr::ggarrange(list.fig.sim1[[1]], list.fig.sim1[[3]], list.fig.sim1[[5]], list.fig.sim1[[7]], ncol=2, nrow=2, common.legend = T, legend="bottom")
fig2 = ggpubr::ggarrange(list.fig.sim1[[2]], list.fig.sim1[[4]], list.fig.sim1[[6]], list.fig.sim1[[8]], ncol=2, nrow=2, common.legend = T, legend="bottom")
ggsave(fig1, file="fig1.eps", device="eps", width=7, height=7)
ggsave(fig2, file="fig2.eps", device="eps", width=7, height=7)



df.rd1 = data.frame(procs=factor((0:99)%/%20+1), inds=rep(seq(0.005,0.1,by=0.005),5), values=c(compMWAS[, c(4,2,3,5,6)]))
fig_rd1 = ggplot(data = df.rd1, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line() + geom_point() + scale_shape_manual(values = c(1,6,4,5,3), labels=c("BH", expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,6,4,5,3), labels=c("BH", expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab('desired FDR level') + ylab('num of rejections') + ylim(0,160) + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(A) MWAS")
#fig_rd1
df.rd2 = data.frame(procs=factor((0:99)%/%20+1), inds=rep(seq(0.005,0.1,by=0.005),5), values=c(compADHD[, c(4,2,3,5,6)]))
fig_rd2 = ggplot(data = df.rd2, mapping = aes(x = inds, y = values, colour = procs, group = procs, shape = procs)) + geom_line() + geom_point() + scale_shape_manual(values = c(1,6,4,5,3), labels=c("BH", expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + scale_color_manual(values = c(1,6,4,5,3), labels=c("BH", expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr")) + xlab('desired FDR level') + ylab('num of rejections') + theme(legend.title=element_blank(), plot.caption=element_text(hjust=0,vjust=0)) + labs(subtitle="(B) ADHD")
#fig_rd2
fig4 = ggpubr::ggarrange(fig_rd1, fig_rd2, ncol=2, nrow=1, common.legend = T, legend="bottom")
ggsave(fig4, file="fig4.eps", device="eps", width=7, height=3.5)






