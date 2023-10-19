# simulation 1

prior.strengths <- c('None', 'Moderate', 'Strong')
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

df.param = as.data.frame(cbind(rep(prior.strengths, each=3),rep(sig.densities,3)))
df.param = as.data.frame(cbind(rep(df.param$V1, each=length(sig.strengths)),rep(df.param$V2, each=length(sig.strengths)),rep(sig.strengths,9)))
np = dim(df.param)[1]

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

for(i in c(5,6,8,9)) {
	plot_label = switch(EXPR=paste(i), '5'='(a) ', '6'='(b) ', '8'='(c) ', '9'='(d) ')
	prior.strength = df.param[length(sig.strengths)*i, 1]
	sig.density = df.param[length(sig.strengths)*i, 2]
	
	matplot(seq(2, 6, len=5), fdr.mthd[(length(sig.strengths)*(i-1)+2):(length(sig.strengths)*i), ], 
		type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main=paste0(plot_label, prior.strength,', ',sig.density), ylim=c(0.001, 0.1),
		xlab=expression(k[t]),  ylab='FDR')

	abline(h=0.05)

	legend("topleft", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)
}

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

for(i in c(5,6,8,9)) {
	plot_label = switch(EXPR=paste(i), '5'='(a) ', '6'='(b) ', '8'='(c) ', '9'='(d) ')
	prior.strength = df.param[length(sig.strengths)*i, 1]
	sig.density = df.param[length(sig.strengths)*i, 2]
	
	matplot(seq(2, 6, len=5), etp.mthd[(length(sig.strengths)*(i-1)+2):(length(sig.strengths)*i), ], 
		type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main=paste0(plot_label, prior.strength,', ',sig.density), ylim=c(0.001, 1),
		xlab=expression(k[t]), ylab='Power')
	legend("bottomright", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)
}








# simulation 2

pi0.vec<-seq(from=0.5, to=0.9, by=0.1)
mu0.vec<-seq(from=2, to=4, by=0.5)

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(mu0.vec, fdr1.mthd[, c(1,6,7,3,5,2)], type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main="(a) FDR Comparison", xlab=expression(mu), ylab="FDR", ylim=c(0, 0.10))
abline(h=0.05)
legend("bottomright", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)

matplot(mu0.vec, etp1.mthd[, c(1,6,7,3,5,2)], type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main="(b) Power Comparison", xlab=expression(mu), ylab="Power", ylim=c(0, 1))
legend("bottomright", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)

matplot(pi0.vec, fdr2.mthd[, c(1,6,7,3,5,2)], type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main="(c) FDR Comparison", xlab=expression(pi[l]), ylab="FDR", ylim=c(0, 0.10))
abline(h=0.05)
legend("bottomright", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)

matplot(pi0.vec, etp2.mthd[, c(1,6,7,3,5,2)], type="o", pch=1:6, col=c(3,2,6,4,5,7), lwd=2, main="(d) Power Comparison", xlab=expression(pi[l]), ylab="Power", ylim=c(0, 1))
legend("bottomright", c("BH", expression(rho-BH.OR), expression(rho-BH.DD), "CAMT", "LAWS", "Clfdr"), pch=1:6, col=c(3,2,6,4,5,7), lwd=2, cex = 0.8)










# real data

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(seq(0.005, 0.1, 0.005), compMWAS[, c(4,2,3,5,6)], 
        type="o", pch=c(1,3,4,5,6), col=c(3,6,4,5,7), lwd=1.5, main='(a) MWAS',
        xlab='desired FDR level', ylab='num of rejections', ylim=c(0,160))
legend("topleft", c('BH',expression(rho-BH), "CAMT", 'LAWS', 'Clfdr'), pch=c(1,3,4,5,6), col=c(3,6,4,5,7), lwd=1.5, cex = 0.8)

matplot(seq(0.005, 0.1, 0.005), compADHD[, c(4,2,3,5,6)], 
        type="o", pch=c(1,3,4,5,6), col=c(3,6,4,5,7), lwd=1.5, main='(b) ADHD',
        xlab='desired FDR level', ylab='num of rejections')
legend("topleft", c('BH',expression(rho-BH), "CAMT", 'LAWS', 'Clfdr'), pch=c(1,3,4,5,6), col=c(3,6,4,5,7), lwd=1.5, cex = 0.8)





