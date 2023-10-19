rm(list = ls())

library(CAMT)
library(qvalue)
library(splines)
library(rhoBH)

registerDoParallel(50)

nrep = 100



## Block patterns

## Simulation 1: Effect sizes

m<-5000
mu0.vec<-seq(from=2, to=4, by=0.5)
pis<-rep(0.01, m)
pis[1001:1200]<-0.9
pis[2001:2200]<-0.9
pis[3001:3200]<-0.6
pis[4001:4200]<-0.6
q<-0.05
np<-length(mu0.vec)
locs<-1:m

# Methods

fdr1.mthd = matrix(0, np, 7)
etp1.mthd = matrix(0, np, 7)

for(i in 1:np) {
	cat("\n", "iteration i= ", i, "\n", "iteration j=")
	mu0<-mu0.vec[i]
		
	forOut = matrix(0, nrep, 14)
	#forOut = foreach(j = 1:nrep, .combine = 'rbind') %dopar% {
	for(j in 1:nrep) {
		cat(j)
    	theta<-rbinom(m, size=1, prob=pis)
    	pii<-sum(theta)/m
		x0<-rnorm(m, mean=0, sd=1)
		x1<-rnorm(m, mean=mu0, sd=1)
		x<-(1-theta)*x0+theta*x1
		pv<-2*pnorm(-abs(x), 0, 1)
		s_aux = matrix(1:m,m,1)
		
		# BH
		bh.res<-bh.func(pv, q)
		bh.de<-bh.res$de
		bh.fdp<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
		bh.ntp<-sum(theta*bh.de)/sum(theta)	

		law.or.fdp = 0
		law.or.ntp = 0
		
		# rhoBH: OR
		rbh.res = rhoBH.OR(x, rep(mu0, m), rep(1, m), pis/(1-pis), q)
		rbh.de<-rbh.res$de
		rbh.fdp<-sum((1-theta)*rbh.de)/max(sum(rbh.de), 1)
		rbh.ntp<-sum(theta*rbh.de)/sum(theta)

		# rhoBH: DD
		rbh.dd.res = rhoBH.prl(x, s_aux, s_aux, q)
		rbh.dd.de = rbh.dd.res$de
		rbh.dd.fdp<-sum((1-theta)*rbh.dd.de)/max(sum(rbh.dd.de), 1)
		rbh.dd.ntp<-sum(theta*rbh.dd.de)/sum(theta)

		# LAWS: DD
		dim_s1 = dim(s_aux)[2]
		bh.th<-bh.func(pv, 0.9)$th
    	h = matrix(Hns(s_aux), dim_s1, dim_s1)
    	if(sum(diag(h)) == 0) {
      		h = diag(dim_s1)
    	}
		pis.hat<-pisLAWS.prl(s1=s_aux, pv2=pv, s2=s_aux, tau=bh.th, h=h)
		law.dd.res<-law.func(pvs=pv, pis.hat, q)
		law.dd.de<-law.dd.res$de
		law.dd.fdp<-sum((1-theta)*law.dd.de)/max(sum(law.dd.de), 1)
		law.dd.ntp<-sum(theta*law.dd.de)/sum(theta)
		
		# CLFDR: DD
		clfdr.dd.de = rep(0, m)
		rhow = rbh.dd.res$rhow1
    	clfdrs = rhow / (1 + rhow)
    	st.clfdrs = sort(clfdrs)
    	fdps = cumsum(st.clfdrs) / (1 : length(clfdrs))
    	if (sum(fdps<=q)==0)
    	    {
    	        k<-0
    	        pk<-1
    	    } else {
    	        k<-max(which(fdps<=q))
    	        pk<-st.clfdrs[k]
    	    }
    	clfdr.dd.de[rbh.dd.res$rdsp1][which(clfdrs <= pk)] = 1
    	rhow = rbh.dd.res$rhow2
    	clfdrs = rhow / (1 + rhow)
    	st.clfdrs = sort(clfdrs)
    	fdps = cumsum(st.clfdrs) / (1 : length(clfdrs))
    	if (sum(fdps<=q)==0)
    	{
    	    k<-0
    	    pk<-1
    	} else {
    	    k<-max(which(fdps<=q))
    	    pk<-st.clfdrs[k]
    	}
    	clfdr.dd.de[rbh.dd.res$rdsp2][which(clfdrs <= pk)] = 1
		#clfdr.dd.res = CLFDR.DD(x, s_aux1, s_aux2, q)
		#clfdr.dd.de = clfdr.dd.res$de
		clfdr.dd.fdp<-sum((1-theta)*clfdr.dd.de)/max(sum(clfdr.dd.de), 1)
		clfdr.dd.ntp<-sum(theta*clfdr.dd.de)/sum(theta)

		# CAMT
		camt.res = camt.fdr(pvals = pv, pi0.var = ns(s_aux,6), f1.var = ns(s_aux,6), alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)
		camt.fdr = camt.res$fdr
		if (sum(camt.fdr <= q) == 0) {
			camt.fdp <- 0
			camt.ntp <- 0
		} else {
			camt.fdp <- mean(theta[camt.fdr <= q] == 0)
			camt.ntp <- sum(theta[camt.fdr <= q]) / sum(theta)
		}

		forOut[j, ] = c(bh.fdp, clfdr.dd.fdp, camt.fdp, law.or.fdp, law.dd.fdp, rbh.fdp, rbh.dd.fdp,
			bh.ntp, clfdr.dd.ntp, camt.ntp, law.or.ntp, law.dd.ntp, rbh.ntp, rbh.dd.ntp)
	}

	fdr1.mthd[i, ] = colMeans(forOut)[1:7]
	etp1.mthd[i, ] = colMeans(forOut)[8:14]

	str = paste('./simu2_prl_mid', i, Sys.time(), '.Rdata', sep='_')
	save(fdr1.mthd, etp1.mthd, file = str)
}

str = paste('./simu2_prl_mid', 'm', m, 'nrep', nrep, Sys.time(), '.Rdata', sep='_')
save(fdr1.mthd, etp1.mthd, file = str)





## Simulation 2: Sparsity levels

pi0.vec<-seq(from=0.5, to=0.9, by=0.1)
pis<-rep(0.01, m)
np<-length(pi0.vec)
mu0<-3

# Methods
fdr2.mthd = matrix(0, np, 7)
etp2.mthd = matrix(0, np, 7)

for(i in 1:np) {
	cat("\n", "iteration i= ", i, "\n", "iteration j=")
	pi0<-pi0.vec[i]
	pis[1001:1200]<-pi0
	pis[2001:2200]<-pi0
	pis[3001:3200]<-pi0
	pis[4001:4200]<-pi0
		
	forOut = matrix(0, nrep, 14)
	#forOut = foreach(j = 1:nrep, .combine = 'rbind') %dopar% {
	for(j in 1:nrep) {
		cat(j)
    	theta<-rbinom(m, size=1, prob=pis)
    	pii<-sum(theta)/m
		x0<-rnorm(m, mean=0, sd=1)
		x1<-rnorm(m, mean=mu0, sd=1)
		x<-(1-theta)*x0+theta*x1
		pv<-2*pnorm(-abs(x), 0, 1)
		s_aux = matrix(1:m,m,1)

		# BH
		bh.res<-bh.func(pv, q)
		bh.de<-bh.res$de
		bh.fdp<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
		bh.ntp<-sum(theta*bh.de)/sum(theta)	

		law.or.fdp = 0
		law.or.ntp = 0
		
		# rhoBH: OR
		rbh.res = rhoBH.OR(x, rep(mu0, m), rep(1, m), pis/(1-pis), q)
		rbh.de<-rbh.res$de
		rbh.fdp<-sum((1-theta)*rbh.de)/max(sum(rbh.de), 1)
		rbh.ntp<-sum(theta*rbh.de)/sum(theta)

		# rhoBH: DD
		rbh.dd.res = rhoBH.prl(x, s_aux, s_aux, q)
		rbh.dd.de = rbh.dd.res$de
		rbh.dd.fdp<-sum((1-theta)*rbh.dd.de)/max(sum(rbh.dd.de), 1)
		rbh.dd.ntp<-sum(theta*rbh.dd.de)/sum(theta)

		# LAWS: DD
		dim_s1 = dim(s_aux)[2]
		bh.th<-bh.func(pv, 0.9)$th
    	h = matrix(Hns(s_aux), dim_s1, dim_s1)
    	if(sum(diag(h)) == 0) {
      		h = diag(dim_s1)
    	}
		pis.hat<-pisLAWS.prl(s1=s_aux, pv2=pv, s2=s_aux, tau=bh.th, h=h)
		law.dd.res<-law.func(pvs=pv, pis.hat, q)
		law.dd.de<-law.dd.res$de
		law.dd.fdp<-sum((1-theta)*law.dd.de)/max(sum(law.dd.de), 1)
		law.dd.ntp<-sum(theta*law.dd.de)/sum(theta)
		
		# CLFDR: DD
		clfdr.dd.de = rep(0, m)
		rhow = rbh.dd.res$rhow1
    	clfdrs = rhow / (1 + rhow)
    	st.clfdrs = sort(clfdrs)
    	fdps = cumsum(st.clfdrs) / (1 : length(clfdrs))
    	if (sum(fdps<=q)==0)
    	    {
    	        k<-0
    	        pk<-1
    	    } else {
    	        k<-max(which(fdps<=q))
    	        pk<-st.clfdrs[k]
    	    }
    	clfdr.dd.de[rbh.dd.res$rdsp1][which(clfdrs <= pk)] = 1
    	rhow = rbh.dd.res$rhow2
    	clfdrs = rhow / (1 + rhow)
    	st.clfdrs = sort(clfdrs)
    	fdps = cumsum(st.clfdrs) / (1 : length(clfdrs))
    	if (sum(fdps<=q)==0)
    	{
    	    k<-0
    	    pk<-1
    	} else {
    	    k<-max(which(fdps<=q))
    	    pk<-st.clfdrs[k]
    	}
    	clfdr.dd.de[rbh.dd.res$rdsp2][which(clfdrs <= pk)] = 1
		#clfdr.dd.res = CLFDR.DD(x, s_aux1, s_aux2, q)
		#clfdr.dd.de = clfdr.dd.res$de
		clfdr.dd.fdp<-sum((1-theta)*clfdr.dd.de)/max(sum(clfdr.dd.de), 1)
		clfdr.dd.ntp<-sum(theta*clfdr.dd.de)/sum(theta)

		# CAMT
		camt.res = camt.fdr(pvals = pv, pi0.var = ns(s_aux,6), f1.var = ns(s_aux,6), alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)
		camt.fdr = camt.res$fdr
		if (sum(camt.fdr <= q) == 0) {
			camt.fdp <- 0
			camt.ntp <- 0
		} else {
			camt.fdp <- mean(theta[camt.fdr <= q] == 0)
			camt.ntp <- sum(theta[camt.fdr <= q]) / sum(theta)
		}

		forOut[j, ] = c(bh.fdp, clfdr.dd.fdp, camt.fdp, law.or.fdp, law.dd.fdp, rbh.fdp, rbh.dd.fdp,
			bh.ntp, clfdr.dd.ntp, camt.ntp, law.or.ntp, law.dd.ntp, rbh.ntp, rbh.dd.ntp)
	}

	fdr2.mthd[i, ] = colMeans(forOut)[1:7]
	etp2.mthd[i, ] = colMeans(forOut)[8:14]

	str = paste('./simu2_prl_both', i, Sys.time(), '.Rdata', sep='_')
	save(fdr1.mthd, etp1.mthd, fdr2.mthd, etp2.mthd, file = str)
}


str = paste('./simu2_prl_both', 'm', m, 'nrep', nrep, Sys.time(), '.Rdata', sep='_')
save(fdr1.mthd, etp1.mthd, fdr2.mthd, etp2.mthd, file = str)









