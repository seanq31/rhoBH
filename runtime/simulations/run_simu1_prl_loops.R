rm(list = ls())

library(CAMT)
library(qvalue)
library(splines)
library(rhoBH)


simulate.data.rbh <- function (
		paras.mapping       = paras.mapping.func(),
		covariate.strength  = c('None', 'Moderate', 'Strong'),  
		covariate.model     = c('pi0', 'f1', 'both'),
		covariate.dist      = c('Normal', 'Uniform', 'T'),
		null.model          = c('Unif', 'Left', 'Right'),
		skewness            = 0.15,
		f1.sd               = 1.0,
		feature.no          = 10000,
		sig.dist            = c('Normal', 'Gamma'), 
		sig.density         = c('Low', 'Medium', 'High'),
		sig.strength        = c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
		cor.struct          = c('None', 'Block', 'AR1'),
		cor.rho             = 0,
		cor.sign            = c('PosCor', 'PosNegCor'),
		cor.nblock          = 500) {

	
	sig.dist <- match.arg(sig.dist)
	sig.density <- match.arg(sig.density)
	sig.strength <- match.arg(sig.strength)
	covariate.strength <- match.arg(covariate.strength)
	covariate.model <- match.arg(covariate.model)
	covariate.dist <- match.arg(covariate.dist)
	null.model <- match.arg(null.model)
	cor.struct <- match.arg(cor.struct)
	cor.sign <- match.arg(cor.sign)
	
	# These could be altered for other simulation purpose

	sig.densities <- paras.mapping[['sig.densities']]
	sig.strengths <- paras.mapping[['sig.strengths']]
	pi.strengths <- paras.mapping[['pi.strengths']]
	f1.strengths <- paras.mapping[['f1.strengths']]
	
	if (covariate.model %in% c('pi0')) {
		pi.k <- pi.strengths[covariate.strength]
		f1.k <- 0
	} 
	if (covariate.model %in% c('f1')) {
		pi.k <- 0
		f1.k <- f1.strengths[covariate.strength]
	} 
	if (covariate.model %in% c('both')) {
		pi.k <- pi.strengths[covariate.strength]
		f1.k <- f1.strengths[covariate.strength]
	} 
	
	# Generate covariate for pi0 and f1
    # We simulate different covariates
	if (covariate.dist == 'Normal') {
		x1 <- rnorm(feature.no)
		x2 <- rnorm(feature.no)
	} 
	
	if (covariate.dist == 'Uniform') {
		x1 <- scale(runif(feature.no))
		x2 <- scale(runif(feature.no))
	}
	
	if (covariate.dist == 'T') {
		x1 <- rt(feature.no, df = 5)
		x2 <- rt(feature.no, df = 5)
	}
	
	eta <- sig.densities[sig.density] + x1 * pi.k
	pi <- exp(eta) / (1 + exp(eta))
	
	truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	ns <- sum(truth)
	
	# Generate the correlated signals under the null
	if (cor.struct != 'None') {
		if (null.model != 'Unif') {
			cat('Correlated signals will be simulated under the uniform null only! "null.model" will be automatically reset to be "Unif". \n')
		} else {
			if (cor.struct == 'Block') {
				obj <- Tmat(feature.no, cor.nblock, cor.rho)
				bs <- feature.no / cor.nblock
				if (cor.sign == 'PosCor') {
					stats <- as.vector(obj$T1 %*% matrix(rnorm(feature.no), bs, cor.nblock))
					
				}
				if (cor.sign == 'PosNegCor') {
					stats <- as.vector(obj$T2 %*% matrix(rnorm(feature.no), bs, cor.nblock))
				}
			} 
			if (cor.struct == 'AR1') {
				if (cor.sign == 'PosCor') {
					stats <- as.vector(arima.sim(n = feature.no, list(ar = c(cor.rho)), n.start = 100000) * sqrt(1 - cor.rho^2))
				}
				if (cor.sign == 'PosNegCor') {
					stats <- as.vector(arima.sim(n = feature.no, list(ar = c(-cor.rho)), n.start = 100000) * sqrt(1 - cor.rho^2))
				}
			}
			skewness <- 0
		}
		
	} else {
		if (null.model == 'Unif') {
			stats <- rnorm(feature.no)
			skewness <- 0
		}
		if (null.model == 'Left') {
			stats <- rnorm(feature.no, -abs(skewness), 1)
			skewness <-  -abs(skewness)
		}
		if (null.model == 'Right') {
			stats <- rnorm(feature.no, abs(skewness), 1)
			skewness <- abs(skewness)
		}
		
	}
	
	if (sig.dist == 'Normal') {
		#f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2
		f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2 
				
		if (null.model == 'Unif') {
			stats[truth == 1] <- stats[truth == 1] * f1.sd + f1.mean[truth == 1]
		} else {
			stats[truth == 1] <- rnorm(ns) * f1.sd + f1.mean[truth == 1]
		}
		
		lfdr <- pi * dnorm(stats, mean = skewness) / (pi * dnorm(stats, mean = skewness) + (1 - pi) * dnorm(stats, f1.mean, sd = f1.sd))
		
		out <- sort(lfdr, index.return = TRUE)
		fdr <- cumsum(out$x) / (1:(length(lfdr)))
		fdr <- fdr[order(out$ix)]
	}
	
	if (sig.dist == 'Gamma') {
		
		# Under gamma distribution (a) the nulls and non-nulls are not correlated (b) we do not shrink the variance in the covariate model 'f1' and 'both'
		#f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2
		f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2
		
		# Match the mean and variance of the normal counterpart
		stats[truth == 1] <- rgamma(ns, shape = 2, scale = 1 / sqrt(2)) - sqrt(2) + f1.mean[truth == 1]
		
		# Calculate local false discovery rate
		lfdr <- pi * dnorm(stats, mean = skewness) / (pi * dnorm(stats, mean = skewness) + (1 - pi) * 
					dgamma(stats - f1.mean + sqrt(2), shape = 2, scale = 1 / sqrt(2))) 
		
		out <- sort(lfdr, index.return = TRUE)
		fdr <- cumsum(out$x) / (1:(length(lfdr)))
		fdr <- fdr[order(out$ix)]
		
	}
	
	# Generate p values - one sided p-value
	#pvals <- (1 - pnorm(stats))
	pvals <- 2*pnorm(-abs(stats), 0, 1)
	
	return(list(pvals = pvals, pi0.var = x1, f1.var = x2, pis = 1 - pi, truth = truth, fdr = fdr, lfdr = lfdr, stats = stats, mu1 = f1.mean, sig1 = f1.sd))
}


paras.mapping.func <- function (
		#sig.densities = c(Low=3.5, Medium=2.5, High=1.5),
		#sig.strengths = log(seq(exp(2), exp(3), len=6)),
		#pi.strengths = c(None=0.0, Moderate=1.0, Strong=1.5),
		#f1.strengths = c(None=0.0, Moderate=0.25, Strong=0.5)
		sig.densities = c(Low=3, Medium=2, High=1),
		sig.strengths = c(1, seq(2, 6, len=5)),
		pi.strengths = c(None=0.0, Moderate=1.5, Strong=2.5),
		f1.strengths = c(None=0.0, Moderate=0.4, Strong=0.6)
) {
	
	names(sig.strengths) <- paste0('L', 1:6)
	names(sig.densities) <- c('Low', 'Medium', 'High')
	names(pi.strengths) <- c('None', 'Moderate', 'Strong')
	names(f1.strengths) <- c('None', 'Moderate', 'Strong')
	return(list(sig.strengths = sig.strengths, sig.densities = sig.densities,
					pi.strengths = pi.strengths, f1.strengths = f1.strengths))
	
}



registerDoParallel(50)

m = 5000
nrep = 100
q = 0.05

prior.strengths <- c('None', 'Moderate', 'Strong')
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

df.param = as.data.frame(cbind(rep(prior.strengths, each=3),rep(sig.densities,3)))
df.param = as.data.frame(cbind(rep(df.param$V1, each=length(sig.strengths)),rep(df.param$V2, each=length(sig.strengths)),rep(sig.strengths,9)))
np = dim(df.param)[1]

bh.fdp.np.nrep = matrix(0, np, nrep)
rbh.fdp.np.nrep = matrix(0, np, nrep)
rbh.dd.fdp.np.nrep = matrix(0, np, nrep)
camt.fdp.np.nrep = matrix(0, np, nrep)

bh.ntp.np.nrep = matrix(0, np, nrep)
rbh.ntp.np.nrep = matrix(0, np, nrep)
rbh.dd.ntp.np.nrep = matrix(0, np, nrep)
camt.ntp.np.nrep = matrix(0, np, nrep)

fdr.mthd = matrix(0, np, 6)
etp.mthd = matrix(0, np, 6)

for(i in c(26:30, 32:36, 44:48, 50:54)) {
	prior.strength = df.param[i, 1]
	sig.density = df.param[i, 2]
	sig.strength = df.param[i, 3]
	print(df.param[i, ])

	#forOut = matrix(0, 12, nrep)
	forOut = foreach(j = 1:nrep, .combine = 'rbind') %dopar% {
	#for(j in 1:nrep) {
		cat(j)
		data <- simulate.data.rbh(covariate.strength = prior.strength, sig.density = sig.density, sig.strength = sig.strength,
			sig.dist = 'Normal', feature.no = m, null.model = 'Unif', covariate.model = 'both', covariate.dist = 'Normal',
			cor.struct = 'None', cor.nblock = 500, cor.sign = 'PosCor', cor.rho = 0)

		x = data$stats
		theta = data$truth
		pis = data$pis
		mu1 = data$mu1
		sig1 = rep(data$sig1, m)
		pv = data$pvals
		s_aux1 = matrix(data$pi0.var,m,1)
		s_aux2 = matrix(data$f1.var,m,1)

		# BH
		bh.res<-bh.func(pv, q)
		bh.de<-bh.res$de
		bh.fdp<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
		bh.ntp<-sum(theta*bh.de)/sum(theta)	

		# CAMT
		camt.res = camt.fdr(pvals = data$pvals, pi0.var = s_aux1, f1.var = s_aux2, alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)
		camt.fdr = camt.res$fdr
		if (sum(camt.fdr <= q) == 0) {
			camt.fdp <- 0
			camt.ntp <- 0
		} else {
			camt.fdp <- mean(theta[camt.fdr <= q] == 0)
			camt.ntp <- sum(theta[camt.fdr <= q]) / sum(theta)
		}

		# rhoBH: OR
		rbh.res = rhoBH.OR(x, mu1, sig1, pis/(1-pis), q)
		rbh.de<-rbh.res$de
		rbh.fdp<-sum((1-theta)*rbh.de)/max(sum(rbh.de), 1)
		rbh.ntp<-sum(theta*rbh.de)/sum(theta)

		# rhoBH: DD
		rbh.dd.res = rhoBH.prl(x, s_aux1, s_aux2, q)
		rbh.dd.de = rbh.dd.res$de
		rbh.dd.fdp<-sum((1-theta)*rbh.dd.de)/max(sum(rbh.dd.de), 1)
		rbh.dd.ntp<-sum(theta*rbh.dd.de)/sum(theta)

		# LAWS: DD
		dim_s1 = dim(s_aux1)[2]
		bh.th<-bh.func(pv, 0.9)$th
    	h = matrix(Hns(s_aux1), dim_s1, dim_s1)
    	if(sum(diag(h)) == 0) {
      		h = diag(dim_s1)
    	}
		pis.hat<-pisLAWS.prl(s1=s_aux1, pv2=pv, s2=s_aux1, tau=bh.th, h=h)
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


		#forOut[, j] = c(bh.fdp, rbh.fdp, rbh.dd.fdp, camt.fdp,
		#	bh.ntp, rbh.ntp, rbh.dd.ntp, camt.ntp,
		#	law.dd.fdp, law.dd.ntp, clfdr.dd.fdp, clfdr.dd.ntp)
		
		c(bh.fdp, rbh.fdp, rbh.dd.fdp, camt.fdp,
			bh.ntp, rbh.ntp, rbh.dd.ntp, camt.ntp,
			law.dd.fdp, law.dd.ntp, clfdr.dd.fdp, clfdr.dd.ntp)
	}

	fdr.mthd[i, ] = colMeans(forOut)[c(1:4,9,11)]
	etp.mthd[i, ] = colMeans(forOut)[c(5:8,10,12)]

	str = paste('./simu1_prl', i, 'm', m, Sys.time(), '.Rdata', sep='_')
	save(fdr.mthd, etp.mthd, file = str)
}

str = paste('./simu1_prl', 'm', m, 'nrep', nrep, Sys.time(), '.Rdata', sep='_')
save(fdr.mthd, etp.mthd, file = str)




