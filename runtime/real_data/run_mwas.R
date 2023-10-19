library(splines)
library(CAMT)
library(rhoBH)

registerDoParallel(50)

datapvx1 = readRDS('./data/mwas.p.value.rds')
pvals = datapvx1[, 1]
x1 = datapvx1[, 2]
x1.ns = ns(x1, df=6)
x1 = matrix(x1, length(x1), 1)


m = length(pvals)
pvals[pvals < 1e-8] = 1e-8
stats_a = qnorm(1 - pvals/2) * (2*rbinom(m, 1, 0.5)-1)
#sum(is.infinite(stats_a))


CAMT.obj = camt.fdr(pvals = pvals, pi0.var = x1.ns, f1.var = x1.ns, 
      alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)

rbhres = rhoBH.prl(stats_a, x1, x1, 0.05)

dim_s1 = dim(x1)[2]
bh.th = bh.func(pvals, 0.9)$th
h = matrix(Hns(x1), dim_s1, dim_s1)
if(sum(diag(h)) == 0) {
    h = diag(dim_s1)
}
pishat = pisLAWS.prl(s1=x1, pv2=pvals, s2=x1, tau=bh.th, h=h)


compCR = matrix(0, 20, 6)
for(i in 1:20) {
    q = i/20*0.1

    rhow = rbhres$rhow1
    fdps = rbhres$fdp1
    st.rhow = sort(rhow)
    if (sum(fdps<=q)==0)
    {
        k<-0
            pk<-0
    } else {
        k<-max(which(fdps<=q))
        pk<-st.rhow[k]
    }
    sum1 = sum(rhow<=pk)
    rhow = rbhres$rhow2
    fdps = rbhres$fdp2
    st.rhow = sort(rhow)
    if (sum(fdps<=q)==0)
    {
        k<-0
        pk<-0
    } else {
        k<-max(which(fdps<=q))
        pk<-st.rhow[k]
    }
    compCR[i, 1] = q
    compCR[i, 2] = sum1+sum(rhow<=pk)
    compCR[i, 3] = sum(CAMT.obj$fdr<=q)
    compCR[i, 4] = sum(bh.func(pvals, q)$de)
    compCR[i, 5] = sum(law.func(pvals, pishat, q)$de)

    rhow = rbhres$rhow1
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
    sum1=sum(clfdrs<=pk)
    rhow = rbhres$rhow2
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
    compCR[i, 6] = sum1 + sum(clfdrs<=pk)
}

compMWAS = compCR

save(compMWAS, rbhres, CAMT.obj, file='mwas_results.Rdata')







