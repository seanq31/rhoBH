library(splines)
library(CAMT)
library(rhoBH)

registerDoParallel(50)

load("./data/ADHD_data.Rdata")

temp_aux = ns(s_aux, 6)
pvals = 2*pnorm(-abs(z_values), 0, 1)

CAMT.obj <- camt.fdr(pvals = pvals, pi0.var = temp_aux, f1.var = temp_aux, 
      alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)

rbhres = rhoBH.prl(z_values, s_aux, s_aux, 0.05)

dim_s1 = dim(s_aux)[2]
bh.th<-bh.func(pvals, 0.9)$th
h = matrix(Hns(s_aux), dim_s1, dim_s1)
if(sum(diag(h)) == 0) {
  h = diag(dim_s1)
}
pishat<-pisLAWS.prl(s1=s_aux, pv2=pvals, s2=s_aux, tau=bh.th, h=h)

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
    sum1 = sum(clfdrs<=pk)
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

compADHD = compCR

save(compADHD, CAMT.obj, rbhres, file='ADHD_results.Rdata')



