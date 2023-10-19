library(splines)
library(CAMT)
library(R.matlab)
library(abind)
library(readr)
library(reticulate)
library(rhoBH)

registerDoParallel(50)

numpy = import('numpy')

dt_train = numpy$load ('./data/dt_train_swssd.npy')
dt_test = numpy$load  ('./data/dt_test_swssd.npy')
label_train = read_csv('./data/label_for_train.csv', col_names = FALSE)
label_test = read_csv ('./data/label_for_test.csv', col_names = FALSE)
dt_all = abind(dt_train, dt_test, along = 1)
label_all = abind(label_train, label_test, along = 1)

x1 = dt_all[which(label_all==0), , , ]
x2 = dt_all[which(label_all==1), , , ]
n1 = dim(x1)[1]
n2 = dim(x2)[1]
var1 = apply(x1, 2:4, var)
x12 = apply(x1, 2:4, mean)
var2 = apply(x2, 2:4, var)
x22 = apply(x2, 2:4, mean)

x_d = x12 - x22
x_t = x12 - x22
x_t[which(x_t!=0, arr.ind=T)] = x_d[which(x_t!=0, arr.ind=T)]/sqrt(var1[which(x_t!=0, arr.ind=T)]/n1 + var2[which(x_t!=0, arr.ind=T)]/n2)

dims = dim(x_t)
m = prod(dims)
z_values = rep(0, m)
s_aux = matrix(0, m, 3)
cnts = 0

for(i in 1:dims[1]) {
    for(j in 1:dims[2]) {
        for(k in 1:dims[3]) {
            cnts = cnts + 1
            z_values[cnts] = x_t[i,j,k]
            s_aux[cnts,] = c(i,j,k)
        }
    }
}

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



