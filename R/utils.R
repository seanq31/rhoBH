#' A correct version of the Hns function in package ks.
#' 
#' A correct version of the Hns function in package ks, all definitions are the same as in the original function.
#'
#' 
#' @export
Hns = function (x, deriv.order=0) 
{
    if (is.vector(x)) {
        d <- 1
        n <- length(x)
    }
    else {
        n <- nrow(x)
        d <- ncol(x)
    }
    r <- deriv.order
    H <- (4/(n * (d + 2 * r + 2)))^(2/(d + 2 * r + 4)) * var(x)
    return(H)
}


calDist = function(row1, row2, h) {
    if(length(row1) == 1) {
        dnorm(row1-row2, 0, sqrt(h))
    } else {
        dmvn(row1-row2, rep(0, length(row1)), h)
    }
}


calPest = function(s1row, s2, h, inds) {
    kht = apply(s2, 1, function(s2row) calDist(s1row, s2row, h) )

    sum(kht[inds])/sum(kht)
}


#' Implementing the pi estimation in the paper.
#' 
#' Implementing the pi estimation in the paper.
#'
#' @param s1 the covariates to estimate non-null probability for.
#' @param pv2 the p-values used for estimation.
#' @param s2 the covariates corresponding to the p-values pv2.
#' @param tau threshold for null p-values.
#' @param h bandwidth matrix H.
#'
#' @return Estimated non-null probability for each covariate in s1.
#' 
#' @export
pisLAWS<- function(s1, pv2, s2, tau=0.1, h=0.3)
{
    m1 <- dim(s1)[1]
    m2 = length(pv2)
    dim_s1 = dim(s1)[2]
    p.est <-rep(0, m1)

    inds = which(pv2>=tau)
    p.est = apply(s1, 1, function(s1row) calPest(s1row, s2, h, inds)) / (1-tau)

    p.est[which(p.est>=1)] <- 1-1e-16
    p.est[which(p.est<=1e-16)] <-1e-16

    return(1-p.est)
}

calDenkdnorm = function(rowxs, s_aux, hs, x, hx, clfdrnest) {
    rowx = rowxs[1]
    rows = rowxs[-1]

    nestw = apply(s_aux, 1, function(row) calDist(row, rows, hs) ) * (1 - clfdrnest)
    nestw[is.na(nestw)] = 0
    nestw[is.infinite(nestw)] = 1e16
    if(sum(abs(nestw)) == 0) {
        nestw = nestw + 1
    }

    sum(dnorm(x - rowx, 0, hx) * nestw) / sum(nestw)
}


calDenkdnormSim = function(rowpiws, s_aux22, hs, x2, hx, clfdrnest, muhat, sighat, st.rhow) {
    pis.hati = rowpiws[1]
    wsi = rowpiws[2]
    rows = rowpiws[c(-1,-2)]

    nestw = apply(s_aux22, 1, function(row) calDist(row, rows, hs) ) * (1 - clfdrnest)
    nestw[is.na(nestw)] = 0
    nestw[is.infinite(nestw)] = 1e16
    if(sum(abs(nestw)) == 0) {
        nestw = nestw + 1
    }

    x_sim = rnorm(1000, muhat, sighat)
    fdenest = density(x2, from = min(x2,x_sim)*1.1-max(x2,x_sim)*0.1 , to = max(x2,x_sim)*1.1-min(x2,x_sim)*0.1, n = 1000, weights = (nestw/sum(nestw)))
    rho_sim = dnorm(x_sim, muhat, sighat) / lin.itp(x_sim, fdenest$x, fdenest$y)
    ecdfrho = ecdf(rho_sim)

    (1 - pis.hati) * ecdfrho(st.rhow * wsi)
}


rhoBHmultif1 = function(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, muhat=0, sighat=1) {    
    m1 = length(x1)
    m2 = length(x2)
    dim_s2 = dim(s_aux21)[2]
    ws = pis.hat1 / (1 - pis.hat1)

    hx = density(x2,from = min(x2)-5,to=max(x2)+5,n=1000)$bw
    hs = matrix(Hns(s_aux22), dim_s2, dim_s2)
    if(sum(diag(hs)) == 0) {
        hs = diag(dim_s2)
    }

    denom = apply(cbind(matrix(x2,m2,1), s_aux22), 1, function(rowxs) calDenkdnorm(rowxs, s_aux22, hs, x2, hx, 0))
    num = (1-pis.hat2f) * dnorm(x2,muhat,sighat)
    clfdrnest = pmin(num/denom,1)
    clfdrnest[is.na(clfdrnest)] = 0
    clfdrnest[is.infinite(clfdrnest)] = 1

    denom = apply(cbind(matrix(x1,m1,1), s_aux21), 1, function(rowxs) calDenkdnorm(rowxs, s_aux22, hs, x2, hx, clfdrnest))
    rho = dnorm(x1,muhat,sighat) / denom
    rho[is.na(rho)] = 0
    rho[is.infinite(rho)] = 1e16
    rhow = rho / ws
    rhow[is.na(rhow)] = 0
    rhow[is.infinite(rhow)] = 1e16
    st.rhow = sort(rhow)

    matCQW = apply(cbind(matrix(pis.hat1,m1,1),matrix(ws,m1,1),s_aux21), 1, function(rowpiws) calDenkdnormSim(rowpiws, s_aux22, hs, x2, hx, clfdrnest, muhat, sighat, st.rhow))
    matCQW[is.na(matCQW)] = 1
    fdps = rowSums(matCQW) / (1 : m1)

    de = rep(0, m1)
    if (sum(fdps <= q) == 0)
        {
            k = 0
            pk = 1
        } else {
            k = max(which(fdps <= q))
            pk = st.rhow[k]
            de[which(rhow <= pk)] = 1
        }
        y<-list(nr=k, th=pk, de=de, rhow=rhow, fdp=fdps)
      
    return (y)
}


CLFDRmultif1 = function(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, muhat=0, sighat=1) {    
    m1 = length(x1)
    m2 = length(x2)
    dim_s2 = dim(s_aux21)[2]
    hx=density(x2,from = min(x2)-5,to=max(x2)+5,n=1000)$bw
    hs = matrix(Hns(s_aux22), dim_s2, dim_s2)
    if(sum(diag(hs)) == 0) {
        hs = diag(dim_s2)
    }
    ws = pis.hat1 / (1 - pis.hat1)

    denom = apply(cbind(matrix(x2,m2,1), s_aux22), 1, function(rowxs) calDenkdnorm(rowxs, s_aux22, hs, x2, hx, 0))
    num = (1-pis.hat2f) * dnorm(x2,muhat,sighat)
    clfdrnest = pmin(num/denom,1)
    clfdrnest[is.na(clfdrnest)] = 0
    clfdrnest[is.infinite(clfdrnest)] = 1

    denom = apply(cbind(matrix(x1,m1,1), s_aux21), 1, function(rowxs) calDenkdnorm(rowxs, s_aux22, hs, x2, hx, clfdrnest))
    rho = dnorm(x1,muhat,sighat) / denom

    rho[is.na(rho)] = 0
    rho[is.infinite(rho)] = 1e16
    rhow = rho / ws
    rhow[is.na(rhow)] = 0
    rhow[is.infinite(rhow)] = 1e16
    
    clfdrs = rhow / (1 + rhow)
    st.clfdrs = sort(clfdrs)
    fdps = cumsum(st.clfdrs) / (1 : m1)

    de = rep(0, m1)
    if(sum(fdps <= q) == 0)
    {
        k = 0
        pk = 1
    } else {
        k = max(which(fdps <= q))
        pk = st.clfdrs[k]
        de[which(clfdrs <= pk)] = 1
    }
    y = list(nr=k, th=pk, de=de, rhow=clfdrs, fdp=fdps)
      
    return (y)
}


#' A parallel computing version of pisLAWS utilizing foreach function.
#' 
#' A parallel computing version of pisLAWS utilizing foreach function.
#'
#' @param s1 the covariates to estimate a non-null probability for.
#' @param pv2 the p-values used for estimation.
#' @param s2 the covariates corresponding to the p-values pv2.
#' @param tau threshold for null p-values.
#' @param h bandwidth matrix H.
#'
#' @return Estimated non-null probability for each covariate in s1.
#' 
#' @export
pisLAWS.prl<- function(s1, pv2, s2, tau=0.1, h=0.3)
{
    m1 <- dim(s1)[1]
    m2 = length(pv2)
    dim_s1 = dim(s1)[2]
    p.est <-rep(0, m1)
    inds = which(pv2>=tau)

    forOut = foreach(j = 1:m1, .combine = 'rbind') %dopar% {
      c(j, calPest(s1[j,], s2, h, inds)/(1-tau))
    }
    orderfor = order(forOut[, 1])
    p.est = forOut[orderfor, 2]

    p.est[which(p.est>=1)] <- 1-1e-16
    p.est[which(p.est<=1e-16)] <-1e-16

    return(1-p.est)
}


rhoBHmultif1.prl = function(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, muhat=0, sighat=1) {    
    m1 = length(x1)
    m2 = length(x2)
    dim_s2 = dim(s_aux21)[2]
    ws = pis.hat1 / (1 - pis.hat1)

    hx = density(x2,from = min(x2)-5,to=max(x2)+5,n=1000)$bw
    hs = matrix(Hns(s_aux22), dim_s2, dim_s2)
    if(sum(diag(hs)) == 0) {
        hs = diag(dim_s2)
    }

    matxs = cbind(matrix(x2,m2,1), s_aux22)
    forOut = foreach(j = 1:m2, .combine = 'rbind') %dopar% {
        c(j, calDenkdnorm(matxs[j,], s_aux22, hs, x2, hx, 0))
    }
    orderfor = order(forOut[, 1])
    denom = forOut[orderfor, 2]

    num = (1-pis.hat2f) * dnorm(x2,muhat,sighat)
    clfdrnest = pmin(num/denom,1)
    clfdrnest[is.na(clfdrnest)] = 0
    clfdrnest[is.infinite(clfdrnest)] = 1

    matxs = cbind(matrix(x1,m1,1), s_aux21)
    forOut = foreach(j = 1:m1, .combine = 'rbind') %dopar% {
        c(j, calDenkdnorm(matxs[j,], s_aux22, hs, x2, hx, clfdrnest))
    }
    orderfor = order(forOut[, 1])
    denom = forOut[orderfor, 2] 

    rho = dnorm(x1,muhat,sighat) / denom
    rho[is.na(rho)] = 0
    rho[is.infinite(rho)] = 1e16
    rhow = rho / ws
    rhow[is.na(rhow)] = 0
    rhow[is.infinite(rhow)] = 1e16
    st.rhow = sort(rhow)

    matxs = cbind(matrix(pis.hat1,m1,1),matrix(ws,m1,1),s_aux21)
    forOut = foreach(j = 1:m1, .combine = 'rbind') %dopar% {
        c(j, calDenkdnormSim(matxs[j,], s_aux22, hs, x2, hx, clfdrnest, muhat, sighat, st.rhow))
    }
    orderfor = order(forOut[, 1])
    matCQW = forOut[orderfor, -1]
    matCQW[is.na(matCQW)] = 1
    fdps = colSums(matCQW) / (1 : m1)

    de = rep(0, m1)
    if (sum(fdps <= q) == 0)
        {
            k = 0
            pk = 1
        } else {
            k = max(which(fdps <= q))
            pk = st.rhow[k]
            de[which(rhow <= pk)] = 1
        }
        y<-list(nr=k, th=pk, de=de, rhow=rhow, fdp=fdps)
      
    return (y)
}








