#' Implementing the oracle rho-BH procedure.
#' 
#' This is a simple implementation where null and alaternative distributions are all assumed to be Gaussian.
#' Null distribution is fixed as N(0,1).
#' Alternative distribution is N(mu1_i, sigma_i), where mu1_i and sigma_i refer to each element of the input mu1 and sigma.
#' Rho-BH with other distributions can be implemented through easy modification of this function.
#'
#' @param x the test statistic, whose null distribution is fixed to be N(0,1).
#' @param mu1 mean vector for the alternative Gaussian distribution.
#' @param sigma sigma vector for the alternative Gaussian distribution.
#' @param ws weight vector.
#' @param q desired FDR level.
#' @param pis non-null probability vector.
#'
#' @return A list with the elements
#' \item{nr}{number of rejections}
#' \item{th}{the threshold for rejection of the procedure}
#' \item{de}{the rejections}
#' \item{rhow}{q_i, weighted rho-values}
#' \item{fdps}{FDP estimate at each rejection threshold}
#' 
#' @export

rhoBH.OR = function(x, mu1, sigma, ws=NULL, q, pis=NULL) {
    if(is.null(pis) & is.null(ws)) {
      print(paste0('Please specify either ', expression(pi[i]), ' or ', expression(w[i]), ' for the oracle procedure.'))
      return(0)
    }
    if(is.null(pis)) {
        pis = ws / (1 + ws)
    }
    if(is.null(ws)) {
        ws = pis / (1 - pis)
    }

    num = length(x)
    rhow = dnorm(x, 0, 1) / dnorm(x, mu1, sigma) / ws
    st.rhow = sort(rhow)
    matCQW = matrix(0, num, num)
    for(i in 1:num) {
            x_sim = rnorm(2000, 0, 1)
            rho_sim = dnorm(x_sim, 0, 1) / dnorm(x_sim, mu1[i], sigma[i])
            ecdfrho = ecdf(rho_sim)
            matCQW[i, ] = (1 - pis[i]) * ecdfrho(st.rhow * ws[i])
        }
    matCQW[is.na(matCQW)] = 1
    fdps = colSums(matCQW) / (1 : num)

    de = rep(0, num)
    if (sum(fdps <= q) == 0)
        {
            k = 0
            pk = 1
        } else {
            k = max(which(fdps <= q))
            pk = st.rhow[k]
            de[which(rhow <= pk)] = 1
        }
        y = list(nr=k, th=pk, de=de, rhow=rhow, fdp=fdps)
    
    return (y)
}


#' Implementing the data-driven rho-BH procedure.
#' 
#' This is a simple implementation where null distributions are all assumed to be Gaussian.
#' Rho-BH with other null distributions can be implemented through easy modification.
#'
#' @param x the test statistic, whose null distribution is Gaussian.
#' @param s_aux1 a m*dim1 matrix for estimation of non-null probability.
#' @param s_aux2 a m*dim2 matrix for estimation of alternative distributions.
#' @param q desired FDR level.
#' @param mu0 mean vector for the null Gaussian distribution.
#' @param sig0 sigma vector for the null Gaussian distribution.
#'
#' @return A list with the elements
#' \item{nr1}{number of rejections, for the first half of data}
#' \item{th1}{the threshold for rejection of the procedure, for the first half of data}
#' \item{nr2}{number of rejections, for the second half of data}
#' \item{th2}{the threshold for rejection of the procedure, for the second half of data}
#' \item{de}{the rejections fot the whole procedure}
#' \item{rdsp1}{the index for the first half of data}
#' \item{rhow1}{q_i, weighted rho-values, for the first half of data}
#' \item{fdps1}{FDP estimate at each rejection threshold, for the first half of data}
#' \item{rdsp2}{the index for the second half of data}
#' \item{rhow2}{q_i, weighted rho-values, for the second half of data}
#' \item{fdps2}{FDP estimate at each rejection threshold, for the second half of data}
#' 
#' @export
#' @importFrom mvnfast dmvn

rhoBH = function(x, s_aux1, s_aux2, q, mu0=0, sig0=1) {
    m = length(x)
    dim_s1 = dim(s_aux1)[2]
    dim_s2 = dim(s_aux2)[2]
    rdsp1 = sample(1:m, m/2)
    rdsp2 = setdiff(1:m, rdsp1)
    m1 = length(rdsp1)
    m2 = length(rdsp2)
    x1 = x[rdsp1]
    x2 = x[rdsp2]
    s_aux11 = matrix(s_aux1[rdsp1,],m1,dim_s1)
    s_aux12 = matrix(s_aux1[rdsp2,],m2,dim_s1)
    s_aux21 = matrix(s_aux2[rdsp1,],m1,dim_s2)
    s_aux22 = matrix(s_aux2[rdsp2,],m2,dim_s2)
    pvx1 = 2 * pnorm(-abs(x1), 0, 1)
    pvx2 = 2 * pnorm(-abs(x2), 0, 1)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
        h = diag(dim_s1)
    }
    pis.hat1<-pisLAWS(s1=s_aux11, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2<-pisLAWS(s1=s_aux12, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat1f<-pisLAWS(s1=s_aux11, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2f<-pisLAWS(s1=s_aux12, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    resTemp1 = rhoBHmultif1(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, mu0, sig0)
    resTemp2 = rhoBHmultif1(x2, x1, pis.hat2, pis.hat1, pis.hat2f, pis.hat1f, s_aux22, s_aux21, q, mu0, sig0)
    rbh.dd.de = rep(0, m)
    rbh.dd.de[rdsp1] = resTemp1$de
    rbh.dd.de[rdsp2] = resTemp2$de
    y = list(nr1=resTemp1$nr, th1=resTemp1$th, nr2=resTemp2$nr, th2=resTemp2$th, de=rbh.dd.de, rdsp1=rdsp1, rhow1=resTemp1$rhow, fdp1=resTemp1$fdp, rdsp2=rdsp2, rhow2=resTemp2$rhow, fdp2=resTemp2$fdp, pis.hat1=pis.hat1, pis.hat2=pis.hat2)
    
    return (y)
}


#' Implementing the procedure based on conditional local FDR.
#' 
#' This is a simple implementation where null distributions are all assumed to be Gaussian.
#' The procedure with other null distributions can be implemented through easy modification.
#'
#' @param x the test statistic, whose null distribution is Gaussian.
#' @param s_aux1 a m*dim1 matrix for estimation of non-null probability.
#' @param s_aux2 a m*dim2 matrix for estimation of alternative distributions.
#' @param q desired FDR level.
#' @param mu0 mean vector for the null Gaussian distribution.
#' @param sig0 sigma vector for the null Gaussian distribution.
#'
#' @return A list with the elements
#' \item{nr1}{number of rejections, for the first half of data}
#' \item{th1}{the threshold for rejection of the procedure, for the first half of data}
#' \item{nr2}{number of rejections, for the second half of data}
#' \item{th2}{the threshold for rejection of the procedure, for the second half of data}
#' \item{de}{the rejections fot the whole procedure}
#' \item{rdsp1}{the index for the first half of data}
#' \item{rhow1}{q_i, weighted rho-values, for the first half of data}
#' \item{fdps1}{FDP estimate at each rejection threshold, for the first half of data}
#' \item{rdsp2}{the index for the second half of data}
#' \item{rhow2}{q_i, weighted rho-values, for the second half of data}
#' \item{fdps2}{FDP estimate at each rejection threshold, for the second half of data}
#' 
#' @export
#' @importFrom mvnfast dmvn

CLFDR.DD = function(x, s_aux1, s_aux2, q, mu0=0, sig0=1) {
    m = length(x)
    dim_s1 = dim(s_aux1)[2]
    dim_s2 = dim(s_aux2)[2]
    rdsp1 = sample(1:m, m/2)
    rdsp2 = setdiff(1:m, rdsp1)
    m1 = length(rdsp1)
    m2 = length(rdsp2)
    x1 = x[rdsp1]
    x2 = x[rdsp2]
    s_aux11 = matrix(s_aux1[rdsp1,],m1,dim_s1)
    s_aux12 = matrix(s_aux1[rdsp2,],m2,dim_s1)
    s_aux21 = matrix(s_aux2[rdsp1,],m1,dim_s2)
    s_aux22 = matrix(s_aux2[rdsp2,],m2,dim_s2)
    pvx1 = 2 * pnorm(-abs(x1), 0, 1)
    pvx2 = 2 * pnorm(-abs(x2), 0, 1)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat1<-pisLAWS(s1=s_aux11, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2<-pisLAWS(s1=s_aux12, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat1f<-pisLAWS(s1=s_aux11, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2f<-pisLAWS(s1=s_aux12, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    resTemp1 = CLFDRmultif1(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, mu0, sig0)
    resTemp2 = CLFDRmultif1(x2, x1, pis.hat2, pis.hat1, pis.hat2f, pis.hat1f, s_aux22, s_aux21, q, mu0, sig0)
    rbh.dd.de = rep(0, m)
    rbh.dd.de[rdsp1] = resTemp1$de
    rbh.dd.de[rdsp2] = resTemp2$de
    y = list(nr1=resTemp1$nr, th1=resTemp1$th, nr2=resTemp2$nr, th2=resTemp2$th, de=rbh.dd.de, rdsp1=rdsp1, rhow1=resTemp1$rhow, fdp1=resTemp1$fdp, rdsp2=rdsp2, rhow2=resTemp2$rhow, fdp2=resTemp2$fdp, pis.hat1=pis.hat1, pis.hat2=pis.hat2)
    
    return (y)
}



#' A parallel computing version of rho.BH utilizing foreach function.
#' 
#' This is a simple implementation where null distributions are all assumed to be Gaussian.
#' Rho-BH with other null distributions can be implemented through easy modification.
#'
#' @param x the test statistic, whose null distribution is Gaussian.
#' @param s_aux1 a m*dim1 matrix for estimation of non-null probability.
#' @param s_aux2 a m*dim2 matrix for estimation of alternative distributions.
#' @param q desired FDR level.
#' @param mu0 mean vector for the null Gaussian distribution.
#' @param sig0 sigma vector for the null Gaussian distribution.
#'
#' @return A list with the elements
#' \item{nr1}{number of rejections, for the first half of data}
#' \item{th1}{the threshold for rejection of the procedure, for the first half of data}
#' \item{nr2}{number of rejections, for the second half of data}
#' \item{th2}{the threshold for rejection of the procedure, for the second half of data}
#' \item{de}{the rejections fot the whole procedure}
#' \item{rdsp1}{the index for the first half of data}
#' \item{rhow1}{q_i, weighted rho-values, for the first half of data}
#' \item{fdps1}{FDP estimate at each rejection threshold, for the first half of data}
#' \item{rdsp2}{the index for the second half of data}
#' \item{rhow2}{q_i, weighted rho-values, for the second half of data}
#' \item{fdps2}{FDP estimate at each rejection threshold, for the second half of data}
#' 
#' @export
#' @import foreach
#' @import parallel
#' @import doParallel
#' @importFrom mvnfast dmvn
#' @importFrom foreach foreach

rhoBH.prl = function(x, s_aux1, s_aux2, q, mu0=0, sig0=1) {
    m = length(x)
    dim_s1 = dim(s_aux1)[2]
    dim_s2 = dim(s_aux2)[2]
    rdsp1 = sample(1:m, m/2)
    rdsp2 = setdiff(1:m, rdsp1)
    m1 = length(rdsp1)
    m2 = length(rdsp2)
    x1 = x[rdsp1]
    x2 = x[rdsp2]
    s_aux11 = matrix(s_aux1[rdsp1,],m1,dim_s1)
    s_aux12 = matrix(s_aux1[rdsp2,],m2,dim_s1)
    s_aux21 = matrix(s_aux2[rdsp1,],m1,dim_s2)
    s_aux22 = matrix(s_aux2[rdsp2,],m2,dim_s2)
    pvx1 = 2 * pnorm(-abs(x1), 0, 1)
    pvx2 = 2 * pnorm(-abs(x2), 0, 1)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
        h = diag(dim_s1)
    }
    pis.hat1<-pisLAWS.prl(s1=s_aux11, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2<-pisLAWS.prl(s1=s_aux12, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx1, 0.9)$th
    h = matrix(Hns(s_aux11), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat1f<-pisLAWS.prl(s1=s_aux11, pv2=pvx1, s2=s_aux11, tau=bh.th, h=h)

    bh.th<-bh.func(pvx2, 0.9)$th
    h = matrix(Hns(s_aux12), dim_s1, dim_s1)
    if(sum(diag(h)) == 0) {
      h = diag(dim_s1)
    }
    pis.hat2f<-pisLAWS.prl(s1=s_aux12, pv2=pvx2, s2=s_aux12, tau=bh.th, h=h)

    resTemp1 = rhoBHmultif1.prl(x1, x2, pis.hat1, pis.hat2, pis.hat1f, pis.hat2f, s_aux21, s_aux22, q, mu0, sig0)
    resTemp2 = rhoBHmultif1.prl(x2, x1, pis.hat2, pis.hat1, pis.hat2f, pis.hat1f, s_aux22, s_aux21, q, mu0, sig0)
    rbh.dd.de = rep(0, m)
    rbh.dd.de[rdsp1] = resTemp1$de
    rbh.dd.de[rdsp2] = resTemp2$de
    y = list(nr1=resTemp1$nr, th1=resTemp1$th, nr2=resTemp2$nr, th2=resTemp2$th, de=rbh.dd.de, rdsp1=rdsp1, rhow1=resTemp1$rhow, fdp1=resTemp1$fdp, rdsp2=rdsp2, rhow2=resTemp2$rhow, fdp2=resTemp2$fdp, pis.hat1=pis.hat1, pis.hat2=pis.hat2)
    
    return (y)
}


