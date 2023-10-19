#' Implementing BH procedure.
#' 
#' Implementing BH procedure.
#'
#' @param pv p-values.
#' @param q desired FDR level.
#' 
#' @return A list with the elements
#' \item{nr}{number of rejections}
#' \item{th}{the threshold for rejection of the procedure}
#' \item{de}{the rejections}
#' 
#' @export
bh.func<-function(pv, q)
{ 
  # the input 
    # pv: the p-values
    # q: the FDR level
  # the output 
    # nr: the number of hypothesis to be rejected
    # th: the p-value threshold
    # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)   
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}


lin.itp<-function(x, X, Y){
  ## x: the coordinates of points where the density needs to be interpolated
  ## X: the coordinates of the estimated densities
  ## Y: the values of the estimated densities
  ## the output is the interpolated densities
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}


EstNull.func<-function (x,gamma=0.1)
{
 # x is a vector of z-values
 # gamma is a parameter, default is 0.1
 # output the estimated mean and standard deviation

 n = length(x)
 t = c(1:1000)/200
 
 gan    = n^(-gamma)
 that   = 0 
 shat   = 0
 uhat   = 0
 epshat = 0

 phiplus   = rep(1,1000)
 phiminus  = rep(1,1000)
 dphiplus  = rep(1,1000)
 dphiminus = rep(1,1000)
 phi       = rep(1,1000)
 dphi      = rep(1,1000)

 for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
 }

 ind = min(c(1:1000)[(phi - gan) <= 0])
 tt = t[ind]
 a  = phiplus[ind]
 b  = phiminus[ind]
 da = dphiplus[ind]
 db = dphiminus[ind]
 c  = phi[ind]

 that   = tt
 shat   = -(a*da + b*db)/(tt*c*c)
 shat   = sqrt(shat) 
 uhat   = -(da*b - db*a)/(c*c)
 epshat = 1 - c*exp((tt*shat)^2/2)

 return(musigma=list(mu=uhat,s=shat))
}


#' Implementing LAWS procedure.
#' 
#' Implementing LAWS procedure.
#'
#' @param pvs p-values.
#' @param pis estimation of non-null propability pi.
#' @param q desired FDR level.
#' 
#' @return A list with the elements
#' \item{nr}{number of rejections}
#' \item{th}{the threshold for rejection of the procedure}
#' \item{de}{the rejections}
#' 
#' @export
law.func<-function(pvs, pis, q)
{
    ## implementing "spatial multiple testing by locally adaptive weighting"
  ## Arguments
   # pvs: p-values
   # pis: conditional probabilities
   # q: FDR level
  ## Values
   # de: the decision
   # th: the threshold for weighted p-values
   
  m<-length(pvs)
  nu<-10e-5
  pis[which(pis<nu)]<-nu # stabilization
  pis[which(pis>1-nu)]<-1-nu # stabilization
  ws<-pis/(1-pis)
    pws<-pvs/ws
    st.pws<-sort(pws)
    fdps<-sum(pis)*st.pws/(1:m)
    de<-rep(0, m)
    if(sum(fdps<=q)==0)
    {
      k<-0
      pwk<-1
    }
    else
    {
      k<-max(which(fdps<=q))
      pwk<-st.pws[k]
      de[which(pws<=pwk)]<-1
    }
    y<-list(nr=k, th=pwk, de=de)
  return (y)
}  


epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  # estimate alternative proportion
  # this implemets the estimator in Jin and Cai
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0, tmax, 0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) { 
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}


