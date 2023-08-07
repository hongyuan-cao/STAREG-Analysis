# ----------- Functions for performing the MaRR procedure ------------------------

# get.SS: a function to produce a vector of sums of squared differences between observed
#         and actual survival functions for k-hat=0,...,n-1
# INPUTS: max.rank = a vector of maximum rank statistics
# OUTPUTS: mySS = a vector of SS(i/k) values as described in manuscript
get.SS = function(max.rank)
{
  n = length(max.rank)
  x = 1:n/n
  mySS=rep(NA,n)
  my.W = sapply(1:n,function(i){sum(max.rank==i)})
  surv.function = ( n-cumsum(my.W))/n
  for(k in 0:(n-1))
  {
    i=k+1; pi1=k/n; pi0=1-pi1
    temp.W =surv.function[i:n]
    Sn = pi0*(1-(x-pi1)^2/pi0^2)[i:n]
    temp.diff = temp.W-Sn
    sq.diff = (temp.diff*temp.diff)/n
    mySS[i] = sum(sq.diff)/pi0
  }
  return(mySS)
}

# est.fdr: a function to produce a vector of estimated false discovery rats based on 
#          a vector of maximum rank statistics and a khat.
# INPUTS: max.rank = a vector of maximum rank statistics
#         khat = a value of khat calculated from the argmin of mySS
# OUTPUTS: a vector of estimated FDR values for each potential threshold Nhat=1,...,n
est.fdr = function(khat,max.rank)
{
  n=length(max.rank)
  Nhat = (khat+1):n
  Q.khat = sum(max.rank<=khat)
  R.Nhat = sapply(Nhat,function(i){sum(max.rank<=i)})
  indicate.R = as.numeric(R.Nhat>0)
  R.Nhat[indicate.R==0] = 1
  temp = Nhat-khat
  numer=temp*temp
  denom = (n-khat)*R.Nhat
  FDR.Nhatkhat = (numer/denom)*indicate.R
  return(c(rep(0,khat),FDR.Nhatkhat))
}

# MaRR: a function that performs the MaRR procedure based on a vector of maximum rank statistics.
# INPUTS: max.rank = a vector of maximum rank statistics
#         cutoff = a value between 0 and 1 that provides the maximum allowed value for pi-hat.
#         alpha = desired level of FDR control
#         khat.to.zero = T/F, whether or not to set k-hat to zero for mFDR calculation (recommended for very small pi1)
# OUTPUTS: khat = n*pi-hat, discrete estimate of where irreproducible signals begin
#          Nhat = estimated cut-off for maximum ranks that will control FDR at level alpha
#          est.fdr = the estimated fdr value for each of potential N-hat
#          SS = vector of values for SS loss function evaluated at i/n =0, 1/n, 2/n,..., 1
#          which.sig = vector of indices for maximum ranks declared to be reproducible
MaRR = function(max.rank,cutoff=.9,alpha=.05,khat.to.zero=F){
  maxx=floor(cutoff*length(max.rank))
  mySS = get.SS(max.rank)
  khat = which(mySS[1:maxx]==min(mySS[1:maxx]))-1
  if(khat.to.zero==T){khat=0}
  temp.fdr = est.fdr(khat,max.rank)
  Nhat = max(which(temp.fdr<=alpha))
  which.sig=which(max.rank<=Nhat)
  return(list(Nhat=Nhat,khat=khat,est.fdr=temp.fdr,SS=mySS,which.sig=which.sig))
}
