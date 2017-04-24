library(mvtnorm)

materncov=function(d,b){
  r=(1+sqrt(5)*d/b+5*d^2/(3*b^2))*exp(-sqrt(5)*d/b)
  return=r
}

gibbssampler=function(y,gene,gr,ti,Cinv,t){
  gene_number=max(gene)
  group_number=max(gr)
  time_span=max(ti)
  total_data=length(y)
  gnenelist=rep(0,gene_number)
  group_size=rep(0,group_number)
  for(i in 1:gene_number){
    genelist[gene[i]]=gr[i]
    group_size[gene[i]]=group_size[gene[i]]+1
  }
  sigma_eps=1
  sigma_tau=rep(1,gene_number)
  sigma_h=rep(1,gene_number)
  sigma_f=rep(1,group_number)
  f=matrix(0,nrow=time_span,ncol=group_number)
  h=matrix(0,nrow=time_span,ncol=gene_number)
  for(i in 1:(1000+t)){
    hmean=matrix(0,nrow=time_span,ncol=gene_number)
    for(j in 1:total_data){
      hmean[ti[j],gene[j]]=hmean[ti[j],gene[j]]+y[j]-f[ti[j],gr[j]]
    }
    for(j in 1:gene_number){
      effsigma=(sigma_eps+sigma_tau[j])/(sigma_eps*sigma_tau[j])
      Vh=solve(3*diag(effsigma,time_span)+Cinv/sigma_h)
      hmean[,j]=effsigma*(Vh%*%hmean[,j])
      h[,j]=rmvnorm(1,mean=hmean[,j],sigma=Vh)
    }
    fmean=matrix(0,nrow=time_span,ncol=group_number)
    for(j in 1:total_data){
      fmean[ti[j],gr[j]]=fmean[ti[j],gr[j]]+y[j]-h[ti[j],gene[j]]
    }
    for(j in 1:gr_number){
      effsigma=
        Vf=solve(3*diag(effsigma,time_span)+Cinv/sigma_f)
      fmean[,j]=effsigma*(Vf%*%fmean[,j])
      f[,j]=rmvnorm(1,mean=fmean[,j],sigma=Vf)
    }
    rms=rep(0,gene_number)
    for(j in 1:total_data){
      dif=y[j]-h[ti[j],gene[j]]-f[ti[j],gr[j]]
      rms[gene[j]]=rms[gene[j]]+crossprod(dif,dif)
    }
    sigma_eps=rgamma(1,shape=(2+total_data)/2,rate=sum(rms)/2)
    sigma_eps=1/sigma_eps
    for(j in 1:gene_number){
      sigma_tau[j]=rgamma(1,shape=2.5,rate=rms[j]/2)
      sigma_h[j]=rgamma(1,shape=2.5,rate=t(h[,j])%*%Cinv%*%h[,j]/2)
    }
    sigma_tau=1/sigma_tau
    sigma_h=1/sigma_h
    for(j in 1:group_number){
      sigma_f[j]=rgamma(1,shape=(2+group_size[j])/2,rate=t(h[,j])%*%Cinv%*%h[,j]/2)
    }
    sigma_f=1/sigma_f
  }
}