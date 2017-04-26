library(mvtnorm)

materncov=function(d,b){
  r=(1+sqrt(5)*d/b+5*d^2/(3*b^2))*exp(-sqrt(5)*d/b)
  return=r
}

gibbssampler=function(y,gene,gr,ti,C,t){
  Cinv=solve(C)
  gene_number=max(gene)
  group_number=max(gr)
  time_span=max(ti)
  total_data=length(y)
  genelist=rep(0,gene_number)
  group_size=rep(0,group_number)
  for(i in 1:gene_number){
    genelist[gene[i]]=gr[i]
    group_size[gr[i]]=group_size[gr[i]]+1
  }
  sigma_eps=1
  sigma_tau=rep(1,gene_number)
  sigma_h=rep(1,gene_number)
  sigma_f=rep(1,group_number)
  f=matrix(mean(y),nrow=time_span,ncol=group_number)
  mf=rowMeans(f)
  h=matrix(0,nrow=time_span,ncol=gene_number)
  mh=rowMeans(h)
  hsample=array(dim=c(t,time_span,gene_number))
  fsample=array(dim=c(t,time_span,group_number))
  sigma_eps_sample=rep(0,t)
  sigma_f_sample=matrix(nrow=group_number,ncol=t)
  sigma_h_sample=matrix(nrow=gene_number,ncol=t)
  sigma_tau_sample=matrix(nrow=gene_number,ncol=t)
  for(i in 1:(1000+t)){
    hmean=matrix(0,nrow=time_span,ncol=gene_number)
    for(j in 1:total_data){
      hmean[ti[j],gene[j]]=hmean[ti[j],gene[j]]+y[j]-f[ti[j],gr[j]]
    }
    for(j in 1:gene_number){
      effsigma=(sigma_eps+sigma_tau[j])/(sigma_eps*sigma_tau[j])
      Vh=solve(3*diag(effsigma,time_span)+Cinv/sigma_h[j])
      hmean[,j]=Vh%*%(effsigma*hmean[,j]+Cinv%*%mh/sigma_h[j])
      h[,j]=rmvnorm(1,mean=hmean[,j],sigma=Vh)
    }
    mh=as.vector(rmvnorm(1,mean=rowMeans(h),sum(1/sigma_h)*C/14))
    fmean=matrix(0,nrow=time_span,ncol=group_number)
    for(j in 1:total_data){
      fmean[ti[j],gr[j]]=fmean[ti[j],gr[j]]+(y[j]-h[ti[j],gene[j]])*(1/sigma_eps+1/sigma_tau[gene[j]])
    }
    for(j in 1:group_number){
      effsigma=group_size[j]/sigma_eps
      for(l in 1:gene_number){
        if(genelist[l]==j)effsigma=effsigma+1/sigma_tau[l]
      }
        Vf=solve(3*diag(effsigma,time_span)+Cinv/sigma_f[j])
      fmean[,j]=Vf%*%(fmean[,j]+Cinv%*%mf/sigma_f[j])
      f[,j]=rmvnorm(1,mean=fmean[,j],sigma=Vf)
    }
    mf=as.vector(rmvnorm(1,mean=rowMeans(f),sum(1/sigma_f)*C/3))
    rms=rep(0,gene_number)
    for(j in 1:total_data){
      dif=y[j]-h[ti[j],gene[j]]-f[ti[j],gr[j]]
      rms[gene[j]]=rms[gene[j]]+dif^2
    }
    sigma_eps=rgamma(1,shape=2+total_data/2,rate=sum(rms)/2)
    sigma_eps=1/sigma_eps
    for(j in 1:gene_number){
      sigma_tau[j]=rgamma(1,shape=3.5,rate=rms[j]/2)
      dh=h[,j]-mh
      sigma_h[j]=rgamma(1,shape=3.5,rate=crossprod(dh,Cinv)%*%dh/2)
    }
    sigma_tau=1/sigma_tau
    sigma_h=1/sigma_h
    for(j in 1:group_number){
      df=f[,j]-mf
      sigma_f[j]=rgamma(1,shape=2+group_size[j]/2,rate=crossprod(df,Cinv)%*%df/2)
    }
    sigma_f=1/sigma_f
    if(i>1000){
      hsample[i-1000,,]=h
      fsample[i-1000,,]=f
      sigma_f_sample[,i-1000]=sigma_f
      sigma_h_sample[,i-1000]=sigma_h
      sigma_tau_sample[,i-1000]=sigma_tau
      sigma_eps_sample[i-1000]=sigma_eps
    }
  }
  return=list(h=hsample,f=fsample,sf=sigma_f_sample,sh=sigma_h_sample,st=sigma_tau_sample,se=sigma_eps_sample)
}