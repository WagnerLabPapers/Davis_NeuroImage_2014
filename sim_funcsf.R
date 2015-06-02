act_gen = function(xp,gamma0,gammap,e0s,eps,e0vs,epvs,etvs){
	
	nstim = length(xp)
	nvox = length(e0vs)
	error = matrix(rnorm(nstim*nvox,0,etvs),nrow=nstim,ncol=nvox)
	
		acttvs=gamma0+e0s+outer(rep(1,nstim),e0vs)+gammap*xp+eps*xp+outer(xp,epvs)+error
	
	
		
	return(list(full = t(acttvs)))
}

within_between_sim=function(cor_matrix,cat){
	
	nstim=length(cor_matrix[1,])
	cclass=rep(cat, each =nstim)
	rclass=rep(cat, nstim)
	flat_vec=c(cor_matrix)
	identity = c(diag(nstim)) 
	all=mean(flat_vec[identity!=1])
	between=mean(flat_vec[rclass!=cclass&identity!=1])
	within=mean(flat_vec[rclass==cclass&identity!=1])
	singlewithin = tapply(flat_vec[rclass==cclass&identity!=1],rclass[rclass==cclass&identity!=1],mean)

	return(c(all, within,between,singlewithin))
}

tstat=function(x){return(t.test(x)$statistic)}
q1=function(x){return(quantile(x,0.25))}
q3=function(x){return(quantile(x,0.75))}
