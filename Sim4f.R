
#SIM 4: Sigmap Varies
library(e1071)

source('~/Dropbox/DLMNP_sims/sim_funcsf.r')

# "experiment" parameters
nsims=100
nstim=60
nvox=50
nsubject=20 # number subjects per simulation
#conditioning variable X
xp=rep(c(0,1),each=30)

#mean coefficients (gamma)
par.gamma0 = 0 
par.gammap = 1


#trial-level
par.etvs = 1


#subject-level
par.sigma0 = 0 
par.sigmap = seq(0,5,0.5)
nlevs=length(par.sigmap)
par.e0s = rnorm(nsubject*nsims*nlevs,0,par.sigma0) 
par.eps = rnorm(nsubject*nsims*nlevs,0,rep(par.sigmap,each=nsubject*nsims))

#voxel-level
par.tau0 = 0
par.taup = 0.5


par.e0vs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,par.tau0),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)
par.epvs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,par.taup),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)


ss_data = matrix(nrow=nsubject*nsims*nlevs,ncol=8)
ss_data = cbind(ss_data,rep(par.sigmap,each=nsubject*nsims),rep(1:nsubject,each=nlevs*nsims),rep(1:nsims,nlevs*nsubject))

	
for(subj in 1:(nsubject*nsims*nlevs)){	
	act_patterns=act_gen(xp, par.gamma0,par.gammap,par.e0s[subj],par.eps[subj],par.e0vs[,subj],par.epvs[,subj],par.etvs)	
	sim_mat=cor(act_patterns$full)
	mean_act=apply(act_patterns$full, 2, mean)
	temp_act = tapply(mean_act,xp,mean)
	
	temp_withinbetween = within_between_sim(sim_mat,xp)
	
	
	#SVM
	xp_short=xp[c(1:15,31:45)]
	
	train=data.frame(t(scale(act_patterns$full[,c(1:15,31:45)])))
	test=data.frame(t(scale(act_patterns$full[,c(16:30,46:60)])))
	
	svm1=svm(as.factor(xp_short)~.,data=train,type='nu-classification',nu=0.01,kernel='linear')
perf_svm1=sum(predict(svm1,test)==xp_short)/30

	
	ss_data[subj,1:8]=c(temp_withinbetween,temp_act,perf_svm1)
	
	
}

ss_data=data.frame(ss_data)
ss_data=cbind(ss_data,rep(1:nsims,nlevs*nsubject))

names(ss_data) = c('all','within','between','withinA','withinB','actA','actB','svm','lev','subj','sim')

ss_act_means=aggregate((actB-actA)~lev+sim,data=ss_data,FUN=mean)
ss_act_ts=aggregate((actB-actA)~lev+sim,data=ss_data,FUN=tstat)
ss_svm=aggregate(svm~lev+sim,data=ss_data,FUN=mean)
ss_cor_means=aggregate((withinB-withinA)~lev+sim,data=ss_data,FUN=mean)
ss_cor_ts=aggregate((withinB-withinA)~lev+sim,data=ss_data,FUN=tstat)

#####################
#Plot Mean Act (Fig 6A)
#####################
sim_act_means=aggregate(ss_act_means[,3]~lev,data=ss_act_means,FUN=mean)
sim_q1_1=aggregate(ss_act_means[,3]~lev,data=ss_act_means,FUN=q1)
sim_q3_1=aggregate(ss_act_means[,3]~lev,data=ss_act_means,FUN=q3)

plot(sim_act_means$lev,sim_act_means[,2],pch=20,ylim=c(-3,3),xlab="Sigmap",ylab="uB-uA")
uppers=sim_q3_1[,2]
lowers=sim_q1_1[,2]
segments(sim_act_means$lev,uppers,sim_act_means$lev,lowers,lwd=1.3)
segments(sim_act_means$lev-.07,uppers,sim_act_means$lev+.07,uppers)
segments(sim_act_means$lev-.07,lowers,sim_act_means$lev+.07,lowers)
#####################
#Plot Act Ts (Fig 6B)
#####################
sim_act_ts=aggregate(ss_act_ts[,3]~lev,data=ss_act_ts,FUN=mean)
sim_q1_ts=aggregate(ss_act_ts[,3]~lev,data=ss_act_ts,FUN=q1)
sim_q3_ts=aggregate(ss_act_ts[,3]~lev,data=ss_act_ts,FUN=q3)

plot(sim_act_ts$lev,sim_act_ts[,2],pch=20,ylim=c(-0.5,70),xlab="Sigmap",ylab="t-statistic")
uppers=sim_q3_ts[,2]
lowers=sim_q1_ts[,2]
segments(sim_act_ts$lev,uppers,sim_act_ts$lev,lowers,lwd=1.3)
segments(sim_act_ts$lev-.07,uppers,sim_act_ts$lev+.07,uppers)
segments(sim_act_ts$lev-.07,lowers,sim_act_ts$lev+.07,lowers)
abline(h=qt(0.975,19),lty=2)

#####################
#Plot Cor Means (6C)
#####################

sim_cor_means=aggregate(ss_cor_means[,3]~lev,data=ss_cor_means,FUN=mean)
sim_q1_1=aggregate(ss_cor_means[,3]~lev,data=ss_cor_means,FUN=q1)
sim_q3_1=aggregate(ss_cor_means[,3]~lev,data=ss_cor_means,FUN=q3)

plot(sim_cor_means$lev,sim_cor_means[,2],pch=20,ylim=c(0,1),xlab="Sigmap",ylab="Pearson's r within B - within A")
uppers=sim_q3_1[,2]
lowers=sim_q1_1[,2]
segments(sim_cor_means$lev,uppers,sim_cor_means$lev,lowers,lwd=1.3)
segments(sim_cor_means$lev-.07,uppers,sim_cor_means$lev+.07,uppers)
segments(sim_cor_means$lev-.07,lowers,sim_cor_means$lev+.07,lowers)

#####################
#Plot Cor Ts (6D)
#####################
sim_cor_ts=aggregate(ss_cor_ts[,3]~lev,data=ss_cor_ts,FUN=mean)
sim_q1_1=aggregate(ss_cor_ts[,3]~lev,data=ss_cor_ts,FUN=q1)
sim_q3_1=aggregate(ss_cor_ts[,3]~lev,data=ss_cor_ts,FUN=q3)

plot(sim_cor_ts$lev,sim_cor_ts[,2],pch=20,ylim=c(-0.5,70),xlab="Tp",ylab="t-stat: within B - within A")
uppers=sim_q3_1[,2]
lowers=sim_q1_1[,2]
segments(sim_cor_ts$lev,uppers,sim_cor_ts$lev,lowers,lwd=1.3)
segments(sim_cor_ts$lev-.07,uppers,sim_cor_ts$lev+.07,uppers)
segments(sim_cor_ts$lev-.07,lowers,sim_cor_ts$lev+.07,lowers)


#####################
#Plot SVM (6E)
#####################

sim_svm_means=aggregate(ss_svm[,3]~lev,data=ss_svm,FUN=mean)
sim_q1_svm=aggregate(ss_svm[,3]~lev,data=ss_svm,FUN=q1)
sim_q3_svm=aggregate(ss_svm[,3]~lev,data=ss_svm,FUN=q3)

plot(sim_svm_means$lev,sim_svm_means[,2],pch=20,ylim=c(0,1),xlab="Tp",ylab="Classifier Accuracy")
uppers=sim_q3_svm[,2]
lowers=sim_q1_svm[,2]
segments(sim_svm_means$lev,uppers,sim_svm_means$lev,lowers,lwd=1.3)
segments(sim_svm_means$lev-.07,uppers,sim_svm_means$lev+.07,uppers)
segments(sim_svm_means$lev-.07,lowers,sim_svm_means$lev+.07,lowers)

