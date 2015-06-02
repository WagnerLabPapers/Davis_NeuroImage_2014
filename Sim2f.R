
#SIM 2: Just taup varies; Dummy coded conditioning variable
library(e1071)

source('~/Dropbox/DLMNP_sims/sim_funcsf.r')

# "experiment" parameters
nsims=100
nstim=60
nvox=50
nsubject=1 # number subjects per simulation
#conditioning variable X
xp=rep(c(0,1),each=nstim/2)

#mean coefficients (gamma)
par.gamma0 = 0 
par.gammap = 0


#trial-level
par.etvs = 1

#voxel-level
par.tau0 = 0
par.taup = seq(0,5,0.5)
nlevs=length(par.taup)

par.e0vs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,rep(par.tau0,each=nsubject*nsims)),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)
par.epvs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,rep(par.taup,each=nsubject*nsims)),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)

#subject-level
par.sigma0 = 0 
par.sigmap = 0
par.e0s = rnorm(nsubject*nsims*nlevs,0,par.sigma0) 
par.eps = rnorm(nsubject*nsims*nlevs,0,par.sigmap)


ss_data = matrix(nrow=nsubject*nsims*nlevs,ncol=8)
ss_data = cbind(ss_data,rep(par.taup,each=nsubject*nsims),rep(1:nsubject,each=nlevs*nsims))

	
for(subj in 1:(nsubject*nsims*nlevs)){	
	act_patterns=act_gen(xp, par.gamma0,par.gammap,par.e0s[subj],par.eps[subj],par.e0vs[,subj],par.epvs[,subj],par.etvs)	
	sim_mat=cor(act_patterns$full)
	mean_act=apply(act_patterns$full, 2, mean)
	
	temp_withinbetween = within_between_sim(sim_mat,xp)
	
	temp_act = tapply(mean_act,xp,mean)
	
	#SVM
	xp_short=xp[c(1:15,31:45)]
	
	train=data.frame(t(act_patterns$full[,c(1:15,31:45)]))
	test=data.frame(t(act_patterns$full[,c(16:30,46:60)]))
	
	svm1=svm(as.factor(xp_short)~.,data=train,type='nu-classification',nu=0.01,kernel='linear')
	perf_svm1=sum(predict(svm1,test)==xp_short)/30
		
	ss_data[subj,1:8]=c(temp_withinbetween,temp_act,perf_svm1)
	
	
}

ss_data=data.frame(ss_data)


names(ss_data) = c('all','within','between','withinA','withinB','actA','actB','svm','lev','subj')

ss_means_a=aggregate(withinA~lev+subj,data=ss_data,FUN=mean)
ss_q1_a=aggregate(withinA~lev+subj,data=ss_data,FUN=q1)
ss_q3_a=aggregate(withinA~lev+subj,data=ss_data,FUN=q3)

ss_means_b=aggregate(withinB~lev+subj,data=ss_data,FUN=mean)
ss_q1_b=aggregate(withinB~lev+subj,data=ss_data,FUN=q1)
ss_q3_b=aggregate(withinB~lev+subj,data=ss_data,FUN=q3)

ss_means_bet=aggregate(between~lev+subj,data=ss_data,FUN=mean)
ss_q1_bet=aggregate(between~lev+subj,data=ss_data,FUN=q1)
ss_q3_bet=aggregate(between~lev+subj,data=ss_data,FUN=q3)

#####
#Plot Mean Correlations (Figure 5B)
#####

plot(ss_means_a$lev,ss_means_a$withinA,pch=20,ylim=c(0,1),xlab="tau_p",ylab="Pearson's r")

uppers=ss_q3_a$withinA
lowers=ss_q1_a$withinA
segments(ss_means_a$lev,uppers,ss_means_a$lev,lowers,lwd=1.3)
segments(ss_means_a$lev-.07,uppers,ss_means_a$lev+.07,uppers)
segments(ss_means_a$lev-.07,lowers,ss_means_a$lev+.07,lowers)

points(ss_means_b$lev,ss_means_b$withinB,pch=18)

lines(ss_means_b$lev,ss_means_b$lev^2/(par.e+ss_means_b$lev^2))
lines(ss_means_b$lev,0*ss_means_b$lev^2/(par.e+ss_means_b$lev^2))


uppers=ss_q3_b$withinB
lowers=ss_q1_b$withinB
segments(ss_means_a$lev,uppers,ss_means_a$lev,lowers,lwd=1.3)
segments(ss_means_a$lev-.07,uppers,ss_means_a$lev+.07,uppers)
segments(ss_means_a$lev-.07,lowers,ss_means_a$lev+.07,lowers)


points(ss_means_bet$lev,ss_means_bet$between,pch="*")

uppers=ss_q3_bet$between
lowers=ss_q1_bet$between
segments(ss_means_a$lev,uppers,ss_means_a$lev,lowers,lwd=1.3)
segments(ss_means_a$lev-.07,uppers,ss_means_a$lev+.07,uppers)
segments(ss_means_a$lev-.07,lowers,ss_means_a$lev+.07,lowers)


#######
#Plot SVM Figure 5C
######

ss_means_svm=aggregate(svm~lev+subj,data=ss_data,FUN=mean)
ss_q1_svm=aggregate(svm~lev+subj,data=ss_data,FUN=q1)
ss_q3_svm=aggregate(svm~lev+subj,data=ss_data,FUN=q3)

plot(ss_means_svm$lev,ss_means_svm$svm,pch=20,ylim=c(0,1),xlab="tau_p",ylab="Classification Accuracy")
uppers=ss_q3_svm$svm
lowers=ss_q1_svm$svm
segments(ss_means_svm$lev,uppers,ss_means_svm$lev,lowers,lwd=1.3)
segments(ss_means_svm$lev-.07,uppers,ss_means_svm$lev+.07,uppers)
segments(ss_means_svm$lev-.07,lowers,ss_means_svm$lev+.07,lowers)
