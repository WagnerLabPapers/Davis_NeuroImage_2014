
#SIM 3: Just taup varies; Continuous X variable
library(e1071)

source('~/Dropbox/DLMNP_sims/sim_funcsf.r')

# "experiment" parameters
nsims=100
nstim=60
nvox=50
nsubject=1 # number subjects per simulation
#conditioning variable X
xp=rep(c(1:6),each=10)

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


ss_data = matrix(nrow=nsubject*nsims*nlevs,ncol=10)
ss_data = cbind(ss_data,rep(par.taup,each=nsubject*nsims),rep(1:nsubject,each=nlevs*nsims))

	
for(subj in 1:(nsubject*nsims*nlevs)){	
	act_patterns=act_gen(xp, par.gamma0,par.gammap,par.e0s[subj],par.eps[subj],par.e0vs[,subj],par.epvs[,subj],par.etvs)	
	sim_mat=cor(act_patterns$full)
	mean_act=apply(act_patterns$full, 2, mean)
	
	temp_withinbetween = within_between_sim(sim_mat,xp)
	
	
	#SVM
	xp_short=xp[c(1:5,11:15,21:25,31:35,41:45,51:55)]
	train=data.frame(t(act_patterns$full[,c(1:5,11:15,21:25,31:35,41:45,51:55)]))
	test=data.frame(t(act_patterns$full[,c(6:10,16:20,26:30,36:40,46:50,56:60)]))
	
	svm1=svm(as.factor(xp_short)~.,data=train,type='nu-classification',nu=0.01,kernel='linear')
perf_svm1=sum(predict(svm1,test)==xp_short)/30

	
	ss_data[subj,1:10]=c(temp_withinbetween,perf_svm1)
	
	
}

ss_data=data.frame(ss_data)


names(ss_data) = c('all','within','between','within1','within2','within3','within4','within5','within6','svm','lev','subj')
#######
#Plot mean correlation (Figure 5D)
#######


ss_means_1=aggregate((within1)~lev+subj,data=ss_data,FUN=mean)
ss_q1_1=aggregate((within1)~lev+subj,data=ss_data,FUN=q1)
ss_q3_1=aggregate((within1)~lev+subj,data=ss_data,FUN=q3)

ss_means_2=aggregate((within2)~lev+subj,data=ss_data,FUN=mean)
ss_q1_2=aggregate((within2)~lev+subj,data=ss_data,FUN=q1)
ss_q3_2=aggregate((within2)~lev+subj,data=ss_data,FUN=q3)

ss_means_3=aggregate((within3)~lev+subj,data=ss_data,FUN=mean)
ss_q1_3=aggregate((within3)~lev+subj,data=ss_data,FUN=q1)
ss_q3_3=aggregate((within3)~lev+subj,data=ss_data,FUN=q3)

ss_means_4=aggregate((within4)~lev+subj,data=ss_data,FUN=mean)
ss_q1_4=aggregate((within4)~lev+subj,data=ss_data,FUN=q1)
ss_q3_4=aggregate((within4)~lev+subj,data=ss_data,FUN=q3)

ss_means_5=aggregate((within5)~lev+subj,data=ss_data,FUN=mean)
ss_q1_5=aggregate((within5)~lev+subj,data=ss_data,FUN=q1)
ss_q3_5=aggregate((within5)~lev+subj,data=ss_data,FUN=q3)

ss_means_6=aggregate((within6)~lev+subj,data=ss_data,FUN=mean)
ss_q1_6=aggregate((within6)~lev+subj,data=ss_data,FUN=q1)
ss_q3_6=aggregate((within6)~lev+subj,data=ss_data,FUN=q3)


ss_means_bet=aggregate(between~lev+subj,data=ss_data,FUN=mean)
ss_q1_bet=aggregate(between~lev+subj,data=ss_data,FUN=q1)
ss_q3_bet=aggregate(between~lev+subj,data=ss_data,FUN=q3)



plot(ss_means_1$lev,ss_means_1[,3],pch=20,ylim=c(0,1),xlab="tau_p",ylab="Pearson's r")

uppers=ss_q3_1[,3]
lowers=ss_q1_1[,3]
segments(ss_means_1$lev,uppers,ss_means_1$lev,lowers,lwd=1.3)
segments(ss_means_1$lev-.07,uppers,ss_means_1$lev+.07,uppers)
segments(ss_means_1$lev-.07,lowers,ss_means_1$lev+.07,lowers)

points(ss_means_2$lev,ss_means_2[,3],pch=18)


uppers=ss_q3_2[,3]
lowers=ss_q1_2[,3]
segments(ss_means_2$lev,uppers,ss_means_2$lev,lowers,lwd=1.3)
segments(ss_means_2$lev-.07,uppers,ss_means_2$lev+.07,uppers)
segments(ss_means_2$lev-.07,lowers,ss_means_2$lev+.07,lowers)

points(ss_means_3$lev,ss_means_3[,3],pch=24)


uppers=ss_q3_3[,3]
lowers=ss_q1_3[,3]
segments(ss_means_3$lev,uppers,ss_means_3$lev,lowers,lwd=1.3)
segments(ss_means_3$lev-.07,uppers,ss_means_3$lev+.07,uppers)
segments(ss_means_3$lev-.07,lowers,ss_means_3$lev+.07,lowers)

points(ss_means_4$lev,ss_means_4[,3],pch=20)


uppers=ss_q3_4[,3]
lowers=ss_q1_4[,3]
segments(ss_means_4$lev,uppers,ss_means_4$lev,lowers,lwd=1.3)
segments(ss_means_4$lev-.07,uppers,ss_means_4$lev+.07,uppers)
segments(ss_means_4$lev-.07,lowers,ss_means_4$lev+.07,lowers)

points(ss_means_5$lev,ss_means_5[,3],pch=21)


uppers=ss_q3_5[,3]
lowers=ss_q1_5[,3]
segments(ss_means_5$lev,uppers,ss_means_5$lev,lowers,lwd=1.3)
segments(ss_means_5$lev-.07,uppers,ss_means_5$lev+.07,uppers)
segments(ss_means_5$lev-.07,lowers,ss_means_5$lev+.07,lowers)

points(ss_means_6$lev,ss_means_6[,3],pch=22)


uppers=ss_q3_6[,3]
lowers=ss_q1_6[,3]
segments(ss_means_6$lev,uppers,ss_means_6$lev,lowers,lwd=1.3)
segments(ss_means_6$lev-.07,uppers,ss_means_6$lev+.07,uppers)
segments(ss_means_6$lev-.07,lowers,ss_means_6$lev+.07,lowers)


points(ss_means_bet$lev,ss_means_bet$between,pch="*")

uppers=ss_q3_bet$between
lowers=ss_q1_bet$between
segments(ss_means_1$lev,uppers,ss_means_1$lev,lowers,lwd=1.3)
segments(ss_means_1$lev-.07,uppers,ss_means_1$lev+.07,uppers)
segments(ss_means_1$lev-.07,lowers,ss_means_1$lev+.07,lowers)

################
#PLOT SVM Figure 5E
################

ss_means_svm=aggregate(svm~lev+subj,data=ss_data,FUN=mean)
ss_q1_svm=aggregate(svm~lev+subj,data=ss_data,FUN=q1)
ss_q3_svm=aggregate(svm~lev+subj,data=ss_data,FUN=q3)

plot(ss_means_svm$lev,ss_means_svm$svm,pch=20,ylim=c(0,1),xlab="tau_p",ylab="Classification Accuracy")
uppers=ss_q3_svm$svm
lowers=ss_q1_svm$svm
segments(ss_means_svm$lev,uppers,ss_means_svm$lev,lowers,lwd=1.3)
segments(ss_means_svm$lev-.07,uppers,ss_means_svm$lev+.07,uppers)
segments(ss_means_svm$lev-.07,lowers,ss_means_svm$lev+.07,lowers)
