
#SIM 1: Just taup varies

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
par.tau0 = seq(0,5,0.5)
nlevs=length(par.tau0)
par.taup = 0


par.e0vs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,rep(par.tau0,each=nsubject*nsims)),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)
par.epvs = matrix(rnorm(nvox*nsubject*nsims*nlevs,0,par.taup),nrow=nvox,ncol=nsubject*nsims*nlevs,byrow = TRUE)

#subject-level
par.sigma0 = 0 
par.sigmap = 0
par.e0s = rnorm(nsubject*nsims*nlevs,0,par.sigma0) 
par.eps = rnorm(nsubject*nsims*nlevs,0,par.sigmap)


ss_data = matrix(nrow=nsubject*nsims*nlevs,ncol=7)
ss_data = cbind(ss_data,rep(par.tau0,each=nsubject*nsims),rep(1:nsubject,each=nlevs*nsims))

	
for(subj in 1:(nsubject*nsims*nlevs)){	
	act_patterns=act_gen(xp, par.gamma0,par.gammap,par.e0s[subj],par.eps[subj],par.e0vs[,subj],par.epvs[,subj],par.etvs)	
	sim_mat=cor(act_patterns$full)
	mean_act=apply(act_patterns$full, 2, mean)
	
	temp_withinbetween = within_between_sim(sim_mat,xp)
	
	temp_act = tapply(mean_act,xp,mean)
	
	ss_data[subj,1:7]=c(temp_withinbetween,temp_act)
	
}

ss_data=data.frame(ss_data)


names(ss_data) = c('all','within','between','withinA','withinB','actA','actB','lev','subj')

ss_means=aggregate(all~lev+subj,data=ss_data,FUN=mean)
ss_q1=aggregate(all~lev+subj,data=ss_data,FUN=q1)
ss_q3=aggregate(all~lev+subj,data=ss_data,FUN=q3)

###########################
#Plot Mean correlation (Figure 5A)
###########################

plot(ss_means$lev,ss_means$all,pch=20,ylim=c(0,1),xlab="tau0",ylab="Pearson's r")

uppers=ss_q3$all
lowers=ss_q1$all
segments(ss_means$lev,uppers,ss_means$lev,lowers,lwd=1.3)
segments(ss_means$lev-.07,uppers,ss_means$lev+.07,uppers)
segments(ss_means$lev-.07,lowers,ss_means$lev+.07,lowers)
lines(ss_means$lev,ss_means$lev^2/(par.etvs+ss_means$lev^2))




