library(reshape2)
library(ggplot2)
library(deSolve)
library(truncnorm)


influenza_england_1978_school<-
  data.frame(time=c(1:14),
             date=as.Date(c("1978-01-22","1978-01-23","1978-01-24","1978-01-25","1978-01-26","1978-01-27","1978-01-28","1978-01-29","1978-01-30","1978-01-31","1978-02-01","1978-02-02","1978-02-03","1978-02-04")),
             in_bed=c(3,8,26,76,225,298,258,233,189,128,68,29,14,4),
             convalescent=c(0,0,0,0,9,17,105,162,176,166,150,85,47,20),
             total=rep(763,14))


#================#
#= Exercise 1.1 =#
#================#
  ggplot(influenza_england_1978_school)+
  geom_point(aes(x=date,y=in_bed))


SIR.model<-function(time,state,parms){
  with(as.list(c(state,parms)), {
    
    beta=parms[1]
    gamma=parms[2]
    
    S=state[1]
    I=state[2]
    R=state[3]
    N=S+I+R
    
    dS= - beta*I*S/N
    dI= beta*I*S/N - gamma*I
    dR= gamma*I
    
    return(list(c(dS,dI,dR)))
  })
}

time.points<-seq(0,nrow(influenza_england_1978_school),0.01)##time unit in days
initial.condition<-c(S=762,I=1,R=0)##start with one infected person
parameters<-c(beta=1.1,gamma=0.5)

as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR

##lets plot the model output on top of the data
  ggplot(merge(solution_SIR,influenza_england_1978_school))+
  geom_line(aes(x=time,y=I),color="darkgreen")+
  geom_point(aes(x=time,y=in_bed))




#=================#
#= Exercise 1.2  =#
#=================#

##residual sum of squares
rss<-function(parameters){
  time.points<-seq(0,20,0.1)
  initial.condition<-c(S=762,I=1,R=0)
  
  as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR
  
  solution_SIR<-merge(solution_SIR,influenza_england_1978_school)
  RSS<-sum((solution_SIR$I-solution_SIR$in_bed)^2)
  return(RSS)
}

##optimize residual sum of squares over parameter space
start<-c(1.2, 0.4)#initial value
##run newton method for optimization
opt.param <- optim(
  par=start,
  rss,
  method = "L-BFGS-B",
  lower = c(0.1, 0),
  upper = c(2,2),
  hessian = TRUE,
  control = list(parscale = c(10^-4,10^-4),
                 factr = 1))$par
##let's plot the optimal solution
as.data.frame(ode(initial.condition,time.points,SIR.model,opt.param))->solution_SIR_optim

##lets plot the model output on top of the data
  ggplot(merge(solution_SIR_optim,influenza_england_1978_school))+
  geom_line(aes(x=time,y=I),color="darkgreen")+
  geom_point(aes(x=time,y=in_bed))



#=================#
#= Exercise 1.3  =#
#=================#

##the SEIR model
SEIR.model<-function(time,state,parms){
  with(as.list(c(state,parms)), {
    
    alpha=parms[1]
    beta=parms[2]
    gamma=parms[3]
    
    S=state[1]
    E=state[2]
    I=state[3]
    R=state[4]
    N=S+I+R
    
    dS= - beta*I*S/N
    dE= beta*I*S/N - alpha*E
    dI= alpha*E - gamma*I
    dR= gamma*I
    
    return(list(c(dS,dE,dI,dR)))
  })
}


##residual sum of squares
rss1<-function(parameters){
  time.points<-seq(0,20,0.1)
  initial.condition<-c(S=762,E=0,I=1,R=0)
  
  as.data.frame(ode(initial.condition,time.points,SEIR.model,parameters))->solution_SEIR
  
  solution_SEIR<-merge(solution_SEIR,influenza_england_1978_school)
  RSS<-sum((solution_SEIR$I-solution_SEIR$in_bed)^2)
  return(RSS)
}

##optimize residual sum of squares over parameter space
start<-c(0.8,1.2, 0.4)#initial value
##run newton method for optimization
opt.param1 <- optim(
  par=start,
  rss1,
  method = "L-BFGS-B",
  lower = c(0.1,0.1, 0),
  upper = c(10,10,2),
  hessian = TRUE,
  control = list(parscale = c(10^-4,10^-4,10^-4),
                 factr = 1))$par
##let's plot the optimal solution
initial.condition<-c(S=762,E=0,I=1,R=0)
as.data.frame(ode(initial.condition,time.points,SEIR.model,opt.param1))->solution_SEIR_optim

##lets plot the model outputs on top of the data
  ggplot(merge(solution_SEIR_optim,influenza_england_1978_school))+
  geom_line(aes(x=time,y=I),color="darkgreen")+
  geom_point(aes(x=time,y=in_bed))+
  geom_line(data=merge(solution_SIR_optim,influenza_england_1978_school),
            aes(x=time,y=I),color="red")


##let's compare the two models, Occam's razor:
rss(opt.param)*2
rss1(opt.param1)*3




#================#
#= Exercise 2.2 =#
#================#

##get peak prevalence from data
max(influenza_england_1978_school$in_bed)->peak_prev_data

##simulate peak prevalence, fix gamma=1/3
peak_prev_model<-function(parameters=c(beta=1.1,gamma=1/3)){
  
  time.points<-seq(0,nrow(influenza_england_1978_school),0.01)##time unit in days
  initial.condition<-c(S=762,I=1,R=0)
  
  as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR
  
  return(floor(max(solution_SIR$I))/763)
  
}

###prior probability
pbeta(q=peak_prev_model(parameters=c(beta=2,gamma=1/2)),
      shape1=3,
      shape2=3)

###posterior probability of peak prevalence data, using conjugate prior
pbeta(q=peak_prev_model(parameters=c(beta=2,gamma=1/2)),
      shape1=3+peak_prev_data,
      shape2=3+763-peak_prev_data)

##compare prior and posterior distribution
compare_prior_posterior<-
  melt(
  data.frame(prevalence=seq(0,1,0.01),
             prior=dbeta(x=seq(0,1,0.01),shape1=3,shape2=3),
             posterior=dbeta(x=seq(0,1,0.01),shape1=3+peak_prev_data,shape2=3+763-peak_prev_data)),
  id.vars="prevalence")


ggplot(compare_prior_posterior)+
  geom_line(aes(x=prevalence,y=value,linetype=variable))+
  geom_vline(xintercept=peak_prev_data/763,color="red")+
  geom_vline(xintercept=peak_prev_model(parameters=c(beta=1.7,gamma=1/2)),color="blue")+
  annotate("text",x=0.2,y=15,label="model output",color="blue")+
  annotate("text",x=0.5,y=11,label="data",color="red")+
  xlab("peak prevalence")




###prior probability
pbeta(q=peak_prev_model(parameters=c(beta=2,gamma=1/2)),
      shape1=3,
      shape2=3)

###posterior probability, with conjugate prior
pbeta(q=peak_prev_model(parameters=c(beta=2,gamma=1/2)),
      shape1=3+peak_prev_data,
      shape2=3+763-peak_prev_data)

parameters=c(beta=2,gamma=1/2)
initial.condition<-c(S=762,I=1,R=0)##start with one infected person
as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR

ggplot(merge(solution_SIR,influenza_england_1978_school))+
  geom_line(aes(x=time,y=I),color="darkgreen")+
  geom_point(aes(x=time,y=in_bed))


#================#
#= Exercise 2.3 =#
#================#
library(truncnorm)

likelihood<-function(parameters=c(beta=1.1,gamma=1/3),log=T){
  
  ##simulate prevalence, per time point, calculate binomial log likelihood
  as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR
  merge(solution_SIR,influenza_england_1978_school)->solution_SIR
  
  result<-dbinom(solution_SIR$in_bed,
                          size=763,prob=solution_SIR$I/(solution_SIR$S+solution_SIR$I+solution_SIR$R),
                          log=T)
  result<-ifelse(log==T,sum(result),prod(exp(result)))
  
  return(result)
}


MonteCarlo_flat<-function(N=100){
  set.seed(123)
  marginal_likelihood<-NULL
  for (omega in runif(N,min=0,max=2)){
    marginal_likelihood=unlist(c(marginal_likelihood,
                          likelihood(parameters=c(beta=omega,gamma=1/3))))
    
  }
  return(list(integral=sum(marginal_likelihood)/N,
              var_integral=var(marginal_likelihood)/N))
}

MonteCarlo_importance<-function(N=100){
  set.seed(123)
  weight<-function(x){dunif(x,min=0,max=2)/dtruncnorm(x,mean=1.5,sd=1,a=0,b=2)}
  marginal_likelihood<-NULL
  for (omega in rtruncnorm(N,mean=1.5,sd=1,a=0,b=2)){
    marginal_likelihood=unlist(c(marginal_likelihood,
                          likelihood(parameters=c(beta=omega,gamma=1/3))*weight(omega)))
    
  }
  return(list(integral=sum(marginal_likelihood)/N,
              var_integral=var(marginal_likelihood)/N))
  
}

##compare variance of estimators
MonteCarlo_flat(100)
MonteCarlo_importance(100)

##calculate posterior with flat prior and importance sampling
posterior_flat<-function(H){likelihood(parameters=c(beta=H,gamma=1/3))*dunif(H,0,2)/MonteCarlo_flat(100)$integral}
posterior_importance<-function(H){likelihood(parameters=c(beta=H,gamma=1/3))*dunif(H,0,2)/MonteCarlo_importance(100)$integral}



###calculate posterior with conjugate prior
posterior_conjugate<-function(H){
  parameters<-c(beta=H,gamma=1/3)
  as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR
  merge(solution_SIR,influenza_england_1978_school)->output
  
  pbeta(q=output$I/763,
        shape1=3+output$in_bed,
        shape2=3+763-output$in_bed)
}




#==============#
#= Exercise 2.4 =#
#==============#


##define prior probability to start with
prior_probability=function(parameters=c(beta=1.1,gamma=1/3),log=T){
  prior_beta=dunif(parameters["beta"],min=0,max=10, log = T)
  prior_gamma=dunif(parameters["gamma"],min=0,max=1, log = T)
  result<-ifelse(log==F,exp(prior_beta)*exp(prior_beta),prior_beta+prior_beta)
  return(result)
}


##define binomial likelihood function for prevalence, logarithmic by default
likelihood<-function(parameters=c(beta=1.1,gamma=1/3),log=T){
  
  ##simulate prevalence, per time point, calculate binomial log likelihood
  as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR
  merge(solution_SIR,influenza_england_1978_school)->solution_SIR
  
  result<-dbinom(solution_SIR$in_bed,
                 size=763,prob=solution_SIR$I/(solution_SIR$S+solution_SIR$I+solution_SIR$R),
                 log=T)
  result<-ifelse(log==T,sum(result),prod(exp(result)))
  
  return(result)
}

##posterior with Bayesian updating
posterior_probability=function(parameters=c(beta=1.1,gamma=1/3),log=T){
  result<-ifelse(log==T,prior_probability(parameters)+likelihood(parameters),
                 prior_probability(parameters,log=F)*likelihood(parameters,log=F))
  return(result)
}


##draw from proposal density function q
q_sample = function(parameters){
  return(c("beta"=rtruncnorm(1,mean=parameters[["beta"]],sd=0.01,a=0,b=10),
           "gamma"=rtruncnorm(1,mean=parameters[["gamma"]],sd=0.01,a=0,b=1)))
}

##proposal density function q, we assume beta and gamma are independent
q_density = function(parameters_i, parameters_j){
  return(dtruncnorm(parameters_j[["beta"]],mean=parameters_i[["beta"]],sd=0.17,a=0,b=10)*
           dtruncnorm(parameters_j[["gamma"]],mean=parameters_i[["gamma"]],sd=0.044,a=0,b=1))
}


##Metropolis-Hastings algorithm
MCMC_MetropolisHastings=function(N=1000, X0=NULL, my_seed=1){
  
  set.seed(my_seed) #set random number seed
  U=runif(N-1,0,1)  # random number to sample from
  
  # initialisation
  X = X0 = c(beta=1.9,gamma=1/3) #always the same prior, could also use random draws for different seeds of Markov Chain
  pi_i = posterior_probability(X0,log=F)
  markov_chain=data.frame(t(X0),posterior=pi_i,posterior_log=log(pi_i),acc=1)
  
  # iterations
  for (k in 2:N){
    # simulate candidate value Y
    Y = q_sample(X)
    pi_j=posterior_probability(Y,log=F)
    
    # accept candidate Y with Metropolis Hastings probability
    a=min(1,pi_j*q_density(X,Y)/(pi_i*q_density(Y,X)))
    acc=0
    if(U[k-1]<a & is.na(a)==FALSE){
      X<-Y
      pi_i<-pi_j
      acc=1
    }
    markov_chain<-rbind(markov_chain,
                        c(X,pi_i,log(pi_i),acc))
  }
  return(markov_chain)
}

markov_chain_result<-MCMC_MetropolisHastings(N=2000)

hist(markov_chain_result$beta)
hist(markov_chain_result$gamma)


#==============#
#= Exercise 2.5 =#
#==============#

##the SIR model
time.points<-seq(0,nrow(influenza_england_1978_school),1)##time unit in days
initial.condition<-c(S=762,I=1,R=0)##start with one infected person

##setup
N<-1000 #number of particles
epsilon<-100 #tolerance
n_par <- 1 #number of parameters
results <- matrix(ncol=n_par+1,nrow=N) #results
summary_statistics_dist<-function(model,data){sqrt(sum((model-data)^2))}

###ABC rejection algorithm
i<-1 #counter of accepted particles
j<-1 #counter of proposed particles
while(i<=N){
  prior_sample_abc<-runif(1,1.2,2) #sample from prior
  model_output<-data.frame(ode(initial.condition,time.points,SIR.model,parms=c(beta=prior_sample_abc,gamma=0.5)))
  model_output<-merge(model_output,influenza_england_1978_school)
  distance<-summary_statistics_dist(influenza_england_1978_school$in_bed,model_output$I)
  
  if(distance <= epsilon){
    results[i,]<-c(prior_sample_abc,distance)
    i<-i+1
  }
  j<-j+1
  acc_rate<-i/j
  cat("current acceptance rate = ",round(acc_rate,2),"\n")
}

###ABC rejection algorithm with adjustment, for peak prevalence
h<-0.8#bandwidth
N<-1000
##generate prior samples with kernel weights
beta_prior<-runif(N,1.2,2)
results<-NULL
for (i in beta_prior)
{
  model_output<-data.frame(ode(initial.condition,time.points,SIR.model,parms=c(beta=i,gamma=0.5)))
  model_output<-merge(model_output,influenza_england_1978_school)
  kernel_weights<-1/(N*h)*sum(dnorm((max(influenza_england_1978_school$in_bed)-max(model_output$I))/h))##kernel weights with bandwidth h
  
  results<-rbind(results,
                 c(beta=i,
                   summary_statistics=max(model_output$I),
                   weight=kernel_weights))
}
results<-data.frame(results)
##now do adjustment
lm(beta~summary_statistics,results,weights=weight)->regression
beta_posterior<-predict(regression)
hist(beta_posterior)





#=================#
#= Exercise 3.1  =#
#=================#
###Gillespie's direct algorithm
Gillespie_direct<-function(N=763,seed=123,parameters=c(beta=2,gamma=0.5)){
  
  set.seed(seed)
  ##intial conditions
  states<-c(time=0,S=N-1,I=1,R=0)
  traj<-rbind(c(time=0,S=N,I=0,R=0),c(time=0,S=N-1,I=1,R=0))
  ##events:
  rates=c(rep(parameters['beta']*states['I']/N,states['S']),
          rep(parameters['gamma'],states['I']))
  
  while(tail(traj[,'I'],1)>0){
    states=traj[nrow(traj),]
    rates=c(rep(parameters['beta']*states['I']/N,states['S']),
            rep(parameters['gamma'],states['I']))
    time_to_event=rexp(n=1,rate=sum(rates))
    which_event=sample(x=names(rates),size=1,prob=rates/sum(rates))
    #stoichometric update
    if(which_event=="beta"){
      print("new infection")
      traj<-rbind(traj,
                  states+c(time_to_event,-1,+1,0))
    }
    if(which_event=="gamma"){
      print("clearance")
      traj<-rbind(traj,
                  states+c(time_to_event,0,-1,+1))
    }
  }
  traj<-as.data.frame(tail(traj,-1))
  traj$seed<-seed
  return(traj)
}

##use parameters obtained from RSS optimisation: c(beta=1.6692258,gamma=0.4434502)
runs<-NULL
for (seed in 1:700){
  runs<-rbind(runs,
              Gillespie_direct(seed=seed,parameters=c(beta=1.6692258,gamma=0.4434502)))
}

###SIR ODE 
time.points<-seq(0,20,0.01)##time unit in days
initial.condition<-c(S=762,I=1,R=0)##start with one infected person
parameters<-c(beta=1.6692258,gamma=0.4434502)
as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))->solution_SIR

##lets plot the models output on top of the data
ggplot(runs)+
  geom_line(aes(x=time,y=I,group=as.factor(seed)),color="grey")+
  geom_line(data=solution_SIR,aes(x=time,y=I),color="darkgreen",size=3)+
  geom_point(data=influenza_england_1978_school,aes(x=time,y=in_bed),color="red")+
  xlim(c(0,20))+
  theme_bw()->gillespie_direct_plot



#=================#
#= Exercise 4.1  =#
#=================#

##sobol index for SIR model for MCMC posterior for f=peak prevalence

posterior_beta<-tail(markov_chain_result$beta,1000)
posterior_gamma<-tail(markov_chain_result$gamma,1000)

f<-function(parameters=c(beta=1.1,gamma=0.5)){
  return(max(as.data.frame(ode(initial.condition,time.points,SIR.model,parameters))$I))
}
##V_beta, V_gamma
##sample from posterior 1000 times, 1000 rows, 2 columns 
N=1000

A<-cbind(sample(posterior_beta,N),
         sample(posterior_gamma,N))
B<-cbind(sample(posterior_beta,N),
         sample(posterior_gamma,N))

A_beta=A_gamma<-A
A_beta[,1]<-B[,1]
A_gamma[,2]<-B[,2]


V_beta=V_gamma<-NULL
for(j in 1:N){
  V_beta<-c(V_beta,
  f(B[j,])*(f(A_beta[j,])-f(A[j,])))

  V_gamma<-c(V_gamma,
  f(B[j,])*(f(A_gamma[j,])-f(A[j,])))
}
V_beta<-1/N*sum(V_beta)
V_gamma<-1/N*sum(V_gamma)

##V_total
samples_A<-NULL
for (j in 1:N){
  samples_A<-c(samples_A,
             f(B[j,]))
}
V_total=var(samples_A)

##Sobol indices
S_beta=V_beta/V_total
S_gamma=V_gamma/V_total

##interaction
S_interaction=(V_total-V_beta-V_gamma)/V_total


