library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)


wdir<-"C:/Users/selich/work/education/ContagiousMathmatics"

###endemic equilibrium of coupled system
x_star_fct<-function(CN,CI,k,beta,sigma,gamma){
  # Define I_crit
  I_crit <- CN / (k*sigma * beta * CI)
  # Calculate x_star
  numerator <- (gamma + (1 - sigma) * beta * I_crit) * (beta - gamma - beta * I_crit)
  denominator <- sigma * beta * gamma
  x_star <- numerator / denominator
  return(x_star)
}

I_star_fct<-function(CN,CI,beta,sigma,gamma,x){
  A=beta^2*(1-sigma)
  B=beta*gamma*(2-sigma)-beta^2*(1-sigma)
  C=gamma^2*(1-beta/gamma*(1-sigma*x))
  I_star=(-B+sqrt(B^2-4*A*C))/(2*A)
  return(I_star)
}


# Parameters, CN=0.1
p_test <- c(beta=0.5, gamma=0.1, CN=0.1, CI=2.0, k=2)
sigmas <- c(0.5,0.9)
init <- c(x = 0.5, IU = 0.05, IN = 0.1)
times <- seq(0, 1000, by = 0.1)

bednet_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    I <- IU + IN
    dx <- k * x * (1 - x) * (beta * I * CI * sigma - CN)
    dIU <- (1 - sigma) * beta * (x - IU) * I - gamma * IU
    dIN <- beta * ((1 - x) - IN) * I - gamma * IN
    return(list(c(dx, dIU, dIN)))
  })
}

# Simulate for each sigma
all_data <- lapply(sigmas, function(s) {
  params <- c(p_test, sigma = s)
  out <- as.data.frame(ode(y = init, times = times, func = bednet_model, parms = params))
  out$sigma <- s
  return(out)
}) %>% bind_rows%>%cbind(as.list(p_test))

# Parameters, CN=0.5
p_test <- c(beta=0.5, gamma=0.1, CN=0.5, CI=2.0, k=2)

# Simulate for each sigma
all_data <- all_data%>%bind_rows(
  lapply(sigmas, function(s) {
  params <- c(p_test, sigma = s)
  out <- as.data.frame(ode(y = init, times = times, func = bednet_model, parms = params))
  out$sigma <- s
  return(out)
}) %>% bind_rows%>%cbind(as.list(p_test))
)

# Parameters, CN=0.01
p_test <- c(beta=0.5, gamma=0.1, CN=0.01, CI=2.0, k=2)

# Simulate for each sigma
all_data <- all_data%>%bind_rows(
  lapply(sigmas, function(s) {
    params <- c(p_test, sigma = s)
    out <- as.data.frame(ode(y = init, times = times, func = bednet_model, parms = params))
    out$sigma <- s
    return(out)
  }) %>% bind_rows%>%cbind(as.list(p_test))
)

all_data%>%
  mutate(I=IU+IN,
         gr=paste0(sigma,CN),
         Ibar=CN/(sigma*beta*CI),
         x_star=1/sigma*((1-gamma/beta)+beta/gamma*(1-sigma)*Ibar*(1-Ibar)-(2-sigma)*Ibar),
         x_e=1/sigma*(1-gamma/beta))->all_data

all_data%>%
  mutate(x_e=if_else(x_e<1,x_e,NA))%>%
  filter(CN!=0.01)%>%
ggplot()+
  geom_path(aes(x = x, y = I,color=as.factor(sigma),group=as.factor(gr),linetype=as.factor(CN)),linewidth = 1)+
  scale_color_manual(name = "efficacy (sigma)",values=c("orange","navyblue","chartreuse4")) +
  geom_vline(aes(xintercept = x_e,color=as.factor(sigma)))+
  annotate("point", y = init["IU"]+init["IN"], x = init["x"], color = "red", size = 3) +
  scale_x_continuous(name="user frequency (x)")+
  scale_y_continuous(name="disase prevalence (I)")+
  scale_linetype_manual(name="nuisance cost (CN)",values = c("solid", "dotted"),labels=c("low","high"))+
  theme_minimal()->plot_x_I

ggsave(file.path(wdir,"dec_sisgmae-x_I.pdf"), plot = plot_x_I, width = 6, height = 4)


all_data%>%
  filter(sigma==0.9)%>%
  mutate(eradication_gap=x_e-x)%>%
  ggplot()+
  geom_line((aes(x=time,y=eradication_gap,linetype=as.factor(CN))),linewidth = 1,color="navyblue")+
  scale_y_continuous(name="eradication gap (x^e-x)")+
  scale_linetype_manual(name="nuisance cost (CN)",
                        values = c("dashed","solid", "dotted"),
                        labels=c("very low","low","high"))+
  theme_minimal()->plot_gap
  

ggsave(file.path(wdir,"dec_sisgmae-gap.pdf"), plot = plot_gap, width = 6, height = 4)


library(cowplot)

plot_grid(plot_x_I+guides(linetype = "none"),
          plot_gap, 
          labels = "AUTO", ncol = 2, nrow = 1)->plt_grid
ggsave(file.path(wdir,"dec_sisgmae-combo.pdf"), plot = plt_grid, width = 12, height = 4)


