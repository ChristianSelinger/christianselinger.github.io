---
layout: post
title: Transmission in a queue
date: 2026-01-02
description: Physical distancing and mask wearing have become ubiquitous non-pharmaceutical prevention measures against respiratory infectious diseases. Starting from a situation where susceptible individuals queue while waiting to pass at a counter, we use elementary probability calculus to demonstrate that disease spreads to everyone unless a certain fraction of individuals comply with distancing or mask wearing.   
tags: probability
categories: disease_modeling
related_posts: false
thumbnail: assets/img/dec_sisgmae-x_I.png
toc:
  sidebar: left
---

Consider a population at risk of an infectious disease, where individuals line up in a queue while waiting to pass
at a counter (e.g. your favourite bakery, security check at the airport). A highly contagious but asymptomatic infectious disease with prevalence $$p$$ is spreading in the population, where infectious hosts only infect their direct neighbours in the queue. To contain the disease, the management of the premise (e.g. bakery or airport) decides to split the individuals into several queues such that no transmission can happen between the queues. Is this measure sufficient to stop the spread of the disease? 

## The transmission model

We consider a population of $N$ hosts standing in line. Since the disease is asymptomatic, we do not know who is infectious and assume that the probability for a host to be infectious is $$p=\tfrac{k}{N}$$ for $0 < k < N$. The disease is highly contagious such that the probability of transmission $$\beta$$ while standing in the queue for at least 15 minutes is non-negligible (e.g. 0.01 for seasonal influenza, 0.03 for pandemic influenza, 0.1 for on omicron-like covid and 0.9 for measles). 
In our model, we seek to calculate the average number of new infections within a queueing period. To do this, we distinguishhosts with regards to their position in the queue. For hosts at the boundary (i.e. head or tail) of the queue, a new infection arises if the boundary host is not infectious and the immediate neighbour is infectious: 
$P_{\text{boundary}}=(1-p)$.
For any host in the interior of the queue, the probability that at least one of the two possible neighbours is infectious is the complement of the fact the neither of the two is infectious:
$P_{\text{interior}}=1-(1-p)^2=p(2-p)$.







$$\begin{align}
    \frac{dS_U}{dt} &= -(1-\sigma) \beta S_U (I_U + I_N) + \gamma I_U \\
    \frac{dI_U}{dt} &= (1-\sigma) \beta S_U (I_U + I_N) - \gamma I_U \\
    \frac{dS_N}{dt} &= -\beta S_N (I_U + I_N) + \gamma I_N \\
    \frac{dI_N}{dt} &= \beta S_N (I_U + I_N) - \gamma I_N
\end{align}$$

Let $$x$$ denote the frequency of intervention users in the population ($$x = S_U + I_U$$). Thus setting $$S_U=x-I_U$$ and $$S_V=(1-x)-I_V$$, the system can be reduced to
$$\begin{align}
    \frac{dI_U}{dt} &= (1-\sigma) \beta (x-I_U) (I_U + I_N) - \gamma I_U \\
    \frac{dI_N}{dt} &= \beta (1-x-I_N) (I_U + I_N) - \gamma I_N
\end{align}$$

Using the next generation matrix approach to calculate the basic reproduction number $R_0$:

$$F = \begin{pmatrix} (1-\sigma)\beta x & (1-\sigma)\beta x \\ \beta(1-x) & \beta(1-x) \end{pmatrix}, \quad V = \begin{pmatrix} \gamma & 0 \\ 0 & \gamma \end{pmatrix}$$ 

The Next-Generation Matrix is 
$K = FV^{-1}$:

$$K = \begin{pmatrix} \frac{(1-\sigma)\beta x}{\gamma} & \frac{(1-\sigma)\beta x}{\gamma} \\ \frac{\beta(1-x)}{\gamma} & \frac{\beta(1-x)}{\gamma} \end{pmatrix}$$ 

$R_0$ is the spectral radius (the largest eigenvalue) of matrix $K$ and reads $$R_0(x,\sigma) = \frac{(1-\sigma)\beta x}{\gamma} + \frac{\beta(1-x)}{\gamma}$$. This gives $$R_0(x,\sigma) = \frac{\beta}{\gamma} ( (1-\sigma)x + (1-x) )=\frac{\beta}{\gamma} ( 1-\sigma x  )$$

The disease dies out if $$R_0(x,\sigma)<1 \Leftrightarrow  x > \frac{1}{\sigma} (1 - \frac{\gamma}{\beta})$$.

We call $$x^e=\frac{1}{\sigma} (1 - \frac{\gamma}{\beta})$$ the user frequency elimination threshold.


## Replicator dynamics for intervention usage

The decision to use an intervention is modeled as a game. Like in a lottery, each host can play two different strategies against the bank, which is here the pool of infectious hosts in the population. If the host decides to not use the intervention, she will suffer from the potential consequences of an infection at cost $$C_I>0$$. If the host decides to play as a user, he will suffer a little less from the consequences of an infection, but incur the nuisance cost $$C_N\geq 0$$ of using the health intervention. For both strategies, we can write the individual payoff as negative values of the cost:

 $$f_U =  -C_N - (1-\sigma) \beta I C_I$$

 $$f_N =  -\beta I C_I$$

The population payoff is

$$\bar{f}(x)= xf_U+ (1-x)f_N$$

We assume that individuals do not interact with other individuals’ strategies to determine their own cost, like it would be in a card game. To play against the entire infectious population is actually in accordance with the mass-action principle for disease transmission.

The evolution of user frequency $$x$$ satisfies the replicator equation:

$$\begin{equation} \dot{x} = x (f_U - \bar{f}(x)) \end{equation}$$

In our case, this equation simplifies to

$$\begin{equation} \dot{x} = x ((1-x) f_U - (1-x) f_N) = x(1-x)(f_U-f_N) \end{equation}$$

also known as **logistic growth** equation. 

The intervention user behavior dynamics is governed by 

$$\begin{equation}
    \dot{x} = x(1-x) ( \sigma\beta I C_I  - C_N )
\end{equation}$$

We distinguish:
- **Adoption:** $$x$$ increases if $$\beta I > \frac{C_N}{C_I \sigma}$$. The risk-reduction benefit outweighs the nuisance cost.
- **Abandonment:** $$x$$ decreases if $$\beta I < \frac{C_N}{C_I \sigma}$$. The nuisance cost is perceived as greater than the expected cost of infection.
- **Evolutionary Stability:** A mixed equilibrium exists when $$\beta I = \frac{C_N}{C_I \sigma}$$.



## Equilibrium for infectious disease dynamics

To find the endemic equilibrium of the reduced system 
$$\begin{align}
    \frac{dI_U}{dt} &= (1-\sigma) \beta (x-I_U) (I_U + I_N) - \gamma I_U \\
    \frac{dI_N}{dt} &= \beta (1-x-I_N) (I_U + I_N) - \gamma I_N
\end{align}$$
We need to consider the case $$R_0(x,\sigma)>1$$ and solve $$dI_U/dt = 0$$ and $$dI_N/dt = 0$$.

Let $$I = I_U + I_N$$ be the total prevalence. 

Then $$I_U^* = \frac{(1-\sigma) \beta x I}{(1-\sigma) \beta I + \gamma}$$ and $$I_N^* = \frac{\beta (1-x) I}{\beta I + \gamma}.$$ 

Since $I = I_U + I_N$, we sum the two equations: $$I =( \frac{(1-\sigma) \beta x}{(1-\sigma) \beta I + \gamma} + \frac{\beta (1-x)}{\beta I + \gamma} ) I$$

If $I > 0$, we can divide both sides by $I$:

$$1 = \frac{(1-\sigma) \beta x}{(1-\sigma) \beta I + \gamma} + \frac{\beta (1-x)}{\beta I + \gamma}$$

This is a quadratic equation in $I$: 
$$\begin{equation}
\beta^2 (1-\sigma) I^2 + \beta ( \gamma(2-\sigma) - \beta(1 - \sigma ) ) I + \gamma^2 ( 1 - R_0(x,\sigma) ) = 0
\end{equation}$$

The positive root of this quadratic gives the unique endemic equilibrium level when $$R_0(x,\sigma ) > 1$$:

$$I^* = \frac{-B + \sqrt{B^2 - 4AC}}{2A}$$ for $$A = \beta^2 (1-\sigma)$$, $$B = \beta \gamma (2-\sigma) - \beta^2(1 - \sigma )$$ and $$C = \gamma^2 (1 - R_0(x,\sigma))$$.


## Equilibrium for the coupled disease–replicator dynamics

We can also solve the quadratic equation in $$I$$ for $$x$$ to obtain the user frequency at the endemic equilibrium. From the replicator equation we know that the user frequency is at equilibrium if the payoffs for users and non-users are equal: $$f_U = f_N$$ which amounts to mainting a specific prevalence:

   $$ \bar{I} = \frac{C_N}{\sigma\beta C_I }$$ in the population.
   
Solving the quadratic equation for $$x$$ by using the expression for $$R_0(x,\sigma)$$ yields

$$\bar{x} = \frac{1}{\sigma}(1-\frac{\gamma}{\beta})- \frac{(\gamma(2-\sigma) - \beta(1 - \sigma ))}{\gamma \sigma} \bar{I}- \frac{\beta (1-\sigma) }{\gamma \sigma}\bar{I}^2$$

The strategy $$\bar{x}$$ is an **evolutionary stable strategy** for the case where the disease will persist and eventually reach an endemic equilbrium $$\bar{I}$$. The coupled disease-replicator dynamics is at equilibrium $$(\bar{I},\bar{x})$$, which depends on both disease or efficacy parameters $$\beta, \gamma, \sigma$$ and intervention usage behavior parameters $$C_N$$ and $$C_I$$.


## Elimination threshold and eradication gap

Let us recall the elimination threshold of user frequency $$x^e = \frac{1}{\sigma} \left( 1 - \frac{\gamma}{\beta} \right)$$

It follows that $$\begin{equation}\bar{x} = x^e + \frac{\beta(1-\sigma) - \gamma(2-\sigma)}{\gamma \sigma} \bar{I} - \frac{\beta (1-\sigma)}{\gamma \sigma} \bar{I}^2\end{equation}$$

The difference $$\Delta=x^e-\bar{x}$$ is called the **eradication gap**. 

$$\Delta = \frac{\bar{I}}{\sigma} (2 - \sigma - \frac{\beta}{\gamma}(1-\sigma) + \frac{\beta}{\gamma}(1-\sigma)\bar{I})>\frac{\bar{I}}{\sigma} (2 - \sigma - \frac{\beta}{\gamma}(1-\sigma))= \frac{C_N}{\beta C_I }(1 + (1-\sigma)(1-\frac{\beta}{\gamma}))$$

The lower bound for the eradication gap 

$$\Delta > \frac{C_N}{\beta C_I }(1 + (1-\sigma)(1-\frac{\beta}{\gamma}))$$

is positive if the basic reproduction number and intervention efficacy satisfy:

$$\frac{\beta}{\gamma}<\frac{2-\sigma}{1-\sigma}$$

As the nuisance cost $$C_N$$ increases, the gap is widening. Equality between the elimination threshold and the evolutionary stable strategy can only hold if the nuisance cost $$C_N=0$$. This means that as long as interventions are effective and have a non-zero nuisance cost, non-users will free ride and jeopardize eradication.

<div class="row justify-content-center">
  <div class="col-sm-10 col-md-8">
    {% include figure.liquid path="assets/img/dec_sisgmae-x_I.png" class="img-fluid rounded" alt="Phase diagram of user frequency versus disease prevalence" zoomable=true %}
  </div>
</div>

The figure depicts the phase diagram of the coupled model of user frequency and disease prevalence for two different choices of efficacy and nuisance. The dynamical systems starts at the red point. At <span class="fig-orange">low efficacy and low nuisance</span>, everyone will end up using the intervention and the equilibrium is at high prevalence. Increasing nuisance at low efficacy will result in no users and the entire population infected. <span class="fig-blue">Increasing efficacy but having low nuisance</span> creates an opportunity for free-riders, resulting in lower endemic equilibrium and intermediate user frequency. Again, increasing nuisance in the high efficacy case will deter users and result in a high endemic equilibrium. The elimination threshold for $$\sigma=0.9$$ is located at the cube on the x-axis, at zero prevalence. The disease parameters are fixed at $\beta=0.5$ and $\gamma=0.1$.
