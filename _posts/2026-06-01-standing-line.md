---
layout: post
title: Transmission in a queue
date: 2026-01-02
description: Physical distancing and mask wearing have become ubiquitous non-pharmaceutical prevention measures against respiratory infectious diseases. Starting from a situation where susceptible individuals queue while waiting to pass at a counter, we use elementary probability calculus to demonstrate that disease spreads to everyone unless a certain fraction of individuals comply with distancing or mask wearing.   
tags: probability
categories: disease_modeling
related_posts: false
thumbnail: assets/img/post_standing_line-fig1.png
toc:
  sidebar: left
---

Consider a population at risk of an infectious disease, where individuals line up in a queue while waiting to pass
at a counter (e.g. your favourite bakery, security check at the airport). A highly contagious but asymptomatic infectious disease with prevalence $$p$$ is spreading in the population, where infectious hosts only infect their direct neighbours in the queue. To contain the disease, the management of the premise (e.g. bakery or airport) decides to split the individuals into several queues such that no transmission can happen between the queues. Under which condition is this measure sufficient to stop the spread of the disease? 

## The transmission model

We consider a population of $N$ hosts standing in line. Since the disease is asymptomatic, we do not know who is infectious and assume that the probability for a host to be infectious is $$p=\tfrac{k}{N}$$ for $0 < k < N$. The disease is highly contagious such that the probability of transmission $$\beta$$ while standing in the queue for at least 15 minutes is non-negligible (e.g. 0.01 for seasonal influenza, 0.03 for pandemic influenza, 0.1 for on omicron-like covid and 0.9 for measles). 
In our model, we seek to calculate the average number of new infections within a queueing period. To do this, we distinguish hosts with regards to their position in the queue. For hosts at the boundary (i.e. head or tail) of the queue, a new infection arises if the boundary host is not infectious and the immediate neighbour is infectious: 

$$P_{\text{boundary}}=(1-p)p$$

For any susceptible host in the interior of the queue, the probability that at least one of the two possible neighbours is infectious is the complement of the fact the neither of the two is infectious:

$$P_{\text{interior}}=(1-p)(1-(1-p)^2)=p(2-p)$$


<div class="row justify-content-center">
  <div class="col-sm-10 col-md-8">
    {% include figure.liquid path="assets/img/post_standing_line-fig1.png" class="img-fluid rounded" alt="Initial queue with three infectious hosts chosen at random." zoomable=true %}
  </div>
</div>


Now we add up over all positions (2 hosts at the boundary, $N-2$ in the interior), under the assumption that all positions are equally likely to carry infectious hosts. The expected number of new infections within a period is:

$$ \lambda=\mathbb{E}(I)= 2 p(1-p) + (N-2)p(1-p)(2-p)$$

But what does the quantity $\lambda$ tell us about the spread of the disease in the population?
Let's step back from our general formula and assume that only one host is infectious, i.e. $p=\tfrac{1}{N}$. In this case, within in one queue period (e.g. 15 minutes), $\lambda$ new infections would have arised. The queue moves forward in a first-in first-out fashion, the head position leaves the queue, and the tail position is populated by a newly arriving host. 
Since we assumed uniformly random initial conditions, we know that the queue movement does not change the dynamics: on average there the number of infectious hosts rises from $1$ to $1+\lambda$ within a queueing period. Moving to the next queueing period and assuming large $N$, we would end up with $(1+\lambda)^T$ after $T$ periods. From this formula we see, that if $\lambda<1$, then the disease would evenutally die out ("subcritical"), and if $\lambda>1$ the disease would spread to the entire population ("supercritical"). 

This holds only for small $T$ at the beginning of the queueing, at some point, the queue would contain mostly infectious hosts, such that it would be difficult to find a susceptible host for new infections. In epidemiology, this phenomenon is called "depletion of susceptibles".   


## Transmission uncertainty

Here we implicitely assumed that $\beta=1$, i.e. that the probability of transmission was certain, given that a susceptible host is located next to an infectious host. It is straight-forward to change the probability of transmission from an infectious host:

$$ P_{\text{boundary}}=(1-p)\beta p $$

and

$$ P_{\text{interior}}=(1-p)(1-(1-\beta p)^2) $$

such that

$$ \lambda(p,\beta)=\mathbb{E}(I)= (1-p)\left[ 2 \beta p + (N-2)(1-(1-\beta p)^2)\right] $$


<div class="row justify-content-center">
  <div class="col-sm-10 col-md-8">
    {% include figure.liquid path="assets/img/post_standing_line-fig2_lambda.png" class="img-fluid rounded" alt="Expected new infections with transmission uncertainty. The cross shows that in our ideal case with certain transmission at lowest prevalence, the diseaes would spread to the entire population. The colored dots show the prevalence levels above which disease spread becomes supercritical for realistic transmission probabilities." zoomable=true %}
  </div>
</div>  



## Splitting the queue

Let's move to our initial question of splitting the initial queue into two independent queues. It is tempting to split the calculations for $\lambda$ into two parts for $N=2M$:

$$\lambda_2=2 p(1-p) + (M-2)p(1-p)(2-p) + 2 p(1-p) + (M-2)p(1-p)(2-p)$$

From our figure below we can see that the cut creates more boundary points for which by definition the probability of infection is lower than for the interior points:
<div class="row justify-content-center">
  <div class="col-sm-10 col-md-8">
    {% include figure.liquid path="assets/img/post_standing_line-fig2.png" class="img-fluid rounded" alt="Expected new infections with transmission uncertainty. The colored dots show the prevalence levels above which disease spread becomes supercritical." zoomable=true %}
  </div>
</div>                                                                                                                      
Without doing the calculations it follows that cutting the queue in two equal parts would decrease the compound number of new infections. From our above interpretation of $\lambda$ we are actually rather interested in the threshold of switching from supercriticality to subcriticality. Elementary calculations show that

$$\lambda_k=\lambda_1-2kp(1-p)^2$$

such that the condition $\lambda_k$ to become subcritical reads:

$$k>\frac{\lambda_1-1}{2p(1-p)^2}

For our ideal case of $\beta=1$ and $p=\tfrac{1}{N}$ this becomes:

$$k>\frac{N-3}{2} - \frac{2N-1}{2(N-1)^2}

but since $k$ needs to be an integer, the smallest number of cuts such that the disease dynamics becomes subcritical is

$$k_{\text{min}}=\lfloor \frac{N}{2} \rfloor $$

which means that only by cutting the queue into pairs, we can achieve subcriticality.



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
