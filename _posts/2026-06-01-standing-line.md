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
    {% include figure.liquid path="assets/img/post_standing_line-fig2-lambda.png" class="img-fluid rounded" alt="Expected new infections with transmission uncertainty. The cross shows that in our ideal case with certain transmission at lowest prevalence, the diseaes would spread to the entire population. The colored dots show the prevalence levels above which disease spread becomes supercritical for realistic transmission probabilities." zoomable=true %}
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

$$k>\frac{\lambda_1-1}{2p(1-p)^2}$$

For our ideal case of $\beta=1$ and $p=\tfrac{1}{N}$ this becomes:

$$k>\frac{N-3}{2} - \frac{2N-1}{2(N-1)^2}$$

but since $k$ needs to be an integer, the smallest number of cuts such that the disease dynamics becomes subcritical is

$$k_{\text{min}}=\left\lfloor \frac{N}{2} \right\rfloor $$

which means that only by cutting the queue into pairs, we can achieve subcriticality. But in our case we had only $1$ initial infectious individuals. By definition it can only infect its neighbour in the queue, and all other queues would remain uninfected. As we increase the prevalence, even more cuts are needed to remain subcritical, resulting in more and more queues with single individuals. From the point of view of public health, this is effectively imposing physical distancing between individuals. But there is another lever in our model. In lowering the transmission probability $\beta$ by wearing masks, we would need less cuts and could allow queues with larger population size. But this is exactly the trade-off which we experience during the beginning of the most recent pandemic.

<div class="row justify-content-center">
  <div class="col-sm-10 col-md-8">
    {% include figure.liquid path="assets/img/post_standing_line-fig3.png" class="img-fluid rounded" alt="Minimum number of cuts to achieve subcriticality" zoomable=true %}
  </div>                                                                                                                    </div>  
