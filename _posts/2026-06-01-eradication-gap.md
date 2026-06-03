---
layout: post
title: Benefits and nuisance of health interventions
date: 2026-01-01
description: Coupling disease transmission dynamics and game theory, we investigate how individual preference can lead to endemic disease accepting some level nuisance individual sacrifice will open
tags: game_theory differential_equations
categories: sample-posts
related_posts: false
---


Consider a population at risk of an infectious disease, where individuals recover
but will eventually become susceptible again. The population is divided into
into two groups. The first group decides to use some health intervention, which comes to the
benefit of preventing infections but also to some nuisance cost. The non-users
will not use the health intervention, they are at higher risk of infection, but do not bear the
nuisance cost.

\paragraph{The SIS Epidemiological Model}
We consider a population divided into two strategy groups: Users ($U$) and Non-users ($N$). Each group is further divided into Susceptible ($S$) and Infectious ($I$) compartments. The total population is normalized to $1$, where:
\[ S_U + I_U + S_N + I_N = 1 \]

The transmission is governed by the force of infection $\lambda = (1-\sigma)\beta (I_U + I_N)$, where $\beta$ is the transmission rate and $\sigma$ is the protective efficacy of the bednet ($0 \leq \sigma \leq 1$). The recovery rate is $\gamma$. The system of differential equations is:

\begin{align}
    \frac{dS_U}{dt} &= -(1-\sigma) \beta S_U (I_U + I_N) + \gamma I_U \\
    \frac{dI_U}{dt} &= (1-\sigma) \beta S_U (I_U + I_N) - \gamma I_U \\
    \frac{dS_N}{dt} &= -\beta S_N (I_U + I_N) + \gamma I_N \\
    \frac{dI_N}{dt} &= \beta S_N (I_U + I_N) - \gamma I_N
\end{align}
Let $x$ denote the frequency of bednet users in the population ($x = S_U + I_U$). Thus setting $S_U=x-I_U$ and $S_V=(1-x)-I_V$, the system can be reduced to
\begin{align*}
    \frac{dI_U}{dt} &= (1-\sigma) \beta (x-I_U) (I_U + I_N) - \gamma I_U \\
    \frac{dI_N}{dt} &= \beta (1-x-I_N) (I_U + I_N) - \gamma I_N
\end{align*}
Using the next generation matrix approach to calculate the basic reproduction number $R_0$:$$F = \begin{pmatrix} (1-\sigma)\beta x & (1-\sigma)\beta x \\ \beta(1-x) & \beta(1-x) \end{pmatrix}, \quad V = \begin{pmatrix} \gamma & 0 \\ 0 & \gamma \end{pmatrix}$$ The Next-Generation Matrix is $K = FV^{-1}$:$$K = \begin{pmatrix} \frac{(1-\sigma)\beta x}{\gamma} & \frac{(1-\sigma)\beta x}{\gamma} \\ \frac{\beta(1-x)}{\gamma} & \frac{\beta(1-x)}{\gamma} \end{pmatrix}$$ $R_0$ is the spectral radius (the largest eigenvalue) of matrix $K$. $$R_0(x,\sigma) = \frac{(1-\sigma)\beta x}{\gamma} + \frac{\beta(1-x)}{\gamma}$$ This gives $$R_0(x,\sigma) = \frac{\beta}{\gamma} \left[ (1-\sigma)x + (1-x) \right]=\frac{\beta}{\gamma} ( 1-\sigma x  )$$
The disease dies out if $$R_0(x,\sigma)<1 \Leftrightarrow  x > \frac{1}{\sigma} (1 - \frac{\gamma}{\beta})$$ We call $x^e=\frac{1}{\sigma} (1 - \frac{\gamma}{\beta})$ the user frequency elimination threshold.

\paragraph{Replicator dynamics of net usage}
The decision to use a bednet is modeled as a game. Like in a lottery, each host can play two different strategies against the bank, which is here the pool of infectious hosts in the population. If the host decides to not use the net, she will suffer from the potential consequences of an infection at cost $C_I>0$. If the host decides to play as a user, he will suffer a little less from the consequences of an infection, but incur the nuisance cost $C_N>0$ of putting and maintaining the net every day. For both strategies we can write the indivdual payoff as negative values of the cost:  
\begin{itemize}
    \item $f_U =  -C_N - (1-\sigma) \beta I C_I$
    \item $f_N =  -\beta I C_I$
\end{itemize}
The population payoff is
\begin{itemize}
    \item $\bar{f}(x)= xf_U+ (1-x)f_N$
\end{itemize}
We assume that individuals do not interact with other individuals’ strategies to determine their own cost, like it would be in a card game. To play against the entire infectious population is actually in accordance with the mass-action principle for disease transmission.
The evolution of user frequency $x$ satisfies the replicator equation:
\[\dot{x} = x [f_U - \bar{f}(x)]\]
In our case, this equation simplifies to
\[ \dot{x} = x [(1-x) f_U - (1-x) f_N] = x(1-x)(f_U-f_N) \]
also known as \HandRight\underline{logistic growth} equation. The net user behavior dynamics is governed by 
\begin{equation*}
    \dot{x} = x(1-x) \left[ \sigma\beta I C_I  - C_N \right]
\end{equation*}
\begin{itemize}
    \item \textbf{Adoption:} $x$ increases if $\beta I > \frac{C_N}{C_I \sigma}$. The risk-reduction benefit outweighs the nuisance cost.
    \item \textbf{Abandonment:} $x$ decreases if $\beta I < \frac{C_N}{C_I \sigma}$. The nuisance cost is perceived as greater than the expected cost of infection.
    \item \textbf{Evolutionary Stability:} A mixed equilibrium exists when $\beta I = \frac{C_N}{C_I \sigma}$.
\end{itemize}


\paragraph{Equilibrium for infectious disease dynamics}
To find the endemic equilibrium of the reduced system 
\begin{align*}
    \frac{dI_U}{dt} &= (1-\sigma) \beta (x-I_U) (I_U + I_N) - \gamma I_U \\
    \frac{dI_N}{dt} &= \beta (1-x-I_N) (I_U + I_N) - \gamma I_N
\end{align*}
we need to consider the case $R_0(x,\sigma)>1$ and solve $dI_U/dt = 0$ and $dI_N/dt = 0$.
Let $I = I_U + I_N$ be the total prevalence. $$I_U^* = \frac{(1-\sigma) \beta x I}{(1-\sigma) \beta I + \gamma}$$$$I_N^* = \frac{\beta (1-x) I}{\beta I + \gamma}$$ Since $I = I_U + I_N$, we sum the two equations:$$I = \left[ \frac{(1-\sigma) \beta x}{(1-\sigma) \beta I + \gamma} + \frac{\beta (1-x)}{\beta I + \gamma} \right] I$$ If $I > 0$, we divide both sides by $I$:$$1 = \frac{(1-\sigma) \beta x}{(1-\sigma) \beta I + \gamma} + \frac{\beta (1-x)}{\beta I + \gamma}$$This is a quadratic equation in $I$: 
$$\beta^2 (1-\sigma) I^2 + \beta \left[ \gamma(2-\sigma) - \beta(1 - \sigma ) \right] I + \gamma^2 \left[ 1 - R_0(x,\sigma) \right] = 0$$

The positive root of this quadratic gives the unique endemic equilibrium level when $R_0(x,\sigma ) > 1$:$$I^* = \frac{-B + \sqrt{B^2 - 4AC}}{2A}$$ for $A = \beta^2 (1-\sigma)$, $B = \beta \gamma (2-\sigma) - \beta^2(1 - \sigma )$ and $C = \gamma^2 (1 - R_0(x,\sigma))$.

\paragraph{Equilibrium for the coupled disease-replicator dynamics}
We can also solve the quadratic equation in $I$ for $x$ to obtain the user frequency at the endemic equilbrium. From the replicator equation we know that the user frequency is at equilibrium if the payoffs for users and non-users are equal: $f_U = f_N$ which amounts to mainting a specific prevalence $\bar{I}$:
   $$ \bar{I} = \frac{C_N}{\sigma\beta C_I }$$ in the population.
Combining both equations yields
$$\bar{x} = \frac{\gamma}{\beta \sigma} \left( \frac{\beta}{\gamma} - 1 \right) - \frac{[ \gamma(2-\sigma) - \beta(1 - \sigma ) ] }{\gamma \sigma} \bar{I}- \frac{\beta (1-\sigma) }{\gamma \sigma}\bar{I}^2$$
The strategy $\bar{x}$ is an \HandRight\underline{evolutionary stable} strategy for the case where the disease will persist and eventually reach an endemic equilbrium $\bar{I}$. The coupled disease-replicator dynamics is at equilibrium $(\bar{I},\bar{x})$, which depends on both disease or efficacy parameters $\beta, \gamma, \sigma$ and net usage behavior parameters $C_N$ and $C_I$.

\paragraph{Elimination threshold and eradication gap}
Let us recall the elimination threshold of user frequency $$x^e = \frac{1}{\sigma} \left( 1 - \frac{\gamma}{\beta} \right)$$
It follows that $$\bar{x} = x^e + \underbrace{\left[ \frac{\beta(1-\sigma) - \gamma(2-\sigma)}{\gamma \sigma} \right] \bar{I}}_{\text{Linear Interaction}} - \underbrace{\left[ \frac{\beta (1-\sigma)}{\gamma \sigma} \right] \bar{I}^2}_{\text{Quadratic Interaction}}$$
The difference
$$\Delta=\bar{x}-x^e$$ is called the eradication gap. As the nuisance cost $C_N$ increases, the gap is widening,  

Equality between the elimination threshold and the evolutionary stable strategy can only hold if $\bar{I}=0$, i.e. the nuisance cost of using a net is $C_N=0$. , and states that as long as nets are effective and have a non-zero nuisance cost, non-users will free-ride at the expense of users.

\paragraph{Heureka}
On the long term, the equilibrium $x^*$ is the optimum based on individual preference leading to endemic disease, whereas $x_e$ is the societal optimum leading to disease eradication. As efficacy $\sigma$ increases, $x_e$ decreases and $x^*$ increases, i.e. the eradication gap narrows. Nevertheless, the gap will never close.
