---
title: 'Simulate multivariate normal distribution with known covariance'
date: 2026-03-25
permalink: /posts/2025/07/blog-post-7-cellregmap-simulation/
tags:
  - statistics
---

-

Background
===

I recently came across this [paper](https://link.springer.com/article/10.15252/msb.202110663) at work, which introduces *CellRegMap*, a statistical method for identifying context‐specific regulatory variants using scRNA‐seq. 

The model can be specified as:

$$
y = g\beta_G + g \odot \beta_{G \times C} + u + c + e,
$$

where

- $$y$$ is the log-transformed single-cell expression for a given gene,  
- $$g$$ is the SNP genotype,  
- $$\beta_G$$ is the persistent genetic effect,  
- $$\beta_{G \times C} \sim \mathcal{N}(0, \sigma^2_{G \times C} \Sigma)$$ is the cell-specific GxC effect,  
- $$u \sim \mathcal{N}(0, \sigma^2_{RC} \, R \odot \Sigma)$$ accounts for repeat samples,  
- $$c \sim \mathcal{N}(0, \sigma^2_C \Sigma)$$ captures effects of cell context,  
- $$e \sim \mathcal{N}(0, \sigma^2_n I)$$ is the noise term,  
- $$ \Sigma$$ is the cellular context covariance, defined as $$ \Sigma = C C^T$$, where matrix of environmental contexts $$C$$ is column-normalized (mean = 0, standard deviation = 1). 

and $$\odot$$ denotes the element-wise Hadamard product. So far, everything is fairly standard. Nothing surprising.  


The authors then used following method to simulate a synthetic dataset to benchmark their method: 

$$
y \sim \text{Poisson}(\lambda), \quad 
\lambda = \exp\left(
y_{\text{base}} + \sum_{i=1}^{k} g \odot c_i (\beta_{G \times C})_i + g \beta_G
\right),
$$

where

- $$y_{\text{base}}$$ is the log-transformed observed (background) gene expression vector for a given gene in the reference dataset,  
- $$g$$ is the SNP genotype vector from the reference dataset,  
- $$c_i$$ denotes the $$i$$-th context variable (MOFA factor),  
- $$(\beta_{G \times C})_i \sim \mathcal{N}(0, \sigma_G^2 \rho_{G \times C})$$ is the interaction effect size for context $$i$$,  
- $$\beta_G \sim \mathcal{N}(0, \sigma_G^2 (1 - \rho_{G \times C}))$$ is the effect size of the persistent genetic effect, and  
- $$\sigma_G^2$$ denotes the total genetic variance and $$\rho_{G \times C}$$ is the fraction of genetic variance explained by GxC.

Why Does This Simulation Strategy Work?
===

What I found interesting is how the author simulated the cell-specific GxC effect $$\beta_{G \times C} $$, which is a n-variate (n=cell number) normal distribution that has a pre-defined covariance matrix, namely 

$$\beta_{G \times C} \sim \mathcal{N}(0, \sigma^2_{G \times C} \Sigma)$$

...where $$ \Sigma = C C^T$$.


To simulate $$\beta_{G \times C} $$, the authors started with univariate normal random variables i.i.d. as 

$$(\beta_{G \times C})_i \sim \mathcal{N}(0, \sigma_G^2 \rho_{G \times C})$$

They then constructed $$\beta_{G \times C} $$ used this formula:

$$\beta_{G \times C} = \sum_{i=1}^{k} c_i (\beta_{G \times C})_i $$


At first glance, this might look interesting. Why does this produce the desired covariance structure?


The Key Idea: Bilinearity of Covariance
===

It turns out that the authors used Bilinearity Property of Covariance. To demonstrate how this trick works, we first expand the expression $$\sum_{i=1}^{k} c_i (\beta_{G \times C})_i$$, and it gives:

$$
\begin{bmatrix}
c_1\beta_1 + c_2\beta_2 + \cdots + c_k\beta_k
\end{bmatrix}
$$

..which is a column vector of length n. 

And now consider two entries *i* and *j*, and calculate the covariance between *i-th* variate and *j-th* variate in this vector. So let:

$$
\begin{aligned}
A = c_{i1}\beta_1 + c_{i2}\beta_2 + \cdots + c_{ik}\beta_k \\
B = c_{j1}\beta_1 + c_{j2}\beta_2 + \cdots + c_{jk}\beta_k
\end{aligned}
$$

And through Bilinearity Property of Covariance in matrix form, we have 

$$
\begin{aligned}
A = \begin{bmatrix}
c_{i1} \quad c_{i2} \quad \cdots \quad c_{ik} \\
c_{j1} \quad c_{j2} \quad \cdots \quad c_{jk} 
\end{bmatrix} 

\\

A\beta = 
\begin{bmatrix}
a^T\beta \\
c^T\beta
\end{bmatrix}
\end{aligned}
$$

And we obtain

$$

Cov(a^T\beta, c^T\beta) = a^T\Sigma_{\beta}c = a_1c_1\sigma^2_G + a_2c_2\sigma^2_G + 
\cdots + a_kc_k\sigma^2_G

$$

...because $$\Sigma_{\beta}$$ is a covariance matrix, and $$\beta$$s are independent. 

Connecting back to $$ \Sigma = C C^T$$
===
Now lets turn back to $$\beta_{G \times C} \sim \mathcal{N}(0, \sigma^2_{G \times C} \Sigma)$$. It is easy to show that with $$ \Sigma = C C^T$$, the covariance of $$\beta_i$$ and $$\beta_j$$ (*i,j-th* entry in $$\sigma^2_{G \times C} \Sigma$$) is also 

$$ a_1c_1\sigma^2_G + a_2c_2\sigma^2_G + 
\cdots + a_kc_k\sigma^2_G $$

And so it works!