---
layout: post
title:  "On cross-fitting with plug-in estimators"
date:   2022-10-24 23:36:56 +0000
categories: semiparametric
---

<style>
.nm_conditions ol {
  type: "i"
}
</style>

<h3> Table of Contents </h3>
* toc
{:toc}

## Introduction
There's a strong correlation between cross-fitting and doubly robust methods in the recent literature on semiparametric statistics that can make it difficult for new researchers to easily understand the purpose of each of these. While extremely powerful when used together, they solve different "problems" that can come up with semiparametric two-stage methods.

In this post, I will discuss the role of cross-fitting, how it affects the asymptotic distribution of an estimator, and why, in some cases, it isn't necessary. Additionally, by comparing estimators with and without cross-fitting to sample-splitting estimators, we will more clearly understand the effect of nuisance parameter estimation on the asymptotic distribution of the estimator.

### Background

Here, we restrict our discussion to a specific subset of two-stage methods, where we first estimate some unknown (nonparametric) nuisance parameter $$ \gamma(\cdot) $$ with a function $$ \hat{\gamma}(x) $$ and then plug it in to a known function $$g(\cdot)$$ that we average,

$$ \hat{\theta} = \frac{1}{n} \sum_{i=1}^n g(Z_i, \hat{\gamma}). $$

Many of the ideas will apply to other semiparametric estimators, but this simplifies a lot of the algebra.

An evident issue with such an estimator is the re-use of data for fitting the nuisance parameter estimator $$ \hat{\gamma} $$ and evaluating it in $$ g(Z_i, \hat{\gamma}) $$. For example, suppose we wanted to estimate the average treatment effect of the treatment $$W$$ on the outcome $$Y$$ with observed confounders $$X$$ using the augmented inverse propensity weighting (AIPW) estimator, where

$$ g_{\text{AIPW}}(Z_i, \gamma) = \mu_1(X_i) - \mu_0(X_i) + \frac{W_i - \pi(X_i)}{\pi(X_i)(1 - \pi(X_i))}(Y_i - \mu_{W_i}(X_i)). $$

The residuals $$Y_i - \mu_{W_i}(X_i)$$ will be _in-sample_ residuals, instead of _out-of-sample_ residuals; It is well-known that in-sample residuals are often systematically smaller than out-of-sample residuals. However, the properties of out-of-sample residuals are important to giving the AIPW estimator the doubly-robust property. One solution to this problem is to use a technique called _cross-fitting_, which is closely related to _cross-validation_ in machine learning.

By cross-fitting, we mean the following procedure:
1. splitting the data into $$k$$ folds, and, 
2. iterating over folds $$ j=1,\dots,k $$, fitting $$ \hat{\gamma}^{j}(\cdot) $$ on the data from all *but* the $$ j $$-th fold, and
3. finally, computing the estimate as the average

$$ \hat{\theta}_{\text{cross-fit}} = \frac{1}{n} \sum_{j=1}^k \sum_{i~\text{in fold}~j} g(Z_i, \hat{\gamma}^{j}). $$

Cross-fitting has gained significant popularity in the semiparametric literature recently, due to it's importance for developing clean black-box master theorems for the use of machine learning in semiparametric estimation methods. While cross-fitting intuitively removes own-observation biases as in AIPW example above, it isn't always clear in these master theorems what other roles the cross-fitting plays in defining the statistical properties of the resulting method.

### Nuisance Parameter Estimates with Kernel Smoothing 

An interesting point of reference is a classical semiparametric method and theoretical analysis from the Handbook of Econometrics, Vol IV, Chp 36, by Newey and McFadden. In Section 8, they discuss a plug-in method as suggested above using kernel smoothing estimators to fit the nuisance parameters nonparametrically, using no cross-fitting. Kernel smoothing is a good nonparametric estimation technique when the function to be estimated has a similar level of smoothness all across the input space. Assuming this is true, and the function is sufficiently smooth (we will formalize this later), then cross-fitting is not necessary to get an asymptotic normality of the plug-in estimator.

<details><summary>Expand for details on kernel smoothing nuisance parameter estimation</summary>

</details>
<br />

The formal assumptions about the nuisance parameters and the bandwidth of the estimator to ensure than the estimated nuisance parameters converge sufficiently quickly are below. With these, we can restate Newey and McFadden's result for plug-in semiparametric estimators with kernel density nuisance parameter estimators. Here, we state a slightly more precise statement (asymptotic linearity) that is shown in the proof in Newey and McFadden, but not stated in the main theorem. The asymptotic linearity directly implies the asymptotic normality due to the central limit theorem.

<details><summary>Expand for assumptions for kernel smoothing (Assumptions 8.1-8.3)</summary>
Placeholder.
</details>
<br />

**Theorem 8.11 *(Newey and McFadden, 1994)***

*Suppose that Assumptions 8.1-8.3 are satisfied, $$ E[g(Z, \gamma_0)] = 0, $$ $$E[ \|g(Z, \gamma_0)\|^2 ] < \infty,$$ $$\mathcal{X}$$ is a compact set, the kernel bandwidth $$ \sigma $$ satisfies $$ n \sigma^{2r + 4d} / (\operatorname{ln} n)^2 \to \infty $$ and $$ n \sigma^{2m} \to 0, $$ and there is a vector of functionals $$ G(z, \gamma) $$ that is linear in $$ \gamma $$ such that*

  1. *(linearization) for all $$ \gamma $$ with $$ \gamma - \gamma_0 $$ small enough, $$ \| g(z, \gamma) - g(z, \gamma_0) - G(z, \gamma - \gamma_0) \| \le b(z) \| \gamma - \gamma_0 \|^2 $$, and $$ E[b(Z)] \sqrt{n} \| \hat{\gamma} - \gamma_0 \|^2 \stackrel{p}{\to} 0 $$;*
  2. *(bounded linear) $$ \vert G(z, \gamma)\vert \le c(z) \| \gamma \|, $$ with $$E[ c^2(Z) ] < \infty$$;*
  3. *(representation) there is $$v(x)$$ with $$\int G(z, \gamma) \text{d}P(z) = \int v(x) \gamma(x) \text{d}{x} $$ for all $$\| \gamma \| < \infty$$;*
  4. *$$ v(x) $$ is continuous almost everywhere, $$\int \|v(x)\| \text{d} x < \infty $$, and there is $$ \epsilon > 0 $$ such that $$ E[\sup_{\| \nu \| \le \epsilon} \| v(x + \nu) \|^4 ]. $$*

*Then, for $$ \delta(z) = v(x) y - E[v(X) Y], $$*

$$ \frac{1}{\sqrt{n}} \sum_{i=1}^n g(Z_i, \hat{\gamma}) = \frac{1}{\sqrt{n}} \sum_{i = 1}^n g(Z_i, \gamma_0) + \delta(Z_i) + o_P\left(1\right). $$

## Same statistical properties with cross-fitting

In the rest of this post, we will show that using cross-fitting versus not does not affect the asymptotic properties of the plug-in estimator with kernel-smoothing estimates of the nuisance parameters, by articulating a cross-fit alternative of this result.

### Cross-fit (and not) master theorems

We start with a more general master theorem (Theorem 8.1 in N+M) that can be used with any nuisance parameter estimator. As before, we state the asymptotic linearity result implied in the proof, but not stated in the main theorem.

**Theorem 8.1 *(Newey and McFadden, 1994)***

*Suppose that $$ E[g(Z, \gamma_0)] = 0,~E[\|g(Z, \gamma_0)\|^2 ] < \infty, $$ and there is $$ \delta(z) $$ with $$ E[\delta(Z)] = 0, $$ $$E[\|\delta(Z)\|^2] < \infty $$ and*
  1. *(linearization) there is a function $$ G(z, \gamma - \gamma_0) $$ that is linear is $$ \gamma - \gamma_0 $$ such that for all $$ \gamma $$ with $$ \gamma - \gamma_0 $$ small enough, $$ \| g(z, \gamma) - g(z, \gamma_0) - G(z, \gamma - \gamma_0) \| \le b(z) \| \gamma - \gamma_0 \|^2 $$, and $$ E[b(Z)] \sqrt{n} \| \hat{\gamma} - \gamma_0 \|^2 \stackrel{p}{\to} 0 $$;*
  2. *(stochastic equicontinuity) $$ \tfrac{1}{\sqrt{n}} \sum_{i=1}^n [G(Z_i, \hat{\gamma} - \gamma_0) - \int G(z, \hat{\gamma} - \gamma_0) \text{d}{P}(z)] \stackrel{p}{\to} 0 $$;*
  3. *(mean-square differentiability) there is a measure $$ \tilde{F} $$ such that for all $$ \| \gamma - \gamma_0 \| $$ small enough, $$ \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) = \int \delta(z) \text{d}\tilde{F}(z) $$;*
  4. *the measure $$ \tilde{F} $$ is approaches the empirical measure, in the sense that $$ \sqrt{n}[ \int \delta(z) \text{d}\tilde{F}(z) - \tfrac{1}{n} \sum_{i=1}^n \delta(Z_i) ] \stackrel{p}{\to} 0. $$*

*Then, $$ \tfrac{1}{\sqrt{n}} \sum_{i=1}^n g(Z_i, \hat{\gamma}) = \tfrac{1}{\sqrt{n}} \sum_{i=1}^n g(Z_i, \gamma_0) + \delta(Z_i) + o_{P}(1). $$*

Essentially, this result provides very general conditions, not specific to the choice of nuisance parameter estimator, under which we can *linearize* the first order effect of plugging in the estimated nuisance parameter $$ \hat{\gamma} $$ instead of the true nuisance parameter $$ \gamma_0 $$, as simply a sum of per-observation components $$ \delta(Z_i). $$

This result can easily be generalized to a cross-fitting estimator, and in the process we will simplify the *stochastic equicontinuity* condition into a simple condition that the linearization $$G(z, \gamma)$$ is has $$L_2$$-bounded coefficients.

The first step in doing so is articulating the result for a sample-splitting estimator. A sample plitting estimator is a simple form of cross-fitting where we split the data into two pieces, $$I_1$$ and $$I_2$$, each a constant fraction of the total sample size, $$ \vert I_1 \vert = a \cdot n$$, for $$ 0 < a < 1 $$, use $$I_2$$ to fit the nuisance parameters $$\hat{\gamma}$$, and use $$ I_1 $$ to construct the average, $$ \hat{\theta}_{\text{sample-split}} = \tfrac{1}{an} \sum_{i \in I_1} g(Z_i, \hat{\gamma}). $$

**Theorem 1**

*Suppose that $$ E[g(Z, \gamma_0)] = 0,~E[\|g(Z, \gamma_0)\|^2 ] < \infty, $$ and there is $$ \delta(z) $$ with $$ E[\delta(Z)] = 0, $$ $$E[\|\delta(Z)\|^2] < \infty $$ and*
  1. *(linearization) there is a function $$ G(z, \gamma - \gamma_0) $$ that is linear is $$ \gamma - \gamma_0 $$ such that for all $$ \gamma $$ with $$ \gamma - \gamma_0 $$ small enough, $$ \| g(z, \gamma) - g(z, \gamma_0) - G(z, \gamma - \gamma_0) \| \le b(z) \| \gamma - \gamma_0 \|^2 $$, and $$ E[b(Z)] \sqrt{n} \| \hat{\gamma} - \gamma_0 \|^2 \stackrel{p}{\to} 0 $$;*
  2. *(bounded linear) $$ \vert G(z, \gamma)\vert \le c(z) \| \gamma \|, $$ with $$E[ c^2(Z) ] < \infty$$;*
  3. *(mean-square differentiability) there is a measure $$ \tilde{F} $$ such that for all $$ \| \gamma - \gamma_0 \| $$ small enough, $$ \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) = \int \delta(z) \text{d}\tilde{F}(z) $$;*
  4. *the measure $$ \tilde{F} $$ is approaches the empirical measure over $$I_2$$, in the sense that $$ \sqrt{n}[ \int \delta(z) \text{d}\tilde{F}(z) - \tfrac{1}{(1-a)n} \sum_{j \in I_2} \delta(Z_j) ] \stackrel{p}{\to} 0. $$*

*Then, $$ \tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \hat{\gamma}) = \tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \gamma_0) + \tfrac{1}{(1-a) \sqrt{n}} \sum_{j \in I_2} \delta(Z_i) + o_{P}(1). $$*

<details><summary>Expand to see proof of Theorem 1.</summary>
<p><b><i>Proof.</i></b></p>
<p>
The proof here is very similar to the proof of Theorem 8.1 by Newey and McFadden (1994). However, we provide a bit more detail on the expansion here. The main idea is to take the terms in the sum and write them as an approximation that will be easier to analyze, and an approximation error term that we will show will be lower order. The first is to add and subtract a linearization of \( g(z, \gamma) \) as \( G(z, \gamma) \).
\[ \begin{aligned} \tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \hat{\gamma}) &= \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \gamma_0)}_{\text{Term 1}} \\
&~~~~~~~~~~~ + \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} (g(Z_i, \hat{\gamma}) - g(Z_i, \gamma_0)) - G(Z_i, \hat{\gamma} - \gamma_0)}_{\text{Term 2}} \\
&~~~~~~~~~~~ + \tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma} - \gamma_0) \end{aligned} \]
Then, we can look at the difference between this empirical average with the population average of the linearization, \( \int G(z, \hat\gamma - \gamma_0) \text{d}P(z), \) 
\[ \begin{aligned} \tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma} - \gamma_0) &= \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma} - \gamma_0) - \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z)}_{\text{Term 3}} \\
&~~~~~~~~~~~ + \sqrt{n} \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) \end{aligned} \]
Finally, using the assumptions 3 and 4 of the theorem, we can linearize the integral as
\[ \begin{aligned}
\sqrt{n} \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) &= \underbrace{\sqrt{n}\left(\int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) - \int \delta(z) \text{d}\tilde{F}(z)\right)}_{\text{Term 4}} \\
&~~~~~~~~~~~ + \underbrace{\sqrt{n} \left[\int \delta(z) \text{d}\tilde{F}(z) - \tfrac{1}{(1-a)n} \sum_{i \in I_2} \delta(Z_i)\right]}_{\text{Term 5}} \\
&~~~~~~~~~~~ + \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} \tfrac{1}{(1-a)n} \sum_{i \in I_2} \delta(Z_i)}_{\text{Term 6}} \end{aligned} \]
Putting all of these pieces together, we have
\[
\begin{aligned}
\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \hat{\gamma}) &= \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} g(Z_i, \gamma_0)}_{\text{Term 1}} \\
&~~~~~~~~~~~ + \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} (g(Z_i, \hat{\gamma}) - g(Z_i, \gamma_0)) - G(Z_i, \hat{\gamma} - \gamma_0)}_{\text{Term 2}} \\
&~~~~~~~~~~~ + \underbrace{\tfrac{1}{a \sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma} - \gamma_0) - \int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z)}_{\text{Term 3}} \\
&~~~~~~~~~~~ + \underbrace{\sqrt{n}\left(\int G(z, \hat{\gamma} - \gamma_0) \text{d}P(z) - \int \delta(z) \text{d}\tilde{F}(z)\right)}_{\text{Term 4}} \\
&~~~~~~~~~~~ + \underbrace{\sqrt{n} \left[\int \delta(z) \text{d}\tilde{F}(z) - \tfrac{1}{(1-a)n} \sum_{i \in I_2} \delta(Z_i)\right]}_{\text{Term 5}} \\
&~~~~~~~~~~~ + \underbrace{\tfrac{1}{(1-a)\sqrt{n}} \sum_{i \in I_2} \delta(Z_i)}_{\text{Term 6}}
\end{aligned}
\]
Terms 1 and 6 will remain as dominant terms, and the rest will be lower order. Terms 4 and 5 are lower order by assumptions 3 and 4 of the theorem, respectively. Therefore, to prove the result, it remains to show that terms 2 and 3 are lower order.
</p>

<p><b>Term 2</b></p>
<p>
Assumption 1 of the theorem and Chebyshev's inequality shows that this term converges to zero in probability--from the assumed convergence of \(\hat\gamma\), for some \(n_0\) large enough, we know that \( n \ge n_0\) implies \(\hat\gamma - \gamma_0\) is small enough such that the linearization holds. Therefore, for all \(n \ge n_0\),
\[ \begin{aligned} \left\vert \tfrac{1}{a\sqrt{n}} \sum_{i \in I_1} (g(Z_i, \hat{\gamma}) - g(Z_i, \gamma_0) - G(Z_i, \hat{\gamma} - \gamma_0)) \right\vert &\le \tfrac{1}{an} \sum_{i \in I_1} \sqrt{n} \vert g(Z_i, \hat{\gamma}) - g(Z_i, \gamma_0) - G(Z_i, \hat{\gamma} - \gamma_0) \vert \\
&\le \tfrac{1}{an} \sum_{i \in I_1}  \sqrt{n} b(Z_i)\| \hat\gamma - \gamma_0 \|^2. \end{aligned} \]
The weak law of large numbers implies that \( \tfrac{1}{an} \sum_{i \in I_1} b(Z_i) \stackrel{p}{\to} E[b(Z)], \) and by assumption, \( E[b(Z)] \sqrt{n} \| \hat\gamma - \gamma_0 \|^2 \stackrel{p}{\to} 0, \) so by the continuous mapping theorem, 
\[ \tfrac{1}{an} \sum_{i \in I_1}  \sqrt{n} b(Z_i)\| \hat\gamma - \gamma_0 \|^2 \stackrel{p}{\to} 0. \]
</p>

<p><b>Term 3</b></p>
<p>
This term is where sample splitting is crucial, and where our proof diverges from Newey and McFadden (1994). Given \(\hat{\gamma},\) the observations \( \{ G(Z_i, \hat{\gamma}-\gamma_0) - \int G(z, \hat{\gamma}-\gamma_0) \}_{i \in I_1} \) are independent and mean zero. Therefore, by Chebyshev's inequality,
\[ \begin{aligned} &P\left( \left\vert \frac{1}{a\sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma}-\gamma_0) - \int G(z, \hat{\gamma}-\gamma_0) \right\vert > \epsilon \middle\vert \hat\gamma\right) \\
&~~~~~~~~\le \frac{\operatorname{var}\left( \frac{1}{a\sqrt{n}} \sum_{i \in I_1} G(Z_i, \hat{\gamma}-\gamma_0) - \int G(z, \hat{\gamma}-\gamma_0) \middle\vert \hat\gamma \right)}{\epsilon^2}  \\
&~~~~~~~~\le \frac{\operatorname{var}\left( G(Z_i, \hat{\gamma}-\gamma_0) - \int G(z, \hat{\gamma}-\gamma_0) \middle\vert \hat\gamma \right)}{a \epsilon^2}  \\
&~~~~~~~~\le \frac{\operatorname{E}\left[ G^2(Z_i, \hat{\gamma}-\gamma_0) \middle\vert \hat\gamma \right]}{a \epsilon^2}.
\end{aligned} \]

By assumption 2, \(\operatorname{E}\left[ G^2(Z_i, \hat{\gamma}-\gamma_0 \mid \hat{\gamma}\right] \le \operatorname{E}\left[ c^2(Z)\| \hat{\gamma} - \gamma_0 \|^2 \mid \hat{\gamma}\right] = \operatorname{E}[ c^2(Z) ] \|\hat\gamma - \gamma_0\|^2. \)
By assumption, \(\|\hat\gamma - \gamma_0\|^2 \stackrel{p}{\to} 0\), so conditionally on \(\hat\gamma\), this sequence converges. Chernozhukov et al. (2018), Lemma 6.1, shows that conditional convergence implies unconditional convergence, and so\(\text{Term 3} \stackrel{p}{\to} 0.\)
<span align="right" style="display:block; clear: both;"> \(\Box\)</span>
</p>

</details>
<br />


This representation more clearly isolates the component of the influence function that comes from estimating the nuisance parameters, in that we can see that the $$ \delta(Z_j) $$ terms come from a sum over $$I_2,$$ where $$ \hat{\gamma} $$ was estimated.

From here, we can think of a cross-fit estimator as an average of sample-split estimators. We will use $$I_{1, j}$$ to denote the two pieces of each of the estimators, and in each one, $$a = 1/k.$$

$$ \begin{aligned} \hat{\theta}_{\text{cross-fit}} &= \frac{1}{k} \sum_{j=1}^k \frac{k}{n} \sum_{i \in I_{1, j}} g(Z_i, \hat{\gamma}^{j}) \\ &= \frac{1}{k} \sum_{j=1}^k \left( \frac{k}{n} \sum_{i \in I_{1, j}} g(Z_i, \gamma_0) + \frac{k}{(k - 1) n} \sum_{i \in I_{2, j}} \delta(Z_i) + o_P\left(\tfrac{1}{\sqrt{n}}\right) \right) \\ &= \frac{1}{n} \sum_{j=1}^k \sum_{i \in I_{1, j}} g(Z_i, \gamma_0) + \sum_{i \in I_{2, j}} \frac{1}{(k - 1)}\delta(Z_i) + o_P\left(\tfrac{1}{\sqrt{n}}\right) \end{aligned} $$

A given example $$i$$ appears in $$I_{1,j}$$ for just 1 value of $$j$$, and in $$I_{2,j}$$ for $$k-1$$ of them. Therefore, the sum $$ \sum_{j=1}^k \sum_{i \in I_{2, j}} \frac{1}{k-1} \delta(Z_i) = \sum_{i=1}^n \delta(Z_i). $$ Altogether, this careful accounting gives us the following corollary:

**Corollary 1**

*Under the assuptions of Theorem 1, the cross-fitting estimator has the following asymptotically linear representation:*

$$ \hat{\theta}_{\text{cross-fit}} = \frac{1}{n} \sum_{i = 1}^n g(Z_i, \gamma_0) + \delta(Z_i) + o_P\left(\tfrac{1}{\sqrt{n}}\right). $$

### Kernel smoothing and cross-fitting (analogue of Theorem 8.11)

With this, we can give an analogue of Theorem 8.11 from Newey and McFadden (1994) for when the nuisance parameters are estimated using kernel smoothing. Interestingly, the assumption replacing the stochastic equicontinuity assumption in our Theorem 1 is exactly the same as used in Assuption 8.? from Newey and McFadden, so there are no additional assumptions or conditions in this theorem.

**Theorem 2 *(cross-fitting analogue of Thm 8.11, Newey and McFadden, 1994)***

*Suppose that Assumptions 8.1-8.3 are satisfied, $$ E[g(Z, \gamma_0)] = 0, $$ $$E[ \|g(Z, \gamma_0)\|^2 ] < \infty,$$ $$\mathcal{X}$$ is a compact set, the kernel bandwidth $$ \sigma $$ satisfies $$ n \sigma^{2r + 4d} / (\operatorname{ln} n)^2 \to \infty $$ and $$ n \sigma^{2m} \to 0, $$ and there is a vector of functionals $$ G(z, \gamma) $$ that is linear in $$ \gamma $$ such that*

  1. *(linearization) for all $$ \gamma $$ with $$ \gamma - \gamma_0 $$ small enough, $$ \| g(z, \gamma) - g(z, \gamma_0) - G(z, \gamma - \gamma_0) \| \le b(z) \| \gamma - \gamma_0 \|^2 $$, and $$ E[b(Z)] \sqrt{n} \| \hat{\gamma} - \gamma_0 \|^2 \stackrel{p}{\to} 0 $$;*
  2. *(bounded linear) $$ \vert G(z, \gamma)\vert \le c(z) \| \gamma \|, $$ with $$E[ c^2(Z) ] < \infty$$;*
  3. *(representation) there is $$v(x)$$ with $$\int G(z, \gamma) \text{d}P(z) = \int v(x) \gamma(x) \text{d}{x} $$ for all $$\| \gamma \| < \infty$$;*
  4. *$$ v(x) $$ is continuous almost everywhere, $$\int \|v(x)\| \text{d} x < \infty $$, and there is $$ \epsilon > 0 $$ such that $$ E[\sup_{\| \nu \| \le \epsilon} \| v(x + \nu) \|^4 ]. $$*

*Then, for $$ \delta(z) = v(x) y - E[v(X) Y], $$*

$$ \frac{1}{\sqrt{n}} \sum_{j=1}^k \sum_{i~\text{in fold}~j} g(Z_i, \hat{\gamma}^{j}) = \frac{1}{\sqrt{n}} \sum_{i = 1}^n g(Z_i, \gamma_0) + \delta(Z_i) + o_P\left(1\right). $$

The proof is almost identical to the original proof of Theorem 8.11 in Newey and McFadden (1994), using Corollary 1 in place of their Theorem 8.1. Note that the second condition is the same between this result and that of Theorem 1 / Corollary 1, making the conditions of Theorem 2 identical to those assumed by Newey and McFadden. Therefore, we can see that cross-fitting is not necessary for kernel smoothing under the assumed conditions, and does not change the result. Therefore, any time a result in the literature cites Theorem 8.11 from N+M, a cross-fit version of the same estimator would work in the same place.

## References

  1. Newey, Whitney K., and Daniel McFadden. "Large sample estimation and hypothesis testing." Handbook of Econometrics, Vol 4 (1994): 2111-2245.
  2. Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, and James Robins. "Double/debiased machine learning for treatment and structural parameters." (2018): C1-C68.
