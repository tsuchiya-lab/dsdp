## (e.g. [@nesterov2000squared; @fushiki2006maximum])

## Let $p(x; \alpha)$ be a univariate k-th degree polynomial with coefficient $\alpha$.

## Then, there exists a unique quadratic form corresponding to $\alpha$

## in each of the following cases.

## 
## 1. $p(x; \alpha) \ge 0$ over $S = (-\infty, \infty)$

##   + $k$ is even and $d = \frac{k}{2}$

##     $$

##     p(x; \alpha) = \textbf{x}_{d}^T Q \textbf{x}_{d} = \mathrm{trace}(X_d Q),

##     $$

##     for some unique symmetric matrix $Q \succeq 0$.

## 2. $p(x; \alpha) \ge 0$ over $S = [0, \infty)$

##   + $k$ is odd and $d = \frac{k-1}{2}$

##     $$

##     p(x;\alpha) = \textbf{x}_d^T Q_1 \textbf{x}_d + x \cdot \textbf{x}_d^T Q_2

##     \textbf{x}_d = \mathrm{trace}(X_d Q_1) + \mathrm{trace}(x X_d Q_2)

##     $$

##     for some unique symmetric matrices $Q_1, Q_2 \succeq 0$.

##   + $k$ is even and $d = \frac{k}{2}$

##     $$

##     p(x;\alpha) = \textbf{x}_d^T Q_1 \textbf{x}_d + x \cdot \textbf{x}_{d-1}^T Q_3

##     \textbf{x}_{d-1} = \mathrm{trace}(X_d Q_1) + \mathrm{trace}(x X_{d-1} Q_3),

##     $$

##     for some unique symmetric matrices $Q_1, Q_3 \succeq 0$.

## 
