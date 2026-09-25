# Monte Carlo Option Pricing in C++

A simple Monte Carlo implementation for pricing a European call option under the Black–Scholes model.

The program:

* simulates stock prices using Geometric Brownian Motion;
* estimates the option price from discounted payoffs;
* computes the Monte Carlo standard error;
* compares the result with the analytical Black–Scholes price;
* tests how the error changes with the number of simulations.

### Model

$$
dS_t = rS_t\,dt + \sigma S_t\,dW_t
$$

with

$$
S_T = S_0\exp\left((r-\frac12\sigma^2)T+\sigma W_T\right).
$$

The project is mainly a small exercise in C++, Monte Carlo methods, and quantitative finance.
