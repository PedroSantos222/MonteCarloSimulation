#include <iostream>
#include <cmath>
#include<random>
#include<algorithm>
#include<vector>

//We assume the stock follows a geometric brownian motion under the risk neutral measure so the solution to the SDE is straightforward
//Estimating the price of the european call option today using the black scholes formulation

double normal_cdf(double x)
{
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double BlackScholesPrice(double S_0, double K, double T, double r, double sigma) {
    const double d1{ (std::log(S_0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T)) };
    const double d2{ d1 - sigma * std::sqrt(T) };

    const double C{ S_0 * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2) };

    return C;
}

std::vector<double> MonteCarloMethod(int N, double S_0, double K, double T, double r, double sigma, std::mt19937& generator ) {

    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> v{};

    std::vector<double> results{};

    double sum{ 0.0 };

    for (int i = 0;i < N;i++) {
        double z{ normal(generator) };
        double Wt{ std::sqrt(T) * z }; //Because we want the normal with variation T
        double sigmaexpWt{ std::exp(sigma * Wt) };
        double S_T{ S_0 * std::exp((r - 0.5 * sigma * sigma) * T) * sigmaexpWt };
        double M{ std::max(S_T - K ,0.0) };
        sum += M;
        v.push_back(M);
    }

    double mean{ sum / N };
    double result{ std::exp(-r * T) * mean };

    double sample_variance{};
    for (std::size_t i{ 0 };i < v.size();i++) {
        sample_variance += (std::exp(-r*T)*v[i] - result) * (std::exp(-r*T)*v[i] - result);
    }
    sample_variance /= (N - 1);

    double standard_deviation{ std::sqrt(sample_variance) };
    double standard_error{ standard_deviation / std::sqrt(N) };

    results.push_back(result);
    results.push_back(sample_variance);
    results.push_back(standard_error);

    return results; //the first value returns the expected value by the monte carlo method, the second value is the standard variance, the third value in the vector is the satndard error
}

int main()
{

    const double T{ 1.0 };
    const double S_0{ 50.0 };
    const double r{ 0.05 };
    const double sigma{ 0.20 };
    const double K{ 50.0 };




    for (int N{ 1000 }; N < 100000000; N *= 10) {
        std::mt19937 generator(50); //generate random numbers with the seed being 50

        std::vector<double> results{ MonteCarloMethod(N,S_0,K,T,r,sigma, generator)};

        double MC{ results[0]};

        //std::cout << "The expected value of the call option with maturity " << T << " and strike price " << K << " today is " << MC << "\n";

        double C{ BlackScholesPrice(S_0,K,T,r,sigma) };

        //std::cout << "The value of the call option given by the Black-Scholes formula is " << C << "\n";

        std::cout << "N=" << N << ", MC=" << MC << " , standard error=" << results[2] << " , error=" << std::abs(C - MC) << "\n";
    }
    return 0;
}

