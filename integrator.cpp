// This program performs numerical integration. Edit the values in main() as desired.
// Trapizoid method based on Professor Zingale's trapezoid integration, but the bug is fixed.
// The bug was that in the integration function, "I" was not passed by reference, missing the &).
// Montecarlo method uses RNG as seed, then mersenee twister to generate random distribution.
// Main evaluates 2 different functions, then prints a table of their values for increasing samples

// By Aaron Burke for PHY 504 HW 7

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>
#include <random>

void trapezoid_integrate(int N, double xmin, double xmax, double &I, std::function<double(double)> f)
{
    // preliminaries
    double dx = (xmax - xmin) / N;
    double xleft{xmin};
    double fleft = f(xleft);
    double xright{0.0}; // to be calculated in loop
    double fright{0.0}; // to be calculated in loop

    // integration loop
    for (int n = 0; n <= N; n++)
    {
        xright = xmin + (n + 1) * dx;
        fright = f(xright);

        I += 0.5 * dx * (fleft + fright);

        xleft = xright;
        fleft = fright;
    }
}

void montecarlo_integrate(int N, double xmin, double xmax, double &I, std::function<double(double)> f)
{
    std::random_device rd;                                      // uniform RNG for seed
    std::mt19937 generator(rd());                               // mersenne twister takes seed
    std::uniform_real_distribution<double> uniform(xmin, xmax); // create random number distribution
    for (int n = 0; n < N; n++)                                 // Evaluate f at N random points
    {
        double x = uniform(generator);
        I += f(x);
    }
    I *= (xmax - xmin) / N; // integration
}

double f1(double x)
{
    return exp(-1.0 * (x * x)); // Function to integrate
}
double f2(double x)
{
    return sin(1.0 / (x * (2 - x) + 1e-12)) * sin(1.0 / (x * (2 - x) + 1e-12)); // Function to integrate
}

int main()
{
    int min_N = 8;
    int max_N = 1024;
    int factor = 2; // multiplication factor
    int precision = 8;
    int width = 12;

    std::cout << "Function 1: exp(-x^2)" << std::endl;
    std::cout << "N" << std::setw(width) << "Trapezoid" << std::setw(width) << "Monte Carlo" << std::endl;
    for (int N = min_N; N <= max_N; N *= factor)
    {
        double trapezoid_result1 = 0.0;
        trapezoid_integrate(N, -5.0, 5.0, trapezoid_result1, f1);

        double monte_carlo_result1 = 0.0;
        montecarlo_integrate(N, -5.0, 5.0, monte_carlo_result1, f1);

        std::cout << N << std::setw(width) << std::setprecision(precision) << trapezoid_result1 << std::setw(width) << monte_carlo_result1 << std::endl;
    }
    std::cout << "Function 2: sin^2((x(2-x)+1e-12)^-1)" << std::endl;
    std::cout << "N" << std::setw(width) << "Trapezoid" << std::setw(width) << "Monte Carlo" << std::endl;
    for (int N = min_N; N <= max_N; N *= factor)
    {
        double trapezoid_result2 = 0.0;
        trapezoid_integrate(N, 0.0, 2.0, trapezoid_result2, f2);

        double monte_carlo_result2 = 0.0;
        montecarlo_integrate(N, 0.0, 2.0, monte_carlo_result2, f2);

        std::cout << N << std::setw(width) << std::setprecision(precision) << trapezoid_result2 << std::setw(width) << monte_carlo_result2 << std::endl;
    }
}