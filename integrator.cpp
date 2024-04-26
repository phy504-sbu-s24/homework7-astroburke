// This program performs numerical integration. Edit the values in main() as desired
// Trapizoid method based on Professor Zingale's trapezoid integration, but the bug is fixed
// The bug was that in the integration function, "I" was not passed by reference, missing the &)

// By Aaron Burke for PHY 504 HW 7

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>

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

double f(double x)
{
    return x * x * x; // Function to integrate
}

int main()
{
    int N{100};       // Number of intervals
    double xmin{0.0}; // Lower Bound
    double xmax{6.0}; // Upper Bound
    double I{0.0};    // Initial value (probably 0)

    assert(N > 0);

    trapezoid_integrate(N, xmin, xmax, I, f);

    std::cout << "Value: " << std::setprecision(15) << I << std::endl;
}