#pragma once

namespace MonteCarloIntegral {
    double integrate_gauss(double f(double), const double& start, const double& end, const unsigned long& samples);

    double gaussian(double mu, double sigma, double x);
}
