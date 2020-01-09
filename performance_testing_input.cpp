#include <vector>

std::vector<double> as{0.1, 0.5, 1, 10, 100};
std::vector<double> bs{0.1, 0.5, 1, 10, 100};
std::vector<double> zs{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

// neg_binomial_2
// E[Y] = mu ~ 50.
// phi[i] = mu[i]^2.
// sqrt[Var[Y]]] = sqrt[mu + mu^2/phi] ~ sqrt[51] ~ 7.
// Reasonable ns are like {30, ..., 70} or something.
std::vector<int> n{37, 40, 50, 50, 68};
std::vector<double> mu{40, 45, 50, 55, 60};
std::vector<double> phi{1600, 2025, 2500, 3025, 3600};
