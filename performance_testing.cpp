#include <benchmark/benchmark.h>
#include <cmath>
#include <stan/math.hpp>

// Add branch name to output.

#ifndef BRANCH
#define BRANCH UNKNOWN_BRANCH
#endif
#define STRINGIFY_(name) #name
#define STRINGIFY(name) STRINGIFY_(name)
#define BRANCHNAME STRINGIFY(BRANCH)

#define BENCHMARK_CAPTURE_(func, test_case_name, ...)                          \
  BENCHMARK_CAPTURE__(func, BRANCH test_case_name, __VA_ARGS__)

#define BENCHMARK_CAPTURE__(func, test_case_name, ...)                         \
  BENCHMARK_CAPTURE(func, test_case_name, __VA_ARGS__)

using dummy = bool;
dummy dummy_;

#define BENCHMARK_(func) BENCHMARK_CAPTURE_(func, x, dummy_)

using stan::math::var;

template <typename T> double maybe_do_a_grad(var &, std::vector<T> &);

double maybe_do_a_grad(var &dependent, std::vector<var> &independents) {
  using stan::math::set_zero_all_adjoints;
  stan::math::set_zero_all_adjoints();
  std::vector<double> grad;
  dependent.grad(independents, grad);
  return stan::math::sum(grad);
}

double maybe_do_a_grad(var &dependent, std::vector<double> &independents) {
  return 0;
}

template <typename T_a, typename T_b, typename T_z>
double do_inc_beta(const std::vector<double> &as, const std::vector<double> &bs,
                   const std::vector<double> &zs) {
  stan::math::start_nested();
  std::vector<T_a> as_T(as.begin(), as.end());
  std::vector<T_b> bs_T(bs.begin(), bs.end());
  std::vector<T_z> zs_T(zs.begin(), zs.end());
  var accum = 0;
  for (T_a &a : as_T) {
    for (T_b &b : bs_T) {
      for (T_z &z : zs_T) {
        accum += stan::math::inc_beta(a, b, z);
      }
    }
  }
  double accum2 = 0;
  accum2 += maybe_do_a_grad(accum, as_T);
  accum2 += maybe_do_a_grad(accum, bs_T);
  accum2 += maybe_do_a_grad(accum, zs_T);
  stan::math::recover_memory_nested();
  return accum2;
}

template <typename T_location, typename T_precision>
double do_neg_binomial_2_cdf(const std::vector<int> &n,
                             const std::vector<double> &mu,
                             const std::vector<double> &phi) {
  stan::math::start_nested();
  std::vector<T_location> mu_as_T_location(mu.begin(), mu.end());
  std::vector<T_precision> phi_as_T_precision(phi.begin(), phi.end());
  var lp =
      stan::math::neg_binomial_2_cdf(n, mu_as_T_location, phi_as_T_precision);
  double accum2 = 0;
  accum2 += maybe_do_a_grad(lp, mu_as_T_location);
  accum2 += maybe_do_a_grad(lp, phi_as_T_precision);
  stan::math::recover_memory_nested();
  return accum2;
}

double do_grad_reg_inc_beta(const std::vector<double> &as,
                            const std::vector<double> &bs,
                            const std::vector<double> &zs) {
  double accum = 0;
  for (const double &a : as) {
    for (const double &b : bs) {
      for (const double &z : zs) {
        double digamma_a = stan::math::digamma(a);
        double digamma_b = stan::math::digamma(b);
        double digamma_ab = stan::math::digamma(a + b);
        double beta_ab = stan::math::beta(a, b);
        double d_a = 0;
        double d_b = 0;
        stan::math::grad_reg_inc_beta(d_a, d_b, a, b, z, digamma_a, digamma_b,
                                      digamma_ab, beta_ab);
        accum += d_a;
        accum += d_b;
      }
    }
  }
  return accum;
}

double do_inc_beta_dda_ddb(const std::vector<double> &as,
                           const std::vector<double> &bs,
                           const std::vector<double> &zs) {
  double accum = 0;
  for (const double &a : as) {
    for (const double &b : bs) {
      for (const double &z : zs) {
        double digamma_a = stan::math::digamma(a);
        double digamma_b = stan::math::digamma(b);
        double digamma_ab = stan::math::digamma(a + b);
        accum += stan::math::inc_beta_dda(a, b, z, digamma_a, digamma_ab);
        accum += stan::math::inc_beta_ddb(a, b, z, digamma_b, digamma_ab);
      }
    }
  }
  return accum;
}

extern std::vector<double> as;
extern std::vector<double> bs;
extern std::vector<double> zs;

extern std::vector<int> n;
extern std::vector<double> mu;
extern std::vector<double> phi;

void test_do_inc_beta(benchmark::State &, dummy &) {
  const double result_ddd = do_inc_beta<double, double, double>(as, bs, zs);
  const double expected_result_ddd = 0;
  assert(expected_result_ddd == result_ddd);
  double expected_result_vdd = 0;
  double expected_result_dvd = 0;
  double expected_result_ddv = 0;
  for (const double &a : as) {
    for (const double &b : bs) {
      for (const double &z : zs) {
        double digamma_a = stan::math::digamma(a);
        double digamma_b = stan::math::digamma(b);
        double digamma_ab = stan::math::digamma(a + b);
        expected_result_vdd +=
            stan::math::inc_beta_dda(a, b, z, digamma_a, digamma_ab);
        expected_result_dvd +=
            stan::math::inc_beta_ddb(a, b, z, digamma_b, digamma_ab);
        expected_result_ddv += stan::math::inc_beta_ddz(a, b, z);
      }
    }
  }
  const double result_vdd = do_inc_beta<var, double, double>(as, bs, zs);
  const double result_dvd = do_inc_beta<double, var, double>(as, bs, zs);
  const double result_ddv = do_inc_beta<double, double, var>(as, bs, zs);
  assert(fabs(expected_result_vdd - result_vdd) < 1e-3);
  assert(fabs(expected_result_dvd - result_dvd) < 1e-3);
  assert(fabs(expected_result_ddv - result_ddv) < 1e-3);
}
BENCHMARK_(test_do_inc_beta);

void test_do_neg_binomial_2_cdf(benchmark::State &, dummy &) {
  const double result_dd = do_neg_binomial_2_cdf<double, double>(n, mu, phi);
  const double expected_result_dd = 0;
  assert(expected_result_dd == result_dd);
  const double result_vd = do_neg_binomial_2_cdf<var, double>(n, mu, phi);
  assert(result_vd != 0);
  assert(std::isfinite(result_vd));
  const double result_dv = do_neg_binomial_2_cdf<double, var>(n, mu, phi);
  assert(result_dv != 0);
  assert(std::isfinite(result_dv));
}
BENCHMARK_(test_do_neg_binomial_2_cdf);

void test_grad_vs_dda_ddb(benchmark::State &, dummy &) {
  const double result1 = do_grad_reg_inc_beta(as, bs, zs);
  const double result2 = do_inc_beta_dda_ddb(as, bs, zs);
  std::cout << result1 << " vs " << result2 << std::endl;
  assert(fabs(result1 - result2) < 1e-3);
}
BENCHMARK_(test_grad_vs_dda_ddb);

double sum_inc_beta_along_lattice(const std::vector<double> &as,
                                  const std::vector<double> &bs,
                                  const std::vector<double> &zs) {
  double accum = 0;
  accum += do_inc_beta<double, double, double>(as, bs, zs);
  accum += do_inc_beta<var, double, double>(as, bs, zs);
  accum += do_inc_beta<double, var, double>(as, bs, zs);
  accum += do_inc_beta<double, double, var>(as, bs, zs);
  accum += do_inc_beta<var, var, double>(as, bs, zs);
  accum += do_inc_beta<double, var, var>(as, bs, zs);
  accum += do_inc_beta<var, double, var>(as, bs, zs);
  accum += do_inc_beta<var, var, var>(as, bs, zs);
  return accum;
}

static bool once = true;
void eval_inc_beta_print_output(benchmark::State &state, dummy &) {
  std::vector<double> as{0.1, 0.5, 1, 10, 100};
  std::vector<double> bs{0.1, 0.5, 1, 10, 100};
  std::vector<double> zs{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  double accum = 0;
  accum += sum_inc_beta_along_lattice(as, bs, zs);
  accum += sum_inc_beta_along_lattice(as, bs, zs);
  // "Benchmark" executed many times, only want 1 copy of the output. XXX.
  if (once)
    std::cout << "On " << BRANCHNAME << " we get " << accum << "." << std::endl;
  accum = 0;
  std::vector<double> as_with_large_a = as;
  as_with_large_a.push_back(15000);
  accum += sum_inc_beta_along_lattice(as_with_large_a, bs, zs);
  std::vector<double> bs_with_large_b = bs;
  bs_with_large_b.push_back(15000);
  accum += sum_inc_beta_along_lattice(as, bs_with_large_b, zs);
  if (once)
    std::cout << "Including large a and b on " << BRANCHNAME << " we get "
              << accum << "." << std::endl;
  once = false;
}
BENCHMARK_(eval_inc_beta_print_output);

template <typename T_a, typename T_b, typename T_z>
void benchmark_inc_beta(benchmark::State &state, const T_a &, const T_b &,
                        const T_z &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_inc_beta<T_a, T_b, T_z>(as, bs, zs));
  }
}

template <typename T_mu, typename T_phi>
void benchmark_neg_binomial_2_cdf(benchmark::State &state, const T_mu &,
                                  const T_phi &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_neg_binomial_2_cdf<T_mu, T_phi>(n, mu, phi));
  }
}

void benchmark_grad_reg_inc_beta(benchmark::State &state, dummy &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_grad_reg_inc_beta(as, bs, zs));
  }
}
BENCHMARK_(benchmark_grad_reg_inc_beta);

void benchmark_inc_beta_dda_ddb(benchmark::State &state, dummy &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_inc_beta_dda_ddb(as, bs, zs));
  }
}
BENCHMARK_(benchmark_inc_beta_dda_ddb);

// Would like to say:
//     BENCHMARK_CAPTURE(benchmark_inc_beta<v, v, v>(as, bs, zs))
// But we can't. So instead we say:
//     BENCHMARK_CAPTURE(benchmark_inc_beta, var_, var_, var_)
// Which delegates to:
//     benchmark_inc_beta<v, v, v>(as, bs, zs);
double double_;
var var_;

const int repetitions = 3;

BENCHMARK_CAPTURE_(benchmark_inc_beta, vvv, var_, var_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vdv, var_, double_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, dvv, double_, var_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vvd, var_, var_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, dvd, double_, var_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vdd, var_, double_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, ddv, double_, double_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, ddd, double_, double_, double_)
    ->Repetitions(repetitions);

BENCHMARK_CAPTURE_(benchmark_neg_binomial_2_cdf, vv, var_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_neg_binomial_2_cdf, dv, double_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_neg_binomial_2_cdf, vd, var_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_neg_binomial_2_cdf, dd, double_, double_)
    ->Repetitions(repetitions);

BENCHMARK_MAIN();
