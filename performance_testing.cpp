#include <benchmark/benchmark.h>
#include <cmath>
#include <stan/math.hpp>

// See run_performnce_testing.sh

#ifndef BRANCH
#define BRANCH UNKNOWN_BRANCH
#endif

// Add branch name to benchmark output.

#define BENCHMARK_CAPTURE_(func, test_case_name, ...)                          \
  BENCHMARK_CAPTURE__(func, BRANCH test_case_name, __VA_ARGS__)

#define BENCHMARK_CAPTURE__(func, test_case_name, ...)                         \
  BENCHMARK_CAPTURE(func, test_case_name, __VA_ARGS__)

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
double do_calculation(const std::vector<double> &as,
                      const std::vector<double> &bs,
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
  double xxx = 0;
  xxx += maybe_do_a_grad(accum, as_T);
  xxx += maybe_do_a_grad(accum, bs_T);
  xxx += maybe_do_a_grad(accum, zs_T);
  stan::math::recover_memory_nested();
  return xxx;
}

template <typename T_location, typename T_precision>
double do_other_calculation(const std::vector<int> &n,
                            const std::vector<double> &mu,
                            const std::vector<double> &phi) {
  stan::math::start_nested();
  std::vector<T_location> mu_as_T_location(mu.begin(), mu.end());
  std::vector<T_precision> phi_as_T_precision(phi.begin(), phi.end());
  var lp =
      stan::math::neg_binomial_2_cdf(n, mu_as_T_location, phi_as_T_precision);
  double xxx = 0;
  xxx += maybe_do_a_grad(lp, mu_as_T_location);
  xxx += maybe_do_a_grad(lp, phi_as_T_precision);
  stan::math::recover_memory_nested();
  return xxx;
}

extern std::vector<double> as;
extern std::vector<double> bs;
extern std::vector<double> zs;

extern std::vector<int> n;
extern std::vector<double> mu;
extern std::vector<double> phi;

void test_do_calculation(benchmark::State &state) {
  const double result_ddd = do_calculation<double, double, double>(as, bs, zs);
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
  const double result_vdd = do_calculation<var, double, double>(as, bs, zs);
  const double result_dvd = do_calculation<double, var, double>(as, bs, zs);
  const double result_ddv = do_calculation<double, double, var>(as, bs, zs);
  assert(fabs(expected_result_vdd - result_vdd) < 1e-3);
  assert(fabs(expected_result_dvd - result_dvd) < 1e-3);
  assert(fabs(expected_result_ddv - result_ddv) < 1e-3);
}
BENCHMARK(test_do_calculation);

void test_do_other_calculation(benchmark::State &state) {
  const double result_dd = do_other_calculation<double, double>(n, mu, phi);
  const double expected_result_dd = 0;
  assert(expected_result_dd == result_dd);
  const double result_vd = do_other_calculation<var, double>(n, mu, phi);
  assert(result_vd != 0);
  assert(std::isfinite(result_vd));
  const double result_dv = do_other_calculation<double, var>(n, mu, phi);
  assert(result_dv != 0);
  assert(std::isfinite(result_dv));
}
BENCHMARK(test_do_other_calculation);

template <typename T_a, typename T_b, typename T_z>
void benchmark_inc_beta(benchmark::State &state, const T_a &, const T_b &,
                        const T_z &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_calculation<T_a, T_b, T_z>(as, bs, zs));
  }
}

template <typename T_mu, typename T_phi>
void benchmark_neg_binomial_2_cdf(benchmark::State &state, const T_mu &,
                                  const T_phi &) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(do_other_calculation<T_mu, T_phi>(n, mu, phi));
  }
}

// Would like to say:
//     BENCHMARK_CAPTURE(benchmark_inc_beta<v, v, v>(as, bs, zs))
// But we can't. So instead we say:
//     BENCHMARK_CAPTURE(benchmark_inc_beta, var_, var_, var_)
// Which delegates to:
//     benchmark_inc_beta<v, v, v>(as, bs, zs);
double double_;
var var_;

const int repetitions = 1;

BENCHMARK_CAPTURE_(benchmark_inc_beta, vvv, var_, var_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vvd, var_, var_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, ddv, double_, double_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vdv, var_, double_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, dvv, double_, var_, var_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, dvd, double_, var_, double_)
    ->Repetitions(repetitions);
BENCHMARK_CAPTURE_(benchmark_inc_beta, vdd, var_, double_, double_)
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
