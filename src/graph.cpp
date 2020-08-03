#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>

#include "graph.hpp"
#include "pcg_random.hpp"

using namespace std;

std::ostream& log_out = std::cout;

void
Graph::load(potts_model* model)
{
  params = model;
}

ostream&
Graph::print_distribution(ostream& os)
{
  vector<size_t> conf(n);

  double norm = 0;
  while (true) {
    double x = 0;
    for (size_t i = 0; i < n; ++i) {
      x += params->h(conf[i], i);
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        x += params->J(i, j)(conf[i], conf[j]);
      }
    }
    norm += exp(x);
    size_t j = 0;
    while (j < n && ++conf[j] == q) {
      conf[j] = 0;
      j++;
    }
    if (j == n) {
      break;
    }
  }
  while (true) {
    double x = 0;
    for (size_t i = 0; i < n; ++i) {
      x += params->h(conf[i], i);
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        x += params->J(i, j)(conf[i], conf[j]);
      }
    }
    os << "G2 " << exp(x) / norm << endl;
    size_t j = 0;
    while (j < n && ++conf[j] == q) {
      conf[j] = 0;
      j++;
    }
    if (j == n) {
      break;
    }
  }
  return os;
};

void
Graph::sample_mcmc_energies(double* ptr,
                            size_t m,
                            size_t mc_iters0,
                            size_t mc_iters,
                            long int seed,
                            double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf(i) = size_t(q * uniform(rng));
  }

  double E = 0;
  for (size_t i = 0; i < n; i++) {
    E -= params->h(conf(i), i);
    for (size_t j = i + 1; j < n; j++) {
      E -= params->J(i, j)(conf(i), conf(j));
    }
  }

  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h(q0, i);
    double e1 = -params->h(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J(j, i)(conf(j), q0);
        e1 -= params->J(j, i)(conf(j), q1);
      } else if (i < j) {
        e0 -= params->J(i, j)(q0, conf(j));
        e1 -= params->J(i, j)(q1, conf(j));
      }
    }
    double de = e1 - e0;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf(i) = q1;
      E += de;
    }
  }

  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * uniform(rng));
      size_t dq = 1 + size_t((q - 1) * uniform(rng));

      size_t q0 = conf(i);
      size_t q1 = (q0 + dq) % q;

      double e0 = -params->h(q0, i);
      double e1 = -params->h(q1, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e0 -= params->J(j, i)(conf(j), q0);
          e1 -= params->J(j, i)(conf(j), q1);
        } else if (i < j) {
          e0 -= params->J(i, j)(q0, conf(j));
          e1 -= params->J(i, j)(q1, conf(j));
        }
      }
      double de = e1 - e0;
      if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
        conf(i) = q1;
        E += de;
      }
    }
    *(ptr+s) = E;
  }
  return;
};

void
Graph::sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed,
                   double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf(i) = size_t(q * uniform(rng));
    assert(conf(i) < q);
  }

  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h(q0, i);
    double e1 = -params->h(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J(j, i)(conf(j), q0);
        e1 -= params->J(j, i)(conf(j), q1);
      } else if (i < j) {
        e0 -= params->J(i, j)(q0, conf(j));
        e1 -= params->J(i, j)(q1, conf(j));
      }
    }
    double de = e1 - e0;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf(i) = q1;
    }
  }

  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * uniform(rng));
      size_t dq = 1 + size_t((q - 1) * uniform(rng));

      size_t q0 = conf(i);
      size_t q1 = (q0 + dq) % q;

      double e0 = -params->h(q0, i);
      double e1 = -params->h(q1, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e0 -= params->J(j, i)(conf(j), q0);
          e1 -= params->J(j, i)(conf(j), q1);
        } else if (i < j) {
          e0 -= params->J(i, j)(q0, conf(j));
          e1 -= params->J(i, j)(q1, conf(j));
        }
      }
      double de = e1 - e0;
      if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
        conf(i) = q1;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr)(s, i) = conf(i);
    }
  }
  return;
};

void
Graph::sample_mcmc_init(arma::Mat<int>* ptr,
                        size_t m,
                        size_t mc_iters0,
                        size_t mc_iters,
                        arma::Col<int>* init_ptr,
                        long int seed,
                        double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf(i) = (*ptr)(i);
    assert(conf(i) < q);
  }

  for (size_t k = 0; k < mc_iters0; ++k) {
    size_t i = size_t(n * uniform(rng));
    size_t dq = 1 + size_t((q - 1) * uniform(rng));

    size_t q0 = conf(i);
    size_t q1 = (q0 + dq) % q;

    double e0 = -params->h(q0, i);
    double e1 = -params->h(q1, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J(j, i)(conf(j), q0);
        e1 -= params->J(j, i)(conf(j), q1);
      } else if (i < j) {
        e0 -= params->J(i, j)(q0, conf(j));
        e1 -= params->J(i, j)(q1, conf(j));
      }
    }
    double de = e1 - e0;
    if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
      conf(i) = q1;
    }
  }

  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      size_t i = size_t(n * uniform(rng));
      size_t dq = 1 + size_t((q - 1) * uniform(rng));

      size_t q0 = conf(i);
      size_t q1 = (q0 + dq) % q;

      double e0 = -params->h(q0, i);
      double e1 = -params->h(q1, i);
      for (size_t j = 0; j < n; ++j) {
        if (i > j) {
          e0 -= params->J(j, i)(conf(j), q0);
          e1 -= params->J(j, i)(conf(j), q1);
        } else if (i < j) {
          e0 -= params->J(i, j)(q0, conf(j));
          e1 -= params->J(i, j)(q1, conf(j));
        }
      }
      double de = e1 - e0;
      if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
        conf(i) = q1;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      (*ptr)(s, i) = conf(i);
    }
  }
  return;
};

void
Graph::sample_mcmc_zanella(arma::Mat<int>* ptr,
                           size_t m,
                           size_t mc_iters0,
                           size_t mc_iters,
                           long int seed,
                           std::string mode,
                           double temperature)
{
  pcg32 rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> uniform(0, 1);

  arma::Mat<double> de = arma::Mat<double>(n, q, arma::fill::zeros);
  arma::Mat<double> g = arma::Mat<double>(n, q, arma::fill::zeros);
  double lambda = 0.0;

  arma::Mat<double> de_old = arma::Mat<double>(n, q, arma::fill::zeros);
  arma::Mat<double> g_old = arma::Mat<double>(n, q, arma::fill::zeros);
  double lambda_old = 0.0;

  arma::Col<size_t> conf = arma::Col<size_t>(n);
  for (size_t i = 0; i < n; ++i) {
    conf(i) = size_t(q * uniform(rng));
    assert(conf(i) < q);
  }

  // Compute initial neighborhood.
  for (size_t i = 0; i < n; i++) {
    size_t q0 = conf(i);
    double e0 = -params->h(q0, i);
    for (size_t j = 0; j < n; ++j) {
      if (i > j) {
        e0 -= params->J(j, i)(conf(j), q0);
      } else if (i < j) {
        e0 -= params->J(i, j)(q0, conf(j));
      }
    }

    for (size_t q1 = 0; q1 < q; q1++) {
      if (q0 == q1) {
        de(i, q1) = 0.0;
      } else {
        double e1 = -params->h(q1, i);
        for (size_t j = 0; j < n; ++j) {
          if (i > j) {
            e1 -= params->J(j, i)(conf(j), q1);
          } else if (i < j) {
            e1 -= params->J(i, j)(q1, conf(j));
          }
        }
        de(i, q1) = e1 - e0;
      }
    }
  }

  if (mode == "sqrt") {
    g = arma::exp(de * -0.5 / temperature);
    lambda = arma::accu(g) - n; // n*exp(0) needs to be subtracted.
  } else if (mode == "barker") {
    g = 1.0 / (1.0 + arma::exp(de / temperature));
    lambda = arma::accu(g) - .5 * n;
  }

  de_old = de;
  g_old = g;
  lambda_old = lambda;

  for (size_t k = 0; k < mc_iters0; ++k) {
    double rand = uniform(rng) * lambda;
    double r_sum = 0.0;
    size_t i = 0;
    size_t q0 = 0;
    size_t q1 = 0;
    bool done = false;
    for (i = 0; i < n; i++) {
      for (q1 = 0; q1 < q; q1++) {
        if (conf(i) == q1) {
          continue;
        } else {
          r_sum += g(i, q1);
        }
        if (r_sum > rand) {
          done = true;
          break;
        }
      }
      if (done) {
        break;
      }
    }

    double tmp = de(i, q1);
    q0 = conf(i);
    conf(i) = q1;
    for (size_t pos = 0; pos < n; pos++) {
      for (size_t aa = 0; aa < q; aa++) {
        if (pos < i) {
          de(pos, aa) += params->J(pos, i)(conf(pos), q1) -
                         params->J(pos, i)(conf(pos), q0) -
                         params->J(pos, i)(aa, q1) + params->J(pos, i)(aa, q0);
        } else if (pos > i) {
          de(pos, aa) += params->J(i, pos)(q1, conf(pos)) -
                         params->J(i, pos)(q0, conf(pos)) -
                         params->J(i, pos)(q1, aa) + params->J(i, pos)(q0, aa);
        } else {
          if (q1 == aa) {
            de(pos, aa) = 0;
          } else if (q0 == aa) {
            de(pos, aa) = -tmp;
          } else {
            de(pos, aa) += params->h(q1, pos) - params->h(q0, pos);
            for (size_t pos2 = 0; pos2 < n; pos2++) {
              if (pos2 < i) {
                de(pos, aa) += params->J(pos2, i)(conf(pos2), q1) -
                               params->J(pos2, i)(conf(pos2), q0);
              } else if (pos2 > i) {
                de(pos, aa) += params->J(i, pos2)(q1, conf(pos2)) -
                               params->J(i, pos2)(q0, conf(pos2));
              }
            }
          }
        }
      }
    }

    if (mode == "sqrt") {
      g = arma::exp(de * -0.5 / temperature);
      lambda = arma::accu(g) - n; // n*exp(0) needs to be subtracted.
    } else if (mode == "barker") {
      g = 1.0 / (1.0 + arma::exp(de / temperature));
      lambda = arma::accu(g) - .5 * n;
    }

    if (uniform(rng) < lambda_old / lambda) {
      conf(i) = q1;
      de_old = de;
      g_old = g;
      lambda_old = lambda;
    } else {
      conf(i) = q0;
      g = g_old;
      lambda = lambda_old;
      de = de_old;
    }
  }
  for (size_t s = 0; s < m; ++s) {
    for (size_t k = 0; k < mc_iters; ++k) {
      double rand = uniform(rng) * lambda;
      double r_sum = 0.0;
      size_t i = 0;
      size_t q0 = 0;
      size_t q1 = 0;
      bool done = false;
      for (i = 0; i < n; i++) {
        for (q1 = 0; q1 < q; q1++) {
          if (conf(i) == q1) {
            continue;
          } else {
            r_sum += g(i, q1);
          }
          if (r_sum > rand) {
            done = true;
            break;
          }
        }
        if (done) {
          break;
        }
      }

      double tmp = de(i, q1);
      q0 = conf(i);
      conf(i) = q1;
      for (size_t pos = 0; pos < n; pos++) {
        for (size_t aa = 0; aa < q; aa++) {
          if (pos < i) {
            de(pos, aa) += params->J(pos, i)(conf(pos), q1) -
                           params->J(pos, i)(conf(pos), q0) -
                           params->J(pos, i)(aa, q1) +
                           params->J(pos, i)(aa, q0);
          } else if (pos > i) {
            de(pos, aa) += params->J(i, pos)(q1, conf(pos)) -
                           params->J(i, pos)(q0, conf(pos)) -
                           params->J(i, pos)(q1, aa) +
                           params->J(i, pos)(q0, aa);
          } else {
            if (q1 == aa) {
              de(pos, aa) = 0;
            } else if (q0 == aa) {
              de(pos, aa) = -tmp;
            } else {
              de(pos, aa) += params->h(q1, pos) - params->h(q0, pos);
              for (size_t pos2 = 0; pos2 < n; pos2++) {
                if (pos2 < i) {
                  de(pos, aa) += params->J(pos2, i)(conf(pos2), q1) -
                                 params->J(pos2, i)(conf(pos2), q0);
                } else if (pos2 > i) {
                  de(pos, aa) += params->J(i, pos2)(q1, conf(pos2)) -
                                 params->J(i, pos2)(q0, conf(pos2));
                }
              }
            }
          }
        }
      }

      if (mode == "sqrt") {
        g = arma::exp(de * -0.5 / temperature);
        lambda = arma::accu(g) - n; // n*exp(0) needs to be subtracted.
      } else if (mode == "barker") {
        g = 1.0 / (1.0 + arma::exp(de / temperature));
        lambda = arma::accu(g) - .5 * n;
      }

      if (uniform(rng) < lambda_old / lambda) {
        conf(i) = q1;
        de_old = de;
        g_old = g;
        lambda_old = lambda;
      } else {
        conf(i) = q0;
        de = de_old;
        g = g_old;
        lambda = lambda_old;
      }
    }

    for (size_t i = 0; i < n; ++i) {
      (*ptr)(s, i) = (int)conf(i);
    }
  }
  return;
};

ostream&
Graph::print_parameters(ostream& os)
{
  log_out << "printing parameters J" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      for (size_t yi = 0; yi < q; ++yi) {
        for (size_t yj = 0; yj < q; ++yj) {
          os << "J " << i << " " << j << " " << yi << " " << yj << " "
             << params->J(i, j)(yi, yj) << endl;
        }
      }
    }
  }
  log_out << "printing parameters h" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t yi = 0; yi < q; ++yi) {
      os << "h " << i << " " << yi << " " << params->h(yi, i) << endl;
    }
  }
  log_out << "done" << endl;
  return os;
};

void
Graph::print_parameters(FILE* of)
{
  log_out << "printing parameters J" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      for (size_t yi = 0; yi < q; ++yi) {
        for (size_t yj = 0; yj < q; ++yj) {
          fprintf(of,
                  "J %lu %lu %lu %lu %g\n",
                  i,
                  j,
                  yi,
                  yj,
                  params->J(i, j)(yi, yj));
        }
      }
    }
  }
  log_out << "printing parameters h" << endl;
  for (size_t i = 0; i < n; ++i) {
    for (size_t yi = 0; yi < q; ++yi) {
      fprintf(of, "h %lu %lu %g\n", i, yi, params->h(yi, i));
    }
  }
  log_out << "done" << endl;
};
