#include "sampler.hpp"

#include <string>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "pcg_random.hpp"

Sampler::Sampler(size_t N, size_t Q)
  : N(N)
  , Q(Q){};

Sampler::Sampler(size_t N, size_t Q, potts_model* m)
  : N(N)
  , Q(Q)
  , model(m){};

void
Sampler::setModel(potts_model* m)
{
  model = m;
};

void
Sampler::sampleEnergies(arma::Mat<double>* p,
                        size_t walkers,
                        size_t samples_per_walk,
                        size_t burn_in,
                        size_t burn_between,
                        long int seed,
                        double temperature)
{
#pragma omp parallel
  {
#pragma omp for
    for (size_t walker = 0; walker < walkers; walker++) {
      double* ptr = ((*p).colptr(walker));

      pcg32 rng;
      rng.seed(seed + walker);
      std::uniform_real_distribution<double> uniform(0, 1);

      arma::Col<size_t> conf = arma::Col<size_t>(N);
      for (size_t i = 0; i < N; ++i) {
        conf(i) = size_t(Q * uniform(rng));
      }

      double E = 0;
      for (size_t i = 0; i < N; i++) {
        E -= model->h(conf(i), i);
        for (size_t j = i + 1; j < N; j++) {
          E -= model->J(i, j)(conf(i), conf(j));
        }
      }

      for (size_t k = 0; k < burn_in; ++k) {
        size_t i = size_t(N * uniform(rng));
        size_t dq = 1 + size_t((Q - 1) * uniform(rng));

        size_t q0 = conf(i);
        size_t q1 = (q0 + dq) % Q;

        double e0 = -model->h(q0, i);
        double e1 = -model->h(q1, i);
        for (size_t j = 0; j < N; ++j) {
          if (i > j) {
            e0 -= model->J(j, i)(conf(j), q0);
            e1 -= model->J(j, i)(conf(j), q1);
          } else if (i < j) {
            e0 -= model->J(i, j)(q0, conf(j));
            e1 -= model->J(i, j)(q1, conf(j));
          }
        }
        double de = e1 - e0;
        if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
          conf(i) = q1;
          E += de;
        }
      }

      for (size_t s = 0; s < samples_per_walk; ++s) {
        for (size_t k = 0; k < burn_between; ++k) {
          size_t i = size_t(N * uniform(rng));
          size_t dq = 1 + size_t((Q - 1) * uniform(rng));

          size_t q0 = conf(i);
          size_t q1 = (q0 + dq) % Q;

          double e0 = -model->h(q0, i);
          double e1 = -model->h(q1, i);
          for (size_t j = 0; j < N; ++j) {
            if (i > j) {
              e0 -= model->J(j, i)(conf(j), q0);
              e1 -= model->J(j, i)(conf(j), q1);
            } else if (i < j) {
              e0 -= model->J(i, j)(q0, conf(j));
              e1 -= model->J(i, j)(q1, conf(j));
            }
          }
          double de = e1 - e0;
          if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
            conf(i) = q1;
            E += de;
          }
        }
        *(ptr + s) = E;
      }
    }
  }
  return;
};

void
Sampler::sampleSequences(arma::Cube<int>* p,
                         size_t walkers,
                         size_t samples_per_walk,
                         size_t burn_in,
                         size_t burn_between,
                         long int seed,
                         double temperature)
{
#pragma omp parallel
  {
#pragma omp for
    for (size_t walker = 0; walker < walkers; walker++){
      arma::Mat<int>* ptr = (arma::Mat<int>*)&((*p).slice(walker));

      pcg32 rng;
      rng.seed(seed + walker);
      std::uniform_real_distribution<double> uniform(0, 1);

      arma::Col<size_t> conf = arma::Col<size_t>(N);
      for (size_t i = 0; i < N; ++i) {
        conf(i) = size_t(Q * uniform(rng));
        assert(conf(i) < Q);
      }

      for (size_t k = 0; k < burn_in; ++k) {
        size_t i = size_t(N * uniform(rng));
        size_t dq = 1 + size_t((Q - 1) * uniform(rng));

        size_t q0 = conf(i);
        size_t q1 = (q0 + dq) % Q;

        double e0 = -model->h(q0, i);
        double e1 = -model->h(q1, i);
        for (size_t j = 0; j < N; ++j) {
          if (i > j) {
            e0 -= model->J(j, i)(conf(j), q0);
            e1 -= model->J(j, i)(conf(j), q1);
          } else if (i < j) {
            e0 -= model->J(i, j)(q0, conf(j));
            e1 -= model->J(i, j)(q1, conf(j));
          }
        }
        double de = e1 - e0;
        if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
          conf(i) = q1;
        }
      }

      for (size_t s = 0; s < samples_per_walk; ++s) {
        for (size_t k = 0; k < burn_between; ++k) {
          size_t i = size_t(N * uniform(rng));
          size_t dq = 1 + size_t((Q - 1) * uniform(rng));

          size_t q0 = conf(i);
          size_t q1 = (q0 + dq) % Q;

          double e0 = -model->h(q0, i);
          double e1 = -model->h(q1, i);
          for (size_t j = 0; j < N; ++j) {
            if (i > j) {
              e0 -= model->J(j, i)(conf(j), q0);
              e1 -= model->J(j, i)(conf(j), q1);
            } else if (i < j) {
              e0 -= model->J(i, j)(q0, conf(j));
              e1 -= model->J(i, j)(q1, conf(j));
            }
          }
          double de = e1 - e0;
          if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
            conf(i) = q1;
          }
        }

        for (size_t i = 0; i < N; ++i) {
          (*ptr)(s, i) = conf(i);
        }
      }
    }
  }
  return;
};

void
Sampler::sampleSequences(arma::Mat<int>* ptr,
                         size_t walkers,
                         size_t burn_in,
                         long int seed,
                         double temperature)
{
#pragma omp parallel
  {
#pragma omp for
    for (size_t walker = 0; walker < walkers; walker++){
      pcg32 rng;
      rng.seed(seed + walker);
      std::uniform_real_distribution<double> uniform(0, 1);

      arma::Col<size_t> conf = arma::Col<size_t>(N);
      for (size_t i = 0; i < N; ++i) {
        conf(i) = size_t(Q * uniform(rng));
        assert(conf(i) < Q);
      }

      for (size_t k = 0; k < burn_in; ++k) {
        size_t i = size_t(N * uniform(rng));
        size_t dq = 1 + size_t((Q - 1) * uniform(rng));

        size_t q0 = conf(i);
        size_t q1 = (q0 + dq) % Q;

        double e0 = -model->h(q0, i);
        double e1 = -model->h(q1, i);
        for (size_t j = 0; j < N; ++j) {
          if (i > j) {
            e0 -= model->J(j, i)(conf(j), q0);
            e1 -= model->J(j, i)(conf(j), q1);
          } else if (i < j) {
            e0 -= model->J(i, j)(q0, conf(j));
            e1 -= model->J(i, j)(q1, conf(j));
          }
        }
        double de = e1 - e0;
        if ((de < 0) || (uniform(rng) < exp(-de / temperature))) {
          conf(i) = q1;
        }
      }
      for (size_t i = 0; i < N; ++i) {
        (*ptr)(walker, i) = conf(i);
      }
    }
  }
  return;
};

void
Sampler::sampleSequencesZanella(arma::Cube<int>* p,
                                size_t walkers,
                                size_t samples_per_walk,
                                size_t burn_in,
                                size_t burn_between,
                                long int seed,
                                std::string mode,
                                double temperature)
{
#pragma omp parallel
  {
#pragma omp for
    for (size_t walker = 0; walker < walkers; walker++){
      arma::Mat<int>* ptr = (arma::Mat<int>*)&((*p).slice(walker));

      pcg32 rng;
      rng.seed(seed + walker);
      std::uniform_real_distribution<double> uniform(0, 1);

      arma::Mat<double> de = arma::Mat<double>(N, Q, arma::fill::zeros);
      arma::Mat<double> g = arma::Mat<double>(N, Q, arma::fill::zeros);
      double lambda = 0.0;

      arma::Mat<double> de_old = arma::Mat<double>(N, Q, arma::fill::zeros);
      arma::Mat<double> g_old = arma::Mat<double>(N, Q, arma::fill::zeros);
      double lambda_old = 0.0;

      arma::Col<size_t> conf = arma::Col<size_t>(N);
      for (size_t i = 0; i < N; ++i) {
        conf(i) = size_t(Q * uniform(rng));
        assert(conf(i) < Q);
      }

      // Compute initial neighborhood.
      for (size_t i = 0; i < N; i++) {
        size_t q0 = conf(i);
        double e0 = -model->h(q0, i);
        for (size_t j = 0; j < N; ++j) {
          if (i > j) {
            e0 -= model->J(j, i)(conf(j), q0);
          } else if (i < j) {
            e0 -= model->J(i, j)(q0, conf(j));
          }
        }

        for (size_t q1 = 0; q1 < Q; q1++) {
          if (q0 == q1) {
            de(i, q1) = 0.0;
          } else {
            double e1 = -model->h(q1, i);
            for (size_t j = 0; j < N; ++j) {
              if (i > j) {
                e1 -= model->J(j, i)(conf(j), q1);
              } else if (i < j) {
                e1 -= model->J(i, j)(q1, conf(j));
              }
            }
            de(i, q1) = e1 - e0;
          }
        }
      }

      if (mode == "sqrt") {
        g = arma::exp(de * -0.5 / temperature);
        lambda = arma::accu(g) - N; // n*exp(0) needs to be subtracted.
      } else if (mode == "barker") {
        g = 1.0 / (1.0 + arma::exp(de / temperature));
        lambda = arma::accu(g) - .5 * N;
      }

      de_old = de;
      g_old = g;
      lambda_old = lambda;

      for (size_t k = 0; k < burn_in; ++k) {
        double rand = uniform(rng) * lambda;
        double r_sum = 0.0;
        size_t i = 0;
        size_t q0 = 0;
        size_t q1 = 0;
        bool done = false;
        for (i = 0; i < N; i++) {
          for (q1 = 0; q1 < Q; q1++) {
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
        for (size_t pos = 0; pos < N; pos++) {
          for (size_t aa = 0; aa < Q; aa++) {
            if (pos < i) {
              de(pos, aa) += model->J(pos, i)(conf(pos), q1) -
                             model->J(pos, i)(conf(pos), q0) -
                             model->J(pos, i)(aa, q1) + model->J(pos, i)(aa, q0);
            } else if (pos > i) {
              de(pos, aa) += model->J(i, pos)(q1, conf(pos)) -
                             model->J(i, pos)(q0, conf(pos)) -
                             model->J(i, pos)(q1, aa) + model->J(i, pos)(q0, aa);
            } else {
              if (q1 == aa) {
                de(pos, aa) = 0;
              } else if (q0 == aa) {
                de(pos, aa) = -tmp;
              } else {
                de(pos, aa) += model->h(q1, pos) - model->h(q0, pos);
                for (size_t pos2 = 0; pos2 < N; pos2++) {
                  if (pos2 < i) {
                    de(pos, aa) += model->J(pos2, i)(conf(pos2), q1) -
                                   model->J(pos2, i)(conf(pos2), q0);
                  } else if (pos2 > i) {
                    de(pos, aa) += model->J(i, pos2)(q1, conf(pos2)) -
                                   model->J(i, pos2)(q0, conf(pos2));
                  }
                }
              }
            }
          }
        }

        if (mode == "sqrt") {
          g = arma::exp(de * -0.5 / temperature);
          lambda = arma::accu(g) - N; // n*exp(0) needs to be subtracted.
        } else if (mode == "barker") {
          g = 1.0 / (1.0 + arma::exp(de / temperature));
          lambda = arma::accu(g) - .5 * N;
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
      for (size_t s = 0; s < samples_per_walk; ++s) {
        for (size_t k = 0; k < burn_between; ++k) {
          double rand = uniform(rng) * lambda;
          double r_sum = 0.0;
          size_t i = 0;
          size_t q0 = 0;
          size_t q1 = 0;
          bool done = false;
          for (i = 0; i < N; i++) {
            for (q1 = 0; q1 < Q; q1++) {
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
          for (size_t pos = 0; pos < N; pos++) {
            for (size_t aa = 0; aa < Q; aa++) {
              if (pos < i) {
                de(pos, aa) += model->J(pos, i)(conf(pos), q1) -
                               model->J(pos, i)(conf(pos), q0) -
                               model->J(pos, i)(aa, q1) +
                               model->J(pos, i)(aa, q0);
              } else if (pos > i) {
                de(pos, aa) += model->J(i, pos)(q1, conf(pos)) -
                               model->J(i, pos)(q0, conf(pos)) -
                               model->J(i, pos)(q1, aa) +
                               model->J(i, pos)(q0, aa);
              } else {
                if (q1 == aa) {
                  de(pos, aa) = 0;
                } else if (q0 == aa) {
                  de(pos, aa) = -tmp;
                } else {
                  de(pos, aa) += model->h(q1, pos) - model->h(q0, pos);
                  for (size_t pos2 = 0; pos2 < N; pos2++) {
                    if (pos2 < i) {
                      de(pos, aa) += model->J(pos2, i)(conf(pos2), q1) -
                                     model->J(pos2, i)(conf(pos2), q0);
                    } else if (pos2 > i) {
                      de(pos, aa) += model->J(i, pos2)(q1, conf(pos2)) -
                                     model->J(i, pos2)(q0, conf(pos2));
                    }
                  }
                }
              }
            }
          }

          if (mode == "sqrt") {
            g = arma::exp(de * -0.5 / temperature);
            lambda = arma::accu(g) - N; // n*exp(0) needs to be subtracted.
          } else if (mode == "barker") {
            g = 1.0 / (1.0 + arma::exp(de / temperature));
            lambda = arma::accu(g) - .5 * N;
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

        for (size_t i = 0; i < N; ++i) {
          (*ptr)(s, i) = (int)conf(i);
        }
      }
    }
  }
  return;
};

void
Sampler::sampleSequencesZanella(arma::Mat<int>* ptr,
                                size_t walkers,
                                size_t burn_in,
                                long int seed,
                                std::string mode,
                                double temperature)
{
#pragma omp parallel
  {
#pragma omp for
    for (size_t walker = 0; walker < walkers; walker++){
      pcg32 rng;
      rng.seed(seed + walker);
      std::uniform_real_distribution<double> uniform(0, 1);

      arma::Mat<double> de = arma::Mat<double>(N, Q, arma::fill::zeros);
      arma::Mat<double> g = arma::Mat<double>(N, Q, arma::fill::zeros);
      double lambda = 0.0;

      arma::Mat<double> de_old = arma::Mat<double>(N, Q, arma::fill::zeros);
      arma::Mat<double> g_old = arma::Mat<double>(N, Q, arma::fill::zeros);
      double lambda_old = 0.0;

      arma::Col<size_t> conf = arma::Col<size_t>(N);
      for (size_t i = 0; i < N; ++i) {
        conf(i) = size_t(Q * uniform(rng));
        assert(conf(i) < Q);
      }

      // Compute initial neighborhood.
      for (size_t i = 0; i < N; i++) {
        size_t q0 = conf(i);
        double e0 = -model->h(q0, i);
        for (size_t j = 0; j < N; ++j) {
          if (i > j) {
            e0 -= model->J(j, i)(conf(j), q0);
          } else if (i < j) {
            e0 -= model->J(i, j)(q0, conf(j));
          }
        }

        for (size_t q1 = 0; q1 < Q; q1++) {
          if (q0 == q1) {
            de(i, q1) = 0.0;
          } else {
            double e1 = -model->h(q1, i);
            for (size_t j = 0; j < N; ++j) {
              if (i > j) {
                e1 -= model->J(j, i)(conf(j), q1);
              } else if (i < j) {
                e1 -= model->J(i, j)(q1, conf(j));
              }
            }
            de(i, q1) = e1 - e0;
          }
        }
      }

      if (mode == "sqrt") {
        g = arma::exp(de * -0.5 / temperature);
        lambda = arma::accu(g) - N; // n*exp(0) needs to be subtracted.
      } else if (mode == "barker") {
        g = 1.0 / (1.0 + arma::exp(de / temperature));
        lambda = arma::accu(g) - .5 * N;
      }

      de_old = de;
      g_old = g;
      lambda_old = lambda;

      for (size_t k = 0; k < burn_in; ++k) {
        double rand = uniform(rng) * lambda;
        double r_sum = 0.0;
        size_t i = 0;
        size_t q0 = 0;
        size_t q1 = 0;
        bool done = false;
        for (i = 0; i < N; i++) {
          for (q1 = 0; q1 < Q; q1++) {
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
        for (size_t pos = 0; pos < N; pos++) {
          for (size_t aa = 0; aa < Q; aa++) {
            if (pos < i) {
              de(pos, aa) += model->J(pos, i)(conf(pos), q1) -
                             model->J(pos, i)(conf(pos), q0) -
                             model->J(pos, i)(aa, q1) + model->J(pos, i)(aa, q0);
            } else if (pos > i) {
              de(pos, aa) += model->J(i, pos)(q1, conf(pos)) -
                             model->J(i, pos)(q0, conf(pos)) -
                             model->J(i, pos)(q1, aa) + model->J(i, pos)(q0, aa);
            } else {
              if (q1 == aa) {
                de(pos, aa) = 0;
              } else if (q0 == aa) {
                de(pos, aa) = -tmp;
              } else {
                de(pos, aa) += model->h(q1, pos) - model->h(q0, pos);
                for (size_t pos2 = 0; pos2 < N; pos2++) {
                  if (pos2 < i) {
                    de(pos, aa) += model->J(pos2, i)(conf(pos2), q1) -
                                   model->J(pos2, i)(conf(pos2), q0);
                  } else if (pos2 > i) {
                    de(pos, aa) += model->J(i, pos2)(q1, conf(pos2)) -
                                   model->J(i, pos2)(q0, conf(pos2));
                  }
                }
              }
            }
          }
        }

        if (mode == "sqrt") {
          g = arma::exp(de * -0.5 / temperature);
          lambda = arma::accu(g) - N; // n*exp(0) needs to be subtracted.
        } else if (mode == "barker") {
          g = 1.0 / (1.0 + arma::exp(de / temperature));
          lambda = arma::accu(g) - .5 * N;
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
      for (size_t i = 0; i < N; ++i) {
        (*ptr)(walker, i) = (int)conf(i);
      }
    }
  }
  return;
};
