#include "mcmc.hpp"

#include <string>

#include "graph.hpp"

void
MCMC::load(potts_model* model)
{
  graph.load(model);
};

MCMC::MCMC(size_t N, size_t Q)
  : graph(N, Q)
{
  n = N;
  q = Q;
};

MCMC::MCMC(size_t N, size_t Q, potts_model* params)
  : graph(N, Q)
{
  n = N;
  q = Q;
  graph.load(params);
};

void
MCMC::sample_energies(arma::Mat<double>* ptr,
                      int reps,
                      int M,
                      int N,
                      int t_wait,
                      int delta_t,
                      long int seed,
                      double temperature){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps;
         rep++){ graph.sample_mcmc_energies(((*ptr).colptr(rep)),
                                            M,
                                            t_wait,
                                            delta_t,
                                            seed + rep,
                                            temperature);
    }
  }
};

void
MCMC::sample(arma::Cube<int>* ptr,
             int reps,
             int M,
             int N,
             int t_wait,
             int delta_t,
             long int seed,
             double temperature){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps;
         rep++){ graph.sample_mcmc((arma::Mat<int>*)&((*ptr).slice(rep)),
                                   M,
                                   t_wait,
                                   delta_t,
                                   seed + rep,
                                   temperature);
    }
  }
};

void
MCMC::sample_zanella(arma::Cube<int>* ptr,
                     int reps,
                     int M,
                     int N,
                     int t_wait,
                     int delta_t,
                     long int seed,
                     std::string mode,
                     double temperature){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps; rep++){
      graph.sample_mcmc_zanella((arma::Mat<int>*)&((*ptr).slice(rep)),
                                M,
                                t_wait,
                                delta_t,
                                seed + rep,
                                mode,
                                temperature);
    }
  }
};

void
MCMC::sample_init(arma::Cube<int>* ptr,
                  int reps,
                  int M,
                  int N,
                  int t_wait,
                  int delta_t,
                  arma::Col<int>* init_ptr,
                  long int seed,
                  double temperature){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps;
         rep++){ graph.sample_mcmc_init((arma::Mat<int>*)&((*ptr).slice(rep)),
                                        M,
                                        t_wait,
                                        delta_t,
                                        init_ptr,
                                        seed + rep,
                                        temperature);
    }
  }
};
