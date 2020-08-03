#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <armadillo>
#include <iostream>
#include <string>

#include "utils.hpp"

class Graph
{
public:
  Graph(size_t n, size_t q, potts_model* p)
    : n(n)
    , q(q)
    , params(p){};

  Graph(size_t n, size_t q)
    : n(n)
    , q(q){};

  void load(potts_model*);

  size_t n, q;

  potts_model* params;

  std::ostream& print_distribution(std::ostream& os);

  std::ostream& print_parameters(std::ostream& os);

  void sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed,
                   double = 1.0);

  void sample_mcmc_energies(double* ptr,
                            size_t m,
                            size_t mc_iters0,
                            size_t mc_iters,
                            long int seed,
                            double = 1.0);

  void sample_mcmc_init(arma::Mat<int>* ptr,
                        size_t m,
                        size_t mc_iters0,
                        size_t mc_iters,
                        arma::Col<int>* init_ptr,
                        long int seed,
                        double = 1.0);

  void sample_mcmc_zanella(arma::Mat<int>* ptr,
                           size_t,
                           size_t,
                           size_t,
                           long int,
                           std::string = "sqrt",
                           double = 1.0);

  void print_parameters(FILE* of);
};

#endif
