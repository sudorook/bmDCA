#ifndef MCMC_HPP
#define MCMC_HPP

#include <string>
#include <unistd.h>

#include "graph.hpp"
#include "utils.hpp"

class MCMC
{

public:
  MCMC(size_t N, size_t Q);
  MCMC(size_t N, size_t Q, potts_model*);
  void load(potts_model*);
  void run(int, int);
  void
  sample(arma::Cube<int>*, int, int, int, int, int, long int, double = 1.0);
  void sample_energies(arma::Mat<double>*,
                       int,
                       int,
                       int,
                       int,
                       int,
                       long int,
                       double = 1.0);
  void sample_zanella(arma::Cube<int>*,
                      int,
                      int,
                      int,
                      int,
                      int,
                      long int,
                      std::string,
                      double = 1.0);
  void sample_init(arma::Cube<int>*,
                   int,
                   int,
                   int,
                   int,
                   int,
                   arma::Col<int>*,
                   long int,
                   double = 1.0);

private:
  size_t n; // number of positions
  size_t q; // number of amino acids (inc. gaps)

  Graph graph;
};

#endif
