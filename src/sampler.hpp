#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <string>
#include <unistd.h>

#include "utils.hpp"

class Sampler
{
public:
  Sampler(size_t, size_t);
  Sampler(size_t, size_t, potts_model*);
  void setModel(potts_model*);

  void sampleEnergies(arma::Mat<double>*,
                      size_t,
                      size_t,
                      size_t,
                      size_t,
                      long int,
                      double = 1.0);
  void sampleSequences(arma::Cube<int>*,
                       size_t,
                       size_t,
                       size_t,
                       size_t,
                       long int,
                       double = 1.0);
  void sampleSequences(arma::Mat<int>*,
                       size_t,
                       size_t,
                       long int,
                       double = 1.0);
  void sampleSequencesZanella(arma::Cube<int>*,
                              size_t,
                              size_t,
                              size_t,
                              size_t,
                              long int,
                              std::string = "sqrt",
                              double = 1.0);
  void sampleSequencesZanella(arma::Mat<int>*,
                              size_t,
                              size_t,
                              long int,
                              std::string = "sqrt",
                              double = 1.0);

private:
  const size_t N; // number of positions
  const size_t Q; // number of amino acids (inc. gaps)

  potts_model *model;
};

#endif
