/* Boltzmann-machine Direct Coupling Analysis (bmDCA)
 * Copyright (C) 2020
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <string>
#include <unistd.h>

#include "utils.hpp"

/**
 * @brief Sampler for Potts models.
 *
 * Sampler for Potts models. Each sample can be made independently or as part
 * of an MCMC trajectory. The sampler can use either Metropolis-Hastings or
 * Zanella update proposals. For reasons of computational efficiency in
 * computing sample statistics, separate sampling functions are implemented for
 * arma::Cube (for MCMC trajectories) and arma::Mat (independent samples).
 */
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
                      unsigned,
                      double = 1.0);
  void sampleSequences(arma::Cube<int>*,
                       size_t,
                       size_t,
                       size_t,
                       size_t,
                       unsigned,
                       double = 1.0);
  void sampleSequences(arma::Mat<int>*, size_t, size_t, unsigned, double = 1.0);
  void sampleSequencesZanellaSqrt(arma::Cube<int>*,
                                  size_t,
                                  size_t,
                                  size_t,
                                  size_t,
                                  unsigned,
                                  double = 1.0);
  void sampleSequencesZanellaSqrt(arma::Mat<int>*,
                                  size_t,
                                  size_t,
                                  unsigned,
                                  double = 1.0);
  void sampleSequencesZanellaBarker(arma::Cube<int>*,
                                    size_t,
                                    size_t,
                                    size_t,
                                    size_t,
                                    unsigned,
                                    double = 1.0);
  void sampleSequencesZanellaBarker(arma::Mat<int>*,
                                    size_t,
                                    size_t,
                                    unsigned,
                                    double = 1.0);

private:
  const size_t N; ///< number of positions
  const size_t Q; ///< number of amino acids (inc. gaps)

  potts_model* model; ///< pointer to Potts parameters
};

#endif
