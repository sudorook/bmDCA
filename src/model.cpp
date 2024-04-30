/*
 * SPDX-FileCopyrightText: 2020 sudorook <daemon@nullcodon.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "model.hpp"

#include <armadillo>

#include "utils.hpp"

/**
 * @brief Model constructor.
 */
Model::Model(void){};

/**
 * @brief Store the address of the MSA statstics.
 *
 * @param msa_train pointer to the statistics for the training MSA
 * @param msa_validate pointer to stats for the validation MSA
 */
void
Model::setMSAStats(MSAStats* msa_train, MSAStats* msa_validate)
{
  training = msa_train;
  validation = msa_validate;

  N = training->getN();
  Q = training->getQ();

  // Set the validate flag to true if a validation MSA is given (msa_validate
  // != nullptr).
  if (msa_validate) {
    validate = true;
  } else {
    validate = false;
  }
};

/**
 * @brief Store the address of the sampled sequence statistics.
 *
 * @param s pointer to the SampleStats sample statistics
 */
void
Model::setSampleStats(SampleStats* s)
{
  samples = s;
};

/**
 * @brief Set the iteration number.
 *
 * @param s iteration number
 *
 * The Model class needs to be updated with the step number, as it doesn't keep
 * track of it on its own.
 */
void
Model::setStep(int s)
{
  step = s;
};

/**
 * @brief Rescale parameters to have 0 marginals for J and h.
 *
 * Sets the 0-mean (Ising) gauge for the parameters. Note that it causes
 * regularization for J and h to blend together and become non-independent.
 */
void
Model::setZeroGauge(void)
{
  potts_model params_zg;

  // initialize
  params_zg.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params_zg.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params_zg.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  // rescale couplings
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i < j) {
        double J_ij_mean = arma::mean(arma::mean(params.J(i, j)));
        arma::Mat<double> J_ija_mean = arma::mean(params.J(i, j), 1);
        arma::Mat<double> J_ijb_mean = arma::mean(params.J(i, j), 0);
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            params_zg.J(i, j)(a, b) =
              params.J(i, j)(a, b) - J_ija_mean(a) - J_ijb_mean(b) + J_ij_mean;
          }
          params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
        }
      } else if (i > j) {
        double J_ij_mean = arma::mean(arma::mean(params.J(j, i)));
        arma::Mat<double> J_ija_mean = arma::mean(params.J(j, i), 0);
        arma::Mat<double> J_ijb_mean = arma::mean(params.J(j, i), 1);
        for (int a = 0; a < Q; a++) {
          params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
        }
      }
    }
  }

  // rescale fields
  arma::Row<double> h_i_mean = arma::mean(params.h, 0);
  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      params_zg.h(a, i) += params.h(a, i) - h_i_mean(i);
    }
  }

  params.h = params_zg.h;
  params.J = params_zg.J;
};
