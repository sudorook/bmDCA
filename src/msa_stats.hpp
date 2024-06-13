/*
 * SPDX-FileCopyrightText: 2020 - 2022 sudorook <daemon@nullcodon.com>
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

#ifndef SRC_MSA_STATS_HPP_
#define SRC_MSA_STATS_HPP_

#include "msa.hpp"

#include <armadillo>
#include <memory>
#include <string>

/**
 * @brief Class for computing the sequence statistics for a MSA.
 *
 * Compute the 1p and 2p statistics, nothing higher.
 */
class MSAStats
{
public:
  explicit MSAStats(std::shared_ptr<MSA>, bool = false);
  void updateMSA(std::shared_ptr<MSA>, bool = false);

  double getEffectiveM();
  double getN();
  double getM();
  double getQ();

  void computeErrorMSA(int = 100, unsigned = 0);

  // void writeRelEntropy(std::string);
  // void writeRelEntropyAscii(std::string);
  // void writeRelEntropyPos(std::string);
  // void writeRelEntropyPosAscii(std::string);
  // void writeRelEntropyGradient(std::string);
  // void writeRelEntropyGradientAscii(std::string);
  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);

  arma::Mat<double> frequency_1p;              ///< single-position frequencies
  arma::field<arma::Mat<double>> frequency_2p; ///< 2-position frequencies
  arma::Mat<double> rel_entropy_1p;            ///< relative entropy (1p)
  arma::Col<double> rel_entropy_pos_1p;        ///< total 1p relative entropy
  arma::Mat<double> rel_entropy_grad_1p;       ///< 1p relative entropy gradient

  double freq_rms;           ///< multinomial standard error estimate of MSA
  arma::Col<double> msa_rms; ///< RMSE of two equal subset of the MSA

private:
  const double pseudocount = 0.03; ///< pseudocount for missing data

  int M;              ///< number of sequences
  int N;              ///< number of positions
  int Q;              ///< amino acid alphabet size
  double M_effective; ///< effect number of sequences

  std::shared_ptr<MSA> msa; ///< address of MSA for which to compute stats

  // void computeMSAStats(MSA*);
  void computeMSAStats();

  arma::Col<double> aa_background_frequencies; ///< background aa frequencies
};

#endif // SRC_MSA_STATS_HPP_
