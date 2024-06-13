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

#include "msa_stats.hpp"

// #include <algorithm>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
// #include <numeric>
#include <unordered_set>
#include <vector>

#include "pcg_random.hpp"

/**
 * @brief MSAStats constructor.
 *
 * @param msa address of MSA instance
 * @param verbose flag to print MSA info
 */
MSAStats::MSAStats(MSA* msa, bool verbose)
  : msa(msa)
{
  // Initialize
  N = msa->N;
  M = msa->M;
  Q = msa->Q;
  M_effective = sum(msa->sequence_weights);

  if (verbose) {
    std::cout << M << " sequences" << std::endl;
    std::cout << N << " positions" << std::endl;
    std::cout << Q << " amino acids (including gaps)" << std::endl;
    std::cout << M_effective << " effective sequences" << std::endl;
  }

  frequency_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
      }
    }
  }
  rel_entropy_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  rel_entropy_pos_1p = arma::Col<double>(N, arma::fill::zeros);
  rel_entropy_grad_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  aa_background_frequencies = arma::Col<double>(Q, arma::fill::ones);

  if (Q == 21) {
    // fixed background frequencies
    aa_background_frequencies = { 0.000, 0.073, 0.025, 0.050, 0.061, 0.042,
                                  0.072, 0.023, 0.053, 0.064, 0.089, 0.023,
                                  0.043, 0.052, 0.040, 0.052, 0.073, 0.056,
                                  0.063, 0.013, 0.033 };
  } else {
    aa_background_frequencies =
      aa_background_frequencies / static_cast<double>(Q);
  }

  computeMSAStats();
};

/**
 * @brief Update MSAStats with new MSA
 *
 * @param new_msa address of new MSA
 * @param verbose print MSA info
 */
void
MSAStats::updateMSA(MSA* new_msa, bool verbose)
{
  msa = new_msa;

  // Initialize
  N = msa->N;
  M = msa->M;
  Q = msa->Q;
  M_effective = sum(msa->sequence_weights);

  if (verbose) {
    std::cout << M << " sequences" << std::endl;
    std::cout << N << " positions" << std::endl;
    std::cout << Q << " amino acids (including gaps)" << std::endl;
    std::cout << M_effective << " effective sequences" << std::endl;
  }

  computeMSAStats();
};

/**
 * @brief Compute the 1p and 2p statistics for an MSA.
 */
void
MSAStats::computeMSAStats()
{
  frequency_1p.zeros();
  rel_entropy_1p.zeros();
  rel_entropy_pos_1p.zeros();
  rel_entropy_grad_1p.zeros();
  aa_background_frequencies.zeros();

  // Compute the frequecies (1p statistics) for amino acids (and gaps) for each
  // position. Use pointers to make things speedier.
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      int* align_ptr = msa->alignment.colptr(i);
      double* freq_ptr = frequency_1p.colptr(i);
      double* weight_ptr = msa->sequence_weights.memptr();
      for (int m = 0; m < M; m++) {
        *(freq_ptr + *(align_ptr + m)) += *(weight_ptr + m);
      }
    }
  }
  frequency_1p = frequency_1p / M_effective;

  double mean_1p_var =
    arma::accu(frequency_1p % (1. - frequency_1p)) / static_cast<double>(N * Q);

  // Compute the 2p statistics
  arma::Col<double> mean_2p_var_vec = arma::Col<double>(N, arma::fill::zeros);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double* weight_ptr = msa->sequence_weights.memptr();
        frequency_2p(i, j).zeros();

        int* align_ptr1 = msa->alignment.colptr(i);
        int* align_ptr2 = msa->alignment.colptr(j);
        for (int m = 0; m < M; m++) {
          frequency_2p(i, j)(*(align_ptr1 + m), *(align_ptr2 + m)) +=
            *(weight_ptr + m);
        }
        frequency_2p(i, j) = frequency_2p(i, j) / M_effective;
        mean_2p_var_vec(i) +=
          arma::accu(frequency_2p(i, j) % (1. - frequency_2p(i, j)));
      }
    }
  }

  double mean_2p_var =
    arma::accu(mean_2p_var_vec) / (N * (N - 1.) / (2. * Q * Q));
  freq_rms = sqrt(mean_1p_var) + sqrt(mean_2p_var);

  // Update the background frequencies based by computing overall gap frequency
  // theta.
  double theta = 0;
  for (int i = 0; i < N; i++) {
    theta += frequency_1p(0, i);
  }
  theta = theta / N;
  aa_background_frequencies[0] = theta;
  for (int i = 1; i < Q; i++) {
    aa_background_frequencies[i] = aa_background_frequencies[i] * (1. - theta);
  }

  // Use the positional and backgrounds frequencies to estimate the relative
  // entropy gradient for each position.
  arma::Mat<double> tmp = frequency_1p * (1. - pseudocount);
  tmp.each_col() += pseudocount * aa_background_frequencies;
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      double pos_freq = tmp(aa, i);
      double background_freq = aa_background_frequencies(aa);
      if (pos_freq < 1. && pos_freq > 0.) {
        rel_entropy_pos_1p(i) += pos_freq * log(pos_freq / background_freq);
        rel_entropy_1p(aa, i) =
          pos_freq * log(pos_freq / background_freq) +
          (1 - pos_freq) * log((1 - pos_freq) / (1 - background_freq));
        rel_entropy_grad_1p(aa, i) = log((pos_freq * (1. - background_freq)) /
                                         ((1. - pos_freq) * background_freq));
      }
    }
  }
};

/**
 * @brief Compute RMSE for multiple subsets of the MSA.
 *
 * @param reps number of RMSE replicates
 * @param seed random seed for RNG
 */
void
MSAStats::computeErrorMSA(int reps, unsigned seed)
{
  msa_rms = arma::Col<double>(reps, arma::fill::zeros);

  // partition MSA into equal subsets
  int M_1 = (M + 1) / 2;
  int M_2 = M / 2;

  arma::Mat<int> alignment_T = (msa->alignment).t();

  pcg32 rng;
  rng.seed(seed);

  for (int rep = 0; rep < reps; rep++) {
    std::unordered_set<int> elems;

    for (int r = M - M_1; r < M; ++r) {
      // int v = static_cast<int>(rng());
      size_t v = rng();
      if (!elems.insert(v).second) {
        elems.insert(r);
      }
    }

    std::vector<int> idx_1(elems.begin(), elems.end());
    std::vector<int> idx_2(M_2);
    int counter = 0;
    for (int i = 0; i < M; i++) {
      if (elems.find(i) == elems.end()) {
        idx_2.at(counter) = i;
        counter++;
      }
    }

    arma::Col<int> idx_1_v2 = arma::conv_to<arma::Col<int>>::from(idx_1);
    arma::Col<int> idx_2_v2 = arma::conv_to<arma::Col<int>>::from(idx_2);

    arma::Mat<int> alignment_1 = arma::Mat<int>(N, M_1, arma::fill::zeros);
    arma::Mat<int> alignment_2 = arma::Mat<int>(N, M_2, arma::fill::zeros);
    arma::Col<double> weights_1 = arma::Col<double>(M_1, arma::fill::zeros);
    arma::Col<double> weights_2 = arma::Col<double>(M_2, arma::fill::zeros);

    // Partition the MSA into two subsets. Use the full MSA sequence weights to
    // estimate effective sequence size for each subset.
#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < M_1; i++) {
        alignment_1.col(i) = alignment_T.col(idx_1.at(i));
        weights_1(i) = msa->sequence_weights(idx_1.at(i));
      }
    }

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < M_2; i++) {
        alignment_2.col(i) = alignment_T.col(idx_2.at(i));
        weights_2(i) = msa->sequence_weights(idx_2.at(i));
      }
    }

    MSA msa_1 = MSA(alignment_1.t(), weights_1, M_1, N, Q);
    MSA msa_2 = MSA(alignment_2.t(), weights_2, M_2, N, Q);

    double M_1_effective = arma::sum(msa_1.sequence_weights);
    double M_2_effective = arma::sum(msa_2.sequence_weights);

    arma::Mat<double> msa_1_frequency_1p =
      arma::Mat<double>(Q, N, arma::fill::zeros);
    arma::field<arma::Mat<double>> msa_1_frequency_2p =
      arma::field<arma::Mat<double>>(N, N);

    arma::Mat<double> msa_2_frequency_1p =
      arma::Mat<double>(Q, N, arma::fill::zeros);
    arma::field<arma::Mat<double>> msa_2_frequency_2p =
      arma::field<arma::Mat<double>>(N, N);

    // Compute the frequecies (1p statistics) for amino acids (and gaps) for
    // each position. Use pointers to make things speedier.
#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < N; i++) {
        int* align_ptr = msa_1.alignment.colptr(i);
        double* freq_ptr = msa_1_frequency_1p.colptr(i);
        double* weight_ptr = msa_1.sequence_weights.memptr();
        for (int m = 0; m < M_1; m++) {
          *(freq_ptr + *(align_ptr + m)) += *(weight_ptr + m);
        }
      }
    }
    msa_1_frequency_1p = msa_1_frequency_1p / M_1_effective;

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < N; i++) {
        int* align_ptr = msa_2.alignment.colptr(i);
        double* freq_ptr = msa_2_frequency_1p.colptr(i);
        double* weight_ptr = msa_2.sequence_weights.memptr();
        for (int m = 0; m < M_2; m++) {
          *(freq_ptr + *(align_ptr + m)) += *(weight_ptr + m);
        }
      }
    }
    msa_2_frequency_1p = msa_2_frequency_1p / M_2_effective;

    // Compute the 2p statistics
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          double* weight_ptr = msa_1.sequence_weights.memptr();
          msa_1_frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);

          int* align_ptr1 = msa_1.alignment.colptr(i);
          int* align_ptr2 = msa_1.alignment.colptr(j);
          for (int m = 0; m < M_1; m++) {
            msa_1_frequency_2p(i, j)(*(align_ptr1 + m), *(align_ptr2 + m)) +=
              *(weight_ptr + m);
          }
          msa_1_frequency_2p(i, j) = msa_1_frequency_2p(i, j) / M_1_effective;
        }
      }
    }
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          double* weight_ptr = msa_2.sequence_weights.memptr();
          msa_2_frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);

          int* align_ptr1 = msa_2.alignment.colptr(i);
          int* align_ptr2 = msa_2.alignment.colptr(j);
          for (int m = 0; m < M_2; m++) {
            msa_2_frequency_2p(i, j)(*(align_ptr1 + m), *(align_ptr2 + m)) +=
              *(weight_ptr + m);
          }
          msa_2_frequency_2p(i, j) = msa_2_frequency_2p(i, j) / M_2_effective;
        }
      }
    }

    // Compute the difference in 1p stats between the two subsets.
    double error_1p =
      arma::accu(arma::pow(msa_1_frequency_1p - msa_2_frequency_1p, 2));
    error_1p = sqrt(error_1p / (N * Q));

    // Compute the difference in 2p stats between the two subsets.
    arma::Col<double> error_2p_vec = arma::Col<double>(N, arma::fill::zeros);
#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          error_2p_vec(i) = arma::accu(
            arma::pow(msa_1_frequency_2p(i, j) - msa_2_frequency_2p(i, j), 2));
        }
      }
    }
    double error_2p = sqrt(arma::accu(error_2p_vec) / (N * (N - 1) * Q * Q));

    // Compute the RMSE of the 1p and 2p stats.
    msa_rms(rep) = (error_1p + error_2p) /
                   sqrt(M_effective / (M_1_effective + M_2_effective) * 2);
  }

  return;
}

/**
 * @brief Getter function for number of states.
 *
 * @return number of states
 */
double
MSAStats::getQ(void)
{
  return Q;
};

/**
 * @brief Getter function for number of sequences.
 *
 * @return number of sequences
 */
double
MSAStats::getM(void)
{
  return M;
};

/**
 * @brief Getter function for number of positions.
 *
 * @return number of positions
 */
double
MSAStats::getN(void)
{
  return N;
};

/**
 * @brief Getter function for effective number of sequences.
 *
 * @return number of effective sequences
 */
double
MSAStats::getEffectiveM(void)
{
  return M_effective;
};

// void
// MSAStats::writeRelEntropyGradient(std::string output_file)
// {
//   rel_entropy_grad_1p.save(output_file, arma::arma_binary);
// };

// void
// MSAStats::writeRelEntropyGradientAscii(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//   for (int i = 0; i < N; i++) {
//     output_stream << i;
//     for (int aa = 0; aa < Q; aa++) {
//       output_stream << " " << rel_entropy_grad_1p(aa, i);
//     }
//     output_stream << std::endl;
//   }
// };

// void
// MSAStats::writeRelEntropyPos(std::string output_file)
// {
//   rel_entropy_pos_1p.save(output_file, arma::arma_binary);
// };

// void
// MSAStats::writeRelEntropyPosAscii(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//   for (int i = 0; i < N; i++) {
//     output_stream << " " << rel_entropy_pos_1p(i) << std::endl;
//   }
// };

// void
// MSAStats::writeRelEntropy(std::string output_file)
// {
//   rel_entropy_1p.save(output_file, arma::arma_binary);
// };

// void
// MSAStats::writeRelEntropyAscii(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//   for (int i = 0; i < N; i++) {
//     output_stream << i;
//     for (int aa = 0; aa < Q; aa++) {
//       output_stream << " " << rel_entropy_1p(aa, i);
//     }
//     output_stream << std::endl;
//   }
// };

/**
 * @brief Write 1p stats (arma::binary)
 *
 * @param output_file file string for output
 */
void
MSAStats::writeFrequency1p(std::string output_file)
{
  frequency_1p.save(output_file, arma::arma_binary);
};

/**
 * @brief Write 1p stats in text format
 *
 * @param output_file file string for output
 */
void
MSAStats::writeFrequency1pAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int i = 0; i < N; i++) {
    output_stream << i;
    for (int aa = 0; aa < Q; aa++) {
      output_stream << " " << frequency_1p(aa, i);
    }
    output_stream << std::endl;
  }
};

/**
 * @brief Write 2p stats (arma::binary)
 *
 * @param output_file file string for output
 */
void
MSAStats::writeFrequency2p(std::string output_file)
{
  frequency_2p.save(output_file, arma::arma_binary);
};

/**
 * @brief Write 2p stats in text format
 *
 * @param output_file file string for output
 */
void
MSAStats::writeFrequency2pAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      output_stream << i << " " << j;
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << " " << frequency_2p(i, j)(aa1, aa2);
        }
      }
      output_stream << std::endl;
    }
  }
};
