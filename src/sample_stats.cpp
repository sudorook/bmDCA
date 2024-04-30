/*
 * SPDX-FileCopyrightText: 2020 - 2021 sudorook <daemon@nullcodon.com>
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

#include "sample_stats.hpp"

#include <armadillo>
#include <iostream>

#include "utils.hpp"

/**
 * @brief Constructor for SampleStats2D.
 *
 * @param s pointer to sequences (arma::Mat)
 * @param p pointer to Potts model
 * @param p_prev pointer to previous Potts model
 */
SampleStats2D::SampleStats2D(arma::Mat<int>* s,
                             potts_model* p,
                             potts_model* p_prev)
{
  M = s->n_rows;
  N = s->n_cols;
  Q = p->h.n_rows;

  samples = s;
  params = p;
  params_prev = p_prev;
};

/**
 * @brief Wrapper for computing 1p/2p statistics.
 */
void
SampleStats2D::computeStats(void)
{
  computeSampleStats();
};

/**
 * @brief Wrapper for computing sequence energies and correlations.
 */
void
SampleStats2D::computeStatsExtra(void)
{
  computeEnergies();
  computeCorrelations();
};

/**
 * @brief Wrapper for importance sampling loop.
 */
void
SampleStats2D::computeStatsImportance(void)
{
  computeSampleStatsImportance();
};

/**
 * @brief Set burn-between time used for sampling.
 *
 * @param m burn-between time
 */
void
SampleStats2D::setMixingTime(int m)
{
  mixing_time = m;
}

/**
 * @brief Write sample statistics for a step.
 *
 * @param step iteration number
 * @param output_binary flag to write arma::binary or text files
 */
void
SampleStats2D::writeStep(int step, bool output_binary)
{
  if (output_binary) {
    writeFrequency1p("samples_stat_1p_" + std::to_string(step) + ".bin");
    writeFrequency2p("samples_stat_2p_" + std::to_string(step) + ".bin");
  } else {
    writeFrequency1pAscii("samples_stat_1p_" + std::to_string(step) + ".txt");
    writeFrequency2pAscii("samples_stat_2p_" + std::to_string(step) + ".txt");
  }

  writeSamples("samples_" + std::to_string(step) + ".txt");
  writeSampleEnergies("energies_" + std::to_string(step) + ".txt");
};

/**
 * @brief Write sample statistics.
 *
 * @param str label for output data
 * @param output_binary flag to write arma::binary or text files
 */
void
SampleStats2D::writeData(std::string str, bool output_binary)
{
  if (output_binary) {
    writeFrequency1p("samples_stat_1p_" + str + ".bin");
    writeFrequency2p("samples_stat_2p_" + str + ".bin");
  } else {
    writeFrequency1pAscii("samples_stat_1p_" + str + ".txt");
    writeFrequency2pAscii("samples_stat_2p_" + str + ".txt");
  }

  writeSamples("samples_" + str + ".txt");
  writeSampleEnergies("energies_" + str + ".txt");
};

/**
 * @brief Return key statistics for the sampled sequences.
 *
 * @return arma::Col containing sequence correlations, correlation variance,
 * and importance sampling stats.
 */
arma::Col<double>
SampleStats2D::getStats(void)
{
  arma::Col<double> stats = arma::Col<double>(5, arma::fill::zeros);
  stats(0) = overlap_inf;
  stats(1) = overlap_inf_sigma;
  stats(2) = Z_ratio;
  stats(3) = sumw_inv;
  stats(4) = dE_av_tot;
  return stats;
};

/**
 * @brief Compute sequence energies.
 */
void
SampleStats2D::computeEnergies(void)
{
  energies = arma::Col<double>(M, arma::fill::zeros);
#pragma omp parallel
  {
#pragma omp for
    for (int seq = 0; seq < M; seq++) {
      double E = 0;
      for (int i = 0; i < N; i++) {
        E -= params->h(samples->at(seq, i), i);
        for (int j = i + 1; j < N; j++) {
          E -= params->J(i, j)(samples->at(seq, i), samples->at(seq, j));
        }
      }
      energies(seq) = E;
    }
  }
};

/**
 * @brief Compute 1p and 2p statistics.
 */
void
SampleStats2D::computeSampleStats(void)
{
  frequency_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      int* align_ptr = samples->colptr(i);
      double* freq_ptr = frequency_1p.colptr(i);
      for (int m = 0; m < M; m++) {
        (*(freq_ptr + *(align_ptr + m)))++;
      }
    }
  }
  frequency_1p = frequency_1p / M;

  // Compute the 2p statistics
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        int* align_ptr1 = samples->colptr(i);
        int* align_ptr2 = samples->colptr(j);
        for (int m = 0; m < M; m++) {
          frequency_2p(i, j)(*(align_ptr1 + m), *(align_ptr2 + m))++;
        }
        frequency_2p(i, j) = frequency_2p(i, j) / M;
      }
    }
  }
};

/**
 * @brief Importance sampling.
 */
void
SampleStats2D::computeSampleStatsImportance(void)
{
  arma::Col<double> p = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<double> dE = arma::Col<double>(M, arma::fill::zeros);
  double Z_tot = 0;
  double Z_inv_tot = 0;
  // double W = 0;
  double w = 0;
  double sumw = 0;

  double Z = 0;
  double Z_inv = 0;
  double sum = 0;
  double dE_av = 0;

  Z_ratio = 0;
  sumw_inv = 0;
  dE_av_tot = 0;

  for (int m = 0; m < M; m++) {
    for (int i = 0; i < N; i++) {
      dE(m) += params->h(samples->at(m, i)) - params_prev->h(samples->at(m, i));
      for (int j = i + 1; j < N; j++) {
        dE(m) += params->J(i, j)(samples->at(m, i), samples->at(m, j)) -
                 params_prev->J(i, j)(samples->at(m, i), samples->at(m, j));
      }
    }
  }
  dE_av = arma::mean(dE);
  dE_av_tot = dE_av;
  for (int m = 0; m < M; m++) {
    p(m) = exp(dE(m) - dE_av);
    Z += p(m);
    Z_inv += 1. / p(m);
  }
  for (int m = 0; m < M; m++) {
    p(m) = p(m) / Z;
    sum += pow(p(m), 2);
  }
  Z_tot = Z;
  Z_inv_tot += Z_inv;
  // W = 1. / sum;
  w = 1. / sum;

  arma::Col<int> n1 = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2 = arma::Mat<int>(Q, Q, arma::fill::zeros);

  Z_ratio = Z_tot / Z_inv_tot;
  sumw_inv = 1.0 / sumw;

  arma::Col<int> n1av = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2av = arma::Mat<int>(Q, Q, arma::fill::zeros);

  arma::Col<int> n1squared = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2squared = arma::Mat<int>(Q, Q, arma::fill::zeros);

  for (int i = 0; i < N; i++) {
    n1.zeros();
    n1av.zeros();
    n1squared.zeros();
    for (int m = 0; m < M; m++) {
      n1(samples->at(m, i)) += p(m);
    }
    for (int aa = 0; aa < Q; aa++) {
      n1av(aa) += w * n1(aa);
      n1squared(aa) += w * pow(n1(aa), 2);
      frequency_1p(aa, i) = (double)n1av(aa);
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      n2.zeros();
      n2av.zeros();
      n2squared.zeros();
      for (int m = 0; m < M; m++) {
        n2(samples->at(m, i), samples->at(m, j)) += p(m);
      }

      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          n2av(aa1, aa2) += w * n2(aa1, aa2);
          n2squared(aa1, aa2) += w * pow(n2(aa1, aa2), 2);
          frequency_2p(i, j)(aa1, aa2) = (double)n2av(aa1, aa2);
        }
      }
    }
  }
};

/**
 * @brief Compute sequence correlations.
 */
void
SampleStats2D::computeCorrelations(void)
{
  arma::Col<double> d = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<double> d2 = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<int> count = arma::Col<int>(M, arma::fill::zeros);

  arma::Mat<int> samples_T = samples->t();
  // Compute distances between sequences
  for (int seq1 = 0; seq1 < M; seq1++) {
    int* seq1_ptr = samples_T.colptr(seq1);
    for (int seq2 = seq1 + 1; seq2 < M; seq2++) {
      int* seq2_ptr = samples_T.colptr(seq2);
      int id = 0;
      for (int i = 0; i < N; i++) {
        if (*(seq1_ptr + i) == *(seq2_ptr + i)) {
          id++;
        }
      }
      d(seq2 - seq1) += (double)id / N;
      d2(seq2 - seq1) += (double)id * id / (N * N);
      count(seq2 - seq1)++;
    }
  }

  double dinf = arma::sum(d);
  double dinf2 = arma::sum(d2);
  overlap_inf = 2.0 * dinf / (double)(M * (M - 1));
  overlap_inf_sigma = sqrt(2.0 / (M * (M - 1))) *
                      sqrt(2.0 * dinf2 / (double)((M - 1) * M) -
                           pow(2.0 * dinf / (double)(M * (M - 1)), 2));
};

/**
 * @brief Write 1p frequencies in binary format.
 *
 * @param output_file file string for output
 */
void
SampleStats2D::writeFrequency1p(std::string output_file)
{
  frequency_1p.save(output_file, arma::arma_binary);
};

/**
 * @brief Write 1p frequencies in text format.
 *
 * @param output_file file string for output
 */
void
SampleStats2D::writeFrequency1pAscii(std::string output_file)
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
 * @brief Write 2p frequencies in binary format.
 *
 * @param output_file file string for output
 */
void
SampleStats2D::writeFrequency2p(std::string output_file)
{
  frequency_2p.save(output_file, arma::arma_binary);
};

/**
 * @brief Write 2p frequencies in text format.
 *
 * @param output_file file string for output
 */
void
SampleStats2D::writeFrequency2pAscii(std::string output_file)
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

/**
 * @brief Write samples to disk.
 *
 * @param output_file file string for output
 */
void
SampleStats2D::writeSamples(std::string output_file)
{
  std::ofstream output_stream(output_file);

  output_stream << M << " " << N << " " << Q << std::endl;
  for (int m = 0; m < M; m++) {
    output_stream << samples->at(m, 0);
    for (int i = 1; i < N; i++) {
      output_stream << " " << samples->at(m, i);
    }
    output_stream << std::endl;
  }
};

/**
 * @brief Write sampled sequence energies to disk.
 *
 * @param output_file file string for output.
 */
void
SampleStats2D::writeSampleEnergies(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int m = 0; m < M; m++) {
    output_stream << energies(m) << std::endl;
  }
};

/**
 * @brief Constructor for SampleStats3D.
 *
 * @param s pointer to sequences (arma::Cube)
 * @param p pointer to Potts model
 * @param p_prev pointer to previous Potts model
 */
SampleStats3D::SampleStats3D(arma::Cube<int>* s,
                             potts_model* p,
                             potts_model* p_prev)
{
  M = s->n_rows;
  N = s->n_cols;
  reps = s->n_slices;
  Q = p->h.n_rows;

  samples = s;
  params = p;
  params_prev = p_prev;
};

/**
 * @brief Wrapper for computing 1p/2p statistics.
 */
void
SampleStats3D::computeStats(void)
{
  computeSampleStats();
};

/**
 * @brief Wrapper for computing sequence energies and correlations.
 */
void
SampleStats3D::computeStatsExtra(void)
{
  computeEnergies();
  computeEnergiesStats();
  computeCorrelations();
};

/**
 * @brief Wrapper for importance sampling loop.
 */
void
SampleStats3D::computeStatsImportance(void)
{
  computeSampleStatsImportance();
};

/**
 * @brief Set burn-between time used for sampling.
 *
 * @param m burn-between time
 */
void
SampleStats3D::setMixingTime(int m)
{
  mixing_time = m;
}

/**
 * @brief Write sample statistics for a step.
 *
 * @param step iteration number
 * @param output_binary flag to write arma::binary or text files
 */
void
SampleStats3D::writeStep(int step, bool output_binary)
{
  if (output_binary) {
    writeFrequency1p("samples_stat_1p_" + std::to_string(step) + ".bin",
                     "samples_stat_1p_sigma_" + std::to_string(step) + ".bin");
    writeFrequency2p("samples_stat_2p_" + std::to_string(step) + ".bin",
                     "samples_stat_2p_sigma_" + std::to_string(step) + ".bin");
  } else {
    writeFrequency1pAscii("samples_stat_1p_" + std::to_string(step) + ".txt",
                          "samples_stat_1p_sigma_" + std::to_string(step) +
                            ".txt");
    writeFrequency2pAscii("samples_stat_2p_" + std::to_string(step) + ".txt",
                          "samples_stat_2p_sigma_" + std::to_string(step) +
                            ".txt");
  }

  writeSamples("samples_" + std::to_string(step) + ".txt");
  writeSampleEnergies("energies_" + std::to_string(step) + ".txt");
  writeSampleEnergiesRelaxation("energies_relax_" + std::to_string(step) +
                                ".txt");
};

/**
 * @brief Write sample statistics.
 *
 * @param str label for output data
 * @param output_binary flag to write arma::binary or text files
 */
void
SampleStats3D::writeData(std::string str, bool output_binary)
{
  if (output_binary) {
    writeFrequency1p("samples_stat_1p_" + str + ".bin",
                     "samples_stat_1p_sigma_" + str + ".bin");
    writeFrequency2p("samples_stat_2p_" + str + ".bin",
                     "samples_stat_2p_sigma_" + str + ".bin");
  } else {
    writeFrequency1pAscii("samples_stat_1p_" + str + ".txt",
                          "samples_stat_1p_sigma_" + str + ".txt");
    writeFrequency2pAscii("samples_stat_2p_" + str + ".txt",
                          "samples_stat_2p_sigma_" + str + ".txt");
  }

  writeSamples("samples_" + str + ".txt");
  writeSampleEnergies("energies_" + str + ".txt");
  writeSampleEnergiesRelaxation("energies_relax_" + str + ".txt");
};

/**
 * @brief Return key statistics for the sampled sequences.
 *
 * @return arma::Col containing sequence correlations, correlation variance,
 * and importance sampling stats.
 */
arma::Col<double>
SampleStats3D::getStats()
{
  arma::Col<double> stats = arma::Col<double>(19, arma::fill::zeros);
  stats(0) = energies_start_avg;
  stats(1) = energies_start_sigma;
  stats(2) = energies_end_avg;
  stats(3) = energies_end_sigma;
  stats(4) = energies_err;

  stats(5) = overlap_inf;
  stats(6) = overlap_inf_sigma;
  stats(7) = overlap_auto;
  stats(8) = overlap_cross;
  stats(9) = overlap_check;
  stats(10) = sigma_auto;
  stats(11) = sigma_cross;
  stats(12) = sigma_check;
  stats(13) = err_cross_auto;
  stats(14) = err_cross_check;
  stats(15) = err_check_auto;

  stats(16) = Z_ratio;
  stats(17) = sumw_inv;
  stats(18) = dE_av_tot;
  return stats;
};

/**
 * @brief Compute sequence energies.
 */
void
SampleStats3D::computeEnergies(void)
{
  energies = arma::Mat<double>(reps, M, arma::fill::zeros);
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps; rep++) {
      for (int seq = 0; seq < M; seq++) {
        double E = 0;
        for (int i = 0; i < N; i++) {
          E -= params->h(samples->at(seq, i, rep), i);
          for (int j = i + 1; j < N; j++) {
            E -= params->J(i, j)(samples->at(seq, i, rep),
                                 samples->at(seq, j, rep));
          }
        }
        energies(rep, seq) = E;
      }
    }
  }
};

/**
 * @brief Compute energy statistics along MCMC trajectories.
 */
void
SampleStats3D::computeEnergiesStats(void)
{
  energies_relax = arma::Row<double>(M);
  energies_relax_sigma = arma::Row<double>(M);

  energies_start_avg = arma::mean(energies.col(0));
  energies_start_sigma = arma::stddev(energies.col(0), 1);
  energies_end_avg = arma::mean(energies.col(M - 1));
  energies_end_sigma = arma::stddev(energies.col(M - 1), 1);
  energies_err =
    sqrt((pow(energies_start_sigma, 2) + pow(energies_end_sigma, 2)) / reps);

  energies_relax = arma::mean(energies, 0);
  energies_relax_sigma = arma::stddev(energies, 1, 0);
};

/**
 * @brief Compute sequence correlations.
 */
void
SampleStats3D::computeCorrelations(void)
{
  arma::Col<double> d = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<double> d2 = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<int> count = arma::Col<int>(M, arma::fill::zeros);
  arma::Mat<double> d_mat = arma::Mat<double>(reps, M, arma::fill::zeros);
  arma::Mat<double> d2_mat = arma::Mat<double>(reps, M, arma::fill::zeros);
  arma::Mat<int> count_mat = arma::Mat<int>(reps, M, arma::fill::zeros);

  // Compute distances within replicates
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int rep = 0; rep < reps; rep++) {
      arma::Mat<int> slice = (samples->slice(rep)).t();
      for (int seq1 = 0; seq1 < M; seq1++) {
        int* seq1_ptr = slice.colptr(seq1);
        for (int seq2 = seq1 + 1; seq2 < M; seq2++) {
          int* seq2_ptr = slice.colptr(seq2);
          int id = 0;
          for (int i = 0; i < N; i++) {
            if (*(seq1_ptr + i) == *(seq2_ptr + i)) {
              id++;
            }
          }
          d_mat(rep, seq2 - seq1) += (double)id / N;
          d2_mat(rep, seq2 - seq1) += (double)id * id / (N * N);
          count_mat(rep, seq2 - seq1)++;
        }
      }
    }
  }

  d = arma::sum(d_mat, 0).t();
  d2 = arma::sum(d2_mat, 0).t();
  count = arma::sum(count_mat, 0).t();

  arma::Col<double> dinf_col = arma::Col<double>(M, arma::fill::zeros);
  arma::Col<double> dinf2_col = arma::Col<double>(M, arma::fill::zeros);

  // Compute distances between replicates
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int seq = 0; seq < M; seq++) {
      for (int rep1 = 0; rep1 < reps; rep1++) {
        for (int rep2 = rep1 + 1; rep2 < reps; rep2++) {
          int id = 0;
          for (int i = 0; i < N; i++) {
            if (samples->at(seq, i, rep1) == samples->at(seq, i, rep2)) {
              id++;
            }
          }
          dinf_col(seq) += (double)id / N;
          dinf2_col(seq) += (double)(id * id) / (N * N);
        }
      }
    }
  }

  double dinf = arma::sum(dinf_col);
  double dinf2 = arma::sum(dinf2_col);

  overlaps = arma::Col<double>(M - 2, arma::fill::zeros);
  overlaps_sigma = arma::Col<double>(M - 2, arma::fill::zeros);
  for (int i = 1; i < M - 1; i++) {
    overlaps(i - 1) = (double)d(i) / (double)count(i);
    overlaps_sigma(i - 1) =
      sqrt(1.0 / count(i)) *
      sqrt(d2(i) / (double)(count(i)) - pow(d(i) / (double)(count(i)), 2));
  }

  overlap_inf = 2.0 * dinf / (double)(reps * (reps - 1) * M);
  overlap_inf_sigma =
    sqrt(2.0 / (reps * (reps - 1) * M)) *
    sqrt(2.0 * dinf2 / (double)(reps * (reps - 1) * M) -
         pow(2.0 * dinf / (double)(reps * (reps - 1) * M), 2));

  int i_auto = 1;
  int i_check = Max((double)M / 10.0, 1.0);

  overlap_cross = (double)2.0 * dinf / (double)(reps * (reps - 1) * M);
  overlap_auto = d(i_auto) / (double)(count(i_auto));
  overlap_check = d(i_check) / (double)(count(i_check));

  sigma_cross = sqrt(2.0 * dinf2 / (double)(reps * (reps - 1) * M) -
                     pow(2.0 * dinf / (double)(reps * (reps - 1) * M), 2));
  sigma_auto = sqrt(d2(i_auto) / (double)(count(i_auto)) -
                    pow(d(i_auto) / (double)(count(i_auto)), 2));
  sigma_check = sqrt(d2(i_check) / (double)(count(i_check)) -
                     pow(d(i_check) / (double)(count(i_check)), 2));

  err_cross_auto = sqrt(pow(sigma_cross, 2) + pow(sigma_auto, 2)) / sqrt(reps);
  err_cross_check =
    sqrt(pow(sigma_cross, 2) + pow(sigma_check, 2)) / sqrt(reps);
  err_check_auto = sqrt(pow(sigma_check, 2) + pow(sigma_auto, 2)) / sqrt(reps);
};

/**
 * @brief Compute 1p and 2p statistics.
 */
void
SampleStats3D::computeSampleStats(void)
{
  frequency_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_1p_sigma = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
  frequency_2p_sigma = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
      frequency_2p_sigma(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      arma::Mat<double> n1 = arma::Mat<double>(Q, reps, arma::fill::zeros);
      arma::Col<double> n1squared = arma::Col<double>(Q, arma::fill::zeros);
      arma::Col<double> n1av = arma::Col<double>(Q, arma::fill::zeros);
      for (int rep = 0; rep < reps; rep++) {
        for (int m = 0; m < M; m++) {
          n1(samples->at(m, i, rep), rep)++;
        }
      }
      for (int aa = 0; aa < Q; aa++) {
        for (int rep = 0; rep < reps; rep++) {
          n1av(aa) += n1(aa, rep);
          n1squared(aa) += pow(n1(aa, rep), 2);
        }
        frequency_1p(aa, i) = n1av(aa) / M / reps;
        frequency_1p_sigma(aa, i) = Max(sqrt((n1squared(aa) / (M * M * reps) -
                                              pow(n1av(aa) / (M * reps), 2)) /
                                             sqrt(reps)),
                                        0);
      }
    }
  }

#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        arma::Cube<double> n2 =
          arma::Cube<double>(reps, Q, Q, arma::fill::zeros);
        arma::Mat<double> n2av = arma::Mat<double>(Q, Q, arma::fill::zeros);
        arma::Mat<double> n2squared =
          arma::Mat<double>(Q, Q, arma::fill::zeros);
        for (int rep = 0; rep < reps; rep++)
          for (int m = 0; m < M; m++) {
            n2(rep, samples->at(m, i, rep), samples->at(m, j, rep))++;
          }

        n2av = arma::sum(n2, 0) / (M * reps);
        n2squared = arma::sum(arma::pow(n2, 2), 0) / (M * reps);
        frequency_2p(i, j) = n2av;
        frequency_2p_sigma(i, j) =
          arma::pow((n2squared / M - arma::pow(n2av, 2)) / sqrt(reps), .5);
      }
    }
  }
};

/**
 * @brief Importance sampling.
 */
void
SampleStats3D::computeSampleStatsImportance(void)
{
  arma::Mat<double> p = arma::Mat<double>(reps, M, arma::fill::zeros);
  arma::Mat<double> dE = arma::Mat<double>(reps, M, arma::fill::zeros);
  arma::Col<double> dE_av = arma::Col<double>(reps, arma::fill::zeros);
  arma::Col<double> Z = arma::Col<double>(reps, arma::fill::zeros);
  arma::Col<double> Z_inv = arma::Col<double>(reps, arma::fill::zeros);
  arma::Col<double> sum = arma::Col<double>(reps, arma::fill::zeros);
  arma::Col<double> w = arma::Col<double>(reps, arma::fill::zeros);
  double Z_tot = 0;
  double Z_inv_tot = 0;
  double W = 0;
  double sumw = 0;
  Z_ratio = 0;
  sumw_inv = 0;
  dE_av_tot = 0;

  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      for (int i = 0; i < N; i++) {
        dE(rep, m) += params->h(samples->at(m, i, rep)) -
                      params_prev->h(samples->at(m, i, rep));
        for (int j = i + 1; j < N; j++) {
          dE(rep, m) +=
            params->J(i, j)(samples->at(m, i, rep), samples->at(m, j, rep)) -
            params_prev->J(i, j)(samples->at(m, i, rep),
                                 samples->at(m, j, rep));
        }
      }
      dE_av(rep) += dE(rep, m);
    }
    dE_av(rep) = dE_av(rep) / M;
    dE_av_tot += dE_av(rep);
  }

  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      p(rep, m) = exp(dE(rep, m) - dE_av(rep));
      Z(rep) += p(rep, m);
      Z_inv(rep) += 1. / p(rep, m);
    }
    for (int m = 0; m < M; m++) {
      p(rep, m) = p(rep, m) / Z(rep);
      sum(rep) += pow(p(rep, m), 2);
    }
    Z_tot += Z(rep);
    Z_inv_tot += Z_inv(rep);
    w(rep) = 1. / sum(rep);
    W += w(rep);
  }

  for (int rep = 0; rep < reps; rep++) {
    w(rep) = w(rep) / W;
    sumw += pow(w(rep), 2);
  }

  arma::Mat<int> n1 = arma::Mat<int>(Q, reps, arma::fill::zeros);
  arma::field<arma::Mat<int>> n2 = arma::field<arma::Mat<int>>(reps);
  for (int rep = 0; rep < reps; rep++) {
    n2(rep) = arma::Mat<int>(Q, Q, arma::fill::zeros);
  }

  Z_ratio = Z_tot / Z_inv_tot;
  sumw_inv = 1.0 / sumw;

  arma::Col<int> n1av = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2av = arma::Mat<int>(Q, Q, arma::fill::zeros);

  arma::Col<int> n1squared = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2squared = arma::Mat<int>(Q, Q, arma::fill::zeros);

  for (int i = 0; i < N; i++) {
    n1.zeros();
    n1av.zeros();
    n1squared.zeros();
    for (int rep = 0; rep < reps; rep++) {
      for (int m = 0; m < M; m++) {
        n1(samples->at(m, i, rep), rep) += p(rep, m);
      }
    }
    for (int aa = 0; aa < Q; aa++) {
      for (int rep = 0; rep < reps; rep++) {
        n1av(aa) += w(rep) * n1(aa, rep);
        n1squared(aa) += w(rep) * pow(n1(aa, rep), 2);
      }
      frequency_1p(aa, i) = (double)n1av(aa);
      frequency_1p_sigma(aa, i) = Max(
        sqrt(((double)n1squared(aa) - pow((double)n1av(aa), 2)) * sqrt(sumw)),
        0);
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int rep = 0; rep < reps; rep++) {
        n2(rep).zeros();
      }
      n2av.zeros();
      n2squared.zeros();
      for (int rep = 0; rep < reps; rep++)
        for (int m = 0; m < M; m++) {
          n2(rep)(samples->at(m, i, rep), samples->at(m, j, rep)) += p(rep, m);
        }

      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          for (int rep = 0; rep < reps; rep++) {
            n2av(aa1, aa2) += w(rep) * n2(rep)(aa1, aa2);
            n2squared(aa1, aa2) += w(rep) * pow(n2(rep)(aa1, aa2), 2);
          }
          frequency_2p(i, j)(aa1, aa2) = (double)n2av(aa1, aa2);
          frequency_2p_sigma(i, j)(aa1, aa2) =
            Max(sqrt(((double)n2squared(aa1, aa2) -
                      pow((double)n2av(aa1, aa2), 2)) *
                     sqrt(sumw)),
                0);
        }
      }
    }
  }
};

/**
 * @brief Write 1p frequencies in binary format.
 *
 * @param output_file file string for 1p output
 * @param output_file_sigma file string for 1p variance output
 */
void
SampleStats3D::writeFrequency1p(std::string output_file,
                                std::string output_file_sigma)
{
  frequency_1p.save(output_file, arma::arma_binary);
  frequency_1p_sigma.save(output_file_sigma, arma::arma_binary);
};

/**
 * @brief Write 1p frequencies in text format.
 *
 * @param output_file file string for 1p output
 * @param output_file_sigma file string for 1p variance output
 */
void
SampleStats3D::writeFrequency1pAscii(std::string output_file,
                                     std::string output_file_sigma)
{
  std::ofstream output_stream(output_file);
  std::ofstream output_stream_sigma(output_file_sigma);

  for (int i = 0; i < N; i++) {
    output_stream << i;
    output_stream_sigma << i;
    for (int aa = 0; aa < Q; aa++) {
      output_stream << " " << frequency_1p(aa, i);
      output_stream_sigma << " " << frequency_1p_sigma(aa, i);
    }
    output_stream << std::endl;
    output_stream_sigma << std::endl;
  }
};

/**
 * @brief Write 2p frequencies in binary format.
 *
 * @param output_file file string for 2p output
 * @param output_file_sigma file string for 2p variance output
 */
void
SampleStats3D::writeFrequency2p(std::string output_file,
                                std::string output_file_sigma)
{
  frequency_2p.save(output_file, arma::arma_binary);
  frequency_2p_sigma.save(output_file_sigma, arma::arma_binary);
};

/**
 * @brief Write 2p frequencies in text format.
 *
 * @param output_file file string for 2p output
 * @param output_file_sigma file string for 2p variance output
 */
void
SampleStats3D::writeFrequency2pAscii(std::string output_file,
                                     std::string output_file_sigma)
{
  std::ofstream output_stream(output_file);
  std::ofstream output_stream_sigma(output_file_sigma);

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      output_stream << i << " " << j;
      output_stream_sigma << i << " " << j;
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << " " << frequency_2p(i, j)(aa1, aa2);
          output_stream_sigma << " " << frequency_2p_sigma(i, j)(aa1, aa2);
        }
      }
      output_stream << std::endl;
      output_stream_sigma << std::endl;
    }
  }
};

/**
 * @brief Write samples to disk.
 *
 * @param output_file file string for output
 */
void
SampleStats3D::writeSamples(std::string output_file)
{
  std::ofstream output_stream(output_file);

  output_stream << reps * M << " " << N << " " << Q << std::endl;
  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      output_stream << samples->at(m, 0, rep);
      for (int i = 1; i < N; i++) {
        output_stream << " " << samples->at(m, i, rep);
      }
      output_stream << std::endl;
    }
  }
};

/**
 * @brief Write sampled sequence energies to disk.
 *
 * @param output_file file string for output.
 */
void
SampleStats3D::writeSampleEnergies(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      output_stream << energies(rep, m) << std::endl;
    }
  }
};

/**
 * @brief Write sequence energy trajectories to disk.
 *
 * @param output_file file string for output.
 */
void
SampleStats3D::writeSampleEnergiesRelaxation(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int m = 0; m < M; m++) {
    output_stream << m * mixing_time << " " << energies_relax(m) << " "
                  << energies_relax_sigma(m) << std::endl;
  }
};
