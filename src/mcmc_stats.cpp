#include "mcmc_stats.hpp"

#include <armadillo>
#include <iostream>

#include "utils.hpp"

MCMCStats::MCMCStats(arma::Cube<int>* s, potts_model* p)
{
  M = s->n_rows;
  N = s->n_cols;
  reps = s->n_slices;
  Q = p->h.n_rows;

  samples = s;
  params = p;

  computeEnergies();
};

void
MCMCStats::updateData(arma::Cube<int>* s, potts_model* p)
{
  samples = s;
  params = p;

  computeEnergies();
};

void
MCMCStats::computeEnergies(void)
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

void
MCMCStats::computeEnergiesStats(void)
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

void
MCMCStats::writeEnergyStats(std::string output_file_start,
                            std::string output_file_end,
                            std::string output_file_cfr,
                            std::string output_file_cfr_err)
{

  std::ofstream output_stream_start(output_file_start);
  std::ofstream output_stream_end(output_file_end);
  std::ofstream output_stream_cfr(output_file_cfr);
  std::ofstream output_stream_cfr_err(output_file_cfr_err);

  for (int rep = 0; rep < reps; rep++) {
    output_stream_start << energies(rep, 0) << std::endl;
    output_stream_end << energies(rep, M - 1) << std::endl;
  }

  output_stream_cfr << reps << " " << energies_start_avg << " "
                    << energies_start_sigma << " " << reps << " "
                    << energies_end_avg << " " << energies_end_sigma
                    << std::endl;
  output_stream_cfr_err << energies_start_avg << " " << energies_end_avg << " "
                        << energies_err << std::endl;
};

void
MCMCStats::computeCorrelations(void)
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

void
MCMCStats::writeCorrelationsStats(std::string overlap_file,
                                  std::string overlap_inf_file,
                                  std::string ergo_file)
{
  std::ofstream output_stream_overlap(overlap_file);
  std::ofstream output_stream_overlap_inf(overlap_inf_file);
  std::ofstream output_stream_ergo(ergo_file);

  for (int i = 0; i < M - 2; i++) {
    output_stream_overlap << i << " " << overlaps(i) << " " << overlaps_sigma(i)
                          << std::endl;
  }

  output_stream_overlap_inf << "0 " << overlap_inf << " " << overlap_inf_sigma
                            << std::endl;
  output_stream_overlap_inf << M << " " << overlap_inf << " "
                            << overlap_inf_sigma << std::endl;

  output_stream_ergo << overlap_auto << " " << overlap_check << " "
                     << overlap_cross << " " << sigma_auto << " " << sigma_check
                     << " " << sigma_cross << " " << err_cross_auto << " "
                     << err_cross_check << " " << err_check_auto << std::endl;
};

std::vector<double>
MCMCStats::getEnergiesStats(void)
{
  std::vector<double> values = { energies_start_avg,
                                 energies_start_sigma,
                                 energies_end_avg,
                                 energies_end_sigma,
                                 energies_err };
  return values;
};

std::vector<double>
MCMCStats::getCorrelationsStats(void)
{
  std::vector<double> values = {
    overlap_inf,    overlap_inf_sigma, overlap_auto,  overlap_cross,
    overlap_check,  sigma_auto,        sigma_cross,   sigma_check,
    err_cross_auto, err_cross_check,   err_check_auto
  };
  return values;
};

void
MCMCStats::computeSampleStats(void)
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

void
MCMCStats::computeSampleStatsImportance(potts_model* cur, potts_model* prev)
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
        dE(rep, m) +=
          cur->h(samples->at(m, i, rep)) - prev->h(samples->at(m, i, rep));
        for (int j = i + 1; j < N; j++) {
          dE(rep, m) +=
            cur->J(i, j)(samples->at(m, i, rep), samples->at(m, j, rep)) -
            prev->J(i, j)(samples->at(m, i, rep), samples->at(m, j, rep));
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

void
MCMCStats::writeFrequency1p(std::string output_file,
                            std::string output_file_sigma)
{
  frequency_1p.save(output_file, arma::arma_binary);
  frequency_1p_sigma.save(output_file_sigma, arma::arma_binary);
};

void
MCMCStats::writeFrequency1pAscii(std::string output_file,
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

void
MCMCStats::writeFrequency2p(std::string output_file,
                            std::string output_file_sigma)
{
  frequency_2p.save(output_file, arma::arma_binary);
  frequency_2p_sigma.save(output_file_sigma, arma::arma_binary);
};

void
MCMCStats::writeFrequency2pAscii(std::string output_file,
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

void
MCMCStats::writeSamples(std::string output_file)
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

void
MCMCStats::writeSampleEnergies(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      output_stream << energies(rep, m) << std::endl;
    }
  }
};

void
MCMCStats::writeSampleEnergiesRelaxation(std::string output_file, int t_wait)
{
  std::ofstream output_stream(output_file);
  for (int m = 0; m < M; m++) {
    output_stream << m * t_wait << " " << energies_relax(m) << " "
                  << energies_relax_sigma(m) << std::endl;
  }
};
