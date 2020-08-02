#include "msa_stats.hpp"

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <vector>
#include <random>

#include "pcg_random.hpp"

MSAStats::MSAStats(MSA* msa)
  : msa(msa)
{
  // Initialize
  N = msa->N;
  M = msa->M;
  Q = msa->Q;
  M_effective = sum(msa->sequence_weights);

  std::cout << M << " sequences" << std::endl;
  std::cout << N << " positions" << std::endl;
  std::cout << Q << " amino acids (including gaps)" << std::endl;
  std::cout << M_effective << " effective sequences" << std::endl;

  frequency_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
  rel_entropy_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  rel_entropy_pos_1p = arma::Col<double>(N, arma::fill::zeros);
  rel_entropy_grad_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  aa_background_frequencies = arma::Col<double>(Q, arma::fill::ones);

  if (Q == 21) {
    aa_background_frequencies = { 0.000, 0.073, 0.025, 0.050, 0.061, 0.042,
                                  0.072, 0.023, 0.053, 0.064, 0.089, 0.023,
                                  0.043, 0.052, 0.040, 0.052, 0.073, 0.056,
                                  0.063, 0.013, 0.033 };
  } else {
    aa_background_frequencies = aa_background_frequencies / (double)Q;
  }
  pseudocount = 0.03;

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
    arma::accu(frequency_1p % (1. - frequency_1p)) / (double)(N*Q);

  // Compute the 2p statistics
  arma::Col<double> mean_2p_var_vec = arma::Col<double>(N, arma::fill::zeros);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1)
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double* weight_ptr = msa->sequence_weights.memptr();
        frequency_2p(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);

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
    arma::accu(mean_2p_var_vec) / (double)(N * (N - 1) / 2 * Q * Q);
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

  // Use the positonal and backgrounds frequencies to estimate the relative
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

void
MSAStats::computeErrorMSA(int reps, long int seed)
{
  const bool reweight = msa->reweight;
  const double threshold = msa->threshold;
  msa_rms = arma::Col<double>(reps, arma::fill::zeros);

  int M_1 = (int)((M + 1) / 2);
  int M_2 = (int)(M / 2);

  arma::Mat<int> alignment_T = (msa->alignment).t();

  pcg32 rng;
  rng.seed(seed);

  for (int rep = 0; rep < reps; rep++) {
    std::unordered_set<int> elems;

    for (int r = M - M_1; r < M; ++r) {
      std::uniform_int_distribution<int> dist(0, r);
      int v = dist(rng);
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
    arma::Mat<int> alignment_2 = arma::Mat<int>(N, M_1, arma::fill::zeros);

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < M_1; i++) {
        alignment_1.col(i) = alignment_T.col(idx_1.at(i));
      }
    }

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < M_2; i++) {
        alignment_2.col(i) = alignment_T.col(idx_2.at(i));
      }
    }

    MSA msa_1 = MSA(alignment_1.t(), M_1, N, Q, reweight, threshold);
    MSA msa_2 = MSA(alignment_2.t(), M_2, N, Q, reweight, threshold);

    double M_1_effective = sum(msa_1.sequence_weights);
    double M_2_effective = sum(msa_2.sequence_weights);

    arma::Mat<double> msa_1_frequency_1p =
      arma::Mat<double>(Q, N, arma::fill::zeros);
    arma::field<arma::Mat<double>> msa_1_frequency_2p =
      arma::field<arma::Mat<double>>(N, N);

    arma::Mat<double> msa_2_frequency_1p =
      arma::Mat<double>(Q, N, arma::fill::zeros);
    arma::field<arma::Mat<double>> msa_2_frequency_2p =
      arma::field<arma::Mat<double>>(N, N);

    // Compute the frequecies (1p statistics) for amino acids (and gaps) for each
    // position. Use pointers to make things speedier.
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
#pragma omp for schedule(dynamic,1)
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
#pragma omp for schedule(dynamic,1)
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

    double error_1p =
      arma::accu(arma::pow(msa_1_frequency_1p - msa_2_frequency_1p, 2));
    error_1p = sqrt(error_1p / (N * Q));

    arma::Col<double> error_2p_vec = arma::Col<double>(N, arma::fill::zeros);
#pragma omp parallel
    {
#pragma omp for schedule(dynamic,1)
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          error_2p_vec(i) = arma::accu(
            arma::pow(msa_1_frequency_2p(i, j) - msa_2_frequency_2p(i, j), 2));
        }
      }
    }
    double error_2p = sqrt(arma::accu(error_2p_vec) / (N * (N - 1) * Q * Q));

    if (reweight) {
      msa_rms(rep) = (error_1p + error_2p) /
                     sqrt(M_effective / (M_1_effective + M_2_effective) * 2);
    } else {
      msa_rms(rep) = (error_1p + error_2p) / sqrt(2);
    }
  }

  return;
}

double
MSAStats::getQ(void)
{
  return Q;
};

double
MSAStats::getM(void)
{
  return M;
};

double
MSAStats::getN(void)
{
  return N;
};

double
MSAStats::getEffectiveM(void)
{
  return M_effective;
};

void
MSAStats::writeRelEntropyGradient(std::string output_file)
{
  rel_entropy_grad_1p.save(output_file, arma::arma_binary);
};

void
MSAStats::writeRelEntropyGradientAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < N; i++) {
    output_stream << i;
    for (int aa = 0; aa < Q; aa++) {
      output_stream << " " << rel_entropy_grad_1p(aa, i);
    }
    output_stream << std::endl;
  }
};

// void
// MSAStats::writeRelEntropyPos(std::string output_file)
// {
//   rel_entropy_pos_1p.save(output_file, arma::arma_binary);
// };

void
MSAStats::writeRelEntropyPosAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < N; i++) {
    output_stream << " " << rel_entropy_pos_1p(i) << std::endl;
  }
};

void
MSAStats::writeRelEntropy(std::string output_file)
{
  rel_entropy_1p.save(output_file, arma::arma_binary);
};

void
MSAStats::writeRelEntropyAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < N; i++) {
    output_stream << i;
    for (int aa = 0; aa < Q; aa++) {
      output_stream << " " << rel_entropy_1p(aa, i);
    }
    output_stream << std::endl;
  }
};

void
MSAStats::writeFrequency1p(std::string output_file)
{
  frequency_1p.save(output_file, arma::arma_binary);
};

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

void
MSAStats::writeFrequency2p(std::string output_file)
{
  frequency_2p.save(output_file, arma::arma_binary);
};

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
