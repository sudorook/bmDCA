#include "model.hpp"

#include <armadillo>

#include "utils.hpp"

Model::Model(std::string parameters_file,
             std::string parameters_prev_file,
             std::string gradient_file,
             std::string gradient_prev_file,
             std::string moment1_file,
             std::string moment2_file)
{
  params = loadPottsModelAscii(parameters_file);
  params_prev = loadPottsModelAscii(parameters_prev_file);
  gradient = loadPottsModelAscii(gradient_file);
  gradient_prev = loadPottsModelAscii(gradient_prev_file);

  N = params.h.n_cols;
  Q = params.h.n_rows;
}

Model::Model(std::string parameters_file_h,
             std::string parameters_file_J,
             std::string parameters_prev_file_h,
             std::string parameters_prev_file_J,
             std::string gradient_file_h,
             std::string gradient_file_J,
             std::string gradient_prev_file_h,
             std::string gradient_prev_file_J,
             std::string moment1_file_h,
             std::string moment1_file_J,
             std::string moment2_file_h,
             std::string moment2_file_J)
{
  params = loadPottsModel(parameters_file_h, parameters_file_J);
  params_prev = loadPottsModel(parameters_prev_file_h, parameters_prev_file_J);
  gradient = loadPottsModel(gradient_file_h, gradient_file_J);
  gradient_prev = loadPottsModel(gradient_prev_file_h, gradient_prev_file_J);
  moment1 = loadPottsModel(moment1_file_h, moment1_file_J);
  moment2 = loadPottsModel(moment2_file_h, moment2_file_J);

  N = params.h.n_cols;
  Q = params.h.n_rows;
}

Model::Model(MSAStats* msa_stats, bool init_gradient)
{
  N = msa_stats->getN();
  Q = msa_stats->getQ();
  double pseudocount = 1. / msa_stats->getEffectiveM();

  // Initialize the parameters J and h
  params.J = arma::field<arma::Mat<double>>(N, N);
  params_prev.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
      params_prev.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  params.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params_prev.h = arma::Mat<double>(Q, N, arma::fill::zeros);

  if (init_gradient) {
    double avg;
    double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      avg = 0;
      freq_ptr = msa_stats->frequency_1p.colptr(i);
      for (int aa = 0; aa < Q; aa++) {
        avg +=
          log((1. - pseudocount) * (*(freq_ptr + aa)) + pseudocount * (1. / Q));
      }
      for (int aa = 0; aa < Q; aa++) {
        params.h(aa, i) = log((1. - pseudocount) * (*(freq_ptr + aa)) +
                              pseudocount * (1. / Q)) -
                          avg / Q;
      }
    }
  }

  // Initialize the moments
  moment1.h = arma::Mat<double>(Q, N);
  moment2.h = arma::Mat<double>(Q, N);
  moment1.h.fill(arma::fill::zeros);
  moment2.h.fill(arma::fill::zeros);

  moment1.J = arma::field<arma::Mat<double>>(N, N);
  moment2.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      moment1.J(i, j) = arma::Mat<double>(Q, Q);
      moment2.J(i, j) = arma::Mat<double>(Q, Q);
      moment1.J(i, j).fill(arma::fill::zeros);
      moment2.J(i, j).fill(arma::fill::zeros);
    }
  }

  // Initialize the gradient
  gradient.J = arma::field<arma::Mat<double>>(N, N);
  gradient_prev.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
      gradient_prev.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }
  gradient.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  gradient_prev.h = arma::Mat<double>(Q, N, arma::fill::zeros);
};

void
Model::resetModel(MSAStats* msa_stats, bool init_gradient)
{
  N = msa_stats->getN();
  Q = msa_stats->getQ();
  double pseudocount = 1. / msa_stats->getEffectiveM();

  // Initialize the parameters J and h
  params.h.zeros();
  params_prev.h.zeros();
  if (init_gradient) {
    double avg;
    double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      avg = 0;
      freq_ptr = msa_stats->frequency_1p.colptr(i);
      for (int aa = 0; aa < Q; aa++) {
        avg +=
          log((1. - pseudocount) * (*(freq_ptr + aa)) + pseudocount * (1. / Q));
      }
      for (int aa = 0; aa < Q; aa++) {
        params.h(aa, i) = log((1. - pseudocount) * (*(freq_ptr + aa)) +
                              pseudocount * (1. / Q)) -
                          avg / Q;
      }
    }
  }
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J(i, j).zeros();
      params_prev.J(i, j).zeros();
    }
  }

  // Initialize the gradient
  gradient.h.zeros();
  gradient_prev.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J(i, j).zeros();
      gradient_prev.J(i, j).zeros();
    }
  }

  // Initialize the moments
  moment1.h.zeros();
  moment2.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      moment1.J(i, j).zeros();
      moment2.J(i, j).zeros();
    }
  }
};

void
Model::writeParams(std::string output_file_h, std::string output_file_J)
{
  params.h.save(output_file_h, arma::arma_binary);
  params.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeParamsAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = params.h.n_cols;
  int Q = params.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << params.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << params.h(aa, i)
                    << std::endl;
    }
  }
};

void
Model::writeParamsPrevious(std::string output_file_h, std::string output_file_J)
{
  params_prev.h.save(output_file_h, arma::arma_binary);
  params_prev.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeParamsPreviousAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = params_prev.h.n_cols;
  int Q = params_prev.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << params_prev.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << params_prev.h(aa, i)
                    << std::endl;
    }
  }
};

void
Model::writeMoment1(std::string output_file_h, std::string output_file_J)
{
  moment1.h.save(output_file_h, arma::arma_binary);
  moment1.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeMoment2(std::string output_file_h, std::string output_file_J)
{
  moment2.h.save(output_file_h, arma::arma_binary);
  moment2.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeMoment1Ascii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = moment1.h.n_cols;
  int Q = moment1.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << moment1.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << moment1.h(aa, i)
                    << std::endl;
    }
  }
};

void
Model::writeMoment2Ascii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = moment1.h.n_cols;
  int Q = moment1.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << moment1.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << moment1.h(aa, i)
                    << std::endl;
    }
  }
};

void
Model::writeGradient(std::string output_file_h, std::string output_file_J)
{
  gradient.h.save(output_file_h, arma::arma_binary);
  gradient.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeGradientAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = gradient.h.n_cols;
  int Q = gradient.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << gradient.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << gradient.h(aa, i)
                    << std::endl;
    }
  }
};

void
Model::writeGradientPrevious(std::string output_file_h,
                             std::string output_file_J)
{
  gradient_prev.h.save(output_file_h, arma::arma_binary);
  gradient_prev.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeGradientPreviousAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = gradient_prev.h.n_cols;
  int Q = gradient_prev.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << gradient_prev.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << gradient_prev.h(aa, i)
                    << std::endl;
    }
  }
};
