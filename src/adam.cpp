// #include "model.hpp"
#include "adam.hpp"

#include <armadillo>

#include "utils.hpp"

Adam::Adam()
  : Model(){};

void
Adam::loadHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_adam_section = false;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
        reading_adam_section = false;
        continue;
      } else if (line[0] == '[') {
        if (line == "[[adam]]") {
          reading_adam_section = true;
          continue;
        } else {
          reading_adam_section = false;
          continue;
        }
      }
      if (reading_adam_section) {
        auto delim_pos = line.find("=");
        auto key = line.substr(0, delim_pos);
        auto value = line.substr(delim_pos + 1);
        setHyperparameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

bool
Adam::compareHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_adam_section = false;
  bool all_same = true;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
        continue;
      } else if (line[0] == '[') {
        if (line == "[[adam]]") {
          reading_adam_section = true;
          continue;
        } else {
          reading_adam_section = false;
          continue;
        }
      }
      if (reading_adam_section) {
        auto delim_pos = line.find("=");
        auto key = line.substr(0, delim_pos);
        auto value = line.substr(delim_pos + 1);
        all_same = all_same & compareHyperparameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return all_same;
};

void
Adam::setHyperparameter(std::string key, std::string value)
{
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_reg_h") {
    lambda_reg_h = std::stod(value);
  } else if (key == "lambda_reg_J") {
    lambda_reg_J = std::stod(value);
  } else if (key == "alpha_reg") {
    alpha_reg = std::stod(value);
  } else if (key == "initial_params") {
    initial_params = value;
  } else if (key == "set_zero_gauge") {
    if (value.size() == 1) {
      set_zero_gauge = (std::stoi(value) == 1);
    } else {
      set_zero_gauge = (value == "true");
    }
  } else if (key == "learn_rate_h") {
    learn_rate_h = std::stod(value);
  } else if (key == "learn_rate_J") {
    learn_rate_J = std::stod(value);
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

void
Adam::checkHyperparameters(void)
{
  // if ((anneal_period < 1) & (anneal_schedule == "cos")) {
  //   std::cerr << "ERROR: period " << anneal_period
  //             << " invalid for 'cos' schedule." << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
}

void
Adam::writeHyperparameters(std::string output_file, bool append)
{
  std::ofstream stream;
  if (append) {
    stream.open(output_file, std::ofstream::out | std::ofstream::app);
  } else {
    stream.open(output_file, std::ofstream::out | std::ofstream::trunc);
  }
  // Header
  stream << std::endl;
  stream << "[[adam]]" << std::endl;

  // BM settings
  stream << "lambda_reg_h=" << lambda_reg_h << std::endl;
  stream << "lambda_reg_J=" << lambda_reg_J << std::endl;
  stream << "alpha_reg=" << alpha_reg << std::endl;
  stream << "initial_params=" << initial_params << std::endl;
  stream << "set_zero_gauge=" << set_zero_gauge << std::endl;
  stream << "learn_rate_h=" << learn_rate_h << std::endl;
  stream << "learn_rate_J=" << learn_rate_J << std::endl;

  stream.close();
};

bool
Adam::compareHyperparameter(std::string key, std::string value)
{
  bool same = true;
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_reg_h") {
    same = same & (lambda_reg_h == std::stod(value));
  } else if (key == "lambda_reg_J") {
    same = same & (lambda_reg_J == std::stod(value));
  } else if (key == "alpha_reg") {
    same = same & (alpha_reg == std::stod(value));
  } else if (key == "initial_params") {
    same = same & (initial_params == value);
  } else if (key == "set_zero_gauge") {
    if (value.size() == 1) {
      same = same & (set_zero_gauge == (std::stoi(value) == 1));
    } else {
      same = same & (set_zero_gauge == (value == "true"));
    }
  } else if (key == "learn_rate_h") {
    same = same & (learn_rate_h == std::stod(value));
  } else if (key == "learn_rate_J") {
    same = same & (learn_rate_J == std::stod(value));
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return same;
};

bool
Adam::isValidStep(int step, bool output_binary)
{
  bool valid = false;
  if (output_binary) {
    if (checkFileExists("parameters_h_" + std::to_string(step) + ".bin") &
        checkFileExists("parameters_J_" + std::to_string(step) + ".bin") &
        checkFileExists("parameters_h_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("parameters_J_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("gradients_h_" + std::to_string(step) + ".bin") &
        checkFileExists("gradients_J_" + std::to_string(step) + ".bin") &
        checkFileExists("moment1_h_" + std::to_string(step) + ".bin") &
        checkFileExists("moment1_J_" + std::to_string(step) + ".bin") &
        checkFileExists("moment2_h_" + std::to_string(step) + ".bin") &
        checkFileExists("moment2_J_" + std::to_string(step) + ".bin")) {
      valid = true;
    }
  } else {
    if (checkFileExists("parameters_" + std::to_string(step) + ".txt") &
        checkFileExists("parameters_" + std::to_string(step - 1) + ".txt") &
        checkFileExists("gradients_" + std::to_string(step) + ".txt") &
        checkFileExists("moment1_" + std::to_string(step) + ".txt") &
        checkFileExists("moment2_" + std::to_string(step) + ".txt")) {
      valid = true;
    }
  }
  return valid;
};

void
Adam::initialize(void)
{
  double pseudocount = 1. / training->getEffectiveM();

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

  if (initial_params == "profile") {
    double avg;
    double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      avg = 0;
      freq_ptr = training->frequency_1p.colptr(i);
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
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }
  gradient.h = arma::Mat<double>(Q, N, arma::fill::zeros);
};

void
Adam::reset()
{
  double pseudocount = 1. / training->getEffectiveM();

  params.h.zeros();
  params_prev.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J(i, j).zeros();
      params_prev.J(i, j).zeros();
    }
  }

  if (initial_params == "profile") {
    double avg;
    double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      avg = 0;
      freq_ptr = training->frequency_1p.colptr(i);
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

  gradient.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J(i, j).zeros();
    }
  }

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
Adam::restore(int step, bool output_binary)
{
  if (output_binary) {
    std::string param_h_file = "parameters_h_" + std::to_string(step) + ".bin";
    std::string param_J_file = "parameters_J_" + std::to_string(step) + ".bin";
    params = loadPottsModel(param_h_file, param_J_file);

    std::string param_prev_h_file =
      "parameters_h_" + std::to_string(step - 1) + ".bin";
    std::string param_prev_J_file =
      "parameters_J_" + std::to_string(step - 1) + ".bin";
    params_prev = loadPottsModel(param_prev_h_file, param_prev_J_file);

    std::string grad_h_file = "gradients_h_" + std::to_string(step) + ".bin";
    std::string grad_J_file = "gradients_J_" + std::to_string(step) + ".bin";
    gradient = loadPottsModel(grad_h_file, grad_J_file);

    std::string moment1_h_file = "moment1_h_" + std::to_string(step) + ".bin";
    std::string moment1_J_file = "moment1_J_" + std::to_string(step) + ".bin";
    moment1 = loadPottsModel(moment1_h_file, moment1_J_file);

    std::string moment2_h_file = "moment2_h_" + std::to_string(step) + ".bin";
    std::string moment2_J_file = "moment2_J_" + std::to_string(step) + ".bin";
    moment2 = loadPottsModel(moment2_h_file, moment2_J_file);
  } else {
    std::string param_file = "parameters_" + std::to_string(step) + ".txt";
    params = loadPottsModelAscii(param_file);

    std::string param_prev_file =
      "parameters_" + std::to_string(step - 1) + ".txt";
    params_prev = loadPottsModelAscii(param_prev_file);

    std::string grad_file = "gradients_" + std::to_string(step) + ".txt";
    gradient = loadPottsModelAscii(grad_file);

    std::string moment1_file = "moment1_" + std::to_string(step) + ".txt";
    moment1 = loadPottsModelAscii(moment1_file);

    std::string moment2_file = "moment2_" + std::to_string(step) + ".txt";
    moment2 = loadPottsModelAscii(moment2_file);
  }
};

void
Adam::update(void)
{
  params_prev = params;
  updateGradients();
  updateMoments();
  updateParameters();
  if (set_zero_gauge) {
    setZeroGauge();
  }
};

void
Adam::updateGradients(void)
{
  train_error_1p = 0;
  train_error_2p = 0;
  validation_error_1p = 0;
  validation_error_2p = 0;

  // Compute gradient
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      double delta =
        samples->frequency_1p(aa, i) - training->frequency_1p(aa, i) +
        lambda_reg_h *
          (alpha_reg * params.h(aa, i) +
           (1. - alpha_reg) * (0.5 - (double)std::signbit(params.h(aa, i))));
      train_error_1p += pow(delta, 2);
      gradient.h(aa, i) = -delta;
    }
  }
  train_error_1p = sqrt(train_error_1p / (N * Q));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          double delta =
            -(training->frequency_2p(i, j)(aa1, aa2) -
              samples->frequency_2p(i, j)(aa1, aa2) -
              lambda_reg_J *
                (alpha_reg * params.J(i, j)(aa1, aa2) +
                 (1. - alpha_reg) *
                   (0.5 - (double)std::signbit(params.J(i, j)(aa1, aa2)))));
          train_error_2p += pow(delta, 2);
          gradient.J(i, j)(aa1, aa2) = -delta;
        }
      }
    }
  }
  train_error_2p = sqrt(train_error_2p / (N * (N - 1) * Q * Q / 2));

  if (validate) {
    for (int i = 0; i < N; i++) {
      for (int aa = 0; aa < Q; aa++) {
        double delta =
          samples->frequency_1p(aa, i) - validation->frequency_1p(aa, i) +
          lambda_reg_h *
            (alpha_reg * params.h(aa, i) +
             (1. - alpha_reg) * (0.5 - (double)std::signbit(params.h(aa, i))));
        validation_error_1p += pow(delta, 2);
      }
    }
    validation_error_1p = sqrt(validation_error_1p / (N * Q));

    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            double delta =
              -(validation->frequency_2p(i, j)(aa1, aa2) -
                samples->frequency_2p(i, j)(aa1, aa2) -
                lambda_reg_J *
                  (alpha_reg * params.J(i, j)(aa1, aa2) +
                   (1. - alpha_reg) *
                     (0.5 - (double)std::signbit(params.J(i, j)(aa1, aa2)))));
            validation_error_2p += pow(delta, 2);
          }
        }
      }
    }
    validation_error_2p = sqrt(validation_error_2p / (N * (N - 1) * Q * Q / 2));
  }
  return;
};

void
Adam::updateMoments(void)
{
  double beta1 = 0.9;
  double beta2 = 0.999;

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          moment1.J(i, j)(a, b) = (1 - beta1) * gradient.J(i, j)(a, b) +
                                  beta1 * moment1.J(i, j)(a, b);
          moment2.J(i, j)(a, b) = (1 - beta2) * pow(gradient.J(i, j)(a, b), 2) +
                                  beta2 * moment2.J(i, j)(a, b);
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      moment1.h(a, i) =
        (1 - beta1) * gradient.h(a, i) + beta1 * moment1.h(a, i);
      moment2.h(a, i) =
        (1 - beta2) * pow(gradient.h(a, i), 2) + beta2 * moment2.h(a, i);
    }
  }
};

void
Adam::updateParameters(void)
{
  double beta1 = 0.9;
  double beta2 = 0.999;
  double beta1_t = pow(beta1, step);
  double beta2_t = pow(beta2, step);

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          params.J(i, j)(a, b) +=
            learn_rate_J *
            (moment1.J(i, j)(a, b) / (1 - beta1_t) /
             (sqrt(moment2.J(i, j)(a, b) / (1 - beta2_t)) + 0.00000001));
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      params.h(a, i) +=
        learn_rate_h * (moment1.h(a, i) / (1 - beta1_t) /
                        (sqrt(moment2.h(a, i) / (1 - beta2_t)) + 0.00000001));
    }
  }
};

void
Adam::writeData(std::string str, bool output_binary)
{
  if (output_binary) {
    std::string param_h_file = "parameters_h_" + str + ".bin";
    std::string param_J_file = "parameters_J_" + str + ".bin";
    writeParams(param_h_file, param_J_file);

    std::string grad_h_file = "gradients_h_" + str + ".bin";
    std::string grad_J_file = "gradients_J_" + str + ".bin";
    writeGradient(grad_h_file, grad_J_file);

    std::string moment1_h_file = "moment1_h_" + str + ".bin";
    std::string moment1_J_file = "moment1_J_" + str + ".bin";
    writeMoment1(moment1_h_file, moment1_J_file);

    std::string moment2_h_file = "moment2_h_" + str + ".bin";
    std::string moment2_J_file = "moment2_J_" + str + ".bin";
    writeMoment2(moment2_h_file, moment2_J_file);
  } else {
    std::string param_file = "parameters_" + str + ".txt";
    writeParamsAscii(param_file);

    std::string grad_file = "gradients_" + str + ".txt";
    writeGradientAscii(grad_file);

    std::string moment1_file = "moment1_" + str + ".txt";
    writeMoment1Ascii(moment1_file);

    std::string moment2_file = "moment2_" + str + ".txt";
    writeMoment2Ascii(moment2_file);
  }
};

void
Adam::deleteStep(int step, bool output_binary)
{
  std::string file;
  if (output_binary) {
    file = "parameters_h_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "parameters_J_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "gradients_h_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "gradients_J_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment1_h_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment1_J_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment2_h_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment2_J_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);
  } else {
    file = "parameters_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "gradients_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment1_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "moment2_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);
  }
};

void
Adam::writeStep(int step, bool output_binary)
{
  if (output_binary) {
    std::string param_h_file = "parameters_h_" + std::to_string(step) + ".bin";
    std::string param_J_file = "parameters_J_" + std::to_string(step) + ".bin";
    writeParams(param_h_file, param_J_file);

    std::string param_prev_h_file =
      "parameters_h_" + std::to_string(step - 1) + ".bin";
    std::string param_prev_J_file =
      "parameters_J_" + std::to_string(step - 1) + ".bin";
    writeParamsPrevious(param_prev_h_file, param_prev_J_file);

    std::string grad_h_file = "gradients_h_" + std::to_string(step) + ".bin";
    std::string grad_J_file = "gradients_J_" + std::to_string(step) + ".bin";
    writeGradient(grad_h_file, grad_J_file);

    std::string moment1_h_file = "moment1_h_" + std::to_string(step) + ".bin";
    std::string moment1_J_file = "moment1_J_" + std::to_string(step) + ".bin";
    writeMoment1(moment1_h_file, moment1_J_file);

    std::string moment2_h_file = "moment2_h_" + std::to_string(step) + ".bin";
    std::string moment2_J_file = "moment2_J_" + std::to_string(step) + ".bin";
    writeMoment2(moment2_h_file, moment2_J_file);
  } else {
    std::string param_file = "parameters_" + std::to_string(step) + ".txt";
    writeParamsAscii(param_file);

    std::string param_prev_file =
      "parameters_prev_" + std::to_string(step) + ".txt";
    writeParamsPreviousAscii(param_prev_file);

    std::string grad_file = "gradients_" + std::to_string(step) + ".txt";
    writeGradientAscii(grad_file);

    std::string moment1_file = "moment1_" + std::to_string(step) + ".txt";
    writeMoment1Ascii(moment1_file);

    std::string moment2_file = "moment2_" + std::to_string(step) + ".txt";
    writeMoment2Ascii(moment2_file);
  }
};

void
Adam::writeParams(std::string output_file_h, std::string output_file_J)
{
  params.h.save(output_file_h, arma::arma_binary);
  params.J.save(output_file_J, arma::arma_binary);
};

void
Adam::writeParamsAscii(std::string output_file)
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
Adam::writeParamsPrevious(std::string output_file_h, std::string output_file_J)
{
  params_prev.h.save(output_file_h, arma::arma_binary);
  params_prev.J.save(output_file_J, arma::arma_binary);
};

void
Adam::writeParamsPreviousAscii(std::string output_file)
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
Adam::writeMoment1(std::string output_file_h, std::string output_file_J)
{
  moment1.h.save(output_file_h, arma::arma_binary);
  moment1.J.save(output_file_J, arma::arma_binary);
};

void
Adam::writeMoment2(std::string output_file_h, std::string output_file_J)
{
  moment2.h.save(output_file_h, arma::arma_binary);
  moment2.J.save(output_file_J, arma::arma_binary);
};

void
Adam::writeMoment1Ascii(std::string output_file)
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
Adam::writeMoment2Ascii(std::string output_file)
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
Adam::writeGradient(std::string output_file_h, std::string output_file_J)
{
  gradient.h.save(output_file_h, arma::arma_binary);
  gradient.J.save(output_file_J, arma::arma_binary);
};

void
Adam::writeGradientAscii(std::string output_file)
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
