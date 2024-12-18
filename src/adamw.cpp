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

#include "adamw.hpp"

#include <armadillo>

#include "utils.hpp"

/**
 * @brief AdamW constructor.
 */
AdamW::AdamW()
  : Model(){};

/**
 * @brief Load model hyperparameters from config file.
 *
 * @param file_name config file string
 */
void
AdamW::loadHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  if (file.is_open()) {
    std::string line;
    bool reading_adamw_section = false;
    while (std::getline(file, line)) {
      if (line.empty()) {
        reading_adamw_section = false;
        continue;
      } else if (line[0] == '#') {
        reading_adamw_section = false;
        continue;
      } else if (line[0] == '[') {
        if (line == "[[adamw]]") {
          reading_adamw_section = true;
          continue;
        } else {
          reading_adamw_section = false;
          continue;
        }
      }
      if (reading_adamw_section) {
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

/**
 * @brief Check that hyperparameters in config file are same as stored values.
 *
 * @param file_name config file string
 *
 * @return (bool) flag that all hyperparameters are equivalent
 */
bool
AdamW::compareHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool all_same = true;
  if (file.is_open()) {
    std::string line;
    bool reading_adamw_section = false;
    while (std::getline(file, line)) {
      if (line.empty()) {
        continue;
      } else if (line[0] == '#') {
        continue;
      } else if (line[0] == '[') {
        if (line == "[[adamw]]") {
          reading_adamw_section = true;
          continue;
        } else {
          reading_adamw_section = false;
          continue;
        }
      }
      if (reading_adamw_section) {
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

/**
 * @brief Set a hyperparameter to a specific value.
 *
 * @param key hyperparameter to set
 * @param value value at which to set hyperparameter
 */
void
AdamW::setHyperparameter(std::string key, std::string value)
{
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_decay_h") {
    lambda_decay_h = std::stod(value);
  } else if (key == "lambda_decay_J") {
    lambda_decay_J = std::stod(value);
  } else if (key == "initial_params") {
    initial_params = value;
  } else if (key == "set_zero_gauge") {
    if (value.size() == 1) {
      set_zero_gauge = (std::stoi(value) == 1);
    } else {
      set_zero_gauge = (value == "true");
    }
  } else if (key == "allow_gap_couplings") {
    if (value.size() == 1) {
      allow_gap_couplings = (std::stoi(value) == 1);
    } else {
      allow_gap_couplings = (value == "true");
    }
  } else if (key == "learn_rate_h") {
    learn_rate_h = std::stod(value);
  } else if (key == "learn_rate_J") {
    learn_rate_J = std::stod(value);
  } else if (key == "eta_min") {
    eta_min = std::stod(value);
  } else if (key == "eta_max") {
    eta_max = std::stod(value);
  } else if (key == "anneal_schedule") {
    anneal_schedule = value;
  } else if (key == "anneal_scale") {
    anneal_scale = std::stod(value);
  } else if (key == "anneal_period") {
    anneal_period = std::stoi(value);
  } else if (key == "anneal_warm") {
    anneal_warm = std::stoi(value);
  } else if (key == "anneal_hot") {
    anneal_hot = std::stoi(value);
  } else if (key == "anneal_cool") {
    anneal_cool = std::stoi(value);
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

/**
 * @brief Check additional constraints for hyperparameter values.
 *
 * Function for checking any constraints on the possible values for the model
 * hyperparameters. Empty by default.
 */
void
AdamW::checkHyperparameters(void)
{
  if ((anneal_period < 1) & (anneal_schedule == "cos")) {
    std::cerr << "ERROR: period " << anneal_period
              << " invalid for 'cos' schedule." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

/**
 * @brief Write stored hyperparameters to file
 *
 * @param output_file config file string
 * @param append flag for whether or not to append the hyperparameters to the
 * file or overwrite it
 */
void
AdamW::writeHyperparameters(std::string output_file, bool append)
{
  std::ofstream stream;
  if (append) {
    stream.open(output_file, std::ofstream::out | std::ofstream::app);
  } else {
    stream.open(output_file, std::ofstream::out | std::ofstream::trunc);
  }
  // Header
  stream << std::endl;
  stream << "[[adamw]]" << std::endl;

  // BM settings
  stream << "lambda_decay_h=" << lambda_decay_h << std::endl;
  stream << "lambda_decay_J=" << lambda_decay_J << std::endl;
  stream << "initial_params=" << initial_params << std::endl;
  stream << "set_zero_gauge=" << set_zero_gauge << std::endl;
  stream << "allow_gap_couplings=" << allow_gap_couplings << std::endl;
  stream << "learn_rate_h=" << learn_rate_h << std::endl;
  stream << "learn_rate_J=" << learn_rate_J << std::endl;
  stream << "eta_min=" << eta_min << std::endl;
  stream << "eta_max=" << eta_max << std::endl;
  stream << "anneal_schedule=" << anneal_schedule << std::endl;
  stream << "anneal_scale=" << anneal_scale << std::endl;
  stream << "anneal_period=" << anneal_period << std::endl;
  stream << "anneal_warm=" << anneal_warm << std::endl;
  stream << "anneal_hot=" << anneal_hot << std::endl;
  stream << "anneal_cool=" << anneal_cool << std::endl;

  stream.close();
};

/**
 * @brief Compare the value of a specific hyperparameter against a given value.
 *
 * @param key name of the hyperparameter to check
 * @param value value against which to check the 'key' hyperparameter.
 *
 * @return (bool) flag for whether the stored and given values are equal
 */
bool
AdamW::compareHyperparameter(std::string key, std::string value)
{
  bool same = true;
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_decay_h") {
    same = same & (lambda_decay_h == std::stod(value));
  } else if (key == "lambda_decay_J") {
    same = same & (lambda_decay_J == std::stod(value));
  } else if (key == "initial_params") {
    same = same & (initial_params == value);
  } else if (key == "set_zero_gauge") {
    if (value.size() == 1) {
      same = same & (set_zero_gauge == (std::stoi(value) == 1));
    } else {
      same = same & (set_zero_gauge == (value == "true"));
    }
  } else if (key == "allow_gap_couplings") {
    if (value.size() == 1) {
      same = same & (allow_gap_couplings == (std::stoi(value) == 1));
    } else {
      same = same & (allow_gap_couplings == (value == "true"));
    }
  } else if (key == "learn_rate_h") {
    same = same & (learn_rate_h == std::stod(value));
  } else if (key == "learn_rate_J") {
    same = same & (learn_rate_J == std::stod(value));
  } else if (key == "eta_min") {
    same = same & (eta_min == std::stod(value));
  } else if (key == "eta_max") {
    same = same & (eta_max == std::stod(value));
  } else if (key == "anneal_schedule") {
    same = same & (anneal_schedule == value);
  } else if (key == "anneal_scale") {
    same = same & (anneal_scale == std::stod(value));
  } else if (key == "anneal_period") {
    same = same & (anneal_period == std::stoi(value));
  } else if (key == "anneal_warm") {
    same = same & (anneal_warm == std::stoi(value));
  } else if (key == "anneal_hot") {
    same = same & (anneal_hot == std::stoi(value));
  } else if (key == "anneal_cool") {
    same = same & (anneal_cool == std::stoi(value));
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return same;
};

/**
 * @brief Check that all the necessary data exists to reload a given step.
 *
 * @param step iteration to check
 * @param output_binary flag for whether to look for binary data (.bin) or text
 * (.txt)
 *
 * @return (bool) flag for whether the necessary files were found
 */
bool
AdamW::isValidStep(int step, bool output_binary)
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

/**
 * @brief Initialize the model parameters.
 */
void
AdamW::initialize(void)
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
    const double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      double avg = 0;
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

/**
 * @brief Re-initialize the model.
 */
void
AdamW::reset()
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
    const double* freq_ptr = nullptr;
    for (int i = 0; i < N; i++) {
      double avg = 0;
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

  moment1.h.zeros();
  moment2.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      moment1.J(i, j).zeros();
      moment2.J(i, j).zeros();
    }
  }

  gradient.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J(i, j).zeros();
    }
  }
};

/**
 * @brief Re-load the model at a given step.
 *
 * @param step iteration to check
 * @param output_binary flag for whether to look for binary data (.bin) or text
 * (.txt)
 */
void
AdamW::restore(int step, bool output_binary)
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

/**
 * @brief Update the model parameters using the sampled sequence statistics.
 */
void
AdamW::update(void)
{
  params_prev = params;
  updateGradients();
  updateMoments();
  updateParameters();
  if (set_zero_gauge) {
    setZeroGauge();
  }
};

/**
 * @brief Update the gradients and compute the 1p/2p RMS error.
 */
void
AdamW::updateGradients(void)
{
  train_error_1p = 0;
  train_error_2p = 0;
  validation_error_1p = 0;
  validation_error_2p = 0;

  // Compute gradient
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      double delta =
        samples->frequency_1p(aa, i) - training->frequency_1p(aa, i);
      train_error_1p += pow(delta, 2);
      gradient.h(aa, i) = -delta;
    }
  }
  train_error_1p = sqrt(train_error_1p / (N * Q));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          double delta = -(training->frequency_2p(i, j)(aa1, aa2) -
                           samples->frequency_2p(i, j)(aa1, aa2));
          train_error_2p += pow(delta, 2);
          gradient.J(i, j)(aa1, aa2) = -delta;
        }
      }
    }
  }
  train_error_2p = sqrt(train_error_2p / (N * (N - 1.) * Q * Q / 2.));

  if (validate) {
    for (int i = 0; i < N; i++) {
      for (int aa = 0; aa < Q; aa++) {
        double delta =
          samples->frequency_1p(aa, i) - validation->frequency_1p(aa, i);
        validation_error_1p += pow(delta, 2);
      }
    }
    validation_error_1p = sqrt(validation_error_1p / (N * Q));

    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            double delta = -(validation->frequency_2p(i, j)(aa1, aa2) -
                             samples->frequency_2p(i, j)(aa1, aa2));
            validation_error_2p += pow(delta, 2);
          }
        }
      }
    }
    validation_error_2p =
      sqrt(validation_error_2p / (N * (N - 1.) * Q * Q / 2.));
  }
  return;
};

/**
 * @brief Update the first moment and second moment estimates.
 */
void
AdamW::updateMoments(void)
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

/**
 * @brief Update the model parameters using the updated gradients and moments.
 */
void
AdamW::updateParameters(void)
{
  double beta1 = 0.9;
  double beta2 = 0.999;
  double beta1_t = pow(beta1, step);
  double beta2_t = pow(beta2, step);

  double eta = 0;
  if (anneal_schedule == "none") {
    eta = 1.0;
  } else if (anneal_schedule == "cos") {
    if (step <= anneal_warm) {
      eta =
        eta_min + (eta_max - eta_min) * static_cast<double>(step) / anneal_warm;
    } else {
      int epoch_length = anneal_period;
      int total_length = anneal_warm;
      while ((step - total_length) > epoch_length) {
        total_length = total_length + epoch_length;
        epoch_length = epoch_length * anneal_scale;
      }
      eta = eta_min +
            0.5 * (eta_max - eta_min) *
              (1 + cos(3.1415926 * static_cast<double>(step - total_length) /
                       static_cast<double>(epoch_length)));
    }
  } else if (anneal_schedule == "trap") {
    if (step <= anneal_warm) {
      eta =
        eta_min + (eta_max - eta_min) * static_cast<double>(step) / anneal_warm;
    } else if (step <= (anneal_hot + anneal_warm)) {
      eta = eta_max;
    } else if (step <= (anneal_cool + anneal_hot + anneal_warm)) {
      eta = eta_max - (eta_max - eta_min) *
                        static_cast<double>(step - anneal_hot - anneal_warm) /
                        anneal_cool;
    } else {
      eta = eta_min;
    }
  }

  if (allow_gap_couplings) {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            params.J(i, j)(a, b) +=
              eta *
              (learn_rate_J *
                 (moment1.J(i, j)(a, b) / (1 - beta1_t) /
                  (sqrt(moment2.J(i, j)(a, b) / (1 - beta2_t)) + 0.00000001)) -
               lambda_decay_J * params.J(i, j)(a, b));
          }
        }
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 1; a < Q; a++) {
          for (int b = 1; b < Q; b++) {
            params.J(i, j)(a, b) +=
              eta *
              (learn_rate_J *
                 (moment1.J(i, j)(a, b) / (1 - beta1_t) /
                  (sqrt(moment2.J(i, j)(a, b) / (1 - beta2_t)) + 0.00000001)) -
               lambda_decay_J * params.J(i, j)(a, b));
          }
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      params.h(a, i) +=
        eta *
        (learn_rate_h * (moment1.h(a, i) / (1 - beta1_t) /
                         (sqrt(moment2.h(a, i) / (1 - beta2_t)) + 0.00000001)) -
         lambda_decay_h * params.h(a, i));
    }
  }
};

/**
 * @brief Write the current data to disk.
 *
 * @param str ID string for the output files.
 * @param output_binary flag for whether to write text or binary files.
 */
void
AdamW::writeData(std::string str, bool output_binary)
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

/**
 * @brief Delete the existing model data files for a given step.
 *
 * @param step iteration to check
 * @param output_binary flag for whether to look for binary data (.bin) or text
 * (.txt)
 *
 * This function is for clearing out steps with missing or incomplete data,
 * such as when the program is terminated in the middle of disk writes.
 */
void
AdamW::deleteStep(int step, bool output_binary)
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

/**
 * @brief Write the necessary data to be able to restart the inference.
 *
 * @param step current inference iteration
 * @param output_binary flag for whether to write text or binary files.
 *
 * This function is used to write all the data for restarting a given step,
 * which for Adam is the current and previous parameters, gradients, and 1st
 * and 2nd moments.
 */
void
AdamW::writeStep(int step, bool output_binary)
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
      "parameters_" + std::to_string(step - 1) + ".txt";
    writeParamsPreviousAscii(param_prev_file);

    std::string grad_file = "gradients_" + std::to_string(step) + ".txt";
    writeGradientAscii(grad_file);

    std::string moment1_file = "moment1_" + std::to_string(step) + ".txt";
    writeMoment1Ascii(moment1_file);

    std::string moment2_file = "moment2_" + std::to_string(step) + ".txt";
    writeMoment2Ascii(moment2_file);
  }
};

/**
 * @brief Write the parameters to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
AdamW::writeParams(std::string output_file_h, std::string output_file_J)
{
  params.h.save(output_file_h, arma::arma_binary);
  params.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the parameters to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
AdamW::writeParamsAscii(std::string output_file)
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

/**
 * @brief Write the previous step parameters to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
AdamW::writeParamsPrevious(std::string output_file_h, std::string output_file_J)
{
  params_prev.h.save(output_file_h, arma::arma_binary);
  params_prev.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the previous step parameters to disk in text format (1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
AdamW::writeParamsPreviousAscii(std::string output_file)
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

/**
 * @brief Write the 1st moment to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
AdamW::writeMoment1(std::string output_file_h, std::string output_file_J)
{
  moment1.h.save(output_file_h, arma::arma_binary);
  moment1.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the 2nd moment to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
AdamW::writeMoment2(std::string output_file_h, std::string output_file_J)
{
  moment2.h.save(output_file_h, arma::arma_binary);
  moment2.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the 1st moment to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
AdamW::writeMoment1Ascii(std::string output_file)
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

/**
 * @brief Write the 2nd moment to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
AdamW::writeMoment2Ascii(std::string output_file)
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

/**
 * @brief Write the gradients to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
AdamW::writeGradient(std::string output_file_h, std::string output_file_J)
{
  gradient.h.save(output_file_h, arma::arma_binary);
  gradient.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the gradients to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
AdamW::writeGradientAscii(std::string output_file)
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
