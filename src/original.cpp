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

#include "original.hpp"

#include <armadillo>
#include <string>

#include "utils.hpp"

/**
 * @brief Original constructor.
 */
Original::Original()
  : Model(){};

/**
 * @brief Load model hyperparameters from config file.
 *
 * @param file_name config file string
 */
void
Original::loadHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  if (file.is_open()) {
    std::string line;
    bool reading_original_section = false;
    while (std::getline(file, line)) {
      if (line.empty()) {
        reading_original_section = false;
        continue;
      } else if (line[0] == '#') {
        reading_original_section = false;
        continue;
      } else if (line[0] == '[') {
        if (line == "[[original]]") {
          reading_original_section = true;
          continue;
        } else {
          reading_original_section = false;
          continue;
        }
      }
      if (reading_original_section) {
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
Original::compareHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool all_same = true;
  if (file.is_open()) {
    bool reading_original_section = false;
    std::string line;
    while (std::getline(file, line)) {
      if (line.empty()) {
        continue;
      } else if (line[0] == '#') {
        continue;
      } else if (line[0] == '[') {
        if (line == "[[original]]") {
          reading_original_section = true;
          continue;
        } else {
          reading_original_section = false;
          continue;
        }
      }
      if (reading_original_section) {
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
Original::setHyperparameter(std::string key, std::string value)
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
  } else if (key == "allow_gap_couplings") {
    if (value.size() == 1) {
      allow_gap_couplings = (std::stoi(value) == 1);
    } else {
      allow_gap_couplings = (value == "true");
    }
  } else if (key == "epsilon_h") {
    epsilon_h = std::stod(value);
  } else if (key == "epsilon_J") {
    epsilon_J = std::stod(value);
  } else if (key == "learn_rate_h_min") {
    learn_rate_h_min = std::stod(value);
  } else if (key == "learn_rate_J_min") {
    learn_rate_J_min = std::stod(value);
  } else if (key == "learn_rate_h_max") {
    learn_rate_h_max = std::stod(value);
  } else if (key == "learn_rate_J_max") {
    learn_rate_J_max = std::stod(value);
  } else if (key == "adapt_up") {
    adapt_up = std::stod(value);
  } else if (key == "adapt_down") {
    adapt_down = std::stod(value);
  } else if (key == "use_pos_reg") {
    if (value.size() == 1) {
      use_pos_reg = (std::stoi(value) == 1);
    } else {
      use_pos_reg = (value == "true");
    }
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
Original::checkHyperparameters(void) {};

/**
 * @brief Write stored hyperparameters to file
 *
 * @param output_file config file string
 * @param append flag for whether or not to append the hyperparameters to the
 * file or overwrite it
 */
void
Original::writeHyperparameters(std::string output_file, bool append)
{
  std::ofstream stream;
  if (append) {
    stream.open(output_file, std::ofstream::out | std::ofstream::app);
  } else {
    stream.open(output_file, std::ofstream::out | std::ofstream::trunc);
  }
  // Header
  stream << std::endl;
  stream << "[[original]]" << std::endl;

  // BM settings
  stream << "lambda_reg_h=" << lambda_reg_h << std::endl;
  stream << "lambda_reg_J=" << lambda_reg_J << std::endl;
  stream << "alpha_reg=" << alpha_reg << std::endl;
  stream << "initial_params=" << initial_params << std::endl;
  stream << "set_zero_gauge=" << set_zero_gauge << std::endl;
  stream << "allow_gap_couplings=" << allow_gap_couplings << std::endl;
  stream << "epsilon_h=" << epsilon_h << std::endl;
  stream << "epsilon_J=" << epsilon_J << std::endl;
  stream << "learn_rate_h_min=" << learn_rate_h_min << std::endl;
  stream << "learn_rate_J_min=" << learn_rate_J_min << std::endl;
  stream << "learn_rate_h_max=" << learn_rate_h_max << std::endl;
  stream << "learn_rate_J_max=" << learn_rate_J_max << std::endl;
  stream << "adapt_up=" << adapt_up << std::endl;
  stream << "adapt_down=" << adapt_down << std::endl;
  stream << "use_pos_reg=" << use_pos_reg << std::endl;

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
Original::compareHyperparameter(std::string key, std::string value)
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
  } else if (key == "allow_gap_couplings") {
    if (value.size() == 1) {
      same = same & (allow_gap_couplings == (std::stoi(value) == 1));
    } else {
      same = same & (allow_gap_couplings == (value == "true"));
    }
  } else if (key == "epsilon_h") {
    same = same & (epsilon_h == std::stod(value));
  } else if (key == "epsilon_J") {
    same = same & (epsilon_J == std::stod(value));
  } else if (key == "learn_rate_h_min") {
    same = same & (learn_rate_h_min == std::stod(value));
  } else if (key == "learn_rate_J_min") {
    same = same & (learn_rate_J_min == std::stod(value));
  } else if (key == "learn_rate_h_max") {
    same = same & (learn_rate_h_max == std::stod(value));
  } else if (key == "learn_rate_J_max") {
    same = same & (learn_rate_J_max == std::stod(value));
  } else if (key == "adapt_up") {
    same = same & (adapt_up == std::stod(value));
  } else if (key == "adapt_down") {
    same = same & (adapt_down == std::stod(value));
  } else if (key == "use_pos_reg") {
    if (value.size() == 1) {
      same = same & (use_pos_reg == (std::stoi(value) == 1));
    } else {
      same = same & (use_pos_reg == (value == "true"));
    }
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
Original::isValidStep(int step, bool output_binary)
{
  bool valid = false;
  if (output_binary) {
    if (checkFileExists("parameters_h_" + std::to_string(step) + ".bin") &
        checkFileExists("parameters_J_" + std::to_string(step) + ".bin") &
        checkFileExists("parameters_h_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("parameters_J_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("gradients_h_" + std::to_string(step) + ".bin") &
        checkFileExists("gradients_J_" + std::to_string(step) + ".bin") &
        checkFileExists("gradients_h_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("gradients_J_" + std::to_string(step - 1) + ".bin") &
        checkFileExists("learning_rates_h_" + std::to_string(step) + ".bin") &
        checkFileExists("learning_rates_J_" + std::to_string(step) + ".bin")) {
      valid = true;
    }
  } else {
    if (checkFileExists("parameters_" + std::to_string(step) + ".txt") &
        checkFileExists("parameters_" + std::to_string(step - 1) + ".txt") &
        checkFileExists("gradients_" + std::to_string(step) + ".txt") &
        checkFileExists("gradients_" + std::to_string(step - 1) + ".txt") &
        checkFileExists("learning_rates_" + std::to_string(step) + ".txt")) {
      valid = true;
    }
  }
  return valid;
};

/**
 * @brief Initialize the model parameters.
 */
void
Original::initialize(void)
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

  // Initialize the learning rates
  learning_rates.h = arma::Mat<double>(Q, N);
  learning_rates.h.fill(epsilon_h);

  learning_rates.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      learning_rates.J(i, j) = arma::Mat<double>(Q, Q);
      learning_rates.J(i, j).fill(epsilon_J);
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
Original::reset()
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
  learning_rates.h.zeros();
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      learning_rates.J(i, j).zeros();
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
Original::restore(int step, bool output_binary)
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

    std::string grad_prev_h_file =
      "gradients_h_" + std::to_string(step - 1) + ".bin";
    std::string grad_prev_J_file =
      "gradients_J_" + std::to_string(step - 1) + ".bin";
    gradient_prev = loadPottsModel(grad_prev_h_file, grad_prev_J_file);

    std::string learning_rates_h_file =
      "learning_rates_h_" + std::to_string(step) + ".bin";
    std::string learning_rates_J_file =
      "learning_rates_J_" + std::to_string(step) + ".bin";
    learning_rates =
      loadPottsModel(learning_rates_h_file, learning_rates_J_file);
  } else {
    std::string param_file = "parameters_" + std::to_string(step) + ".txt";
    params = loadPottsModelAscii(param_file);

    std::string param_prev_file =
      "parameters_" + std::to_string(step - 1) + ".txt";
    params_prev = loadPottsModelAscii(param_prev_file);

    std::string grad_file = "gradients_" + std::to_string(step) + ".txt";
    gradient = loadPottsModelAscii(grad_file);

    std::string grad_prev_file =
      "gradients_" + std::to_string(step - 1) + ".txt";
    gradient_prev = loadPottsModelAscii(grad_prev_file);

    std::string learning_rates_file =
      "learning_rates_" + std::to_string(step) + ".txt";
    learning_rates = loadPottsModelAscii(learning_rates_file);
  }
};

/**
 * @brief Update the model parameters using the sampled sequence statistics.
 */
void
Original::update(void)
{
  params_prev = params;
  gradient_prev = gradient;
  updateGradients();
  updateLearningRates();
  updateParameters();
  if (set_zero_gauge) {
    setZeroGauge();
  }
};

/**
 * @brief Update the gradients and compute the 1p/2p RMS error.
 */
void
Original::updateGradients(void)
{
  train_error_1p = 0;
  train_error_2p = 0;
  validation_error_1p = 0;
  validation_error_2p = 0;

  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      double delta =
        samples->frequency_1p(aa, i) - training->frequency_1p(aa, i) +
        lambda_reg_h *
          (alpha_reg * params.h(aa, i) +
           (1. - alpha_reg) *
             (0.5 - static_cast<double>(std::signbit(params.h(aa, i)))));
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
                 (1. - alpha_reg) * (0.5 - static_cast<double>(std::signbit(
                                             params.J(i, j)(aa1, aa2))))));
          if (use_pos_reg) {
            delta = delta * 1. /
                    (1. + fabs(training->rel_entropy_grad_1p(aa1, i))) /
                    (1. + fabs(training->rel_entropy_grad_1p(aa2, j)));
          }
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
          samples->frequency_1p(aa, i) - validation->frequency_1p(aa, i) +
          lambda_reg_h *
            (alpha_reg * params.h(aa, i) +
             (1. - alpha_reg) *
               (0.5 - static_cast<double>(std::signbit(params.h(aa, i)))));
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
                   (1. - alpha_reg) * (0.5 - static_cast<double>(std::signbit(
                                               params.J(i, j)(aa1, aa2))))));
            if (use_pos_reg) {
              delta = delta * 1. /
                      (1. + fabs(validation->rel_entropy_grad_1p(aa1, i))) /
                      (1. + fabs(validation->rel_entropy_grad_1p(aa2, j)));
            }
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
Original::updateLearningRates(void)
{
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          double alfa = gradient.J(i, j)(a, b) * gradient_prev.J(i, j)(a, b);
          alfa = Theta<double, double>(alfa) * adapt_up +
                 Theta<double, double>(-alfa) * adapt_down +
                 Delta<double, double>(alfa);

          learning_rates.J(i, j)(a, b) = Min<double>(
            learn_rate_J_max,
            Max<double>(learn_rate_J_min, alfa * learning_rates.J(i, j)(a, b)));
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      double alfa = gradient.h(a, i) * gradient_prev.h(a, i);
      alfa = Theta<double, double>(alfa) * adapt_up +
             Theta<double, double>(-alfa) * adapt_down +
             Delta<double, double>(alfa);
      learning_rates.h(a, i) = Min<double>(
        learn_rate_h_max,
        Max<double>(learn_rate_h_min, alfa * learning_rates.h(a, i)));
    }
  }
};

/**
 * @brief Update the model parameters using the updated gradients and moments.
 */
void
Original::updateParameters(void)
{
  if (allow_gap_couplings) {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int a = 0; a < Q; a++) {
          for (int b = 0; b < Q; b++) {
            params.J(i, j)(a, b) +=
              learning_rates.J(i, j)(a, b) * gradient.J(i, j)(a, b);
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
              learning_rates.J(i, j)(a, b) * gradient.J(i, j)(a, b);
          }
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      params.h(a, i) += learning_rates.h(a, i) * gradient.h(a, i);
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
Original::writeData(std::string str, bool output_binary)
{
  if (output_binary) {
    std::string param_h_file = "parameters_h_" + str + ".bin";
    std::string param_J_file = "parameters_J_" + str + ".bin";
    writeParams(param_h_file, param_J_file);

    std::string grad_h_file = "gradients_h_" + str + ".bin";
    std::string grad_J_file = "gradients_J_" + str + ".bin";
    writeGradient(grad_h_file, grad_J_file);

    std::string learning_rates_h_file = "learning_rates_h_" + str + ".bin";
    std::string learning_rates_J_file = "learning_rates_J_" + str + ".bin";
    writeLearningRates(learning_rates_h_file, learning_rates_J_file);
  } else {
    std::string param_file = "parameters_" + str + ".txt";
    writeParamsAscii(param_file);

    std::string grad_file = "gradients_" + str + ".txt";
    writeGradientAscii(grad_file);

    std::string learning_rates_file = "learning_rates_" + str + ".txt";
    writeLearningRatesAscii(learning_rates_file);
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
Original::deleteStep(int step, bool output_binary)
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

    file = "learning_rates_h_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "learning_rates_J_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

  } else {
    file = "parameters_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "gradients_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "learning_rates_" + std::to_string(step) + ".txt";
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
Original::writeStep(int step, bool output_binary)
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

    std::string gradient_prev_h_file =
      "gradients_h_" + std::to_string(step - 1) + ".bin";
    std::string gradient_prev_J_file =
      "gradients_J_" + std::to_string(step - 1) + ".bin";
    writeParamsPrevious(gradient_prev_h_file, gradient_prev_J_file);

    std::string learning_rates_h_file =
      "learning_rates_h_" + std::to_string(step) + ".bin";
    std::string learning_rates_J_file =
      "learning_rates_J_" + std::to_string(step) + ".bin";
    writeLearningRates(learning_rates_h_file, learning_rates_J_file);
  } else {
    std::string param_file = "parameters_" + std::to_string(step) + ".txt";
    writeParamsAscii(param_file);

    std::string param_prev_file =
      "parameters_" + std::to_string(step - 1) + ".txt";
    writeParamsPreviousAscii(param_prev_file);

    std::string grad_file = "gradients_" + std::to_string(step) + ".txt";
    writeGradientAscii(grad_file);

    std::string grad_prev_file =
      "gradients_" + std::to_string(step - 1) + ".txt";
    writeGradientPreviousAscii(grad_prev_file);

    std::string learning_rates_file =
      "learning_rates_" + std::to_string(step) + ".txt";
    writeLearningRatesAscii(learning_rates_file);
  }
};

/**
 * @brief Write the parameters to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
Original::writeParams(std::string output_file_h, std::string output_file_J)
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
Original::writeParamsAscii(std::string output_file)
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
Original::writeParamsPrevious(std::string output_file_h,
                              std::string output_file_J)
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
Original::writeParamsPreviousAscii(std::string output_file)
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
 * @brief Write the gradients to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
Original::writeGradient(std::string output_file_h, std::string output_file_J)
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
Original::writeGradientAscii(std::string output_file)
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

/**
 * @brief Write the previous gradients to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
Original::writeGradientPrevious(std::string output_file_h,
                                std::string output_file_J)
{
  gradient_prev.h.save(output_file_h, arma::arma_binary);
  gradient_prev.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the previous gradients to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
Original::writeGradientPreviousAscii(std::string output_file)
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

/**
 * @brief Write the learning rates to disk in arma::binary format.
 *
 * @param output_file_h file string for fields
 * @param output_file_J file string for couplings
 */
void
Original::writeLearningRates(std::string output_file_h,
                             std::string output_file_J)
{
  learning_rates.h.save(output_file_h, arma::arma_binary);
  learning_rates.J.save(output_file_J, arma::arma_binary);
};

/**
 * @brief Write the learning rates to disk in text format (all 1 file).
 *
 * @param output_file file string for fields and couplings
 */
void
Original::writeLearningRatesAscii(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = learning_rates.h.n_cols;
  int Q = learning_rates.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << learning_rates.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << learning_rates.h(aa, i)
                    << std::endl;
    }
  }
};
