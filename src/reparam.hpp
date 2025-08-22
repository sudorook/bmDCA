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

#ifndef SRC_REPARAM_HPP_
#define SRC_REPARAM_HPP_

#include "model.hpp"
#include "utils.hpp"

#include <string>

/**
 * @brief 'Original' gradient descent algorithm using a reparameterized model.
 *
 * Class for using a reparameterized Potts model for updating the parameters at
 * each iteration, using the 'original' technique for learning parameters. Each
 * parameter has its own learning rate, which is scaled up or down at runtime
 * if the sign of the gradient did not or did change between the current and
 * previous step, respectively. Regularizations for fields and couplings are
 * not independent.
 *
 * Abstracts from the Model class and implements its virtual functions.
 */
class Reparam : public Model
{
public:
  Reparam(void);

  void update(void) override;
  void initialize(void) override;
  void reset(void) override;
  void restore(int, bool = true) override;

  void writeData(const std::string, bool = true) override;
  void writeStep(int, bool = true) override;
  void deleteStep(int, bool = true) override;

  void loadHyperparameters(const std::string) override;
  void writeHyperparameters(const std::string, bool = true) override;
  bool isValidStep(int, bool = true) override;

private:
  potts_model gradient;       ///< model gradient
  potts_model gradient_prev;  ///< previous step model gradient
  potts_model learning_rates; ///< learning rates for each parameter

  double lambda_reg_phi_h = 0.01; ///< regularization strength for fields
  double lambda_reg_phi_J = 0.01; ///< regularization strength for couplings
  std::string initial_params = "profile"; ///< initialization for parameters
  bool set_zero_gauge = false; ///< re-scale parameter for 0-mean J marginals
  double epsilon_h = 0.01;     ///< initial learning rate for fields
  double epsilon_J = 0.001;    ///< initial learning rate for couplings
  double learn_rate_h_min = 1e-04; ///< minimum learning rate for fields
  double learn_rate_J_min = 1e-04; ///< maximum learning rate for fields
  double learn_rate_h_max = 0.5;   ///< minimum learning rate for couplings
  double learn_rate_J_max = 0.5;   ///< maximum learning rate for couplings
  double adapt_up = 1.5;   ///< scaling factor for increasing learning rate
  double adapt_down = 0.6; ///< scaling factor for decreasing learning rate

  void updateGradients(void);
  void updateLearningRates(void);
  void updateParameters(void);

  bool compareHyperparameters(const std::string) override;
  bool compareHyperparameter(const std::string, const std::string);
  void checkHyperparameters(void) override;
  void setHyperparameter(const std::string, const std::string);

  void writeParams(const std::string, const std::string);
  void writeParamsPrevious(const std::string, const std::string);
  void writeGradient(const std::string, const std::string);
  void writeGradientPrevious(const std::string, const std::string);
  void writeLearningRates(const std::string, const std::string);

  void writeParamsAscii(const std::string);
  void writeParamsPreviousAscii(const std::string);
  void writeGradientAscii(const std::string);
  void writeGradientPreviousAscii(const std::string);
  void writeLearningRatesAscii(const std::string);
};

#endif // SRC_REPARAM_HPP_
