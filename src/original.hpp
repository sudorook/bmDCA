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

#ifndef SRC_ORIGINAL_HPP_
#define SRC_ORIGINAL_HPP_

#include "model.hpp"
#include "utils.hpp"

#include <string>

/**
 * @brief 'Original' gradient descent algorithm.
 *
 * Class for using, for lack of a better name, the original technique for
 * learning parameters. Each parameter has its own learning rate, which is
 * scaled up or down at runtime if the sign of the gradient did not or did
 * change between the current and previous step, respectively.
 *
 * Abstracts from the Model class and implements its virtual functions.
 */
class Original : public Model
{
public:
  Original(void);

  void update(void) override;
  void initialize(void) override;
  void reset(void) override;
  void restore(int, bool = true) override;

  void writeData(std::string, bool = true) override;
  void writeStep(int, bool = true) override;
  void deleteStep(int, bool = true) override;

  void loadHyperparameters(std::string) override;
  void writeHyperparameters(std::string, bool = true) override;
  bool isValidStep(int, bool = true) override;

private:
  potts_model gradient;       ///< model gradient
  potts_model gradient_prev;  ///< previous step model gradient
  potts_model learning_rates; ///< learning rates for each parameter

  double lambda_reg_h = 0.01; ///< regularization strength for fields
  double lambda_reg_J = 0.01; ///< regularization strength for couplings
  double alpha_reg = 1.0;     ///< elastic-net L1 / L2 scaling parameter
  std::string initial_params = "zero"; ///< initialization for parameters
  bool set_zero_gauge = false; ///< re-scale parameter for 0-mean J marginals
  bool allow_gap_couplings = true; ///< flag to regularize the gaps
  double epsilon_h = 0.01;         ///< initial learning rate for fields
  double epsilon_J = 0.001;        ///< initial learning rate for couplings
  double learn_rate_h_min = 1e-04; ///< minimum learning rate for fields
  double learn_rate_h_max = 0.5;   ///< maximum learning rate for fields
  double learn_rate_J_min = 1e-04; ///< minimum learning rate for couplings
  double learn_rate_J_max = 0.5;   ///< maximum learning rate for couplings
  double adapt_up = 1.5;    ///< scaling factor for increasing learning rate
  double adapt_down = 0.6;  ///< scaling factor for decreasing learning rate
  bool use_pos_reg = false; ///< flag for position-specific regularization

  void updateGradients(void);
  void updateLearningRates(void);
  void updateParameters(void);

  bool compareHyperparameters(std::string) override;
  bool compareHyperparameter(std::string, std::string);
  void checkHyperparameters(void) override;
  void setHyperparameter(std::string, std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeGradient(std::string, std::string);
  void writeGradientPrevious(std::string, std::string);
  void writeLearningRates(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeGradientAscii(std::string);
  void writeGradientPreviousAscii(std::string);
  void writeLearningRatesAscii(std::string);
};

#endif // SRC_ORIGINAL_HPP_
