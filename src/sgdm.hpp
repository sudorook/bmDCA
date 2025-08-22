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

#ifndef SRC_SGDM_HPP_
#define SRC_SGDM_HPP_

#include "model.hpp"
#include "utils.hpp"

#include <string>

/**
 * @brief Stochastic gradient descent + momentum algorithm.
 *
 * Class for using stochastic gradient descent with momentum to learn
 * parameters. Moments are computed as an exponentially-decaying running
 * average of moment estimates from all steps during the run.
 *
 * Abstracts from the Model class and implements its virtual functions.
 */
class SGDM : public Model
{
public:
  SGDM(void);

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
  potts_model gradient; ///< model gradient
  potts_model moment1;  ///< 1st moment estimate

  double lambda_reg_h = 0.01; ///< regularization strength for fields
  double lambda_reg_J = 0.01; ///< regularization strength for couplings
  double alpha_reg = 1.0;     ///< elastic-net L1 / L2 scaling parameter
  std::string initial_params = "zero"; ///< initialization for parameters
  bool set_zero_gauge = false; ///< re-scale parameter for 0-mean J marginals
  bool allow_gap_couplings = true; ///< flag to regularize the gaps
  double learn_rate_h = 0.01;      ///< base learning rate for fields
  double learn_rate_J = 0.01;      ///< base learning rate for couplings
  double beta_h = 0.9;             ///< momentum average decay rate for fields
  double beta_J = 0.9; ///< momentum average decay rate for couplings

  void updateGradients(void);
  void updateMoments(void);
  void updateParameters(void);

  bool compareHyperparameters(const std::string) override;
  bool compareHyperparameter(const std::string, const std::string);
  void checkHyperparameters(void) override;
  void setHyperparameter(const std::string, const std::string);

  void writeParams(const std::string, const std::string);
  void writeParamsPrevious(const std::string, const std::string);
  void writeMoment1(const std::string, const std::string);
  void writeGradient(const std::string, const std::string);

  void writeParamsAscii(const std::string);
  void writeParamsPreviousAscii(const std::string);
  void writeMoment1Ascii(const std::string);
  void writeGradientAscii(const std::string);
};

#endif // SRC_SGDM_HPP_
