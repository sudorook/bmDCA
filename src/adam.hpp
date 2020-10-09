/* Boltzmann-machine Direct Coupling Analysis (bmDCA)
 * Copyright (C) 2020
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ADAM_HPP
#define ADAM_HPP

#include "model.hpp"
#include "utils.hpp"

/**
 * @brief Adam gradient descent algorithm.
 *
 * Class for using the Adam (_Ada_ptive _m_oment) stochastic gradient descent
 * algorithm for learning parameters. Adam uses a momentum-based update rule
 * that scales each estimate of the first moment by an estimate of the second
 * moment.
 *
 * Abstracts from the Model class and implements its virtual functions.
 */
class Adam : public Model
{
public:
  Adam(void);

  void update(void);
  void initialize(void);
  void reset(void);
  void restore(int, bool = true);

  void writeData(std::string, bool = true);
  void writeStep(int, bool = true);
  void deleteStep(int, bool = true);

  void loadHyperparameters(std::string);
  void writeHyperparameters(std::string, bool = true);
  bool isValidStep(int, bool = true);

private:
  potts_model gradient; ///< model gradient
  potts_model moment1;  ///< 1st moment estimate
  potts_model moment2;  ///< 2nd moment estimate

  double lambda_reg_h = 0.01;          ///< regularization strength for fields
  double lambda_reg_J = 0.01;          ///< regularization strength for couplings
  double alpha_reg = 1.0;              ///< elastic-net L1 / L2 scaling parameter
  std::string initial_params = "zero"; ///< initialization for parameters
  bool set_zero_gauge = false;         ///< re-scale parameter for 0-mean J marginals
  bool allow_gap_couplings = true;     ///< flag to regularize the gaps
  double learn_rate_h = 0.01;          ///< base learning rate for fields
  double learn_rate_J = 0.01;          ///< base learning rate for couplings

  void updateGradients(void);
  void updateMoments(void);
  void updateParameters(void);

  bool compareHyperparameters(std::string);
  bool compareHyperparameter(std::string, std::string);
  void checkHyperparameters(void);
  void setHyperparameter(std::string, std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeMoment1(std::string, std::string);
  void writeMoment2(std::string, std::string);
  void writeGradient(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeMoment1Ascii(std::string);
  void writeMoment2Ascii(std::string);
  void writeGradientAscii(std::string);
};

#endif
