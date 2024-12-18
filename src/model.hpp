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

#ifndef SRC_MODEL_HPP_
#define SRC_MODEL_HPP_

#include "msa_stats.hpp"
#include "sample_stats.hpp"
#include "utils.hpp"
#include <string>

/**
 * @brief Abstract Model class.
 *
 * The inference scheme for learning parameters is not known at compile time.
 * Instead, it's set via the config file at runtime. This Model class defines a
 * set of virtual functions to be shared amongst all classes derived from this
 * one. A common function table allows an instance of any derived class
 * specified in the config file, recast to a Model, to be loaded at runtime.
 */
class Model
{
public:
  Model();
  virtual ~Model(){};

  void setMSAStats(std::shared_ptr<MSAStats>, std::shared_ptr<MSAStats>);
  void setSampleStats(std::shared_ptr<SampleStats>);
  void setStep(int);

  double train_error_1p = 1000;      ///< training RMSE for fields
  double train_error_2p = 1000;      ///< training RMSE for couplings
  double validation_error_1p = 1000; ///< validation RMSE for fields
  double validation_error_2p = 1000; ///< validation RMSE for couplings

  // Define virtual functions that form the Model 'API'
  virtual void update(void) = 0;       ///< update parameters, gradients, etc.
  virtual void initialize(void) = 0;   ///< initialize the model
  virtual void reset(void) = 0;        ///< reset the model to initial values
  virtual void restore(int, bool) = 0; ///< re-load parameters from disk

  virtual void writeData(std::string,
                         bool = true) = 0; ///< write current model to disk
  virtual void writeStep(
    int,
    bool = true) = 0; ///< write all data to reload current step
  virtual void deleteStep(int, bool = true) = 0; ///< delete data for a step
  virtual bool isValidStep(
    int,
    bool = true) = 0; ///< check if data exists to reload a step

  virtual void loadHyperparameters(
    std::string) = 0; ///< read parameters from config file
  virtual void checkHyperparameters(
    void) = 0; ///< check config parameters against saved state
  virtual bool compareHyperparameters(
    std::string) = 0; ///< check that config parameters are valid
  virtual void writeHyperparameters(std::string,
                                    bool = false) = 0; ///< write config file

  potts_model params;      ///< current step parameters
  potts_model params_prev; ///< previous step parameters

protected:
  std::shared_ptr<MSAStats> training;   ///< pointer to training MSA stats
  std::shared_ptr<MSAStats> validation; ///< pointer to validation MSA stats
  std::shared_ptr<SampleStats> samples; ///< pointer to sampled sequence stats

  std::string hyperparameter_file; ///< file string for storing hyperparams

  bool validate; ///< flag for whether to compute validation error

  void setZeroGauge(void);

  int N; ///< number of positions
  int Q; ///< number of states

  double temperature = 1.0; ///< inference temperature (unused)

  int step = 1; ///< iteration number
};

#endif // SRC_MODEL_HPP_
