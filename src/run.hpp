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

#ifndef SRC_RUN_HPP_
#define SRC_RUN_HPP_

#include "model.hpp"
#include "msa.hpp"
#include "msa_stats.hpp"
#include "pcg_random.hpp"
#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"

#include <cstdint>
#include <memory>
#include <string>

/**
 * @brief Class for performing Boltzmann-machine inference.
 */
class Sim
{
public:
  Sim(std::shared_ptr<MSA>,
      std::shared_ptr<MSA>,
      std::string,
      std::string,
      bool);
  void run(void);

private:
  // Functions for processing hyperparameter files.
  void loadHyperparameters(std::string);
  void writeHyperparameters(std::string, bool = false);
  void setHyperparameter(std::string, std::string);
  bool compareHyperparameters(std::string);
  bool compareHyperparameter(std::string, std::string);
  void checkHyperparameters(void);

  void initialize(void);

  // Functions for reloading run states
  void setStepOffset(void);
  bool isValidStep(int);
  void deleteStep(int);
  void restoreRunState(void);

  void writeData(std::string);
  void writeStep(int);
  // void clearFiles(std::string);

  // BM settings
  int step_max = 500;                  ///< maximum number of iterations
  int save_period = 40;                ///< period for saving parameters to disk
  bool save_best_steps = false;        ///< flag to save steps with min rmse
  std::string stop_mode = "threshold"; ///< stopping criterion
  double stop_threshold = 0.00001;     ///< threshold for stopping
  std::string train_mode = "sgdm";     ///< training model
  bool cross_validate = true;          ///< flag to cross-validate the MSA
  int validation_seqs = 100;           ///< effective M for validatino set
  unsigned random_seed = 142857;       ///< random seed for RNG
  bool output_binary = true; ///< flag to output in arma::binary format

  // Sampler settings
  std::string update_rule = "mh"; ///< update proposal rule for sampler
  int burn_in_start = 1000;       ///< initial burn-in
  int burn_between_start = 100;   ///< initial burn-between
  bool update_burn_time = true;   ///< flag to tune burn times
  double adapt_up_time = 1.2;     ///< scaling factor for increasing burn times
  double adapt_down_time = 0.9;   ///< scaling factorfor decreasing burn times
  int step_importance_max = 0;    ///< max number of importance sampling steps
  double coherence_min =
    0.9999;            ///< model coherence to maintain for importance sampling
  bool use_ss = false; ///< flag to sample effective M sequences
  int walkers = 1000;  ///< number of trajectories
  int samples_per_walk = 1; ///< samples per trajectory

  // Run settings
  int step = 1; ///< current step
  int step_offset =
    0; ///< placeholder for restoring step number when restarting
  int burn_in = burn_in_start;           ///< burn-in time
  int burn_between = burn_between_start; ///< burn-between time
  double train_err_tot_min = 1000;       ///< min RMSE for training MSA
  double validate_err_tot_min = 1000;    ///< min RMSE validation MSA
  double diff_avg_energy =
    0; ///< difference in training and validation energies

  std::string hyperparameter_file =
    "bmdca_params.conf"; ///< file string for runtime config file
  std::string run_log_file = "bmdca_run.log"; ///< file string for run log

  bool checkErgodicity(void);

  // Buffers
  arma::Mat<double> run_buffer; ///< stores values to write to run log
  void initializeRunLog();
  void writeRunLog(int = -1, int = 0, bool = false);

  // Sample data
  arma::Cube<int> samples_3d; ///< stores samples when samples_per_walk > 1
  arma::Mat<int> samples_2d;  ///< stores samples when samples_per_walk == 1

  // Stats from original MSA
  std::shared_ptr<MSA> msa_train;            ///< address of training MSA
  std::shared_ptr<MSA> msa_validate;         ///< address of validation MSA
  std::shared_ptr<MSAStats> msa_train_stats; ///< address of training MSA stats
  std::shared_ptr<MSAStats>
    msa_validate_stats; ///< address of validation MSA stats

  arma::Col<double> msa_train_energies;    ///< energies for training MSA
  arma::Col<double> msa_validate_energies; ///< energies for validation MSA
  void computeMSAEnergies(arma::Col<double>*,
                          std::shared_ptr<MSA>,
                          std::shared_ptr<potts_model>);
  void writeMSAEnergies(int);
  void writeMSAEnergies(std::string);

  // TODO: should these be std::unique_ptr?
  std::shared_ptr<Model> model; ///< training model

  std::shared_ptr<Sampler> sampler; ///< Potts model sampler

  std::shared_ptr<SampleStats>
    sample_stats; ///< address of sample sequence statistics

  pcg32 rng;                      ///< PCG32 random number generator
  arma::Mat<uint64_t> rng_buffer; ///< stores rng values to write to run log
};

#endif // SRC_RUN_HPP_
