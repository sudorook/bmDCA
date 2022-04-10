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

#ifndef BMDCA_RUN_HPP
#define BMDCA_RUN_HPP

#include "model.hpp"
#include "msa.hpp"
#include "msa_stats.hpp"
#include "pcg_random.hpp"
#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"

/**
 * @brief Class for performing Boltzmann-machine inference.
 */
class Sim
{
public:
  Sim(MSA*, MSA*, std::string, std::string, bool);
  ~Sim(void);
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
  std::string train_mode = "original"; ///< training model
  bool cross_validate = false;         ///< flag to cross-validate the MSA
  int validation_seqs = 100;           ///< effective M for validatino set
  unsigned random_seed = 142857;       ///< random seed for RNG
  bool output_binary = true; ///< flag to output in arma::binary format

  // Sampler settings
  std::string update_rule = "mh"; ///< update proposal rule for sampler
  int burn_in_start = 10000;      ///< initial burn-in
  int burn_between_start = 100;   ///< initial burn-between
  bool update_burn_time = true;   ///< flag to tune burn times
  double adapt_up_time = 1.2;     ///< scaling factor for increasing burn times
  double adapt_down_time = 0.9;   ///< scaling factorfor decreasing burn times
  int step_importance_max = 1;    ///< max number of importance sampling steps
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
  MSA* msa_train = nullptr;               ///< address of training MSA
  MSA* msa_validate = nullptr;            ///< address of validation MSA
  MSAStats* msa_train_stats = nullptr;    ///< address of training MSA stats
  MSAStats* msa_validate_stats = nullptr; ///< address of validation MSA stats

  arma::Col<double> msa_train_energies;    ///< energies for training MSA
  arma::Col<double> msa_validate_energies; ///< energies for validation MSA
  void computeMSAEnergies(arma::Col<double>*, MSA*, potts_model*);
  void writeMSAEnergies(int);
  void writeMSAEnergies(std::string);

  Model* model; ///< training model

  Sampler* sampler; ///< Potts model sampler

  SampleStats* sample_stats; ///< address of sample sequence statistics

  pcg32 rng; ///< PCG32 random number generator
  arma::Mat<unsigned long>
    rng_buffer; ///< stores rng values to write to run log
};

#endif
