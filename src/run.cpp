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

#include "run.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
// #include <random>
#include <regex>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include "adam.hpp"
#include "adamw.hpp"
#include "original.hpp"
#include "radam.hpp"
#include "reparam.hpp"
#include "sgdm.hpp"

#define EPSILON 0.00000001
#define PI 3.1415926

/**
 * @brief Constructor for inferencce loop.
 *
 * @param msa_train pointer to training MSA
 * @param msa_validate pointer to validation MSA (optional)
 * @param config_file file string for hyperparameter file
 * @param dest_dir directory to output bmDCA results
 * @param force_restart force the inference to start from the beginning
 */
Sim::Sim(std::shared_ptr<MSA> msa_train,
         std::shared_ptr<MSA> msa_validate,
         std::string config_file,
         std::string dest_dir,
         bool force_restart)
  : msa_train(msa_train)
  , msa_validate(msa_validate)
{
  std::cout << "--------------------------- bmDCA ---------------------------"
            << std::endl;

  // Load bmDCA hyperparameters.
  if (!config_file.empty()) {
    std::cout << "loading bmDCA hyperparameters... " << std::flush;
    loadHyperparameters(config_file);
    std::cout << "done." << std::endl;
  } else if ((!force_restart) &
             (checkFileExists(dest_dir + "/" + hyperparameter_file))) {
    std::cout << "loading previous bmDCA hyperparameters... " << std::flush;
    loadHyperparameters(dest_dir + "/" + hyperparameter_file);
    std::cout << "done." << std::endl;
  }
  checkHyperparameters();

  // Load the training model.
  if (train_mode == "adam") {
    model = std::make_shared<Adam>();
  } else if (train_mode == "adamw") {
    model = std::make_shared<AdamW>();
  } else if (train_mode == "original") {
    model = std::make_shared<Original>();
  } else if (train_mode == "radam") {
    model = std::make_shared<RAdam>();
  } else if (train_mode == "reparametrization") {
    model = std::make_shared<Reparam>();
  } else if (train_mode == "sgdm") {
    model = std::make_shared<SGDM>();
  } else {
    std::cerr << "ERROR: unrecognised training model '" << train_mode
              << "' given. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Load model hyperparameters.
  if (!config_file.empty()) {
    std::cout << "loading model hyperparameters... " << std::flush;
    model->loadHyperparameters(config_file);
    std::cout << "done." << std::endl;
  } else if ((!force_restart) &
             (checkFileExists(dest_dir + "/" + hyperparameter_file))) {
    std::cout << "loading previous model hyperparameters... " << std::flush;
    model->loadHyperparameters(dest_dir + "/" + hyperparameter_file);
    std::cout << "done." << std::endl;
  }
  model->checkHyperparameters();

  if (!dest_dir.empty()) {
    chdir(dest_dir.c_str());
  }

  if ((!force_restart) & (checkFileExists(hyperparameter_file))) {
    if ((!compareHyperparameters(hyperparameter_file)) |
        (!model->compareHyperparameters(hyperparameter_file))) {
      std::cerr << "ERROR: current and previous hyperparameters mismatched."
                << std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      setStepOffset();
    }
  }
  writeHyperparameters(hyperparameter_file, false);
  model->writeHyperparameters(hyperparameter_file, true);

  std::cout << std::endl;

  initialize();
};

/**
 * @brief Initialize the run.
 *
 * Set up the training and validation MSA, initialize the sample and energy
 * data structures, create instances of the MSAStats, Model, SampleStats, and
 * Sampler classes, and compute the MSA statistics.
 */
void
Sim::initialize(void)
{
  // Subset the MSA if cross-validation is enabled but not validation MSA was
  // provided.
  if (!msa_validate) {
    if (cross_validate) {
      std::tie(msa_train, msa_validate) =
        msa_train->partitionAlignment(validation_seqs, random_seed);
      msa_validate->writeSequenceWeights("msa_weights_validate.txt");
      msa_validate->writeMatrix("msa_numerical_validate.txt");
    }
  }

  msa_train_energies = arma::Col<double>(msa_train->M, arma::fill::zeros);
  msa_train->writeMatrix("msa_numerical.txt");
  msa_train->writeSequenceWeights("msa_weights.txt");

  // Compute stats for training MSA
  msa_train_stats = std::make_shared<MSAStats>(msa_train, true);
  msa_train_stats->writeFrequency1p("msa_stat_1p.bin");
  msa_train_stats->writeFrequency2p("msa_stat_2p.bin");
  std::cout << std::endl;

  if (msa_validate) {
    msa_validate_energies =
      arma::Col<double>(msa_validate->M, arma::fill::zeros);
    msa_validate->writeMatrix("msa_numerical_validate.txt");
    msa_validate->writeSequenceWeights("msa_weights_validate.txt");

    // Compute stats for validation MSA
    msa_validate_stats = std::make_shared<MSAStats>(msa_validate, true);
    msa_validate_stats->writeFrequency1p("msa_stat_1p_validate.bin");
    msa_validate_stats->writeFrequency2p("msa_stat_2p_validate.bin");
    std::cout << std::endl;

    if (!cross_validate) {
      std::cout
        << "NOTE: Ignoring cross_validate=false due to provided validation MSA."
        << std::endl;
    }
  }

  // If using stochastic sampling, set M to match the number of sequences.
  if (use_ss) {
    samples_per_walk = 1;
    burn_between_start = 0;
    burn_between = 0;
    walkers = static_cast<int>((round(msa_train_stats->getM())));
  }

  int N = msa_train_stats->getN();
  int Q = msa_train_stats->getQ();

  // Set stat objects and initialize model
  model->setMSAStats(msa_train_stats, msa_validate_stats);
  if (step_offset == 0) {
    model->initialize();
  } else {
    model->restore(step_offset, output_binary);
  }

  // Initialize sample data structure and statistic class
  if (samples_per_walk > 1) {
    samples_3d =
      arma::Cube<int>(samples_per_walk, N, walkers, arma::fill::zeros);
    sample_stats = std::make_shared<SampleStats3D>(
      &samples_3d, &(model->params), &(model->params_prev));
    sample_stats->setMixingTime(burn_between);
  } else {
    samples_2d = arma::Mat<int>(walkers, N, arma::fill::zeros);
    sample_stats = std::make_shared<SampleStats2D>(
      &samples_2d, &(model->params), &(model->params_prev));
  }

  model->setSampleStats(sample_stats);
  model->setStep(step_offset);

  // Initialize the sampler
  sampler = std::make_shared<Sampler>(N, Q, &(model->params));
};

/**
 * @brief Load sampler hyperparameters from config file.
 *
 * @param file_name config file string
 */
void
Sim::loadHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  if (file.is_open()) {
    std::string line;
    bool reading_bmdca_section = false;
    while (std::getline(file, line)) {
      if (line.empty()) {
        reading_bmdca_section = false;
        continue;
      } else if (line[0] == '#') {
        reading_bmdca_section = false;
        continue;
      } else if (line[0] == '[') {
        if (line == "[bmDCA]") {
          reading_bmdca_section = true;
          continue;
        } else {
          reading_bmdca_section = false;
          continue;
        }
      }
      if (reading_bmdca_section) {
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
Sim::compareHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool all_same = true;
  if (file.is_open()) {
    std::string line;
    bool reading_bmdca_section = false;
    while (std::getline(file, line)) {
      if (line.empty()) {
        continue;
      } else if (line[0] == '#') {
        continue;
      } else if (line[0] == '[') {
        if (line == "[bmDCA]") {
          reading_bmdca_section = true;
          continue;
        } else {
          reading_bmdca_section = false;
          continue;
        }
      }
      if (reading_bmdca_section) {
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
Sim::setHyperparameter(const std::string key, const std::string value)
{
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "step_max") {
    step_max = std::stoi(value);
  } else if (key == "save_period") {
    save_period = std::stoi(value);
  } else if (key == "save_best_steps") {
    if (value.size() == 1) {
      save_best_steps = (std::stoi(value) == 1);
    } else {
      save_best_steps = (value == "true");
    }
  } else if (key == "stop_mode") {
    stop_mode = value;
  } else if (key == "stop_threshold") {
    stop_threshold = std::stod(value);
  } else if (key == "train_mode") {
    train_mode = value;
  } else if (key == "cross_validate") {
    if (value.size() == 1) {
      cross_validate = (std::stoi(value) == 1);
    } else {
      cross_validate = (value == "true");
    }
  } else if (key == "validation_seqs") {
    validation_seqs = std::stoi(value);
  } else if (key == "random_seed") {
    random_seed = std::stol(value);
  } else if (key == "output_binary") {
    if (value.size() == 1) {
      output_binary = (std::stoi(value) == 1);
    } else {
      output_binary = (value == "true");
    }
  } else if (key == "update_rule") {
    update_rule = value;
  } else if (key == "burn_in_start") {
    burn_in_start = std::stoi(value);
  } else if (key == "burn_between_start") {
    burn_between_start = std::stoi(value);
  } else if (key == "update_burn_time") {
    if (value.size() == 1) {
      update_burn_time = (std::stoi(value) == 1);
    } else {
      update_burn_time = (value == "true");
    }
  } else if (key == "adapt_up_time") {
    adapt_up_time = std::stod(value);
  } else if (key == "adapt_down_time") {
    adapt_down_time = std::stod(value);
  } else if (key == "step_importance_max") {
    step_importance_max = std::stoi(value);
  } else if (key == "coherence_min") {
    coherence_min = std::stod(value);
  } else if (key == "use_ss") {
    if (value.size() == 1) {
      use_ss = (std::stoi(value) == 1);
    } else {
      use_ss = (value == "true");
    }
  } else if (key == "walkers") {
    walkers = std::stoi(value);
  } else if (key == "samples_per_walk") {
    samples_per_walk = std::stoi(value);
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

/**
 * @brief Check additional constraints for hyperparameter values.
 *
 * Check any constraints on the possible values for the run hyperparameters.
 */
void
Sim::checkHyperparameters(void)
{
  if ((validation_seqs < 1) & cross_validate) {
    std::cerr << "ERROR: Incompatible number of sequecnes for cross-validation."
              << std::endl;
  }
};

/**
 * @brief Write stored hyperparameters to file
 *
 * @param output_file config file string
 * @param append flag for whether or not to append the hyperparameters to the
 * file or overwrite it
 */
void
Sim::writeHyperparameters(std::string output_file, bool append)
{
  std::ofstream stream;
  if (append) {
    stream.open(output_file, std::ofstream::out | std::ofstream::app);
  } else {
    stream.open(output_file, std::ofstream::out | std::ofstream::trunc);
  }

  // Header
  stream << "[bmDCA]" << std::endl;

  // BM settings
  stream << "step_max=" << step_max << std::endl;
  stream << "save_period=" << save_period << std::endl;
  stream << "save_best_steps=" << save_best_steps << std::endl;
  stream << "stop_mode=" << stop_mode << std::endl;
  stream << "stop_threshold=" << stop_threshold << std::endl;
  stream << "train_mode=" << train_mode << std::endl;
  stream << "cross_validate=" << cross_validate << std::endl;
  stream << "validation_seqs=" << validation_seqs << std::endl;
  stream << "random_seed=" << random_seed << std::endl;
  stream << "output_binary=" << output_binary << std::endl;
  stream << "update_rule=" << update_rule << std::endl;
  stream << "burn_in_start=" << burn_in_start << std::endl;
  stream << "burn_between_start=" << burn_between_start << std::endl;
  stream << "update_burn_time=" << update_burn_time << std::endl;
  stream << "adapt_up_time=" << adapt_up_time << std::endl;
  stream << "adapt_down_time=" << adapt_down_time << std::endl;
  stream << "step_importance_max=" << step_importance_max << std::endl;
  stream << "coherence_min=" << coherence_min << std::endl;
  stream << "use_ss=" << use_ss << std::endl;
  stream << "walkers=" << walkers << std::endl;
  stream << "samples_per_walk=" << samples_per_walk << std::endl;

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
Sim::compareHyperparameter(std::string key, std::string value)
{
  bool same = true;
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "step_max") {
  } else if (key == "save_period") {
  } else if (key == "save_best_steps") {
  } else if (key == "stop_mode") {
  } else if (key == "stop_threshold") {
  } else if (key == "train_mode") {
    same = same & (train_mode == value);
  } else if (key == "cross_validate") {
    if (value.size() == 1) {
      same = same & (cross_validate == (std::stoi(value) == 1));
    } else {
      same = same & (cross_validate == (value == "true"));
    }
  } else if (key == "validation_seqs") {
    same = same & (validation_seqs == std::stoi(value));
  } else if (key == "random_seed") {
    same = same & (random_seed == std::stol(value));
  } else if (key == "output_binary") {
    if (value.size() == 1) {
      same = same & (output_binary == (std::stoi(value) == 1));
    } else {
      same = same & (output_binary == (value == "true"));
    }
  } else if (key == "update_rule") {
    same = same & (update_rule == value);
  } else if (key == "burn_in_start") {
    same = same & (burn_in_start == std::stoi(value));
  } else if (key == "burn_between_start") {
    same = same & (burn_between_start == std::stoi(value));
  } else if (key == "update_burn_time") {
    if (value.size() == 1) {
      same = same & (update_burn_time == (std::stoi(value) == 1));
    } else {
      same = same & (update_burn_time == (value == "true"));
    }
  } else if (key == "adapt_up_time") {
    same = same & (adapt_up_time == std::stod(value));
  } else if (key == "adapt_down_time") {
    same = same & (adapt_down_time == std::stod(value));
  } else if (key == "step_importance_max") {
    same = same & (step_importance_max == std::stoi(value));
  } else if (key == "coherence_min") {
    same = same & (coherence_min == std::stod(value));
  } else if (key == "use_ss") {
    if (value.size() == 1) {
      same = same & (use_ss == (std::stoi(value) == 1));
    } else {
      same = same & (use_ss == (value == "true"));
    }
  } else if (key == "walkers") {
    same = same & (walkers == std::stoi(value));
  } else if (key == "samples_per_walk") {
    same = same & (samples_per_walk == std::stoi(value));
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
 *
 * @return (bool) flag for whether the necessary files were found
 */
bool
Sim::isValidStep(int step)
{
  bool valid = false;
  if (checkFileExists("samples_" + std::to_string(step) + ".txt") &
      checkFileExists("samples_stat_1p_" + std::to_string(step) + ".bin") &
      checkFileExists("samples_stat_2p_" + std::to_string(step) + ".bin") &
      checkFileExists("energies_" + std::to_string(step) + ".txt")) {
    valid = true;
  } else if (checkFileExists("samples_" + std::to_string(step) + ".txt") &
             checkFileExists("samples_stat_1p_" + std::to_string(step) +
                             ".txt") &
             checkFileExists("samples_stat_2p_" + std::to_string(step) +
                             ".txt") &
             checkFileExists("energies_" + std::to_string(step) + ".txt")) {
    valid = true;
  }
  return valid;
};

/**
 * @brief Find the last valid step in order to restore run state.
 *
 * Search the output directory and parse file names for step numbers, and then
 * check that the requisite data exists for that step number.
 */
void
Sim::setStepOffset(void)
{
  DIR* dp;
  struct dirent* dirp;

  dp = opendir(".");

  std::vector<int> steps;
  std::vector<int> invalid_steps;
  while ((dirp = readdir(dp)) != nullptr) {
    std::string fname = dirp->d_name;
    if (output_binary) {
      if (fname.find("samples_") == std::string::npos)
        continue;

      const std::regex re_samples("samples_([0-9]+)\\.txt");
      std::smatch match_samples;

      int idx = -1;
      if (std::regex_match(fname, match_samples, re_samples)) {
        if (match_samples.size() == 2) {
          idx = std::stoi(match_samples[1].str());
        }
      }

      if (isValidStep(idx) & model->isValidStep(idx)) {
        steps.push_back(idx);
      } else {
        if (idx > -1) {
          invalid_steps.push_back(idx);
        }
      }
    }
  }
  closedir(dp);

  // Clear out steps with missing files.
  for (auto it_bad = invalid_steps.begin(); it_bad != invalid_steps.end();
       ++it_bad) {
    bool delete_files = true;
    for (auto it_good = steps.begin(); it_good != steps.end(); ++it_good) {
      if (*it_good == *it_bad - 1) {
        delete_files = false;
        break;
      }
      if (*it_good == *it_bad + 1) {
        delete_files = false;
        break;
      }
    }

    if (delete_files) {
      std::cout << "missing data --- clearing step " << *it_bad << std::endl;
      deleteStep(*it_bad);
      model->deleteStep(*it_bad);
    }
  }

  // Pick out the most recent valid step.
  int max = -1;
  for (auto it = steps.begin(); it != steps.end(); ++it) {
    if (*it > max) {
      max = *it;
    }
  }
  if (max < 0) {
    step_offset = 0;
  } else {
    step_offset = max;
  }
};

/**
 * @brief Delete all data associated with a particular step.
 *
 * @param step iteration to delete
 */
void
Sim::deleteStep(int step)
{
  std::string file;

  file = "samples_" + std::to_string(step) + ".txt";
  if (checkFileExists(file))
    deleteFile(file);

  file = "energies_relax" + std::to_string(step) + ".txt";
  if (checkFileExists(file))
    deleteFile(file);

  file = "energies_" + std::to_string(step) + ".txt";
  if (checkFileExists(file))
    deleteFile(file);

  if (output_binary) {
    file = "samples_stat_1p_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_2p_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_1p_sigma_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_2p_sigma_" + std::to_string(step) + ".bin";
    if (checkFileExists(file))
      deleteFile(file);
  } else {
    file = "samples_stat_1p_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_2p_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_1p_sigma_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);

    file = "samples_stat_2p_sigma_" + std::to_string(step) + ".txt";
    if (checkFileExists(file))
      deleteFile(file);
  }
};

/**
 * @brief Restore a run from the last valid position.
 *
 * Parse the run log to re-load burn-times, etc. and reset the RNG to the state
 * it was in when the program terminated. Ensures results are deterministic.
 */
void
Sim::restoreRunState(void)
{
  std::ifstream stream(run_log_file);
  std::string line;
  std::getline(stream, line);
  std::string rng_state;

  while (!stream.eof()) {
    std::getline(stream, line);
    std::stringstream buffer(line);
    std::string field;

    std::getline(buffer, field, '\t');
    if (step_offset == std::stoi(field)) {
      std::vector<std::string> fields;
      std::stringstream ss;
      ss.str(line);
      std::string field;
      while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
      }

      burn_in = std::stoi(fields.at(3));
      if (samples_per_walk > 1) {
        burn_between = std::stoi(fields.at(4));
      }

      // column number for some fields depend on which flags are set
      if (samples_per_walk > 1) {
        train_err_tot_min = std::stod(fields.at(17));
        if (cross_validate) {
          validate_err_tot_min = std::stod(fields.at(22));
          rng_state = fields.at(23) + " " + fields.at(24) + " " + fields.at(25);
        } else {
          rng_state = fields.at(18) + " " + fields.at(19) + " " + fields.at(20);
        }
      } else {
        if (cross_validate) {
          train_err_tot_min = std::stod(fields.at(8));
          validate_err_tot_min = std::stod(fields.at(12));
          rng_state = fields.at(14) + " " + fields.at(15) + " " + fields.at(16);
        } else {
          train_err_tot_min = std::stod(fields.at(8));
          rng_state = fields.at(9) + " " + fields.at(10) + " " + fields.at(11);
        }
      }
      break;
    }
  }

  std::istringstream is(rng_state);
  is >> rng;
  if (is.fail()) {
    std::cerr << "ERROR: failed to restore RNG state." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

/**
 * @brief Run the inference loop.
 *
 * All heavy lifting is done here.
 */
void
Sim::run(void)
{
  std::cout << "initializing run... " << std::flush;

  arma::wall_clock step_timer;

  arma::wall_clock timer;
  timer.tic();

  // Instantiate the PCG random number generator and uniform random
  // distribution.
  rng.seed(random_seed);

  // Initialize the buffer.
  run_buffer = arma::Mat<double>(save_period, 24, arma::fill::zeros);
  rng_buffer = arma::Mat<uint64_t>(save_period, 3, arma::fill::zeros);
  int buffer_offset = 0;

  if (step_offset == 0) {
    initializeRunLog();
    burn_in = burn_in_start;
    burn_between = burn_between_start;
  } else if (step_offset >= step_max) {
    std::cout << "step " << step_max << " already reached... done."
              << std::endl;
    return;
  } else {
    restoreRunState();
    buffer_offset = (step_offset % save_period);
  }
  std::cout << timer.toc() << " sec" << std::endl;
  std::cout << std::endl;

  // Bootstrap MSA to get cutoff error for 'msaerr' stop mode.
  if (stop_mode == "msaerr") {
    std::cout << "bootstrapping msa... " << std::flush;
    timer.tic();
    msa_train_stats->computeErrorMSA(10);
    std::cout << timer.toc() << " sec" << std::endl;

    double error_avg = arma::mean(msa_train_stats->msa_rms);
    // double error_stddev = arma::stddev(msa_train_stats->msa_rms, 1);
    // stop_threshold = (error_avg - 2 * error_stddev) /
    //                   sqrt(M * count_max / msa_train_stats->getEffectiveM());
    stop_threshold = error_avg / sqrt(samples_per_walk * walkers /
                                      msa_train_stats->getEffectiveM());
    std::cout << "convergence threshold is " << stop_threshold << std::endl;
  } else if ((stop_mode == "stderr") | (stop_mode == "stderr_adj")) {
    stop_threshold =
      msa_train_stats->freq_rms / sqrt(samples_per_walk * walkers);
    std::cout << "convergence threshold is " << stop_threshold << std::endl;
  }

  // BM sampling loop
  for (step = 1 + step_offset; step <= step_max; step++) {
    step_timer.tic();
    std::cout << "Step: " << step << std::endl;
    model->setStep(step);

    // Sampling from MCMC (keep trying until correct properties found)
    bool flag_mc = true;
    while (flag_mc) {
      if (update_burn_time & (samples_per_walk == 1)) {
        std::cout << "setting burn time to... " << std::flush;
        timer.tic();
        bool flag_burn = true;
        while (flag_burn) {
          double burn_reps = 48;
          double burn_count = 2;
          arma::Mat<double> energy_burn =
            arma::Mat<double>(burn_count, burn_reps, arma::fill::zeros);

          sampler->sampleEnergies(
            &energy_burn, burn_reps, burn_count, burn_in, burn_in, rng());

          double e_start = arma::mean(energy_burn.row(0));
          double e_start_sigma = arma::stddev(energy_burn.row(0), 1);
          double e_end = arma::mean(energy_burn.row(burn_count - 1));
          double e_end_sigma = arma::stddev(energy_burn.row(burn_count - 1), 1);
          double e_err =
            sqrt((pow(e_start_sigma, 2) + pow(e_end_sigma, 2)) / burn_reps);

          bool flag_twaiting_up = true;
          bool flag_twaiting_down = true;

          // If the difference in average starting (burn-in) and ending
          // (burn-in + n_sequences x burn-between) energies along a trajectory
          // is less than twice the variance, don't increase burn-in time.
          if (e_start - e_end <= 2 * e_err) {
            flag_twaiting_up = false;
          }

          // If the difference in average ending (burn-in + n_sequences x
          // burn-between) and starting (burn-in) energies along a trajectory
          // is less than twice the variance, don't decrease burn-in time.
          if (e_start - e_end >= -2 * e_err) {
            flag_twaiting_down = false;
          }

          if (flag_twaiting_up) {
            burn_in = static_cast<int>(
              round(static_cast<double>(burn_in) * adapt_up_time));
          } else if (flag_twaiting_down) {
            burn_in = Max(static_cast<int>(round(static_cast<double>(burn_in) *
                                                 adapt_down_time)),
                          1);
          }
          if (!flag_twaiting_up) {
            flag_burn = false;
          }
        }
        std::cout << burn_in << "... " << timer.toc() << " sec" << std::endl;
      }

      // Draw from MCMC
      std::cout << "sampling model... " << std::flush;
      timer.tic();
      unsigned seed = rng();
      {
        std::stringstream ss;
        uint64_t rng_mult = 0, rng_inc = 0, rng_state = 0;
        ss << rng;
        ss >> rng_mult >> rng_inc >> rng_state;
        rng_buffer((step - 1) % save_period, 0) = rng_mult;
        rng_buffer((step - 1) % save_period, 1) = rng_inc;
        rng_buffer((step - 1) % save_period, 2) = rng_state;
      }
      if (samples_per_walk > 1) {
        if (update_rule == "mh") {
          sampler->sampleSequences(&samples_3d,
                                   walkers,
                                   samples_per_walk,
                                   burn_in,
                                   burn_between,
                                   seed);
        } else if (update_rule == "z-sqrt") {
          sampler->sampleSequencesZanellaSqrt(&samples_3d,
                                              walkers,
                                              samples_per_walk,
                                              burn_in,
                                              burn_between,
                                              seed);
        } else if (update_rule == "z-barker") {
          sampler->sampleSequencesZanellaBarker(&samples_3d,
                                                walkers,
                                                samples_per_walk,
                                                burn_in,
                                                burn_between,
                                                seed);
        } else {
          std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        if (update_rule == "mh") {
          sampler->sampleSequences(&samples_2d, walkers, burn_in, seed);
        } else if (update_rule == "z-sqrt") {
          sampler->sampleSequencesZanellaSqrt(
            &samples_2d, walkers, burn_in, seed);
        } else if (update_rule == "z-barker") {
          sampler->sampleSequencesZanellaBarker(
            &samples_2d, walkers, burn_in, seed);
        } else {
          std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      std::cout << timer.toc() << " sec" << std::endl;

      std::cout << "computing sample energies and correlations... "
                << std::flush;
      timer.tic();
      sample_stats->computeStatsExtra();
      std::cout << timer.toc() << " sec" << std::endl;

      // Run checks and alter burn-in and wait times
      if (update_burn_time & (samples_per_walk > 1)) {
        flag_mc = checkErgodicity();
        if (flag_mc) {
          sample_stats->setMixingTime(burn_between);
        }
      } else {
        flag_mc = false;
      }
    }

    if (step_importance_max > 1) {
      // Importance sampling loop
      std::cout << "starting importance sampling loop" << std::endl;
      for (int step_importance = 0; step_importance < step_importance_max;
           step_importance++) {
        std::cout << "importance sampling step " << step_importance << "... "
                  << std::flush;
        timer.tic();
        sample_stats->computeStatsImportance();
        std::cout << timer.toc() << " sec" << std::endl;

        std::cout << "updating model with 1p and 2p frequencies... "
                  << std::flush;
        timer.tic();
        model->update();
        std::cout << timer.toc() << " sec" << std::endl;

        double coherence;
        if (samples_per_walk > 1) {
          arma::Col<double> stats = sample_stats->getStats();
          coherence = stats(16); // Z_ratio
        } else {
          arma::Col<double> stats = sample_stats->getStats();
          coherence = stats(2); // Z_ratio
        }
        if (coherence > coherence_min && 1.0 / coherence > coherence_min) {
          break;
        }
      }
    } else {
      std::cout << "computing sample 1p and 2p statistics... " << std::flush;
      timer.tic();
      sample_stats->computeStats();
      std::cout << timer.toc() << " sec" << std::endl;

      std::cout << "updating model with 1p and 2p frequencies... "
                << std::flush;
      timer.tic();
      model->update();
      std::cout << timer.toc() << " sec" << std::endl;
    }

    computeMSAEnergies(&msa_train_energies,
                       msa_train,
                       std::make_shared<potts_model>(model->params));
    if (msa_validate) {
      computeMSAEnergies(&msa_validate_energies,
                         msa_validate,
                         std::make_shared<potts_model>(model->params));
      diff_avg_energy =
        arma::mean(msa_train_energies) - arma::mean(msa_validate_energies);
    }

    double train_err_1p = model->train_error_1p;
    double train_err_2p = model->train_error_2p;
    double train_err_tot = train_err_1p + train_err_2p;
    double validate_err_1p = model->validation_error_1p;
    double validate_err_2p = model->validation_error_2p;
    double validate_err_tot = validate_err_1p + validate_err_2p;

    bool new_min_found = true;
    if (train_err_tot >= train_err_tot_min) {
      new_min_found = new_min_found & false;
    } else {
      train_err_tot_min = train_err_tot;
    }
    if (cross_validate) {
      if (validate_err_tot >= validate_err_tot_min) {
        new_min_found = new_min_found & false;
      } else {
        validate_err_tot_min = validate_err_tot;
      }
    }

    run_buffer((step - 1) % save_period, 0) = step;
    run_buffer((step - 1) % save_period, 1) = walkers;
    run_buffer((step - 1) % save_period, 2) = samples_per_walk;
    run_buffer((step - 1) % save_period, 3) = burn_in;
    run_buffer((step - 1) % save_period, 4) = burn_between;

    double total_corr;

    // Get sample statistics
    if (samples_per_walk > 1) {
      arma::Col<double> stats = sample_stats->getStats();

      double e_start = stats.at(0);
      double e_start_sigma = stats.at(1);
      double e_end = stats.at(2);
      double e_end_sigma = stats.at(3);
      double e_err = stats.at(4);

      total_corr = stats.at(5);
      double auto_corr = stats.at(7);
      double cross_corr = stats.at(8);
      double auto_cross_err = stats.at(13);

      run_buffer((step - 1) % save_period, 5) = total_corr;
      run_buffer((step - 1) % save_period, 6) = auto_corr;
      run_buffer((step - 1) % save_period, 7) = cross_corr;
      run_buffer((step - 1) % save_period, 8) = auto_cross_err;
      run_buffer((step - 1) % save_period, 9) = e_start;
      run_buffer((step - 1) % save_period, 10) = e_start_sigma;
      run_buffer((step - 1) % save_period, 11) = e_end;
      run_buffer((step - 1) % save_period, 12) = e_end_sigma;
      run_buffer((step - 1) % save_period, 13) = e_err;
    } else {
      arma::Col<double> stats = sample_stats->getStats();
      total_corr = stats.at(0);
      run_buffer((step - 1) % save_period, 5) = total_corr;
    }

    run_buffer((step - 1) % save_period, 14) = train_err_1p;
    run_buffer((step - 1) % save_period, 15) = train_err_2p;
    run_buffer((step - 1) % save_period, 16) = train_err_tot;
    run_buffer((step - 1) % save_period, 17) = train_err_tot_min;

    run_buffer((step - 1) % save_period, 18) = validate_err_1p;
    run_buffer((step - 1) % save_period, 19) = validate_err_2p;
    run_buffer((step - 1) % save_period, 20) = validate_err_tot;
    run_buffer((step - 1) % save_period, 21) = validate_err_tot_min;
    run_buffer((step - 1) % save_period, 22) = diff_avg_energy;

    bool converged = false;
    if (train_err_tot < stop_threshold) {
      converged = true;
    } else if ((stop_mode == "stderr") | (stop_mode == "stderr_adj")) {
      double std_err = stop_threshold;
      if (stop_mode == "stderr_adj") {
        std_err = sqrt((1 + total_corr) / (1 - total_corr)) * std_err;
      }
      if (train_err_tot < std_err) {
        converged = true;
      }
    }

    if (converged) {
      run_buffer((step - 1) % save_period, 23) = step_timer.toc();
      std::cout << "converged! writing final results... " << std::flush;
      writeRunLog(step % save_period, buffer_offset);
      writeStep(step);
      std::cout << "done" << std::endl;
      return;
    }

    run_buffer((step - 1) % save_period, 23) = step_timer.toc();

    // Save parameters
    if (step % save_period == 0) {
      std::cout << "writing step " << step << "... " << std::flush;
      timer.tic();
      writeRunLog(step % save_period, buffer_offset);
      buffer_offset = 0;
      writeStep(step);
      std::cout << timer.toc() << " sec" << std::endl;
    } else if (new_min_found & (step > save_period) & save_best_steps) {
      std::cout << "close... writing step... " << std::flush;
      run_buffer((step - 1) % save_period, 23) = step_timer.toc();
      writeRunLog(step % save_period, buffer_offset, true);
      buffer_offset = step % save_period;
      writeStep(step);
      std::cout << "done" << std::endl;
    }

    std::cout << std::endl;
  }

  if (step_offset != step_max) {
    std::cout << "writing final results... " << std::flush;
    if ((step_max % save_period) != 0) {
      writeRunLog(step_max % save_period);
    }
    writeStep(step_max);
  } else {
    std::cout << "all " << step_offset << " steps already... " << std::flush;
  }

  std::cout << "done" << std::endl;
  return;
};

/**
 * @brief Update the burn-in and burn-between times.
 *
 * @return (bool) flag for whether to resample the sequences.
 *
 * This function is called when M > 1 sequences are sampled per MC trajectory.
 * If sequences along a trajectory are too correlated, burn-between time is
 * increased. If energy after burn-in and after burn-in + M x burn-between is
 * too much lower, then burn-in time is increased.
 */
bool
Sim::checkErgodicity(void)
{
  arma::Col<double> stats = sample_stats->getStats();

  double e_start = stats.at(0); // average starting energy for MC trajectory
  double e_end = stats.at(2);   // average ending energy
  double e_err = stats.at(4);   // combined starting and ending variance

  double auto_corr =
    stats.at(7); // sequence correlations for trajectory start/end
  double cross_corr = stats.at(8); // correlations between trajectories
  double check_corr =
    stats.at(9); // correlations for mid-trajectory (length/10) and end
  double cross_check_err = stats.at(14); // combined auto+cross variance
  double auto_cross_err = stats.at(13);  // combined check+cross variance

  bool flag_deltat_up = true;     // flag to increase burn-between time
  bool flag_deltat_down = true;   // flag to decrease burn-between time
  bool flag_twaiting_up = true;   // flag to increase burn-in time
  bool flag_twaiting_down = true; // flag to decrease burn-in time

  // If the difference in sequence correlations within and between trajectories
  // is less than the variance, don't increase waiting time.
  if (check_corr - cross_corr <= cross_check_err) {
    flag_deltat_up = false;
  }

  // If the difference in sequence correlations within and between trajectories
  // is grater than the variance, don't decrease waiting time.
  if (auto_corr - cross_corr >= auto_cross_err) {
    flag_deltat_down = false;
  }

  // If the difference in average starting (burn-in) and ending (burn-in +
  // n_sequences x burn-between) energies along a trajectory is less than twice
  // the variance, don't increase burn-in time.
  if (e_start - e_end <= 2 * e_err) {
    flag_twaiting_up = false;
  }

  // If the difference in average ending (burn-in + n_sequences x burn-between)
  // and starting (burn-in) energies along a trajectory is less than twice the
  // variance, don't decrease burn-in time.
  if (e_start - e_end >= -2 * e_err) {
    flag_twaiting_down = false;
  }

  if (flag_deltat_up) {
    burn_between = static_cast<int>(
      round(static_cast<double>(burn_between) * adapt_up_time));
    std::cout << "increasing burn-between time to " << burn_between
              << std::endl;
  } else if (flag_deltat_down) {
    burn_between = Max(static_cast<int>(round(
                         static_cast<double>(burn_between) * adapt_down_time)),
                       1);
    std::cout << "decreasing burn-between time to " << burn_between
              << std::endl;
  }

  if (flag_twaiting_up) {
    burn_in =
      static_cast<int>(round(static_cast<double>(burn_in) * adapt_up_time));
    std::cout << "increasing burn-in time to " << burn_in << std::endl;
  } else if (flag_twaiting_down) {
    burn_in = Max(
      static_cast<int>(round(static_cast<double>(burn_in) * adapt_down_time)),
      1);
    std::cout << "decreasing burn-in time to " << burn_in << std::endl;
  }

  // If burn-in or burn-between times were increased, then the sequences need
  // to be resampled to ensure proper mixing.
  bool flag_mc = true;
  if (!flag_deltat_up && !flag_twaiting_up) {
    flag_mc = false;
  }
  return flag_mc;
};

/**
 * @brief Compute MSA sequence energies using the current parameters.
 *
 * @param energies pointer to energy data structure
 * @param msa pointer to the MSA
 * @param params pointer to the Potts parameters
 */
void
Sim::computeMSAEnergies(arma::Col<double>* energies,
                        std::shared_ptr<MSA> msa,
                        std::shared_ptr<potts_model> params)
{
  int N = msa->N;
  int M = msa->M;
#pragma omp parallel
  {
#pragma omp for
    for (int seq = 0; seq < M; seq++) {
      double E = 0;
      for (int i = 0; i < N; i++) {
        E -= params->h(msa->alignment(seq, i), i);
        for (int j = i + 1; j < N; j++) {
          E -= params->J(i, j)(msa->alignment(seq, i), msa->alignment(seq, j));
        }
      }
      (*energies)(seq) = E;
    }
  }
};

/**
 * @brief Write the data for the current step.
 *
 * @param step current inference iteration
 *
 * Writes the model parameters necessary for restarting runs, along with the
 * sequence energies and statistics for the MSA and sampled sequences.
 */
void
Sim::writeStep(int step)
{
  model->writeStep(step, output_binary);
  sample_stats->writeStep(step, output_binary);
  writeMSAEnergies(step);
};

/**
 * @brief Write the current data.
 *
 * @param id label for the output data
 *
 * Writes the model parameters along with the sequence energies and statistics
 * for the MSA and sampled sequences.
 */
void
Sim::writeData(const std::string id)
{
  model->writeData(id, output_binary);
  sample_stats->writeData(id, output_binary);
  writeMSAEnergies(id);
};

/**
 * @brief Write the sequence energies for the MSA using current parameters.
 *
 * @param step current iteration
 *
 * Writes both the training and validation (if applicable) sequence energies as
 * computed using the currently inferred Potts model.
 */
void
Sim::writeMSAEnergies(int step)
{
  {
    std::ofstream output_stream("msa_energies_" + std::to_string(step) +
                                ".txt");

    int M = msa_train->M;
    for (int m = 0; m < M; m++) {
      output_stream << msa_train_energies(m) << std::endl;
    }

    output_stream.close();
  }
  if (msa_validate) {
    std::ofstream output_stream("msa_energies_validate_" +
                                std::to_string(step) + ".txt");

    int M = msa_validate->M;
    for (int m = 0; m < M; m++) {
      output_stream << msa_validate_energies(m) << std::endl;
    }

    output_stream.close();
  }
};

/**
 * @brief Write the sequence energies for the MSA using current parameters.
 *
 * @param id output label
 *
 * Writes both the training and validation (if applicable) sequence energies as
 * computed using the currently inferred Potts model.
 */
void
Sim::writeMSAEnergies(const std::string id)
{
  {
    std::ofstream output_stream("msa_energies_" + id + ".txt");

    int M = msa_train->M;
    for (int m = 0; m < M; m++) {
      output_stream << msa_train_energies(m) << std::endl;
    }

    output_stream.close();
  }
  if (msa_validate) {
    std::ofstream output_stream("msa_energies_validate_" + id + ".txt");

    int M = msa_validate->M;
    for (int m = 0; m < M; m++) {
      output_stream << msa_validate_energies(m) << std::endl;
    }

    output_stream.close();
  }
};

/**
 * @brief Write the column headers for the run log.
 */
void
Sim::initializeRunLog()
{
  std::ofstream stream{ run_log_file, std::ios_base::out };
  stream << "step" << "\t" << "walkers" << "\t" << "samples-per-walk" << "\t"
         << "burn-in" << "\t";
  if (samples_per_walk > 1) {
    stream << "burn-between" << "\t";
  }
  stream << "total-corr" << "\t";
  if (samples_per_walk > 1) {
    stream << "auto-corr" << "\t" << "cross-corr" << "\t" << "auto-cross-err"
           << "\t" << "energy-start-avg" << "\t" << "energy-start-sigma" << "\t"
           << "energy-end-avg" << "\t" << "energy-end-sigma" << "\t"
           << "energy-err" << "\t";
  }
  stream << "train-err-1p" << "\t" << "train-err-2p" << "\t" << "train-err-tot"
         << "\t" << "train-err-tot-min" << "\t";
  if (cross_validate) {
    stream << "validate-err-1p" << "\t" << "validate-err-2p" << "\t"
           << "validate-err-tot" << "\t" << "validate-err-tot-min" << "\t"
           << "diff-avg-energy" << "\t";
  }
  stream << "rng-mult" << "\t" << "rng-inc" << "\t" << "rng-state" << "\t"
         << "duration" << std::endl;
  stream.close();
};

/**
 * @brief Flush the run buffer data to the run log.
 *
 * @param current_step current iteration
 * @param offset step offset (if the buffer isn't full)
 * @param keep (bool) flag to keep the data in the buffer
 */
void
Sim::writeRunLog(int current_step, int offset, bool keep)
{
  std::ofstream stream{ run_log_file, std::ios_base::app };

  int n_entries;
  if (current_step == 0) {
    n_entries = save_period;
  } else {
    n_entries = current_step;
  }
  for (int i = 0 + offset; i < n_entries; i++) {
    stream << static_cast<int>(run_buffer(i, 0)) << "\t";
    stream << static_cast<int>(run_buffer(i, 1)) << "\t";
    stream << static_cast<int>(run_buffer(i, 2)) << "\t";
    stream << static_cast<int>(run_buffer(i, 3)) << "\t";
    if (samples_per_walk > 1) {
      stream << run_buffer(i, 4) << "\t";
    }
    stream << run_buffer(i, 5) << "\t";
    if (samples_per_walk > 1) {
      stream << run_buffer(i, 6) << "\t";
      stream << run_buffer(i, 7) << "\t";
      stream << run_buffer(i, 8) << "\t";
      stream << run_buffer(i, 9) << "\t";
      stream << run_buffer(i, 10) << "\t";
      stream << run_buffer(i, 11) << "\t";
      stream << run_buffer(i, 12) << "\t";
      stream << run_buffer(i, 13) << "\t";
    }
    stream << run_buffer(i, 14) << "\t";
    stream << run_buffer(i, 15) << "\t";
    stream << run_buffer(i, 16) << "\t";
    stream << run_buffer(i, 17) << "\t";
    if (cross_validate) {
      stream << run_buffer(i, 18) << "\t";
      stream << run_buffer(i, 19) << "\t";
      stream << run_buffer(i, 20) << "\t";
      stream << run_buffer(i, 21) << "\t";
      stream << run_buffer(i, 22) << "\t";
    }
    stream << rng_buffer(i, 0) << "\t";
    stream << rng_buffer(i, 1) << "\t";
    stream << rng_buffer(i, 2) << "\t";
    stream << run_buffer(i, 23) << std::endl;
  }
  if (!keep) {
    run_buffer.zeros();
    rng_buffer.zeros();
  }
  stream.close();
};
