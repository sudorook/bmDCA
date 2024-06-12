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

#include "generator.hpp"

#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"

/**
 * @brief Generator constructor.
 *
 * @param params input parameters
 * @param n number of positions
 * @param q number of states
 * @param config_file file string for hyperparameter file
 */
Generator::Generator(potts_model params, int n, int q, std::string config_file)
  : N(n)
  , Q(q)
  , model(params)
{
  if (config_file.length() != 0) {
    loadParameters(config_file);
  }
};

/**
 * @brief Generator destructor.
 */
Generator::~Generator(void)
{
  delete sample_stats;
  delete sampler;
};

/**
 * @brief Function for where hyperparameter checks can go...
 *
 * Empty for now, but if checks for whether the input hyperparameters are
 * sensible are desired in future, this is where they go.
 */
void
Generator::checkParameters(void){};

/**
 * @brief Load sampler hyperparameters from config file.
 *
 * @param file_name config file string
 */
void
Generator::loadParameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
        continue;
      } else if (line[0] == '[') {
        if (line == "[bmDCA_sample]") {
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
        setParameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

/**
 * @brief Set a hyperparameter to a specific value.
 *
 * @param key hyperparameter to set
 * @param value value at which to set hyperparameter
 */
void
Generator::setParameter(std::string key, std::string value)
{
  if (key == "resample_max") {
    resample_max = std::stoi(value);
  } else if (key == "random_seed") {
    random_seed = std::stoi(value);
  } else if (key == "save_interim_samples") {
    if (value.size() == 1) {
      save_interim_samples = (std::stoi(value) == 1);
    } else {
      save_interim_samples = (value == "true");
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
  } else if (key == "temperature") {
    temperature = std::stod(value);
  }
};

/**
 * @brief Write sampled sequences in a FASTA files.
 *
 * @param output_file file string for output file
 *
 * Converts numerical sequences to letters and stores the output in FASTA
 * format.
 */
void
Generator::writeAASequences(std::string output_file)
{
  // Label each sequence such that each sequence has a unique number in the
  // FASTA header.  If more than 1 sequence was sampled in each trajectory, use
  // replicate x count_per_replicate + count as the label. Otherwise, just use
  // the array index number.
  if (samples_per_walk > 1) {
    int M = samples_3d.n_rows;
    int N = samples_3d.n_cols;
    int reps = samples_3d.n_slices;

    std::ofstream output_stream(output_file);

    char aa;
    for (int rep = 0; rep < reps; rep++) {
      for (int m = 0; m < M; m++) {
        output_stream << ">sample" << m * rep + m << std::endl;
        for (int n = 0; n < N; n++) {
          aa = convertAA(samples_3d(m, n, rep));
          if (aa != '\0') {
            output_stream << aa;
          }
        }
        output_stream << std::endl << std::endl;
      }
    }
  } else {
    int M = samples_2d.n_rows;
    int N = samples_2d.n_cols;

    std::ofstream output_stream(output_file);

    char aa;
    for (int m = 0; m < M; m++) {
      output_stream << ">sample" << m << std::endl;
      for (int n = 0; n < N; n++) {
        aa = convertAA(samples_2d(m, n));
        if (aa != '\0') {
          output_stream << aa;
        }
      }
      output_stream << std::endl << std::endl;
    }
  }
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
Generator::checkErgodicity(void)
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
    burn_between = (int)(round((double)burn_between * adapt_up_time));
    std::cout << "increasing burn-between time to " << burn_between
              << std::endl;
  } else if (flag_deltat_down) {
    burn_between = Max((int)(round((double)burn_between * adapt_down_time)), 1);
    std::cout << "decreasing burn-between time to " << burn_between
              << std::endl;
  }

  if (flag_twaiting_up) {
    burn_in = (int)(round((double)burn_in * adapt_up_time));
    std::cout << "increasing burn-in time to " << burn_in << std::endl;
  } else if (flag_twaiting_down) {
    burn_in = Max((int)(round((double)burn_in * adapt_down_time)), 1);
    std::cout << "decreasing burn-in time to " << burn_in << std::endl;
  }

  // If burn-in or burn-between times were increased, then the sequences need
  // to be resampled to ensure proper mixing.
  bool flag_mc = true;
  if (not flag_deltat_up and not flag_twaiting_up) {
    flag_mc = false;
  }
  return flag_mc;
};

/**
 * @brief Estimate the appropriate burn-in time.
 *
 * Function is called when M==1 sequences are sampled per MC trajectory. It
 * samples a small set of dummy sequences, checking the value of burn-in where
 * mean energy after one burn-in length stable relative to multiples of that
 * burn-in lengthn.
 */
void
Generator::estimateBurnTime(void)
{
  bool flag_burn = true;
  while (flag_burn) {
    double burn_reps = 24;
    double burn_count = 4;
    arma::Mat<double> energy_burn =
      arma::Mat<double>(burn_count, burn_reps, arma::fill::zeros);

    sampler->sampleEnergies(
      &energy_burn, burn_reps, burn_count, burn_in, burn_in, rng());

    double e_start = arma::mean(energy_burn.row(0));
    double e_start_sigma = arma::stddev(energy_burn.row(0), 1);
    double e_end =
      (arma::mean(energy_burn.row(burn_count - 1)) +
       arma::mean(energy_burn.row(burn_count - 2))) /
      2; // lazy way to double the number of sequences from the end step
    double e_end_sigma =
      sqrt(pow(arma::stddev(energy_burn.row(burn_count - 1), 1), 2) +
           pow(arma::stddev(energy_burn.row(burn_count - 2), 1), 2));
    double e_err =
      sqrt((pow(e_start_sigma, 2) / burn_reps + pow(e_end_sigma, 2)) / 2 /
           burn_reps);

    bool flag_twaiting_up = true;
    bool flag_twaiting_down = true;

    // If the difference in average starting (burn-in) and ending (burn-in +
    // n_sequences x burn-between) energies along a trajectory is less than
    // twice the variance, don't increase burn-in time.
    if (e_start - e_end <= 2 * e_err) {
      flag_twaiting_up = false;
    }

    // If the difference in average ending (burn-in + n_sequences x
    // burn-between) and starting (burn-in) energies along a trajectory is less
    // than twice the variance, don't decrease burn-in time.
    if (e_start - e_end >= -2 * e_err) {
      flag_twaiting_down = false;
    }

    if (flag_twaiting_up) {
      burn_in = (int)(round((double)burn_in * adapt_up_time));
    } else if (flag_twaiting_down) {
      burn_in = Max((int)(round((double)burn_in * adapt_down_time)), 1);
    }

    // If burn-in was increased, then keep resample at the longer time.
    if (!flag_twaiting_up) {
      flag_burn = false;
    }
  }
};

/**
 * @brief Write the output sequences in numerical format.
 *
 * @param output_file file string for output file
 */
void
Generator::writeNumericalSequences(std::string output_file)
{
  int idx = output_file.find_last_of(".");
  std::string raw_file = output_file.substr(0, idx);

  sample_stats->writeSamples(raw_file + "_numerical.txt");
  sample_stats->writeSampleEnergies(raw_file + "_energies.txt");
}

/**
 * @brief Run the sampler.
 *
 * @param n_indep_runs number of independent trajectories
 * @param n_per_run number of samples per trajectory
 * @param output_file file string for the output file for the sequences
 */
void
Generator::run(int n_indep_runs, int n_per_run, std::string output_file)
{
  std::cout << "initializing sampler... " << std::flush;

  arma::wall_clock timer;
  timer.tic();

  walkers = n_indep_runs;
  samples_per_walk = n_per_run;

  checkParameters();

  // Use the 2d (arma::Mat) to store sequences if only one sample per walk. Use
  // the 3d (arma::Cube) if more than 1 sampled.
  if (samples_per_walk == 1) {
    samples_2d = arma::Mat<int>(walkers, N, arma::fill::zeros);
    sample_stats = new SampleStats2D(&samples_2d, &(model));
  } else {
    samples_3d =
      arma::Cube<int>(samples_per_walk, N, walkers, arma::fill::zeros);
    sample_stats = new SampleStats3D(&samples_3d, &(model));
  }

  // Instantiate the PCG random number generator and uniform random
  // distribution.
  rng.seed(random_seed);

  std::cout << timer.toc() << " sec" << std::endl << std::endl;

  burn_in = burn_in_start;
  burn_between = burn_between_start;

  sampler = new Sampler(N, Q, &model);

  if ((samples_per_walk == 1) & update_burn_time) {
    std::cout << "setting burn time to... " << std::flush;
    timer.tic();
    estimateBurnTime();
    std::cout << burn_in << "... " << timer.toc() << " sec" << std::endl;
  }

  bool flag_mc = true;
  int resample_counter = 0;
  while (flag_mc) {
    std::cout << "sampling the model... " << std::flush;
    timer.tic();
    if (samples_per_walk > 1) {
      if (update_rule == "mh") {
        sampler->sampleSequences(&samples_3d,
                                 walkers,
                                 samples_per_walk,
                                 burn_in,
                                 burn_between,
                                 rng(),
                                 temperature);
      } else if (update_rule == "z-sqrt") {
        sampler->sampleSequencesZanellaSqrt(&samples_3d,
                                            walkers,
                                            samples_per_walk,
                                            burn_in,
                                            burn_between,
                                            rng(),
                                            temperature);
      } else if (update_rule == "z-barker") {
        sampler->sampleSequencesZanellaBarker(&samples_3d,
                                              walkers,
                                              samples_per_walk,
                                              burn_in,
                                              burn_between,
                                              rng(),
                                              temperature);
      } else {
        std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      if (update_rule == "mh") {
        sampler->sampleSequences(
          &samples_2d, walkers, burn_in, rng(), temperature);
      } else if (update_rule == "z-sqrt") {
        sampler->sampleSequencesZanellaSqrt(
          &samples_2d, walkers, burn_in, rng(), temperature);
      } else if (update_rule == "z-barker") {
        sampler->sampleSequencesZanellaBarker(
          &samples_2d, walkers, burn_in, rng(), temperature);
      } else {
        std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    std::cout << timer.toc() << " sec" << std::endl;

    std::cout << "computing sequence energies and correlations... "
              << std::flush;
    sample_stats->computeStatsExtra();
    timer.tic();
    std::cout << timer.toc() << " sec" << std::endl;

    if (update_burn_time & (samples_per_walk > 1)) {
      flag_mc = checkErgodicity();
      if (flag_mc) {
        if (resample_counter >= resample_max) {
          std::cout << "maximum number of resamplings (" << resample_counter
                    << ") reached. stopping..." << std::endl;
          flag_mc = false;
        } else {
          std::cout << "resampling..." << std::endl;
          resample_counter++;

          if (save_interim_samples) {
            std::cout << "writing temporary files" << std::endl;
            writeAASequences("temp_" + output_file);
            writeNumericalSequences("temp_" + output_file);
          }
        }
      }
    } else {
      flag_mc = false;
    }
  }

  int idx = output_file.find_last_of(".");
  std::string output_name = output_file.substr(0, idx);
  if (save_interim_samples) {
    if (deleteFile("temp_" + output_file) != 0)
      std::cerr << "temporary file deletion failed!" << std::endl;
    else if (deleteFile("temp_" + output_name + "_numerical.txt") != 0)
      std::cerr << "temporary file deletion failed!" << std::endl;
    else if (deleteFile("temp_" + output_name + "_energies.txt") != 0)
      std::cerr << "temporary file deletion failed!" << std::endl;
  }

  std::cout << "writing final sequences... " << std::flush;
  writeAASequences(output_file);
  writeNumericalSequences(output_file);
  std::cout << "done" << std::endl;

  return;
};
