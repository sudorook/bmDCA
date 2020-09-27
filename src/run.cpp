#include "run.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>
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

Sim::Sim(MSA* msa_train,
         MSA* msa_validate,
         std::string config_file,
         std::string dest_dir,
         bool force_restart)
  : msa_train(msa_train)
  , msa_validate(msa_validate)
{
  std::cout << "----------------------- bmDCA -----------------------"
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

  // Load model hyperparameters.
  if (train_mode == "adam") {
    model = new Adam();
  } else if (train_mode == "adamw") {
    model = new AdamW();
  } else if (train_mode == "original") {
    model = new Original();
  } else if (train_mode == "radam") {
    model = new RAdam();
  } else if (train_mode == "reparametrization") {
    model = new Reparam();
  } else if (train_mode == "sgdm") {
    model = new SGDM();
  } else {
    std::cerr << "ERROR: unrecognised training model '" << train_mode
              << "' given. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }

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

void
Sim::initialize(void)
{

  if (!msa_validate) {
    if (cross_validate) {
      std::vector<MSA*> tmp =
        msa_train->partitionAlignment(validation_seqs, random_seed);
      msa_train = tmp[0];
      msa_validate = tmp[1];
      msa_validate->writeSequenceWeights("msa_weights_validate.txt");
      msa_validate->writeMatrix("msa_numerical_validate.txt");
    }
  }

  msa_train_energies = arma::Col<double>(msa_train->M, arma::fill::zeros);
  msa_train->writeMatrix("msa_numerical.txt");
  msa_train->writeSequenceWeights("msa_weights.txt");

  // Compute stats for training MSA
  msa_train_stats = new MSAStats(msa_train, true);
  msa_train_stats->writeFrequency1p("msa_stat_1p.bin");
  msa_train_stats->writeFrequency2p("msa_stat_2p.bin");
  std::cout << std::endl;

  if (msa_validate) {
    msa_validate_energies =
      arma::Col<double>(msa_validate->M, arma::fill::zeros);
    msa_validate->writeMatrix("msa_numerical_validate.txt");
    msa_validate->writeSequenceWeights("msa_weights_validate.txt");

    // Compute stats for validation MSA
    msa_validate_stats = new MSAStats(msa_validate, true);
    msa_validate_stats->writeFrequency1p("msa_stat_1p_validate.bin");
    msa_validate_stats->writeFrequency2p("msa_stat_2p_validate.bin");
    std::cout << std::endl;

    if (!cross_validate) {
      std::cout
        << "NOTE: Ignoring cross_validate=false due to provided validation MSA."
        << std::endl;
    }
  }

  // If using stochastic sampling, set M to match the effective number of
  // sequences.
  if (use_ss) {
    samples_per_walk = 1;
    burn_between_start = 0;
    burn_between = 0;
    walkers = (int)(round(msa_train_stats->getEffectiveM()));
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
    sample_stats =
      new SampleStats3D(&samples_3d, &(model->params), &(model->params_prev));
    sample_stats->setMixingTime(burn_between);
  } else {
    samples_2d = arma::Mat<int>(walkers, N, arma::fill::zeros);
    sample_stats =
      new SampleStats2D(&samples_2d, &(model->params), &(model->params_prev));
  }

  model->setSampleStats(sample_stats);
  model->setStep(step_offset);

  // Initialize the sampler
  sampler = new Sampler(N, Q, &(model->params));
};

Sim::~Sim(void)
{
  delete msa_train_stats;
  delete msa_validate_stats;
  delete model;
  delete sampler;
  delete sample_stats;
};

void
Sim::loadHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
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

bool
Sim::compareHyperparameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
  bool all_same = true;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
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

void
Sim::setHyperparameter(std::string key, std::string value)
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

void
Sim::checkHyperparameters(void)
{
  if ((validation_seqs < 1) & cross_validate) {
    std::cerr << "ERROR: Incompatible number of sequecnes for cross-validation."
              << std::endl;
  }
};

void
Sim::writeHyperparameters(std::string output_file, bool append)
{
  // std::ofstream stream(output_file);
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

void
Sim::setStepOffset(void)
{
  DIR* dp;
  struct dirent* dirp;

  dp = opendir(".");

  std::vector<int> steps;
  std::vector<int> invalid_steps;
  while ((dirp = readdir(dp)) != NULL) {
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

// void
// Sim::clearFiles(std::string dest_dir)
// {
//   DIR* dp;
//   struct dirent* dirp;
//
//   std::vector<std::string> files;
//   dp = opendir(dest_dir.c_str());
//
//   while ((dirp = readdir(dp)) != NULL) {
//     std::string fname = dirp->d_name;
//
//     if (fname.find("parameters_"))
//       files.push_back(fname);
//     // else if (fname.find("parameters_avg_"))
//     //   files.push_back(fname);
//     else if (fname.find("gradients_"))
//       files.push_back(fname);
//     else if (fname.find("bmdca_"))
//       files.push_back(fname);
//     else if (fname.find("moment1_"))
//       files.push_back(fname);
//     else if (fname.find("moment2_"))
//       files.push_back(fname);
//     else if (fname.find("MC_energies_"))
//       files.push_back(fname);
//     else if (fname.find("MC_samples_"))
//       files.push_back(fname);
//     else if (fname.find("msa_numerical"))
//       files.push_back(fname);
//     else if (fname.find("overlap_"))
//       files.push_back(fname);
//     else if (fname.find("ergo_"))
//       files.push_back(fname);
//     else if (fname.find("stat_MC_"))
//       files.push_back(fname);
//     else if (fname.find("stat_align_"))
//       files.push_back(fname);
//     else if (fname.find("rel_ent_grad_"))
//       files.push_back(fname);
//   }
//   closedir(dp);
//
//   for (auto it = files.begin(); it != files.end(); ++it) {
//     std::remove((*it).c_str());
//   }
// };

void
Sim::restoreRunState(void)
{
  long int prev_seed = -1;
  std::ifstream stream(run_log_file);
  std::string line;
  std::getline(stream, line);
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
          prev_seed = std::stol(fields.at(23));
        } else {
          prev_seed = std::stol(fields.at(18));
        }
      } else {
        if (cross_validate) {
          train_err_tot_min = std::stod(fields.at(8));
          validate_err_tot_min = std::stod(fields.at(12));
          prev_seed = std::stol(fields.at(14));
        } else {
          train_err_tot_min = std::stod(fields.at(8));
          prev_seed = std::stol(fields.at(9));
        }
      }
      break;
    }
  }

  std::uniform_int_distribution<long int> dist(0, RAND_MAX - step_max);
  int counter = 1;
  while (dist(rng) != prev_seed) {
    if (counter > 1000000 * step_max) {
      std::cerr << "WARNING: cannot restore RNG state." << std::endl;
      break;
    }
    counter++;
  }
};

void
Sim::run(void)
{
  std::cout << "initializing run... " << std::flush;

  arma::wall_clock step_timer;

  arma::wall_clock timer;
  timer.tic();

  // Instantiate the PCG random number generator and unifrom random
  // distribution.
  rng.seed(random_seed);
  std::uniform_int_distribution<long int> dist(0, RAND_MAX - step_max);

  // Initialize the buffer.
  run_buffer = arma::Mat<double>(save_period, 25, arma::fill::zeros);
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
    double error_stddev = arma::stddev(msa_train_stats->msa_rms, 1);
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

  int N = msa_train_stats->getN();
  int Q = msa_train_stats->getQ();

  // BM sampling loop
  for (step = 1 + step_offset; step <= step_max; step++) {
    step_timer.tic();
    std::cout << "Step: " << step << std::endl;
    model->setStep(step);

    // Sampling from MCMC (keep trying until correct properties found)
    bool flag_mc = true;
    long int seed;
    while (flag_mc) {
      if (update_burn_time & (samples_per_walk == 1)) {
        std::cout << "setting burn time to... " << std::flush;
        timer.tic();
        bool flag_burn = true;
        while (flag_burn) {
          double burn_reps = 24;
          double burn_count = 4;
          arma::Mat<double> energy_burn =
            arma::Mat<double>(burn_count, burn_reps, arma::fill::zeros);

          sampler->sampleEnergies(
            &energy_burn, burn_reps, burn_count, burn_in, burn_in, dist(rng));

          double e_start = arma::mean(energy_burn.row(0));
          double e_start_sigma = arma::stddev(energy_burn.row(0), 1);
          double e_end = (arma::mean(energy_burn.row(burn_count - 1)) +
                          arma::mean(energy_burn.row(burn_count - 2))) /
                         2;
          double e_end_sigma =
            sqrt(pow(arma::stddev(energy_burn.row(burn_count - 1), 1), 2) +
                 pow(arma::stddev(energy_burn.row(burn_count - 2), 1), 2));
          double e_err =
            sqrt((pow(e_start_sigma, 2) / burn_reps + pow(e_end_sigma, 2)) / 2 /
                 burn_reps);

          bool flag_twaiting_up = true;
          bool flag_twaiting_down = true;
          if (e_start - e_end <= 2 * e_err) {
            flag_twaiting_up = false;
          }
          if (e_start - e_end >= -2 * e_err) {
            flag_twaiting_down = false;
          }
          if (flag_twaiting_up) {
            burn_in = (int)(round((double)burn_in * adapt_up_time));
          } else if (flag_twaiting_down) {
            burn_in = Max((int)(round((double)burn_in * adapt_down_time)), 1);
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
      seed = dist(rng);
      run_buffer((step - 1) % save_period, 23) = seed;
      if (samples_per_walk > 1) {
        if (update_rule == "mh") {
          sampler->sampleSequences(&samples_3d,
                                   walkers,
                                   samples_per_walk,
                                   burn_in,
                                   burn_between,
                                   seed);
        } else if (update_rule == "z-sqrt") {
          sampler->sampleSequencesZanella(&samples_3d,
                                          walkers,
                                          samples_per_walk,
                                          burn_in,
                                          burn_between,
                                          seed,
                                          "sqrt");
        } else if (update_rule == "z-barker") {
          sampler->sampleSequencesZanella(&samples_3d,
                                          walkers,
                                          samples_per_walk,
                                          burn_in,
                                          burn_between,
                                          seed,
                                          "barker");
        } else {
          std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        if (update_rule == "mh") {
          sampler->sampleSequences(&samples_2d, walkers, burn_in, seed);
        } else if (update_rule == "z-sqrt") {
          sampler->sampleSequencesZanella(
            &samples_2d, walkers, burn_in, seed, "sqrt");
        } else if (update_rule == "z-barker") {
          sampler->sampleSequencesZanella(
            &samples_2d, walkers, burn_in, seed, "barker");
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
          std::cout << "resampling..." << std::endl;
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

    computeMSAEnergies(&msa_train_energies, msa_train, &(model->params));
    if (msa_validate) {
      computeMSAEnergies(
        &msa_validate_energies, msa_validate, &(model->params));
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

    double e_start = 0;
    double e_start_sigma = 0;
    double e_end = 0;
    double e_end_sigma = 0;
    double e_err = 0;

    double total_corr = 0;
    double auto_corr = 0;
    double cross_corr = 0;
    double auto_cross_err = 0;

    // Get sample statistics
    if (samples_per_walk > 1) {
      arma::Col<double> stats = sample_stats->getStats();

      e_start = stats.at(0);
      e_start_sigma = stats.at(1);
      e_end = stats.at(2);
      e_end_sigma = stats.at(3);
      e_err = stats.at(4);

      total_corr = stats.at(5);
      auto_corr = stats.at(7);
      cross_corr = stats.at(8);
      auto_cross_err = stats.at(13);

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
      run_buffer((step - 1) % save_period, 24) = step_timer.toc();
      std::cout << "converged! writing final results... " << std::flush;
      writeRunLog(step % save_period, buffer_offset);
      writeStep(step);
      // writeData(std::to_string(step) + "_final");
      std::cout << "done" << std::endl;
      return;
    }

    run_buffer((step - 1) % save_period, 24) = step_timer.toc();

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
      run_buffer((step - 1) % save_period, 24) = step_timer.toc();
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

bool
Sim::checkErgodicity(void)
{
  arma::Col<double> stats = sample_stats->getStats();

  double e_start = stats.at(0);
  double e_end = stats.at(2);
  double e_err = stats.at(4);

  double auto_corr = stats.at(7);
  double cross_corr = stats.at(8);
  double check_corr = stats.at(9);
  double cross_check_err = stats.at(14);
  double auto_cross_err = stats.at(13);

  bool flag_deltat_up = true;
  bool flag_deltat_down = true;
  bool flag_twaiting_up = true;
  bool flag_twaiting_down = true;

  if (check_corr - cross_corr <= cross_check_err) {
    flag_deltat_up = false;
  }
  if (auto_corr - cross_corr >= auto_cross_err) {
    flag_deltat_down = false;
  }

  if (e_start - e_end <= 2 * e_err) {
    flag_twaiting_up = false;
  }
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

  bool flag_mc = true;
  if (not flag_deltat_up and not flag_twaiting_up) {
    flag_mc = false;
  }
  return flag_mc;
};

void
Sim::computeMSAEnergies(arma::Col<double>* energies,
                        MSA* msa,
                        potts_model* params)
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

void
Sim::writeStep(int step)
{
  model->writeStep(step, output_binary);
  sample_stats->writeStep(step, output_binary);
  writeMSAEnergies(step);
};

void
Sim::writeData(std::string id)
{
  model->writeData(id, output_binary);
  sample_stats->writeData(id, output_binary);
  writeMSAEnergies(id);
};

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

void
Sim::writeMSAEnergies(std::string id)
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

void
Sim::initializeRunLog()
{
  std::ofstream stream{ run_log_file, std::ios_base::out };
  stream << "step"
         << "\t"
         << "walkers"
         << "\t"
         << "samples-per-walk"
         << "\t"
         << "burn-in"
         << "\t";
  if (samples_per_walk > 1) {
    stream << "burn-between"
           << "\t";
  }
  stream << "total-corr"
         << "\t";
  if (samples_per_walk > 1) {
    stream << "auto-corr"
           << "\t"
           << "cross-corr"
           << "\t"
           << "auto-cross-err"
           << "\t"
           << "energy-start-avg"
           << "\t"
           << "energy-start-sigma"
           << "\t"
           << "energy-end-avg"
           << "\t"
           << "energy-end-sigma"
           << "\t"
           << "energy-err"
           << "\t";
  }
  stream << "train-err-1p"
         << "\t"
         << "train-err-2p"
         << "\t"
         << "train-err-tot"
         << "\t"
         << "train-err-tot-min"
         << "\t";
  if (cross_validate) {
    stream << "validate-err-1p"
           << "\t"
           << "validate-err-2p"
           << "\t"
           << "validate-err-tot"
           << "\t"
           << "validate-err-tot-min"
           << "\t"
           << "diff-avg-energy"
           << "\t";
  }
  stream << "seed"
         << "\t"
         << "duration" << std::endl;
  stream.close();
};

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
    stream << (int)run_buffer(i, 0) << "\t";
    stream << (int)run_buffer(i, 1) << "\t";
    stream << (int)run_buffer(i, 2) << "\t";
    stream << (int)run_buffer(i, 3) << "\t";
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
    stream << (long int)run_buffer(i, 23) << "\t";
    stream << run_buffer(i, 24) << std::endl;
  }
  if (!keep) {
    run_buffer.zeros();
  }
  stream.close();
};
