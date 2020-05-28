#include "run.hpp"

#include <armadillo>
#include <cassert>
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

#define EPSILON 0.00000001

void
Sim::initializeParameters(void)
{
  // BM settings
  lambda_reg1 = 0.01;
  lambda_reg2 = 0.01;
  step_max = 2000;
  error_max = 0.00001;
  save_parameters = 20;
  random_seed = 1;
  use_reparametrization = true;

  // Learning rate settings
  epsilon_0_h = 0.01;
  epsilon_0_J = 0.001;
  adapt_up = 1.5;
  adapt_down = 0.6;
  min_step_h = 0.000001;
  max_step_h_N = 25; // divide by N later
  min_step_J = 0.0000001;
  max_step_J_N = 2.5; // divide by N later
  error_min_update = -1;

  // sampling time settings
  t_wait_0 = 10000;
  delta_t_0 = 100;
  check_ergo = true;
  adapt_up_time = 1.5;
  adapt_down_time = 0.600;

  output_binary = true;

  // importance sampling settings
  step_importance_max = 1;
  coherence_min = 0.9999;

  // mcmc settings
  M = 1000;            // importance sampling max iterations
  count_max = 10;      // number of independent MCMC runs
  init_sample = false; // flag to load first position for mcmc seqs
};

void
Sim::checkParameters(void)
{
  // If using stochastic sampling, set M to match the effective number of
  // sequences.
  if (use_ss) {
    M = 1;
    count_max = (int)(round(msa_stats.getEffectiveM()));
  }

  // Ensure that the set of ergodiciy checks is disabled if M=1
  if ((M == 1) && check_ergo) {
    check_ergo = false;
    std::cerr << "WARNING: disabling 'check_ergo' when M=1." << std::endl;
  }
}

void
Sim::writeParameters(std::string output_file)
{
  std::ofstream stream(output_file);

  // Header
  stream << "[bmDCA]" << std::endl;

  // BM settings
  stream << "lambda_reg1=" << lambda_reg1 << std::endl;
  stream << "lambda_reg2=" << lambda_reg2 << std::endl;
  stream << "step_max=" << step_max << std::endl;
  stream << "error_max=" << error_max << std::endl;
  stream << "save_parameters=" << save_parameters << std::endl;
  stream << "save_best_steps=" << save_best_steps << std::endl;
  stream << "random_seed=" << random_seed << std::endl;
  stream << "use_reparametrization=" << use_reparametrization << std::endl;

  // Learning rate settings
  stream << "epsilon_0_h=" << epsilon_0_h << std::endl;
  stream << "epsilon_0_J=" << epsilon_0_J << std::endl;
  stream << "adapt_up=" << adapt_up << std::endl;
  stream << "adapt_down=" << adapt_down << std::endl;
  stream << "min_step_h=" << min_step_h << std::endl;
  stream << "max_step_h_N=" << max_step_h_N << std::endl;
  stream << "min_step_J=" << min_step_J << std::endl;
  stream << "max_step_J_N=" << max_step_J_N << std::endl;
  stream << "error_min_update=" << error_min_update << std::endl;

  // sampling time settings
  stream << "t_wait_0=" << t_wait_0 << std::endl;
  stream << "delta_t_0=" << delta_t_0 << std::endl;
  stream << "check_ergo=" << check_ergo << std::endl;
  stream << "adapt_up_time=" << adapt_up_time << std::endl;
  stream << "adapt_down_time=" << adapt_down_time << std::endl;

  // importance sampling settings
  stream << "step_importance_max=" << step_importance_max << std::endl;
  stream << "coherence_min=" << coherence_min << std::endl;

  // mcmc settings
  stream << "use_ss=" << use_ss << std::endl;
  stream << "M=" << M << std::endl;
  stream << "count_max=" << count_max << std::endl;
  stream << "init_sample=" << init_sample << std::endl;
  stream << "init_sample_file=" << init_sample_file << std::endl;
  stream << "sampler=" << sampler << std::endl;
  stream << "use_pos_reg=" << use_pos_reg << std::endl;

  stream << "output_binary=" << output_binary << std::endl;
};

bool
Sim::compareParameters(std::string file_name)
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
        all_same = all_same & compareParameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return all_same;
};

void
Sim::loadParameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
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
        setParameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

bool
Sim::compareParameter(std::string key, std::string value)
{
  bool same = true;
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_reg1") {
    same = same & (lambda_reg1 == std::stod(value));
  } else if (key == "lambda_reg2") {
    same = same & (lambda_reg2 == std::stod(value));
  } else if (key == "step_max") {
  } else if (key == "error_max") {
  } else if (key == "save_parameters") {
  } else if (key == "save_best_steps") {
  } else if (key == "random_seed") {
    same = same & (random_seed == std::stoi(value));
  } else if (key == "use_reparametrization") {
    if (value.size() == 1) {
      same = same & (use_reparametrization == (std::stoi(value) == 1));
    } else {
      same = same & (use_reparametrization == (value == "true"));
    }
  } else if (key == "epsilon_0_h") {
    same = same & (epsilon_0_h == std::stod(value));
  } else if (key == "epsilon_0_J") {
    same = same & (epsilon_0_J == std::stod(value));
  } else if (key == "adapt_up") {
    same = same & (adapt_up == std::stod(value));
  } else if (key == "adapt_down") {
    same = same & (adapt_down == std::stod(value));
  } else if (key == "min_step_h") {
    same = same & (min_step_h == std::stod(value));
  } else if (key == "max_step_h_N") {
    same = same & (max_step_h_N == std::stod(value));
  } else if (key == "min_step_J") {
    same = same & (min_step_J == std::stod(value));
  } else if (key == "max_step_J_N") {
    same = same & (max_step_J_N == std::stod(value));
  } else if (key == "error_min_update") {
    same = same & (error_min_update == std::stod(value));
  } else if (key == "t_wait_0") {
    same = same & (t_wait_0 == std::stoi(value));
  } else if (key == "delta_t_0") {
    same = same & (delta_t_0 == std::stoi(value));
  } else if (key == "check_ergo") {
    if (value.size() == 1) {
      same = same & (check_ergo == (std::stoi(value) == 1));
    } else {
      same = same & (check_ergo == (value == "true"));
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
  } else if (key == "M") {
    same = same & (M == std::stoi(value));
  } else if (key == "count_max") {
    same = same & (count_max == std::stoi(value));
  } else if (key == "init_sample") {
    if (value.size() == 1) {
      same = same & (init_sample == (std::stoi(value) == 1));
    } else {
      same = same & (init_sample == (value == "true"));
    }
  } else if (key == "init_sample_file") {
    same = same & (init_sample_file == value);
  } else if (key == "sampler") {
    same = same & (sampler == value);
  } else if (key == "use_pos_reg") {
    if (value.size() == 1) {
      same = same & (use_pos_reg == (std::stoi(value) == 1));
    } else {
      same = same & (use_pos_reg == (value == "true"));
    }
  } else if (key == "output_binary") {
    if (value.size() == 1) {
      same = same & (output_binary == (std::stoi(value) == 1));
    } else {
      same = same & (output_binary == (value == "true"));
    }
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return same;
};

void
Sim::setParameter(std::string key, std::string value)
{
  // It's not possible to use switch blocks on strings because they are char*
  // arrays, not actual types.
  if (key == "lambda_reg1") {
    lambda_reg1 = std::stod(value);
  } else if (key == "lambda_reg2") {
    lambda_reg2 = std::stod(value);
  } else if (key == "step_max") {
    step_max = std::stoi(value);
  } else if (key == "error_max") {
    error_max = std::stod(value);
  } else if (key == "save_parameters") {
    save_parameters = std::stoi(value);
  } else if (key == "save_best_steps") {
    if (value.size() == 1) {
      save_best_steps = (std::stoi(value) == 1);
    } else {
      save_best_steps = (value == "true");
    }
  } else if (key == "random_seed") {
    random_seed = std::stoi(value);
  } else if (key == "use_reparametrization") {
    if (value.size() == 1) {
      use_reparametrization = (std::stoi(value) == 1);
    } else {
      use_reparametrization = (value == "true");
    }
  } else if (key == "epsilon_0_h") {
    epsilon_0_h = std::stod(value);
  } else if (key == "epsilon_0_J") {
    epsilon_0_J = std::stod(value);
  } else if (key == "adapt_up") {
    adapt_up = std::stod(value);
  } else if (key == "adapt_down") {
    adapt_down = std::stod(value);
  } else if (key == "min_step_h") {
    min_step_h = std::stod(value);
  } else if (key == "max_step_h_N") {
    max_step_h_N = std::stod(value);
  } else if (key == "min_step_J") {
    min_step_J = std::stod(value);
  } else if (key == "max_step_J_N") {
    max_step_J_N = std::stod(value);
  } else if (key == "error_min_update") {
    error_min_update = std::stod(value);
  } else if (key == "t_wait_0") {
    t_wait_0 = std::stoi(value);
  } else if (key == "delta_t_0") {
    delta_t_0 = std::stoi(value);
  } else if (key == "check_ergo") {
    if (value.size() == 1) {
      check_ergo = (std::stoi(value) == 1);
    } else {
      check_ergo = (value == "true");
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
  } else if (key == "M") {
    M = std::stoi(value);
  } else if (key == "count_max") {
    count_max = std::stoi(value);
  } else if (key == "init_sample") {
    if (value.size() == 1) {
      init_sample = (std::stoi(value) == 1);
    } else {
      init_sample = (value == "true");
    }
  } else if (key == "init_sample_file") {
    init_sample_file = value;
  } else if (key == "sampler") {
    sampler = value;
  } else if (key == "use_pos_reg") {
    if (value.size() == 1) {
      use_pos_reg = (std::stoi(value) == 1);
    } else {
      use_pos_reg = (value == "true");
    }
  } else if (key == "output_binary") {
    if (value.size() == 1) {
      output_binary = (std::stoi(value) == 1);
    } else {
      output_binary = (value == "true");
    }
  } else {
    std::cerr << "ERROR: unknown parameter '" << key << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

Sim::Sim(MSAStats msa_stats,
         std::string config_file,
         std::string dest_dir,
         bool force_restart)
  : msa_stats(msa_stats)
{

  initializeParameters();
  if (!config_file.empty()) {
    std::cout << "loading hyperparameters... " << std::flush;
    loadParameters(config_file);
    std::cout << "done." << std::endl;
  } else if ((!force_restart) &
             (checkFileExists(dest_dir + "/" + hyperparameter_file))) {
    std::cout << "loading previously used hyperparameters... " << std::flush;
    loadParameters(dest_dir + "/" + hyperparameter_file);
    std::cout << "done." << std::endl;
  }
  checkParameters();

  if (!dest_dir.empty()) {
    chdir(dest_dir.c_str());
  }

  if ((!force_restart) & (checkFileExists(hyperparameter_file))) {
    if (!compareParameters(hyperparameter_file)) {
      std::cerr << "ERROR: current and previous hyperparameters mismatched. "
                << std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      setStepOffset();
    }
  }

  if (step_offset == 0) {
    model = new Model(msa_stats, epsilon_0_h, epsilon_0_J);
  } else {
    if (output_binary) {
      model =
        new Model("parameters_h_" + std::to_string(step_offset) + ".bin",
                  "parameters_J_" + std::to_string(step_offset) + ".bin",
                  "parameters_h_" + std::to_string(step_offset - 1) + ".bin",
                  "parameters_J_" + std::to_string(step_offset - 1) + ".bin",
                  "gradients_h_" + std::to_string(step_offset) + ".bin",
                  "gradients_J_" + std::to_string(step_offset) + ".bin",
                  "gradients_h_" + std::to_string(step_offset - 1) + ".bin",
                  "gradients_J_" + std::to_string(step_offset - 1) + ".bin",
                  "learning_rates_h_" + std::to_string(step_offset) + ".bin",
                  "learning_rates_J_" + std::to_string(step_offset) + ".bin");
    } else {
      model =
        new Model("parameters_" + std::to_string(step_offset) + ".txt",
                  "parameters_" + std::to_string(step_offset - 1) + ".txt",
                  "gradients_" + std::to_string(step_offset) + ".txt",
                  "gradients_" + std::to_string(step_offset - 1) + ".txt",
                  "learning_rates_" + std::to_string(step_offset) + ".txt");
    }
  }
  mcmc = new MCMC(msa_stats.getN(), msa_stats.getQ());
};

void
Sim::clearFiles(std::string dest_dir)
{
  DIR* dp;
  struct dirent* dirp;

  std::vector<std::string> files;
  dp = opendir(dest_dir.c_str());

  while ((dirp = readdir(dp)) != NULL) {
    std::string fname = dirp->d_name;

    if (fname.find("parameters_"))
      files.push_back(fname);
    else if (fname.find("gradients_"))
      files.push_back(fname);
    else if (fname.find("bmdca_"))
      files.push_back(fname);
    else if (fname.find("learning_rates_"))
      files.push_back(fname);
    else if (fname.find("MC_energies_"))
      files.push_back(fname);
    else if (fname.find("MC_samples_"))
      files.push_back(fname);
    else if (fname.find("msa_numerical"))
      files.push_back(fname);
    else if (fname.find("overlap_"))
      files.push_back(fname);
    else if (fname.find("ergo_"))
      files.push_back(fname);
    else if (fname.find("stat_MC_"))
      files.push_back(fname);
    else if (fname.find("stat_align_"))
      files.push_back(fname);
    else if (fname.find("rel_ent_grad_"))
      files.push_back(fname);
  }
  closedir(dp);

  for (auto it = files.begin(); it != files.end(); ++it) {
    std::remove((*it).c_str());
  }
};

void
Sim::setStepOffset(void)
{
  DIR* dp;
  struct dirent* dirp;

  dp = opendir(".");

  int run_count = 0;
  std::ifstream stream(run_log_file);
  if (stream) {
    std::string line;
    std::getline(stream, line); // skip the file header
    while (std::getline(stream, line)) {
      run_count++;
    }
  }
  stream.close();

  std::vector<int> steps;
  std::vector<int> invalid_steps;
  while ((dirp = readdir(dp)) != NULL) {
    std::string fname = dirp->d_name;
    if (output_binary) {
      if (fname.find("parameters_h") == std::string::npos)
        continue;

      const std::regex re_param_h("parameters_h_([0-9]+)\\.bin");
      std::smatch match_param_h;

      int idx = -1;
      if (std::regex_match(fname, match_param_h, re_param_h)) {
        if (match_param_h.size() == 2) {
          idx = std::stoi(match_param_h[1].str());
        }
      }

      if (checkFileExists("parameters_J_" + std::to_string(idx) + ".bin") &
          checkFileExists("parameters_h_" + std::to_string(idx - 1) + ".bin") &
          checkFileExists("parameters_J_" + std::to_string(idx - 1) + ".bin") &
          checkFileExists("gradients_h_" + std::to_string(idx) + ".bin") &
          checkFileExists("gradients_J_" + std::to_string(idx) + ".bin") &
          checkFileExists("gradients_h_" + std::to_string(idx - 1) + ".bin") &
          checkFileExists("gradients_J_" + std::to_string(idx - 1) + ".bin") &
          checkFileExists("learning_rates_h_" + std::to_string(idx) + ".bin") &
          checkFileExists("learning_rates_J_" + std::to_string(idx) + ".bin")) {
        steps.push_back(idx);
      } else {
        if (idx > -1) {
          invalid_steps.push_back(idx);
        }
      }

    } else {
      if (fname.find("parameters") == std::string::npos)
        continue;

      const std::regex re_param("parameters_([0-9]+)\\.txt");
      std::smatch match_param;

      int idx = -1;
      if (std::regex_match(fname, match_param, re_param)) {
        if (match_param.size() == 2) {
          idx = std::stoi(match_param[1].str());
        }
      }

      if (checkFileExists("parameters_" + std::to_string(idx - 1) + ".txt") &
          checkFileExists("gradients_" + std::to_string(idx) + ".txt") &
          checkFileExists("gradients_" + std::to_string(idx - 1) + ".txt") &
          checkFileExists("learning_rates_" + std::to_string(idx) + ".txt")) {
        steps.push_back(idx);
      } else {
        invalid_steps.push_back(idx);
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
      if (output_binary) {
        std::string file;
        file = "parameters_h_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "parameters_J_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "parameters_h_" + std::to_string(*it_bad - 1) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "parameters_J_" + std::to_string(*it_bad - 1) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_h_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_J_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_h_" + std::to_string(*it_bad - 1) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_J_" + std::to_string(*it_bad - 1) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "learning_rates_h_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);

        file = "learning_rates_J_" + std::to_string(*it_bad) + ".bin";
        if (checkFileExists(file))
          deleteFile(file);
      } else {
        std::string file;
        file = "parameters_" + std::to_string(*it_bad) + ".txt";
        if (checkFileExists(file))
          deleteFile(file);

        file = "parameters_" + std::to_string(*it_bad - 1) + ".txt";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_" + std::to_string(*it_bad) + ".txt";
        if (checkFileExists(file))
          deleteFile(file);

        file = "gradients_" + std::to_string(*it_bad - 1) + ".txt";
        if (checkFileExists(file))
          deleteFile(file);

        file = "learning_rates_" + std::to_string(*it_bad) + ".txt";
        if (checkFileExists(file))
          deleteFile(file);
      }
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

Sim::~Sim(void)
{
  delete model;
  delete mcmc;
  delete mcmc_stats;
};

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

      t_wait = std::stoi(fields.at(2));
      delta_t = std::stoi(fields.at(3));

      // position of some variables depend on whether some flags are set
      if (check_ergo) {
        error_tot_min = std::stod(fields.at(17));
        prev_seed = std::stol(fields.at(18));
      } else {
        error_tot_min = std::stod(fields.at(7));
        prev_seed = std::stol(fields.at(8));
      }
      break;
    }
  }

  std::uniform_int_distribution<long int> dist(0, RAND_MAX - count_max);
  int counter = 1;
  while (dist(rng) != prev_seed) {
    if (counter > 1000 * step_max * step_importance_max) {
      std::cerr << "WARNING: cannot restore RNG state." << std::endl;
      break;
    }
    counter++;
  }
};

void
Sim::readInitialSample(int N, int Q)
{
  std::ifstream input_stream(init_sample_file);

  if (!input_stream) {
    std::cerr << "ERROR: cannot read '" << init_sample_file << "'."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string line;
  int aa;
  std::getline(input_stream, line);
  for (int n = 0; n < N; n++) {
    std::istringstream iss(line);
    iss >> aa;
    assert(aa < Q);
    initial_sample(n) = aa;
  }
  input_stream.close();
};

void
Sim::run(void)
{
  std::cout << "initializing run... " << std::flush;

  arma::wall_clock step_timer;

  arma::wall_clock timer;
  timer.tic();

  int N = model->N;
  int Q = model->Q;

  // Initialize sample data structure
  samples = arma::Cube<int>(M, N, count_max, arma::fill::zeros);

  mcmc->load(&(model->params));
  mcmc_stats = new MCMCStats(&samples, &(model->params));

  if (init_sample) {
    initial_sample = arma::Col<int>(N, arma::fill::zeros);
    readInitialSample(N, Q);
  }

  // Instantiate the PCG random number generator and unifrom random
  // distribution.
  rng.seed(random_seed);
  std::uniform_int_distribution<long int> dist(0, RAND_MAX - count_max);

  // Initialize the buffer.
  run_buffer = arma::Mat<double>(save_parameters, 20, arma::fill::zeros);
  int buffer_offset = 0;

  if (step_offset == 0) {
    initializeRunLog();
    t_wait = t_wait_0;
    delta_t = delta_t_0;
  } else if (step_offset >= step_max) {
    std::cout << "step " << step_max << " already reached... done."
              << std::endl;
    return;
  } else {
    restoreRunState();
    buffer_offset = (step_offset % save_parameters);
  }

  std::cout << timer.toc() << " sec" << std::endl << std::endl;

  // BM sampling loop
  for (step = 1 + step_offset; step <= step_max; step++) {
    step_timer.tic();
    std::cout << "Step: " << step << std::endl;

    // Sampling from MCMC (keep trying until correct properties found)
    bool flag_mc = true;
    long int seed;
    while (flag_mc) {
      // Draw from MCMC
      std::cout << "sampling model with mcmc... " << std::flush;
      timer.tic();
      seed = dist(rng);
      run_buffer((step - 1) % save_parameters, 18) = seed;
      if (sampler == "mh") {
        if (init_sample) {
          mcmc->sample_init(
            &samples, count_max, M, N, t_wait, delta_t, &initial_sample, seed);
        } else {
          mcmc->sample(&samples, count_max, M, N, t_wait, delta_t, seed);
        }
      } else if (sampler == "z-sqrt") {
        mcmc->sample_zanella(
          &samples, count_max, M, N, t_wait, delta_t, seed, "sqrt");
      } else if (sampler == "z-barker") {
        mcmc->sample_zanella(
          &samples, count_max, M, N, t_wait, delta_t, seed, "barker");
      } else {
        std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::cout << timer.toc() << " sec" << std::endl;

      std::cout << "updating mcmc with samples... " << std::flush;
      timer.tic();
      mcmc_stats->updateData(&samples, &(model->params));
      std::cout << timer.toc() << " sec" << std::endl;

      // Run checks and alter burn-in and wait times
      if (check_ergo) {
        std::cout << "computing sequence energies and correlations... "
                  << std::flush;
        timer.tic();
        mcmc_stats->computeEnergiesStats();
        mcmc_stats->computeCorrelations();
        std::cout << timer.toc() << " sec" << std::endl;

        std::vector<double> energy_stats = mcmc_stats->getEnergiesStats();
        std::vector<double> corr_stats = mcmc_stats->getCorrelationsStats();

        double auto_corr = corr_stats.at(2);
        double cross_corr = corr_stats.at(3);
        double check_corr = corr_stats.at(4);
        double cross_check_err = corr_stats.at(9);
        double auto_cross_err = corr_stats.at(8);

        double e_start = energy_stats.at(0);
        double e_start_sigma = energy_stats.at(1);
        double e_end = energy_stats.at(2);
        double e_end_sigma = energy_stats.at(3);
        double e_err = energy_stats.at(4);

        run_buffer((step - 1) % save_parameters, 4) = auto_corr;
        run_buffer((step - 1) % save_parameters, 5) = cross_corr;
        run_buffer((step - 1) % save_parameters, 6) = check_corr;
        run_buffer((step - 1) % save_parameters, 7) = auto_cross_err;
        run_buffer((step - 1) % save_parameters, 8) = cross_check_err;

        run_buffer((step - 1) % save_parameters, 9) = e_start;
        run_buffer((step - 1) % save_parameters, 10) = e_start_sigma;
        run_buffer((step - 1) % save_parameters, 11) = e_end;
        run_buffer((step - 1) % save_parameters, 12) = e_end_sigma;
        run_buffer((step - 1) % save_parameters, 13) = e_err;

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
          delta_t = (int)(round((double)delta_t * adapt_up_time));
          std::cout << "increasing wait time to " << delta_t << std::endl;
        } else if (flag_deltat_down) {
          delta_t = Max((int)(round((double)delta_t * adapt_down_time)), 1);
          std::cout << "decreasing wait time to " << delta_t << std::endl;
        }

        if (flag_twaiting_up) {
          t_wait = (int)(round((double)t_wait * adapt_up_time));
          std::cout << "increasing burn-in time to " << t_wait << std::endl;
        } else if (flag_twaiting_down) {
          t_wait = Max((int)(round((double)t_wait * adapt_down_time)), 1);
          std::cout << "decreasing burn-in time to " << t_wait << std::endl;
        }

        if (not flag_deltat_up and not flag_twaiting_up) {
          flag_mc = false;
        } else {
          std::cout << "resampling..." << std::endl;
        }
      } else {
        flag_mc = false;
      }
    }

    run_buffer((step - 1) % save_parameters, 0) = step;
    run_buffer((step - 1) % save_parameters, 1) = count_max;
    run_buffer((step - 1) % save_parameters, 2) = t_wait;
    run_buffer((step - 1) % save_parameters, 3) = delta_t;

    // Importance sampling loop
    int step_importance = 0;
    bool flag_coherence = true;
    bool new_min_found = false;
    while (step_importance < step_importance_max and flag_coherence == true) {
      step_importance++;
      if (step_importance > 1) {
        std::cout << "importance sampling step " << step_importance << "... "
                  << std::flush;
        timer.tic();
        mcmc_stats->computeSampleStatsImportance(&(model->params),
                                                 &(model->params_prev));
        std::cout << timer.toc() << " sec" << std::endl;

        double coherence = mcmc_stats->Z_ratio;
        if (coherence > coherence_min && 1.0 / coherence > coherence_min) {
          flag_coherence = true;
        } else {
          flag_coherence = false;
        }
      } else {
        std::cout << "computing mcmc 1p and 2p statistics... " << std::flush;
        timer.tic();
        mcmc_stats->computeSampleStats();
        std::cout << timer.toc() << " sec" << std::endl;
      }

      // Compute error reparametrization
      model->gradient_prev.h = model->gradient.h;
      model->gradient_prev.J = model->gradient.J;

      std::cout << "computing error and updating gradient... " << std::flush;
      timer.tic();
      computeErrorReparametrization();
      std::cout << timer.toc() << " sec" << std::endl;

      run_buffer((step - 1) % save_parameters, 14) = error_1p;
      run_buffer((step - 1) % save_parameters, 15) = error_2p;
      run_buffer((step - 1) % save_parameters, 16) = error_tot;
      if (error_tot < error_tot_min) {
        error_tot_min = error_tot;
        error_stat_tot_min = error_stat_tot;
        new_min_found = true;
      }
      run_buffer((step - 1) % save_parameters, 17) = error_tot_min;

      bool converged = false;
      if (error_tot < error_max) {
        converged = true;
      }

      if (converged) {
        run_buffer((step - 1) % save_parameters, 19) = step_timer.toc();
        std::cout << "converged! writing final results... " << std::flush;
        writeRunLog(step % save_parameters, buffer_offset);
        writeData(step);
        writeData(std::to_string(step) + "_final");
        std::cout << "done" << std::endl;
        return;
      }

      // Update learning rate
      std::cout << "update learning rate... " << std::flush;
      timer.tic();
      updateLearningRate();
      std::cout << timer.toc() << " sec" << std::endl;

      // Update parameters
      model->params_prev.h = model->params.h;
      model->params_prev.J = model->params.J;
      std::cout << "updating parameters... " << std::flush;
      timer.tic();
      updateReparameterization();
      std::cout << timer.toc() << " sec" << std::endl;

      run_buffer((step - 1) % save_parameters, 19) = step_timer.toc();

      // Save parameters
      if (step % save_parameters == 0 &&
          (step_importance == step_importance_max || flag_coherence == false)) {
        std::cout << "writing step " << step << "... " << std::flush;
        timer.tic();
        writeRunLog(step % save_parameters,  buffer_offset);
        buffer_offset = 0;
        writeData(step);
        std::cout << timer.toc() << " sec" << std::endl;
      } else if (new_min_found & (step > save_parameters) & save_best_steps) {
        new_min_found = false;
        std::cout << "close... writing step... " << std::flush;
        run_buffer((step - 1) % save_parameters, 19) = step_timer.toc();
        writeRunLog(step % save_parameters, buffer_offset, true);
        buffer_offset = step % save_parameters;
        writeData(step);
        std::cout << "done" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  if (step_offset != step_max) {
    std::cout << "writing final results... " << std::flush;
    if ((step_max % save_parameters) != 0) {
      writeRunLog(step_max % save_parameters);
    }
    writeData(step_max);
    writeData("final");
  } else {
    std::cout << "all " << step_offset << " steps already... " << std::flush;
  }

  std::cout << "done" << std::endl;
  return;
};

void
Sim::computeErrorReparametrization(void)
{
  double M_eff = msa_stats.getEffectiveM();
  int N = msa_stats.getN();
  int Q = msa_stats.getQ();

  error_1p = 0;
  error_2p = 0;
  error_stat_1p = 0;
  error_stat_2p = 0;

  double delta = 0;
  double delta_stat = 0;
  double deltamax_1 = 0;
  double deltamax_2 = 0;
  double rho = 0, beta = 0, den_beta = 0, num_beta = 0, num_rho = 0,
         den_stat = 0, den_mc = 0, c_mc_av = 0, c_stat_av = 0, rho_1p = 0,
         num_rho_1p = 0, den_stat_1p = 0, den_mc_1p = 0;

  double phi = 0;
  double lambda_h = lambda_reg1;
  double lambda_j = lambda_reg2;

  int count1 = 0;
  int count2 = 0;

  // Compute gradient
  if (use_reparametrization) {
    for (int i = 0; i < N; i++) {
      for (int aa = 0; aa < Q; aa++) {
        phi = 2 * model->params.h(aa, i);
        for (int j = 0; j < N; j++) {
          if (i < j) {
            for (int bb = 0; bb < Q; bb++) {
              phi +=
                model->params.J(i, j)(aa, bb) * msa_stats.frequency_1p(bb, j);
            }
          } else if (i > j) {
            for (int bb = 0; bb < Q; bb++) {
              phi +=
                model->params.J(j, i)(bb, aa) * msa_stats.frequency_1p(bb, j);
            }
          }
        }
        delta = mcmc_stats->frequency_1p(aa, i) -
                msa_stats.frequency_1p(aa, i) + lambda_h * 0.5 * phi;
        delta_stat =
          (mcmc_stats->frequency_1p(aa, i) - msa_stats.frequency_1p(aa, i)) /
          (pow(msa_stats.frequency_1p(aa, i) *
                   (1. - msa_stats.frequency_1p(aa, i)) / M_eff +
                 pow(mcmc_stats->frequency_1p_sigma(aa, i), 2) + EPSILON,
               0.5));
        error_1p += pow(delta, 2);
        error_stat_1p += pow(delta_stat, 2);
        if (pow(delta, 2) > pow(deltamax_1, 2))
          deltamax_1 = sqrt(pow(delta, 2));

        if (sqrt(pow(delta_stat, 2)) > error_min_update) {
          model->gradient.h(aa, i) = -delta;
          count1++;
        }
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      for (int aa = 0; aa < Q; aa++) {
        delta = mcmc_stats->frequency_1p(aa, i) -
                msa_stats.frequency_1p(aa, i) +
                lambda_h * model->params.h(aa, i);
        delta_stat =
          (mcmc_stats->frequency_1p(aa, i) - msa_stats.frequency_1p(aa, i)) /
          (pow(msa_stats.frequency_1p(aa, i) *
                   (1. - msa_stats.frequency_1p(aa, i)) / M_eff +
                 pow(mcmc_stats->frequency_1p_sigma(aa, i), 2) + EPSILON,
               0.5));
        error_1p += pow(delta, 2);
        error_stat_1p += pow(delta_stat, 2);
        if (pow(delta, 2) > pow(deltamax_1, 2))
          deltamax_1 = sqrt(pow(delta, 2));

        if (sqrt(pow(delta_stat, 2)) > error_min_update) {
          model->gradient.h(aa, i) = -delta;
          count1++;
        }
      }
    }
  }

  double error_c = 0;
  double c_mc = 0, c_stat = 0;

  if (use_reparametrization) {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            if (use_pos_reg) {
              delta = -(msa_stats.frequency_2p(i, j)(aa1, aa2) -
                        mcmc_stats->frequency_2p(i, j)(aa1, aa2) +
                        (mcmc_stats->frequency_1p(aa1, i) -
                         msa_stats.frequency_1p(aa1, i)) *
                          msa_stats.frequency_1p(aa2, j) +
                        (mcmc_stats->frequency_1p(aa2, j) -
                         msa_stats.frequency_1p(aa2, j)) *
                          msa_stats.frequency_1p(aa1, i) -
                        lambda_j * model->params.J(i, j)(aa1, aa2) * 1. /
                          (1. + fabs(msa_stats.rel_entropy_grad_1p(aa1, i))) /
                          (1. + fabs(msa_stats.rel_entropy_grad_1p(aa2, j))));
            } else {
              delta = -(msa_stats.frequency_2p(i, j)(aa1, aa2) -
                        mcmc_stats->frequency_2p(i, j)(aa1, aa2) +
                        (mcmc_stats->frequency_1p(aa1, i) -
                         msa_stats.frequency_1p(aa1, i)) *
                          msa_stats.frequency_1p(aa2, j) +
                        (mcmc_stats->frequency_1p(aa2, j) -
                         msa_stats.frequency_1p(aa2, j)) *
                          msa_stats.frequency_1p(aa1, i) -
                        lambda_j * model->params.J(i, j)(aa1, aa2));
            }
            delta_stat =
              (mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
               msa_stats.frequency_2p(i, j)(aa1, aa2)) /
              (pow(msa_stats.frequency_2p(i, j)(aa1, aa2) *
                       (1.0 - msa_stats.frequency_2p(i, j)(aa1, aa2)) / M_eff +
                     pow(mcmc_stats->frequency_2p(i, j)(aa1, aa2), 2) + EPSILON,
                   0.5));

            c_mc = mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
                   mcmc_stats->frequency_1p(aa1, i) *
                     mcmc_stats->frequency_1p(aa2, j);
            c_stat =
              msa_stats.frequency_2p(i, j)(aa1, aa2) -
              msa_stats.frequency_1p(aa1, i) * msa_stats.frequency_1p(aa2, j);
            c_mc_av += c_mc;
            c_stat_av += c_stat;
            error_c += pow(c_mc - c_stat, 2);
            error_2p += pow(delta, 2);
            error_stat_2p += pow(delta_stat, 2);

            if (pow(delta, 2) > pow(deltamax_2, 2)) {
              deltamax_2 = sqrt(pow(delta, 2));
            }
            if (sqrt(pow(delta_stat, 2)) > error_min_update) {
              model->gradient.J(i, j)(aa1, aa2) = -delta;
              count2++;
            }
          }
        }
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            if (use_pos_reg) {
              delta = -(msa_stats.frequency_2p(i, j)(aa1, aa2) -
                        mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
                        lambda_j * model->params.J(i, j)(aa1, aa2) * 1. /
                          (1. + fabs(msa_stats.rel_entropy_grad_1p(aa1, i))) /
                          (1. + fabs(msa_stats.rel_entropy_grad_1p(aa2, j))));
            } else {
              delta = -(msa_stats.frequency_2p(i, j)(aa1, aa2) -
                        mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
                        lambda_j * model->params.J(i, j)(aa1, aa2));
            }
            delta_stat =
              (mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
               msa_stats.frequency_2p(i, j)(aa1, aa2)) /
              (pow(msa_stats.frequency_2p(i, j)(aa1, aa2) *
                       (1.0 - msa_stats.frequency_2p(i, j)(aa1, aa2)) / M_eff +
                     pow(mcmc_stats->frequency_2p(i, j)(aa1, aa2), 2) + EPSILON,
                   0.5));

            c_mc = mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
                   mcmc_stats->frequency_1p(aa1, i) *
                     mcmc_stats->frequency_1p(aa2, j);
            c_stat =
              msa_stats.frequency_2p(i, j)(aa1, aa2) -
              msa_stats.frequency_1p(aa1, i) * msa_stats.frequency_1p(aa2, j);
            c_mc_av += c_mc;
            c_stat_av += c_stat;
            error_c += pow(c_mc - c_stat, 2);
            error_2p += pow(delta, 2);
            error_stat_2p += pow(delta_stat, 2);

            if (pow(delta, 2) > pow(deltamax_2, 2)) {
              deltamax_2 = sqrt(pow(delta, 2));
            }
            if (sqrt(pow(delta_stat, 2)) > error_min_update) {
              model->gradient.J(i, j)(aa1, aa2) = -delta;
              count2++;
            }
          }
        }
      }
    }
  }

  c_stat_av /= ((N * (N - 1) * Q * Q) / 2);
  c_mc_av /= ((N * (N - 1) * Q * Q) / 2);

  num_rho = num_beta = den_stat = den_mc = den_beta = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          c_mc =
            mcmc_stats->frequency_2p(i, j)(aa1, aa2) -
            mcmc_stats->frequency_1p(aa1, i) * mcmc_stats->frequency_1p(aa2, j);
          c_stat =
            msa_stats.frequency_2p(i, j)(aa1, aa2) -
            msa_stats.frequency_1p(aa1, i) * msa_stats.frequency_1p(aa2, j);
          num_rho += (c_mc - c_mc_av) * (c_stat - c_stat_av);
          num_beta += (c_mc) * (c_stat);
          den_stat += pow(c_stat - c_stat_av, 2);
          den_mc += pow(c_mc - c_mc_av, 2);
          den_beta += pow(c_stat, 2);
        }
      }
    }
  }

  num_rho_1p = den_stat_1p = den_mc_1p = 0;
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      num_rho_1p += (mcmc_stats->frequency_1p(aa, i) - 1.0 / Q) *
                    (msa_stats.frequency_1p(aa, i) - 1.0 / Q);
      den_stat_1p += pow(msa_stats.frequency_1p(aa, i) - 1.0 / Q, 2);
      den_mc_1p += pow(mcmc_stats->frequency_1p(aa, i) - 1.0 / Q, 2);
    }
  }

  beta = num_beta / den_beta;
  rho = num_rho / sqrt(den_mc * den_stat);
  rho_1p = num_rho_1p / sqrt(den_mc_1p * den_stat_1p);

  error_1p = sqrt(error_1p / (N * Q));
  error_2p = sqrt(error_2p / ((N * (N - 1) * Q * Q) / 2));

  error_stat_1p = sqrt(error_stat_1p / (N * Q));
  error_stat_2p = sqrt(error_stat_2p / (N * (N - 1) * Q * Q) / 2);

  error_c = sqrt(error_c / (N * (N - 1) * Q * Q) / 2);

  error_tot = error_1p + error_2p;
  error_stat_tot = error_stat_1p + error_stat_2p;

  return;
};

void
Sim::updateLearningRate(void)
{
  int N = msa_stats.getN();
  int Q = msa_stats.getQ();
  double max_step_h = max_step_h_N / N;
  double max_step_J = max_step_J_N / N;

  double alfa = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          alfa =
            model->gradient.J(i, j)(a, b) * model->gradient_prev.J(i, j)(a, b);
          alfa =
            Theta(alfa) * adapt_up + Theta(-alfa) * adapt_down + Delta(alfa);

          model->learning_rates.J(i, j)(a, b) =
            Min(max_step_J,
                Max(min_step_J, alfa * model->learning_rates.J(i, j)(a, b)));
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      alfa = model->gradient.h(a, i) * model->gradient_prev.h(a, i);
      alfa = Theta(alfa) * adapt_up + Theta(-alfa) * adapt_down + Delta(alfa);
      model->learning_rates.h(a, i) =
        Min(max_step_h, Max(min_step_h, alfa * model->learning_rates.h(a, i)));
    }
  }
};

void
Sim::updateReparameterization(void)
{
  int N = msa_stats.getN();
  int Q = msa_stats.getQ();

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          model->params.J(i, j)(a, b) +=
            model->learning_rates.J(i, j)(a, b) * model->gradient.J(i, j)(a, b);
        }
      }
    }
  }

  if (use_reparametrization) {
    arma::Mat<double> Dh = arma::Mat<double>(Q, N, arma::fill::zeros);
    for (int i = 0; i < N; i++) {
      for (int a = 0; a < Q; a++) {
        for (int j = 0; j < N; j++) {
          if (i < j) {
            for (int b = 0; b < Q; b++) {
              Dh(a, i) += msa_stats.frequency_1p(b, j) *
                          model->learning_rates.J(i, j)(a, b) *
                          model->gradient.J(i, j)(a, b);
            }
          }
          if (i > j) {
            for (int b = 0; b < Q; b++) {
              Dh(a, i) += msa_stats.frequency_1p(b, j) *
                          model->learning_rates.J(j, i)(b, a) *
                          model->gradient.J(j, i)(b, a);
            }
          }
        }
      }
    }

    for (int i = 0; i < N; i++) {
      for (int a = 0; a < Q; a++) {
        model->params.h(a, i) +=
          model->learning_rates.h(a, i) * model->gradient.h(a, i) +
          0.5 * Dh(a, i);
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      for (int a = 0; a < Q; a++) {
        model->params.h(a, i) +=
          model->learning_rates.h(a, i) * model->gradient.h(a, i);
      }
    }
  }
};

void
Sim::writeData(int step)
{
  if (output_binary) {
    model->writeParamsPrevious(
      "parameters_h_" + std::to_string(step - 1) + ".bin",
      "parameters_J_" + std::to_string(step - 1) + ".bin");
    model->writeGradientPrevious(
      "gradients_h_" + std::to_string(step - 1) + ".bin",
      "gradients_J_" + std::to_string(step - 1) + ".bin");
    model->writeParams("parameters_h_" + std::to_string(step) + ".bin",
                       "parameters_J_" + std::to_string(step) + ".bin");
    model->writeGradient("gradients_h_" + std::to_string(step) + ".bin",
                         "gradients_J_" + std::to_string(step) + ".bin");
    model->writeLearningRates(
      "learning_rates_h_" + std::to_string(step) + ".bin",
      "learning_rates_J_" + std::to_string(step) + ".bin");

    mcmc_stats->writeFrequency1p("stat_MC_1p_" + std::to_string(step) + ".bin",
                                 "stat_MC_1p_sigma_" + std::to_string(step) +
                                   ".bin");
    mcmc_stats->writeFrequency2p("stat_MC_2p_" + std::to_string(step) + ".bin",
                                 "stat_MC_2p_sigma_" + std::to_string(step) +
                                   ".bin");
  } else {
    model->writeParamsPreviousAscii("parameters_" + std::to_string(step - 1) +
                                    ".txt");
    model->writeGradientPreviousAscii("gradients_" + std::to_string(step - 1) +
                                      ".txt");

    model->writeParamsAscii("parameters_" + std::to_string(step) + ".txt");
    model->writeGradientAscii("gradients_" + std::to_string(step) + ".txt");
    model->writeLearningRatesAscii("learning_rates_" + std::to_string(step) +
                                   ".txt");

    mcmc_stats->writeFrequency1pAscii(
      "stat_MC_1p_" + std::to_string(step) + ".txt",
      "stat_MC_1p_sigma_" + std::to_string(step) + ".txt");
    mcmc_stats->writeFrequency2pAscii(
      "stat_MC_2p_" + std::to_string(step) + ".txt",
      "stat_MC_2p_sigma_" + std::to_string(step) + ".txt");
  }
  mcmc_stats->writeSamples("MC_samples_" + std::to_string(step) + ".txt");
  mcmc_stats->writeSampleEnergies("MC_energies_" + std::to_string(step) +
                                  ".txt");

  if (check_ergo) {
    mcmc_stats->writeCorrelationsStats(
      "overlap_" + std::to_string(step) + ".txt",
      "overlap_inf_" + std::to_string(step) + ".txt",
      "ergo_" + std::to_string(step) + ".txt");
  }
};

void
Sim::writeData(std::string id)
{
  if (output_binary) {
    model->writeParams("parameters_h_" + id + ".bin",
                       "parameters_J_" + id + ".bin");
    model->writeGradient("gradients_h_" + id + ".bin",
                         "gradients_J_" + id + ".bin");
    model->writeLearningRates("learning_rates_h_" + id + ".bin",
                              "learning_rates_J_" + id + ".bin");

    mcmc_stats->writeFrequency1p("stat_MC_1p_" + id + ".bin",
                                 "stat_MC_1p_sigma_" + id + ".bin");
    mcmc_stats->writeFrequency2p("stat_MC_2p_" + id + ".bin",
                                 "stat_MC_2p_sigma_" + id + ".bin");
  } else {
    model->writeParamsAscii("parameters_" + id + ".txt");
    model->writeGradientAscii("gradients_" + id + ".txt");
    model->writeLearningRatesAscii("learning_rates_" + id + ".txt");

    mcmc_stats->writeFrequency1pAscii("stat_MC_1p_" + id + ".txt",
                                      "stat_MC_1p_sigma_" + id + ".txt");
    mcmc_stats->writeFrequency2pAscii("stat_MC_2p_" + id + ".txt",
                                      "stat_MC_2p_sigma_" + id + ".txt");
  }
  mcmc_stats->writeSamples("MC_samples_" + id + ".txt");
  mcmc_stats->writeSampleEnergies("MC_energies_" + id + ".txt");

  if (check_ergo) {
    mcmc_stats->writeCorrelationsStats("overlap_" + id + ".txt",
                                       "overlap_inf_" + id + ".txt",
                                       "ergo_" + id + ".txt");
  }
};

void
Sim::initializeRunLog()
{
  std::ofstream stream{ run_log_file, std::ios_base::out };
  stream << "step"
         << "\t"
         << "reps"
         << "\t"
         << "burn-in"
         << "\t"
         << "burn-between"
         << "\t";
  if (check_ergo) {
    stream << "auto-corr"
           << "\t"
           << "cross-corr"
           << "\t"
           << "check-corr"
           << "\t"
           << "auto-cross-err"
           << "\t"
           << "cross-check-err"
           << "\t"
           << "energy-start-avg"
           << "\t"
           << "sigma-energy-start-sigma"
           << "\t"
           << "energy-end-avg"
           << "\t"
           << "energy-end-sigma"
           << "\t"
           << "energy-err"
           << "\t";
  }
  stream << "error-h"
         << "\t"
         << "error-J"
         << "\t"
         << "error-tot"
         << "\t"
         << "error-tot-min"
         << "\t"
         << "seed"
         << "\t"
         << "step-time" << std::endl;
  stream.close();
};

void
Sim::writeRunLog(int current_step, int offset, bool keep)
{
  std::ofstream stream{ run_log_file, std::ios_base::app };

  int n_entries;
  if (current_step == 0) {
    n_entries = save_parameters;
  } else {
    n_entries = current_step;
  }
  for (int i = 0 + offset; i < n_entries; i++) {
    stream << (int)run_buffer(i, 0) << "\t";
    stream << (int)run_buffer(i, 1) << "\t";
    stream << (int)run_buffer(i, 2) << "\t";
    stream << (int)run_buffer(i, 3) << "\t";
    if (check_ergo) {
      stream << run_buffer(i, 4) << "\t";
      stream << run_buffer(i, 5) << "\t";
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
    stream << (long int)run_buffer(i, 18) << "\t";
    stream << run_buffer(i, 19) << std::endl;
  }
  if (!keep) {
    run_buffer.zeros();
  }
  stream.close();
};
