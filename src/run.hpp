#ifndef BMDCA_RUN_HPP
#define BMDCA_RUN_HPP

#include "model.hpp"
#include "msa.hpp"
#include "msa_stats.hpp"
#include "pcg_random.hpp"
#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"

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
  int step_max = 500;
  int save_period = 40;
  bool save_best_steps = false;
  std::string stop_mode = "threshold";
  double stop_threshold = 0.00001;
  std::string train_mode = "original";
  bool cross_validate = false;
  int validation_seqs = 100;
  long int random_seed = 1;
  bool output_binary = true;

  // Sampler settings
  std::string update_rule = "mh";
  int burn_in_start = 10000;
  int burn_between_start = 100;
  bool update_burn_time = true;
  double adapt_up_time = 1.2;
  double adapt_down_time = 0.9;
  int step_importance_max = 1;
  double coherence_min = 0.9999;
  bool use_ss = false;
  int walkers = 1000;
  int samples_per_walk = 1;

  // Run settings
  int step = 1;
  int step_offset = 0;
  int burn_in = burn_in_start;
  int burn_between = burn_between_start;
  double train_err_tot_min = 1000;
  double validate_err_tot_min = 1000;

  std::string hyperparameter_file = "bmdca_params.conf";
  std::string run_log_file = "bmdca_run.log";

  bool checkErgodicity(void);

  // Buffers
  arma::Mat<double> run_buffer;
  void initializeRunLog();
  void writeRunLog(int = -1, int = 0, bool = false);

  // Sample data
  arma::Cube<int> samples_3d;
  arma::Mat<int> samples_2d;

  // Stats from original MSA
  MSA* msa_train = nullptr;
  MSA* msa_validate = nullptr;
  MSAStats* msa_train_stats = nullptr;
  MSAStats* msa_validate_stats = nullptr;

  // Training model
  Model* model;

  Sampler* sampler;

  SampleStats* sample_stats;

  // RNG
  pcg32 rng;
};

#endif
