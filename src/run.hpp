#ifndef BMDCA_RUN_HPP
#define BMDCA_RUN_HPP

#include "mcmc.hpp"
#include "mcmc_stats.hpp"
#include "model.hpp"
#include "msa.hpp"
#include "msa_stats.hpp"
#include "pcg_random.hpp"
#include "utils.hpp"

class Sim
{
public:
  Sim(MSAStats, std::string, std::string, bool);
  ~Sim(void);
  void run(void);
  void loadParameters(std::string);
  void writeParameters(std::string);

private:
  void initializeParameters(void);
  bool compareParameters(std::string);
  void checkParameters(void);
  void setStepOffset(void);
  void restoreRunState(void);
  void readInitialSample(int, int);
  void computeErrorReparametrization(void);
  void updateLearningRate(void);
  void updateReparameterization(void);
  void writeData(std::string);
  void writeData(int);
  void clearFiles(std::string);

  // BM settings
  double lambda_reg1;  // L2 regularization strength for 1p statistics (fields)
  double lambda_reg2;  // L2 regularization strength for 2p statistics (cpling)
  int step_max;        // max number of BM steps
  double error_max;    // exit error
  int save_parameters; // multiple of iterations at which to save parameters
  int save_best_steps; // multiple of iterations at which to save parameters
  int random_seed;
  bool use_reparametrization = true;

  // Learning parameters
  double epsilon_0_h;      // starting learning rate for fields
  double epsilon_0_J;      // starting learning rate for couplings
  double adapt_up;         // positive adaptive step for learning rate
  double adapt_down;       // negative sdaptive step for learning rate
  double min_step_h;       // min learning rate for fields
  double max_step_h_N;     // maximum learning rate for fields
  double min_step_J;       // minimum learning rate for couplings
  double max_step_J_N;     // maximum learning rate for couplings (to be
                           // divided by N)
  double error_min_update; // minimal number of standard deviation s of z
                           // variable for having parameter update (if
                           // negative or zero all parameters are updated)

  // Sampling times
  int t_wait_0;           // staring thermalization time for MCMC
  int delta_t_0;          // starging samplign time for MCMC
  bool check_ergo;        // flag to check if MC is well thermalized and
                          // decorrelated
  double adapt_up_time;   // negative adaptive step for sampling/
                          // thermalization times
  double adapt_down_time; // positive adaptive step for sampling/
                          // thermalization times

  int t_wait;
  int delta_t;

  // Importance sampling settings
  int step_importance_max; // importance sampling maximum iterations
  double coherence_min;    // coherence importance sampling

  // MCMC settings
  bool use_ss = false;          // flag to use stochastic sampling mode
  int step;                     // current step number
  int step_offset = 0;
  int M;                        // number of samples for each MCMC run
  int count_max;                // number of independent MCMC runs
  bool init_sample = false;     // flag for loading the first positions when
                                // initializing the mcmc from a file
  std::string init_sample_file; // name of file with mcmc initial sample
  std::string sampler = "mh";   // MC sampler type ('mh' or 'z')
  bool use_pos_reg = false;     // enable for position-specific regularizetion

  bool output_binary = true;

  double error_1p;
  double error_2p;
  double error_tot;
  double error_tot_min = 1;
  double error_stat_1p;
  double error_stat_2p;
  double error_stat_tot;
  double error_stat_tot_min = 1;

  std::string hyperparameter_file = "bmdca_params.conf";
  std::string run_log_file = "bmdca_run.log";

  // Key-value wrapper for loading parameters from a file.
  void setParameter(std::string, std::string);
  bool compareParameter(std::string, std::string);

  // Buffers
  arma::Mat<double> run_buffer;
  void initializeRunLog();
  void writeRunLog(int = -1, int = 0, bool = false);

  // Sample data
  arma::Cube<int> samples;
  arma::Col<int> initial_sample;

  // Stats from original MSA
  MSAStats msa_stats;

  Model* model;

  // MCMC
  MCMC* mcmc;

  // Stats from MCMC samples
  MCMCStats* mcmc_stats;

  // RNG
  pcg32 rng;
};

#endif
