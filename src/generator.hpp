#include "pcg_random.hpp"
#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"

class Generator
{
public:
  Generator(potts_model, int, int, std::string);
  ~Generator(void);
  void run(int, int, std::string);
  void writeAASequences(std::string);
  void writeNumericalSequences(std::string);

private:
  int N;                // number of positions
  int Q;                // number of amino acids
  int samples_per_walk; // number of independent sampling runs
  int walkers;          // number of sequences sampled from independent runs

  int resample_max = 20;
  long int random_seed = 1;
  double adapt_up_time = 1.5;
  double adapt_down_time = 0.6;
  int burn_in_start = 100000;
  int burn_between_start = 1000;
  bool update_burn_time = true;
  bool save_interim_samples = false;
  std::string update_rule = "mh";
  double temperature = 1.0;

  int burn_in;
  int burn_between;

  pcg32 rng;

  arma::Mat<int> samples_2d;
  arma::Cube<int> samples_3d;
  potts_model model;

  Sampler* sampler;
  SampleStats* sample_stats;

  void estimateBurnTime(void);
  bool checkErgodicity(void);
  void loadParameters(std::string);
  void checkParameters(void);
  void setParameter(std::string, std::string);
};
