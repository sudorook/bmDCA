#include "generator.hpp"

#include "sampler.hpp"
#include "sample_stats.hpp"
#include "utils.hpp"

Generator::Generator(potts_model params, int n, int q, std::string config_file)
  : N(n)
  , Q(q)
  , model(params)
{
  if (config_file.length() != 0) {
    loadParameters(config_file);
  }
};

Generator::~Generator(void)
{
  delete sample_stats;
  delete sampler;
};

void
Generator::checkParameters(void){};

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

void
Generator::writeAASequences(std::string output_file)
{
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

bool
Generator::checkErgodicity(void)
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
    std::cout << "increasing burn-between time to " << burn_between << std::endl;
  } else if (flag_deltat_down) {
    burn_between = Max((int)(round((double)burn_between * adapt_down_time)), 1);
    std::cout << "decreasing burn-between time to " << burn_between << std::endl;
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
Generator::estimateBurnTime(void)
{ 
  std::uniform_int_distribution<long int> dist(0, RAND_MAX - walkers);
  
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
    // sampler->setBurnTime(burn_in, burn_in);
  }
};


void
Generator::writeNumericalSequences(std::string output_file)
{
  int idx = output_file.find_last_of(".");
  std::string raw_file = output_file.substr(0, idx);

  sample_stats->writeSamples(raw_file + "_numerical.txt");
  sample_stats->writeSampleEnergies(raw_file + "_energies.txt");
}

void
Generator::run(int n_indep_runs, int n_per_run, std::string output_file)
{
  std::cout << "initializing sampler... " << std::flush;

  arma::wall_clock timer;
  timer.tic();

  walkers = n_indep_runs;
  samples_per_walk = n_per_run;

  checkParameters();

  if (samples_per_walk == 1) {
    samples_2d = arma::Mat<int>(walkers, N, arma::fill::zeros);
    sample_stats = new SampleStats2D(&samples_2d, &(model));
  } else {
    samples_3d = arma::Cube<int>(samples_per_walk, N, walkers, arma::fill::zeros);
    sample_stats = new SampleStats3D(&samples_3d, &(model));
  }

  // Instantiate the PCG random number generator and unifrom random
  // distribution.
  rng.seed(random_seed);
  std::uniform_int_distribution<long int> dist(0, RAND_MAX - walkers);

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
    if(samples_per_walk > 1) {
      if (update_rule == "mh") {
        sampler->sampleSequences(&samples_3d,
                                 walkers,
                                 samples_per_walk,
                                 burn_in,
                                 burn_between,
                                 dist(rng),
                                 temperature);
      } else if (update_rule == "z-sqrt") {
        sampler->sampleSequencesZanella(
          &samples_3d, walkers, samples_per_walk, burn_in, burn_between, dist(rng), "sqrt", temperature);
      } else if (update_rule == "z-barker") {
        sampler->sampleSequencesZanella(
          &samples_3d, walkers, samples_per_walk, burn_in, burn_between, dist(rng), "barker", temperature);
      } else {
        std::cerr << "ERROR: sampler '" << sampler << "' not recognized."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      if (update_rule == "mh") {
        sampler->sampleSequences(
          &samples_2d, walkers, burn_in, dist(rng), temperature);
      } else if (update_rule == "z-sqrt") {
        sampler->sampleSequencesZanella(
          &samples_2d, walkers, burn_in, dist(rng), "sqrt", temperature);
      } else if (update_rule == "z-barker") {
        sampler->sampleSequencesZanella(
          &samples_2d, walkers, burn_in, dist(rng), "barker", temperature);
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
