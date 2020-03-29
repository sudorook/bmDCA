#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "generator.hpp"
#include "utils.hpp"

int
main(int argc, char* argv[])
{
  int num_sequences = 1000;
  int num_replicates = 10;

  std::string parameters_file, J_file;
  std::string dest_dir = ".";
  std::string config_file;
  std::string output_file = "MC_samples.fasta";

  bool dest_dir_given = false;
  bool compat_mode = true;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "p:P:d:n:c:o:r:")) != -1) {
    switch (c) {
      case 'p':
        parameters_file = optarg;
        break;
      case 'P':
        J_file = optarg;
        compat_mode = false;
        break;
      case 'd':
        dest_dir = optarg;
        {
          struct stat st = { 0 };
          if (stat(dest_dir.c_str(), &st) == -1) {
#if __unix__
            mkdir(dest_dir.c_str(), 0700);
#elif _WIN32
            mkdir(dest_dir.c_str());
#endif
          }
        }
        dest_dir_given = true;
        break;
      case 'o':
        output_file = optarg;
        break;
      case 'n':
        num_sequences = std::stoi(optarg);
        break;
      case 'r':
        num_replicates = std::stoi(optarg);
        break;
      case 'c':
        config_file = optarg;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  // Check inputs
  if (compat_mode == true) {
    if (parameters_file.size() == 0) {
      std::cerr << "ERROR: Parameters file not given." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    if ((parameters_file.size() == 0) || (J_file.size() == 0)) {
      std::cerr << "ERROR: Parameters files not given." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Load Potts model
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelCompat(parameters_file);
  } else {
    params = loadPottsModel(parameters_file, J_file);
  }

  int N = params.h.n_cols;
  int Q = params.h.n_rows;

  Generator generator = Generator(params, N, Q, config_file);

  if (dest_dir_given == true) {
    chdir(dest_dir.c_str());
  }

  generator.run(num_replicates, num_sequences, output_file);

  return 0;
};
