#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>
#include <iostream>

#include "msa.hpp"
// #include "msa_stats.hpp"
#include "run.hpp"

void
print_usage(void)
{
  std::cout << "bmdca usage:" << std::endl;
  std::cout << "(e.g. bmdca -i <training MSA> -w <training sequence weights> -d <output directory> -c <config file>)"
            << std::endl;
  std::cout << "  -i: training MSA (FASTA format)" << std::endl;
  std::cout << "  -I: validation MSA (FASTA format)" << std::endl;
  std::cout << "  -d: destination directory" << std::endl;
  std::cout << "  -n: training numerical MSA" << std::endl;
  std::cout << "  -N: validation numerical MSA" << std::endl;
  std::cout << "  -w: training sequence weights" << std::endl;
  std::cout << "  -W: validation sequence weights" << std::endl;
  std::cout << "  -c: config file" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
  std::cout << "  -f: force a restart of the inference loop" << std::endl;
};

int
main(int argc, char* argv[])
{
  std::string train_file;
  std::string train_weight_file;
  
  std::string validate_file;
  std::string validate_weight_file;

  std::string config_file;
  std::string dest_dir = ".";

  bool is_numeric = false;
  bool force_restart = false;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "c:d:fhi:I:n:N:w:W:")) != -1) {
    switch (c) {
      case 'c':
        config_file = optarg;
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
        break;
      case 'f':
        force_restart = true;
        break;
      case 'h':
        print_usage();
        return 0;
        break;
      case 'i':
        train_file = optarg;
        break;
      case 'I':
        validate_file = optarg;
        break;
      case 'n':
        train_file = optarg;
        is_numeric = true;
        break;
      case 'N':
        validate_file = optarg;
        is_numeric = true;
        break;
      case 'w':
        train_weight_file = optarg;
        break;
      case 'W':
        validate_weight_file = optarg;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
  }

  // Exit if no input file was provided.
  if (train_file.empty()) {
    std::cerr << "ERROR: Missing MSA input file." << std::endl;
    print_usage();
    std::exit(EXIT_FAILURE);
  }

  MSA *msa_train = nullptr;
  MSA *msa_validate = nullptr;

  // Parse the multiple sequence alignment.
  msa_train = new MSA(train_file, train_weight_file, is_numeric);
  msa_train->writeSequenceWeights(dest_dir + "/msa_weights.txt");
  msa_train->writeMatrix(dest_dir + "/msa_numerical.txt");
  
  if (!validate_file.empty()) {
    msa_validate = new MSA(validate_file, validate_weight_file, is_numeric);
    msa_validate->writeSequenceWeights(dest_dir + "/msa_validate_weights.txt");
    msa_validate->writeMatrix(dest_dir + "/msa_validate_numerical.txt");
  }

  Sim sim = Sim(msa_train, msa_validate, config_file, dest_dir, force_restart);
  sim.run();

  delete msa_train;
  delete msa_validate;
  
  return 0;
};
