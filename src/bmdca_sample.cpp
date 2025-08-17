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

#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "generator.hpp"
#include "utils.hpp"

/**
 * @brief Print command line usage.
 */
void
print_usage(void)
{
  std::cout << "bmdca_sample usage:" << std::endl;
  std::cout << "(e.g. bmdca_sample -p <params h> -P <params|params J> -c "
               "<config file> \\"
            << std::endl;
  std::cout << "                   -n <# samples per thread> -r <# threads> \\"
            << std::endl;
  std::cout << "                   -d <directory> -o <output flle>)"
            << std::endl;
  std::cout << "  -p: parameters (txt) _or_ fields h (dat)" << std::endl;
  std::cout << "  -P: couplings J (dat), required if fields h given"
            << std::endl;
  std::cout << "  -d: destination directory" << std::endl;
  std::cout << "  -o: output file name" << std::endl;
  std::cout << "  -r: # independent sampling runs" << std::endl;
  std::cout << "  -n: # samples obtained per sampling run" << std::endl;
  std::cout << "  -c: config file" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
};

/**
 * @brief Sample sequences from a given Potts model.
 *
 * @return Return status.
 */
int
main(int argc, char* argv[])
{
  int num_sequences = 1000; // number of sequence to sample per replicate
  int num_replicates = 10;  // number of independent MC trajectories

  std::string parameters_file, J_file; // file strings for input parameters
  std::string dest_dir = ".";          // directory for storing output
  std::string config_file; // config file for sampler hyperparameters
  std::string output_file = "MC_samples.fasta"; // output file name

  bool dest_dir_given = false;
  bool compat_mode =
    true; // flag for whether the input parameters are text or binary

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "p:P:d:n:c:o:r:h")) != -1) {
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
      case 'h':
        print_usage();
        return 0;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
  }

  // Check inputs
  if (compat_mode == true) {
    if (parameters_file.size() == 0) {
      std::cerr << "ERROR: Parameters file not given." << std::endl;
      print_usage();
      std::exit(EXIT_FAILURE);
    }
  } else {
    if ((parameters_file.size() == 0) || (J_file.size() == 0)) {
      std::cerr << "ERROR: Both parameters files not given." << std::endl;
      print_usage();
      std::exit(EXIT_FAILURE);
    }
  }

  // Load Potts model
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelAscii(parameters_file);
  } else {
    params = loadPottsModel(parameters_file, J_file);
  }

  int N = params.h.n_cols; // number of positions
  int Q = params.h.n_rows; // number of states

  // Initialize the sampler
  Generator generator = Generator(params, N, Q, config_file);

  if (dest_dir_given == true) {
    chdir(dest_dir.c_str());
  }

  // Sample num_replicates x num_sequences sequences
  generator.run(num_replicates, num_sequences, output_file);

  return 0;
};
