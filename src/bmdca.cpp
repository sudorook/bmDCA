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

#include <iostream>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "msa.hpp"
#include "run.hpp"

/**
 * @brief Print command line usage.
 */
void
print_usage(void)
{
  std::cout << "bmdca usage:" << std::endl;
  std::cout << "(e.g. bmdca -i <training MSA> -w <training sequence weights> "
               "-d <output directory> -c <config file>)"
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

/**
 * @brief Run inference loop with input MSA(s).
 *
 * @return return status
 */
int
main(int argc, char* argv[])
{
  std::string train_file;        // file string for training MSA
  std::string train_weight_file; // file string for training sequence weights

  std::string validate_file;        // file string for validation MSA
  std::string validate_weight_file; // file string for validation weights

  std::string config_file;    // file string for bmDCA config file
  std::string dest_dir = "."; // output destination directory

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

  std::shared_ptr<MSA> msa_train =
    std::make_shared<MSA>(train_file, train_weight_file, is_numeric);
  std::shared_ptr<MSA> msa_validate(nullptr);

  // Only load the validation MSA if a validation MSA file is given.
  if (!validate_file.empty()) {
    msa_validate =
      std::make_shared<MSA>(validate_file, validate_weight_file, is_numeric);
  }

  Sim sim = Sim(msa_train, msa_validate, config_file, dest_dir, force_restart);
  sim.run(); // run inference loop

  return 0;
};
