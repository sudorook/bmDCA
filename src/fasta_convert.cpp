/* Boltzmann-machine Direct Coupling Analysis (bmDCA)
 * Copyright (C) 2020
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "utils.hpp"

/**
 * @brief Print command line usage.
 */
void
print_usage(void)
{
  std::cout << "mumeric2fasta usage:" << std::endl;
  std::cout << "(e.g. numeric2fasta -n <numeric MSA> -o <output FASTA file>"
            << std::endl;
  std::cout << "  -n: numeric MSA (only 21 states supported)" << std::endl;
  std::cout << "  -o: output FASTA file" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
};

/**
 * @brief Convert numerical MSA to FASTA format.
 *
 * @return return status
 */
int
main(int argc, char* argv[])
{
  std::string numeric_file;
  std::string fasta_file = "output.fasta";

  bool valid_input = false;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "n:o:h")) != -1) {
    switch (c) {
      case 'n':
        numeric_file = optarg;
        valid_input = true;
        break;
      case 'o':
        fasta_file = optarg;
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

  if (!valid_input) {
    std::cerr << "ERROR: Input file not given." << std::endl;
    print_usage();
    std::exit(EXIT_FAILURE);
  }

  std::ifstream input_stream(numeric_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << numeric_file << "' for reading."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int N;
  int M;
  int Q;

  input_stream >> M >> N >> Q;

  if (Q != 21) {
    std::cerr << "ERROR: only 21 amino acids (states) supported." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Mat<int> alignment = arma::Mat<int>(M, N);

  int counter = 0;
  int i = 0;
  std::string line;
  std::getline(input_stream, line);
  while (std::getline(input_stream, line)) {
    std::istringstream iss(line);
    int n;
    i = 0;

    while (iss >> n) {
      alignment(counter, i) = n;
      i++;
    }
    counter++;
  }

  std::ofstream output_stream(fasta_file);

  char aa;
  for (int m = 0; m < M; m++) {
    output_stream << ">sample" << m << std::endl;
    for (int n = 0; n < N; n++) {
      aa = convertAA(alignment(m, n));
      if (aa != '\0') {
        output_stream << aa;
      }
    }
    output_stream << std::endl << std::endl;
  }
  output_stream.close();

  return 0;
}
