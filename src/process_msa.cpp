/*
 * SPDX-FileCopyrightText: 2020 - 2021 sudorook <daemon@nullcodon.com>
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
#include <cstdlib>
// #include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.hpp"

/**
 * @brief Print command line usage.
 */
void
print_usage(void)
{
  std::cout << "process_msa usage:" << std::endl;
  std::cout << "  -i: input MSA (FASTA format)" << std::endl;
  std::cout << "  -n: numeric MSA" << std::endl;
  std::cout << "  -w: pre-computed sequence weights (optional)" << std::endl;
  std::cout << "  -g: maximum gap frequency for sequences" << std::endl;
  std::cout << "  -G: maximum gap frequency for positions" << std::endl;
  std::cout << "  -s: maximum sequence similarity" << std::endl;
  std::cout << "  -t: sequence reweighting threshold" << std::endl;
  std::cout << "  -r: flag to reweight sequences before filtering" << std::endl;
  std::cout << "  -p: position (by index number) to protect from filtering"
            << std::endl;
  std::cout << "  -q: sequences (by index number) to protect from filtering"
            << std::endl;
  std::cout << "  -o: output file name (numeric format)" << std::endl;
  std::cout << "  -O: output sequence weights file name" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
}

/**
 * @brief Process a multiple sequence alignment
 *
 * @return return status
 *
 * Compute sequence weights, filter similar sequences, and remove gapped
 * sequences and positions.
 */
int
main(int argc, char* argv[])
{
  std::string infile;
  std::string weight_file;

  std::string outfile;
  std::string outfile_weights;

  bool is_verbose = false;
  bool is_numeric = false;
  bool reweight_first = false;

  double reweighting_threshold = 1.0;
  double similarity_threshold = 1.0;
  double position_gap_threshold = 1.0;
  double sequence_gap_threshold = 1.0;

  std::vector<size_t> keep_seq;
  std::vector<size_t> keep_pos;

  char c;
  while ((c = getopt(argc, argv, "i:n:w:rp:q:g:G:s:t:o:O:hv")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        is_numeric = false;
        break;
      case 'n':
        infile = optarg;
        is_numeric = true;
        break;
      case 'w':
        weight_file = optarg;
        break;
      case 'r':
        reweight_first = true;
        break;
      case 'q':
        keep_seq.push_back(std::stoi(optarg));
        break;
      case 'p':
        keep_pos.push_back(std::stoi(optarg));
        break;
      case 't':
        reweighting_threshold = std::stod(optarg);
        break;
      case 's':
        similarity_threshold = std::stod(optarg);
        break;
      case 'g':
        position_gap_threshold = std::stod(optarg);
        break;
      case 'G':
        sequence_gap_threshold = std::stod(optarg);
        break;
      case 'o':
        outfile = optarg;
        break;
      case 'O':
        outfile_weights = optarg;
        break;
      case 'v':
        is_verbose = true;
        break;
      case 'h':
        print_usage();
        return 0;
        break;
      case '?':
        std::cerr << "ERROR: Invalid input." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
  }

  std::cout << "reading sequences... " << std::flush;
  MSA msa = MSA(infile, weight_file, is_numeric);
  std::cout << "done" << std::endl;

  if (keep_seq.size() > 0) {
    msa.setKeepSequences(keep_seq);
  }

  if (keep_pos.size() > 0) {
    msa.setKeepPositions(keep_pos);
  }

  if (reweight_first) {
    if ((reweighting_threshold < 1) & (similarity_threshold >= 0)) {
      std::cout << "reweighting sequences before filtering... " << std::flush;
      msa.computeSequenceWeights(reweighting_threshold);
      std::cout << "done" << std::endl;
    };
  }

  if ((position_gap_threshold < 1) & (position_gap_threshold >= 0)) {
    std::cout << "filter gapped positions... " << std::flush;
    msa.filterPositionGaps(position_gap_threshold, is_verbose);
    std::cout << "done" << std::endl;
  }

  if ((sequence_gap_threshold < 1) & (sequence_gap_threshold >= 0)) {
    std::cout << "filter gapped sequences... " << std::flush;
    msa.filterSequenceGaps(sequence_gap_threshold, is_verbose);
    std::cout << "done" << std::endl;
  }

  if ((similarity_threshold < 1) & (similarity_threshold >= 0)) {
    std::cout << "filtering similar sequences... " << std::flush;
    msa.filterSimilarSequences(similarity_threshold, is_verbose);
    std::cout << "done" << std::endl;
  }

  if ((reweighting_threshold < 1) & (similarity_threshold >= 0)) {
    std::cout << "computing sequence weights... " << std::flush;
    msa.computeSequenceWeights(reweighting_threshold);
    std::cout << "done" << std::endl;
  }

  if (outfile.size() == 0) {
    int idx = infile.find_last_of(".");
    std::string outfile_prefix = infile.substr(0, idx);
    outfile = outfile_prefix + "_processed.txt";
  }
  if (outfile_weights.size() == 0) {
    if (outfile.size() != 0) {
      int idx = outfile.find_last_of(".");
      std::string outfile_prefix = outfile.substr(0, idx);
      outfile_weights = outfile_prefix + "_weights.txt";
    } else {
      int idx = infile.find_last_of(".");
      std::string outfile_prefix = infile.substr(0, idx);
      outfile_weights = outfile_prefix + "_processed_weights.txt";
    }
  }

  msa.writeMatrix(outfile);
  msa.writeSequenceWeights(outfile_weights);

  return 0;
};
