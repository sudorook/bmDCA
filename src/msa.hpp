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

#ifndef SRC_MSA_HPP_
#define SRC_MSA_HPP_

#include <armadillo>
#include <memory>
#include <string>
#include <vector>

#include "utils.hpp"

/**
 * @brief Class for loading, storing, and processing a MSA.
 *
 * The MSA class can read in sequences in FASTA or numerical format, process
 * alignment based on similarity and gap frequencies, compute sequence weights,
 * and sub-sample alignments.
 */
class MSA
{
public:
  explicit MSA(std::string = " ", std::string = " ", bool = true);
  MSA(arma::Mat<int>, arma::Col<double>, int, int, int);
  MSA(arma::Mat<int>, int, int, int);

  arma::Mat<int> alignment; ///< numerical multiple sequence alignment
  int M;                    ///< number of sequences
  int N;                    ///< number of positions
  int Q;                    ///< number of amino acids

  arma::Col<double> sequence_weights;  ///< weights for each sequence
  arma::Col<double> hamming_distances; ///< hamming distances among sequences

  void setKeepPositions(std::vector<size_t>);
  void setKeepSequences(std::vector<size_t>);

  void filterSequenceGaps(double = 0.2, bool = false);
  void filterPositionGaps(double = 0.2, bool = false);
  void filterSimilarSequences(double = 0.8, bool = false);
  void computeSequenceWeights(double = 0.8);
  void computeHammingDistances(void);

  void printAlignment();
  void writeMatrix(std::string);
  void writeSequenceWeights(std::string);
  void writeHammingDistances(std::string);

  MSA sampleAlignment(int = 0, unsigned = 0);
  std::pair<std::shared_ptr<MSA>, std::shared_ptr<MSA>> partitionAlignment(
    int = 0,
    unsigned = 0);

private:
  std::vector<SeqRecord>
    seq_records; ///< vector of sequences loaded from a FASTA file
  int getAASequenceLength(const std::string&);
  int getNTSequenceLength(const std::string&);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void readSequenceWeights(std::string);
  void makeAANumericalMatrix(void);
  void makeNTNumericalMatrix(void);

  std::vector<size_t>
    seq_to_keep; ///< indices for sequences to protect from pruning
  std::vector<size_t>
    pos_to_keep; ///< indices for positions to protect from pruning
};

#endif // SRC_MSA_HPP_
