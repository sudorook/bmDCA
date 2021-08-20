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

#include "msa.hpp"

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "pcg_random.hpp"

#ifndef AA_ALPHABET_SIZE
#define AA_ALPHABET_SIZE 21
#endif

/**
 * @brief Constructor for reading MSA (and weights) from file.
 *
 * @param msa_file file string for input MSA
 * @param weight_file file string for input MSA weights
 * @param is_numeric_msa (bool) flag for input MSA format
 */
MSA::MSA(std::string msa_file, std::string weight_file, bool is_numeric_msa)
{
  if (is_numeric_msa) {
    readInputNumericMSA(msa_file);
  } else {
    readInputMSA(msa_file);
    M = seq_records.size();
    N = getSequenceLength(seq_records.begin()->getSequence());
    Q = AA_ALPHABET_SIZE;
    makeNumericalMatrix();
  }
  if (!weight_file.empty()) {
    readSequenceWeights(weight_file);
  } else {
    sequence_weights = arma::Col<double>(M, arma::fill::ones);
  }
};

/**
 * @brief Constructor for existing MSA matrix (arma::Mat).
 *
 * @param alignment arma::Mat containing numerical alignment
 * @param M number of sequences
 * @param N number of positions
 * @param Q number of states
 */
MSA::MSA(arma::Mat<int> alignment, int M, int N, int Q)
  : alignment(alignment)
  , M(M)
  , N(N)
  , Q(Q)
{
  sequence_weights = arma::Col<double>(M, arma::fill::ones);
};

/**
 * @brief Constructor for existing MSA matrix and weights vector.
 *
 * @param alignment arma::Mat containing numerical alignment
 * @param weights arma::Col containing sequence weights
 * @param M number of sequences
 * @param N number of positions
 * @param Q number of states
 */
MSA::MSA(arma::Mat<int> alignment,
         arma::Col<double> weights,
         int M,
         int N,
         int Q)
  : alignment(alignment)
  , M(M)
  , N(N)
  , Q(Q)
  , sequence_weights(weights){};

/**
 * @brief Read numerical MSA from file.
 *
 * @param numeric_msa_file input numeric msa file string
 */
void
MSA::readInputNumericMSA(std::string numeric_msa_file)
{
  std::ifstream input_stream(numeric_msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << numeric_msa_file
              << "' for reading." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  input_stream >> M >> N >> Q; // read header
  alignment = arma::Mat<int>(M, N);

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
}

/**
 * @brief Read sequence weights from file.
 *
 * @param weights_file file string for input sequence weights.
 */
void
MSA::readSequenceWeights(std::string weights_file)
{
  std::ifstream input_stream(weights_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << weights_file << "' for reading."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  sequence_weights = arma::Col<double>(M, arma::fill::zeros);

  std::string line;
  int counter = 0;
  while (std::getline(input_stream, line)) {
    std::istringstream iss(line);
    double n;

    while (iss >> n) {
      sequence_weights(counter) = n;
    }
    counter++;
  }
}

/**
 * @brief Read FASTA-format input MSA.
 *
 * @param msa_file file string for input MSA.
 */
void
MSA::readInputMSA(std::string msa_file)
{
  std::ifstream input_stream(msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: cannot read from '" << msa_file << "'." << std::endl;
    exit(2);
  }

  /*
   * Read a FASTA-formatted multiple sequence alignment. Each record from the
   * file is stored as a SeqRecord object and appended to the seq_records
   * vector.
   */
  std::string header, sequence, line;
  while (input_stream) {
    std::getline(input_stream, line);
    if (line[0] == '>') {
      if (sequence.length() > 0) {
        seq_records.push_back(SeqRecord(header, sequence));
        sequence.clear();
        header.clear();
      }
      header = line;
      header.erase(0, 1);
    } else {
      sequence += line;
      line.clear();
    }
  };
  seq_records.push_back(SeqRecord(header, sequence));
  input_stream.close();
};

/**
 * @brief Convert a FASTA-formatted alignment into a numerical matrix.
 */
void
MSA::makeNumericalMatrix(void)
{
  alignment = arma::Mat<int>(M, N);

  int row_idx = 0;
  for (auto seq = seq_records.begin(); seq != seq_records.end(); seq++) {
    std::string sequence = seq->getSequence();
    int col_idx = 0;
    for (auto aa = sequence.begin(); aa != sequence.end(); aa++) {
      switch (*aa) {
        case '-':
        case '.':
        case 'B':
        case 'J':
        case 'O':
        case 'U':
        case 'X':
        case 'Z':
          alignment(row_idx, col_idx) = 0;
          col_idx++;
          break;
        case 'A':
          alignment(row_idx, col_idx) = 1;
          col_idx++;
          break;
        case 'C':
          alignment(row_idx, col_idx) = 2;
          col_idx++;
          break;
        case 'D':
          alignment(row_idx, col_idx) = 3;
          col_idx++;
          break;
        case 'E':
          alignment(row_idx, col_idx) = 4;
          col_idx++;
          break;
        case 'F':
          alignment(row_idx, col_idx) = 5;
          col_idx++;
          break;
        case 'G':
          alignment(row_idx, col_idx) = 6;
          col_idx++;
          break;
        case 'H':
          alignment(row_idx, col_idx) = 7;
          col_idx++;
          break;
        case 'I':
          alignment(row_idx, col_idx) = 8;
          col_idx++;
          break;
        case 'K':
          alignment(row_idx, col_idx) = 9;
          col_idx++;
          break;
        case 'L':
          alignment(row_idx, col_idx) = 10;
          col_idx++;
          break;
        case 'M':
          alignment(row_idx, col_idx) = 11;
          col_idx++;
          break;
        case 'N':
          alignment(row_idx, col_idx) = 12;
          col_idx++;
          break;
        case 'P':
          alignment(row_idx, col_idx) = 13;
          col_idx++;
          break;
        case 'Q':
          alignment(row_idx, col_idx) = 14;
          col_idx++;
          break;
        case 'R':
          alignment(row_idx, col_idx) = 15;
          col_idx++;
          break;
        case 'S':
          alignment(row_idx, col_idx) = 16;
          col_idx++;
          break;
        case 'T':
          alignment(row_idx, col_idx) = 17;
          col_idx++;
          break;
        case 'V':
          alignment(row_idx, col_idx) = 18;
          col_idx++;
          break;
        case 'W':
          alignment(row_idx, col_idx) = 19;
          col_idx++;
          break;
        case 'Y':
          alignment(row_idx, col_idx) = 20;
          col_idx++;
          break;
        default:
          std::cerr << *aa << std::endl;
          std::exit(EXIT_FAILURE);
          break;
      }
    }
    row_idx++;
  }
};

/**
 * @brief Write numerical alignment to disk.
 *
 * @param output_file file string for output file
 */
void
MSA::writeMatrix(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << M << " " << N << " " << Q << std::endl;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      if (j + 1 == N) {
        output_stream << alignment(i, j) << std::endl;
      } else {
        output_stream << alignment(i, j) << " ";
      }
    }
  }
};

/**
 * @brief Print the sequences.
 */
void
MSA::printAlignment(void)
{
  for (std::vector<SeqRecord>::iterator it = seq_records.begin();
       it != seq_records.end();
       ++it) {
    it->print();
  }
};

/**
 * @brief Count the number of valid positions in a sequence.
 *
 * @param sequence sequence string
 *
 * @return number of positions
 */
int
MSA::getSequenceLength(std::string sequence)
{
  int valid_aa_count = 0;
  for (std::string::iterator it = sequence.begin(); it != sequence.end();
       ++it) {
    switch (*it) {
      case '-':
      case '.':
      case 'B':
      case 'J':
      case 'O':
      case 'U':
      case 'X':
      case 'Z':
      case 'A':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
      case 'K':
      case 'L':
      case 'M':
      case 'N':
      case 'P':
      case 'Q':
      case 'R':
      case 'S':
      case 'T':
      case 'V':
      case 'W':
      case 'Y':
        valid_aa_count += 1;
        break;
      default:
        std::cerr << *it << std::endl;
        std::exit(EXIT_FAILURE);
        break;
    }
  }
  return valid_aa_count;
};

/**
 * @brief Compute sequence weights based on a similarity threshold.
 *
 * @param threshold for the fraction of allowable sequence similarity
 */
void
MSA::computeSequenceWeights(double threshold)
{
  sequence_weights = arma::Col<double>(M, arma::fill::ones);
  double cutoff = N * threshold;
  arma::Mat<int> alignment_T = alignment.t();

  // arma::Mat is stored in column-major format, but sequences are stored in
  // rows. Run operations on transposed alignment matrix to keep things in the
  // same cache lines.
#pragma omp parallel
  {
#pragma omp for
    for (int m1 = 0; m1 < M; ++m1) {
      int* m1_ptr = alignment_T.colptr(m1);
      for (int m2 = 0; m2 < M; ++m2) {
        if (m1 != m2) {
          int* m2_ptr = alignment_T.colptr(m2);
          double id = 0;
          for (int i = 0; i < N; ++i) {
            if (*(m1_ptr + i) == *(m2_ptr + i)) {
              id += 1;
            }
          }
          if (id > cutoff) {
            sequence_weights(m1) += 1;
          }
        }
      }
    }
  }

  sequence_weights = 1. / sequence_weights;
};

/**
 * @brief Remove sequences above a threshold sequence similarity.
 *
 * @param threshold similarity threshold
 * @param verbose flag to print which sequences are pruned
 *
 * If two (or more) sequences are too similary, _both_ (or _all_) are removed
 * from the alignment.
 */
void
MSA::filterSimilarSequences(double threshold, bool verbose)
{
  arma::Col<int> sequence_status = arma::Col<int>(M, arma::fill::zeros);
  // int sim_cutoff = (int)(N * threshold);
  double sim_cutoff = N * threshold;
  arma::Mat<int> alignment_T = alignment.t();
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 1)
    for (int m1 = 0; m1 < M; ++m1) {
      int* m1_ptr = alignment_T.colptr(m1);
      for (int m2 = 0; m2 < M; ++m2) {
        if (m1 != m2) {
          int* m2_ptr = alignment_T.colptr(m2);
          double id = 0;
          for (int i = 0; i < N; ++i) {
            if (*(m1_ptr + i) == *(m2_ptr + i)) {
              id += 1;
            }
          }
          if (id > sim_cutoff) {
            sequence_status(m1) = 1;
            break;
          }
        }
      }
    }
  }

  // Collect row indices for sequences to remove, but first make sure to keep
  // sequences that are explicitly protected.
  arma::uvec bad_sequences = arma::find(sequence_status == 1);
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == seq_to_keep[i]));
  }
  if (verbose) {
    for (size_t i = 0; i < bad_sequences.size(); i++) {
      if (i == bad_sequences.size() - 1) {
        std::cout << bad_sequences[i] << "... ";
      } else {
        std::cout << bad_sequences[i] << " ";
      }
    }
    std::cout << std::flush;
  }
  alignment.shed_rows(bad_sequences); // drop too-similar sequences
  M = alignment.n_rows; // re-set the number of sequences

  // Once sequences are removed, the indices for sequences to keep need to be
  // updated.
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    int offset = 0;
    for (size_t j = 0; j < bad_sequences.size(); j++) {
      if (bad_sequences.at(j) < seq_to_keep.at(i)) {
        offset++;
      }
    }
    seq_to_keep.at(i) -= offset;
  }
};

/**
 * @brief Set position indicies to protect from processing.
 *
 * @param keep vector of indices
 */
void
MSA::setKeepPositions(std::vector<size_t> keep) {
  pos_to_keep = keep;
};

/**
 * @brief Set sequence indicies to protect from processing.
 *
 * @param keep vector of indices
 */
void
MSA::setKeepSequences(std::vector<size_t> keep) {
  seq_to_keep = keep;
};

/**
 * @brief Remove sequences that have too many gaps.
 *
 * @param seq_threshold maximum tolerable fraction of gaps
 * @param verbose flag to print which sequences are removed
 */
void
MSA::filterSequenceGaps(double seq_threshold, bool verbose)
{
  arma::Col<int> seq_gap_counts = arma::Col<int>(M, arma::fill::zeros);

  int seq_gap_cutoff = (int)(seq_threshold * N);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        if (alignment(i, j) == 0) {
          seq_gap_counts(i)++;
        }
      }
    }
  }

  // Collect row indices for sequences to remove, but first make sure to keep
  // sequences that are explicitly protected.
  arma::uvec bad_sequences = arma::find(seq_gap_counts > seq_gap_cutoff);
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == seq_to_keep[i]));
  }
  if (verbose) {
    for (size_t i = 0; i < bad_sequences.size(); i++) {
      if (i == bad_sequences.size() - 1) {
        std::cout << bad_sequences[i] << "... ";
      } else {
        std::cout << bad_sequences[i] << " ";
      }
    }
    std::cout << std::flush;
  }
  alignment.shed_rows(bad_sequences); // drop overly-gapped sequences
  M = alignment.n_rows; // re-set the number of sequences

  // Once sequences are removed, the indices for sequences to keep need to be
  // updated.
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    int offset = 0;
    for (size_t j = 0; j < bad_sequences.size(); j++) {
      if (bad_sequences.at(j) < seq_to_keep.at(i)) {
        offset++;
      }
    }
    seq_to_keep.at(i) -= offset;
  }
};

/**
 * @brief Remove positions with too many gaps.
 *
 * @param pos_threshold maximum tolerable fraction of gaps
 * @param verbose flag to print the positions removed
 */
void
MSA::filterPositionGaps(double pos_threshold, bool verbose)
{
  arma::Col<int> pos_gap_counts = arma::Col<int>(N, arma::fill::zeros);

  int pos_gap_cutoff = (int)(pos_threshold * arma::sum(sequence_weights));

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        if (alignment(j, i) == 0) {
          pos_gap_counts(i) += sequence_weights(j);
        }
      }
    }
  }

  // Collect column indices for positions to remove, but first make sure to
  // keep positions that are explicitly protected.
  arma::uvec bad_positions = arma::find(pos_gap_counts > pos_gap_cutoff);
  for (size_t i = 0; i < pos_to_keep.size(); i++) {
    bad_positions.shed_rows(arma::find(bad_positions == pos_to_keep[i]));
  }
  if (verbose) {
    for (size_t i = 0; i < bad_positions.size(); i++) {
      if (i == bad_positions.size() - 1) {
        std::cout << bad_positions[i] << "... ";
      } else {
        std::cout << bad_positions[i] << " ";
      }
    }
    std::cout << std::flush;
  }
  alignment.shed_cols(bad_positions); // drop overly-gapped positions
  N = (int)alignment.n_cols; // re-set the number of positions

  // Once positions are removed, the indices for positions to keep need to be
  // updated.
  for (size_t i = 0; i < pos_to_keep.size(); i++) {
    int offset = 0;
    for (size_t j = 0; j < bad_positions.size(); j++) {
      if (bad_positions.at(j) < pos_to_keep.at(i)) {
        offset++;
      }
    }
    pos_to_keep.at(i) -= offset;
  }
};

/**
 * @brief Compute maximum sequence similarity for each sequence.
 */
void
MSA::computeHammingDistances(void)
{
  hamming_distances = arma::Col<double>(M, arma::fill::zeros);
  arma::Mat<int> alignment_T = alignment.t();

  int* i_ptr = nullptr;
  int* j_ptr = nullptr;
  int count = 0;
  double id = 0;
  for (int i = 0; i < M; i++) {
    i_ptr = alignment_T.colptr(i);
    for (int j = i + 1; j < M; j++) {
      count = 0;
      j_ptr = alignment_T.colptr(j);
      for (int n = 0; n < N; n++) {
        if (*(i_ptr + n) == *(j_ptr + n)) {
          count++;
        }
      }
      id = (double)count / N;
      if (id > hamming_distances(i)) {
        hamming_distances(i) = id;
      }
    }
  }
};

/**
 * @brief Write sequence weights to disk.
 *
 * @param output_file file string for output
 */
void
MSA::writeSequenceWeights(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << sequence_weights(i) << std::endl;
  }
};

/**
 * @brief Write max sequence similarities to file.
 *
 * @param output_file file string for output
 */
void
MSA::writeHammingDistances(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << hamming_distances(i) << std::endl;
  }
};

/**
 * @brief Sample a random set of sequences from the alignment.
 *
 * @param size number of sequences to sample
 * @param seed random seed for the RNG
 *
 * @return MSA instance with subset of sequences
 */
MSA
MSA::sampleAlignment(int size, long int seed)
{
  pcg32 rng;
  rng.seed(seed);

  arma::Mat<int> alignment_T = alignment.t();
  arma::Mat<int> sub_alignment = arma::Mat<int>(N, size, arma::fill::zeros);

  std::unordered_set<int> elems;
  for (int r = M - size; r < M; ++r) {
    std::uniform_int_distribution<int> dist(0, r);
    int v = dist(rng);
    if (!elems.insert(v).second) {
      elems.insert(r);
    }
  }
  std::vector<int> idx(elems.begin(), elems.end());

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < size; i++) {
      sub_alignment.col(i) = alignment_T.col(idx.at(i));
    }
  }

  return MSA(sub_alignment.t(), size, N, Q);
};

/**
 * @brief Split an alignment into two random sets
 *
 * @param validation_size number of sequences in second set
 * @param seed random seed for RNG
 *
 * @return vector of two MSA objects for each subset
 *
 * This function is used for splitting an MSA into a training set and a
 * validation set. The validation_size corresponds to the effective number of
 * sequences in the validation set. Weights in the two subsets correspond to
 * the weights in the original full MSA.
 */
std::vector<MSA*>
MSA::partitionAlignment(int validation_size, long int seed)
{
  arma::arma_rng::set_seed(seed);

  if (validation_size >= arma::sum(sequence_weights)) {
    std::cerr << "ERROR: cannot sample " << validation_size
              << " from alignment size " << M << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Mat<int> alignment_T = alignment.t();

  // Select sequences for the validation alignment such that the cumulative
  // effective sequence size is greater than the given threshold.
  arma::uvec idx = arma::randperm(M);
  double accu = 0;
  int size = 0;
  for (int i = 0; i < M; i++) {
    size++;
    accu += sequence_weights(idx(i));
    if (accu > (double)validation_size)
      break;
  }

  // don't use more than a quarter of the MSA for cross-validation
  size = (int)Min((double)size, (double)M / 4);

  arma::Mat<int> sub_alignment_1 =
    arma::Mat<int>(N, M - size, arma::fill::zeros);
  arma::Col<double> weights_1 = arma::Col<double>(M - size, arma::fill::ones);

  arma::Mat<int> sub_alignment_2 = arma::Mat<int>(N, size, arma::fill::zeros);
  arma::Col<double> weights_2 = arma::Col<double>(size, arma::fill::ones);

#pragma omp parallel
  {
#pragma omp for
    for (int i = size; i < M; i++) {
      sub_alignment_1.col(i - size) = alignment_T.col(idx(i));
      weights_1(i - size) = sequence_weights(idx(i));
    }
  }

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < size; i++) {
      sub_alignment_2.col(i) = alignment_T.col(idx(i));
      weights_2(i) = sequence_weights(idx(i));
    }
  }

  std::vector<MSA*> return_vector;
  MSA* tmp;

  tmp = new MSA(sub_alignment_1.t(), weights_1, M - size, N, Q);
  return_vector.push_back(tmp);
  tmp = new MSA(sub_alignment_2.t(), weights_2, size, N, Q);
  return_vector.push_back(tmp);

  return return_vector;
};
