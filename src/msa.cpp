#include "msa.hpp"

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>

#include "pcg_random.hpp"

#ifndef AA_ALPHABET_SIZE
#define AA_ALPHABET_SIZE 21
#endif

MSA::MSA(std::string msa_file,
         std::string weight_file,
         bool reweight,
         bool is_numeric_msa,
         double threshold)
  : reweight(reweight)
  , threshold(threshold)
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
  if (reweight) {
    computeSequenceWeights(threshold);
  } else if (!weight_file.empty()) {
    readSequenceWeights(weight_file);
  } else {
    sequence_weights = arma::Col<double>(M, arma::fill::ones);
  }
};

MSA::MSA(arma::Mat<int> alignment,
         int M,
         int N,
         int Q,
         bool reweight,
         double threshold)
  : alignment(alignment)
  , M(M)
  , N(N)
  , Q(Q)
  , reweight(reweight)
  , threshold(threshold)
{
  if (reweight) {
    computeSequenceWeights(threshold);
  } else {
    sequence_weights = arma::Col<double>(M, arma::fill::ones);
  }
};

void
MSA::readInputNumericMSA(std::string numeric_msa_file)
{
  std::ifstream input_stream(numeric_msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << numeric_msa_file
              << "' for reading." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  input_stream >> M >> N >> Q;
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
      }
    }
    row_idx++;
  }
};

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

void
MSA::printAlignment(void)
{
  for (std::vector<SeqRecord>::iterator it = seq_records.begin();
       it != seq_records.end();
       ++it) {
    it->print();
  }
};

int
MSA::getSequenceLength(std::string sequence)
{
  int valid_aa_count = 0;
  for (std::string::iterator it = sequence.begin(); it != sequence.end();
       ++it) {
    switch (*it) {
      case '-':
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
    }
  }
  return valid_aa_count;
};

void
MSA::computeSequenceWeights(double threshold)
{
  sequence_weights = arma::Col<double>(M, arma::fill::ones);
  arma::Mat<int> alignment_T = alignment.t();

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
          if (id > threshold * N) {
            sequence_weights(m1) += 1;
          }
        }
      }
    }
  }

  sequence_weights = 1. / sequence_weights;
};

void
MSA::writeSequenceWeights(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << sequence_weights(i) << std::endl;
  }
};

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

void
MSA::writeHammingDistances(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << hamming_distances(i) << std::endl;
  }
};

MSA
MSA::subsampleAlignment(int size, long int seed)
{
  pcg32 rng;
  rng.seed(seed);

  arma::Mat<int> alignment_T = alignment.t();
  arma::Mat<int> sub_alignment =
    arma::Mat<int>(N, size, arma::fill::zeros);

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

  return MSA(sub_alignment.t(), size, N, Q, reweight, threshold);
};
