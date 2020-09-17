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

MSA::MSA(arma::Mat<int> alignment,
         int M,
         int N,
         int Q)
  : alignment(alignment)
  , M(M)
  , N(N)
  , Q(Q)
{
  sequence_weights = arma::Col<double>(M, arma::fill::ones);
};

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
MSA::filterSimilarSequences(double threshold)
{
  arma::Col<int> sequence_status = arma::Col<int>(M, arma::fill::zeros);
  int sim_cutoff = (int)(N * threshold);
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

  arma::uvec bad_sequences = arma::find(sequence_status == 1);
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == seq_to_keep[i]));
  }
  alignment.shed_rows(bad_sequences);
  M = alignment.n_rows;
};

void
MSA::filterSequenceGaps(double seq_threshold)
{
  arma::Col<int> seq_gap_counts = arma::Col<int>(M, arma::fill::zeros);

  int seq_gap_cutoff = (int)(seq_threshold * N);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++ ) {
        // seq_gap_counts(i) += alignment(i, Q * j);
        if (alignment(i, j) == 0) {
          seq_gap_counts(i)++;
        }
      }
    }
  }

  arma::uvec bad_sequences = arma::find(seq_gap_counts > seq_gap_cutoff);
  for (size_t i = 0; i < seq_to_keep.size(); i++) {
    bad_sequences.shed_rows(arma::find(bad_sequences == seq_to_keep[i]));
  }
  alignment.shed_rows(bad_sequences);
  M = alignment.n_rows;
};

void
MSA::filterPositionGaps(double pos_threshold)
{
  arma::Col<int> pos_gap_counts =
    arma::Col<int>(N, arma::fill::zeros);

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

  arma::uvec bad_positions = arma::find(pos_gap_counts > pos_gap_cutoff);
  for (size_t i = 0; i < pos_to_keep.size(); i++) {
    bad_positions.shed_rows(arma::find(bad_positions == pos_to_keep[i]));
  }
  alignment.shed_cols(bad_positions);
  N = (int)alignment.n_cols;
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
MSA::writeSequenceWeights(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << sequence_weights(i) << std::endl;
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
MSA::sampleAlignment(int size, long int seed)
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

  return MSA(sub_alignment.t(), size, N, Q);
};

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
  size = Min(size, (int) M / 4);

  arma::Mat<int> sub_alignment_1 =
    arma::Mat<int>(N, M - size, arma::fill::zeros);
  arma::Col<double> weights_1 =
    arma::Col<double>(M - size, arma::fill::ones);

  arma::Mat<int> sub_alignment_2 =
    arma::Mat<int>(N, size, arma::fill::zeros);
  arma::Col<double> weights_2 =
    arma::Col<double>(size, arma::fill::ones);

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
