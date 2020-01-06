#ifndef MSA_HPP
#define MSA_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "utils.hpp"

class MSA
{
public:
  arma::Mat<int> alignment;           // numerical multiple sequence alignment
  arma::Col<double> sequence_weights; // weights for each sequence
  int M;                              // number of sequences
  int N;                              // number of positions

  MSA(std::string, bool = true, double = 0.8);
  void printAlignment();
  void writeMatrix(std::string);
  void writeSequenceWeights(std::string);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void makeNumericalMatrix(void);
  void computeSequenceWeights(double);
};

#endif