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
  int Q;                              // number of amino acids
  const bool reweight;                // whether reweighting performed
  const double threshold;             // reweighting threshold

  arma::Col<double> hamming_distances;

  MSA(std::string, std::string = "", bool = true, bool = false, double = 0.8);
  MSA(arma::Mat<int>, int, int, int, bool = true, double = 0.8);
  void printAlignment();
  void writeMatrix(std::string);
  void writeSequenceWeights(std::string);
  void writeHammingDistances(std::string);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void readSequenceWeights(std::string);
  void makeNumericalMatrix(void);
  void computeSequenceWeights(double);
  void computeHammingDistances(void);
};

#endif
