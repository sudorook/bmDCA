#ifndef MSA_HPP
#define MSA_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "utils.hpp"

class MSA
{
public:
  MSA(std::string, std::string = "", bool = true);
  MSA(arma::Mat<int>, arma::Col<double>, int, int, int);
  MSA(arma::Mat<int>, int, int, int);

  arma::Mat<int> alignment; // numerical multiple sequence alignment
  int M;                    // number of sequences
  int N;                    // number of positions
  int Q;                    // number of amino acids

  arma::Col<double> sequence_weights; // weights for each sequence
  arma::Col<double> hamming_distances;

  void setKeepPositions(std::vector<int>);
  void setKeepSequences(std::vector<int>);

  void filterSequenceGaps(double = 0.2, bool = false);
  void filterPositionGaps(double = 0.2, bool = false);
  void filterSimilarSequences(double = 0.8, bool = false);
  void computeSequenceWeights(double = 0.8);
  void computeHammingDistances(void);

  void printAlignment();
  void writeMatrix(std::string);
  void writeSequenceWeights(std::string);
  void writeHammingDistances(std::string);

  MSA sampleAlignment(int = 0, long int = 0);
  std::vector<MSA*> partitionAlignment(int = 0, long int = 0);

private:
  std::vector<SeqRecord> seq_records;
  int getSequenceLength(std::string);
  void readInputMSA(std::string);
  void readInputNumericMSA(std::string);
  void readSequenceWeights(std::string);
  void makeNumericalMatrix(void);

  std::vector<int> seq_to_keep;
  std::vector<int> pos_to_keep;
};

#endif
