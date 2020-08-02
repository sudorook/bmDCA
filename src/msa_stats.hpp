#ifndef MSA_STATS_HPP
#define MSA_STATS_HPP

#include "msa.hpp"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA*);

  double getEffectiveM();
  double getN();
  double getM();
  double getQ();
  void writeRelEntropy(std::string);
  void writeRelEntropyAscii(std::string);
  // void writeRelEntropyPos(std::string);
  void writeRelEntropyPosAscii(std::string);
  void writeRelEntropyGradient(std::string);
  void writeRelEntropyGradientAscii(std::string);
  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);

  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;
  arma::Mat<double> rel_entropy_1p;
  arma::Col<double> rel_entropy_pos_1p;
  arma::Mat<double> rel_entropy_grad_1p;

  double freq_rms;

private:
  double pseudocount;
  int M;              // number of sequences
  int N;              // number of positions
  int Q;              // amino acid alphabet size
  double M_effective; // effect number of sequences

  MSA* msa;

  arma::Col<double> aa_background_frequencies;
};

#endif
