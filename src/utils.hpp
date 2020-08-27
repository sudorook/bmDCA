#ifndef UTILS_HPP
#define UTILS_HPP

#include <armadillo>
#include <string>

class SeqRecord
{
private:
  const std::string header;
  const std::string sequence;

public:
  SeqRecord(std::string, std::string);
  void print();
  std::string getRecord();
  std::string getHeader();
  std::string getSequence();
};

typedef struct
{
  arma::field<arma::Mat<double>> J;
  arma::Mat<double> h;
} potts_model;

potts_model loadPottsModel(std::string, std::string);

potts_model loadPottsModelAscii(std::string);

void convertFrequencyToAscii(std::string);

void convertParametersToAscii(std::string, std::string);

int
Theta(double);

int
Delta(double);

double
Max(double, double);

double
Min(double, double);

int deleteFile(std::string);

bool checkFileExists(std::string);

void deleteAllFiles(std::string);

char
convertAA(int);

#endif
