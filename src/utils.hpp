#ifndef UTILS_HPP
#define UTILS_HPP

#include <armadillo>
#include <string>

#define BINS 201

typedef struct
{
  arma::Mat<unsigned long long int> range;
  int bins = BINS;
  double bin_width;
  double min = 0;
  double max = 1;
} histogram1d;

typedef struct
{
  arma::Mat<unsigned long long int> grid;
  int bins = BINS;
  double bin_width;
  double min = 0;
  double max = 1;
} histogram2d;

typedef struct
{
  double a;
  double b;
  double R2;
} linear_model;

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

void writeHistogram1D(std::string, histogram1d);

void writeHistogram2D(std::string, histogram2d);

void writeLinearModel(std::string, linear_model);

#endif
