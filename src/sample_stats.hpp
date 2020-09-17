#ifndef SAMPLE_STATS_HPP
#define SAMPLE_STATS_HPP

#include <armadillo>

#include "utils.hpp"

class SampleStats
{
public:
  SampleStats(void) {};
  virtual ~SampleStats(void){};

  virtual void writeStep(int, bool = true) = 0;
  virtual void writeData(std::string, bool = true) = 0;
  virtual void writeSamples(std::string) = 0;
  virtual void writeSampleEnergies(std::string) = 0;
  virtual void computeStats(void) = 0;
  virtual void computeStatsExtra(void) = 0;
  virtual void computeStatsImportance(void) = 0;
  virtual void setMixingTime(int) = 0;
  virtual arma::Col<double> getStats(void) = 0;

  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;

protected:
  int mixing_time = 1;
};

class SampleStats2D : public SampleStats
{
public:
  SampleStats2D(arma::Mat<int>*, potts_model*, potts_model* = nullptr);

  void computeStats(void);
  void computeStatsExtra(void);
  void writeStep(int, bool=true);
  void writeData(std::string, bool=true);
  void computeStatsImportance(void);
  void setMixingTime(int);
  arma::Col<double> getStats(void);

private:
  void computeEnergies(void);
  void computeCorrelations(void);
  void computeSampleStats(void);
  void computeSampleStatsImportance(void);

  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);

  void writeSamples(std::string);
  void writeSampleEnergies(std::string);

  double Z_ratio;
  double sumw_inv;
  double dE_av_tot;

  double overlap_inf;
  double overlap_inf_sigma;

  potts_model *params;
  potts_model *params_prev;
  arma::Mat<int>* samples;
  arma::Col<double> energies;

  int N;
  int Q;
  int M;
};

class SampleStats3D : public SampleStats
{
public:
  SampleStats3D(arma::Cube<int>*, potts_model*, potts_model* = nullptr);

  void computeStats(void);
  void computeStatsExtra(void);
  void writeStep(int, bool=true);
  void writeData(std::string, bool=true);
  void computeStatsImportance(void);
  void setMixingTime(int);
  arma::Col<double> getStats(void);

private:
  void computeEnergies(void);
  void computeEnergiesStats(void);
  void computeCorrelations(void);
  void computeSampleStats(void);
  void computeSampleStatsImportance(void);

  void writeFrequency1p(std::string, std::string);
  void writeFrequency2p(std::string, std::string);
  void writeFrequency1pAscii(std::string, std::string);
  void writeFrequency2pAscii(std::string, std::string);

  void writeSamples(std::string);
  void writeSampleEnergies(std::string);
  void writeSampleEnergiesRelaxation(std::string);
  
  arma::Mat<double> frequency_1p_sigma;
  arma::field<arma::Mat<double>> frequency_2p_sigma;

  arma::Row<double> energies_relax;
  arma::Row<double> energies_relax_sigma;

  arma::Col<double> overlaps;
  arma::Col<double> overlaps_sigma;

  double energies_start_avg;
  double energies_start_sigma;
  double energies_end_avg;
  double energies_end_sigma;
  double energies_err;

  double Z_ratio;
  double sumw_inv;
  double dE_av_tot;

  double overlap_inf;
  double overlap_inf_sigma;
  double overlap_auto;
  double overlap_cross;
  double overlap_check;
  double sigma_auto;
  double sigma_cross;
  double sigma_check;
  double err_cross_auto;
  double err_cross_check;
  double err_check_auto;

  potts_model *params;
  potts_model *params_prev;
  arma::Cube<int>* samples;
  arma::Mat<double> energies;

  int reps;
  int N;
  int Q;
  int M;
};

#endif
