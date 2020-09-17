#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "sample_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  Model(void);
  virtual ~Model(){};

  void setMSAStats(MSAStats*, MSAStats*);
  void setSampleStats(SampleStats*);
  void setStep(int);

  double train_error_1p = 1000;
  double train_error_2p = 1000;
  double validation_error_1p = 1000;
  double validation_error_2p = 1000;

  virtual void update(void) = 0;
  virtual void initialize(void) = 0;
  virtual void reset(void) = 0;
  virtual void restore(int, bool) = 0;

  virtual void writeData(std::string, bool = true) = 0;
  virtual void writeStep(int, bool = true) = 0;
  virtual void deleteStep(int, bool = true) = 0;
  virtual bool isValidStep(int, bool = true) = 0;

  virtual void loadHyperparameters(std::string) = 0;
  virtual void checkHyperparameters(void) = 0;
  virtual bool compareHyperparameters(std::string) = 0;
  virtual void writeHyperparameters(std::string, bool = false) = 0;

  potts_model params;
  potts_model params_prev;

protected:
  MSAStats* training = nullptr;
  MSAStats* validation = nullptr;
  SampleStats* samples = nullptr;

  std::string hyperparameter_file;

  bool validate;

  void setZeroGauge(void);

  int N;
  int Q;

  double temperature = 1.0;

  int step = 1;
};

#endif
