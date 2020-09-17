#ifndef ADAMW_HPP
#define ADAMW_HPP

#include "model.hpp"
#include "utils.hpp"

class AdamW : public Model
{
public:
  AdamW(void);
  // virtual ~AdamW(void){};

  void update(void);
  void initialize(void);
  void reset(void);
  void restore(int, bool = true);

  void writeData(std::string, bool = true);
  void writeStep(int, bool = true);
  void deleteStep(int, bool = true);
  
  void loadHyperparameters(std::string);
  void writeHyperparameters(std::string, bool = true);
  bool isValidStep(int, bool = true);

private:
  potts_model gradient;
  potts_model moment1;
  potts_model moment2;

  double lambda_decay_h = 0.01;
  double lambda_decay_J = 0.01;
  std::string initial_params = "zero";
  bool set_zero_gauge = false;
  double learn_rate_h = 0.01;
  double learn_rate_J = 0.01;

  double eta_min = 0.1;
  double eta_max = 0.1;
  std::string anneal_schedule = "none";
  double anneal_scale = 2.0;
  int anneal_period = 40;
  int anneal_warm = 20;
  int anneal_hot = 0;
  int anneal_cool = 0;

  void updateGradients(void);
  void updateMoments(void);
  void updateParameters(void);

  bool compareHyperparameters(std::string);
  bool compareHyperparameter(std::string, std::string);
  void checkHyperparameters(void);
  void setHyperparameter(std::string, std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeMoment1(std::string, std::string);
  void writeMoment2(std::string, std::string);
  void writeGradient(std::string, std::string);
  
  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeMoment1Ascii(std::string);
  void writeMoment2Ascii(std::string);
  void writeGradientAscii(std::string);
};

#endif
