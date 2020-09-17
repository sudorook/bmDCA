#ifndef SGDM_HPP
#define SGDM_HPP

#include "model.hpp"
#include "utils.hpp"

class SGDM : public Model
{
public:
  SGDM(void);
  // virtual ~SGDM(void){};

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

  double lambda_reg_h = 0.01;
  double lambda_reg_J = 0.01;
  double alpha_reg = 1.0;
  std::string initial_params = "zero";
  bool set_zero_gauge = false;
  double learn_rate_h = 0.01;
  double learn_rate_J = 0.01;
  double beta_h = 0.9;
  double beta_J = 0.9;

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
  void writeGradient(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeMoment1Ascii(std::string);
  void writeGradientAscii(std::string);
};

#endif
