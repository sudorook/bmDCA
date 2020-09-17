#ifndef ORIGINAL_HPP
#define ORIGINAL_HPP

#include "model.hpp"
#include "utils.hpp"

class Original : public Model
{
public:
  Original(void);
  // virtual ~Original(void){};

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
  potts_model gradient_prev;
  potts_model learning_rates;

  double lambda_reg_h = 0.01;
  double lambda_reg_J = 0.01;
  double alpha_reg = 1.0;
  std::string initial_params = "zero";
  bool set_zero_gauge = false;
  double epsilon_h = 0.01;
  double epsilon_J = 0.001;
  double learn_rate_h_min = 1e-04;
  double learn_rate_h_max = 0.5;
  double learn_rate_J_min = 1e-04;
  double learn_rate_J_max = 0.5;
  double adapt_up = 1.5;
  double adapt_down = 0.6;
  bool use_pos_reg = false;

  void updateGradients(void);
  void updateLearningRates(void);
  void updateParameters(void);

  bool compareHyperparameters(std::string);
  bool compareHyperparameter(std::string, std::string);
  void checkHyperparameters(void);
  void setHyperparameter(std::string, std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeGradient(std::string, std::string);
  void writeGradientPrevious(std::string, std::string);
  void writeLearningRates(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeGradientAscii(std::string);
  void writeGradientPreviousAscii(std::string);
  void writeLearningRatesAscii(std::string);
};

#endif
