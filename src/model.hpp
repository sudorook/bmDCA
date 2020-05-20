#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  potts_model params;
  potts_model params_prev;
  potts_model learning_rates;
  potts_model gradient;
  potts_model gradient_prev;
  int N;
  int Q;

  Model(MSAStats, double, double, bool=true);
  Model(std::string, std::string, std::string, std::string, std::string);
  Model(std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeLearningRates(std::string, std::string);
  void writeGradient(std::string, std::string);
  void writeGradientPrevious(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeLearningRatesAscii(std::string);
  void writeGradientAscii(std::string);
  void writeGradientPreviousAscii(std::string);
};

#endif
