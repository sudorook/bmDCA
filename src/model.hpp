#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  potts_model params;
  potts_model params_prev;
  potts_model gradient;
  potts_model gradient_prev;
  potts_model moment1;
  potts_model moment2;
  int N;
  int Q;

  Model(MSAStats*, bool = true);
  Model(std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string);
  Model(std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string);

  void resetModel(MSAStats*, bool = true);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeMoment1(std::string, std::string);
  void writeMoment2(std::string, std::string);
  void writeGradient(std::string, std::string);
  void writeGradientPrevious(std::string, std::string);

  void writeParamsAscii(std::string);
  void writeParamsPreviousAscii(std::string);
  void writeMoment1Ascii(std::string);
  void writeMoment2Ascii(std::string);
  void writeGradientAscii(std::string);
  void writeGradientPreviousAscii(std::string);
};

#endif
