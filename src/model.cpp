#include "model.hpp"

#include <armadillo>

#include "utils.hpp"

Model::Model(void){};

void
Model::setMSAStats(MSAStats* msa_train, MSAStats* msa_validate)
{
  training = msa_train;
  validation = msa_validate;

  N = training->getN();
  Q = training->getQ();

  if (msa_validate) {
    validate = true;
  } else {
    validate = false;
  }
};

void
Model::setSampleStats(SampleStats* s)
{
  samples = s;
};

void
Model::setStep(int s)
{
  step = s;
};

void
Model::setZeroGauge(void){};
