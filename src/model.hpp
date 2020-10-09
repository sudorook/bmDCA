#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "sample_stats.hpp"
#include "utils.hpp"

/**
 * @brief Abstract Model class.
 *
 * The inference scheme for learning parameters is not known at compile time.
 * Instead, it's set via the config file at runtime. This Model class defines a
 * set of vitrual functions to be shared amongst all classes derived from this
 * one. A common function table allows an instance of any derived class
 * specified in the config file, recast to a Model, to be loaded at runtime.
 */
class Model
{
public:
  Model(void);
  virtual ~Model(){};

  void setMSAStats(MSAStats*, MSAStats*);
  void setSampleStats(SampleStats*);
  void setStep(int);

  double train_error_1p = 1000;      ///< training RMSE for fields
  double train_error_2p = 1000;      ///< training RMSE for couplings
  double validation_error_1p = 1000; ///< validation RMSE for fields
  double validation_error_2p = 1000; ///< validation RMSE for couplings

  // Define virtual functions that form the Model 'API'
  virtual void update(void) = 0;       ///< update parameters, gradients, etc.
  virtual void initialize(void) = 0;   ///< initialize the model
  virtual void reset(void) = 0;        ///< reset the model to initial values
  virtual void restore(int, bool) = 0; ///< re-load parameters from disk

  virtual void writeData(std::string, bool = true) = 0; ///< write current model to disk
  virtual void writeStep(int, bool = true) = 0;         ///< write all data to reload current step
  virtual void deleteStep(int, bool = true) = 0;        ///< delete data for a step
  virtual bool isValidStep(int, bool = true) = 0;       ///< check if data exists to reload a step

  virtual void loadHyperparameters(std::string) = 0;                ///< read parameters from config file
  virtual void checkHyperparameters(void) = 0;                      ///< check config parameters against saved state
  virtual bool compareHyperparameters(std::string) = 0;             ///< check that config parameters are valid
  virtual void writeHyperparameters(std::string, bool = false) = 0; ///< write config file

  potts_model params;      ///< current step parameters
  potts_model params_prev; ///< previous step parameters

protected:
  MSAStats* training = nullptr;   ///< pointer to training MSA stats
  MSAStats* validation = nullptr; ///< pointer to validation MSA stats
  SampleStats* samples = nullptr; ///< pointer to sampled sequence stats

  std::string hyperparameter_file; ///< file string for storing hyperparams

  bool validate; ///< flag for whether to compute validation error

  void setZeroGauge(void);

  int N; ///< number of positions
  int Q; ///< number of states

  double temperature = 1.0; ///< inference temperature (unused)

  int step = 1; ///< iteration number
};

#endif
