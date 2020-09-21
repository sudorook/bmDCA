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
Model::setZeroGauge(void)
{
  potts_model params_zg;

  // initialize
  params_zg.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params_zg.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params_zg.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

   // rescale couplings
   for (int i = 0; i < N; i++) {
     for (int j = 0; j < N; j++) {
       if (i < j) {
         double J_ij_mean = arma::mean(arma::mean(params.J(i, j)));
         arma::Mat<double> J_ija_mean = arma::mean(params.J(i, j), 1);
         arma::Mat<double> J_ijb_mean = arma::mean(params.J(i, j), 0);
         for (int a = 0; a < Q; a++) {
           for (int b = 0; b < Q; b++) {
             params_zg.J(i, j)(a, b) =
               params.J(i, j)(a, b) - J_ija_mean(a) - J_ijb_mean(b) + J_ij_mean;
           }
           params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
         }
       } else if (i > j) {
         double J_ij_mean = arma::mean(arma::mean(params.J(j, i)));
         arma::Mat<double> J_ija_mean = arma::mean(params.J(j, i), 0);
         arma::Mat<double> J_ijb_mean = arma::mean(params.J(j, i), 1);
         for (int a = 0; a < Q; a++) {
           params_zg.h(a, i) += J_ija_mean(a) - J_ij_mean;
         }
       }
     }
   }

   // rescale fields
   arma::Row<double> h_i_mean = arma::mean(params.h, 0);
   for (int i = 0; i < N; i++) {
     for (int a = 0; a < Q; a++) {
       params_zg.h(a, i) += params.h(a, i) - h_i_mean(i);
     }
   }

   params.h = params_zg.h;
   params.J = params_zg.J;
};
