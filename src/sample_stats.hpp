/*
 * SPDX-FileCopyrightText: 2020 - 2022 sudorook <daemon@nullcodon.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SRC_SAMPLE_STATS_HPP_
#define SRC_SAMPLE_STATS_HPP_

#include <armadillo>
#include <string>

#include "utils.hpp"

/**
 * @brief Abstract class for computing statistics for sampled sequences.
 *
 * Because 1p/2p statistics are more efficiently computed with different data
 * structures depending on whether or not multiple sequences are sampled along
 * trajectories, the SampleStats abstract class acts as a wrapper around the
 * key functions. SampleStats2D or SampleStats3D can be recasted to
 * SampleStats.
 */
class SampleStats
{
public:
  SampleStats(void){};
  virtual ~SampleStats(void){};

  virtual void writeStep(int,
                         bool = true) = 0; ///< write all stats for the step
  virtual void writeData(std::string,
                         bool = true) = 0;    ///< write all stats for the step
  virtual void writeSamples(std::string) = 0; ///< write the samples to disk
  virtual void writeSampleEnergies(
    std::string) = 0;                  ///< write the sample energies to disk
  virtual void computeStats(void) = 0; ///< compute the 1p/2p statistics
  virtual void computeStatsExtra(
    void) = 0; ///< compute the energies and sequence correlations
  virtual void computeStatsImportance(
    void) = 0; ///< run an importance sampling to re-estimate parameters
  virtual void setMixingTime(
    int) = 0; ///< set the burn-between time for the samples
  virtual arma::Col<double> getStats(void) = 0; ///< return important metrics

  arma::Mat<double> frequency_1p;              ///< 1p frequencies
  arma::field<arma::Mat<double>> frequency_2p; ///< 2p sample frequencies

protected:
  int mixing_time = 1; ///< burn-between time used by sampler
};

/**
 * @brief Class to compute sample statistics (independent).
 */
class SampleStats2D : public SampleStats
{
public:
  SampleStats2D(arma::Mat<int>*, potts_model*, potts_model* = nullptr);

  void computeStats(void) override;
  void computeStatsExtra(void) override;
  void writeStep(int, bool = true) override;
  void writeData(std::string, bool = true) override;
  void computeStatsImportance(void) override;
  void setMixingTime(int) override;
  arma::Col<double> getStats(void) override;

private:
  void computeEnergies(void);
  void computeCorrelations(void);
  void computeSampleStats(void);
  void computeSampleStatsImportance(void);

  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  void writeFrequency1pAscii(std::string);
  void writeFrequency2pAscii(std::string);

  void writeSamples(std::string) override;
  void writeSampleEnergies(std::string) override;

  double Z_ratio;   ///< for importance sampling...
  double sumw_inv;  ///< for importance sampling...
  double dE_av_tot; ///< for importance sampling...

  double overlap_inf;       ///< total sequence correlations
  double overlap_inf_sigma; ///< variance of correlations

  potts_model* params; ///< Potts parameters
  potts_model*
    params_prev; ///< previous Potts parameters (for importance sampling)
  arma::Mat<int>* samples;    ///< pointer to sampled sequences
  arma::Col<double> energies; ///< vector of energies for sampled sequences

  int N; ///< number of positions
  int Q; ///< number of states
  int M; ///< number of sampled sequences
};

/**
 * @brief Class to compute sample statistics (MCMC).
 */
class SampleStats3D : public SampleStats
{
public:
  SampleStats3D(arma::Cube<int>*, potts_model*, potts_model* = nullptr);

  void computeStats(void) override;
  void computeStatsExtra(void) override;
  void writeStep(int, bool = true) override;
  void writeData(std::string, bool = true) override;
  void computeStatsImportance(void) override;
  void setMixingTime(int) override;
  arma::Col<double> getStats(void) override;

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

  void writeSamples(std::string) override;
  void writeSampleEnergies(std::string) override;
  void writeSampleEnergiesRelaxation(std::string);

  arma::Mat<double>
    frequency_1p_sigma; ///< variance of 1p estimates across trajectories
  arma::field<arma::Mat<double>>
    frequency_2p_sigma; ///< variance of 1p estimates across trajectories

  arma::Row<double>
    energies_relax; ///< mean energies at each step along trajectories
  arma::Row<double> energies_relax_sigma; ///< variance of energies at each step
                                          ///< along trajectories

  arma::Col<double> overlaps;       ///< sequence correlations
  arma::Col<double> overlaps_sigma; ///< variance in sequence correlations

  double energies_start_avg;   ///< average energy after burn-in
  double energies_start_sigma; ///< variance of average energy after burn-in
  double energies_end_avg;     ///< average energy at trajectory end
  double energies_end_sigma;   ///< variance of average energy at trajectory end
  double energies_err;         ///< combined variance of energy at start and end

  double Z_ratio;   ///< for importance sampling...
  double sumw_inv;  ///< for importance sampling...
  double dE_av_tot; ///< for importance sampling...

  double overlap_inf;       ///< total sequence correlations
  double overlap_inf_sigma; ///< variance in total correlations
  double overlap_auto;      ///< correlation at start and end of trajectories
  double overlap_cross;     ///< correlations between trajectories
  double overlap_check;     ///< correlatino at middle and end of trajectories
  double sigma_auto;        ///< variance of auto-correlation
  double sigma_cross;       ///< variance of cross-correlation
  double sigma_check;       ///< variance of mid-trajectory and end correlation
  double err_cross_auto;  ///< combined between-trajectory and within-trajectory
                          ///< variance (start)
  double err_cross_check; ///< combined between-trajectory and within-trajectory
                          ///< variance (middle)
  double err_check_auto;  ///< combined within-trajectory (start) and
                          ///< within-trajectory variance (middke)

  potts_model* params; ///< Potts parameters
  potts_model*
    params_prev; ///< previous Potts parameters (for importance sampling)
  arma::Cube<int>* samples;   ///< pointer to sampled sequences
  arma::Col<double> energies; ///< vector of energies for sampled sequences

  int reps; ///< number of trajectories
  int N;    ///< number of positions
  int Q;    ///< number of states
  int M;    ///< number of sampled sequences per trajectory
};

#endif // SRC_SAMPLE_STATS_HPP_
