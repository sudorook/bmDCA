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

#ifndef SRC_GENERATOR_HPP_
#define SRC_GENERATOR_HPP_

#include "pcg_random.hpp"
#include "sample_stats.hpp"
#include "sampler.hpp"
#include "utils.hpp"
#include <memory>
#include <string>

/**
 * @brief Generator class for samplign sequences from given Potts model.
 *
 * Samples pre-defined number of sequences along pre-defined number of
 * trajectories. Burn times can be fixed or adaptive, and sampling temperature
 * can also be specified.
 */
class Generator
{
public:
  Generator(potts_model, int, int, const std::string&);
  void run(int, int, const std::string&);
  void writeAASequences(const std::string&);
  void writeNumericalSequences(const std::string&);

private:
  int N = 0;                ///< number of positions
  int Q = 0;                ///< number of amino acids
  int samples_per_walk = 1; ///< number of independent sampling runs
  int walkers = 1000; ///< number of sequences sampled from independent runs

  int resample_max = 20;         ///< max number of times to resample sequences
  unsigned random_seed = 1;      ///< random seed for the PCG RNG
  double adapt_up_time = 1.5;    ///< scaling factor for increasing burn times
  double adapt_down_time = 0.6;  ///< scaling factor for decreasing burn times
  int burn_in_start = 100000;    ///< starting burn-in time
  int burn_between_start = 1000; ///< starting burn-between time
  bool update_burn_time =
    true; ///< flag to update burn-time or stick to starting values
  bool save_interim_samples = false; ///< write sequence between resamplings
  std::string update_rule =
    "mh";                   ///< update rule for accepting proposed samples
  double temperature = 1.0; ///< temperature for sampling

  int burn_in = 100000;    ///< stores current burn-in time
  int burn_between = 1000; ///< stores current burn-between time

  pcg32 rng; ///< the RNG for the sampler

  arma::Mat<int>
    samples_2d; ///< efficiently store sequences when samples_per_run=1
  arma::Cube<int>
    samples_3d;      ///< efficiently store sequences when samples_per_run>1
  potts_model model; ///< structure for Potts model parameters

  std::shared_ptr<Sampler> sampler; ///< sampler
  std::shared_ptr<SampleStats>
    sample_stats; ///< pointer for sampled sequence statistics

  void estimateBurnTime(void);
  bool checkErgodicity(void);
  void loadParameters(const std::string&);
  void checkParameters(void);
  void setParameter(const std::string&, const std::string&);
};

#endif // SRC_GENERATOR_HPP_
