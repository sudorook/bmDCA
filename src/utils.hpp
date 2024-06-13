/*
 * SPDX-FileCopyrightText: 2020 sudorook <daemon@nullcodon.com>
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

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

#include <armadillo>
#include <cstdint>
#include <string>

constexpr int BINS = 201;

/**
 * @brief Structure for storing 1D histograms.
 */
typedef struct histogram1d
{
  arma::Mat<uint64_t> range; ///< histogram counts
  int bins = BINS;           ///< number of bins in histogram
  double bin_width;          ///< width of bins
  double min = 0;            ///< minimum value
  double max = 1;            ///< maximum value
} histogram1d;

/**
 * @brief Structure for storing 2D histograms.
 */
typedef struct histogram2d
{
  arma::Mat<uint64_t> grid; ///< histogram counts
  int bins = BINS;          ///< number of bins in histogram
  double bin_width;         ///< width of bins
  double min = 0;           ///< minimum value
  double max = 1;           ///< maximum value
} histogram2d;

/**
 * @brief Structure for parameters of linear fit.
 */
typedef struct
{
  double a;  ///< intercept
  double b;  ///< slope
  double R2; ///< r^2
} linear_model;

/**
 * @brief Wrapper class for storing a FASTA sequence.
 */
class SeqRecord
{
private:
  const std::string header;   ///< sequence header
  const std::string sequence; ///< sequence

public:
  SeqRecord(std::string, std::string);
  void print();
  std::string getRecord();
  std::string getHeader();
  std::string getSequence();
};

/**
 * @brief Structure for storing Potts models.
 */
typedef struct
{
  arma::field<arma::Mat<double>> J; ///< couplings
  arma::Mat<double> h;              ///< fields
} potts_model;

potts_model loadPottsModel(std::string, std::string);

potts_model loadPottsModelAscii(std::string);

void convertFrequencyToAscii(std::string);

void convertParametersToAscii(std::string, std::string);

/**
 * @brief Step function.
 *
 * @param x input value
 *
 * @return 1 or 0
 */
template<typename T, typename U>
U
Theta(T x)
{
  if (x > 0)
    return 1;
  return 0;
}

/**
 * @brief Delta function.
 *
 * @param x input value
 *
 * @return 1 or 0
 */
template<typename T, typename U>
U
Delta(T x)
{
  if (x == 0)
    return 1;
  return 0;
}

/**
 * @brief Return larger of two values
 *
 * @param a first value
 * @param b second value
 *
 * @return larger of a or b
 */
template<typename T>
T
Max(T a, T b)
{
  return (a > b) ? a : b;
}

/**
 * @brief Return smaller of two values
 *
 * @param a first value
 * @param b second value
 *
 * @return smaller of a or b
 */
template<typename T>
T
Min(T a, T b)
{
  return (a < b) ? a : b;
}

int deleteFile(std::string);

bool checkFileExists(std::string);

void deleteAllFiles(std::string);

char
convertAA(int);

void
writeHistogram1D(std::string, histogram1d&);

void
writeHistogram2D(std::string, histogram2d&);

void
writeLinearModel(std::string, linear_model&);

#endif // SRC_UTILS_HPP_
