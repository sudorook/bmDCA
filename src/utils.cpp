#include "utils.hpp"

#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>

SeqRecord::SeqRecord(std::string h, std::string s)
  : header(h)
  , sequence(s){};

void
SeqRecord::print(void)
{
  std::cout << ">" << header << std::endl;
  std::cout << sequence << std::endl;
};

std::string
SeqRecord::getRecord(void)
{
  std::string record_string = ">" + header + "\n" + sequence + "\n";
  return record_string;
};

std::string
SeqRecord::getHeader(void)
{
  return header;
};

std::string
SeqRecord::getSequence(void)
{
  return sequence;
};

potts_model
loadPottsModel(std::string h_file, std::string J_file)
{
  potts_model params;
  params.h.load(h_file);
  params.J.load(J_file);
  return params;
};

potts_model
loadPottsModelAscii(std::string parameters_file)
{
  std::ifstream input_stream(parameters_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << parameters_file << "' for reading."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int N = 0;
  int Q = 0;

  int count = 1;
  int n1, n2, aa1, aa2;
  double value;
  std::string tmp = "";
  std::getline(input_stream, tmp);
  while (std::getline(input_stream, tmp)) {
    input_stream >> tmp;
    input_stream >> n1 >> n2 >> aa1 >> aa2;
    input_stream >> value;
    count++;

    if ((n1 == 1) & (N == 0)) {
      N = count;
    }
    if ((aa2 == 0) & (Q == 0)) {
      Q = count;
    }
    if ((N != 0) & (Q != 0))
      break;
  }
  N = (int)((double)N / (double)Q / double(Q)) + 1;

  input_stream.clear();
  input_stream.seekg(0);

  potts_model params;
  params.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  for (int count = 0; count < (int)N * (N - 1) / 2 * Q * Q; count++) {
    input_stream >> tmp;
    input_stream >> n1 >> n2 >> aa1 >> aa2;
    input_stream >> value;
    params.J(n1, n2)(aa1, aa2) = value;
  }

  for (int count = 0; count < N * Q; count++) {
    input_stream >> tmp;
    input_stream >> n1 >> aa1;
    input_stream >> value;
    params.h(aa1, n1) = value;
  }
  return params;
};

void
convertFrequencyToAscii(std::string stats_file)
{
  int idx = stats_file.find_last_of(".");
  std::string stats_name = stats_file.substr(0, idx);
  std::string stats_ext = stats_file.substr(idx + 1);

  if (stats_ext != "bin") {
    std::cerr << "ERROR: input file does not have 'bin' extension."
              << std::endl;
    return;
  }

  // Guess if 1p vs 2p statistics
  bool is_1p = false;
  bool is_2p = false;
  if (stats_name.find("_1p") != std::string::npos) {
    is_1p = true;
  } else if (stats_name.find("_2p") != std::string::npos) {
    is_2p = true;
  }

  std::string output_file = stats_name + ".txt";
  std::ofstream output_stream(output_file);

  if (is_1p) {
    arma::Mat<double> frequency_1p;
    frequency_1p.load(stats_file, arma::arma_binary);

    int N = frequency_1p.n_cols;
    int Q = frequency_1p.n_rows;

    for (int i = 0; i < N; i++) {
      output_stream << i;
      for (int aa = 0; aa < Q; aa++) {
        output_stream << " " << frequency_1p(aa, i);
      }
      output_stream << std::endl;
    }
  } else if (is_2p) {
    arma::field<arma::Mat<double>> frequency_2p;
    frequency_2p.load(stats_file, arma::arma_binary);

    int N = frequency_2p.n_rows;
    int Q = frequency_2p(0, 1).n_rows;

    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        output_stream << i << " " << j;
        for (int aa1 = 0; aa1 < Q; aa1++) {
          for (int aa2 = 0; aa2 < Q; aa2++) {
            output_stream << " " << frequency_2p(i, j)(aa1, aa2);
          }
        }
        output_stream << std::endl;
      }
    }
  } else { // if name doesn't say if 1p or 2p, load files and guess again
    arma::Mat<double> frequency_1p;
    frequency_1p.load(stats_file, arma::arma_binary);

    int N = frequency_1p.n_cols;
    int Q = frequency_1p.n_rows;

    if ((Q != 0) & (N != 0)) { // 1p
      for (int i = 0; i < N; i++) {
        output_stream << i;
        for (int aa = 0; aa < Q; aa++) {
          output_stream << " " << frequency_1p(aa, i);
        }
        output_stream << std::endl;
      }
    } else { // 2p
      arma::field<arma::Mat<double>> frequency_2p;
      frequency_2p.load(stats_file, arma::arma_binary);

      int N = frequency_2p.n_rows;
      int Q = frequency_2p(0, 1).n_rows;

      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          output_stream << i << " " << j;
          for (int aa1 = 0; aa1 < Q; aa1++) {
            for (int aa2 = 0; aa2 < Q; aa2++) {
              output_stream << " " << frequency_2p(i, j)(aa1, aa2);
            }
          }
          output_stream << std::endl;
        }
      }
    }
  }
};

void
convertParametersToAscii(std::string h_file, std::string J_file)
{

  // Check file extensions and parse out file names.
  int idx = h_file.find_last_of(".");
  std::string h_name = h_file.substr(0, idx);
  std::string h_ext = h_file.substr(idx + 1);

  idx = J_file.find_last_of(".");
  std::string J_name = J_file.substr(0, idx);
  std::string J_ext = J_file.substr(idx + 1);

  if ((J_ext != "bin") & (h_ext != "bin")) {
    std::cerr << "ERROR: input parameters do not have 'bin' extension."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  arma::Mat<double> h;

  h.load(h_file, arma::arma_binary);
  int N = h.n_cols;
  int Q = h.n_rows;

  arma::field<arma::Mat<double>> J(N, N);
  J.load(J_file, arma::arma_binary);

  if ((N != (int)J.n_rows) & (N != (int)J.n_cols)) {
    std::cerr << "ERROR: parameters N dimension mismatch." << std::endl;
    return;
  }
  if ((Q != (int)J(0, 1).n_cols) & (Q != (int)J(0, 1).n_rows)) {
    std::cerr << "ERROR: parameters Q dimension mismatch." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Generate an output file name.
  std::string output_file;
  for (int i = 0; i < Min(h_name.size(), J_name.size()); i++) {
    if (h_name[i] == J_name[i]) {
      if ((output_file.back() == '_') && (h_name[i] == '_'))
        continue;
      output_file += h_name[i];
    }
  }
  if (output_file.back() == '_')
    output_file.pop_back();
  std::ofstream output_stream(output_file + ".txt");

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << h(aa, i) << std::endl;
    }
  }
  return;
};

int
Theta(double x)
{
  if (x > 0)
    return 1;
  return 0;
};

int
Delta(double x)
{
  if (x == 0)
    return 1;
  return 0;
};

double
Max(double a, double b)
{
  if (a > b)
    return a;
  return b;
};

double
Min(double a, double b)
{
  if (a < b)
    return a;
  return b;
};

int
deleteFile(std::string filename)
{
  std::fstream fs;
  fs.open(filename);
  if (!fs.fail()) {
    if (std::remove(filename.c_str()) != 0) {
      fs.close();
      return 1;
    }
  }
  fs.close();
  return 0;
};

bool
checkFileExists(std::string filename)
{
  std::fstream fs;
  fs.open(filename);
  if (fs.fail()) {
    fs.close();
    return false;
  } else {
    fs.close();
    return true;
  }
};

void
deleteAllFiles(std::string directory)
{
  DIR* dp;
  struct dirent* dirp;

  dp = opendir(".");
  std::vector<int> steps;
  while ((dirp = readdir(dp)) != NULL) {
    std::string fname = dirp->d_name;
    if (fname == ".")
      continue;
    if (fname == "..")
      continue;
    if (std::remove(fname.c_str()) != 0) {
      std::cerr << "ERROR: deletion of '" << fname << "' failed." << std::endl;
    }
  }
  closedir(dp);
  return;
};

char
convertAA(int n)
{
  char aa = '\0';
  switch (n) {
    case 0:
      aa = '-';
      break;
    case 1:
      aa = 'A';
      break;
    case 2:
      aa = 'C';
      break;
    case 3:
      aa = 'D';
      break;
    case 4:
      aa = 'E';
      break;
    case 5:
      aa = 'F';
      break;
    case 6:
      aa = 'G';
      break;
    case 7:
      aa = 'H';
      break;
    case 8:
      aa = 'I';
      break;
    case 9:
      aa = 'K';
      break;
    case 10:
      aa = 'L';
      break;
    case 11:
      aa = 'M';
      break;
    case 12:
      aa = 'N';
      break;
    case 13:
      aa = 'P';
      break;
    case 14:
      aa = 'Q';
      break;
    case 15:
      aa = 'R';
      break;
    case 16:
      aa = 'S';
      break;
    case 17:
      aa = 'T';
      break;
    case 18:
      aa = 'V';
      break;
    case 19:
      aa = 'W';
      break;
    case 20:
      aa = 'Y';
      break;
  }
  return aa;
};

void
writeHistogram1D(std::string file, histogram1d hist)
{
  std::ofstream output_stream(file);
  int N = hist.range.n_elem;
  for (int i = 0; i < N; i++) {
    output_stream << hist.min + i * hist.bin_width << "\t" << hist.range(i)
                  << std::endl;
  }
  output_stream.close();
  return;
};

void
writeHistogram2D(std::string file, histogram2d hist)
{
  std::ofstream output_stream(file);
  int N1 = hist.grid.n_rows;
  int N2 = hist.grid.n_cols;
  for (int i = 0; i < N1; i++) {
    for (int j = 0; j < N2; j++) {
      output_stream << hist.min + i * hist.bin_width << "\t"
                    << hist.min + j * hist.bin_width << "\t" << hist.grid(i, j)
                    << std::endl;
    }
  }
  output_stream.close();
  return;
};

void
writeLinearModel(std::string file, linear_model model)
{
  std::ofstream output_stream(file);
  output_stream << "a"
                << "\t" << model.a << std::endl;
  output_stream << "b"
                << "\t" << model.b << std::endl;
  output_stream << "R2"
                << "\t" << model.R2 << std::endl;
  output_stream.close();
  return;
};
