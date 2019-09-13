/*
 ApproxMC

 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#ifndef AppMC_H_
#define AppMC_H_

#include "approxmcconfig.h"
#include <cryptominisat5/cryptominisat.h>
#include <cstdint>
#include <fstream>
#include <map>
#include <random>

using std::string;
using std::vector;
using namespace CMSat;

class AppMC {
public:
  AppMC(AppMCConfig conf, SATSolver *solver);

  ~AppMC() {}

  static void printVersionInfo(SATSolver *solver);
  void printVersionInfo() const;
  std::pair<uint64_t, uint32_t> approxCount(double epsilon, double delta);
  // The number of solutions can be very large thus we return it in the form
  // cell_size*2^num_hash
  std::pair<uint64_t, uint32_t>
  approxCountWithAssumptions(double epsilon, double delta,
                             std::vector<Lit> const &assumps);

private:
  void readInAFile(SATSolver *solver2, const string &filename);
  void readInStandardInput(SATSolver *solver2);
  void openLogFile();
  void addHash(vector<Lit> &assumps, double xor_density);
  void setHash(uint32_t num_hashes, std::vector<Lit> &hash_vars,
               vector<Lit> &assumps, double pivot);
  int correctReturnValue(const lbool ret) const;

  uint64_t boundedSolCount(uint32_t max_solutions, const vector<Lit> &assumps,
                           const uint32_t hash_count);

  double start_time_;
  double total_runtime_;
  std::ofstream logfile_;
  std::mt19937 random_engine_;
  AppMCConfig conf_;
  SATSolver *solver_;
};

#endif // AppMC_H_
