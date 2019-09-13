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

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <ctime>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string.h>
#include <sys/stat.h>

#include "GitSHA1.h"
#include "approxmc.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "time_mem.h"

using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::map;

AppMC::AppMC(AppMCConfig conf, SATSolver *solver)
    : conf_(conf), solver_(solver) {
  if (solver_ == NULL) {
    cerr << "No solver was defined" << endl;
    exit(1);
  }
  random_engine_.seed(conf_.seed);
}

void printXor(const vector<uint32_t> &vars, const uint32_t rhs) {
  cout << "[appmc] Added XOR ";
  for (size_t i = 0; i < vars.size(); i++) {
    cout << vars[i] + 1;
    if (i < vars.size() - 1) {
      cout << " + ";
    }
  }
  cout << " = " << (rhs ? "True" : "False") << endl;
}

void AppMC::openLogFile() {
  if (!conf_.logfilename.empty()) {
    logfile_.open(conf_.logfilename.c_str());
    if (!logfile_.is_open()) {
      cout << "[appmc] Cannot open AppMC log file '" << conf_.logfilename
           << "' for writing." << endl;
      exit(1);
    }
  }
}

template <class T> inline T findMedian(vector<T> &num_list) {
  size_t med_index = num_list.size() / 2;
  std::nth_element(num_list.begin(), num_list.begin() + med_index,
                   num_list.end());
  return num_list[med_index];
}

template <class T> inline T findMin(vector<T> &num_list) {
  T min = num_list[0];
  for (auto it = num_list.begin() + 1; it != num_list.end(); ++it) {
    if ((*it) < min) {
      min = *it;
    }
  }
  return min;
}

double inverseBinEntropy(double x, double precision = 0.01) {
  double a = 0;
  double b = 0.5;
  double p = (b + a) / 2;
  while ((b - a) > precision) {
    p = (b + a) / 2;
    if ((-p * log2(p) - (1 - p) * log2(1 - p)) < x) {
      a = p;
    } else {
      b = p;
    }
  }
  assert(p != 0);
  return p;
}

void AppMC::addHash(vector<Lit> &assumps, double xor_density) {
  std::uniform_real_distribution<double> dist{0.0, 1.0};
  vector<uint32_t> vars;
  // new activation variable, use to control the number of active hashes
  solver_->new_var();
  uint32_t act_var = solver_->nVars() - 1;
  assumps.push_back(Lit(act_var, true));
  vars.clear();
  vars.push_back(act_var);
  for (uint32_t j = 0; j < conf_.sampling_set.size(); j++) {
    if (dist(random_engine_) < xor_density) {
      vars.push_back(conf_.sampling_set[j]);
    }
  }
  bool rhs = dist(random_engine_) < 0.5;
  solver_->add_xor_clause(vars, rhs);
  if (conf_.verb_appmc_cls) {
    printXor(vars, rhs);
  }
}

void AppMC::setHash(uint32_t num_hashes, vector<Lit> &hash_vars,
                    vector<Lit> &assumps, double pivot) {
  if (num_hashes < assumps.size()) {
    uint64_t numberToRemove = assumps.size() - num_hashes;
    for (uint64_t i = 0; i < numberToRemove; i++) {
      assumps.pop_back();
    }
  } else {
    if (num_hashes > assumps.size() && assumps.size() < hash_vars.size()) {
      for (size_t i = assumps.size(); i < hash_vars.size() && i < num_hashes;
           i++) {
        assumps.push_back(hash_vars[i]);
      }
    }
    if (num_hashes > hash_vars.size()) {
      for (size_t i = hash_vars.size(); i < num_hashes; i++) {
        double xor_density = 0.5;
        if (conf_.sparse) {
          // i+1 because we start at 1 and not 0
          double entropy = (i + 1) / (i + 1 + log2(pivot));
          xor_density = std::min(0.5, 16 / inverseBinEntropy(entropy) *
                                          log2(i + 1) / (i + 1));
        }
        assert(xor_density <= 0.5);
        addHash(assumps, xor_density);
        hash_vars.push_back(assumps[i]);
      }
    }
  }
}

uint64_t AppMC::boundedSolCount(uint32_t max_solutions,
                                const vector<Lit> &assumps,
                                const uint32_t hash_count) {
  if (conf_.verb >= 2) {
    cout << "[appmc] "
            "[ "
         << std::setw(7) << std::setprecision(2) << std::fixed
         << (cpuTimeTotal() - total_runtime_) << " ]"
         << " bounded_sol_count looking for " << std::setw(4) << max_solutions
         << " solutions"
         << " -- hashes active: " << hash_count << endl;
  }
  // Set up things for adding clauses that can later be removed
  vector<lbool> model;
  vector<Lit> new_assumps(assumps);
  solver_->new_var();
  uint32_t act_var = solver_->nVars() - 1;
  new_assumps.push_back(Lit(act_var, true));
  if (hash_count > 2) {
    solver_->simplify(&new_assumps);
  }

  uint64_t solutions = 0;
  lbool ret = l_True;
  double last_found_time = cpuTimeTotal();
  while (solutions < max_solutions && ret == l_True) {
    ret = solver_->solve(&new_assumps);
    assert(ret == l_False || ret == l_True);

    if (conf_.verb >= 2) {
      cout << "[appmc] bounded_sol_count ret: " << std::setw(7) << ret;
      if (ret == l_True) {
        cout << " sol no.  " << std::setw(3) << solutions;
      } else {
        cout << " No more. " << std::setw(3) << "";
      }
      cout << " T: " << std::setw(7) << std::setprecision(2) << std::fixed
           << (cpuTimeTotal() - total_runtime_)
           << " -- hashes act: " << hash_count
           << " -- T since last: " << std::setw(7) << std::setprecision(2)
           << std::fixed << (cpuTimeTotal() - last_found_time) << endl;
      last_found_time = cpuTimeTotal();
    }

    if (ret == l_True) {
      model = solver_->get_model();
      solutions++;

      vector<Lit> lits;
      lits.push_back(Lit(act_var, false));
      for (const uint32_t var : conf_.sampling_set) {
        assert(model[var] != l_Undef);
        lits.push_back(
            Lit(var, model[var] == l_True)); // NOTE: adding Lit(x,true) means x
                                             // is inverted
      }
      if (conf_.verb_appmc_cls) {
        cout << "[appmc] Adding banning clause: " << lits << endl;
      }
      solver_->add_clause(lits);
    }
  }

  // Remove clauses added
  vector<Lit> cl_that_removes;
  cl_that_removes.push_back(Lit(act_var, false));
  solver_->add_clause(cl_that_removes);

  assert(ret != l_Undef);
  return solutions;
}

std::pair<uint64_t, uint32_t>
AppMC::approxCountWithAssumptions(double epsilon, double delta,
                                  std::vector<Lit> const &assumps) {
  vector<uint64_t> num_count_list;
  vector<uint32_t> num_hash_list;
  vector<Lit> xor_assumps;
  vector<Lit> full_assumps(assumps);
  vector<Lit> hash_vars; // assumption var to XOR hash
  assert(epsilon > 0);
  assert(delta > 0 && delta < 1);
  uint32_t threshold =
      uint32_t(1 + 9.84 * (1 + (1 / epsilon)) * (1 + (1 / epsilon)) *
                       (1 + (epsilon / (1 + epsilon))));
  uint32_t measurements = (int)std::ceil(std::log2(3.0 / delta) * 17);
  uint32_t hash_count = conf_.start_iter;

  double myTime = cpuTimeTotal();

  // Note, the rank of a random NxN matrix is not N of course. It has an
  // expected rank that is of course lower than N. So we need to shoot
  // higher.
  // https://math.stackexchange.com/questions/324150/expected-rank-of-a-random-binary-matrix
  // Apparently this question is analyzed in Kolchin's book Random Graphs
  // in sect. 3.2. Thanks to Yash Pote to digging this one out. Very
  // helpful.
  uint64_t total_max_xors =
      std::ceil((double)conf_.sampling_set.size() * 1.2) + 5;
  double pivot = 78.72 * conf_.ro * (1 + 1 / epsilon) *
                 (1 + 1 / epsilon); // used only if conf_.sparse
  for (uint32_t j = 0; j < measurements; j++) {
    hash_vars.clear();
    xor_assumps.clear();

    uint64_t lower_bound = 0;
    uint64_t upper_bound = total_max_xors;
    uint64_t best_num_solution = 0;
    uint32_t best_num_hash = 0;
    // NOTE: we don"t reset hash_count : we start from the previous optimal
    // hash_count
    setHash(hash_count, hash_vars, xor_assumps, pivot);
    full_assumps = assumps;
    full_assumps.insert(full_assumps.end(), xor_assumps.begin(),
                        xor_assumps.end());
    uint64_t current_num_solutions =
        boundedSolCount(threshold + 1, full_assumps, hash_count);
    uint64_t jump = 1;
    // exponential search in O(ln(|last_optimal-new_optimal|))
    if (current_num_solutions <= threshold) {
      upper_bound = hash_count;
      best_num_solution = current_num_solutions;
      best_num_hash = hash_count;
      while (current_num_solutions <= threshold &&
             upper_bound - lower_bound > jump) {
        best_num_solution = current_num_solutions;
        best_num_hash = hash_count;
        upper_bound = hash_count;
        hash_count -= jump;
        jump *= 2;
        setHash(hash_count, hash_vars, xor_assumps, pivot);
        full_assumps = assumps;
        full_assumps.insert(full_assumps.end(), xor_assumps.begin(),
                            xor_assumps.end());
        current_num_solutions =
            boundedSolCount(threshold + 1, full_assumps, hash_count);
      }
      if (current_num_solutions > threshold) {
        lower_bound = hash_count;
      }
    } else {
      lower_bound = hash_count;
      while (current_num_solutions > threshold &&
             upper_bound - lower_bound > jump) {
        lower_bound = hash_count;
        hash_count += jump;
        jump *= 2;
        setHash(hash_count, hash_vars, xor_assumps, pivot);
        full_assumps = assumps;
        full_assumps.insert(full_assumps.end(), xor_assumps.begin(),
                            xor_assumps.end());
        current_num_solutions =
            boundedSolCount(threshold + 1, full_assumps, hash_count);
      }
      if (current_num_solutions <= threshold) {
        upper_bound = hash_count;
        best_num_solution = current_num_solutions;
        best_num_hash = hash_count;
      }
    }

    while (upper_bound - lower_bound > 1) {
      myTime = cpuTimeTotal();
      hash_count = (lower_bound + upper_bound) / 2;
      if (!conf_.logfilename.empty()) {
        logfile_ << "appmc:" << j << ":" << hash_count << ":" << std::fixed
                 << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                 << (int)(current_num_solutions == (threshold + 1)) << ":"
                 << current_num_solutions << endl;
      }
      setHash(hash_count, hash_vars, xor_assumps, pivot);
      full_assumps = assumps;
      full_assumps.insert(full_assumps.end(), xor_assumps.begin(),
                          xor_assumps.end());
      current_num_solutions =
          boundedSolCount(threshold + 1, full_assumps, hash_count);
      if (current_num_solutions <= threshold) {
        best_num_solution = current_num_solutions;
        best_num_hash = hash_count;
        upper_bound = hash_count;
      } else {
        lower_bound = hash_count;
      }
    }
    num_count_list.push_back(best_num_solution);
    num_hash_list.push_back(best_num_hash);
  }
  uint32_t hash_min = findMin(num_hash_list);
  auto num_hash_it = num_hash_list.begin();
  for (auto num_sol_it = num_count_list.begin();
       num_sol_it != num_count_list.end() && num_hash_it != num_hash_list.end();
       ++num_sol_it, ++num_hash_it) {
    *num_sol_it *= pow(2, (*num_hash_it) - hash_min);
  }

  return std::make_pair(findMedian(num_count_list), hash_min);
}

std::pair<uint64_t, uint32_t> AppMC::approxCount(double epsilon, double delta) {
  vector<Lit> assumps;
  return approxCountWithAssumptions(epsilon, delta, assumps);
}

///////////
// Useful helper functions
///////////

static void printVersionInfoAppMC() {
  cout << "c AppMC SHA revision " << ::get_version_sha1() << endl;
  cout << "c AppMC version " << ::get_version_tag() << endl;
  cout << "c AppMC compilation env " << ::get_compilation_env() << endl;
#ifdef __GNUC__
  cout << "c AppMC compiled with gcc version " << __VERSION__ << endl;
#else
  cout << "c AppMC compiled with non-gcc compiler" << endl;
#endif
}

void AppMC::printVersionInfo(SATSolver *solver) {
  ::printVersionInfoAppMC();
  cout << solver->get_text_version_info();
}
void AppMC::printVersionInfo() const {
  ::printVersionInfoAppMC();
  cout << solver_->get_text_version_info();
}

int AppMC::correctReturnValue(const lbool ret) const {
  int retval = -1;
  if (ret == l_True) {
    retval = 10;
  } else if (ret == l_False) {
    retval = 20;
  } else if (ret == l_Undef) {
    retval = 15;
  } else {
    std::cerr << "Something is very wrong, output is neither l_Undef, nor "
                 "l_False, nor l_True"
              << endl;
    exit(-1);
  }

  return retval;
}
