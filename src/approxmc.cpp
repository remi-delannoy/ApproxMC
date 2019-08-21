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

AppMC::AppMC(AppMCConfig conf, SATSolver *solver) : conf(conf), solver(solver) {
  if (solver == NULL) {
    cerr << "No solver was defined" << endl;
    exit(1);
  }
  random_engine.seed(conf.seed);
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
  if (!conf.logfilename.empty()) {
    logfile.open(conf.logfilename.c_str());
    if (!logfile.is_open()) {
      cout << "[appmc] Cannot open AppMC log file '" << conf.logfilename
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
  // NOTE: is useless atm
  T min = num_list[0];
  for (auto it = num_list.begin() + 1; it != num_list.end(); ++it) {
    if ((*it) < min) {
      min = *it;
    }
  }
  return min;
}

void AppMC::addHash(uint32_t num_hashes, vector<Lit> &assumps,
                    uint32_t total_num_hashes) {
  std::uniform_real_distribution<uint32_t> dist{0.0, 1.0};
  double cutoff = 0.5;
  if (conf.sparse && total_num_hashes > 132) {
    // NOTE: magic numbers are related to yash's work
    cutoff = 13.46 * std::log(total_num_hashes) / total_num_hashes;
    assert(cutoff < 0.5);

    cout << "[appmc] sparse hashing used, cutoff: " << cutoff << endl;
  }
  vector<uint32_t> vars;
  for (uint32_t i = 0; i < num_hashes; i++) {
    // new activation variable, use to control the number of active hashes
    solver->new_var();
    uint32_t act_var = solver->nVars() - 1;
    assumps.push_back(Lit(act_var, false));

    vars.clear();
    vars.push_back(act_var);

    for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
      if (dist(random_engine) < cutoff) {
        vars.push_back(conf.sampling_set[j]);
      }
    }
    bool rhs = dist(random_engine) < 0.5;
    solver->add_xor_clause(vars, rhs);
    if (conf.verb_appmc_cls) {
      printXor(vars, rhs);
    }
  }
}

uint64_t AppMC::boundedSolCount(uint32_t max_solutions,
                                const vector<Lit> &assumps,
                                const uint32_t hash_count) {
  cout << "[appmc] "
          "[ "
       << std::setw(7) << std::setprecision(2) << std::fixed
       << (cpuTimeTotal() - total_runtime) << " ]"
       << " bounded_sol_count looking for " << std::setw(4) << max_solutions
       << " solutions"
       << " -- hashes active: " << hash_count << endl;

  // Set up things for adding clauses that can later be removed
  vector<lbool> model;
  vector<Lit> new_assumps(assumps);
  solver->new_var();
  uint32_t act_var = solver->nVars() - 1;
  new_assumps.push_back(Lit(act_var, false));
  if (hash_count > 2) {
    solver->simplify(&new_assumps);
  }

  uint64_t solutions = 0;
  lbool ret = l_True;
  double last_found_time = cpuTimeTotal();
  while (solutions < max_solutions && ret == l_True) {
    ret = solver->solve(&new_assumps);
    assert(ret == l_False || ret == l_True);

    if (conf.verb >= 2) {
      cout << "[appmc] bounded_sol_count ret: " << std::setw(7) << ret;
      if (ret == l_True) {
        cout << " sol no.  " << std::setw(3) << solutions;
      } else {
        cout << " No more. " << std::setw(3) << "";
      }
      cout << " T: " << std::setw(7) << std::setprecision(2) << std::fixed
           << (cpuTimeTotal() - total_runtime)
           << " -- hashes act: " << hash_count
           << " -- T since last: " << std::setw(7) << std::setprecision(2)
           << std::fixed << (cpuTimeTotal() - last_found_time) << endl;
      last_found_time = cpuTimeTotal();
    }

    if (ret == l_True) {
      model = solver->get_model();
      solutions++;

      vector<Lit> lits;
      lits.push_back(Lit(act_var, true));
      for (const uint32_t var : conf.sampling_set) {
        assert(model[var] != l_Undef);
        lits.push_back(
            Lit(var, model[var] == l_True)); // NOTE: shouldn't it be != true ?
      }
      if (conf.verb_appmc_cls) {
        cout << "[appmc] Adding banning clause: " << lits << endl;
      }
      solver->add_clause(lits);
    }
  }

  // Remove clauses added
  vector<Lit> cl_that_removes;
  cl_that_removes.push_back(Lit(act_var, true));
  solver->add_clause(cl_that_removes);

  assert(ret != l_Undef);
  return solutions;
}

// int AppMC::solve() {
//   openLogFile();
//   random_engine.seed(conf.seed);
//   total_runtime = cpuTimeTotal();
//   cout << "[appmc] Using start iteration " << conf.start_iter << endl;
//
//   bool finished = approxCount(solCount);
//   assert(finished);
//
//   cout << "[appmc] FINISHED AppMC T: " << (cpuTimeTotal() - start_time) << "
//   s"
//        << endl;
//   if (solCount.hash_count == 0 && solCount.cell_sol_count == 0) {
//     cout << "[appmc] Formula was UNSAT " << endl;
//   }
//
//   if (conf.verb > 2) {
//     solver->print_stats();
//   }
//
//   cout << "[appmc] Number of solutions is: " << solCount.cell_sol_count
//        << " x 2^" << solCount.hash_count << endl;
//
//   return correctReturnValue(l_True);
// }

void AppMC::setHash(uint32_t num_hashes, vector<Lit> &hash_vars,
                    vector<Lit> &assumps) {
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
      addHash(num_hashes - hash_vars.size(), assumps, num_hashes);
      for (size_t i = hash_vars.size(); i < num_hashes; i++) {
        hash_vars.push_back(assumps[i]);
      }
    }
  }
}

uint64_t AppMC::approxCount(double epsilon, double delta) {
  count.clear();
  vector<uint64_t> num_count_list;
  vector<Lit> assumps;
  vector<Lit> hash_vars; // assumption var to XOR hash
  assert(epsilon > 0);
  assert(delta > 0 && delta < 1);
  uint32_t threshold =
      uint32_t(1 + 9.84 * (1 + (1 / epsilon)) * (1 + (1 / epsilon)) *
                       (1 + (epsilon / (1 + epsilon))));
  uint32_t measurements = (int)std::ceil(std::log2(3.0 / conf.delta) * 17);
  uint32_t hash_count = conf.start_iter;

  double myTime = cpuTimeTotal();

  // Note, the rank of a random NxN matrix is not N of course. It has an
  // expected rank that is of course lower than N. So we need to shoot
  // higher.
  // https://math.stackexchange.com/questions/324150/expected-rank-of-a-random-binary-matrix
  // Apparently this question is analyzed in Kolchin's book Random Graphs
  // in sect. 3.2. Thanks to Yash Pote to digging this one out. Very
  // helpful.
  uint64_t total_max_xors =
      std::ceil((double)conf.sampling_set.size() * 1.2) + 5;

  for (uint32_t j = 0; j < measurements; j++) {
    hash_vars.clear();
    assumps.clear();
    if (conf.sparse) {
      // NOTE: for the sparse approach we need to add all the hashes first to
      // have the correct probability
      setHash(total_max_xors, hash_vars, assumps);
    }

    uint64_t lower_bound = 0;
    uint64_t upper_bound = total_max_xors;
    uint64_t best_num_solution;
    // NOTE: we don"t reset hash_count : we start from the previous optimal
    // hash_count
    setHash(hash_count, hash_vars, assumps);
    uint64_t current_num_solutions =
        boundedSolCount(threshold + 1, assumps, hash_count);
    uint64_t jump = 1;
    // exponential search in O(ln(|last_optimal-new_optimal|))
    if (current_num_solutions <= threshold) {
      while (current_num_solutions <= treshold &&
             upper_bound - jump > lower_bound) {
        best_num_solution = current_num_solutions;
        upper_bound = hash_count;
        hash_count -= jump;
        jump *= 2;
        setHash(hash_count, hash_vars, assumps);
        uint64_t current_num_solutions =
            boundedSolCount(threshold + 1, assumps, hash_count);
      }
      if (current_num_solutions > treshold) {
        lower_bound = hash_count;
      }
    } else {
      while (current_num_solutions > treshold &&
             lower_bound + jump < upper_bound) {
        lower_bound = hash_count;
        hash_count += jump;
        jump *= 2;
        setHash(hash_count, hash_vars, assumps);
        uint64_t current_num_solutions =
            boundedSolCount(threshold + 1, assumps, hash_count);
      }
      if (current_num_solutions <= treshold) {
        upper_bound = hash_count;
        best_num_solution = current_num_solutions;
      }
    }

    while (lower_bound < upper_bound - 1) {
      cout << "[appmc] gap to explore: " << std::setw(4)
           << upper_bound - lower_bound << " ind set size: " << std::setw(6)
           << conf.sampling_set.size() << endl;
      myTime = cpuTimeTotal();
      hash_count = (lower_bound + upper_bound) / 2;
      cout << "[appmc] hashes active: " << std::setw(6) << hash_count << endl;
      if (!conf.logfilename.empty()) {
        logfile << "appmc:" << j << ":" << hash_count << ":" << std::fixed
                << std::setprecision(2) << (cpuTimeTotal() - myTime) << ":"
                << (int)(current_num_solutions == (threshold + 1)) << ":"
                << current_num_solutions << endl;
      }
      setHash(hash_count, hash_vars, assumps);
      uint64_t current_num_solutions =
          boundedSolCount(threshold + 1, assumps, hash_count);
      if (current_num_solutions <= treshold) {
        best_num_solution = current_num_solutions;
        upper_bound = hash_count;
      } else {
        lower_bound = hash_count;
      }
    }
    num_count_list.push_back(best_num_solution);
  }

  return findMedian(num_count_list);
}

///////////
// Useful helper functions
///////////

void printVersionInfoAppMC() {
  cout << "c AppMC SHA revision " << ::get_version_sha1() << endl;
  cout << "c AppMC version " << ::get_version_tag() << endl;
  cout << "c AppMC compilation env " << ::get_compilation_env() << endl;
#ifdef __GNUC__
  cout << "c AppMC compiled with gcc version " << __VERSION__ << endl;
#else
  cout << "c AppMC compiled with non-gcc compiler" << endl;
#endif
}

void AppMC::printVersionInfo() const {
  ::printVersionInfoAppMC();
  cout << solver->get_text_version_info();
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
