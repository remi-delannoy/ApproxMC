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

#include <boost/program_options.hpp>
using boost::lexical_cast;
namespace po = boost::program_options;
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

#include "approxmc.h"
#include "approxmcconfig.h"
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"
#include "time_mem.h"
#include <cryptominisat5/cryptominisat.h>
#include <memory>

using namespace CMSat;
using std::cerr;
using std::cout;
using std::endl;

po::options_description appmc_options =
    po::options_description("ApproxMC options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

void add_appmc_options(double epsilon, double delta, AppMCConfig &conf) {
  std::ostringstream my_epsilon;
  std::ostringstream my_delta;
  my_epsilon << std::setprecision(8) << epsilon;
  my_delta << std::setprecision(8) << delta;

  appmc_options.add_options()("help,h", "Prints help")(
      "version", "Print version info")("input", po::value<vector<string>>(),
                                       "file(s) to read")(
      "verb,v", po::value(&conf.verb)->default_value(conf.verb), "verbosity")(
      "seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")(
      "epsilon", po::value(&epsilon)->default_value(epsilon, my_epsilon.str()),
      "epsilon parameter as per PAC guarantees")(
      "delta", po::value(&delta)->default_value(delta, my_delta.str()),
      "delta parameter as per PAC guarantees; 1-delta is the confidence")(
      "start", po::value(&conf.start_iter)->default_value(conf.start_iter),
      "Start at this many XORs")(
      "log", po::value(&conf.logfilename)->default_value(conf.logfilename),
      "Logs of ApproxMC execution")(
      "th", po::value(&conf.num_threads)->default_value(conf.num_threads),
      "How many solving threads to use per solver call")(
      "vcl",
      po::value(&conf.verb_appmc_cls)->default_value(conf.verb_appmc_cls),
      "Print banning clause + xor clauses. Highly verbose.")(
      "sparse", po::value(&conf.sparse)->default_value(conf.sparse),
      "Generate sparse XORs when possible");

  help_options.add(appmc_options);
}

void add_supported_options(int argc, char **argv, AppMCConfig &conf,
                           SATSolver *solver) {
  add_appmc_options(0.8, 0.2, conf);
  p.add("input", 1);

  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(help_options)
                  .positional(p)
                  .run(),
              vm);
    if (vm.count("help")) {
      cout << "Probably Approximate counter" << endl;

      cout << "approxmc [options] inputfile" << endl << endl;

      cout << help_options << endl;
      std::exit(0);
    }

    if (vm.count("version")) {
      AppMC::printVersionInfo(solver);
      std::exit(0);
    }

    po::notify(vm);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<po::unknown_option>>
               &c) {
    cerr << "ERROR: Some option you gave was wrong. Please give '--help' to "
            "get help"
         << endl
         << "       Unkown option: " << c.what() << endl;
    std::exit(-1);
  } catch (boost::bad_any_cast &e) {
    std::cerr << "ERROR! You probably gave a wrong argument type" << endl
              << "       Bad cast: " << e.what() << endl;

    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<
               po::invalid_option_value>> &what) {
    cerr << "ERROR: Invalid value '" << what.what() << "'" << endl
         << "       given to option '" << what.get_option_name() << "'" << endl;

    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<
               po::multiple_occurrences>> &what) {
    cerr << "ERROR: " << what.what() << " of option '" << what.get_option_name()
         << "'" << endl;

    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<po::required_option>>
               &what) {
    cerr << "ERROR: You forgot to give a required option '"
         << what.get_option_name() << "'" << endl;

    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<
               po::too_many_positional_options_error>> &what) {
    cerr << "ERROR: You gave too many positional arguments. Only the input CNF "
            "can be given as a positional option."
         << endl;
    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<po::ambiguous_option>>
               &what) {
    cerr << "ERROR: The option you gave was not fully written and matches"
         << endl
         << "       more than one option. Please give the full option name."
         << endl
         << "       The option you gave: '" << what.get_option_name() << "'"
         << endl
         << "       The alternatives are: ";
    for (size_t i = 0; i < what.alternatives().size(); i++) {
      cout << what.alternatives()[i];
      if (i + 1 < what.alternatives().size()) {
        cout << ", ";
      }
    }
    cout << endl;

    std::exit(-1);
  } catch (boost::exception_detail::clone_impl<
           boost::exception_detail::error_info_injector<
               po::invalid_command_line_syntax>> &what) {
    cerr << "ERROR: The option you gave is missing the argument or the" << endl
         << "       argument is given with space between the equal sign."
         << endl
         << "       detailed error message: " << what.what() << endl;
    std::exit(-1);
  }
}

void readInAFile(SATSolver *solver, const string &filename, AppMCConfig &conf) {
  solver->add_sql_tag("filename", filename);
#ifndef USE_ZLIB
  FILE *in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE *, FN>> parser(solver, NULL, 2);
#else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, GZ>> parser(solver, NULL, 2);
#endif

  if (in == NULL) {
    std::cerr << "ERROR! Could not open file '" << filename
              << "' for reading: " << strerror(errno) << endl;

    std::exit(-1);
  }

  if (!parser.parse_DIMACS(in, false)) {
    exit(-1);
  }

  conf.sampling_set.swap(parser.sampling_vars);

#ifndef USE_ZLIB
  fclose(in);
#else
  gzclose(in);
#endif
}

void readInStandardInput(SATSolver *solver, AppMCConfig &conf) {
  cout << "c Reading from standard input... Use '-h' or '--help' for help."
       << endl;

#ifndef USE_ZLIB
  FILE *in = stdin;
#else
  gzFile in = gzdopen(0, "rb"); // opens stdin, which is 0
#endif

  if (in == NULL) {
    std::cerr << "ERROR! Could not open standard input for reading" << endl;
    std::exit(1);
  }

#ifndef USE_ZLIB
  DimacsParser<StreamBuffer<FILE *, FN>> parser(solver, NULL, 2);
#else
  DimacsParser<StreamBuffer<gzFile, GZ>> parser(solver, NULL, 2);
#endif

  if (!parser.parse_DIMACS(in, false)) {
    exit(-1);
  }

  conf.sampling_set.swap(parser.sampling_vars);

#ifdef USE_ZLIB
  gzclose(in);
#endif
}

void set_sampling_vars(SATSolver *solver, AppMCConfig &conf) {
  if (conf.sampling_set.empty()) {
    cout << "[appmc] WARNING! Sampling set was not declared with 'c ind var1 "
            "[var2 var3 ..] 0'"
            " notation in the CNF."
         << endl
         << "[appmc] we may work substantially worse!" << endl;
    for (size_t i = 0; i < solver->nVars(); i++) {
      conf.sampling_set.push_back(i);
    }
  }
  cout << "[appmc] Sampling set size: " << conf.sampling_set.size() << endl;
  if (conf.sampling_set.size() > 100) {
    cout << "[appmc] Sampling var set contains over 100 variables, not "
            "displaying"
         << endl;
  } else {
    cout << "[appmc] Sampling set: ";
    for (auto v : conf.sampling_set) {
      cout << v + 1 << ", ";
    }
    cout << endl;
  }
  solver->set_sampling_vars(&conf.sampling_set);
}

int main(int argc, char **argv) {
#if defined(__GNUC__) && defined(__linux__)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  std::unique_ptr<SATSolver> solver(new SATSolver());
  AppMCConfig conf;
  add_supported_options(argc, argv, conf, solver.get());
  cout << "[appmc] using seed: " << conf.seed << endl;

  if (vm.count("log") == 0) {
    if (vm.count("input") != 0) {
      conf.logfilename = vm["input"].as<vector<string>>()[0] + ".log";
      cout << "[appmc] Logfile name not given, assumed to be "
           << conf.logfilename << endl;
    } else {
      std::cerr << "[appmc] ERROR: You must provide the logfile name" << endl;
      exit(-1);
    }
  }

  // startTime = cpuTimeTotal();

  solver->set_up_for_scalmc();

  if (conf.verb > 2) {
    solver->set_verbosity(conf.verb - 2);
  }
  solver->set_allow_otf_gauss();

  if (conf.num_threads > 1) {
    solver->set_num_threads(conf.num_threads);
  }

  // parsing the input
  if (vm.count("input") != 0) {
    vector<string> inp = vm["input"].as<vector<string>>();
    if (inp.size() > 1) {
      cout << "[appmc] ERROR: can only parse in one file" << endl;
    }
    readInAFile(solver.get(), inp[0].c_str(), conf);
  } else {
    readInStandardInput(solver.get(), conf);
  }
  set_sampling_vars(solver.get(), conf);

  if (conf.start_iter > conf.sampling_set.size()) {
    cout << "[appmc] ERROR: Manually-specified start_iter"
            "is larger than the size of the sampling set.\n"
         << endl;
    exit(-1);
  }

  AppMC appmc(conf, solver.get());
  appmc.printVersionInfo();
  auto count = appmc.approxCount(0.8, 0.2);
  cout << "[appmc] Number of solutions is: " << count.first << "*2^"
       << count.second << endl;
  return 0;
}
