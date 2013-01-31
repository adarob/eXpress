/**
 *  lengthdistribution.cpp
 *  express
 *
 *  Created by Adam Roberts on 1/30/13.
 *  Copyright 2013 Adam Roberts. All rights reserved.
 */

#include "lengthdistribution.h"
#include "main.h"
#include <numeric>
#include <boost/assign.hpp>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace std;

LengthDistribution::LengthDistribution(double alpha, size_t max_val,
                                       size_t prior_mu, size_t prior_sigma,
                                       size_t kernel_n, double kernel_p,
                                       size_t bin_size)
    : _hist(max_val/bin_size+1),
      _tot_mass(LOG_0),
      _sum(LOG_0),
      _min(max_val/bin_size),
      _bin_size(bin_size) {
  
  max_val = max_val/bin_size;
  kernel_n = kernel_n/bin_size;
  assert(kernel_n % 2 == 0);
        
  double tot = log(alpha);

  // Set to prior distribution
  if (prior_mu) {
    boost::math::normal norm(prior_mu/bin_size,
                             prior_sigma/(bin_size*bin_size));

    for (size_t i = 0; i <= max_val; ++i) {
      double norm_mass = boost::math::cdf(norm, i+0.5) -
                         boost::math::cdf(norm, i-0.5);
      double mass = LOG_EPSILON;
      if (norm_mass != 0) {
        mass = tot + log(norm_mass);
      }
      _hist[i] = mass;
      _sum = log_add(_sum, log((double)i)+mass);
      _tot_mass = log_add(_tot_mass, mass);
    }
  } else {
    _hist = vector<double>(max_val + 1, tot - log(max_val));
    _hist[0] = LOG_0;
    _sum = _hist[1] + log(max_val * (max_val + 1)) - log(2);
    _tot_mass = tot;
  }

  // Define kernel
  boost::math::binomial_distribution<double> binom(kernel_n, kernel_p);
  _kernel = vector<double>(kernel_n + 1);
  for (size_t i = 0; i <= kernel_n; i++) {
    _kernel[i] = log(boost::math::pdf(binom, i));
  }
}

LengthDistribution::LengthDistribution(string param_file_name,
                                       string length_type) :
    _bin_size(1) {
  ifstream infile (param_file_name.c_str());
  const size_t BUFF_SIZE = 99999;
  char line_buff[BUFF_SIZE];
  
  if (!infile.is_open()) {
    cerr << "ERROR: Unable to open paramater file '" << param_file_name
    << "'.\n";
    exit(1);
  }

  do {
    infile.getline (line_buff, BUFF_SIZE, '\n');
  } while (strncmp(line_buff + 1, length_type.c_str(), length_type.size()));

  infile.getline (line_buff, BUFF_SIZE, '\n');
  char *p = strtok(line_buff, "\t");
  size_t i = 0;
  
  _tot_mass = 0;
  _sum = 0;
  do {
    double val = strtod(p,NULL);
    _hist.push_back(log(val));
    _tot_mass += val;
    _sum += i*val;
    i++;
    p = strtok(NULL, "\t");
  } while (p);
  
  _tot_mass = log(_tot_mass);
  _sum = log(_sum);
  _min = max_val();;
}

size_t LengthDistribution::max_val() const {
  return (_hist.size()-1) * _bin_size;
}

size_t LengthDistribution::min_val() const {
  if (_min == _hist.size() - 1) {
    return 1;
  }
  return _min;
}

void LengthDistribution::add_val(size_t len, double mass) {
  assert(!isnan(mass));
  assert(_kernel.size());
  
  len /= _bin_size;

  if (len > max_val()) {
      len = max_val();
  }
  if (len < _min) {
    _min = len;
  }

  size_t offset = len - _kernel.size()/2;

  for (size_t i = 0; i < _kernel.size(); i++) {
    if (offset > 0 && offset <= max_val()) {
      double k_mass = mass + _kernel[i];
      _hist[offset] = log_add(_hist[offset], k_mass);
      _sum = log_add(_sum, log((double)offset)+k_mass);
      _tot_mass = log_add(_tot_mass, k_mass);
    }
    offset++;
  }
}

double LengthDistribution::pmf(size_t len) const {
  len /= _bin_size;
  if (len > max_val()) {
    len = max_val();
  }
  return _hist[len]-_tot_mass;
}

vector<double> LengthDistribution::cmf() const {
  double cum = LOG_0;
  vector<double> cdf(_hist.size());
  for (size_t i = 0; i < _hist.size(); ++i) {
    cum = log_add(cum, _hist[i]);
    cdf[i] = cum - _tot_mass;
  }
  assert(approx_eq(cum, _tot_mass));

  return cdf;
}

double LengthDistribution::tot_mass() const {
  return _tot_mass;
}

double LengthDistribution::mean() const {
  return _sum - tot_mass();
}

string LengthDistribution::to_string() const {
  string s = "";
  char buffer[50];
  for(size_t i = 0; i <= max_val()*_bin_size; i++) {
    sprintf(buffer, "%e\t",sexp(pmf(i/_bin_size)/_bin_size));
    s += buffer;
  }
  s.erase(s.length()-1,1);
  return s;
}

void LengthDistribution::append_output(ofstream& outfile,
                                       string length_type) const {
  outfile << ">" << length_type << " Length Distribution (0-" << max_val()*_bin_size;
  outfile << ")\n" << to_string() << endl;
}
