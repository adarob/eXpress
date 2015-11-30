/**
 *  logger.cpp
 *  express
 *
 *  Created by Adam Roberts on 6/22/13.
 *  Copyright 2013 Adam Roberts. All rights reserved.
 **/

#ifndef __express__logger__
#define __express__logger__

#include <iostream>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <cstdarg>
#include <cstdio>

namespace pt = boost::posix_time;

class Logger {
  const static size_t BUFF_SIZE = 4096;
  
  std::ostream* _info_out;
  std::ostream* _warn_out;
  std::ostream* _severe_out;
  mutable boost::mutex _mut;
 
  std::string get_time() const {
    return pt::to_simple_string(pt::second_clock::local_time());
  }
  
public:
  Logger()
    : _info_out(&std::cerr), _warn_out(&std::cerr), _severe_out(&std::cerr) {}

  void info_out(std::ostream* out) {
    boost::unique_lock<boost::mutex>(_mut);
    _info_out = out;
  }

  void warn_out(std::ostream* out) {
    boost::unique_lock<boost::mutex>(_mut);
    _warn_out = out;
  }
  
  void severe_out(std::ostream* out) {
    boost::unique_lock<boost::mutex>(_mut);
    _severe_out = out;
  }
  
  void info(const char* msg, ...) const {
    boost::unique_lock<boost::mutex>(_mut);
    char buffer[BUFF_SIZE];
    std::va_list arg;
    va_start(arg, msg);
    vsnprintf(buffer, BUFF_SIZE, msg, arg);
    va_end(arg);
    *_info_out << get_time() << " - " << buffer << std::endl;
  }

  void warn(const char* msg, ...) const {
    boost::unique_lock<boost::mutex>(_mut);
    char buffer[BUFF_SIZE];
    std::va_list arg;
    va_start(arg, msg);
    vsnprintf(buffer, BUFF_SIZE, msg, arg);
    va_end(arg);
    *_warn_out << get_time() << " - WARNING: " << buffer << std::endl;
  }

  void severe(const char* msg, ...) const {
    boost::unique_lock<boost::mutex>(_mut);
    char buffer[BUFF_SIZE];
    std::va_list arg;
    va_start(arg, msg);
    vsnprintf(buffer, BUFF_SIZE, msg, arg);
    va_end(arg);
    *_severe_out << get_time() << " - SEVERE: " << buffer << std::endl;
    exit(1);
  }
};

#endif /* defined(__express__logger__) */
