#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include "msg.h"

namespace {
  enum LogLevel log_level;
  char prefix[8]= "# ";
}


void msg_set_loglevel(const enum LogLevel lv)
{
  log_level= lv;
}

LogLevel msg_get_loglevel()
{
  return log_level;
}


void msg_set_prefix(const char prefix_[])
{
  // Start messages with given prefix such as '#'
  strncpy(prefix, prefix_, 7);
}


void msg_printf(const enum LogLevel msg_level, char const * const fmt, ...)
{
  if(msg_level >= msg_error || msg_level >= log_level) {
    va_list argp;

    va_start(argp, fmt);

    fprintf(stdout, "%s", prefix);
    vfprintf(stdout, fmt, argp);
    fflush(stdout);

    va_end(argp);
  }
}


void msg_abort(char const * const fmt, ...)
{
  va_list argp;

  if(log_level <= msg_fatal) {
    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp); fflush(stdout);
    va_end(argp);
  }

  abort();
}
