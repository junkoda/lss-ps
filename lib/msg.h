#ifndef MSG_H
#define MSG_H 1

enum LogLevel {msg_debug=0, msg_verbose, msg_info, msg_warn, msg_error, msg_fatal, msg_silent};


void msg_set_loglevel(const LogLevel level);
void msg_set_prefix(const char prefix[]);

void msg_printf(const LogLevel level, char const * const fmt, ...);
void msg_abort(char const * const fmt, ...);

#endif
