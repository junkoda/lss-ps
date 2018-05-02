#ifndef ERROR_H
#define ERROR_H

class Error {};

class RuntimeError: Error {};

class AssertionError: RuntimeError{};

class IOError: RuntimeError {};

class FileNotFoundError: IOError {};

class MemoryError: RuntimeError {};

class ValError: RuntimeError {}; // check Python error name

class TypeError: RuntimeError {};

#endif
