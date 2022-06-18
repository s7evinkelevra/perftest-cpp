//
// Created by Jan on 25.05.2022.
//

#include "CSVWriter.h"

#include <utility>

void CSVWriter::flush() {
  _file.flush();
}