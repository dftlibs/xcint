#pragma once

#include <vector>
#include "balboa.h"

void compress(balboa_context_t *balboa_context,
              const bool use_gradient,
              int &aoc_num,
              int aoc_index[],
              double aoc[],
              double ao[],
              const std::vector<int> &coor);
