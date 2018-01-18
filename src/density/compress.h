#pragma once

#include <vector>

void compress(const bool use_gradient,
              const int block_length,
              int &num_compressed_aos,
              int compressed_aos_indices[],
              double compressed_aos[],
              const int num_aos,
              const double aos[],
              const int ao_centers[],
              const std::vector<int> &coor,
              const int slice_offsets[4]);
