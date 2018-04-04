#include "compress.h"

// size_t
#include <cstddef>

// abs
#include <cmath>

bool is_same_center(const int c, const std::vector<int> &carray)
{
    for (unsigned int i = 0; i < carray.size(); i++)
    {
        if (c != carray[i])
            return false;
    }
    // returns true if carray is empty (no derivatives)
    return true;
}

void compress(const bool use_gradient,
              const int block_length,
              int &num_compressed_aos,
              int compressed_aos_indices[],
              double compressed_aos[],
              const int num_aos,
              const double aos[],
              const int ao_centers[],
              const std::vector<int> &coor,
              const int slice_offsets[4])
{
    std::vector<int> cent;
    for (size_t j = 0; j < coor.size(); j++)
    {
        cent.push_back((coor[j] - 1) / 3);
    }

    int num_slices;
    (use_gradient) ? (num_slices = 4) : (num_slices = 1);

    int n = 0;
    for (int i = 0; i < num_aos; i++)
    {
        if (is_same_center(ao_centers[i], cent))
        {
            double tmax = 0.0;
            for (int ib = 0; ib < block_length; ib++)
            {
                double t = std::abs(aos[i * block_length + ib]);
                if (t > tmax)
                    tmax = t;
            }
            if (tmax > 1.0e-15)
            {
                compressed_aos_indices[n] = i;
                n++;
            }
        }
    }
    num_compressed_aos = n;

    if (num_compressed_aos == 0)
        return;

    for (int i = 0; i < num_compressed_aos; i++)
    {
        for (int islice = 0; islice < num_slices; islice++)
        {
            int iuoff = slice_offsets[islice];
            int icoff = islice * num_aos;

            int iu = block_length * (iuoff + compressed_aos_indices[i]);
            int ic = block_length * (icoff + i);

            std::copy(&aos[iu], &aos[iu + block_length], &compressed_aos[ic]);
        }
    }
}
