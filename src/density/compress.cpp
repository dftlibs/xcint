#include <math.h> // fabs

#include "compress.h"
#include "balboa.h"
#include "density_parameters.h" // FIXME

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

void compress(balboa_context_t *balboa_context,
              const bool use_gradient,
              int &aoc_num,
              int *(&aoc_index),
              double *(&aoc),
              double ao[],
              const std::vector<int> &coor)
{
    std::vector<int> cent;
    for (size_t j = 0; j < coor.size(); j++)
    {
        cent.push_back((coor[j] - 1) / 3);
    }

    int num_slices;
    (use_gradient) ? (num_slices = 4) : (num_slices = 1);

    int n = 0;
    for (int i = 0; i < balboa_get_num_aos(balboa_context); i++)
    {
        if (is_same_center(balboa_get_ao_center(balboa_context, i), cent))
        {
            double tmax = 0.0;
            for (int ib = 0; ib < AO_BLOCK_LENGTH; ib++)
            {
                double t = fabs(ao[i * AO_BLOCK_LENGTH + ib]);
                if (t > tmax)
                    tmax = t;
            }
            if (tmax > 1.0e-15)
            {
                aoc_index[n] = i;
                n++;
            }
        }
    }
    aoc_num = n;

    if (aoc_num == 0)
        return;

    int kp[3] = {0, 0, 0};
    for (size_t j = 0; j < coor.size(); j++)
    {
        kp[(coor[j] - 1) % 3]++;
    }

    int off[4];
    off[0] = balboa_get_geo_offset(balboa_context, kp[0], kp[1], kp[2]);
    off[1] = balboa_get_geo_offset(balboa_context, kp[0] + 1, kp[1], kp[2]);
    off[2] = balboa_get_geo_offset(balboa_context, kp[0], kp[1] + 1, kp[2]);
    off[3] = balboa_get_geo_offset(balboa_context, kp[0], kp[1], kp[2] + 1);

    for (int i = 0; i < aoc_num; i++)
    {
        for (int islice = 0; islice < num_slices; islice++)
        {
            int iuoff = off[islice];
            int icoff = islice * balboa_get_num_aos(balboa_context);

            int iu = AO_BLOCK_LENGTH * (iuoff + aoc_index[i]);
            int ic = AO_BLOCK_LENGTH * (icoff + i);

            std::copy(&ao[iu], &ao[iu + AO_BLOCK_LENGTH], &aoc[ic]);
        }
    }
}
