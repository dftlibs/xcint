#include "Functional.h"
#include "xcfun.h"

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <cstdio>
#include <cstring>


Functional::Functional() { nullify(); }

Functional::~Functional()
{
    delete[] functional_line;
    nullify();
}

void Functional::nullify()
{
    is_gga = false;
    is_tau_mgga = false;
    is_synced = false;
    keys.clear();
    weights.clear();
    functional_line = NULL;
}

void Functional::set_functional(const char *line)
{
    int ierr;

    parse(line);

    delete[] functional_line;
    functional_line = NULL;
    functional_line = new char[strlen(line) + 1];
    for (int i = 0; i < strlen(line); i++)
        functional_line[i] = line[i];
    functional_line[strlen(line)] = '\0';

    xc_functional fun;
    fun = xc_new_functional();
    for (int i = 0; i < keys.size(); i++)
    {
        ierr = xc_set(fun, keys[i].c_str(), weights[i]);
        if (ierr != 0)
        {
            fprintf(stderr,
                    "ERROR in fun init: \"%s\" not recognized, quitting.\n",
                    keys[i].c_str());
            exit(-1);
        }
    }
    xc_free_functional(fun);
}

void Functional::parse(const char *line)
{
    int pos;
    double w;
    std::string key;

    keys.clear();
    weights.clear();

    std::istringstream iss(line);

    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string> >(tokens));

    is_gga = false;
    is_tau_mgga = false;

    for (int i = 0; i < tokens.size(); i++)
    {
        pos = tokens[i].find("=");
        if (pos != std::string::npos)
        {
            key = tokens[i].substr(0, pos);
            w = atof(tokens[i].substr(pos + 1).c_str());
        }
        else
        {
            key = tokens[i];
            w = 1.0;
        }

        // convert to lowercase
        std::transform(key.begin(), key.end(), key.begin(), ::tolower);

        if (key == "lda")
        {
            keys.push_back("slaterx");
            weights.push_back(w);
            keys.push_back("vwn5c");
            weights.push_back(w);
        }

        if (key == "blyp")
        {
            keys.push_back("beckex");
            weights.push_back(w);
            keys.push_back("lypc");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "camb3lyp")
        {
            keys.push_back("camb3lyp");
            weights.push_back(w * 1.0);
            is_gga = true;
        }

        if (key == "b3lyp")
        {
            keys.push_back("slaterx");
            weights.push_back(w * 0.8);
            keys.push_back("beckecorrx");
            weights.push_back(w * 0.72);
            keys.push_back("vwn5c");
            weights.push_back(w * 0.19);
            keys.push_back("lypc");
            weights.push_back(w * 0.81);
            is_gga = true;
        }

        if (key == "pbe")
        {
            keys.push_back("pbex");
            weights.push_back(w);
            keys.push_back("pbec");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "pbe0")
        {
            keys.push_back("pbex");
            weights.push_back(w * 0.75);
            keys.push_back("pbec");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "m06")
        {
            keys.push_back("m06x");
            weights.push_back(w);
            keys.push_back("m06c");
            weights.push_back(w);
            is_tau_mgga = true;
        }

        if (key == "slaterx")
        {
            keys.push_back("slaterx");
            weights.push_back(w);
        }

        if (key == "pw86x")
        {
            keys.push_back("pw86x");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "vwn5c")
        {
            keys.push_back("vwn5c");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "pbex")
        {
            keys.push_back("pbex");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "lypc")
        {
            keys.push_back("lypc");
            weights.push_back(w);
            is_gga = true;
        }

        if (key == "beckecorrx")
        {
            keys.push_back("beckecorrx");
            weights.push_back(w);
            is_gga = true;
        }
    }

    if (keys.size() == 0)
    {
        fprintf(stderr, "ERROR: functional '%s' not recognized\n", line);
        fprintf(stderr, "       list of implemented functionals:\n");
        fprintf(stderr, "       b3lyp     \n");
        fprintf(stderr, "       beckecorrx\n");
        fprintf(stderr, "       blyp      \n");
        fprintf(stderr, "       camb3lyp  \n");
        fprintf(stderr, "       lda       \n");
        fprintf(stderr, "       lypc      \n");
        fprintf(stderr, "       m06       \n");
        fprintf(stderr, "       pbe       \n");
        fprintf(stderr, "       pbe0      \n");
        fprintf(stderr, "       pbex      \n");
        fprintf(stderr, "       pw86x     \n");
        fprintf(stderr, "       slaterx   \n");
        fprintf(stderr, "       vwn5c     \n");
        exit(-1);
    }
}

int Functional::set_order(const int order, xc_functional fun) const
{
    int ierr = -1;

    int dens_offset = (int)pow(2, order);

    if (is_tau_mgga)
    {
        ierr = xc_eval_setup(fun, XC_N_NX_NY_NZ_TAUN, XC_CONTRACTED, order);
    }
    else if (is_gga)
    {
        ierr = xc_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    }
    else
    {
        ierr = xc_eval_setup(fun, XC_N, XC_CONTRACTED, order);
    }

    if (ierr != 0)
    {
        fprintf(stderr, "ERROR in set_order (called with order %i).\n", order);
        exit(-1);
    }

    return dens_offset;
}
