
// has to be first include
#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include <math.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#include "io.h"
#include "Functional.h"


Functional::Functional()
{
    fun = xc_new_functional();

    nullify();
}


Functional::~Functional()
{
    xc_free_functional(fun);
    delete [] functional_line;
    nullify();
}


void Functional::nullify()
{
    is_gga = false;
    is_tau_mgga = false;
    is_synced = false;
    dens_offset = -1;
    keys.clear();
    weights.clear();
    functional_line = NULL;
}


void Functional::set_functional(const int verbosity, const char *line, double &hfx, double &mu, double &beta)
{
    int ierr;

    parse(line, hfx, mu, beta);

    functional_line = new char[strlen(line)+1];
    for (int i = 0; i < strlen(line); i++) functional_line[i] = line[i];
    functional_line[strlen(line)] = '\0';

    if (verbosity > 0)
    {
        io::speak_your_mind("\n\nXCint functional\n");
        io::speak_your_mind("----------------\n\n");
        io::speak_your_mind("XC functional: %s\n", line);
        io::speak_your_mind("Composition (XC part):\n");
    }

    for (int i = 0; i < keys.size(); i++)
    {
        ierr = xc_set(fun, keys[i].c_str(), weights[i]);
        if (ierr != 0)
        {
            fprintf(stderr, "ERROR in fun init: \"%s\" not recognized, quitting.\n", keys[i].c_str());
            exit(-1);
        }
        if (verbosity > 0)
        {
            io::speak_your_mind("        %s=%.4f\n", keys[i].c_str(), weights[i]);
        }
    }
}


#ifdef ENABLE_MPI
void Functional::sync_functional(const MPI_Comm &comm)
{
    if (is_synced) return;

    int rank = 0;
    int num_proc = 1;

    MPI_Comm_size(comm, &num_proc);
    MPI_Comm_rank(comm, &rank);

    if (num_proc > 1)
    {
        int l;
        if (rank == 0) l = strlen(functional_line);

        MPI_Bcast(&l, 1, MPI_INT, 0, comm);

        if (rank > 0)
        {
            functional_line = new char[l+1];
        }

        MPI_Bcast(functional_line, l, MPI_CHAR, 0, comm);

        if (rank > 0)
        {
            double hfx, mu, beta;
            functional_line[l] = '\0';
            set_functional(0, functional_line, hfx, mu, beta);
        }

        is_synced = true;
    }
}
#endif


void Functional::parse(const char *line,
                       double &hfx,
                       double &mu,
                       double &beta)
{
    int pos;
    double w;
    std::string key;

    std::istringstream iss(line);

    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string> >(tokens));

    is_gga = false;
    is_tau_mgga = false;

    hfx  = 0.0;
    mu   = 0.0;
    beta = 0.0;

    for (int i = 0; i < tokens.size(); i++)
    {
        pos = tokens[i].find("=");
        if (pos != std::string::npos)
        {
            key = tokens[i].substr(0, pos);
            w = atof(tokens[i].substr(pos+1).c_str());
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
            keys.push_back("beckesrx");
            weights.push_back(w);
            keys.push_back("vwn5c");
            weights.push_back(w*0.19);
            keys.push_back("lypc");
            weights.push_back(w*0.81);
            is_gga = true;
            hfx = w*0.19;
            mu = w*0.33;
            beta = w*0.46;
        }

        if (key == "b3lyp")
        {
            keys.push_back("slaterx");
            weights.push_back(w*0.8);
            keys.push_back("beckecorrx");
            weights.push_back(w*0.72);
            keys.push_back("vwn5c");
            weights.push_back(w*0.19);
            keys.push_back("lypc");
            weights.push_back(w*0.81);
            is_gga = true;
            hfx = w*0.2;
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
            weights.push_back(w*0.75);
            keys.push_back("pbec");
            weights.push_back(w);
            is_gga = true;
            hfx = w*0.25;
        }

        if (key == "m06")
        {
            keys.push_back("m06x");
            weights.push_back(w);
            keys.push_back("m06c");
            weights.push_back(w);
            is_tau_mgga = true;
            hfx = w*0.27;
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


void Functional::set_order(const int order)
{
    int ierr = -1;

    dens_offset = (int)pow(2, order);

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
}
