#ifndef Functional_h_
#define Functional_h_

#include <xcfun.h>
#include <string>
#include <vector>

class Functional
{
    public:

        Functional();
        ~Functional();

        void set_functional(const int verbosity, const char *line, double &hfx, double &mu, double &beta);
        void set_order(const int order);

        bool is_gga;      // FIXME make private
        bool is_tau_mgga; // FIXME make private
        int  dens_offset; // FIXME make private
        xc_functional fun;

    private:

        Functional(const Functional &rhs);            // not implemented
        Functional &operator=(const Functional &rhs); // not implemented

        char *functional_line;

        void parse(const char *line,
                   double &hfx,
                   double &mu,
                   double &beta);
        void nullify();

        std::vector<std::string> keys;
        std::vector<double>      weights;

        bool is_synced;
};

#endif // Functional_h_
