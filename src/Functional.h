#pragma once

#include <string>
#include <vector>

#include "xcfun.h"

class Functional
{
  public:
    Functional();
    ~Functional();

    void set_functional(const char *line);
    int set_order(const int order, xc_functional fun) const;

    bool is_gga;                   // FIXME make private
    bool is_tau_mgga;              // FIXME make private
    std::vector<std::string> keys; // FIXME make private
    std::vector<double> weights;   // FIXME make private

  private:
    Functional(const Functional &rhs);            // not implemented
    Functional &operator=(const Functional &rhs); // not implemented

    char *functional_line;

    void parse(const char *line);
    void nullify();

    bool is_synced;
};
