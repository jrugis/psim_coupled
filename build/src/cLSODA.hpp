#ifndef _CLSODA_H_
#define _CLSODA_H_

#include <fstream>
#include <string>

#include "libsoda/LSODA.h"

#include "global_defs.hpp"

class cCell_flow;

class cLSODA {
  public:
  cLSODA(cCell_flow* flow_, std::ofstream& out_, double abstol_, double reltol_);
  void init(Array1IC& y);
  void run(double t, double tout, Array1IC& y);

  private:
  cCell_flow* flow;
  std::ofstream& out;
  std::vector<double> yin, yout;
  double abstol, reltol; // tolerances for the solver
};

#endif
