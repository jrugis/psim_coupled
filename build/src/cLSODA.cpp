#include <fstream>
#include <string>

#include "libsoda/LSODA.h"

#include "cLSODA.hpp"
#include "cCell_flow.hpp"
#include "global_defs.hpp"
#include "utils.hpp"

/*
 * f routine. Compute function f(t,y).
 */

static void f(double t, double* y, double* ydot, void* data)
{
  cCell_flow* pt_flow = static_cast<cCell_flow*>(data);

  Array1IC ymat;
  for (int i = 0; i < IONCOUNT; i++) { ymat(i) = y[i]; }
  Array1IC ydotmat;

  pt_flow->secretion(t, ymat, ydotmat);

  for (int i = 0; i < IONCOUNT; i++) { ydot[i] = ydotmat(i); }
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

cLSODA::cLSODA(cCell_flow* flow_, std::ofstream& out_, double abstol_, double reltol_)
  : flow(flow_), out(out_), abstol(abstol_), reltol(reltol_)
{
  out << std::scientific;
  out << "<LSODA>: creating LSODA solver" << std::endl;
  out << " tolerances are " << abstol << " (absolute) and ";
  out << reltol << " (relative)" << std::endl;
  out << std::fixed;
}

void cLSODA::init(Array1IC& y)
{
  yin.resize(IONCOUNT);
  yout.resize(IONCOUNT);
}

void cLSODA::run(double t, double tout, Array1IC& y)
{
  // transfer input values into required format
  for (int i = 0; i < IONCOUNT; i++) { yin[i] = y(i); }

  // call the solver
  int istate = 1;
  LSODA lsoda;
  lsoda.lsoda_update(f, IONCOUNT, yin, yout, &t, tout, &istate, static_cast<void*>(flow), reltol, abstol);
  if (istate < 1) { utils::fatal_error("lsoda failed to compute the solution", out); }

  //  extract the results
  for (int i = 0; i < IONCOUNT; i++) { y(i) = yout[i + 1]; }
}
