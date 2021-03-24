/*
 * cCell_flow.hpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#ifndef CCELL_FLOW_H_
#define CCELL_FLOW_H_

#include <fstream>
#include <string>
#include <unordered_map>
#include "global_defs.hpp"
#include "cCVode.hpp"
#include "cLSODA.hpp"


class cCell_calcium;

class constant_values {                                                           // invariant cell properties
  public:
  double aNaK, aNkcc1, GtNa, GtK, GCl, GK, G1, G4, GB, St, Sb, Sa, V0, wl;
};

class cCell_flow {
  public:
  cCell_flow(cCell_calcium* parent);
  ~cCell_flow();
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)
  Array1IC solvec, prev_solvec, dsolvec, prev_dsolvec;   // solution vectors for ions
  void step(double t, double dt);
  void secretion(double t, Array1IC& x_ion, Array1IC& dx_ion);
  void save_results();

  private:
  cCell_calcium* parent;
  std::unordered_map<std::string, double> p;
  //int cell_number;
  constant_values s;              // secretion constants vector
  bool solver_initialised;
  cCVode* cvode_solver;
  cLSODA* lsoda_solver;
  std::ofstream ion_file;
  void init_solvec();
  void init_const();
  void init_solver();
  void compute_osmolarities(Array1IC& x_ion, double& Qa, double& Qb, double& Qt);
};

#endif /* CCELL_FLOW_H_ */
