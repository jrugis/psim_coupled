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
  void step();
  void secretion(double t, Array1IC& x_ion, Array1IC& dx_ion);

  private:
  cCell_calcium* parent;
  std::unordered_map<std::string, double> p;
  //int cell_number;
  Array1IC solvec, prev_solvec;   // solution vectors for ions
  constant_values s;              // secretion constants vector
  void init_solvec();
  void init_const();
};

#endif /* CCELL_FLOW_H_ */
