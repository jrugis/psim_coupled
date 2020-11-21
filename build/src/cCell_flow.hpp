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

class cCell_calcium;

#include "global_defs.hpp"

class cCell_flow {
  public:
  cCell_flow(cCell_calcium* parent);
  ~cCell_flow();
  void run();

  private:
  cCell_calcium* parent;
  std::unordered_map<std::string, double>* p;
  int my_index;

  // double omega; // volume
  // double na;    // sodium
  // double nk;    // potassium
  // double cl;    // chloride
  // double hc03;  // bicarbonate
  // double h;     // hydrogen

  // double snd_recv(double t, double dt);
};

#endif /* CCELL_FLOW_H_ */
