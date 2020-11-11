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

class cCell_calcium;

#include "global_defs.hpp"

class cCell_flow {
public:
  //cCell_flow(int my_index, tCalcs p[], std::ofstream& out);
  cCell_flow(cCell_calcium* parent);
  ~cCell_flow();
  void run();

private:
  cCell_calcium* parent;
  tCalcs *p; // the fluid flow parameters array
  int my_index;
  tCalcs omega; // volume
  tCalcs na;    // sodium
  tCalcs nk;    // potassium
  tCalcs cl;    // chloride
  tCalcs hc03;  // bicarbonate
  tCalcs h;     // hydrogen
  
  //tCalcs snd_recv(tCalcs t, tCalcs dt);
};

#endif /* CCELL_FLOW_H_ */

