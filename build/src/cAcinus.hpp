/*
 * cAcinus.hpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#ifndef CACINUS_H_
#define CACINUS_H_

#include <fstream>
#include <string>

class cLumen;

#include "global_defs.hpp"

class cAcinus {
friend class cLumen;
public:
  //cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count, int lumen_rank);
  cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count);
  ~cAcinus();
  void run();

private:
  std::string id;
  std::ofstream out;
  double p[PCOUNT]; // the calcium model parameters array
  //int my_rank, cell_rank, cell_count, lumen_rank;
  int my_rank, cell_rank, cell_count;
  cLumen* lumen;

  double snd_recv(double t, double dt);
};

#endif /* CACINUS_H_ */

