/*
 * cLumen.hpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#ifndef CLUMEN_H_
#define CLUMEN_H_

#include <fstream>
#include <string>
#include <vector>

class cAcinus;

#include "global_defs.hpp"

class cLumenTree {
  public:
  cLumenTree(std::ofstream& out);
  ~cLumenTree();

  double get_dnl(const Eigen::Vector3d p);

  private:
  std::string id;
  std::ofstream* out;
  int points_count, segments_count; // the number of points and segments in the lumen tree
  MatrixN3d points;                 // 3x coordinate
  MatrixN2i segments;               // line segment indices, 2x points

  void get_segments();
  void print_info();
};

class cLumen {
  public:
  // cLumen(std::string host_name, int my_rank, int cell_count);
  cLumen(cAcinus* parent, std::string id);
  ~cLumen();
  void step();

  private:
  cAcinus* parent;
  std::string id;
  // tCalcs p[FPCOUNT];  // the fluid flow parameters array
  // cCell_flow *cells[CELLS_COUNT];  // the cells connected to this lumen
  // int my_rank, cell_count;
  // tCalcs na; // sodium
  // tCalcs k;  // patasium
};

#endif /* CLUMEN_ */
