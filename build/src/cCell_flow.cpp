/*
 * cCell_flow.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include "cCell_flow.hpp"
#include "cCell_calcium.hpp"

// cCell_flow::cCell_flow(int index, double parms[], std::ofstream& out) {
cCell_flow::cCell_flow(cCell_calcium* _parent)
{
  parent = _parent;
  parent->out << "<Cell_flow> instantiated" << std::endl;
  p = &parent->p;
}

cCell_flow::~cCell_flow() {}
