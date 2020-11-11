/*
 * cLumen.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "cAcinus.hpp"
#include "cLumen.hpp"

//cLumen::cLumen(std::string host_name, int rank, int c_count, int a_rank) {
cLumen::cLumen(cAcinus* _parent,   std::string _id) {
  parent = _parent;
  id = _id;

  parent->out << "<Lumen> instantiated" << std::endl;
}

cLumen::~cLumen() {
}

void cLumen::step() {
}


