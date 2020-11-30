/*
 * cCell_flow.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <iostream>

#include "cCell_flow.hpp"
#include "cCell_calcium.hpp"

// cCell_flow::cCell_flow(int index, double parms[], std::ofstream& out) {
cCell_flow::cCell_flow(cCell_calcium* _parent)
{
  parent = _parent;
  init_solvec();   // initialise solution buffer
  init_constvec();  // initialise the constant values
}

cCell_flow::~cCell_flow() {
}

void cCell_flow::init_constvec()
{
  scv(aNaK) = 0.0;
  scv(GtNa) = 0.0;
  scv(GtK) = 0.0;
  scv(GCl) = 0.0;
  scv(GK) = 0.0;
  scv(G1) = 0.0;
  scv(G4) = 0.0;
  scv(GB) = 0.0;
  scv(St) = 0.0;
  scv(Sb) = 0.0;
  scv(Sa) = 0.0;
}	

void cCell_flow::init_solvec()
{
  parent->out << "<Cell_flow> initialising solution vector" << std::endl;
  prev_solvec(Nal) = parent->p.at("Nal0");
  prev_solvec(Kl) = parent->p.at("Kl0");
  prev_solvec(Cll) = parent->p.at("Nal0") + parent->p.at("Kl0");
  prev_solvec(VOL) = parent->volume;
  prev_solvec(Na) = parent->p.at("Na0");
  prev_solvec(K) = parent->p.at("K0");
  prev_solvec(Cl) = parent->p.at("Cl0");
  prev_solvec(HCO3) = parent->p.at("HCO30");
  prev_solvec(H) = 1e3 * pow(10, -(parent->p.at("pHi")));
  prev_solvec(Va) = parent->p.at("Va0");
  prev_solvec(Vb) = parent->p.at("Vb0");
  parent->out << "<Cell_flow> H: " << prev_solvec(H) << std::endl;
  //std::cout << "<Cell_flow> H: " << prev_solvec(H) << std::endl;
  solvec = prev_solvec;
}
	  
void cCell_flow::step(){
  //// (3Na+)/(2K+) ATP-ase pump (NaK)
  //double NaKbasalfactor = 0.7;   // Fraction of NaK ATPase in the basal membrane
  //double JNaKb = NaKbasalfactor*Sb*aNaK*( param.r*param.Ke^2 * Na^3 / ( param.Ke^2 + param.alpha1 * Na^3 ) );
  //double JNaKa = (1-NaKbasalfactor)*Sa*aNaK*(param.r * Kl^2 * Na^3 / ( Kl^2 + param.alpha1 * Na^3 ) ); 
  
  parent->volume = solvec(VOL); // set the cell volume
}
