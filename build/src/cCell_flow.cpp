/*
 * cCell_flow.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <iostream>

#include "global_defs.hpp"
#include "cCellMesh.hpp"
#include "cCell_flow.hpp"
#include "cCell_calcium.hpp"

// cCell_flow::cCell_flow(int index, double parms[], std::ofstream& out) {
cCell_flow::cCell_flow(cCell_calcium* _parent)
{
  parent = _parent;
  p = &parent->p;
  parent->out << "<Cell_flow> initialising solution vector and constant values" << std::endl;
  init_const();    // first initialise the constant values
  init_solvec();   // initialise solution buffer
}

cCell_flow::~cCell_flow() {
}

void cCell_flow::init_const()
{
    // initial cell volume
	scv(V0) = parent->element_data.col(VOL_e).sum();  

    // apical surface area
	scv(Sa) = 0.0;   
    for (int n = 0; n < parent->mesh->apical_triangles_count; n++) {
      int apical_tri = parent->mesh->apical_triangles(n);
      scv(Sa) += parent->surface_data(apical_tri, AREA_s);
    }

	// basal surface area
    scv(Sb) = 0.0;   
    for (int n = 0; n < parent->mesh->basal_triangles_count; n++) {
      int basal_tri = parent->mesh->basal_triangles(n);
      scv(Sb) += parent->surface_data(basal_tri, AREA_s);
    }

	scv(St) = scv(Sa) / 0.943551157250391;
    scv(aNkcc1) = ( 0.0063812 ) * scv(Sb);

    // Bicarbonate Buffer GB
    double vBB = scv(V0) * ( p->at("kp") * p->at("CO20") - p->at("kn") * p->at("HCO30") * p->at("H0") );
    scv(GB) = 4251.79 / vBB;

    // Tight Junction Na current GtNa
    double vtNa = scv(St) * ( p->at("Vt0") - p->at("VtNa0" )) / F;
    scv(GtNa) = p->at("Qtot0") * p->at("Nal0") / vtNa;
    
    // Tight junction K current GtK
    double vtK = scv(St) * ( p->at("Vt0") - p->at("VtK0") ) / F; 
    scv(GtK) = p->at("Qtot0") * p->at("Kl0") / vtK;
    
    // Apical Ca2+ activated Cl channels GCl
    double PCl = 1.0 / ( 1.0 + pow( p->at("KCaCC") / p->at("c0"), p->at("eta1")) );
    double vCl = PCl * ( p->at("Va0") + p->at("VCl0")) / F;
    scv(GCl) = -p->at("Qtot0") * p->at("Cll0") / vCl;

    // Sodium Proton Antiporter G1
    double JBB0 = scv(GB) * vBB;
    double vNhe1 = scv(Sb) * ( ( p->at("Nae") / ( p->at("Nae") + p->at("KNa") ) ) * ( p->at("H0") / ( p->at("KH") + p->at("H0") ) ) -
		( p->at("Na0") / ( p->at("Na0") + p->at("KNa") ) ) * ( p->at("He") / ( p->at("KH") + p->at("He") ) ) );
    scv(G1) = JBB0 / vNhe1;

    // Sodium Potassium ATPase aNaK
    double JtNa0 = scv(GtNa) * vtNa;
    double JtK0 = scv(GtK) * vtK;
    double JNkcc10 = scv(aNkcc1) * ( ( p->at("a1") - p->at("a2") * p->at("Na0") * p->at("K0") * pow(p->at("Cl0"), 2) )
        / pow( p->at("a3") + p->at("a4") * p->at("Na0") * p->at("K0") * p->at("Cl0"), 2) );
    double vNaK = scv(Sb) * ( p->at("r") * pow(p->at("Ke"), 2) * pow(p->at("Na0"), 3)
        / ( pow(p->at("Ke"), 2) + p->at("alpha1") * pow(p->at("Na0"), 3) ) );                  
	scv(aNaK) = ( ( JtNa0 + JtK0 ) - JNkcc10 ) / ( 3 * vNaK);

    // Ca2+ activated K+ channel GK
    double PKb = 1.0 / ( 1.0 + pow( p->at("KCaKC") / p->at("c0"), p->at("eta2") ));
    double vK = PKb * ( p->at("Vb0") - p->at("VK0") ) / F;
    scv(GK) = ( JNkcc10 + 2 * ( JtNa0 + JtK0 ) ) / ( 3 * vK );

    // Anion exchanger 4 G4
    double JNaK0 = scv(aNaK) * vNaK;
    double JNhe10 = scv(G1) * vNhe1;
    double vAe4 = scv(Sb) *  ( ( p->at("Cle") / ( p->at("Cle") + p->at("KCl") ) ) *
		( p->at("Na0") / ( p->at("Na0") + p->at("KNa") ) ) * pow( p->at("HCO30") / ( p->at("HCO30") + p->at("KB") ), 2) );
    scv(G4) = ( JNkcc10 - 3 * JNaK0 + JNhe10 ) / vAe4;
  
}	

void cCell_flow::init_solvec()
{
  prev_solvec(Nal) = p->at("Nal0");
  prev_solvec(Kl) = p->at("Kl0");
  prev_solvec(Cll) = p->at("Nal0") + p->at("Kl0");
  prev_solvec(VOL) = scv(V0);
  prev_solvec(Na) = p->at("Na0");
  prev_solvec(K) = p->at("K0");
  prev_solvec(Cl) = p->at("Cl0");
  prev_solvec(HCO3) = p->at("HCO30");
  prev_solvec(H) = 1e3 * pow(10, -(p->at("pHi")));
  prev_solvec(Va) = p->at("Va0");
  prev_solvec(Vb) = p->at("Vb0");
  std::cerr << "Cell:" << parent->id << " <Cell_flow> H: " << prev_solvec(H) << std::endl;
}
	  
void cCell_flow::step(){
  //// (3Na+)/(2K+) ATP-ase pump (NaK)
  //double NaKbasalfactor = 0.7;   // Fraction of NaK ATPase in the basal membrane
  //double JNaKb = NaKbasalfactor*Sb*aNaK*( param.r*param.Ke^2 * Na^3 / ( param.Ke^2 + param.alpha1 * Na^3 ) );
  //double JNaKa = (1-NaKbasalfactor)*Sa*aNaK*(param.r * Kl^2 * Na^3 / ( Kl^2 + param.alpha1 * Na^3 ) ); 

  solvec = prev_solvec;

  prev_solvec = solvec;
}

//out << "<Cell_calcium> initial volume: " << volume << std::endl;
