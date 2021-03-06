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
  p = parent->p;  // make a "local" copy of the parameters

  // add dependant parameters to the map
  p["He"] = 1e3 * pow(10, -p.at("pHe"));
  p["Cll0"] = p.at("Nal0") + p.at("Kl0");  
  p["Hl"] = 1e3 * pow(10, -p.at("pHl"));
  p["HCO3l"] = p.at("Kl0") + p.at("Nal0") - p.at("Cll0") + p.at("Hl");
  p["H0"] = 1e3 * pow(10, -p.at("pHi"));
  p["CO20"] = (0.197e4 * (p.at("CO2l") + p.at("CO2e")) - p.at("kn") * p.at("HCO30") * p.at("H0")) / (2 * 0.197e4 - p.at("kp"));
  p["Ul"] = (p.at("B2") / p.at("B1")) * (2 * (p.at("Na0") + p.at("K0") + p.at("H0")) + p.at("CO20") - 
	(p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e"))) -
    (2 * (p.at("Nal0") + p.at("Kl0") - p.at("Na0") - p.at("K0") - p.at("H0")) - p.at("CO20"));
  p["Vt0"] = p.at("Va0") - p.at("Vb0");
  p["VtNa0"] = RTF * log(p.at("Nal0") / p.at("Nae"));
  p["VtK0"] = RTF * log(p.at("Kl0") / p.at("Ke"));
  p["VCl0"] = RTF * log(p.at("Cll0") / p.at("Cl0"));
  p["VK0"] = RTF * log(p.at("Ke") / p.at("K0"));
  
  parent->out << "<Cell_flow> initialising solution vector and constant values" << std::endl;
  init_const();    // first initialise the constant values
  init_solvec();   // initialise solution buffer
}

cCell_flow::~cCell_flow() {
}

void cCell_flow::init_const()
{
    // initial cell volume
	s.V0 = parent->element_data.col(VOL_e).sum();  

    // apical surface area
	s.Sa = 0.0;   
    for (int n = 0; n < parent->mesh->apical_triangles_count; n++) {
      int apical_tri = parent->mesh->apical_triangles(n);
      s.Sa += parent->surface_data(apical_tri, AREA_s);
    }

	// basal surface area
    s.Sb = 0.0;   
    for (int n = 0; n < parent->mesh->basal_triangles_count; n++) {
      int basal_tri = parent->mesh->basal_triangles(n);
      s.Sb += parent->surface_data(basal_tri, AREA_s);
    }

	s.St = s.Sa / 0.943551157250391;
    s.aNkcc1 = ( 0.0063812 ) * s.Sb;

    // Bicarbonate Buffer GB
    double vBB = s.V0 * ( p.at("kp") * p.at("CO20") - p.at("kn") * p.at("HCO30") * p.at("H0") );
    s.GB = 4251.79 / vBB;

    // Tight Junction Na current GtNa
    double vtNa = s.St * ( p.at("Vt0") - p.at("VtNa0" )) / F;
    s.GtNa = p.at("Qtot0") * p.at("Nal0") / vtNa;
    
    // Tight junction K current GtK
    double vtK = s.St * ( p.at("Vt0") - p.at("VtK0") ) / F; 
    s.GtK = p.at("Qtot0") * p.at("Kl0") / vtK;
    
    // Apical Ca2+ activated Cl channels GCl
    double PCl = 1.0 / ( 1.0 + pow( p.at("KCaCC") / p.at("c0"), p.at("eta1")) );
    double vCl = PCl * ( p.at("Va0") + p.at("VCl0")) / F;
    s.GCl = -p.at("Qtot0") * p.at("Cll0") / vCl;

    // Sodium Proton Antiporter G1
    double JBB0 = s.GB * vBB;
    double vNhe1 = s.Sb * ( ( p.at("Nae") / ( p.at("Nae") + p.at("KNa") ) ) * ( p.at("H0") / ( p.at("KH") + p.at("H0") ) ) -
		( p.at("Na0") / ( p.at("Na0") + p.at("KNa") ) ) * ( p.at("He") / ( p.at("KH") + p.at("He") ) ) );
    s.G1 = JBB0 / vNhe1;

    // Sodium Potassium ATPase aNaK
    double JtNa0 = s.GtNa * vtNa;
    double JtK0 = s.GtK * vtK;
    double JNkcc10 = s.aNkcc1 * ( ( p.at("a1") - p.at("a2") * p.at("Na0") * p.at("K0") * pow(p.at("Cl0"), 2) )
        / pow( p.at("a3") + p.at("a4") * p.at("Na0") * p.at("K0") * p.at("Cl0"), 2) );
    double vNaK = s.Sb * ( p.at("r") * pow(p.at("Ke"), 2) * pow(p.at("Na0"), 3)
        / ( pow(p.at("Ke"), 2) + p.at("alpha1") * pow(p.at("Na0"), 3) ) );                  
	s.aNaK = ( ( JtNa0 + JtK0 ) - JNkcc10 ) / ( 3 * vNaK);

    // Ca2+ activated K+ channel GK
    double PKb = 1.0 / ( 1.0 + pow( p.at("KCaKC") / p.at("c0"), p.at("eta2") ));
    double vK = PKb * ( p.at("Vb0") - p.at("VK0") ) / F;
    s.GK = ( JNkcc10 + 2 * ( JtNa0 + JtK0 ) ) / ( 3 * vK );

    // Anion exchanger 4 G4
    double JNaK0 = s.aNaK * vNaK;
    double JNhe10 = s.G1 * vNhe1;
    double vAe4 = s.Sb *  ( ( p.at("Cle") / ( p.at("Cle") + p.at("KCl") ) ) *
		( p.at("Na0") / ( p.at("Na0") + p.at("KNa") ) ) * pow( p.at("HCO30") / ( p.at("HCO30") + p.at("KB") ), 2) );
    s.G4 = ( JNkcc10 - 3 * JNaK0 + JNhe10 ) / vAe4;
}	

void cCell_flow::init_solvec()
{
	// Nal = x(1);
	// Kl = x(2);
	// Cll = x(3);
	// w = x(4);
	// Na = x(5);
	// K = x(6);
	// Cl = x(7);
	// HCO3 = x(8);
	// H = x(9);
	// Va = x(10);
	// Vb = x(11);
  solvec(Nal) = p.at("Nal0");
  solvec(Kl) = p.at("Kl0");
  solvec(Cll) = p.at("Nal0") + p.at("Kl0");
  solvec(VOL) = s.V0;
  solvec(Na) = p.at("Na0");
  solvec(K) = p.at("K0");
  solvec(Cl) = p.at("Cl0");
  solvec(HCO3) = p.at("HCO30");
  solvec(H) = 1e3 * pow(10, -(p.at("pHi")));
  solvec(Va) = p.at("Va0");
  solvec(Vb) = p.at("Vb0");
  //std::cerr << "Cell:" << parent->id << " <Cell_flow> H: " << prev_solvec(H) << std::endl;
}
  
void cCell_flow::step(){
	// %% Currents and fluxes
  	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// % (3Na+)/(2K+) ATP-ase pump (NaK)
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// NaKbasalfactor = 0.7;        % Fraction of NaK ATPase in the basal membrane
  double NaKbasalfactor = 0.7;   // Fraction of NaK ATPase in the basal membrane

	//JNaKb = NaKbasalfactor*Sb*aNaK*( param.r*param.Ke^2 * Na^3 ... 
	//  / ( param.Ke^2 + param.alpha1 * Na^3 ) );   
  double JNaKb = NaKbasalfactor*s.Sb*s.aNaK*(p.at("r")*pow(p.at("Ke"),2) * pow(solvec(Na),3) / 
	  ( pow(p.at("Ke"),2) + p.at("alpha1") * pow(solvec(Na),3) ) );

	//JNaKa = (1-NaKbasalfactor)*Sa*aNaK*(param.r * Kl^2 * Na^3 ...
	//  / ( Kl^2 + param.alpha1 * Na^3 ) ); 
  double JNaKa = (1-NaKbasalfactor)*s.Sa*s.aNaK*(p.at("r") * pow(solvec(Kl),2) * pow(solvec(Na),3) /
	  ( pow(solvec(Kl),2) + p.at("alpha1") * pow(solvec(Na),3) ) ); 

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Nernst potentials (J/C)1e3 = mV
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//VCl = param.RTF * log( Cll / Cl );        
  double VCl = RTF * log( solvec(Cll) / solvec(Cl) );

	//VKb = param.RTF * log( param.Ke / K ); 
  double VKb = RTF * log( p.at("Ke") / solvec(K));    

	//VKa = param.RTF * log( Kl / K );
  double VKa = RTF * log( solvec(Kl) / solvec(K));         

	//VtNa = param.RTF * log( Nal / param.Nae );
  double VtNa = RTF * log( solvec(Nal) / p.at("Nae"));         

	//VtK = param.RTF * log( Kl / param.Ke );   
  double VtK = RTF * log( solvec(Kl) / p.at("Ke"));         

	//Vt = Va - Vb;
  double Vt = solvec(Va) - solvec(Vb);


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Ca2+ activated channels.
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Cl channels
	//PCl = 1 ./ ( 1 + ( param.KCaCC ./ (cav_tri.*(w_IPR>0)) ).^param.eta1 ); 
	//PrCl = sum( PCl .* surf_tri )/Sa;

	//% K Channels
	//PKb = 1 ./ ( 1 + ( param.KCaKC ./ (cav_tri.*(w_basal>0)) ).^param.eta2 );
	//PKa = 1 ./ ( 1 + ( param.KCaKC ./ (cav_tri.*(w_IPR>0)) ).^param.eta2 ); 
	//PrKb = sum( PKb .* surf_tri )/Sb;
	//PrKa = sum( PKa .* surf_tri )/Sa;
  double PCl = 0.0;
  double PKa = 0.0;
  for (int n = 0; n < parent->mesh->apical_triangles_count; n++) { // apical tris
    int this_tri = parent->mesh->apical_triangles(n);
	double area_tri = parent->surface_data(this_tri, AREA_s);

	double ca_tri = 0.0;  // average calcium at this triangle
    for (int i=0; i<3; i++){
      int vertex_index = parent->mesh->mesh_vals.surface_triangles(this_tri, i);
      ca_tri += parent->solvec(vertex_index);  // note: Ca is first
    }
	ca_tri /= 3.0;  

    PCl += (1.0 / (1.0 + pow(p.at("KCaCC") / ca_tri, p.at("eta1")))) * area_tri;
    PKa += (1.0 / (1.0 + pow(p.at("KCaKC") / ca_tri, p.at("eta1")))) * area_tri;
  }
  double PrCl = PCl / s.Sa; 
  double PrKa = PKa / s.Sa; 

  double PKb = 0.0;
  for (int n = 0; n < parent->mesh->basal_triangles_count; n++) { // basal tris
    int this_tri = parent->mesh->basal_triangles(n);
	double area_tri = parent->surface_data(this_tri, AREA_s);

	double ca_tri = 0.0;  // average calcium at this triangle
    for (int i=0; i<3; i++){
      int vertex_index = parent->mesh->mesh_vals.surface_triangles(this_tri, i);
      ca_tri += parent->solvec(vertex_index);  // note: Ca is first
    }
	ca_tri /= 3.0;  
    PKb += (1.0 / (1.0 + pow(p.at("KCaKC") / ca_tri, p.at("eta1")))) * area_tri;
  }
  double PrKb = PKb / s.Sb; 

	//JCl = GCl * PrCl * ( Va + VCl ) / param.F;  % fS.micro-metres^2.mV.mol.C^-1
	//JKa = GK * PrKa * ( Va - VKa ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
	//JKb = GK * PrKb * ( Vb - VKb ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
  double JCl = s.GCl * PrCl * ( solvec(Va) + VCl) / F;
  double JKa = s.GK * PrKa * ( solvec(Va) - VKa) / F;
  double JKb = s.GK * PrKb * ( solvec(Vb) - VKb) / F;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Tight Junction Na+ and K+ currents
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JtNa = GtNa * St * ( Vt - VtNa ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
	//JtK = GtK * St * ( Vt - VtK ) / param.F;      % fS.micro-metres^2.mV.mol.C^-1
  double JtNa = s.GtNa * s.St * ( Vt - VtNa ) / F;
  double JtK = s.GtK * s.St * ( Vt - VtK ) / F;
  
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Osmolarities 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//Qa = param.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - param.CO20 + param.Ul );     % micro-metres^3.s^-1
	//Qb = param.B2 * ( 2 * ( Na + K + H ) + param.CO20 - ...
	//                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );
	//Qt = param.B3 * ( 2 * ( Nal + Kl ) + param.Ul - ....
	//                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); % micro-metres^3.s^-1
	//Qtot=(Qa+Qt);                                     % micro-metres^3.s^-1
  double Qa = p.at("B1") * ( 2 * ( solvec(Nal) + solvec(Kl) - solvec(Na) - solvec(K) - solvec(H) ) - p.at("CO20") + p.at("Ul") );
  double Qb = p.at("B2") * ( 2 * ( solvec(Na) + solvec(K) + solvec(H) ) + p.at("CO20") - ( p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e") ) );
  double t = p.at("B3") * ( 2 * ( solvec(Nal) + solvec(Kl) ) + p.at("Ul") - ( p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e") ) );
  double Qt = p.at("B3") * ( 2 * ( solvec(Nal) + solvec(Kl) ) + p.at("Ul") - ( p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e") ) );
  double Qtot = Qa + Qt;
  
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Na+ K+ 2Cl- co-transporter (Nkcc1)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JNkcc1 = aNkcc1 * Sb * ( param.a1 - param.a2 * Na * K * Cl^2 ) ...
	//                                             / ( param.a3 + param.a4 * Na * K * Cl^2 );
  double JNkcc1 = s.aNkcc1 * s.Sb * ( p.at("a1") - p.at("a2") * solvec(Na) * solvec(K) * pow(solvec(Cl),2) ) / 
	  ( p.at("a3") + p.at("a4") * solvec(Na) * solvec(K) * pow(solvec(Cl),2) );
       
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JAe4 = Sb * G4 * ( ( param.Cle / ( param.Cle + param.KCl ) ) * ( Na / ( Na + param.KNa ) ) ...
	//             * ( HCO3 / ( HCO3 + param.KB ) )^2 );       
  double JAe4 = s.Sb * s.G4 * ( ( p.at("Cle") / ( p.at("Cle") + p.at("KCl") ) ) * ( solvec(Na) / 
	  ( solvec(Na) + p.at("KNa") ) ) * pow( solvec(HCO3) / ( solvec(HCO3) + p.at("KB") ),2) );

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Na+ / H+ Anion exchanger (Nhe1)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JNhe1 = Sb * G1 * ( ( param.Nae / ( param.Nae + param.KNa ) ) * ( H / ( param.KH + H ) )...
	//                          - ( Na / ( Na + param.KNa ) ) * ( param.He / ( param.KH + param.He ) ) ); 
  double JNhe1 = s.Sb * s.G1 * ( ( p.at("Nae") / ( p.at("Nae") + p.at("KNa") ) ) * ( solvec(H) / 
	  ( p.at("KH") + solvec(H) ) ) - ( solvec(Na) / ( solvec(Na) + p.at("KNa") ) ) * ( p.at("He") / ( p.at("KH") + p.at("He") ) ) );
  
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Bicarbonate Buffer (Reaction Term)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% This equation is a reaction inside the cell, note how it depends on the
	//% cellular volume
	//JBB = w * GB * ( param.kp * param.CO20 - param.kn * HCO3 * H ); 
  double JBB = solvec(VOL) * s.GB * ( p.at("kp") * p.at("CO20") - p.at("kn") * solvec(HCO3) * solvec(H) );                 
 
  // invoke the solver here...
  
  prev_solvec = solvec;
}

//out << "<Cell_calcium> initial volume: " << volume << std::endl;
