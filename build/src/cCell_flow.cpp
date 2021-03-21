/*
 * cCell_flow.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <iostream>
#include <string>
#include <iomanip>

#include "global_defs.hpp"
#include "cCellMesh.hpp"
#include "cCell_flow.hpp"
#include "cCell_calcium.hpp"
#include "cCVode.hpp"
#include "cLSODA.hpp"
#include "utils.hpp"

#define DEBUGFLOW 1

// cCell_flow::cCell_flow(int index, double parms[], std::ofstream& out) {
cCell_flow::cCell_flow(cCell_calcium* _parent) : solver_initialised(false)
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

#ifdef DEBUGFLOW
  parent->out << "DEBUG: He = " << p["He"] << std::endl;
  parent->out << "DEBUG: Cll0 = " << p["Cll0"] << std::endl;
  parent->out << "DEBUG: Hl = " << p["Hl"] << std::endl;
  parent->out << "DEBUG: HCO3l = " << p["HCO3l"] << std::endl;
  parent->out << "DEBUG: H0 = " << p["H0"] << std::endl;
  parent->out << "DEBUG: CO20 = " << p["CO20"] << std::endl;
  parent->out << "DEBUG: Ul = " << p["Ul"] << std::endl;
  parent->out << "DEBUG: Nae = " << p["Nae"] << std::endl;
  parent->out << "DEBUG: Vt0 = " << p["Vt0"] << std::endl;
  parent->out << "DEBUG: VtNa0 = " << p["VtNa0"] << std::endl;
  parent->out << "DEBUG: VtK0 = " << p["VtK0"] << std::endl;
  parent->out << "DEBUG: VCl0 = " << p["Vcl0"] << std::endl;
  parent->out << "DEBUG: VK0 = " << p["VK0"] << std::endl;
#endif

  parent->out << "<Cell_flow> initialising solution vector and constant values" << std::endl;
  init_const();    // first initialise the constant values
  init_solvec();   // initialise solution buffer
  init_solver();   // initialise the solver
}

cCell_flow::~cCell_flow() {
  if (solver_initialised) {
    if (p.at("odeSolver") == 0) {
      delete cvode_solver;
    }
    else if (p.at("odeSolver") == 1) {
      delete lsoda_solver;
    }
  }
}

void cCell_flow::init_const()
{
    // initial cell volume
    s.V0 = parent->element_data.col(VOL_e).sum();  
    s.wl = s.V0 * 0.02;

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
    double vtNa = s.St * ( p.at("Vt0") - p.at("VtNa0" )) / F_CONST;
    s.GtNa = p.at("Qtot0") * p.at("Nal0") / vtNa;
    
    // Tight junction K current GtK
    double vtK = s.St * ( p.at("Vt0") - p.at("VtK0") ) / F_CONST;
    s.GtK = p.at("Qtot0") * p.at("Kl0") / vtK;
    
    // Apical Ca2+ activated Cl channels GCl
    double PCl = 1.0 / ( 1.0 + pow( p.at("KCaCC") / p.at("c0"), p.at("eta1")) );
    double vCl = PCl * ( p.at("Va0") + p.at("VCl0")) / F_CONST;
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
#ifdef DEBUGFLOW
    parent->out << std::fixed << std::setprecision(16);
    parent->out << "DEBUGFLOW: s.aNaK = " << s.aNaK << std::endl;
    parent->out << "DEBUGFLOW:   JtNa0 = " << JtNa0 << std::endl;
    parent->out << "DEBUGFLOW:   JtK0 = " << JtK0 << std::endl;
    parent->out << "DEBUGFLOW:   JNkcc10 = " << JNkcc10 << std::endl;
    parent->out << "DEBUGFLOW:   vNaK = " << vNaK << std::endl;
    parent->out << "DEBUGFLOW:     s.Sb = " << s.Sb << std::endl;
    parent->out << "DEBUGFLOW:     p.Ke = " << p.at("Ke") << std::endl;
    parent->out << "DEBUGFLOW:     p.Na0 = " << p.at("Na0") << std::endl;
    parent->out << "DEBUGFLOW:     p.alpha1 = " << p.at("alpha1") << std::endl;
    
#endif

    // Ca2+ activated K+ channel GK
    double PKb = 1.0 / ( 1.0 + pow( p.at("KCaKC") / p.at("c0"), p.at("eta2") ));
    double vK = PKb * ( p.at("Vb0") - p.at("VK0") ) / F_CONST;
    s.GK = ( JNkcc10 + 2 * ( JtNa0 + JtK0 ) ) / ( 3 * vK );

    // Anion exchanger 4 G4
    double JNaK0 = s.aNaK * vNaK;
    double JNhe10 = s.G1 * vNhe1;
    double vAe4 = s.Sb *  ( ( p.at("Cle") / ( p.at("Cle") + p.at("KCl") ) ) *
		( p.at("Na0") / ( p.at("Na0") + p.at("KNa") ) ) * pow( p.at("HCO30") / ( p.at("HCO30") + p.at("KB") ), 2) );
    s.G4 = ( JNkcc10 - 3 * JNaK0 + JNhe10 ) / vAe4;

    s.GtNa = 10.0 * s.GtNa;
    s.GtK = 10.0 * s.GtK;
    s.GCl = 0.7 * s.GCl;
    s.aNaK = 1.1 * s.aNaK;
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

void cCell_flow::init_solver()
{
  // which solver to use
  parent->out << "DEBUG: solver flag = " << p.at("odeSolver") << std::endl;

  if (p.at("odeSolver") == 0) {
    cvode_solver = new cCVode(this, parent->out, p.at("odeSolverAbsTol"), p.at("odeSolverRelTol"));
    cvode_solver->init(solvec);
  }
  else if (p.at("odeSolver") == 1) {
    lsoda_solver = new cLSODA(this, parent->out, p.at("odeSolverAbsTol"), p.at("odeSolverRelTol"));
    lsoda_solver->init(solvec);
  }
  else {
    utils::fatal_error("Unrecognised value for odeSolver fluid flow parameter", parent->out);
  }

  solver_initialised = true;
}
  
void cCell_flow::step(double t, double dt){
  // TODO: precompute values here that don't depend on x_ion: PrCl, PrKa, PrKb



 
#ifdef DEBUGFLOW
  // DEBUGGING: run secretion function with initial conditions passed in and store output
  std::string xfilename = parent->id + "_debugsecretion_xion.txt";
  std::ofstream xfile;
  xfile.open(xfilename);
  xfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < IONCOUNT; i++) {
    xfile << solvec(i) << std::endl;
  }
  xfile.close();

  Array1IC dx;
  secretion(t, solvec, dx);

  std::string dxfilename = parent->id + "_debugsecretion_dxion.txt";
  std::ofstream dxfile;
  dxfile.open(dxfilename);
  dxfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < IONCOUNT; i++) {
    dxfile << dx(i) << std::endl;
  }
  dxfile.close();

  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD, 1);
  // END DEBUGGING
#endif

  // invoke the solver here
  if (p.at("odeSolver") == 0) {
    cvode_solver->run(t, t + dt, solvec);
  }
  else if (p.at("odeSolver") == 1) {
    lsoda_solver->run(t, t + dt, solvec);
  }

#ifdef DEBUGFLOW
  std::string solfilename = parent->id + "_debugflowsol.txt";
  std::ofstream solfile;
  solfile.open(solfilename);
  solfile << std::fixed << std::setprecision(16);
  for (int i = 0; i < IONCOUNT; i++) {
    solfile << solvec(i) << std::endl;
  }
  solfile.close();

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  
  // store solution
  prev_solvec = solvec;
}

void cCell_flow::secretion(double t, Array1IC& x_ion, Array1IC& dx_ion){
#ifdef DEBUGFLOW
  std::ofstream debugf;
  debugf.open(parent->id + "_debugsecretion.txt");
  debugf << std::fixed << std::setprecision(16);
#endif

	// %% Currents and fluxes
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// % (3Na+)/(2K+) ATP-ase pump (NaK)
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// NaKbasalfactor = 0.7;        % Fraction of NaK ATPase in the basal membrane
  double NaKbasalfactor = 0.7;   // Fraction of NaK ATPase in the basal membrane

	//JNaKb = NaKbasalfactor*Sb*aNaK*( param.r*param.Ke^2 * Na^3 ... 
	//  / ( param.Ke^2 + param.alpha1 * Na^3 ) );   
  double JNaKb = NaKbasalfactor*s.Sb*s.aNaK*(p.at("r")*pow(p.at("Ke"),2) * pow(x_ion(Na),3) / 
	  ( pow(p.at("Ke"),2) + p.at("alpha1") * pow(x_ion(Na),3) ) );
#ifdef DEBUGFLOW
  debugf << "JNaKb = " << JNaKb << std::endl;
  debugf << "  s.Sb = " << s.Sb << std::endl;
  debugf << "  s.aNaK = " << s.aNaK << std::endl;
  debugf << "  p.r = " << p.at("r") << std::endl;
  debugf << "  p.Ke = " << p.at("Ke") << std::endl;
  debugf << "  p.alpha1 = " << p.at("alpha1") << std::endl;
  debugf << "  Na = " << x_ion(Na) << std::endl;
#endif

	//JNaKa = (1-NaKbasalfactor)*Sa*aNaK*(param.r * Kl^2 * Na^3 ...
	//  / ( Kl^2 + param.alpha1 * Na^3 ) ); 
  double JNaKa = (1-NaKbasalfactor)*s.Sa*s.aNaK*(p.at("r") * pow(x_ion(Kl),2) * pow(x_ion(Na),3) /
	  ( pow(x_ion(Kl),2) + p.at("alpha1") * pow(x_ion(Na),3) ) ); 
#ifdef DEBUGFLOW
  debugf << "JNaKa = " << JNaKa << std::endl;
  debugf << "  s.Sa = " << s.Sa << std::endl;
  debugf << "  s.aNaK = " << s.aNaK << std::endl;
  debugf << "  p.r = " << p.at("r") << std::endl;
  debugf << "  p.alpha1 = " << p.at("alpha1") << std::endl;
  debugf << "  Kl = " << x_ion(Kl) << std::endl;
  debugf << "  Na = " << x_ion(Na) << std::endl;
#endif

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Nernst potentials (J/C)1e3 = mV
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//VCl = param.RTF * log( Cll / Cl );        
  double VCl = RTF * log( x_ion(Cll) / x_ion(Cl) );
#ifdef DEBUGFLOW
  debugf << "VCl = " << VCl << std::endl;
#endif

	//VKb = param.RTF * log( param.Ke / K ); 
  double VKb = RTF * log( p.at("Ke") / x_ion(K));    
#ifdef DEBUGFLOW
  debugf << "VKb = " << VKb << std::endl;
#endif

	//VKa = param.RTF * log( Kl / K );
  double VKa = RTF * log( x_ion(Kl) / x_ion(K));         
#ifdef DEBUGFLOW
  debugf << "VKa = " << VKa << std::endl;
#endif

	//VtNa = param.RTF * log( Nal / param.Nae );
  double VtNa = RTF * log( x_ion(Nal) / p.at("Nae"));         
#ifdef DEBUGFLOW
  debugf << "VtNa = " << VtNa << std::endl;
#endif

	//VtK = param.RTF * log( Kl / param.Ke );   
  double VtK = RTF * log( x_ion(Kl) / p.at("Ke"));         
#ifdef DEBUGFLOW
  debugf << "VtK = " << VtK << std::endl;
#endif

	//Vt = Va - Vb;
  double Vt = x_ion(Va) - x_ion(Vb);
#ifdef DEBUGFLOW
  debugf << "Vt = " << Vt << std::endl;
#endif

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
    PKa += (1.0 / (1.0 + pow(p.at("KCaKC") / ca_tri, p.at("eta2")))) * area_tri;
  }
  double PrCl = PCl / s.Sa; 
  double PrKa = PKa / s.Sa; 
#ifdef DEBUGFLOW
  debugf << "s.Sa = " << s.Sa << std::endl;
  debugf << "PrCl = " << PrCl << std::endl;
  debugf << "p.KCaKC = " << p.at("KCaKC") << std::endl;
  debugf << "PrKa = " << PrKa << std::endl;
#endif

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
    PKb += (1.0 / (1.0 + pow(p.at("KCaKC") / ca_tri, p.at("eta2")))) * area_tri;
  }
  double PrKb = PKb / s.Sb; 
#ifdef DEBUGFLOW
  debugf << "PrKb = " << PrKb << std::endl;
#endif

	//JCl = GCl * PrCl * ( Va + VCl ) / param.F;  % fS.micro-metres^2.mV.mol.C^-1
	//JKa = GK * PrKa * ( Va - VKa ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
	//JKb = GK * PrKb * ( Vb - VKb ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
  double JCl = s.GCl * PrCl * ( x_ion(Va) + VCl) / F_CONST;
  double JKa = s.GK * PrKa * ( x_ion(Va) - VKa) / F_CONST;
  double JKb = s.GK * PrKb * ( x_ion(Vb) - VKb) / F_CONST;

#ifdef DEBUGFLOW
  debugf << "JCl = " << JCl << std::endl;
  debugf << "JKa = " << JKa << std::endl;
  debugf << "JKb = " << JKb << std::endl;
#endif


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Tight Junction Na+ and K+ currents
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JtNa = GtNa * St * ( Vt - VtNa ) / param.F;   % fS.micro-metres^2.mV.mol.C^-1
	//JtK = GtK * St * ( Vt - VtK ) / param.F;      % fS.micro-metres^2.mV.mol.C^-1
  double JtNa = s.GtNa * s.St * ( Vt - VtNa ) / F_CONST;
  double JtK = s.GtK * s.St * ( Vt - VtK ) / F_CONST;
#ifdef DEBUGFLOW
  debugf << "JtNa = " << JtNa << std::endl;
  debugf << "JtK = " << JtK << std::endl;
#endif

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Osmolarities 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//Qa = param.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - param.CO20 + param.Ul );     % micro-metres^3.s^-1
	//Qb = param.B2 * ( 2 * ( Na + K + H ) + param.CO20 - ...
	//                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );
	//Qt = param.B3 * ( 2 * ( Nal + Kl ) + param.Ul - ....
	//                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); % micro-metres^3.s^-1
	//Qtot=(Qa+Qt);                                     % micro-metres^3.s^-1
  double Qa = p.at("B1") * ( 2 * ( x_ion(Nal) + x_ion(Kl) - x_ion(Na) - x_ion(K) - x_ion(H) ) - p.at("CO20") + p.at("Ul") );
  double Qb = p.at("B2") * ( 2 * ( x_ion(Na) + x_ion(K) + x_ion(H) ) + p.at("CO20") - ( p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e") ) );
  double Qt = p.at("B3") * ( 2 * ( x_ion(Nal) + x_ion(Kl) ) + p.at("Ul") - ( p.at("Nae") + p.at("Ke") + p.at("Cle") + p.at("HCO3e") ) );
  double Qtot = Qa + Qt;

#ifdef DEBUGFLOW
  debugf << "Qa = " << Qa << std::endl;
  debugf << "Qb = " << Qb << std::endl;
  debugf << "Qt = " << Qt << std::endl;
  debugf << "Qtot = " << Qtot << std::endl;
#endif

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Na+ K+ 2Cl- co-transporter (Nkcc1)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JNkcc1 = aNkcc1 * Sb * ( param.a1 - param.a2 * Na * K * Cl^2 ) ...
	//                                             / ( param.a3 + param.a4 * Na * K * Cl^2 );
  double JNkcc1 = s.aNkcc1 * s.Sb * ( p.at("a1") - p.at("a2") * x_ion(Na) * x_ion(K) * pow(x_ion(Cl),2) ) / 
	  ( p.at("a3") + p.at("a4") * x_ion(Na) * x_ion(K) * pow(x_ion(Cl),2) );
#ifdef DEBUGFLOW
  debugf << "JNkcc1 = " << JNkcc1 << std::endl;
  debugf << "  s.aNkcc1 = " << s.aNkcc1 << std::endl;
#endif
       
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JAe4 = Sb * G4 * ( ( param.Cle / ( param.Cle + param.KCl ) ) * ( Na / ( Na + param.KNa ) ) ...
	//             * ( HCO3 / ( HCO3 + param.KB ) )^2 );       
  double JAe4 = s.Sb * s.G4 * ( ( p.at("Cle") / ( p.at("Cle") + p.at("KCl") ) ) * ( x_ion(Na) / 
	  ( x_ion(Na) + p.at("KNa") ) ) * pow( x_ion(HCO3) / ( x_ion(HCO3) + p.at("KB") ),2) );
#ifdef DEBUGFLOW
  debugf << "JAe4 = " << JAe4 << std::endl;
#endif

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Na+ / H+ Anion exchanger (Nhe1)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//JNhe1 = Sb * G1 * ( ( param.Nae / ( param.Nae + param.KNa ) ) * ( H / ( param.KH + H ) )...
	//                          - ( Na / ( Na + param.KNa ) ) * ( param.He / ( param.KH + param.He ) ) ); 
  double JNhe1 = s.Sb * s.G1 * ( ( p.at("Nae") / ( p.at("Nae") + p.at("KNa") ) ) * ( x_ion(H) / 
	  ( p.at("KH") + x_ion(H) ) ) - ( x_ion(Na) / ( x_ion(Na) + p.at("KNa") ) ) * ( p.at("He") / ( p.at("KH") + p.at("He") ) ) );
#ifdef DEBUGFLOW
  debugf << "JNhe1 = " << JNhe1 << std::endl;
#endif

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Bicarbonate Buffer (Reaction Term)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% This equation is a reaction inside the cell, note how it depends on the
	//% cellular volume
	//JBB = w * GB * ( param.kp * param.CO20 - param.kn * HCO3 * H ); 
  double JBB = x_ion(VOL) * s.GB * ( p.at("kp") * p.at("CO20") - p.at("kn") * x_ion(HCO3) * x_ion(H) );                 
#ifdef DEBUGFLOW
  debugf << "JBB = " << JBB << std::endl;
#endif

  // Equations
//  dx(1) = ( JtNa - Qtot*Nal + 3*JNaKa )/param.wl(i);
//  dx(2) = ( JtK  - Qtot*Kl + JKa - 2*JNaKa)/param.wl(i);
//  dx(3) = ( -JCl - Qtot*Cll     )/param.wl(i);
//  dx(4) = Qb - Qa;
//  dx(5) = ( JNkcc1 - 3*(JNaKb+JNaKa) + JNhe1 - JAe4 - dx(4) * Na ) / w;
//  dx(6) = ( JNkcc1 + 2*(JNaKb+JNaKa) - JKb - JKa - dx(4) * K ) / w;
//  dx(7) = ( 2 * JNkcc1 + JAe4 + JCl - dx(4) * Cl ) / w;
//  dx(8) = ( JBB - 2 * JAe4 - dx(4) * HCO3 ) / w;
//  dx(9) = ( JBB - JNhe1 - dx(4) * H ) / w;
//  dx(10) = 100*(-JCl - JNaKa - JKa - JtK - JtNa);      % Note the arbitrary factor of 100, just to make sure Va is fast.
//  dx(11) = 100*(     - JNaKb - JKb + JtK + JtNa);
//  enum solution_values { Nal, Kl, Cll, VOL, Na, K, Cl, HCO3, H, Va, Vb, IONCOUNT }; // solution vector components
  dx_ion(Nal) = (JtNa - Qtot * x_ion(Nal) + 3.0 * JNaKa) / s.wl;
  dx_ion(Kl) = (JtK - Qtot * x_ion(Kl) + JKa - 2.0 * JNaKa) / s.wl;
  dx_ion(Cll) = (-JCl - Qtot * x_ion(Cll)) / s.wl;
  dx_ion(VOL) = Qb - Qa;
  dx_ion(Na) = (JNkcc1 - 3.0 * (JNaKb + JNaKa) + JNhe1 - JAe4 - dx_ion(VOL) * x_ion(Na)) / x_ion(VOL);
  dx_ion(K) = (JNkcc1 + 2.0 * (JNaKb + JNaKa) - JKb - JKa - dx_ion(VOL) * x_ion(K)) / x_ion(VOL);
  dx_ion(Cl) = (2.0 * JNkcc1 + JAe4 + JCl - dx_ion(VOL) * x_ion(Cl)) / x_ion(VOL);
  dx_ion(HCO3) = (JBB - 2.0 * JAe4 - dx_ion(VOL) * x_ion(HCO3)) / x_ion(VOL);
  dx_ion(H) = (JBB - JNhe1 - dx_ion(VOL) * x_ion(H)) / x_ion(VOL);
  dx_ion(Va) = 100.0 * (-JCl - JNaKa - JKa - JtK - JtNa);  // Note the arbitrary factor of 100, just to make sure Va is fast.
  dx_ion(Vb) = 100.0 * (     - JNaKb - JKb + JtK + JtNa);

#ifdef DEBUGFLOW
  debugf << std::endl;
  debugf << "dx(1) = " << dx_ion(Nal) << std::endl;
  debugf << "dx(2) = " << dx_ion(Kl) << std::endl;
  debugf << "dx(3) = " << dx_ion(Cll) << std::endl;
  debugf << "dx(4) = " << dx_ion(VOL) << std::endl;
  debugf << "dx(5) = " << dx_ion(Na) << std::endl;
  debugf << "dx(6) = " << dx_ion(K) << std::endl;
  debugf << "dx(7) = " << dx_ion(Cl) << std::endl;
  debugf << "dx(8) = " << dx_ion(HCO3) << std::endl;
  debugf << "dx(9) = " << dx_ion(H) << std::endl;
  debugf << "dx(10) = " << dx_ion(Va) << std::endl;
  debugf << "dx(11) = " << dx_ion(Vb) << std::endl;

  debugf.close();
#endif
}

//out << "<Cell_calcium> initial volume: " << volume << std::endl;
