/*
 * cCell_calcium.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <vector>

#include "cCellMesh.hpp"
#include "cCell_calcium.hpp"
#include "cCell_flow.hpp"
#include "global_defs.hpp"
#include "utils.hpp"

cCell_calcium::cCell_calcium(std::string host_name, int my_rank, int a_rank)
{
  cell_number = my_rank;
  acinus_rank = a_rank;
  acinus_id = "a" + std::to_string(acinus_rank + 1);
  id = acinus_id + "c" + std::to_string(my_rank);
  out.open(id + ".out");
  out << "<Cell_calcium> host_name: " << host_name << std::endl;

  utils::get_parameters(acinus_id, cell_number, p, out);
  // add dependant parameters to the parameters map
  p["ce0"] = (p.at("ct") - p.at("c0")) / p.at("Gamma");
  p["h0"] = pow(p.at("K_h"), 4) / (pow(p.at("K_h"), 4) + pow(p.at("c0"), 4));

  // add dependant default parameters to the parameters map if not already defined
  if (not p.count("ip0"))
    p["ip"] = 0.5 * (2e-4 / 0.1) * ((pow(p.at("K3K"), 2) + pow(p.at("c0"), 2)) / pow(p.at("c0"), 2));

  mesh = new cCellMesh(id, this);
  mesh->print_info();

  // common cells for inter-cellular exchange
  out << "<Cell_calcium> common faces with cells:";
  int other_cell = -1;
  int face_count = 0;
  for (int r = 0; r < mesh->common_triangles.rows(); r++) {
    if (mesh->common_triangles(r, oCell) != other_cell) {
      if (other_cell != -1) {
        cells.push_back({other_cell, face_count});
        face_count = 0;
      }
      other_cell = mesh->common_triangles(r, oCell);
      out << " " << other_cell + 1; // cells are zero indexed
    }
    face_count++;
  }
  cells.push_back({other_cell, face_count}); // one more time
  out << std::endl;

  make_matrices(); // create the constant matrices
  init_solvec();   // initialise solution buffer
  ca_file.open(id + "_ca.bin", std::ios::binary);
  ip3_file.open(id + "_ip3.bin", std::ios::binary);
  cer_file.open(id + "_cer.bin", std::ios::binary);

  flow = new cCell_flow(this);
}

cCell_calcium::~cCell_calcium()
{
  ca_file.close();
  ip3_file.close();
  cer_file.close();
  out.close();
  delete mesh;
  delete flow;
}

void cCell_calcium::init_solvec()
{
  out << "<Cell_calcium> initialising solution vector..." << std::endl;
  int np = mesh->mesh_vals.vertices_count;

  solvec.resize(DIFVARS * np, 1); // NOTE: the variable ordering is c, ip, ce
  prev_solvec.resize(DIFVARS * np, 1);
  prev_solvec.block(0, 0, np, 1) = MatrixN1d().Constant(np, 1, p.at("c0"));
  prev_solvec.block(np, 0, np, 1) = MatrixN1d().Constant(np, 1, p.at("ip0"));
  prev_solvec.block(2 * np, 0, np, 1) = MatrixN1d().Constant(np, 1, p.at("ce0"));
  solvec = prev_solvec;

  //nd_solvec.resize(NONDIFVARS * np, 1); // NOTE: the variable ordering is g, h
  //prev_nd_solvec.resize(NONDIFVARS * np, 1);
  //prev_nd_solvec.block(0, 0, np, 1) = MatrixN1d().Constant(np, 1, p.at("g0"));
  //prev_nd_solvec.block(np, 0, np, 1) = MatrixN1d().Constant(np, 1, 0.0); // default to 0.0
  nd_solvec.resize(NONDIFVARS * np, 1); // variable h
  prev_nd_solvec.resize(NONDIFVARS * np, 1);
  prev_nd_solvec.block(0, 0, np, 1) = MatrixN1d().Constant(np, 1, p.at("h0")); 
  for (int n = 0; n < np; n++) {
    if (node_data(n, BOOL_apical) == 1.0) { // apical nodes only
      prev_nd_solvec(n) = p.at("h0");
    }
  }
  nd_solvec = prev_nd_solvec;
}

ArrayRefMass cCell_calcium::make_ref_mass()
{
  ArrayRefMass ref_mass;
  double v = (1.0 / 6.0) * 0.25 * 0.25;
  for (int i = 0; i < REF_MASS_SIZE; i++) {
    for (int j = 0; j < REF_MASS_SIZE; j++) { ref_mass(i, j) = v; }
  }
  return ref_mass;
}

void cCell_calcium::make_matrices()
{
  out << "<Cell_calcium> id:" << id << " calculating the spatial factors..." << std::endl;

  int nt = mesh->mesh_vals.tetrahedrons_count;
  element_data.resize(nt, Eigen::NoChange);
  int ns = mesh->mesh_vals.surface_triangles_count;
  surface_data.resize(ns, Eigen::NoChange);
  int np = mesh->mesh_vals.vertices_count;
  node_data.resize(np, Eigen::NoChange);

  // make the reference mass matrix
  ArrayRefMass ref_mass;
  ref_mass = make_ref_mass();

  // make the mass and stiffness matrices and the element constants matrix
  out << "<Cell_calcium> calculating the constant matrices..." << std::endl;
  MatrixNNd stiffc, stiffp, small_mass;
  MatrixNNd stiffce;
  stiffc = stiffc.Zero(np, np);
  stiffp = stiffp.Zero(np, np);
  stiffce = stiffce.Zero(np, np);
  small_mass = small_mass.Zero(np, np);

  // --------------------------------
  // for each volume element...
  // --------------------------------
  for (int n = 0; n < (mesh->mesh_vals.tetrahedrons_count); n++) {
    Eigen::Matrix<int, 1, 4> vi; // tetrahedron vertex indices
    vi = mesh->mesh_vals.tetrahedrons.block<1, 4>(n, 0);

    Eigen::Matrix<double, 4, 3> vert; // tetrahedron vertex coordinates
    for (int i = 0; i < 4; i++)
      vert.block<1, 3>(i, 0) = mesh->mesh_vals.vertices.block<1, 3>(int(vi(i)), 0); // why is typecast needed???

    Eigen::Matrix<double, 3, 3> J; // tetrahedron edge vectors
    for (int i = 0; i < 3; i++) J.block<1, 3>(i, 0) = vert.block<1, 3>(i + 1, 0) - vert.block<1, 3>(0, 0);
    double V, Vx6; // tetrahedron volume, (6x) volume
    Vx6 = J.determinant();
    V = Vx6 / 6.0;
    element_data(n, VOL_e) = V; // save the tetrahedron volume

    // RyR and PLC spatial factors per element
    //element_data(n, RYR_e) = ((mesh->e_dfa[n] < p.at("d_RyR")) ? (mesh->e_dfa[n] / p.at("d_RyR")) : 1.0);
    element_data(n, PLC_e) = (mesh->e_dfb[n] < p.at("PLCds") && mesh->e_dfa[n] > p.at("PLCdl")) ? 1.0 : 0.0;

    double Ic = V * p.at("Dc"); // diffusion coefficients
    double Ip = V * p.at("Dp");
    double Ice = V * p.at("De");

    Eigen::Matrix<double, 4, 4> M, C, G;
    M.col(0) << 1, 1, 1, 1;
    M.block<4, 3>(0, 1) = vert;
    C = M.inverse();
    G = C.block<3, 4>(1, 0).transpose() * C.block<3, 4>(1, 0); // gradients of the basis functions

    // construct the mass and stiffness matrix components
    for (int i = 0; i < 4; i++) {
      stiffc(vi(i), vi(i)) += G(i, i) * Ic;
      stiffp(vi(i), vi(i)) += G(i, i) * Ip;
      stiffce(vi(i), vi(i)) += G(i, i) * Ice;
      small_mass(vi(i), vi(i)) += ref_mass(i, i) * Vx6;
    }
    for (int i = 0; i < 3; i++) {
      stiffc(vi(0), vi(i + 1)) += G(0, i + 1) * Ic;
      stiffp(vi(0), vi(i + 1)) += G(0, i + 1) * Ip;
      stiffce(vi(0), vi(i + 1)) += G(0, i + 1) * Ice;
      small_mass(vi(0), vi(i + 1)) += ref_mass(0, i + 1) * Vx6;
      stiffc(vi(i + 1), vi(0)) = stiffc(vi(0), vi(i + 1));
      stiffp(vi(i + 1), vi(0)) = stiffp(vi(0), vi(i + 1));
      stiffce(vi(i + 1), vi(0)) = stiffce(vi(0), vi(i + 1));
      small_mass(vi(i + 1), vi(0)) = small_mass(vi(0), vi(i + 1));
    }
    for (int i = 0; i < 2; i++) {
      stiffc(vi(1), vi(i + 2)) += G(1, i + 2) * Ic;
      stiffp(vi(1), vi(i + 2)) += G(1, i + 2) * Ip;
      stiffce(vi(1), vi(i + 2)) += G(1, i + 2) * Ice;
      small_mass(vi(1), vi(i + 2)) += ref_mass(1, i + 2) * Vx6;
      stiffc(vi(i + 2), vi(1)) = stiffc(vi(1), vi(i + 2));
      stiffp(vi(i + 2), vi(1)) = stiffp(vi(1), vi(i + 2));
      stiffce(vi(i + 2), vi(1)) = stiffce(vi(1), vi(i + 2));
      small_mass(vi(i + 2), vi(1)) = small_mass(vi(1), vi(i + 2));
    }
    stiffc(vi(2), vi(3)) += G(2, 3) * Ic;
    stiffp(vi(2), vi(3)) += G(2, 3) * Ip;
    stiffce(vi(2), vi(3)) += G(2, 3) * Ice;
    small_mass(vi(2), vi(3)) += ref_mass(2, 3) * Vx6;
    stiffc(vi(3), vi(2)) = stiffc(vi(2), vi(3));
    stiffp(vi(3), vi(2)) = stiffp(vi(2), vi(3));
    stiffce(vi(3), vi(2)) = stiffce(vi(2), vi(3));
    small_mass(vi(3), vi(2)) = small_mass(vi(2), vi(3));
  }

  // construct sparse mass matrix from a list of triplets (non zero elements)
  std::vector<Triplet> triplet_list;
  int np2 = np * 2;
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      double v_ij = small_mass(i, j);
      if (v_ij != 0) { // add non zeros in first, second and third blocks
        triplet_list.push_back(Triplet(i, j, v_ij));
        triplet_list.push_back(Triplet(np + i, np + j, v_ij));
        triplet_list.push_back(Triplet(np2 + i, np2 + j, v_ij));
      }
    }
  }
  sparseMass.resize(DIFVARS * np, DIFVARS * np);
  sparseMass.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // construct sparse stiffness matrix from list of triplets
  triplet_list.clear();
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      double v_stiffc = stiffc(i, j); // first block
      if (v_stiffc != 0) { triplet_list.push_back(Triplet(i, j, v_stiffc)); }
      double v_stiffp = stiffp(i, j); // second block
      if (v_stiffp != 0) { triplet_list.push_back(Triplet(np + i, np + j, v_stiffp)); }
      double v_stiffce = stiffce(i, j); // third block
      if (v_stiffce != 0) { triplet_list.push_back(Triplet(np2 + i, np2 + j, v_stiffce)); }
    }
  }
  sparseStiff.resize(DIFVARS * np, DIFVARS * np);
  sparseStiff.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // set the A matrix size
  sparseA.resize(DIFVARS * np, DIFVARS * np);

  // --------------------------------
  // for each surface triangle...
  // --------------------------------
  for (int n = 0; n < ns; n++) {
    Eigen::Array<int, 1, 3> vi; // surface triangle vertex indices
    vi = mesh->mesh_vals.surface_triangles.block<1, 3>(n, 0);
    Eigen::Array<double, 3, 3> vert; // triangle vertex coordinates
    for (int i = 0; i < 3; i++)
      vert.block<1, 3>(i, 0) = mesh->mesh_vals.vertices.block<1, 3>(int(vi(i)), 0); // why is typecast needed???
    Eigen::Matrix<double, 1, 3> side1 = vert.block<1, 3>(0, 0) - vert.block<1, 3>(1, 0);
    Eigen::Matrix<double, 1, 3> side2 = vert.block<1, 3>(0, 0) - vert.block<1, 3>(2, 0);
    double triangle_area = 0.5 * (side1.cross(side2)).norm();
    surface_data(n, AREA_s) = triangle_area; // save the triangle area
  }
  // --------------------------------
  // for each node...
  // --------------------------------
  for (int n = 0; n < np; n++) {
    node_data(n, BOOL_apical) = 0.0; // default
  }
  for (int n = 0; n < mesh->apical_triangles_count; n++) { // for each apical (surface) element...
    Eigen::Array<int, 1, 3> vi;                            // apical triangle vertex indices
    vi = mesh->mesh_vals.surface_triangles.block<1, 3>(mesh->apical_triangles(n), 0);
    for (int i = 0; i < 3; i++) {          // for each apical triangle vertex
      node_data(vi(i), BOOL_apical) = 1.0; // flag it as apical
    }
  }
}

Array1VC cCell_calcium::get_apical_reactions(double c, double ip, double ce, double h)
{
  double phi_c = pow(c, 4) / (pow(p.at("K_c"), 4) + pow(c, 4));
  double phi_p = pow(ip, 2) / (pow(p.at("K_p"), 2) + pow(ip, 2));
  double h_alpha = pow(p.at("K_h"), 4) / (pow(p.at("K_h"), 4) + pow(c, 4));
  double beta = phi_c * phi_p * h;
  double alpha = (1 - phi_p) * (1 - (phi_c * h_alpha));
  double po = beta / (beta + p.at("k_beta") * (beta + alpha));

  Array1VC reactions;
  reactions(0) = p.at("k_IPR") * po * (ce - c);
  reactions(1) = ip;
  reactions(2) = -reactions(0) / p.at("Gamma");
  return reactions;
}

//double cCell_calcium::get_g_reaction(double c, double g)
//{ // RYR dynamics
//  double ginf = pow(p.at("K_hRyR"), 2) / (pow(p.at("K_hRyR"), 2) + pow(c, 2));
//  return (ginf - g) / p.at("tau");
//}

double cCell_calcium::get_h_reaction(double c, double h)
{ // IPR dynamics
  double hinf = pow(p.at("K_h"), 4) / (pow(p.at("K_h"), 4) + pow(c, 4));
  double htau = p.at("tau_max") * pow(p.at("K_tau"), 4) / (pow(p.at("K_tau"), 4) + pow(c, 4));
  return (hinf - h) / htau;
}

//Array1VC cCell_calcium::get_body_reactions(double c, double ip, double ce, double g, double ryr_f, double plc_f)
Array1VC cCell_calcium::get_body_reactions(double c, double ip, double ce, double plc_f)
{
  double J_SERCA = p.at("V_p") * (pow(c, 2) - p.at("K_bar") * pow(ce, 2)) / (pow(p.at("k_p"), 2) + pow(c, 2));
  //double J_RYR = ryr_f * p.at("V_RyR") * (pow(c, p.at("n_RyR")) / (pow(c, p.at("n_RyR")) + pow(p.at("K_RyR"), p.at("n_RyR")))) *
  //               (pow(ce, p.at("m_RyR")) / (pow(ce, p.at("m_RyR")) + pow(p.at("K_RyR2"), p.at("m_RyR")))) * g;
  double J_RYR = 0.0;
  double vplc = plc_f * p.at("V_PLC") * pow(c, 2) / (pow(p.at("K_PLC"), 2) + pow(c, 2));
  double vdeg = (p.at("V_5K") + (p.at("V_3K") * pow(c, 2)) / (pow(p.at("K3K"), 2) + pow(c, 2))) * ip;

  Array1VC reactions;
  reactions(0) = (J_RYR * (ce - c)) - J_SERCA;
  reactions(1) = vplc - vdeg;
  reactions(2) = -reactions(0) / p.at("Gamma");

  // scale reaction terms by (at_rest_volume / new_volume)
  reactions *= volume_term1;

  return reactions;
}

//MatrixN1d cCell_calcium::make_load(double dt, bool plc)
MatrixN1d cCell_calcium::make_load(bool plc)
{
  int np = mesh->mesh_vals.vertices_count;
  ArrayX1C c, ip, g, h;
  ArrayX1C load_c, load_ip;
  ArrayX1C ce;
  ArrayX1C load_ce;
  MatrixN1d load;

  c = prev_solvec.block(0, 0, np, 1);
  ip = prev_solvec.block(np, 0, np, 1);
  ce = prev_solvec.block(2 * np, 0, np, 1);
  //g = prev_nd_solvec.block(0, 0, np, 1);
  //h = prev_nd_solvec.block(np, 0, np, 1); // note: only apical (surface) nodes used
  h = prev_nd_solvec.block(0, 0, np, 1); // note: only apical (surface) nodes used

  load_c = load_c.Zero(np, 1);
  load_ip = load_ip.Zero(np, 1);
  load_ce = load_ce.Zero(np, 1);
  load.resize(DIFVARS * np, Eigen::NoChange);

  // volume reaction terms
  for (int n = 0; n < (mesh->mesh_vals.tetrahedrons_count); n++) { // for each volume element...
    Eigen::Array<int, 1, 4> vi;                                    // tetrahedron vertex indices
    vi = mesh->mesh_vals.tetrahedrons.block<1, 4>(n, 0);

    double cav = 0.25 * (c(vi(0)) + c(vi(1)) + c(vi(2)) + c(vi(3)));
    double ipav = 0.25 * (ip(vi(0)) + ip(vi(1)) + ip(vi(2)) + ip(vi(3)));
    double ceav = 0.25 * (ce(vi(0)) + ce(vi(1)) + ce(vi(2)) + ce(vi(3)));
    //double gav = 0.25 * (g(vi(0)) + g(vi(1)) + g(vi(2)) + g(vi(3)));

    //Array1VC reactions =
    //  get_body_reactions(cav, ipav, ceav, gav, double(element_data(n, RYR_e)), double(plc ? element_data(n, PLC_e) : 0.0));
    Array1VC reactions =
      get_body_reactions(cav, ipav, ceav, double(plc ? element_data(n, PLC_e) : 0.0));

    for (int i = 0; i < 4; i++) {                                    // for each tetrahedron vertex
      // reaction terms, scaled by 1/4 volume
      load_c(vi(i)) += element_data(n, VOL_e) * 0.25 * (reactions(0) - volume_term2 * c(vi(i)));
      load_ip(vi(i)) += element_data(n, VOL_e) * 0.25 * (reactions(1) - volume_term2 * ip(vi(i)));
      load_ce(vi(i)) += element_data(n, VOL_e) * 0.25 * (reactions(2) - volume_term2 * ce(vi(i)));
    }
  }

  // apical (surface) reaction terms
  for (int n = 0; n < (mesh->apical_triangles_count); n++) { // for each apical (surface) triangle...
    Eigen::Array<int, 1, 3> vi;                              // apical triangle vertex indices
    vi = mesh->mesh_vals.surface_triangles.block<1, 3>(mesh->apical_triangles(n), 0);

    double cav = (c(vi(0)) + c(vi(1)) + c(vi(2))) / 3.0;
    double ipav = (ip(vi(0)) + ip(vi(1)) + ip(vi(2))) / 3.0;
    double ceav = (ce(vi(0)) + ce(vi(1)) + ce(vi(2))) / 3.0;
    double hav = (h(vi(0)) + h(vi(1)) + h(vi(2))) / 3.0;

    Array1VC reactions = get_apical_reactions(cav, ipav, ceav, hav);
    for (int i = 0; i < 3; i++) { // for each apical triangle vertex
      load_c(vi(i)) +=            // reaction term scaled by 1/3 area
        (surface_data(mesh->apical_triangles(n), AREA_s) / 3.0) * reactions(0);
      load_ce(vi(i)) += // reaction term scaled by 1/3 area
        (surface_data(mesh->apical_triangles(n), AREA_s) / 3.0) * reactions(2);
    }
  }

  // the diffusing variables
  load.block(0, 0, np, 1) = load_c;
  load.block(np, 0, np, 1) = load_ip;
  // load.block(np, 0, np, 1) = load_ip + exchange_load_ip;
  load.block(2 * np, 0, np, 1) = load_ce;

  return load;
}

MatrixN1d cCell_calcium::solve_nd(double dt)
{ // the non-diffusing variables
  int np = mesh->mesh_vals.vertices_count;
  ArrayX1C c, g, h;
  MatrixN1d svec;

  c = prev_solvec.block(0, 0, np, 1);
  //g = prev_nd_solvec.block(0, 0, np, 1);
  //h = prev_nd_solvec.block(np, 0, np, 1); // note: only apical (surface) nodes used
  h = prev_nd_solvec.block(0, 0, np, 1); // note: only apical (surface) nodes used
  svec.resize(NONDIFVARS * np, Eigen::NoChange);

  for (int n = 0; n < np; n++) {                               // for each node...
    //svec(n) = g(n) + (dt * get_g_reaction(c(n), g(n)));        // g
    if (node_data(n, BOOL_apical) == 1.0) {                    // only the apical nodes
      //svec(np + n) = h(n) + (dt * get_h_reaction(c(n), h(n))); // h
	  svec(n) = h(n) + (dt * get_h_reaction(c(n), h(n))); // h
    }
  }
  return svec;
}

void cCell_calcium::exchange()
{
  float* msg;
  MPI_Status stat;
  MPI_Request request;
  // send common face values to other cells (non-blocking)
  for (std::vector<cfc>::iterator it = cells.begin(); it != cells.end(); ++it) {
    int mlength = C2CCOUNT * it->fcount;
    msg = new float[mlength];
    //...
    MPI_CHECK(MPI_Isend(msg, mlength, MPI_FLOAT, it->cell + 1, CELL_CELL_TAG, MPI_COMM_WORLD, &request));
    delete[] msg;
  }
  // receive common face values back from other cells
  for (std::vector<cfc>::iterator it = cells.begin(); it != cells.end(); ++it) {
    int mlength = C2CCOUNT * it->fcount;
    msg = new float[mlength];
    MPI_CHECK(MPI_Recv(msg, mlength, MPI_FLOAT, it->cell + 1, CELL_CELL_TAG, MPI_COMM_WORLD, &stat));
    //...
    delete[] msg;
  }
}

void cCell_calcium::run()
{
  float msg[ACCOUNT]; // mpi message from and back to acinus
  double delta_time, current_time;
  double prev_delta_time = 0.0;
  MPI_Status stat;
  struct timespec start, end, flow_start, flow_end;
  double elapsed, flow_elapsed;
  int np = mesh->mesh_vals.vertices_count;
  double V0 = element_data.col(VOL_e).sum();  // initial volume

  MatrixN1d rhs; // the right-hand-side vector
  rhs.resize(DIFVARS * np, Eigen::NoChange);
  bool plc;

  int step = 0;
  while (true) {
    step++;
    prev_solvec = solvec;       // c, ip, ce
    //prev_nd_solvec = nd_solvec; // g, h
    prev_nd_solvec = nd_solvec; // h

    // get current time and time step value from acinus
    MPI_CHECK(MPI_Recv(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD, &stat));
    delta_time = msg[dTime];
    current_time = msg[cTime];
    if (delta_time == 0.0) { // done?
      MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD));
      break;
    }

    // not done
    out << std::fixed << std::setprecision(3);
    out << "<Cell_calcium> step: " << step << " current_time: " << current_time << "s";
    out << " delta_time: " << delta_time << "s" << std::endl;
    plc = ((current_time >= p.at("PLCsrt")) and (current_time <= p.at("PLCfin"))); // PLC on or off?

    // coupling flow back to calcium model
    volume_term1 = V0 / flow->prev_solvec(VOL); // at rest volume divided by current volume
    volume_term2 = flow->prev_dsolvec(VOL) / flow->prev_solvec(VOL);  // volume derivative divided by current volume

    // calculate the fluid flow
    clock_gettime(CLOCK_REALTIME, &flow_start);
    if(p.at("fluidFlow")) flow->step(current_time, delta_time);
    clock_gettime(CLOCK_REALTIME, &flow_end);
    flow_elapsed = (flow_end.tv_sec - flow_start.tv_sec) + ((flow_end.tv_nsec - flow_start.tv_nsec) / 1e9);
	
    if (delta_time != prev_delta_time) { // recalculate A matrix if time step changed
      sparseA = sparseMass + (delta_time * sparseStiff);
      solver.compute(sparseA);
      if (solver.info() != Eigen::Success) { // decomposition failed?
        utils::fatal_error("matrix decomposition failed", out);
      }
      prev_delta_time = delta_time;
    }

    // exchange common face data
    exchange();

    // calculate solution for diffusing variables
    //rhs = (sparseMass * prev_solvec) + (delta_time * make_load(delta_time, plc));
    rhs = (sparseMass * prev_solvec) + (delta_time * make_load(plc));
    clock_gettime(CLOCK_REALTIME, &start);
    solvec = solver.solve(rhs); // Eigen solver
    if (solver.info() != Eigen::Success) {
      utils::fatal_error("solver failed", out);
      ;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1e9);
    out << std::fixed << std::setprecision(3);
    out << "<Cell_calcium> flow duration: " << flow_elapsed << "s ; solver duration: " << elapsed << "s" << std::endl;

    // check solver error and send it to acinus
    // ...
    msg[sError] = 0.0;
    MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD));

    // calculate solution for non-diffusing variables
    nd_solvec = solve_nd(delta_time);

    // save results
    if (step % int(p.at("Tstride")) == 0) {
      save_results(ca_file, 0);  // 0 = calcium
      save_results(ip3_file, 1); // 1 = ip3
      save_results(cer_file, 2); // 2 = cer
      flow->save_results();
    }
  }
}

void cCell_calcium::save_results(std::ofstream& data_file, int var)
{
  int np = mesh->mesh_vals.vertices_count;
  float* fbuf = new float[np];
  for (int n = 0; n < np; n++) fbuf[n] = solvec[var * np + n]; // convert to float for reduced file size
  data_file.write(reinterpret_cast<char*>(fbuf), np * sizeof(float));
}
