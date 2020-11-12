/*
 * cAcinus.cpp
 *
 *  Created on: 11/11/2020
 *      Author: jrugis
 */

#include <time.h>
#include <iomanip>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cAcinus.hpp"
#include "cLumen.hpp"

cAcinus::cAcinus(std::string host_name, int rank, int c_rank, int c_count) {
  my_rank = rank;
  cell_rank = c_rank;
  cell_count = c_count;

  id = "a1";             // "a1" hard coded for now

  out.open(id + ".out");
  out << "<Acinus> id: " << id << std::endl;
  out << "<Acinus> host_name: " << host_name << std::endl;

  utils::get_parameters(id, calciumParms, 1, p, out);

  lumen = new cLumen(this, "l1"); // "l1" hard coded for now
}

cAcinus::~cAcinus() {
  out.close();
  delete lumen;
}

// NOTE: mpi send to all first, then receive from all
double cAcinus::snd_recv(double t, double dt) {
  float msg[ACCOUNT]; 
  MPI_Status stat;

  out << "<Acinus> t: " << t << std::endl;
  msg[dTime] = dt;
  msg[cTime] = t;
  msg[sError] = 0.0;
  for(int n = cell_rank; n < (cell_rank + cell_count); n++){
    MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, n, ACINUS_CELL_TAG, MPI_COMM_WORLD));
  }
  for(int n = cell_rank; n < (cell_rank + cell_count); n++){
    MPI_CHECK(MPI_Recv(&msg, ACCOUNT, MPI_FLOAT, n, ACINUS_CELL_TAG, MPI_COMM_WORLD, &stat));
  }
  return(msg[sError]);
}

void cAcinus::run() {
  double t = 0.0;
  double solver_dt = p[delT];
  double prev_dt = solver_dt;
  double error;
  struct timespec start, end;
  double elapsed;

  // simulation time stepping and synchronization
  clock_gettime(CLOCK_REALTIME, &start);
  while((p[totalT] - t) > 0.000001 ) {  // HARD CODED: assumes solver_dt always > 1us
    error = snd_recv(t, solver_dt);  // invoke the calcium solver
    if(error != 0.0) { // change time step?
      // ...
    } 
    // ... // invoke the fluid flow solver
    clock_gettime(CLOCK_REALTIME, &end);
	elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
	out << std::fixed << std::setprecision(3);
	out << "<Acinus> step duration: " << elapsed << "s"<< std::endl;
    start = end;
    t += solver_dt;
  }  

  // instruct cells to finish
  solver_dt = 0.0;
  snd_recv(t, solver_dt);
}


