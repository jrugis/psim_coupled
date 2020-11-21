/*
 * global_defs.hpp
 *
 *  Created on: 20/11/20
 *      Author: jrugis
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <mpi.h>
#include <string>

//************************************************************************
// model sizes
//  NOTE: a single acinus with lumen (for now)
#define CELLS_COUNT 7

//************************************************************************
// MPI defs, macros

#define MPI_ERROR_ABORT -101
#define MPI_NODES_ABORT -102
#define ACINUS_CELL_TAG 400
#define CELL_CELL_TAG 401

// mpi shutdown on error
#define mpi_abort(err)                                                \
  {                                                                   \
    std::cout << "Program shutdown error code: " << err << std::endl; \
    MPI_Abort(MPI_COMM_WORLD, err);                                   \
  }

// mpi error checking macro
#define MPI_CHECK(call)                                      \
  if ((call) != MPI_SUCCESS) {                               \
    std::cout << "MPI error calling: " << call << std::endl; \
    mpi_abort(MPI_ERROR_ABORT);                              \
  }

//************************************************************************
// acinus to cell mpi message
// delta time, current time, solver error
enum acinus2cell { dTime, cTime, sError, ACCOUNT };

//************************************************************************
// acinus to lumen mpi message
// delta time, current time, apical calcium, basal calcium, cell volume
// enum acinus2lumen { dlTime, clTime, aC, Bc, cV, ALCOUNT };

//************************************************************************
// cell to cell mpi message
// triangle index, value
enum cell2cell { tInd, tVal, C2CCOUNT };

//************************************************************************
// cell to cell connectivity
// this_triangle, other_cell, other_triangle
enum cell_conn { tTri, oCell, oTri, CCONNCOUNT };

//************************************************************************
// some convenience typedefs
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixN1d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixNNd;
typedef Eigen::SparseMatrix<double> SparceMatrixd;

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> MatrixN1i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajorBit> MatrixN2i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> MatrixN4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, CCONNCOUNT, Eigen::RowMajorBit> MatrixNCi;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixNNi;

typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector4i Vector4i;

//************************************************************************
//************************************************************************
// cell mesh values
struct sMeshVals {
  int vertices_count;
  int surface_triangles_count;
  int tetrahedrons_count;
  MatrixN3d vertices;          // 3x coordinate
  MatrixN3i surface_triangles; // 3x vertex
  MatrixN4i tetrahedrons;      // 4x vertex
};

//************************************************************************
// thermodynamic constants
#define R 8.314462100000000
#define T 310
#define F 9.645833650000000e4
const double RTF = 1000.0 * R * T / F;

//************************************************************************
//************************************************************************

#endif /* DEFS_H_ */
