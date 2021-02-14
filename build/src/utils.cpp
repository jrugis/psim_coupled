/*
 * utils.cpp
 *
 *	Created on: 27/04/2018
 *	Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <fstream>
#include <iostream>

#include "global_defs.hpp"
#include "utils.hpp"

// get the center points of triangles
void utils::calc_tri_centers(MatrixN3d& centers, const MatrixN3d& vertices, const MatrixN3i& triangles)
{
  centers.resize(triangles.rows(), Eigen::NoChange);
  for (int n = 0; n < triangles.rows(); n++) {
    Vector3d v1 = Vector3d(vertices.row((triangles)(n, 0)));
    Vector3d v2 = Vector3d(vertices.row((triangles)(n, 1)));
    Vector3d v3 = Vector3d(vertices.row((triangles)(n, 2)));
    centers.row(n) = (v1 + v2 + v3) / 3.0;
  }
}

// NOTE: outputs error message to stderr and the process "out" file
void utils::fatal_error(const std::string msg, std::ofstream& out)
{
  std::string m = "ERROR: " + msg;
  out << m << std::endl;
  out.close();
  std::cerr << m << std::endl;
  MPI_CHECK(MPI_Abort(MPI_COMM_WORLD, 1)); // make sure the whole program fails
  exit(1);
}

// NOTE: used by each of the acinus, lumen and cell objects
void utils::get_parameters(const std::string file_id, int cell_num, std::unordered_map<std::string, double>& p, std::ofstream& out)
{
  std::string file_name = file_id + ".dat";
  std::ifstream model_file(file_name); // open the model parameters file
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line

  if (not model_file.is_open()) { fatal_error("the model parameters file " + file_name + " could not be opened", out); }
  out << "<utils> reading model parameters:";
  while (getline(model_file, line)) {
    if (line.data()[0] == '#') continue;
    int ci = line.find_first_of("#");      // remove comment, if any
    if (ci > 0) line = line.substr(0, ci); //
    line = boost::trim_right_copy(line);   // remove trailing whitespace
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    p[tokens[0]] = atof(tokens[((tokens.size() == 2) ? 1 : cell_num)].c_str());
    out << " " << tokens[0] << std::flush;
  }
  out << std::endl;
  model_file.close();
};

double utils::get_distance(const Vector3d& p, const Vector3d& v, const Vector3d& w)
{
  // Return minimum distance between line segment vw and point p
  double l2 = (w - v).squaredNorm();      // |w-v|^2		avoid a sqrt
  if (l2 == 0.0) return ((v - p).norm()); // v == w case, return distance(p, v)
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // Find projection of point p onto the line. It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // Clamp t from [0,1] to handle points outside the segment vw.
  double t = std::max(0.0, std::min(1.0, (p - v).dot(w - v) / l2)); // max(0, min(1, dot(p - v, w - v) / l2));
  const Vector3d projection = v + (t * (w - v));                    // Projection falls on the segment
  return ((projection - p).norm());                                 // return distance(p, projection)
}

void utils::read_mesh(const std::string file_name, sMeshVals& mesh_vals, std::ofstream& out)
{
  std::ifstream cell_file(file_name + ".bmsh", std::ios::in | std::ios::binary); // open the mesh file
  if (not cell_file.is_open()) { fatal_error("mesh file " + file_name + " could not be opened", out); }

  // get the mesh vertices (int count, 3x-double vertices)
  cell_file.read(reinterpret_cast<char*>(&(mesh_vals.vertices_count)), sizeof(int));
  mesh_vals.vertices.resize(mesh_vals.vertices_count, Eigen::NoChange);
  cell_file.read(reinterpret_cast<char*>(mesh_vals.vertices.data()), 3 * mesh_vals.vertices_count * sizeof(double));

  // get the surface triangles (int count, 3x-int vertex indices)
  cell_file.read(reinterpret_cast<char*>(&(mesh_vals.surface_triangles_count)), sizeof(int));
  mesh_vals.surface_triangles.resize(mesh_vals.surface_triangles_count, Eigen::NoChange);
  cell_file.read(reinterpret_cast<char*>(mesh_vals.surface_triangles.data()), 3 * mesh_vals.surface_triangles_count * sizeof(int));

  // get the element tetrahedrons (int count, 4x-int vertex indices)
  cell_file.read(reinterpret_cast<char*>(&(mesh_vals.tetrahedrons_count)), sizeof(int));
  mesh_vals.tetrahedrons.resize(mesh_vals.tetrahedrons_count, Eigen::NoChange);
  cell_file.read(reinterpret_cast<char*>(mesh_vals.tetrahedrons.data()), 4 * mesh_vals.tetrahedrons_count * sizeof(int));

  cell_file.close();
}

void utils::save_matrix(const std::string file_name, int bytes, char* data)
{
  std::ofstream data_file;
  data_file.open(file_name, std::ios::binary);
  data_file.write(data, bytes);
  data_file.close();
}
