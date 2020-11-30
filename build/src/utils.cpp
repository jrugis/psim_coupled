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
    out << " " << tokens[0];
  }
  out << std::endl;
  model_file.close();
  
  // add dependant parameters to the map
  p["ce0"] = (p.at("ct") - p.at("c0")) / p.at("Gamma");
  p["h0"] = pow(p.at("K_h"), 4) / (pow(p.at("K_h"), 4) + pow(p.at("c0"), 4));
  p["He"] = 1e3 * pow(10, -p.at("pHe"));
  p["Cll0"] = p.at("Nal0") = p.at("Kl0");
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

  // add dependant default parameters to the map if not already defined
  if (not p.count("ip0"))
    p["ip"] = 0.5 * (2e-4 / 0.1) * ((pow(p.at("K3K"), 2) + pow(p.at("c0"), 2)) / pow(p.at("c0"), 2));
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
