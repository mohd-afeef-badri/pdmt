/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     VTK writers for explicit three-dimensional polyhedron topology.

*****************************************************************************/

#ifndef PDMT_POLY_MESH_3D_WRITE_HPP
#define PDMT_POLY_MESH_3D_WRITE_HPP

#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace Pdmt3DWriter {

inline void validateTopology(const KNM<double> &nodes,
                             const KN<KN<long> > &cells,
                             const KN<KN<long> > &faces) {
  for (long f = 0; f < faces.N(); ++f) {
    if (faces[f].N() < 3)
      ExecError("PdmtPolyMeshWrite: a polygonal face must contain at least three nodes");
    for (long i = 0; i < faces[f].N(); ++i)
      if (faces[f][i] < 0 || faces[f][i] >= nodes.N())
        ExecError("PdmtPolyMeshWrite: a face references an out-of-range node");
  }
  for (long c = 0; c < cells.N(); ++c) {
    if (cells[c].N() < 4)
      ExecError("PdmtPolyMeshWrite: a polyhedron must contain at least four faces");
    for (long f = 0; f < cells[c].N(); ++f) {
      const long id = std::labs(cells[c][f]) - 1;
      if (cells[c][f] == 0 || id < 0 || id >= faces.N())
        ExecError("PdmtPolyMeshWrite: a 3D cell references an out-of-range face");
    }
  }
}

inline long faceIndex(long encoded) {
  if (encoded == 0)
    ExecError("PdmtPolyMeshWrite: a 3D cell contains invalid face id 0");
  return std::labs(encoded) - 1;
}

inline void orientedFace(long encoded, const KN<KN<long> > &faces,
                         std::vector<long> &vertices) {
  const long id = faceIndex(encoded);
  if (id < 0 || id >= faces.N())
    ExecError("PdmtPolyMeshWrite: a 3D cell references an out-of-range face");
  vertices.resize(faces[id].N());
  if (encoded > 0) {
    for (long i = 0; i < faces[id].N(); ++i)
      vertices[i] = faces[id][i];
  } else {
    for (long i = 0; i < faces[id].N(); ++i)
      vertices[i] = faces[id][faces[id].N() - 1 - i];
  }
}

inline std::vector<long> cellVertices(const KN<long> &cell,
                                      const KN<KN<long> > &faces) {
  std::set<long> unique;
  for (long i = 0; i < cell.N(); ++i) {
    const long id = faceIndex(cell[i]);
    if (id < 0 || id >= faces.N())
      ExecError("PdmtPolyMeshWrite: a 3D cell references an out-of-range face");
    for (long j = 0; j < faces[id].N(); ++j)
      unique.insert(faces[id][j]);
  }
  return std::vector<long>(unique.begin(), unique.end());
}

inline long faceStreamSize(const KN<long> &cell,
                           const KN<KN<long> > &faces) {
  long size = 1;
  for (long i = 0; i < cell.N(); ++i) {
    const long id = faceIndex(cell[i]);
    if (id < 0 || id >= faces.N())
      ExecError("PdmtPolyMeshWrite: a 3D cell references an out-of-range face");
    size += 1 + faces[id].N();
  }
  return size;
}

inline void writeFaceFieldVtu(std::ofstream &out, const KN<KN<long> > &faces,
                              const KN<long> *faceLabels) {
  out << "    <FieldData>\n"
      << "      <DataArray type=\"Int64\" Name=\"pdmt_face_connectivity\" format=\"ascii\">\n        ";
  for (long i = 0; i < faces.N(); ++i) {
    for (long j = 0; j < faces[i].N(); ++j)
      out << faces[i][j] << " ";
  }
  out << "\n      </DataArray>\n";
  out << "      <DataArray type=\"Int64\" Name=\"pdmt_face_offsets\" format=\"ascii\">\n        ";
  long offset = 0;
  for (long i = 0; i < faces.N(); ++i) {
    offset += faces[i].N();
    out << offset << " ";
  }
  out << "\n      </DataArray>\n";
  if (faceLabels) {
    if (faceLabels->N() != faces.N())
      ExecError("PdmtPolyMeshWrite: faceLabels must have one value per face");
    out << "      <DataArray type=\"Int64\" Name=\"pdmt_face_labels\" format=\"ascii\">\n        ";
    for (long i = 0; i < faceLabels->N(); ++i)
      out << (*faceLabels)[i] << " ";
    out << "\n      </DataArray>\n";
  }
  out << "    </FieldData>\n";
}

} // namespace Pdmt3DWriter

inline void writePolyVtu3D(const std::string *fileName, KNM<double> *nodes,
                           KN<KN<long> > *cells, KN<KN<long> > *faces,
                           KN<long> *labels, KN<long> *faceLabels) {
  using namespace Pdmt3DWriter;
  if (!nodes || !cells || !faces || nodes->M() < 3)
    ExecError("PdmtPolyMeshWrite: invalid 3D polyhedron arrays");
  validateTopology(*nodes, *cells, *faces);
  if (labels && labels->N() != cells->N())
    ExecError("PdmtPolyMeshWrite: labels must have one value per polyhedron");

  std::ofstream out(fileName->c_str());
  if (!out)
    ExecError("PdmtPolyMeshWrite: cannot open the requested VTU file");

  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
      << "  <UnstructuredGrid>\n";

  writeFaceFieldVtu(out, *faces, faceLabels);
  out << "    <Piece NumberOfPoints=\"" << nodes->N() << "\" NumberOfCells=\"" << cells->N() << "\">\n";
  if (labels) {
    out << "      <CellData Scalars=\"label\">\n"
        << "        <DataArray type=\"Int64\" Name=\"label\" format=\"ascii\">\n          ";
    for (long i = 0; i < labels->N(); ++i)
      out << (*labels)[i] << " ";
    out << "\n        </DataArray>\n      </CellData>\n";
  }

  out << "      <Points>\n"
      << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (long i = 0; i < nodes->N(); ++i)
    out << "          " << (*nodes)(i, 0L) << " " << (*nodes)(i, 1L) << " " << (*nodes)(i, 2L) << "\n";
  out << "        </DataArray>\n      </Points>\n      <Cells>\n";

  out << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
  std::vector<std::vector<long> > uniqueVertices(cells->N());
  for (long c = 0; c < cells->N(); ++c) {
    uniqueVertices[c] = cellVertices((*cells)[c], *faces);
    for (std::vector<long>::const_iterator it = uniqueVertices[c].begin();
         it != uniqueVertices[c].end(); ++it)
      out << *it << " ";
    out << "\n          ";
  }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
  long offset = 0;
  for (long c = 0; c < cells->N(); ++c) {
    offset += uniqueVertices[c].size();
    out << offset << " ";
  }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
  for (long c = 0; c < cells->N(); ++c)
    out << "42 ";
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";

  std::vector<long> oriented;
  for (long c = 0; c < cells->N(); ++c) {
    out << "          " << (*cells)[c].N() << " ";
    for (long f = 0; f < (*cells)[c].N(); ++f) {
      orientedFace((*cells)[c][f], *faces, oriented);
      out << oriented.size() << " ";
      for (std::vector<long>::const_iterator it = oriented.begin(); it != oriented.end(); ++it)
        out << *it << " ";
    }
    out << "\n";
  }
  out << "        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n          ";
  offset = 0;
  for (long c = 0; c < cells->N(); ++c) {
    offset += faceStreamSize((*cells)[c], *faces);
    out << offset << " ";
  }
  out << "\n        </DataArray>\n"
      << "      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
}

inline void writePolyVtk3D(const std::string *fileName, KNM<double> *nodes,
                           KN<KN<long> > *cells, KN<KN<long> > *faces,
                           KN<long> *labels, KN<long> *faceLabels) {
  using namespace Pdmt3DWriter;
  if (!nodes || !cells || !faces || nodes->M() < 3)
    ExecError("PdmtPolyMeshWrite: invalid 3D polyhedron arrays");
  validateTopology(*nodes, *cells, *faces);
  if (labels && labels->N() != cells->N())
    ExecError("PdmtPolyMeshWrite: labels must have one value per polyhedron");
  if (faceLabels && faceLabels->N() != faces->N())
    ExecError("PdmtPolyMeshWrite: faceLabels must have one value per face");

  std::ofstream out(fileName->c_str());
  if (!out)
    ExecError("PdmtPolyMeshWrite: cannot open the requested VTK file");

  out << "# vtk DataFile Version 2.0\nPDMT 3D polyhedral mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n\n"
      << "FIELD FieldData " << (faceLabels ? 3 : 2) << "\n";
  long faceConnectivitySize = 0;
  for (long f = 0; f < faces->N(); ++f)
    faceConnectivitySize += (*faces)[f].N();
  out << "pdmt_face_connectivity 1 " << faceConnectivitySize << " long\n";
  for (long f = 0; f < faces->N(); ++f)
    for (long i = 0; i < (*faces)[f].N(); ++i)
      out << (*faces)[f][i] << "\n";
  out << "pdmt_face_offsets 1 " << faces->N() << " long\n";
  long faceOffset = 0;
  for (long f = 0; f < faces->N(); ++f) {
    faceOffset += (*faces)[f].N();
    out << faceOffset << "\n";
  }
  if (faceLabels) {
    out << "pdmt_face_labels 1 " << faceLabels->N() << " long\n";
    for (long f = 0; f < faceLabels->N(); ++f)
      out << (*faceLabels)[f] << "\n";
  }

  out << "\nPOINTS " << nodes->N() << " double\n";
  for (long i = 0; i < nodes->N(); ++i)
    out << (*nodes)(i, 0L) << " " << (*nodes)(i, 1L) << " " << (*nodes)(i, 2L) << "\n";

  long connectivitySize = 0;
  for (long c = 0; c < cells->N(); ++c)
    connectivitySize += 1 + faceStreamSize((*cells)[c], *faces);
  out << "\nCELLS " << cells->N() << " " << connectivitySize << "\n";

  std::vector<long> oriented;
  for (long c = 0; c < cells->N(); ++c) {
    out << faceStreamSize((*cells)[c], *faces) << " " << (*cells)[c].N() << " ";
    for (long f = 0; f < (*cells)[c].N(); ++f) {
      orientedFace((*cells)[c][f], *faces, oriented);
      out << oriented.size() << " ";
      for (std::vector<long>::const_iterator it = oriented.begin(); it != oriented.end(); ++it)
        out << *it << " ";
    }
    out << "\n";
  }

  out << "\nCELL_TYPES " << cells->N() << "\n";
  for (long c = 0; c < cells->N(); ++c)
    out << "42\n";
  if (labels) {
    out << "\nCELL_DATA " << cells->N() << "\nSCALARS label long 1\nLOOKUP_TABLE default\n";
    for (long c = 0; c < cells->N(); ++c)
      out << (*labels)[c] << "\n";
  }

}

#endif
