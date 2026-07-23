#ifndef PDMT_MED_WRITER_HPP
#define PDMT_MED_WRITER_HPP

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

void writePolyMed(std::string const * fineName, KNM < double > * nodesPoly, KN < KN < long >> * CellsPoly, KN < KN < long >> * EdgesPoly, KN < long > * LabelsPoly)
{
      //  get nodes of the mesh  //
      int TotalNodes = nodesPoly -> N();

     //double medNodeCoords[TotalNodes * 2];
      double* medNodeCoords = new double[TotalNodes * 2];
      for (int i = 0; i < TotalNodes; i++) {
        medNodeCoords[i * 2] = ( * nodesPoly)(i, 0);
        medNodeCoords[i * 2 + 1] = ( * nodesPoly)(i, 1);
      }

      //  get cells | cells + edges //
      int TotalCells = CellsPoly -> N();
      int TotalCellConnectivity = 0;
      int TotalElements = TotalCells;

      for (int i = 0; i < TotalCells; i++)
        TotalCellConnectivity += ( * CellsPoly)(i).N() + 1;

      if (EdgesPoly) {
        int TotalEdges = EdgesPoly -> N();

        for (int i = 0; i < TotalEdges; i++)
          TotalCellConnectivity += ( * EdgesPoly)(i).N() + 1;

        TotalElements += TotalEdges;
      }

      //mcIdType medCellConn[TotalCellConnectivity - TotalElements];
      mcIdType *medCellConn = new mcIdType[TotalCellConnectivity - TotalElements];

      int count = 0;
      for (int i = 0; i < CellsPoly -> N(); i++) {
        for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
          medCellConn[count] = ( * CellsPoly)(i)(j);
          count++;
        }
      }

      if (EdgesPoly) {
        for (int i = 0; i < EdgesPoly -> N(); i++) {
          for (int j = 0; j < ( * EdgesPoly)(i).N(); j++) {
            medCellConn[count] = ( * EdgesPoly)(i)(j);
            count++;
          }
        }
      }

      // --- fill med meshes ---	//
      MEDCouplingUMesh * medMesh2d = MEDCouplingUMesh::New();
      MEDCouplingUMesh * medMesh1d = MEDCouplingUMesh::New();

      medMesh2d -> setMeshDimension(2);
      medMesh2d -> allocateCells(CellsPoly -> N());
      medMesh2d -> setName("PolyMesh");

      count = 0;
      for (int i = 0; i < CellsPoly -> N(); i++) {
        medMesh2d -> insertNextCell(INTERP_KERNEL::NORM_POLYGON, ( * CellsPoly)(i).N(), medCellConn + count);
        count += ( * CellsPoly)(i).N();
      }

      medMesh2d -> finishInsertingCells();

      if (EdgesPoly) {
        medMesh1d -> setMeshDimension(1);
        medMesh1d -> allocateCells(EdgesPoly -> N());
        medMesh1d -> setName("PolyMesh");

        for (int i = 0; i < EdgesPoly -> N(); i++) {
          medMesh1d -> insertNextCell(INTERP_KERNEL::NORM_SEG2, 2, medCellConn + count);
          count += ( * EdgesPoly)(i).N();
        }
      }

      DataArrayDouble * myCoords = DataArrayDouble::New();
      myCoords -> alloc(TotalNodes, 2);
      myCoords -> setInfoOnComponent(0, "x");
      myCoords -> setInfoOnComponent(1, "y");
      std::copy(medNodeCoords, medNodeCoords + (TotalNodes * 2), myCoords -> getPointer());

      medMesh2d -> setCoords(myCoords);

      if (EdgesPoly)
        medMesh1d -> setCoords(myCoords);

      myCoords -> decrRef();

      if (!LabelsPoly) {
        std::vector < const MEDCouplingUMesh * > finalMesh;
        finalMesh.push_back(medMesh2d);
        if (EdgesPoly)
          finalMesh.push_back(medMesh1d);
        WriteUMeshes( * fineName, finalMesh, true);
      }

      if (LabelsPoly) {

        MCAuto < MEDFileUMesh > finalMeshWithLabel = MEDFileUMesh::New();

        finalMeshWithLabel -> setMeshAtLevel(0, medMesh2d);
        finalMeshWithLabel -> setMeshAtLevel(-1, medMesh1d);

        MCAuto < DataArrayIdType > fam2d = DataArrayIdType::New();
        MCAuto < DataArrayIdType > fam1d = DataArrayIdType::New();

        fam2d -> alloc(CellsPoly -> N(), 1);
        fam1d -> alloc(EdgesPoly -> N(), 1);

        mcIdType elemsFams[CellsPoly -> N() + EdgesPoly -> N()];

        std::set < int > poly2DUniqueLabels;
        std::set < int > poly1DUniqueLabels;

        for (int i = 0; i < CellsPoly -> N(); i++) {
          poly2DUniqueLabels.insert(( * LabelsPoly)(i));
          elemsFams[i] = ( * LabelsPoly)(i) + 1000; // Adding 1000 because med does not like tag zero
        }

        for (int i = CellsPoly -> N(); i < EdgesPoly -> N() + CellsPoly -> N(); i++) {
          poly1DUniqueLabels.insert(( * LabelsPoly)(i));
          elemsFams[i] = ( * LabelsPoly)(i);
        }

#ifdef DEBUG
        // Iterate through all the elements in a set and display the value.
        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it)
          std::cout << " Volume tag_i " << * it << endl;

        for (std::set < int > ::iterator it = poly1DUniqueLabels.begin(); it != poly1DUniqueLabels.end(); ++it)
          std::cout << " Surface tag_i " << * it << endl;
#endif
        std::copy(elemsFams, elemsFams + CellsPoly -> N(), fam2d -> getPointer());
        std::copy(elemsFams + CellsPoly -> N(), elemsFams + int(EdgesPoly -> N() + CellsPoly -> N()), fam1d -> getPointer());

        finalMeshWithLabel -> setFamilyFieldArr(-1, fam1d);
        finalMeshWithLabel -> setFamilyFieldArr(0, fam2d);

        std::map < std::string, std::vector < std::string >> theGroups;
        std::map < std::string, mcIdType > theFamilies;

        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it) {
          theFamilies["cell_family_" + to_string( * it) + ""] = * it + 1000;
          theGroups["cell_group_" + to_string( * it) + ""].push_back("cell_family_" + to_string( * it) + "");
        }

        for (std::set < int > ::iterator it = poly1DUniqueLabels.begin(); it != poly1DUniqueLabels.end(); ++it) {
          theFamilies["boundary_family_" + to_string( * it) + ""] = * it;
          theGroups["boundary_group_" + to_string( * it) + ""].push_back("boundary_family_" + to_string( * it) + "");
        }

    /*
            theFamilies["cells" ]=0;
            theFamilies["border1"]=1;
            theFamilies["border2"]=2;
            theFamilies["border3"]=3;
            theFamilies["border4"]=4;

            theGroups["Face_group"].push_back("cells");
            theGroups["boundary1"].push_back("border1");
            theGroups["boundary2"].push_back("border2");
            theGroups["boundary3"].push_back("border3");
            theGroups["boundary4"].push_back("border4");
    */

        // write mesh information useful for setting boundary conditions //

        ofstream infoWrite;
        infoWrite.open(""+* fineName+"_INFO.txt");

        infoWrite << "-------------------------------------------- \n"
                  << " info on mesh file " << * fineName << "\n"
                  << "-------------------------------------------- \n\n"
                  << "Information on nodes: \n"
                  << "  # nodes " << TotalNodes << "\n"
                  << "  List of families , tags , groups  \n\n"
                  << "Information on cells: \n"
                  << "  # polygons " << CellsPoly -> N() << "\n"
                  << "  List of families , tags , groups  \n\n";
        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it) {
          infoWrite << "    'cell_family_" + to_string( * it) + "'  has tag '" <<   * it + 1000 << "' belongs to group 'cell_group_" + to_string( * it) + "' " <<  endl;
        }

        infoWrite << "\n\n"
                  << "Information on boundary: \n"
                  << "  # edges " << EdgesPoly -> N() << "\n"
                  << "  List of families , tags , groups  \n\n";
        for (std::set < int > ::iterator it = poly1DUniqueLabels.begin(); it != poly1DUniqueLabels.end(); ++it) {
          infoWrite << "    'boundary_family_" + to_string( * it) + "'  has tag '" << * it << "' belongs to group 'boundary_group_" + to_string( * it) + "' " <<endl;
        }

        infoWrite.close();


        finalMeshWithLabel -> setFamilyInfo(theFamilies);
        finalMeshWithLabel -> setGroupInfo(theGroups);

        finalMeshWithLabel -> write( * fineName, 2); // med

      }
      medMesh1d -> decrRef();
      medMesh2d -> decrRef();
      delete[]     medNodeCoords;
      delete[]     medCellConn;
}

namespace PdmtMedWriter {

inline mcIdType cellFamily(long label) {
  // Keep lower-dimensional boundary labels unchanged, as the historical 2D
  // writer does, while putting cell families in a disjoint range.
  return static_cast<mcIdType>(label + 1000);
}

inline mcIdType boundaryFamily(long label) {
  // MED family zero means "no family" and therefore cannot form a group.
  return static_cast<mcIdType>(label == 0 ? 2000 : label);
}

inline void setCoordinates(MEDCouplingUMesh *mesh, KNM<double> *nodes,
                           int spaceDimension) {
  MCAuto<DataArrayDouble> coordinates = DataArrayDouble::New();
  coordinates->alloc(nodes->N(), spaceDimension);
  coordinates->setInfoOnComponent(0, "x");
  coordinates->setInfoOnComponent(1, "y");
  if (spaceDimension == 3)
    coordinates->setInfoOnComponent(2, "z");
  double *out = coordinates->getPointer();
  for (long node = 0; node < nodes->N(); ++node)
    for (int component = 0; component < spaceDimension; ++component)
      out[node * spaceDimension + component] = (*nodes)(node, static_cast<long>(component));
  mesh->setCoords(coordinates);
}

inline void addCellFamilies(MEDFileUMesh *fileMesh, int level,
                            KN<long> *labels, long count,
                            const std::string &kind,
                            std::map<std::string, mcIdType> &families,
                            std::map<std::string, std::vector<std::string> > &groups,
                            bool topLevel) {
  if (!labels || count == 0)
    return;
  MCAuto<DataArrayIdType> field = DataArrayIdType::New();
  field->alloc(count, 1);
  std::set<long> uniqueLabels;
  for (long i = 0; i < count; ++i) {
    const long label = (*labels)[i];
    field->setIJ(i, 0, topLevel ? cellFamily(label) : boundaryFamily(label));
    uniqueLabels.insert(label);
  }
  fileMesh->setFamilyFieldArr(level, field);
  for (std::set<long>::const_iterator label = uniqueLabels.begin();
       label != uniqueLabels.end(); ++label) {
    const mcIdType familyId = topLevel ? cellFamily(*label) : boundaryFamily(*label);
    const std::string familyName = kind + "_family_" + to_string(*label);
    const std::string groupName = kind + "_group_" + to_string(*label);
    families[familyName] = familyId;
    groups[groupName].push_back(familyName);
  }
}

} // namespace PdmtMedWriter

// Write a polygonal surface embedded in 3D.  Boundary edges are reconstructed
// from polygon incidence, so open surfaces retain a valid MED level -1 mesh.
inline void writePolySurfaceMed(std::string const *fileName,
                                KNM<double> *nodes,
                                KN<KN<long> > *polygons,
                                KN<long> *labels) {
  using namespace PdmtMedWriter;
  if (!nodes || !polygons || nodes->M() < 3)
    ExecError("PdmtPolyMeshWrite: invalid 3S polygon arrays for MED output");
  if (labels && labels->N() != polygons->N())
    ExecError("PdmtPolyMeshWrite: labels must have one value per surface polygon");

  MCAuto<MEDCouplingUMesh> surface = MEDCouplingUMesh::New();
  surface->setMeshDimension(2);
  surface->setName("PolySurfaceMesh");
  surface->allocateCells(polygons->N());
  for (long polygon = 0; polygon < polygons->N(); ++polygon) {
    if ((*polygons)[polygon].N() < 3)
      ExecError("PdmtPolyMeshWrite: a surface polygon must contain at least three nodes");
    std::vector<mcIdType> connectivity((*polygons)[polygon].N());
    for (long vertex = 0; vertex < (*polygons)[polygon].N(); ++vertex) {
      const long node = (*polygons)[polygon][vertex];
      if (node < 0 || node >= nodes->N())
        ExecError("PdmtPolyMeshWrite: a surface polygon references an out-of-range node");
      connectivity[vertex] = node;
    }
    surface->insertNextCell(INTERP_KERNEL::NORM_POLYGON,
                            connectivity.size(), &connectivity[0]);
  }
  surface->finishInsertingCells();
  setCoordinates(surface, nodes, 3);
  surface->checkConsistencyLight();

  typedef std::pair<long, long> Edge;
  std::map<Edge, long> edgeUse;
  for (long polygon = 0; polygon < polygons->N(); ++polygon)
    for (long vertex = 0; vertex < (*polygons)[polygon].N(); ++vertex) {
      long a = (*polygons)[polygon][vertex];
      long b = (*polygons)[polygon][(vertex + 1) % (*polygons)[polygon].N()];
      if (a > b)
        std::swap(a, b);
      ++edgeUse[Edge(a, b)];
    }
  std::vector<Edge> boundaryEdges;
  for (std::map<Edge, long>::const_iterator edge = edgeUse.begin();
       edge != edgeUse.end(); ++edge)
    if (edge->second == 1)
      boundaryEdges.push_back(edge->first);

  MCAuto<MEDCouplingUMesh> boundary;
  if (!boundaryEdges.empty()) {
    boundary = MEDCouplingUMesh::New();
    boundary->setMeshDimension(1);
    boundary->setName("PolySurfaceMesh");
    boundary->allocateCells(boundaryEdges.size());
    for (std::vector<Edge>::const_iterator edge = boundaryEdges.begin();
         edge != boundaryEdges.end(); ++edge) {
      mcIdType connectivity[2] = {edge->first, edge->second};
      boundary->insertNextCell(INTERP_KERNEL::NORM_SEG2, 2, connectivity);
    }
    boundary->finishInsertingCells();
    boundary->setCoords(surface->getCoords());
    boundary->checkConsistencyLight();
  }

  MCAuto<MEDFileUMesh> output = MEDFileUMesh::New();
  output->setMeshAtLevel(0, surface);
  if (boundary)
    output->setMeshAtLevel(-1, boundary);

  std::map<std::string, mcIdType> families;
  std::map<std::string, std::vector<std::string> > groups;
  addCellFamilies(output, 0, labels, polygons->N(), "polygon", families, groups, true);
  if (!families.empty()) {
    output->setFamilyInfo(families);
    output->setGroupInfo(groups);
  }
  output->write(*fileName, 2);
}

// Write explicit polyhedra in MED's NORM_POLYHED representation.  Each cell
// contains its oriented polygonal faces separated by -1.  Exterior faces are
// also written at level -1 so MED/SALOME retains the boundary labels.
inline void writePolyMed3D(std::string const *fileName,
                           KNM<double> *nodes,
                           KN<KN<long> > *cells,
                           KN<KN<long> > *faces,
                           KN<long> *labels,
                           KN<long> *faceLabels) {
  using namespace PdmtMedWriter;
  using namespace Pdmt3DWriter;
  if (!nodes || !cells || !faces || nodes->M() < 3)
    ExecError("PdmtPolyMeshWrite: invalid 3D polyhedron arrays for MED output");
  validateTopology(*nodes, *cells, *faces);
  if (labels && labels->N() != cells->N())
    ExecError("PdmtPolyMeshWrite: labels must have one value per polyhedron");
  if (faceLabels && faceLabels->N() != faces->N())
    ExecError("PdmtPolyMeshWrite: faceLabels must have one value per face");

  MCAuto<MEDCouplingUMesh> volume = MEDCouplingUMesh::New();
  volume->setMeshDimension(3);
  volume->setName("PolyMesh");
  volume->allocateCells(cells->N());
  std::vector<long> oriented;
  for (long cell = 0; cell < cells->N(); ++cell) {
    std::vector<mcIdType> connectivity;
    for (long localFace = 0; localFace < (*cells)[cell].N(); ++localFace) {
      if (localFace)
        connectivity.push_back(-1);
      orientedFace((*cells)[cell][localFace], *faces, oriented);
      connectivity.insert(connectivity.end(), oriented.begin(), oriented.end());
    }
    volume->insertNextCell(INTERP_KERNEL::NORM_POLYHED,
                           connectivity.size(), &connectivity[0]);
  }
  volume->finishInsertingCells();
  setCoordinates(volume, nodes, 3);
  volume->checkConsistencyLight();

  std::vector<long> faceUse(faces->N(), 0);
  for (long cell = 0; cell < cells->N(); ++cell)
    for (long localFace = 0; localFace < (*cells)[cell].N(); ++localFace)
      ++faceUse[faceIndex((*cells)[cell][localFace])];
  std::vector<long> boundaryFaceIds;
  for (long face = 0; face < faces->N(); ++face)
    if (faceUse[face] == 1)
      boundaryFaceIds.push_back(face);

  MCAuto<MEDCouplingUMesh> boundary;
  if (!boundaryFaceIds.empty()) {
    boundary = MEDCouplingUMesh::New();
    boundary->setMeshDimension(2);
    boundary->setName("PolyMesh");
    boundary->allocateCells(boundaryFaceIds.size());
    for (std::vector<long>::const_iterator faceId = boundaryFaceIds.begin();
         faceId != boundaryFaceIds.end(); ++faceId) {
      std::vector<mcIdType> connectivity((*faces)[*faceId].N());
      for (long vertex = 0; vertex < (*faces)[*faceId].N(); ++vertex)
        connectivity[vertex] = (*faces)[*faceId][vertex];
      boundary->insertNextCell(INTERP_KERNEL::NORM_POLYGON,
                               connectivity.size(), &connectivity[0]);
    }
    boundary->finishInsertingCells();
    boundary->setCoords(volume->getCoords());
    boundary->checkConsistencyLight();
  }

  MCAuto<MEDFileUMesh> output = MEDFileUMesh::New();
  output->setMeshAtLevel(0, volume);
  if (boundary)
    output->setMeshAtLevel(-1, boundary);

  std::map<std::string, mcIdType> families;
  std::map<std::string, std::vector<std::string> > groups;
  addCellFamilies(output, 0, labels, cells->N(), "cell", families, groups, true);
  if (faceLabels && !boundaryFaceIds.empty()) {
    KN<long> boundaryLabels(boundaryFaceIds.size());
    for (long face = 0; face < static_cast<long>(boundaryFaceIds.size()); ++face)
      boundaryLabels[face] = (*faceLabels)[boundaryFaceIds[face]];
    addCellFamilies(output, -1, &boundaryLabels, boundaryFaceIds.size(),
                    "boundary", families, groups, false);
  }
  if (!families.empty()) {
    output->setFamilyInfo(families);
    output->setGroupInfo(groups);
  }
  output->write(*fileName, 2);
}

#endif
