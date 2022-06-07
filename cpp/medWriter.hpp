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
                  << "  List of familes , tags , groups  \n\n"
                  << "Information on cells: \n"
                  << "  # polygons " << CellsPoly -> N() << "\n"
                  << "  List of familes , tags , groups  \n\n";
        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it) {
          infoWrite << "    'cell_family_" + to_string( * it) + "'  has tag '" <<   * it + 1000 << "' belongs to group 'cell_group_" + to_string( * it) + "' " <<  endl;
        }

        infoWrite << "\n\n"
                  << "Information on boundary: \n"
                  << "  # edges " << EdgesPoly -> N() << "\n"
                  << "  List of familes , tags , groups  \n\n";
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
