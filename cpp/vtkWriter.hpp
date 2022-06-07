
void writePolyVtu(std::string const * fineName, KNM < double > * nodesPoly, KN < KN < long >> * CellsPoly, KN < KN < long >> * EdgesPoly, KN < long > * LabelsPoly)
{
  ofstream polyWrite;
  polyWrite.open( * fineName);

  //------------ Write header ----------------//

  polyWrite << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n" <<
    "  <UnstructuredGrid>\n";

  // calculate NumberOfPoints and  NumberOfCells
  int NumberOfPoints, NumberOfCells; {
    NumberOfPoints = nodesPoly -> N();

    NumberOfCells = CellsPoly -> N();
    if (EdgesPoly)
      NumberOfCells += EdgesPoly -> N();
  }

  if (verbosity) {
    std::cout << "------------------------------------------------------ " << std::endl;
    std::cout << "PMDT NODES IN POLYMESH " << NumberOfPoints << std::endl;
    std::cout << "PMDT CELLS IN POLYMESH " << NumberOfCells << std::endl;
    std::cout << "------------------------------------------------------ " << std::endl;
  }

  //------------ Write Piece info ----------------//
  {
    polyWrite << "    <Piece NumberOfPoints=\"" << NumberOfPoints << "\" NumberOfCells=\"" << NumberOfCells << "\">\n";
  }

  //------------ Write CellData for labels ----------------//
  {
    if (LabelsPoly) {
      polyWrite << "      <CellData Scalars=\"label\">\n" <<
        "        <DataArray type=\"Float32\" Name=\"label\" format=\"ascii\">\n\t";

      for (int i = 0; i < NumberOfCells; i++)
        polyWrite << ( * LabelsPoly)(i) << "\t";

      polyWrite << "\n        </DataArray>\n" <<
        "      </CellData>\n";
    }
  }

  //------------ Write PointsData ----------------//
  {
    polyWrite << "      <Points>\n" <<
      "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n\t";

    for (int i = 0; i < NumberOfPoints; i++)
      polyWrite << ( * nodesPoly)(i, 0) << "\t" << ( * nodesPoly)(i, 1) << "\t 0\n\t";

    polyWrite << "\n        </DataArray>\n" <<
      "      </Points>\n";
  }

  //------------ Write CellData for labels ----------------//
  {
    polyWrite << "      <Cells>\n" <<
      "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n\t";

    for (int i = 0; i < CellsPoly -> N(); i++) {
      for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
        polyWrite << ( * CellsPoly)(i)(j) << " ";
      }
      polyWrite << "\n\t";
    }

    if (EdgesPoly) {
      for (int i = 0; i < EdgesPoly -> N(); i++) {
        for (int j = 0; j < ( * EdgesPoly)(i).N(); j++) {
          polyWrite << ( * EdgesPoly)(i)(j) << " ";
        }
        polyWrite << "\n\t";
      }
    }

    polyWrite << "\n        </DataArray>\n";
    polyWrite << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n\t";

    int offset = 0;
    for (int i = 0; i < CellsPoly -> N(); i++) {
      offset += ( * CellsPoly)(i).N();
      polyWrite << offset << "\t";
    }

    if (EdgesPoly) {
      for (int i = 0; i < EdgesPoly -> N(); i++) {
        offset += ( * EdgesPoly)(i).N();
        polyWrite << offset << "\t";

      }
    }

    polyWrite << "\n        </DataArray>\n";
    polyWrite << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t";

    for (int i = 0; i < CellsPoly -> N(); i++) {
      offset += ( * CellsPoly)(i).N();
      polyWrite << "7\t";
    }

    if (EdgesPoly) {
      for (int i = 0; i < EdgesPoly -> N(); i++) {
        offset += ( * EdgesPoly)(i).N();
        polyWrite << "3\t";

      }
    }

    polyWrite << "\n        </DataArray>\n" <<
      "      </Cells>\n";
  }

  polyWrite << "    </Piece>\n" <<
    "  </UnstructuredGrid>\n" <<
    "</VTKFile>\n";
}


void writePolyVtk(std::string const * fineName, KNM < double > * nodesPoly, KN < KN < long >> * CellsPoly, KN < KN < long >> * EdgesPoly, KN < long > * LabelsPoly)
{

  ofstream polyWrite;
  polyWrite.open( * fineName);

  //------------ Write header ----------------//

  polyWrite << "# vtk DataFile Version 2.0\n" <<
    "Unstructured Grid PDMT\n" <<
    "ASCII\n" <<
    "DATASET UNSTRUCTURED_GRID\n\n";

  //------------ Write nodes ----------------//

  if (verbosity) {
    std::cout << "------------------------------------------------------ " << std::endl;
    std::cout << "PMDT NODES IN POLYMESH " << nodesPoly -> N() << std::endl;
    std::cout << "------------------------------------------------------ " << std::endl;
  }

  {
    int TotalNodes = nodesPoly -> N();

    polyWrite << "POINTS " << TotalNodes << " float\n";
    for (int i = 0; i < TotalNodes; i++)
      polyWrite << ( * nodesPoly)(i, 0) << "\t" << ( * nodesPoly)(i, 1) << "\t 0\n";
  }

  //------------ Write cells ----------------//
  if (verbosity) {
    std::cout << "------------------------------------------------------ " << std::endl;
    std::cout << "PMDT CELS IN POLYMESH " << CellsPoly -> N() << std::endl;
    //std::cout << "PMDT CELS IN POLYMESH " << (*CellsPoly)(0).N()  << std::endl;
    //std::cout << "PMDT CELS IN POLYMESH " << (*CellsPoly)(0)(0)  << std::endl;
    std::cout << "------------------------------------------------------ " << std::endl;
  }

  {
    int TotalCells = CellsPoly -> N();
    int TotalCellConnectivity = 0;

    int TotalVTKConnectionList = TotalCells;

    for (int i = 0; i < TotalCells; i++)
      TotalCellConnectivity += ( * CellsPoly)(i).N() + 1;

    if (EdgesPoly) {
      int TotalEdges = EdgesPoly -> N();
      for (int i = 0; i < TotalEdges; i++)
        TotalCellConnectivity += ( * EdgesPoly)(i).N() + 1;

      TotalVTKConnectionList += TotalEdges;
    }

    polyWrite << "\n";
    polyWrite << "CELLS " << TotalVTKConnectionList << "\t" << TotalCellConnectivity;
    polyWrite << "\n";

    for (int i = 0; i < TotalCells; i++) {
      polyWrite << ( * CellsPoly)(i).N() << " ";
      for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
        polyWrite << ( * CellsPoly)(i)(j) << " ";
      }
      polyWrite << "\n";
    }

    if (EdgesPoly) {
      for (int i = 0; i < EdgesPoly -> N(); i++) {
        polyWrite << ( * EdgesPoly)(i).N() << " ";
        for (int j = 0; j < ( * EdgesPoly)(i).N(); j++) {
          polyWrite << ( * EdgesPoly)(i)(j) << " ";
        }
        polyWrite << "\n";
      }
    }

  }

  //------------ Write cell types -----------------------//
  {
    int TotalCells = CellsPoly -> N();

    if (EdgesPoly)
      TotalCells += EdgesPoly -> N();

    polyWrite << "\n";
    polyWrite << "\n";
    polyWrite << "CELL_TYPES " << TotalCells << " ";
    polyWrite << "\n";
    for (int i = 0; i < CellsPoly -> N(); i++)
      polyWrite << "7" << "\n";

    if (EdgesPoly)
      for (int i = 0; i < EdgesPoly -> N(); i++)
        polyWrite << "3" << "\n";
  }

  //------------ Write cell labels -----------------------//
  {
    int TotalCells = CellsPoly -> N();

    if (EdgesPoly)
      TotalCells += EdgesPoly -> N();

    polyWrite << "\n";
    polyWrite << "\n";
    polyWrite << "CELL_DATA " << TotalCells << " \n";
    polyWrite << "SCALARS label float 1 \n";
    polyWrite << "LOOKUP_TABLE CellColors ";
    polyWrite << "\n";

    if (!LabelsPoly)
      for (int i = 0; i < CellsPoly -> N(); i++)
        polyWrite << "0" << "\n";
    else
      for (int i = 0; i < CellsPoly -> N(); i++)
        polyWrite << ( * LabelsPoly)(i) << "\n";

    if (EdgesPoly && !LabelsPoly)
      for (int i = 0; i < EdgesPoly -> N(); i++)
        polyWrite << "1" << "\n";

    if (EdgesPoly && LabelsPoly)
      for (int i = CellsPoly -> N(); i < EdgesPoly -> N() + CellsPoly -> N(); i++)
        polyWrite << ( * LabelsPoly)(i) << "\n";
  }

}
/**/
