
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


void parallelIO (string *&name, MPI_Comm *const &comm, int *tPvd, int *mpiSize, std::string *baseFilename)
{

  std::string base_filename (*name);
  std::string::size_type p (base_filename.find_last_of ('.'));
  std::string file_without_extension = base_filename.substr (0, p);
  std::string extension;

  if (p == std::string::npos)
    extension = "vtu";
  else
    extension = base_filename.substr (p + 1, std::string::npos);

  p = base_filename.find_last_of ("/\\");

  if (p == std::string::npos)
    base_filename = file_without_extension;
  else
    base_filename = file_without_extension.substr (p + 1, std::string::npos);

  int rank;
  int size;
  MPI_Comm_rank (comm ? *comm : MPI_COMM_WORLD, &rank);
  MPI_Comm_size (comm ? *comm : MPI_COMM_WORLD, &size);

  std::ostringstream str[3];
  str[2] << size;
  str[1] << std::setw (str[2].str ().length ()) << std::setfill ('0') << rank;

  ofstream pvd, pvtu;

  int T = 0;

  if (rank == 0)
    {
      ifstream input;
      input.open (file_without_extension + (size > 1 ? "_" + str[2].str () : "") + ".pvd");
      if (input.peek () != std::ifstream::traits_type::eof ())
        {
          std::string line;
          std::getline (input, line);
          std::getline (input, line);
          std::string delimiter = "\"";
          p = line.find (delimiter);
          line = line.substr (p + 1, std::string::npos);
          p = line.find (delimiter);
          T = std::stoi (line.substr (0, p)) + 1;
        }
    }
  MPI_Bcast (&T, 1, MPI_INT, 0, comm ? *comm : MPI_COMM_WORLD);

  str[0] << std::setw (4) << std::setfill ('0') << T;
  *name = file_without_extension + "_" + (size > 1 ? str[2].str () + "_" : "") + str[0].str ()
          + (size > 1 ? "_" + /*str[1].str()*/ to_string (rank) : "") + "." + extension;
  *tPvd = T;                     // Time for pvd file
  *mpiSize = size;               // MPI size
  *baseFilename = base_filename; // File name without extension
}


void PvtuWriter( string*& pffname , const int mpiSize, const int timePvd, const string basVtuFileName)
{

  std::string ProcZero;

  if (mpiSize > 1)
    ProcZero = "_0.vtu";                            // This will indicate if Proc is 0
  else
    ProcZero = ".vtu";

  std::string fullFileName(*pffname);
  std::size_t foundProcZero = (fullFileName).find(ProcZero);

  if (foundProcZero!=std::string::npos)           // Only proc 0 does the pvtu/pvd writing
  {

    std::string file_without_extension;

    if (mpiSize > 1)
      file_without_extension = fullFileName.substr(0, (foundProcZero - 6 - string(to_string(mpiSize)).length()));
    else
      file_without_extension = fullFileName.substr(0, (foundProcZero - 5));


    ofstream pvtu;

    std::string timeStamp = std::string(4 - string(to_string(timePvd)).length(), '0') + to_string(timePvd);
    pvtu.open(file_without_extension + (mpiSize > 1 ? "_" + to_string(mpiSize) : "") + "_" + timeStamp + ".pvtu");

    pvtu << "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <PUnstructuredGrid GhostLevel=\"0\">\n"
            "    <PCellData>\n"
            "      <PDataArray type=\"Int32\" Name=\"Label\"/>\n"
            "    </PCellData>\n"
            "    <PPoints>\n"
            "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Points\"/>\n"
            "    </PPoints>\n";

    for(int i = 0; i < mpiSize; ++i) {
      pvtu << "    <Piece Source=\"";
      pvtu << basVtuFileName << "_";
      if(mpiSize > 1)
       pvtu << to_string(mpiSize) + "_";
      pvtu << std::setw(4) << std::setfill('0') << timePvd;
      if(mpiSize > 1)
        pvtu << "_" << std::setfill('0') << i;
      pvtu << ".vtu\"/>\n";
    }

    pvtu << "  </PUnstructuredGrid>\n"
            "</VTKFile>\n";
  }
}
