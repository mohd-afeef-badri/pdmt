// Write typ mesh format
void writePolyTyp(std::string const * fineName, KNM < double > * nodesPoly, KN < KN < long >> * CellsPoly)
{
  ofstream polyWrite;
  polyWrite.open( * fineName);

  //------------ Write nodes ----------------//

  if (verbosity) {
    std::cout << "------------------------------------------------------ " << std::endl;
    std::cout << "PMDT NODES IN POLYMESH " << nodesPoly -> N() << std::endl;
    std::cout << "------------------------------------------------------ " << std::endl;
  }

  {
    int TotalNodes = nodesPoly -> N();

    polyWrite << " Vertices\n          " << TotalNodes << "\n";
    for (int i = 0; i < TotalNodes; i++)
      polyWrite << "\t" << ( * nodesPoly)(i, 0) << "\t" << ( * nodesPoly)(i, 1) << "\n";
  }

  //------------ Write cells ----------------//
  if (verbosity) {
    std::cout << "------------------------------------------------------ " << std::endl;
    std::cout << "PMDT CELS IN POLYMESH " << CellsPoly -> N() << std::endl;
    std::cout << "------------------------------------------------------ " << std::endl;
  }

  {
    int TotalCells = CellsPoly -> N();
    int TotalCellConnectivity = 0;

    int TotalVTKConnectionList = TotalCells;

    for (int i = 0; i < TotalCells; i++)
      TotalCellConnectivity += ( * CellsPoly)(i).N() + 1;

    polyWrite << " cells\n          " << TotalCells << "\n";

    for (int i = 0; i < TotalCells; i++) {
      polyWrite << ( * CellsPoly)(i).N() << " ";
      for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
        polyWrite << ( * CellsPoly)(i)(j) << " ";
      }
      polyWrite << "\n";
    }
  }

}
