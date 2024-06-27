#include <iomanip>

// Write typ mesh format
void writePolyTyp(std::string const * fineName, KNM < double > * nodesPoly, KN < KN < long >> * CellsPoly)
{

  bool highPrecision = true;

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

    if (highPrecision) {
        polyWrite << std::fixed << std::setprecision(10);
    }

    for (int i = 0; i < TotalNodes; i++)
      polyWrite << "    " << ( * nodesPoly)(i, 0) << "    " << ( * nodesPoly)(i, 1) << "\n";
  }

    if (highPrecision) {
        // Reset the stream to default precision (optional)
        polyWrite << std::defaultfloat;
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
      polyWrite << std::setw(8) << ( * CellsPoly)(i).N();

      for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
        polyWrite << std::setw(8) << ( * CellsPoly)(i)(j);
      }
      polyWrite << "\n";
    }
  }

}
