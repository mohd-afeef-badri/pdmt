/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 15/01/2022
     Comment  : The program finds convex hull of a set of  points.
                Original  program  found  online  was  adapted for
                PDMT.  For explanation of the orientation()
                www.geeksforgeeks.org/orientation-3-ordered-points

     -------------------------------------------------------------------

     PDMT a parallel  dual meshing tool uses   finite  element framework
     to convert a triangular / tetrahedral mesh into a  polyhedral  mesh.
     PDMT is distributed  in  the  hope that it  will be useful, HOWEVER
     WITHOUT ANY WARRANTY; or without  even  implied warranty of FITNESS
     FOR A PARTICULAR PURPOSE.

*******************************************************************************/

#include <iostream>
#include <stack>
#include <stdlib.h>

using namespace std;

struct Point2D
{
  double x , y;   // x and y coordinates
  long int nn;    // node number in the mesh
};

// A global point needed for  sorting points with reference
// to  the first  point Used in compare function of qsort()
Point2D p0;

// A utility function to find next to top in a stack
Point2D nextToTop(stack<Point2D> &S)
{
  Point2D p = S.top();
  S.pop();
  Point2D res = S.top();
  S.push(p);

  return res;
}

// A utility function to swap two points
void swap(Point2D &p1, Point2D &p2)
{
  Point2D temp = p1;
  p1 = p2;
  p2 = temp;
}

// A utility function to return square of distance
// between p1 and p2
double distSq(Point2D p1, Point2D p2)
{
  return (p1.x - p2.x)*(p1.x - p2.x) +
         (p1.y - p2.y)*(p1.y - p2.y) ;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point2D p, Point2D q, Point2D r)
{
  double val = (q.y - p.y) * (r.x - q.x) -
               (q.x - p.x) * (r.y - q.y);

  if (val == 0) return 0;   // collinear
  return (val > 0)? 1: 2;   // clock or counterclock wise
}

// A function used by library function qsort() to sort an
// array of points with respect to the first point
int compare(const void *vp1, const void *vp2)
{
  Point2D *p1 = (Point2D *)vp1;
  Point2D *p2 = (Point2D *)vp2;

  // Find orientation
  int o = orientation(p0, *p1, *p2);

  if (o == 0)
    return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

  return (o == 2)? -1: 1;
}

// Prints convex hull of a set of n points.
void convexHull(Point2D points[], int n)
{
  // Find the bottommost point
  double ymin = points[0].y; int min = 0;
  for (int i = 1; i < n; i++)
  {
    double y = points[i].y;

    // Pick the bottom-most or chose the left
    // most point in case of tie
    if ((y < ymin) || (ymin == y && points[i].x < points[min].x))
       ymin = points[i].y, min = i;
  }

  // Place the bottom-most point at first position
  swap(points[0], points[min]);

  // Sort n-1 points with respect to the first point.
  // A point p1 comes before p2 in sorted output if p2
  // has larger polar angle (in counter-clockwise
  // direction) than p1

  p0 = points[0];
  qsort(&points[1], n-1, sizeof(Point2D), compare);

  // If two or more points make same angle with p0,
  // Remove all but the one that is farthest from p0
  // Remember that, in above sorting, our criteria was
  // to keep the farthest point at the end when more than
  // one points have same angle.

  int m = 1; // Initialize size of modified array
  for (int i=1; i<n; i++)
  {
    // Keep removing i while angle of i and i+1 is same
    // with respect to p0
    while (i < n-1 && orientation(p0, points[i], points[i+1]) == 0)
      i++;

    points[m] = points[i];
    m++;                      // Update size of modified array
  }

/*
   int m = 1; // Initialize size of modified array
   m = n;
*/

  // If modified array of points has less than 3 points,
  // convex hull is not possible
  if (m < 3) return;

  // Create an empty stack and push first three points
  // to it.
  stack<Point2D> S;
  S.push(points[0]);
  S.push(points[1]);
  S.push(points[2]);

  // Process remaining n-3 points
  for (int i = 3; i < m; i++)
  {
    // Keep removing top while the angle formed by
    // points next-to-top, top, and points[i] makes
    // a non-left turn
    while (S.size()>1 && orientation(nextToTop(S), S.top(), points[i]) != 2)
       S.pop();
    S.push(points[i]);
  }

  // Now stack has the output points
  // print contents of stack
  while (!S.empty())
  {
    Point2D p = S.top();
#ifdef DEBUG
    cout << "(" << p.x << ", " << p.y <<")" << endl;
#endif
    S.pop();
  }
}


int PdmtConvexHull(KN<double> *const & px, KN<double> *const & py, KN<long> *const & nI)
{
  long int nn = px->N();

  Point2D points[nn];

  for(long int i=0; i<nn; i++)
  {
     points[i].x  = *(px[0]+i);
     points[i].y  = *(py[0]+i);
     points[i].nn = *(nI[0]+i);
  }

  int n = sizeof(points)/sizeof(points[0]);

#ifdef DEBUG
  cout << "" << endl;
  cout << "--------------------------------------" << endl;
  cout << " PdmtConvexHull BEFORE C-HULL FUNCTION" << endl;
  cout << "--------------------------------------" << endl;
  cout << "x" << "\t" << "y" << "\t" << "n" << "\n\n" << endl;

  for(long int i=0; i<nn; i++)
    cout << points[i].x << "\t" << points[i].y << "\t" << points[i].nn << "\n" << endl;

  cout << "--------------------------------------\n" << endl;
#endif

  convexHull(points, n);

#ifdef DEBUG
  cout << "" << endl;
  cout << "--------------------------------------" << endl;
  cout << " PdmtConvexHull AFTER C-HULL FUNCTION " << endl;
  cout << "--------------------------------------" << endl;
  cout << "x" << "\t" << "y" << "\t" << "n" << "\n\n" << endl;

  for(long int i=0; i<nn; i++)
    cout << points[i].x << "\t" << points[i].y << "\t" << points[i].nn << "\n" << endl;

  cout << "--------------------------------------\n" << endl;
#endif

  for(long int i=0; i<nn; i++)
      *(nI[0]+i) = points[i].nn ;

  return 0;
}

int PdmtTest()
{
  Point2D points[6];

/*
0.761755	0.761755	6
0.332496	1.02236	5
0.638909	0.995436	1
0.5	        0.866025	14
0.683013	0.683013	23
0.25	        0.933013	22

  points[0].x = 0;     points[0].y = 0;      points[0].nn = 0;
  points[1].x = 0.25;  points[1].y = 0.75;   points[1].nn = 1;
  points[2].x = 0.5;   points[2].y = 0.5;    points[2].nn = 2;
  points[3].x = 0;     points[3].y = 0.75;   points[3].nn = 3;
  points[4].x = 0.5;   points[4].y = 0;      points[4].nn = 4;
*/

  points[0].x = 0.761755;  points[0].y = 0.761755;   points[0].nn = 6 ;
  points[1].x = 0.332496;  points[1].y = 1.022360;   points[1].nn = 5 ;
  points[2].x = 0.638909;  points[2].y = 0.995436;   points[2].nn = 1 ;
  points[3].x = 0.500000;  points[3].y = 0.866025;   points[3].nn = 7 ;
  points[4].x = 0.683013;  points[4].y = 0.683013;   points[4].nn = 23;
  points[5].x = 0.250000;  points[5].y = 0.933013;   points[5].nn = 22;

  int n = sizeof(points)/sizeof(points[0]);

  cout << "" << endl;
  cout << "--------------------------------------" << endl;
  cout << "BEFORE C-HULL FUNCTION     " << endl;
  cout << "--------------------------------------" << endl;
  cout << "x" << "\t" << "y" << "\t" << "n" << "\n\n" <<
  points[0].x << "\t" << points[0].y << "\t" << points[0].nn << "\n" <<
  points[1].x << "\t" << points[1].y << "\t" << points[1].nn << "\n" <<
  points[2].x << "\t" << points[2].y << "\t" << points[2].nn << "\n" <<
  points[3].x << "\t" << points[3].y << "\t" << points[3].nn << "\n" <<
  points[4].x << "\t" << points[4].y << "\t" << points[4].nn << "\n" <<
  points[5].x << "\t" << points[5].y << "\t" << points[5].nn << "\n"   ;
  cout << "--------------------------------------" << endl;
  cout << "" << endl;

  convexHull(points, n);

  cout << "" << endl;
  cout << "--------------------------------------" << endl;
  cout << "AFTER C-HULL FUNCTION     " << endl;
  cout << "--------------------------------------" << endl;
  cout << "x" << "\t" << "y" << "\t" << "n" << "\n\n" <<
  points[0].x << "\t" << points[0].y << "\t" << points[0].nn << "\n" <<
  points[1].x << "\t" << points[1].y << "\t" << points[1].nn << "\n" <<
  points[2].x << "\t" << points[2].y << "\t" << points[2].nn << "\n" <<
  points[3].x << "\t" << points[3].y << "\t" << points[3].nn << "\n" <<
  points[4].x << "\t" << points[4].y << "\t" << points[4].nn << "\n" <<
  points[5].x << "\t" << points[5].y << "\t" << points[5].nn << "\n"   ;
  cout << "--------------------------------------" << endl;
  cout << "" << endl;

 return 0;
}
