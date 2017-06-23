#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <time.h>
#include <sstream>

// define the kernel
typedef double numberType;
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<numberType>  Kernel;

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Apollonius_graph_adaptation_traits_2.h>
#include <CGAL/Apollonius_graph_adaptation_policies_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

// define the Voronoi diagram adaptor
typedef CGAL::Apollonius_graph_traits_2<Kernel> Traits;
typedef CGAL::Apollonius_graph_2<Traits> Apo_graph;
typedef CGAL::Apollonius_graph_adaptation_traits_2<Apo_graph> AT;
typedef CGAL::Apollonius_graph_degeneracy_removal_policy_2<Apo_graph> AP;
typedef CGAL::Voronoi_diagram_2<Apo_graph,AT,AP> Diagram;

typedef CGAL::Apollonius_site_2<Kernel> Site;
typedef Site::Point_2 Point;
typedef Site::Weight Weight;
typedef Diagram::Face Face;
typedef Diagram::Halfedge Halfedge;
typedef Diagram::Ccb_halfedge_circulator Ccb_halfedge_circulator;
typedef Diagram::Face_iterator Face_iterator;



void insertSite(Apo_graph* ag, FILE* f, double x, double y, double weight) {
  Point p(x, y);
  Weight w(weight);
  Site s(p, w);
  ag->insert(s);
  if (f != NULL) {
    fprintf(f, "%.30f %.30f %.30f\n", p.x(), p.y(), w);
  }
}

// inserts sites outside of the visible area so that each hyperbola inside has
// finite endpoints
void insertBoundingSites(Apo_graph* ag, double maxX, double maxY, double maxWeight) {
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      // I don't use x and y in the header of the for loop directly because of
      // to avoid the comparison of unprecise double values.
      double x = -maxWeight - maxX + i * (2*maxWeight + 3*maxX);
      double y = -maxWeight - maxY + j * (2*maxWeight + 3*maxY);
      Point p(x, y);
      Weight w(0);
      Site s(p, w);
      ag->insert(s);
    }
  }
}

// generates the Apollonius graph of n random points with random weights in the
// rectangle width x height
// the generated data is saved to f
Apo_graph* generateApo_graphRandom(int n, double width, double height, std::string sitesPath) {
  Apo_graph* ag = new Apo_graph();
  FILE* sitesFile = fopen(sitesPath.c_str(), "w");

  double maxWeight = 0;

  for (int i = 0; i < n; i++) {
    double x = width * rand()/RAND_MAX;
    double y = height * rand()/RAND_MAX;
    double weight = std::min(width, height)/6 * rand()/RAND_MAX;
    maxWeight = std::max(maxWeight, weight);
    insertSite(ag, sitesFile, x, y, weight);
  }

  insertBoundingSites(ag, width, height, maxWeight);

  fclose(sitesFile);

	return ag;
}

// takes a sites file as input and thus doesn't need to create one
Apo_graph* generateApo_graphCustom(std::string sitesPath) {
  std::ifstream sitesStream(sitesPath.c_str());

  std::vector<std::vector<double> > sites(3);
  for (int i = 0; i < 3; ++i) {
    sites[i] = std::vector<double>();
  }

  double maxWeight = 0;
  double maxX = 0;
  double maxY = 0;

  while (sitesStream) {
    std::string currLine;
    getline(sitesStream, currLine);

    std::stringstream currLineStream;
    currLineStream << currLine;
    double currValue;

    for (int i = 0; i < 3; ++i) {
      if (!(currLineStream >> currValue)) {
        break;
      }
      sites[i].push_back(currValue);
    }
    maxX = std::max(maxX, sites[0].back());
    maxY = std::max(maxY, sites[1].back());
    maxWeight = std::max(maxWeight, sites[2].back());
  }

  Apo_graph* ag = new Apo_graph();

  insertBoundingSites(ag, maxX, maxY, maxWeight);

  int numSites = sites[0].size();

  for (int i = 0; i < numSites; ++i) {
    insertSite(ag, NULL, sites[0][i], sites[1][i], sites[2][i]);
  }

  return ag;
}

Apo_graph* generateApo_graphGrid(std::string weightsPath, std::string sitesPath) {
  FILE* sitesFile = fopen(sitesPath.c_str(), "w");

  std::ifstream weightsStream(weightsPath.c_str());
  int rows = -1, cols = 0;
  std::vector<double> weights;
  double maxWeight = 0;

  while (weightsStream) {
    rows++;
    std::string currLine;
    getline(weightsStream, currLine);

    std::stringstream currLineStream;
    currLineStream << currLine;
    double currWeight;

    while (currLineStream >> currWeight) {
      weights.push_back(currWeight);
      maxWeight = std::max(maxWeight, currWeight);
      if (rows == 0) {
        cols++;
      }
    }
  }

  Apo_graph* ag = new Apo_graph();

  insertBoundingSites(ag, 1, 1, maxWeight);

  // the coordinates of the points are divided by normalizingFactor to make them
  // be contained in the square [0,1]x[0,1]
  int normalizingFactor = std::max(rows, cols);

  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      double x = (col + 0.5) / normalizingFactor;
      double y = (rows - row - 0.5) / normalizingFactor;
      double weight = weights[row*cols + col];
      insertSite(ag, sitesFile, x, y, weight);
    }
  }

  fclose(sitesFile);

  return ag;
}

void printIntersections(Apo_graph* ag, std::string intersectionsPath) {
  FILE* intersectionsFile = fopen(intersectionsPath.c_str(), "w");

  Diagram diagram(*ag);

  // To prevent edges from being considered from both sides, we collect the
  // already visited ones in a set.
  std::set<Halfedge> edges;
  Halfedge currHalfedge, currOppositeHalfedge;
  // the two sites adjacent to a halfedge
  Site origin, neighbor;
  // the two endpoints of a halfedge
  Point source, target;

  int numEdges = 0, numHalfedges = 0;

  for (Face_iterator face_it = diagram.faces_begin();
    face_it != diagram.faces_end(); face_it++) {
    Face currFace = *face_it;
    origin = currFace.dual()->site();

    Ccb_halfedge_circulator edge_ci, edge_ci_end;
    edge_ci = edge_ci_end = currFace.ccb();
    do {
      numHalfedges++;
      currHalfedge = *edge_ci;
      currOppositeHalfedge = *currHalfedge.opposite();
      // the current edge hasn't been visited yet
      if (edges.find(currOppositeHalfedge) == edges.end()) {
        numEdges++;

        neighbor = currOppositeHalfedge.face()->dual()->site();

        if (!currHalfedge.has_source() || !currHalfedge.has_target()) {
          continue;
        }

        fprintf(intersectionsFile, "%.30f %.30f %.30f %.30f %.30f %.30f ",
          origin.point().x(), origin.point().y(), origin.weight(),
          neighbor.point().x(), neighbor.point().y(), neighbor.weight());

        // the first endpoint of the hyperbola
        if (currHalfedge.has_source()) {
          source = (*currHalfedge.source()).point();
          fprintf(intersectionsFile, "%.30f %.30f ", source.x(), source.y());
        } else {
          fprintf(intersectionsFile, "NaN NaN ");
        }

        // the second endpoint of the hyperbola
        if (currHalfedge.has_target()) {
          target = (*currHalfedge.target()).point();
          fprintf(intersectionsFile, "%.30f %.30f\n", target.x(), target.y());
        } else {
          fprintf(intersectionsFile, "NaN NaN\n");
        }
      }
      edges.insert(currHalfedge);
    } while(++edge_ci != edge_ci_end);
  }
  fclose(intersectionsFile);
}

int main(int argc, char* argv[]) {
  srand(time(NULL));

  Apo_graph* ag;
  std::string intersectionsPath;

  switch (argc) {
    default:
      {
      const char* helpMessage =
        "\nThis tool computes the Voronoi diagram for either a random set of sites, "
        "for points arranged on a grid with prescribed weights or for arbitrary "
        "points with associated weights.\n\n"
        "It creates two files, one containing the sites and one containing for each "
        "hyperbola segment between two sites the corresponding sites and the "
        "endpoints of the segment, i.e., the intersection points of the hyperbolas.\n\n"
        "sites file: <x-coordinate> <y-coordinate> <weight>\n"
        "intersections file: <x-coordinate 1> <y-coordinate 1> <weight 1> "
        "<x-coordinate 2> <y-coordinate 2> <weight 2> <1st endpoint> <2nd endpoint>\n"
        "<1st/2nd endpoint> == NaN means that the hyperbola doesn't intersect any "
        "other hyperbola on this side\n\n"
        "Usage:\n\n"
        "RANDOM: create_diagram <number of sites> <width> <height> <sites file> "
        "<intersections file>\n"
        "<number of sites> - the number of sites the diagram should contain\n"
        "<width>, <height> - the sites will be contained in the rectangle "
        "width x height\n"
        "<sites file> - the file to which the sites should be printed\n"
        "<intersections file> - the file to which the intersection points of the "
        "hyperbolas should be printed\n\n"
        "GRID: create_diagram <weight file> <sites file> <intersections file>\n"
        "<weight file> - a file containing the weights of sites arranged on a "
        "pixel grid, the former arranged exactly as the sites they belong to\n"
        "<sites file> - the file to which the sites should be printed\n"
        "<intersections file> - the file to which the intersection points of the "
        "hyperbolas should be printed\n\n"
        "ARBITRARY: create_diagram <sites file> <intersections file>\n"
        "<sites file> - the INPUT file with the sites\n"
        "               three columns: <x-coordinate> <y-coordinate> <weight>\n"
        "<intersections file> - the file to which the intersection points of the "
        "hyperbolas should be printed\n";
        std::cout << helpMessage;
      return -1;
      }
    case 6:
      {
      int n = std::atoi(argv[1]);
      int width = std::atoi(argv[2]);
      int height = std::atoi(argv[3]);
      std::string sitesPath = argv[4];
      intersectionsPath = argv[5];
      ag = generateApo_graphRandom(n, width, height, sitesPath);
      break;
      }
    case 4:
      {
      std::string weightsPath = argv[1];
      std::string sitesPath = argv[2];
      intersectionsPath = argv[3];
      ag = generateApo_graphGrid(weightsPath, sitesPath);
      break;
      }
    case 3:
      {
      std::string sitesPath = argv[1];
      intersectionsPath = argv[2];
      ag = generateApo_graphCustom(sitesPath);
      break;
      }
  }

  printIntersections(ag, intersectionsPath);
}
