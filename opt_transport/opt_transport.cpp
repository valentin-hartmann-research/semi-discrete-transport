/*
 * Author: Valentin Hartmann
 */

#include "Source.h"
#include "Target.h"

#include <lbfgs.h>

#include <vector>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>

/*
 * prints the computed Wasserstein distance
 */
#define PRINT_WASSERSTEIN_DISTANCE

// ***** debug output *****

/*
 * prints the current level of the coarsening of \nu
 */
// #define PRINT_LEVELS

/*
 * prints the return code of L-BFGS for each coarsening of \nu including the
 * original \nu
 * 1 = minimization was successful
 */
// #define PRINT_LBFGS_RETURN

/*
 * prints the return code of L-BFGS only for the original \nu
 * 1 = minimization was successful
 */
#define PRINT_FINAL_LBFGS_RETURN

/*
 * prints the error return codes of lbfgs
 */
// #define LBFGS_ERRORS

// ************************

struct MeasureData {
  std::vector<double> masses;
  double accMass;
};

struct GridMeasureData : MeasureData {
  int rows;
  int cols;
};

struct GeneralMeasureData : MeasureData {
  std::vector<std::pair<double, double> > points;
};

/*
 * Does the actual work for the functions below.
 * @param transportPlanPath:  "" means it is not being provided
 * @param normalizingFactor:  negative values means it is not being provided
 */
void _computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath,
  std::string transportPlanPath,
  double normalizingFactor);

/*
 * @param refinement:         the refinement ratio for the source measure
 * @param targetReduction:    |supp(\nu_l)| = |supp(\nu_{l-1})| / targetReduction
 * @param sourcePath:         the path of the file containing the source masses
 * @param targetPath:         the path of the file containing the target masses
 * @param normalizingFactor:  in the case that targetPath contains a list of
 *                            points and corresponding masses this is the factor
 *                            by which all coordinates are divided to make them
 *                            being contained in [0,1]x[0,1]
 *                            The caller must make sure that this is actually
 *                            the case.
 * @param weightsPath:        the file the weights vector should be printed to
 * @param transportPlanPath:  the file the transport plan should be printed to;
 *                            if omitted, the transport plan isn't printed
 */
void computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath,
  std::string transportPlanPath,
  double normalizingFactor) {
    _computeWeightVector(
      refinement,
      targetReduction,
      sourcePath,
      targetPath,
      weightsPath,
      transportPlanPath,
      normalizingFactor);
  }
void computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath,
  double normalizingFactor) {
    _computeWeightVector(
      refinement,
      targetReduction,
      sourcePath,
      targetPath,
      weightsPath,
      std::string(),
      normalizingFactor);
  }
void computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath,
  std::string transportPlanPath) {
    _computeWeightVector(
      refinement,
      targetReduction,
      sourcePath,
      targetPath,
      weightsPath,
      transportPlanPath,
      -1);
  }
void computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath) {
    _computeWeightVector(
      refinement,
      targetReduction,
      sourcePath,
      targetPath,
      weightsPath,
      std::string(),
      -1);
  }

/*
 * @param filepath: a file containing the masses of a measure on a grid
 * @return: a corresponding GridMeasureData object
 */
GridMeasureData* readGridFile(std::string filepath);

/*
 * The points are assumed to have non-negative coordinates.
 * Note that this function normalizes the points to be contained in [0,1]x[0,1]
 * in such a way that at least one point lies on the line (x,1) or (1,y).
 * @param filepath:           a file containing the points and corresponding
 *                            masses of a measure in a list
 * @param normalizingFactor:  the factor by which the coordinates should be
 *                            divided
 *                            After division they need to be contained in
 *                            [0,1]x[0,1].
 * @return: a corresponding GeneralMeasureData object
 */
GeneralMeasureData* readListFile(std::string filepath,
  double normalizingFactor);

/*
 * Normalizes a vector such that its smallest entry is 0 by adding a constant to
 * all entries.PRINT_LBFGS_RETURN
 * @param vector:     the vector to be normalized
 * @param numEntries: the number of entries of vector
 */
void normalizeVector(lbfgsfloatval_t* vector, int numEntries);

/*
 * Prints a hint on how to use the command line version of the program.
 */
void printUsageHint();


int main(int argc, char* argv[]) {
  srand(time(NULL));


  #ifdef LBFGS_ERRORS
  int codes[] ={  LBFGS_SUCCESS, LBFGS_CONVERGENCE,
                  LBFGS_STOP, LBFGS_ALREADY_MINIMIZED,
                  LBFGSERR_UNKNOWNERROR, LBFGSERR_LOGICERROR,
                  LBFGSERR_OUTOFMEMORY, LBFGSERR_CANCELED,
                  LBFGSERR_INVALID_N, LBFGSERR_INVALID_N_SSE,
                  LBFGSERR_INVALID_X_SSE, LBFGSERR_INVALID_EPSILON,
                  LBFGSERR_INVALID_TESTPERIOD, LBFGSERR_INVALID_DELTA,
                  LBFGSERR_INVALID_LINESEARCH, LBFGSERR_INVALID_MINSTEP,
                  LBFGSERR_INVALID_MAXSTEP, LBFGSERR_INVALID_FTOL,
                  LBFGSERR_INVALID_WOLFE, LBFGSERR_INVALID_GTOL,
                  LBFGSERR_INVALID_XTOL, LBFGSERR_INVALID_MAXLINESEARCH,
                  LBFGSERR_INVALID_ORTHANTWISE, LBFGSERR_INVALID_ORTHANTWISE_START,
                  LBFGSERR_INVALID_ORTHANTWISE_END, LBFGSERR_OUTOFINTERVAL,
                  LBFGSERR_INCORRECT_TMINMAX, LBFGSERR_ROUNDING_ERROR,
                  LBFGSERR_MINIMUMSTEP, LBFGSERR_MAXIMUMSTEP,
                  LBFGSERR_MAXIMUMLINESEARCH, LBFGSERR_MAXIMUMITERATION,
                  LBFGSERR_WIDTHTOOSMALL, LBFGSERR_INVALIDPARAMETERS,
                  LBFGSERR_INCREASEGRADIENT};
  std::string codeNames[] ={"LBFGS_SUCCESS", "LBFGS_CONVERGENCE",
                  "LBFGS_STOP", "LBFGS_ALREADY_MINIMIZED",
                  "LBFGSERR_UNKNOWNERROR", "LBFGSERR_LOGICERROR",
                  "LBFGSERR_OUTOFMEMORY", "LBFGSERR_CANCELED",
                  "LBFGSERR_INVALID_N", "LBFGSERR_INVALID_N_SSE",
                  "LBFGSERR_INVALID_X_SSE", "LBFGSERR_INVALID_EPSILON",
                  "LBFGSERR_INVALID_TESTPERIOD", "LBFGSERR_INVALID_DELTA",
                  "LBFGSERR_INVALID_LINESEARCH", "LBFGSERR_INVALID_MINSTEP",
                  "LBFGSERR_INVALID_MAXSTEP", "LBFGSERR_INVALID_FTOL",
                  "LBFGSERR_INVALID_WOLFE", "LBFGSERR_INVALID_GTOL",
                  "LBFGSERR_INVALID_XTOL", "LBFGSERR_INVALID_MAXLINESEARCH",
                  "LBFGSERR_INVALID_ORTHANTWISE", "LBFGSERR_INVALID_ORTHANTWISE_START",
                  "LBFGSERR_INVALID_ORTHANTWISE_END", "LBFGSERR_OUTOFINTERVAL",
                  "LBFGSERR_INCORRECT_TMINMAX", "LBFGSERR_ROUNDING_ERROR",
                  "LBFGSERR_MINIMUMSTEP", "LBFGSERR_MAXIMUMSTEP",
                  "LBFGSERR_MAXIMUMLINESEARCH", "LBFGSERR_MAXIMUMITERATION",
                  "LBFGSERR_WIDTHTOOSMALL", "LBFGSERR_INVALIDPARAMETERS",
                  "LBFGSERR_INCREASEGRADIENT"};
  for (int i = 0; i < 35; ++i) {
    printf("%d - %s\n", codes[i], codeNames[i].c_str());
  }
  #endif

  // ***************************************************************************

  bool targetAsList;

  if (argc < 5) {
    printUsageHint();
    return -1;
  } else {
    std::string modifier = argv[1];
    if (modifier == "-l") {
      targetAsList = true;
    } else if (modifier == "-g") {
      targetAsList = false;
    } else {
      printUsageHint();
      return -1;
    }
  }

  int refinement = 1000;
  int targetReduction = 5;

  std::string sourcePath = argv[2];
  std::string targetPath = argv[3];
  if (targetAsList) {
    if (argc < 6) {
      printUsageHint();
      return -1;
    } else {
      std::istringstream normalizingFactorStream(argv[4]);
      double normalizingFactor;
      if (!(normalizingFactorStream >> normalizingFactor)) {
        printUsageHint();
        return -1;
      }
      std::string weightsPath = argv[5];
      if (argc == 6) {
        computeWeightVector(refinement, targetReduction, sourcePath, targetPath,
          weightsPath, normalizingFactor);
      } else {
        std::string transportPlanPath = argv[6];
        computeWeightVector(refinement, targetReduction, sourcePath, targetPath,
          weightsPath, transportPlanPath, normalizingFactor);
      }
    }
  } else {
    std::string weightsPath = argv[4];
    if (argc == 5) {
      computeWeightVector(refinement, targetReduction, sourcePath, targetPath,
        weightsPath);
    } else {
      std::string transportPlanPath = argv[5];
      computeWeightVector(refinement, targetReduction, sourcePath, targetPath,
        weightsPath, transportPlanPath);
    }
  }
}


void printUsageHint() {
  printf("Usage: opt_transport <-l/-g> <source> <target> [<normalizing factor for list target>] <weights> [<transport plan>]\n");
}

void _computeWeightVector(
  int refinement,
  int targetReduction,
  std::string sourcePath,
  std::string targetPath,
  std::string weightsPath,
  std::string transportPlanPath,
  double normalizingFactor)
  {

  bool targetAsList;
  if (normalizingFactor < 0) {
    targetAsList = false;
  } else {
    targetAsList = true;
  }

  GridMeasureData* sourceData = readGridFile(sourcePath);
  Source source(sourceData->masses, sourceData->rows, sourceData->cols,
    sourceData->accMass);

  MeasureData* targetData;
  Target* target;

  if (targetAsList) {
    targetData = readListFile(targetPath, normalizingFactor);
    GeneralMeasureData* generalTargetData = static_cast<GeneralMeasureData*>(
      targetData);
    target = new Target(generalTargetData->points, generalTargetData->masses,
      generalTargetData->accMass);
  } else {
    targetData = readGridFile(targetPath);
    GridMeasureData* gridTargetData = static_cast<GridMeasureData*>(targetData);
    target = new Target(gridTargetData->rows, gridTargetData->cols,
      gridTargetData->masses, gridTargetData->accMass);
  }


  // ***** the coarsening of the target measure *****

  int numberOfMeasures = log(targetData->masses.size()) / log(targetReduction);
  // if we increased numberOfMeasures also for <= 2,
  // we would have |supp(\nu_L)| == 1
  if (targetData->masses.size() / pow(5, numberOfMeasures) >= 2) {
    ++numberOfMeasures;
  }
  typedef std::vector<std::pair<Target*, std::vector<std::vector<int> > > > Coarsening;
  Coarsening coarsening(numberOfMeasures);
  coarsening[0] = std::make_pair(target, std::vector<std::vector<int> >());
  for (int l = 1; l < coarsening.size(); ++l) {
    coarsening[l] = coarsening[l - 1].first->coarsen(targetReduction);
  }

  // the start weights of \nu_L
  lbfgsfloatval_t* startWeights = lbfgs_malloc(coarsening.back().second.size());
  for (int i = 0; i < coarsening.back().second.size(); ++i) {
    startWeights[i] = 0;
  }

  // ***************************************************


  // ***** the minimization of \Phi *****

  int lbfgs_return;

  for (int l = coarsening.size() - 1; l > 0; --l) {
    #ifdef PRINT_LEVELS
    printf("***********\nLevel of Coarsening: %d\n***********\n\n", l);
    #endif

    lbfgs_return = source.optimize(coarsening[l].first, startWeights, refinement);

    #ifdef PRINT_LBFGS_RETURN
    printf("L-BFGS return code: %d\n\n", lbfgs_return);
    #endif

    int lengthOldStartWeights = coarsening[l].first->getPoints().size();
    lbfgsfloatval_t* oldStartWeights = lbfgs_malloc(lengthOldStartWeights);
    for (int i = 0; i < lengthOldStartWeights; ++i) {
      oldStartWeights[i] = startWeights[i];
    }

    lbfgs_free(startWeights);
    startWeights = lbfgs_malloc(coarsening[l - 1].first->getPoints().size());

    for (int i = 0; i < coarsening[l].second.size(); ++i) {
      // the points of \nu_{l-1} that belong to the i-th point of \nu_l
      std::vector<int> currAssignedPoints = coarsening[l].second[i];
      // the weight of the i-th point of \nu_l
      double currWeight = oldStartWeights[i];
      for (int j = 0; j < currAssignedPoints.size(); ++j) {
        startWeights[currAssignedPoints[j]] = currWeight;
      }
    }
    lbfgs_free(oldStartWeights);
  }

  #ifdef PRINT_LEVELS
  printf("***********\nLevel of Coarsening: 0\n***********\n\n");
  #endif

  lbfgs_return =  source.optimize(coarsening[0].first, startWeights, refinement);

  #if defined PRINT_LBFGS_RETURN || defined PRINT_FINAL_LBFGS_RETURN
  printf("L-BFGS return code: %d\n", lbfgs_return);
  #endif

  #ifdef PRINT_WASSERSTEIN_DISTANCE
  printf("Wasserstein distance: %.20f\n", source.getWasserstein());
  #endif

  // ************************************


  // ***** print the weight vector *****

  int precision = 10;

  FILE* weightsFile = fopen(weightsPath.c_str(), "w");

  normalizeVector(startWeights, target->getPoints().size());

  int rows, cols;
  if (targetAsList) {
    rows = target->getPoints().size();
    cols = 1;
  } else {
    GridMeasureData* gridTargetData = static_cast<GridMeasureData*>(targetData);
    rows = gridTargetData->rows;
    cols = gridTargetData->cols;
  }

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      fprintf(weightsFile, "%.*f ", precision, startWeights[cols*i + j]);
    }
    fprintf(weightsFile, "\n");
  }

  fclose(weightsFile);

  if (!transportPlanPath.empty()) {
    source.createTransportPlan(coarsening[0].first, startWeights, transportPlanPath);
  }

  lbfgs_free(startWeights);

  for (Coarsening::iterator it = coarsening.begin(); it != coarsening.end(); ++it) {
    delete (*it).first;
  }
  delete targetData;
}

GridMeasureData* readGridFile(std::string filepath) {
  std::ifstream fileStream(filepath.c_str());
  GridMeasureData* m = new GridMeasureData;
  m->accMass = 0;
  m->rows = -1, m->cols = 0;

  while (fileStream) {
    m->rows++;
    std::string currLine;
    getline(fileStream, currLine);

    std::stringstream currLineStream;
    currLineStream << currLine;
    double currMass;

    while (currLineStream >> currMass) {
      m->masses.push_back(currMass);
      m->accMass += currMass;
      if (m->rows == 0) {
        m->cols++;
      }
    }
  }

  return m;
}

GeneralMeasureData* readListFile(std::string filepath,
  double normalizingFactor) {
  std::ifstream fileStream(filepath.c_str());
  GeneralMeasureData* m = new GeneralMeasureData;
  m->accMass = 0;

  while (fileStream) {
    std::string currLine;
    getline(fileStream, currLine);
    if (currLine.empty()) {
      break;
    }

    std::stringstream currLineStream;
    currLineStream << currLine;
    double currX, currY, currMass;
    currLineStream >> currX;
    currLineStream >> currY;
    currLineStream >> currMass;

    m->points.push_back(std::make_pair(currX, currY));
    m->masses.push_back(currMass);
    m->accMass += currMass;
  }

  // normalize the points to be contained in [0,1]x[0,1]
  for (std::vector<std::pair<double, double> >::iterator it = m->points.begin();
    it != m->points.end(); ++it) {
    *it = std::make_pair(it->first / normalizingFactor,
      it->second / normalizingFactor);
  }

  return m;
}

void normalizeVector(lbfgsfloatval_t* vector, int numEntries) {
  lbfgsfloatval_t min = vector[0];
  for (int i = 1; i < numEntries; ++i) {
    min = std::min(min, vector[i]);
  }
  for (int i = 0; i < numEntries; ++i) {
    vector[i] -= min;
  }
}
