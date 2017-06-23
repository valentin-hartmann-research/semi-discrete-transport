#include <fstream>
#include <sstream>
#include <string>

struct Site{
  double x;
  double y;
  std::string weight;
};

void printUsageHint() {
  printf("Usage: create_sites_file <target measure> <normalizing factor> <weights> <sites>");
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    printUsageHint();
    return -1;
  }

  std::string measurePath = argv[1];
  std::ifstream measureStream(measurePath.c_str());

  std::istringstream normalizingFactorStream(argv[2]);
  double normalizingFactor;
  normalizingFactorStream >> normalizingFactor;

  std::string weightsPath = argv[3];
  std::ifstream weightsStream(weightsPath.c_str());
  std::string targetPath = argv[4];
  FILE* targetFile = fopen(targetPath.c_str(), "w");

  std::string currMeasureLine;
  while (getline(measureStream, currMeasureLine)) {
    std::string currWeightsLine;
    getline(weightsStream, currWeightsLine);

    Site currSite;
    currSite.weight = currWeightsLine;

    std::stringstream currMeasureLineStream;
    currMeasureLineStream << currMeasureLine;
    double currX, currY;
    currMeasureLineStream >> currX;
    currMeasureLineStream >> currY;
    currSite.x = currX / normalizingFactor;
    currSite.y = currY / normalizingFactor;

    // const char* currMeasureLinePointer = currMeasureLine.c_str();
    // const char* begin = currMeasureLinePointer;
    // while(*currMeasureLinePointer != ' ') {
    //   ++currMeasureLinePointer;
    // }
    // currSite.x = std::string(begin, currMeasureLinePointer);
    // ++currMeasureLinePointer;
    // begin = currMeasureLinePointer;
    // while(*currMeasureLinePointer != ' ') {
    //   ++currMeasureLinePointer;
    // }
    // currSite.y = std::string(begin, currMeasureLinePointer);


    // int endX = currMeasureLine.find(" ", 0);
    // int endY = currMeasureLine.find(" ", endX);
    // currSite.x = currMeasureLine.substr(0, endX);
    // currSite.y = currMeasureLine.substr(endX + 1, endY);

    // fprintf(targetFile, "%s %s %s\n", currSite.x.c_str(), currSite.y.c_str(), currSite.weight.c_str());
    fprintf(targetFile, "%.30f %.30f %s\n", currSite.x, currSite.y, currSite.weight.c_str());
  }

  fclose(targetFile);
}
