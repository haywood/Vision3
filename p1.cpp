/**
 * p1.cpp
 */

#include "vision_utilities.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        cerr << "usage: " << argv[0] << " <input original image> <output parameters file>\n";
        exit(0);
    }

    char * input = argv[1], * output = argv[2];
    Image im;

    if (readImage(&im, input) == -1) {
        cerr << "error reading " << input << ". exiting...\n";
        exit(1);
    }

    int rows = getNRows(&im), cols = getNCols(&im),
        left = numeric_limits<int>::max(), 
        right = 0, top = -1;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (getPixel(&im, i, j)) {
                if (top < 0) top = i;
                if (j < left) left = j;
                if (j > right) right = j;
            }
        }
    }

    int radius = (right - left)/2, x0 = left + radius, y0 = top + radius;
    ofstream parameter_file;
    parameter_file.open(output);
    parameter_file << x0 << " " << y0 << " " << radius;
    parameter_file.close();

    free(im.data);

    return 0;
}
