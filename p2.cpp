/**
 * p2.cpp
 */

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "vision_utilities.h"

#define NUM_IMG 5

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 8) {
        cerr << "usage: " << *argv << "<input parameters file> <image 1> <image 2> <image 3> <image 4> <image 5> <output directions file>\n";
        exit(0);
    }

    ifstream parameter_file;
    parameter_file.open(argv[1]);
    int x0, y0, radius;
    parameter_file >> x0;
    parameter_file >> y0;
    parameter_file >> radius;
    parameter_file.close();

    int rows, cols, x, y, z, max_intensity, intensity;
    float scale;
    ofstream directions_file;
    directions_file.open(argv[7]);
    Image image;

    for (int l = 0; l < NUM_IMG; ++l) {
        if (readImage(&image, argv[l+2]) == -1) {
            cerr << "error reading " << argv[l+2] << ". exiting...\n";
            exit(1);
        }

        rows = getNRows(&image);
        cols = getNCols(&image);

        max_intensity = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                intensity = getPixel(&image, i, j);
                if (intensity > max_intensity) {
                    max_intensity = intensity;
                    x = j;
                    y = i;
                }
            }
        }

        x -= x0;
        y -= y0;
        z = -sqrt(pow(radius, 2) - (pow(x, 2) + pow(y, 2)));
        scale = -(float)max_intensity/radius;
        directions_file << -scale*x << " " << -scale*y << " " << scale*z;
        if (l+1 < NUM_IMG) directions_file << "\n";

        free(image.data);
    }
    directions_file.close();


    return 0;
}
