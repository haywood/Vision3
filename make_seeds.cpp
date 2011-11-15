/**
 * make_seeds.cpp
 * program for generating seed points for p5.
 */

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "vision_utilities.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 4) {
        cerr << "usage: " << *argv << " <image mask> <seed file> <number of points>\n";
        exit(0);
    }

    Image mask;
    readImage(&mask, argv[1]);
    ofstream seed_file(argv[2]);
    int n_points = 0, rows = getNRows(&mask), cols = getNCols(&mask);
    int i, j, N = atoi(argv[3]);

    while (n_points < N) {
        i = rand() % rows; j = rand() % cols;
        if (getPixel(&mask, i, j)) {
            seed_file << j << " " << i << "\n";
            n_points++;
        }
    }

    seed_file.close();
    free(mask.data);
}
