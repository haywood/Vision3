/**
 * p3.cpp
 */

#include "vision_utilities.h"
#include <iostream>
#include <cstdlib>

#define NUM_IMG 5
#define NUM_MIN 3

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 7) {
        cerr << "usage: " << *argv << "<image 1> <image 2> <image 3> <image 4> <image 5> <output mask>\n";
        exit(0);
    }

    Image images[NUM_IMG], mask;
    for (int i = 0; i < NUM_IMG; ++i) {
        if (readImage(images+i, argv[i+1]) == -1) {
            cerr << "error reading " << argv[i+1] << ". exiting...\n";
            exit(1);
        }
    }

    char * output_file = argv[NUM_IMG+1];
    int rows = getNRows(images), cols = getNCols(images), count;
    setSize(&mask, rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            count = 0;
            for (int k = 0; k < NUM_IMG && count < NUM_MIN; ++k) {
                if (getPixel(images+k, i, j))
                    count++;
            }
            setPixel(&mask, i, j, count >= NUM_MIN ? 1 : 0);
        }
    }
    setColors(&mask, 1);

    writeImage(&mask, output_file);
    for (int i = 0; i < NUM_IMG; ++i) {
        free(images[i].data);
    }
    free(mask.data);

    return 0;
}
