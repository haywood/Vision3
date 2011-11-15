/**
 * p5.cpp
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <limits>

#include "vision_utilities.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 5) {
        cerr << "usage: " << *argv << " <input gradinet file> <seed point file> <binary mask> <output depth map>\n";
        exit(0);
    }

    Image mask, depth_image;
    if (readImage(&mask, argv[3]) == -1) {
        cerr << "error reading " << argv[3] << ". exiting...\n";
        exit(1);
    }

    int rows = getNRows(&mask), cols = getNCols(&mask);
    Gradient *gradients = (Gradient *)calloc(sizeof(Gradient), rows*cols);
    ifstream gradient_file(argv[1]);
    if (!gradient_file) {
        cerr << "error reading " << argv[1] << ". exiting...\n";
        exit(1);
    }

    gradient_file.read((char *)gradients, rows*cols*sizeof(Gradient));
    gradient_file.close();

    float p, q, depth, depth_max, depth_min, *depth_map, *avg_depth, avg;
    int x0, y0, zx, zy, num_seeds = 0;

    ifstream seed_file(argv[2]);
    if (!seed_file) {
        cerr << "error reading " << argv[2] << ". exiting...\n";
        exit(1);
    }
    
    avg_depth = (float *)calloc(sizeof(float), rows*cols);
    depth_map = (float *)calloc(sizeof(float), rows*cols);
    if (!depth_map || !avg_depth) {
        cerr << "out of memory\n";
        exit(1);
    }
    
    while ((seed_file >> x0) && (seed_file >> y0)) {

        depth_max = depth_min = depth_map[y0*cols + x0] = 0.0f;

        for (int dy = -1; dy < 2; dy+=2) {
            // column of the seed point
            for (int y = y0 + dy; y >= 0 && y < rows; y+=dy) {
                if (getPixel(&mask, y, x0)) {
                    q = gradients[y*cols + x0].q;
                    zy = depth_map[(y-dy)*cols + x0];
                    depth = depth_map[y*cols + x0] = zy - dy*q;
                    if (depth > depth_max) depth_max = depth;
                    if (depth < depth_min) depth_min = depth;
                }
            }
        }
        
        for (int dx = -1; dx < 2; dx+=2) {
            // row of the seed point
            for (int x = x0 + dx; x >= 0 && x < cols; x+=dx) {
                if (getPixel(&mask, y0, x)) {
                    p = gradients[y0*cols + x].p;
                    zx = depth_map[y0*cols + x - dx];
                    depth = depth_map[y0*cols + x] = zx - dx*p;
                    if (depth > depth_max) depth_max = depth;
                    if (depth < depth_min) depth_min = depth;
                }
            }
        }


        //integrate quadrants
        for (int dy = -1; dy < 2; dy+=2) {
            for (int dx = -1; dx < 2; dx+=2) {                
                // current quadrant
                for (int y = y0 + dy; y >= 0 && y < rows; y+=dy) {
                    for (int x = x0 + dx; x >= 0 && x < cols; x+=dx) {
                    if (getPixel(&mask, y, x)) {
                            p = gradients[y*cols + x].p;
                            q = gradients[y*cols + x].q;
                            zx = depth_map[y*cols + (x-dx)];
                            zy = depth_map[(y-dy)*cols + x];
                            depth = 0.5*(zx - dx*p + zy - dy*q);
                            depth_map[y*cols + x] = depth;
                            if (depth > depth_max) depth_max = depth;
                            if (depth < depth_min) depth_min = depth;
                        }
                    }
                }
            }
        }
        
        for (int i = 0; i < rows*cols; ++i) {
            if (depth_map[i])
                avg_depth[i] += (depth_map[i] - depth_min)/(depth_max - depth_min);
        }

        num_seeds++;
    }

    seed_file.close();

    if (setSize(&depth_image, rows, cols) == -1) {
        cerr << "out of memory.\n";
        exit(1);
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            avg = avg_depth[i*cols + j]/num_seeds;
            setPixel(&depth_image, i, j, 255*avg);
        }
    }

    setColors(&depth_image, 255);
    writeImage(&depth_image, argv[4]);
    free(depth_image.data);
    free(avg_depth);
    free(depth_map);
    free(mask.data);

    return 0;
}