/**
 * p4.cpp
 */

#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "vision_utilities.h"

#define NUM_IMG 5
#define DIM 3

using namespace std;

Image _images[NUM_IMG];
int _i, _j;

bool intensity_at(int img1, int img2) { return getPixel(_images+img1, _i, _j) > getPixel(_images+img2, _i, _j); }

bool zero_at(int img) { return 0 == getPixel(_images+img, _i, _j); }

int map_normal(float x) { return 255*(1.0f + x)/2; }

void matrix_inverse(float M[DIM][DIM], float M_inv[DIM][DIM]);

int main(int argc, char *argv[])
{
    if (argc < 11) {
        cerr << "usage: " << *argv << " <input directions> <image 1> <image 2> <image 3> <image 4> <image 5> <mask> <normal map> <albedo map> <output gradient>\n";
        exit(0);
    }

    int was_error = 0;

    ifstream direction_file(argv[1]);
    float directions[NUM_IMG][DIM];
    vector<int> image_indices;
    ImageColor normals;
    Image mask, albedo;


    for (int i = 0; i < NUM_IMG; ++i) {
        if (readImage(_images+i, argv[i+2]) == -1) {
            cerr << "error reading " << argv[i+2] << ". exiting...\n";
            exit(1);
        }
        for (int j = 0; j < DIM; ++j) {
            direction_file >> directions[i][j];
        }
        image_indices.push_back(i);
    }

    direction_file.close();

    vector<int>::iterator first_zero;
    readImage(&mask, argv[7]);
    int rows = getNRows(&mask), cols = getNCols(&mask), k;
    float x, y, z, magnitude, albedo_max = 0.0f, *albedo_store=(float *)calloc(sizeof(float), rows*cols),
          S[DIM][DIM], I[DIM], M[DIM], S_inv[DIM][DIM];
    Gradient *gradients=(Gradient *)calloc(sizeof(Gradient), rows*cols);
    
    if (setSizeColor(&normals, rows, cols) == -1) {
        cerr << "out of memory. exiting...\n";
        exit(1);
    }

    if (setSize(&albedo, rows, cols) == -1) {
        cerr << "out of memory. exiting...\n";
        exit(1);
    }

    for (_i = 0; _i < rows; ++_i) {
        for (_j = 0; _j < cols; ++_j) {
            sort(image_indices.begin(), image_indices.end(), intensity_at);
            if (getPixel(&mask, _i, _j)) {
                for (int i = 0; i < DIM; ++i) {
                    k = image_indices[i];
                    memcpy(S[i], directions[k], sizeof(directions[k]));
                    I[i] = getPixel(_images+k, _i, _j);
                }

                matrix_inverse(S, S_inv);
                magnitude = 0;

                for (int i = 0; i < DIM; ++i) {
                    M[i] = 0.0f;
                    for (int j = 0; j < DIM; ++j) {
                        M[i] += S_inv[i][j]*I[j];
                    }
                    magnitude += pow(M[i], 2);
                }

                magnitude = sqrt(magnitude);
                for (int i = 0; i < DIM; ++i)
                    M[i] /= magnitude;

                x = map_normal(M[0]);
                y = map_normal(M[1]);
                z = map_normal(M[2]);

                if (magnitude > albedo_max) albedo_max = magnitude;

                setPixelColor(&normals, _i, _j, x, y, z);
                albedo_store[_i*cols + _j] = magnitude;
                gradients[_i*cols + _j].p = M[0];
                gradients[_i*cols + _j].q = M[1];

            } else {
                setPixelColor(&normals, _i, _j, 0, 0, 255);
            }
        }
    }

    for (_i = 0; _i < rows; ++_i) {
        for (_j = 0; _j < cols; ++_j) {
            setPixel(&albedo, _i, _j, 255*(albedo_store[_i*cols + _j]/albedo_max));
        }
    }

    setColorsColor(&normals, 255);
    setColors(&albedo, 255);
    if (writeImageColor(&normals, argv[8]) == -1) {
        cerr << "error writing normal map.\n";
        was_error = 1;
    }

    if (writeImage(&albedo, argv[9]) == -1) {
        cerr << "error writing albedo map.\n";
        was_error = 1;
    }

    ofstream gradient_file(argv[10], fstream::binary);
    gradient_file.write((char *)gradients, rows*cols*sizeof(Gradient));
    gradient_file.close();

    for (int i = 0; i < NUM_IMG; ++i) {
        free(_images[i].data);
    }

    free(normals.dataR);
    free(normals.dataG);
    free(normals.dataB);
    free(albedo.data);
    free(mask.data);
    free(albedo_store);
    free(gradients);

    return was_error;
}

void matrix_inverse(float M[DIM][DIM], float M_inv[DIM][DIM])
{
    int i, j, k, l, row_max;
    float m, tmp_row[DIM];
    for (i = 0; i < DIM; ++i) {
        for (j = 0; j < DIM; ++j) {
            M_inv[i][j] = i == j ? 1.0 : 0.0;
        }
    }
    i = j = 0;
    while (i < DIM && j < DIM) {
        row_max = i;
        for (k = i+1; k < DIM; ++k) {
            if (fabs(M[k][j]) > fabs(M[row_max][j])) {
                row_max = k;
            }
        }
        if (M[row_max][j]) {
            if (row_max != i) {
                memcpy(tmp_row, M_inv[row_max], sizeof(tmp_row));
                memcpy(M_inv[row_max], M_inv[i], sizeof(M_inv[i]));
                memcpy(M_inv[i], tmp_row, sizeof(tmp_row));

                memcpy(tmp_row, M[row_max], sizeof(tmp_row));
                memcpy(M[row_max], M[i], sizeof(M[i]));
                memcpy(M[i], tmp_row, sizeof(tmp_row));
            }
            m = M[i][j];
            for (k = 0; k < DIM; ++k) {
                M_inv[i][k] /= m;
                M[i][k] /= m;
            }
            for (k = 0; k < DIM; ++k) {
                if (k != i) {
                    m = M[k][j];
                    for (l = 0; l < DIM; ++l) {
                        M_inv[k][l] -= m*M_inv[i][l];
                        M[k][l] -= m*M[i][l];
                    }
                }
            }
            i++;
        }
        j++;
    }
}
