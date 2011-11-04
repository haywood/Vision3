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

void matrix_inverse(float **M, float **M_inv, int n);
void matrix_transpose(float **M, float **M_T, int n, int m);
void naive_product(float **A, float **B, float **C, int n, int m, int p);

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

    readImage(&mask, argv[7]);
    int rows = getNRows(&mask), cols = getNCols(&mask), k, nonzero, N;
    float x, y, z, magnitude, albedo_max = 0.0f, *albedo_store=(float *)calloc(sizeof(float), rows*cols);
    float **S, *I, M[DIM], **S_T, **A, **A_inv, **S_inv;
    float p, q, f, g;
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
                nonzero = 1;
                N = 0;
                for (int i = 0; i < NUM_IMG && nonzero; ++i) {
                    k = image_indices[i];
                    if (getPixel(_images+k, _i, _j)) N++;
                    else nonzero = 0;
                }

                S = (float **)malloc(N*sizeof(float *));
                I = (float *)malloc(N*sizeof(float));

                S_T = (float **)malloc(DIM*sizeof(float *));
                for (int i = 0; i < DIM; ++i)
                    S_T[i] = (float *)malloc(N*sizeof(float));

                A = (float **)malloc(DIM*sizeof(float *));
                A_inv = (float **)malloc(DIM*sizeof(float *));
                S_inv = (float **)malloc(DIM*sizeof(float *));
                for (int i = 0; i < DIM; ++i) {
                    A[i] = (float *)malloc(DIM*sizeof(float));
                    A_inv[i] = (float *)malloc(DIM*sizeof(float));
                    S_inv[i] = (float *)malloc(N*sizeof(float));
                }

                for (int i = 0; i < N; ++i) {
                    S[i] = (float *)malloc(DIM*sizeof(float));
                    k = image_indices[i];
                    memcpy(S[i], directions[k], sizeof(directions[k]));
                    I[i] = getPixel(_images+k, _i, _j);
                }

                matrix_transpose(S, S_T, N, DIM);
                naive_product(S_T, S, A, DIM, N, DIM);
                matrix_inverse(A, A_inv, DIM);
                naive_product(A_inv, S_T, S_inv, DIM, DIM, N);
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

                setPixelColor(&normals, _i, _j, x, y, z);

                if (magnitude > albedo_max) albedo_max = magnitude;
                albedo_store[_i*cols + _j] = magnitude;

                p = M[0]/M[2];
                q = M[1]/M[2];

                magnitude = sqrt(1 + pow(p, 2) + pow(q, 2));

                f = 2*p/(1 + magnitude);
                g = 2*q/(1 + magnitude);

                gradients[_i*cols + _j].p = f;
                gradients[_i*cols + _j].q = g;

                for (int i = 0; i < N; ++i)
                    free(S[i]);
                free(S);
                
                free(I);

                for (int i = 0; i < DIM; ++i)
                    free(A[i]);
                free(A);

                for (int i = 0; i < DIM; ++i)
                    free(S_T[i]);
                free(S_T);

                for (int i = 0; i < DIM; ++i)
                    free(A_inv[i]);
                free(A_inv);

                for (int i = 0; i < DIM; ++i)
                    free(S_inv[i]);
                free(S_inv);

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

void matrix_inverse(float **M, float **M_inv, int n)
{
    int i, j, k, l, row_max;
    float m, *tmp_row =(float *)malloc(n*sizeof(float));
    for (i = 0; i < DIM; ++i) {
        for (j = 0; j < DIM; ++j) {
            M_inv[i][j] = i == j ? 1.0 : 0.0;
        }
    }
    i = j = 0;
    while (i < n && j < n) {
        row_max = i;
        for (k = i+1; k < n; ++k) {
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
            for (k = 0; k < n; ++k) {
                M_inv[i][k] /= m;
                M[i][k] /= m;
            }
            for (k = 0; k < n; ++k) {
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
    free(tmp_row);
}

void matrix_transpose(float **M, float **M_T, int n, int m)
{
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            M_T[i][j] = M[j][i];
        }
    }
}

void naive_product(float **A, float **B, float **C, int n, int m, int p)
{

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i][j] = 0.0f;
            for (int k = 0; k < m; ++k) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}
