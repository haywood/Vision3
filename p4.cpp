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

int map_normal(float x) { return 255*(2.0f + x)/4.0f; }

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
    float p, q, f, g, h;
    Gradient *gradients=(Gradient *)calloc(sizeof(Gradient), rows*cols);
    
    if (setSizeColor(&normals, rows, cols) == -1) {
        cerr << "out of memory. exiting...\n";
        exit(1);
    }

    if (setSize(&albedo, rows, cols) == -1) {
        cerr << "out of memory. exiting...\n";
        exit(1);
    }

    S = (float **)malloc(NUM_IMG*sizeof(float *));
    for (int i = 0; i < NUM_IMG; ++i)
        S[i] = (float *)malloc(DIM*sizeof(float));

    S_T = (float **)malloc(DIM*sizeof(float *));
    for (int i = 0; i < DIM; ++i)
        S_T[i] = (float *)malloc(NUM_IMG*sizeof(float));

    I = (float *)malloc(NUM_IMG*sizeof(float));

    A = (float **)malloc(DIM*sizeof(float *));
    A_inv = (float **)malloc(DIM*sizeof(float *));
    S_inv = (float **)malloc(DIM*sizeof(float *));
    for (int i = 0; i < DIM; ++i) {
        A[i] = (float *)malloc(DIM*sizeof(float));
        A_inv[i] = (float *)malloc(DIM*sizeof(float));
        S_inv[i] = (float *)malloc(NUM_IMG*sizeof(float));
    }


    for (_i = 0; _i < rows; ++_i) {
        for (_j = 0; _j < cols; ++_j) {
            if (getPixel(&mask, _i, _j)) {
                sort(image_indices.begin(), image_indices.end(), intensity_at);
                nonzero = 1;
                N = 0;
                for (int i = 0; i < NUM_IMG && nonzero; ++i) {
                    k = image_indices[i];
                    if (getPixel(_images+k, _i, _j)) {
                        memcpy(S[i], directions[k], sizeof(directions[k]));
                        I[i] = getPixel(_images+k, _i, _j);
                        N++;
                    }
                    else nonzero = 0;
                }

                /* pseudo inverse */
                matrix_transpose(S, S_T, N, DIM);
                naive_product(S_T, S, A, DIM, N, DIM);
                matrix_inverse(A, A_inv, DIM);
                naive_product(A_inv, S_T, S_inv, DIM, DIM, N);
                magnitude = 0.0f;

                for (int i = 0; i < DIM; ++i) {
                    M[i] = 0.0f;
                    for (int j = 0; j < N; ++j) {
                        M[i] += S_inv[i][j]*I[j];
                    }
                    magnitude += pow(M[i], 2);
                }

                magnitude = sqrt(magnitude);

                /* save albedo */
                albedo_store[_i*cols + _j] = magnitude;
                if (magnitude > albedo_max) 
                    albedo_max = magnitude;

                for (int i = 0; i < DIM; ++i)
                    M[i] /= magnitude;

                p = M[0]/M[2];
                q = M[1]/M[2];

                magnitude = sqrt(1 + pow(p, 2) + pow(q, 2));

                f = 2.0f*p/(1 + magnitude);
                g = 2.0f*q/(1 + magnitude);
                h = 2.0f/(1 + magnitude);

                x = map_normal(f);
                y = map_normal(g);
                z = map_normal(h);

                /* save normals */
                setPixelColor(&normals, _i, _j, x, y, z);

                /* save gradients */
                gradients[_i*cols + _j].p = f;
                gradients[_i*cols + _j].q = g;

            } else {
                setPixelColor(&normals, _i, _j, 0, 0, 255);
            }
        }
    }

    for (int i = 0; i < NUM_IMG; ++i)
        free(S[i]);
    free(S);

    for (int i = 0; i < DIM; ++i) {
        free(S_inv[i]);
        free(A_inv[i]);
        free(S_T[i]);
        free(A[i]);
    }
    free(A_inv);
    free(S_inv);
    free(S_T);
    free(A);
    free(I);

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
    float *pivrow = (float *)malloc(n*sizeof(float)),
          m1, m2;
    int imax;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M_inv[i][j] = i == j ? 1 : 0;
    for (int i = 0; i < n; ++i) {
        imax = i;
        for (int j = i + 1; j < n; ++j)
            if (fabs(M[j][i]) > fabs(M[imax][i]))
                imax = j;
        if (imax != i) {
            memcpy(pivrow, M[imax], n*sizeof(float));
            memcpy(M[imax], M[i], n*sizeof(float));
            memcpy(M[i], pivrow, n*sizeof(float));

            memcpy(pivrow, M_inv[imax], n*sizeof(float));
            memcpy(M_inv[imax], M_inv[i], n*sizeof(float));
            memcpy(M_inv[i], pivrow, n*sizeof(float));
        }

        m1 = M[i][i];
        for (int j = 0; j < n; ++j) {
            M_inv[i][j] /= m1;
            M[i][j] /= m1;
        }
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                m2 = M[j][i];
                for (int k = 0; k < n; ++k) {
                    M[j][k] -= m2*M[i][k];
                    M_inv[j][k] -= m2*M_inv[i][k];
                }
            }
        }
    }
    free(pivrow);
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
