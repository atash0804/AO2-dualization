#ifndef DUALIZATION_MATRIX_UTILS
#define DUALIZATION_MATRIX_UTILS

#include <random>
#include <fstream>
#include <iostream>
#include <cerrno>
#include <vector>
#include <cstring>
#include <bitset>

#include "AO2_trajectory.h"

using namespace std;

#ifndef DUALIZATION_CHUNK_SIZE
#define DUALIZATION_CHUNK_SIZE
enum {
    CHUNK_SIZE = sizeof(ull) * 8,
};
#endif

void generate_matrix(size_t height, size_t width, const string &filename,
                     double density = 0.5, int seed = -1) {
    if (seed != -1) srand(seed);
    ofstream out;
    out.open(filename);
    //cout << "Opened file:" << out.is_open() << endl;
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            out << 1 - int(rand() > RAND_MAX * density) << ' ';
        }
        out << endl;
    }
    out.close();
}

ull** read_matrix(const string &filename, size_t height, size_t width) {
    ifstream in;
    in.open(filename);
    if (!in.is_open()) {
        cerr << "Error: " << strerror(errno) << endl;
        cerr << "Failed to open input file" << endl;
        throw bad_exception();
    }
    ull bit;
    size_t chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

    ull** matrix = new ull*[height];
    for (size_t i = 0; i < height; i++) {
        ull* row = new ull[chunks];
        for (size_t j = 0; j < chunks; j++) {
            ull chunk = 0;
            for (size_t k = 1;
                 (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                if (!(in >> bit)) {
                    cerr << "Incorrect size of input matrix" << endl;
                    throw length_error("");
                }
                if (bit > 1) {
                    cerr << "Incorrect value encountered in input stream: "
                         << bit << endl;
                    throw out_of_range("");
                }
                chunk |= (bit << (CHUNK_SIZE - k));
            }
            row[j] = chunk;
        }
        matrix[i] = row;
    }
    in.close();
    return matrix;
}

void print_matrix(ostream &os, ull** bm,
                    size_t height, size_t width) {
    size_t chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

    cout << chunks << ' ' << CHUNK_SIZE << endl;
    os << "Bit matrix with shape (" << height << ", " << width << ")"
         << endl << endl;
    for (size_t i = 0; i < height; i++) {
        if (!bm[i]) {
            os << "NULL" << endl;
            continue;
        }
        for (size_t j = 0; j < chunks - 1; j++) {
            os << bitset<CHUNK_SIZE>(bm[i][j]);
        }
        size_t bits_left = width-(chunks-1)*CHUNK_SIZE;
        for (size_t j = 0; j < bits_left; j++) {
            os << ((bm[i][chunks-1] >> (CHUNK_SIZE - j - 1)) & ull(1));
        }
        os << endl;
    }
    os << endl << flush;
}

#endif
