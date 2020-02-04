#include <iostream>
#include <string>
#include <map>
#include <set>
#include <ctime>

#include "matrix_utils.h"
#include "runc_m.h"

using namespace std;

int main(int argc, char *argv[]) {
    const string filename = "matrix.txt";
    size_t width = 30;
    size_t height = 26;

    // generate_matrix(height, width, filename, 0.1);
    PartialBitMatrix matr = PartialBitMatrix("matrix.txt", height, width);
    cout << "MATRIX" << matr << endl;

    set<customset> cov;
    map<size_t, set<size_t>> supporting_rows;

    cout << "EXP" << endl << flush;
    clock_t start = clock();
    dualization(matr, supporting_rows, true, true, cov);
    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    cout << "Dualization time: " << elapsed << endl;
    cout << "N coverages: " << cov.size() << endl;
    for (auto coverage: cov){
        for (uint64_t j = 0; j < width; j++) {
            if (coverage.in(j)) {
                cout << j << ' ';
            }
        }
        cout << endl;
    }
}
