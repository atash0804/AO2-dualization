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
    size_t width = 10;
    size_t height = 5;

    generate_matrix(height, width, filename, 0.5, 35);
    PartialBitMatrix matr = PartialBitMatrix("matrix.txt", height, width);
    // cout << "MATRIX" << matr << endl;

    set<customset> cov;
    map<size_t, set<size_t>> supporting_rows;

    // cout << "EXP" << endl << flush;
    clock_t start = clock();
    dualization(matr, supporting_rows, true, true, cov);
    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    cout << "Dualization time: " << elapsed << endl;
    cout << "N coverages: " << cov.size() << endl;
    set<customset> should_delete;
    for (auto cov1: cov){
        for (auto cov2: cov){
            bool is_subset1 = true;
            bool is_subset2 = true;
            for (size_t j = 0; j < cov1.chunks; j++) {
                if ((cov1.data[j] & cov2.data[j]) != cov1.data[j]) is_subset1 = false;
                if ((cov1.data[j] & cov2.data[j]) != cov2.data[j]) is_subset2 = false;
            }
            if (is_subset1 && is_subset2) continue;
            if (is_subset1) should_delete.insert(cov2);
            if (is_subset2) should_delete.insert(cov1);
        }
    }
    for (auto tmp: should_delete) {
        cov.erase(tmp);
    }
    cout << "N coverages: " << cov.size() << endl;
    // for (auto coverage: cov){
    //     for (uint64_t j = 0; j < width; j++) {
    //         if (coverage.in(j)) {
    //             cout << j << ' ';
    //         }
    //     }
    //     cout << endl;
    // }
}
