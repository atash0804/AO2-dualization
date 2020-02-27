#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

#include "matrix_utils.h"
#include "AO2_trajectory.h"
#include "AO2.h"

/************************************************************
*   This file contains basic implementation of AO2
*   algorithm
************************************************************/

int main(int argc, char *argv[]) {
    c_int HEIGHT = 20;
    c_int WIDTH = 20;
    double SPARSITY = 0.5;

    srand(time(NULL));
    double elapsed = 0;
    uint64_t n_cov = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;
    // generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
    ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

    if (has_zero_rows(R, HEIGHT, WIDTH)) {
        std::cout << "Matrix contains zero rows \n" << std::endl;
        return 1;
    }

    clock_t start = clock();

    AO2(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);

    clock_t stop = clock();
    elapsed += (double) (stop - start) / CLOCKS_PER_SEC;

    for (c_int i = 0; i < HEIGHT; i++) {
        delete [] R[i];
    }
    delete [] R;

    std::cout << elapsed << " & ";
    std::cout << uint64_t(n_cov) << " & ";
    std::cout << uint64_t(n_extra) << " & ";
    std::cout << uint64_t(n_steps) << " \\\\ \n";
}
