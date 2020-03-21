#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>

#include "matrix_utils.h"
#include "AO2.h"

/************************************************************
*   This file contains experimental implementation of AO2
*   algorithm
************************************************************/

void print_stats(std::string name, double elapsed, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    std::cout << name << ": ";
    std::cout << elapsed << " & ";
    std::cout << uint64_t(n_cov) << " & ";
    std::cout << uint64_t(n_extra) << " & ";
    std::cout << uint64_t(n_steps) << " \\\\ \n";
}

int main(int argc, char *argv[]) {
    double SPARSITY = 0.5;
    clock_t start, stop;
    srand(time(NULL));
    for (coord HEIGHT: std::vector<int>{10, 20, 30}) {
        for (coord WIDTH: std::vector<int>{30, 50, 70}) {
            double elapsed1 = 0, elapsed2 = 0, elapsed3 = 0;
            uint64_t n_cov1 = 0, n_cov2 = 0, n_cov3 = 0;
            uint64_t n_extra1 = 0, n_extra2 = 0, n_extra3 = 0;
            uint64_t n_steps1 = 0, n_steps2 = 0, n_steps3 = 0;
            for (int i = 0; i < 10; i++) {
                generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
                ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

                if (has_zero_rows(R, HEIGHT, WIDTH)) {
                    std::cout << "Matrix contains zero rows \n" << std::endl;
                    i--;
                    continue;
                }

                start = clock();
                AO2(HEIGHT, WIDTH, R, n_cov1, n_extra1, n_steps1);
                stop = clock();
                elapsed1 += (double) (stop - start) / CLOCKS_PER_SEC;

                start = clock();
                AO2Moptimized(HEIGHT, WIDTH, R, n_cov2, n_extra2, n_steps2);
                stop = clock();
                elapsed2 += (double) (stop - start) / CLOCKS_PER_SEC;

                start = clock();
                AO2Best(HEIGHT, WIDTH, R, n_cov3, n_extra3, n_steps3);
                stop = clock();
                elapsed3 += (double) (stop - start) / CLOCKS_PER_SEC;

                for (coord i = 0; i < HEIGHT; i++) {
                    delete [] R[i];
                }
                delete [] R;
            }
            std::cout << HEIGHT << " \\times " << WIDTH << " & ";
            // std::cout << elapsed1 / 10 << " & ";
            std::cout << elapsed2 / 10 << " & ";
            std::cout << elapsed3 / 10 << " & ";
            std::cout << uint64_t(n_cov2 / 10) << " & ";
            // std::cout << uint64_t(n_extra1 / 10) << " & ";
            std::cout << uint64_t(n_extra2 / 10) << " & ";
            std::cout << uint64_t(n_extra3 / 10) << " & ";
            // std::cout << uint64_t(n_steps1 / 10) << " & ";
            std::cout << uint64_t(n_steps2 / 10) << " & ";
            std::cout << uint64_t(n_steps3 / 10) << " \\\\ \n";
        }
    }
    return 0;
}
