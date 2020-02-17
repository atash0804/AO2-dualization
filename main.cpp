#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

#include "matrix_utils.h"
#include "AO2_trajectory.h"

/************************************************************
*   This file contains basic implementation of AO2
*   algorithm
************************************************************/

//! Implementation of basic AO2 algorithm
/*! Cycles through all possible trajectories and collects coverages.*/
void AO2(c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    n_cov += coverages.size();
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    c_int HEIGHT = 20;
    c_int WIDTH = 20;
    double SPARSITY = 0.5;

    double elapsed = 0;
    uint64_t n_cov = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;

    // generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
    ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

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
