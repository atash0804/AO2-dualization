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
*   This file contains experimental implementation of AO2
*   algorithm
************************************************************/


class AO2Trajectory_lightest: public AO2Trajectory {
public:
    AO2Trajectory_lightest(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {}

    //! Find non-zero element of B with the least index
    Coord find_the_least() {
        c_int least_d1 = -1;
        c_int least_d2 = -1;
        for (c_int j = 0; j < col_chunks; j++) {
            for (c_int i = 0; i < n; i++) {
                if (!B[i] || !B[i][j]) continue;
                if (B[i][j] && (least_d2 > (j+1)*CH_SIZE - c_int(log2(B[i][j])) - 1)) {
                    least_d1 = i;
                    least_d2 = (j+1)*CH_SIZE - c_int(log2(B[i][j])) - 1;
                }
            }
            if (least_d1 != c_int(-1)) return Coord(least_d1, least_d2);
        }
        return Coord(least_d1, least_d2);
    }
};


void AO2_lightest_row (c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_lightest traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            // for (auto elem: traj.Q) {
            //     std::cout << "(" << elem.first << "," << elem.second << ")" << std::flush;
            // }
            // std::cout << '\n';
            coverages.push_back(traj.get_coverage());
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    n_cov += coverages.size();
}

int main(int argc, char *argv[]) {
    c_int n = 20;
    c_int m = 20;
    double elapsed = 0;
    uint64_t n_cov = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;

    // generate_matrix(n, m, "matrix.txt", 0.5);
    ull** R = read_matrix("matrix.txt", n, m);

    clock_t start = clock();

    AO2_lightest_row(n, m, R, n_cov, n_extra, n_steps);

    clock_t stop = clock();
    elapsed += (double) (stop - start) / CLOCKS_PER_SEC;

    for (c_int i = 0; i < n; i++) {
        delete [] R[i];
    }
    delete [] R;

    std::cout << elapsed << " & ";
    std::cout << uint64_t(n_cov) << " & ";
    std::cout << uint64_t(n_extra) << " & ";
    std::cout << uint64_t(n_steps) << " \\\\ \n";
}
