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

class AO2TrajectoryExp: public AO2MoptimizedTrajectory {
protected:

    //! Find non-zero element of B with the least index
    virtual Element find_the_least() {
        if (latest_element.first != coord(-1)) {
            for (coord i = latest_element.first + 1; i < n; i++) {
                if (!B[i]) continue;
                if (B[i][latest_element.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - latest_element.second % CH_SIZE))) {
                    latest_element = Element(-1, latest_element.second);
                    return Element(i, latest_element.second);
                }
            }
            latest_element = Element(-1, latest_element.second);
        }
        coord least_d1 = -1;
        coord least_d2 = -1;
        coord least_Er = n*m+1, curr_Er;
        coord *columns = new coord[m]();
        coord best_col = -1;
        bool is_set;

        for (coord i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (coord j = 0; j < m; j++) {
                if (B[i][j/CH_SIZE] & (ull(1) << (CH_SIZE_1 - j % CH_SIZE))) columns[j]++;
            }
        }

        for (coord i = 0; i < n; i++) {
            // if row is not covered or is competing
            if (!(states.top()[i] & ST_IS_COV) || states.top()[i] & ST_IS_COMP) {
                curr_Er = 0;
                is_set = false;
                for (coord j = 0; j < m; j++) {
                    if (M[i][j/CH_SIZE] & (ull(1) << (CH_SIZE_1 - j % CH_SIZE))) {
                        curr_Er += columns[j];
                        if (columns[j] && (!is_set || (columns[best_col] > columns[j]))) {
                            best_col = j;
                            is_set = true;
                        }
                    }
                }
                // cout << curr_Er << "*" << min_weight_col << " \n";
                if (curr_Er < least_Er) {
                    least_Er = curr_Er;
                    least_d2 = best_col;
                    for (coord i1 = 0; i1 < n; i1++) {
                        if (!B[i1]) continue;
                        if (B[i1][least_d2/CH_SIZE] && (B[i1][least_d2/CH_SIZE]&M[i][least_d2/CH_SIZE]) && (B[i1][least_d2/CH_SIZE] & (ull(1) << (CH_SIZE_1 - least_d2 % CH_SIZE)))) {
                            least_d1 = i1;
                            break;
                        }
                    }
                }
            }
        }
        // cout << endl;
        delete [] columns;
        return Element(least_d1, least_d2);
    }

public:
    AO2TrajectoryExp(coord d1, coord d2, ull** Matrix) :
        AO2MoptimizedTrajectory(d1, d2, Matrix) {}
};

void AO2Exp(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2TrajectoryExp traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    n_cov += coverages.size();
}



void print_stats(std::string name, double elapsed, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    std::cout << name << ": ";
    std::cout << elapsed << " & ";
    std::cout << uint64_t(n_cov) << " & ";
    std::cout << uint64_t(n_extra) << " & ";
    std::cout << uint64_t(n_steps) << " \\\\ \n";
}


int main(int argc, char *argv[]) {
    coord HEIGHT = 20;
    coord WIDTH = 50;
    double SPARSITY = 0.5;

    double elapsed = 0;
    uint64_t n_cov = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;
    clock_t start, stop;

    srand(time(NULL));
    int i = 0;
    while (true) {
        if (i > 10) break;
        cout << i++ << '\n';
        generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
        ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

        if (has_zero_rows(R, HEIGHT, WIDTH)) {
            std::cout << "Matrix contains zero rows \n" << std::endl;
            continue;
        }

        start = clock();
        n_cov = n_extra = n_steps = 0;
        AO2(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);
        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        print_stats("AO2  :", elapsed, n_cov, n_extra, n_steps);

        start = clock();
        n_cov = n_extra = n_steps = 0;
        AO2Moptimized(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);
        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        print_stats("AO2M+  :", elapsed, n_cov, n_extra, n_steps);

        start = clock();
        n_cov = n_extra = n_steps = 0;
        AO2Best(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);
        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        print_stats("Best :", elapsed, n_cov, n_extra, n_steps);

        for (coord i = 0; i < HEIGHT; i++) {
            delete [] R[i];
        }
        delete [] R;
        // if ((n_cov1 == n_cov2)) continue;
        // break;
    }
}
