#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

#include "matrix_utils.h"
#include "AO2.h"

/************************************************************
*   This file contains experimental implementation of AO2
*   algorithm
************************************************************/

class AO2ZeroTrajectory_exp: public AO2ZeroTrajectory {
public:
    AO2ZeroTrajectory_exp(coord d1, coord d2, ull** Matrix) :
        AO2ZeroTrajectory(d1, d2, Matrix) {}
protected:
    virtual Element find_the_least() {
        coord least_d1 = -1;
        coord least_d2 = -1;
        coord least_Er = n*m+1, curr_Er;
        coord *columns = new coord[m]();
        bool found;

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
                for (coord j = 0; j < m; j++) {
                    if (M[i][j/CH_SIZE] & (ull(1) << (CH_SIZE_1 - j % CH_SIZE))) curr_Er += columns[j];
                }
                // cout << curr_Er << "***********\n";
                if (curr_Er < least_Er) {
                    least_d1 = -1;
                    least_d2 = -1;
                    for (coord i1 = 0; i1 < n; i1++) {
                        if (!B[i1]) continue;
                        for (coord j1 = 0; j1 < col_chunks; j1++) {
                            // cout << bitset<32>(B[i1][j1]) << ' ' << bitset<32>(B[i1][j1]&M[i][j1]) << '\n';
                            // cout << least_d2 << ' ' << (j1+1)*CH_SIZE - coord(log2(B[i1][j1])) - 1;
                            // cout << "_________\n";
                            if (B[i1][j1] && (B[i1][j1]&M[i][j1]) && (least_d2 > (j1+1)*CH_SIZE - coord(log2(B[i1][j1]&M[i][j1])) - 1)) {

                                least_d1 = i1;
                                least_d2 = (j1+1)*CH_SIZE - coord(log2(B[i1][j1]&M[i][j1])) - 1;
                                least_Er = curr_Er;
                                // cout << "NEW " << least_d1 << ' ' << least_d2 << endl;
                            }
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
    bool complete_trajectory() {
        while (!check_empty()) {
            if (!check_covers_comp_rows()) return false;
            eliminate_dominating_rows();
            if (!check_covers_comp_rows() || check_empty()) {
                delete [] states.top();
                states.pop();
                return false;
            }
            Element candidate = find_the_least();
            Q.push_back(candidate);
            if (!eliminate_incompatible(candidate)) {
                return false;
            }
        }
        return true;
    }
};


void AO2Zero_exp(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2ZeroTrajectory_exp traj(n, m, R);
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


class AO2ZeroTrajectory_exp2: public AO2ZeroTrajectory {
public:
    AO2ZeroTrajectory_exp2(coord d1, coord d2, ull** Matrix) :
        AO2ZeroTrajectory(d1, d2, Matrix) {}
protected:
    virtual Element find_the_least() {
        coord least_d1 = -1;
        coord least_d2 = -1;
        coord least_Er = 0, curr_Er;
        coord *columns = new coord[m]();
        bool found;

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
                for (coord j = 0; j < m; j++) {
                    if (M[i][j/CH_SIZE] & (ull(1) << (CH_SIZE_1 - j % CH_SIZE))) curr_Er += columns[j];
                }
                // cout << curr_Er << "***********\n";
                if (curr_Er > least_Er) {
                    least_d1 = -1;
                    least_d2 = -1;
                    for (coord i1 = 0; i1 < n; i1++) {
                        if (!B[i1]) continue;
                        for (coord j1 = 0; j1 < col_chunks; j1++) {
                            // cout << bitset<32>(B[i1][j1]) << ' ' << bitset<32>(B[i1][j1]&M[i][j1]) << '\n';
                            // cout << least_d2 << ' ' << (j1+1)*CH_SIZE - coord(log2(B[i1][j1])) - 1;
                            // cout << "_________\n";
                            if (B[i1][j1] && (B[i1][j1]&M[i][j1]) && (least_d2 > (j1+1)*CH_SIZE - coord(log2(B[i1][j1]&M[i][j1])) - 1)) {

                                least_d1 = i1;
                                least_d2 = (j1+1)*CH_SIZE - coord(log2(B[i1][j1]&M[i][j1])) - 1;
                                least_Er = curr_Er;
                                // cout << "NEW " << least_d1 << ' ' << least_d2 << endl;
                            }
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
    bool complete_trajectory() {
        // std::cout << "INITIAL:\n";
        // print_B();
        while (!check_empty()) {
            if (!check_covers_comp_rows()) return false;
            eliminate_dominating_rows();
            // std::cout << "ELIM DM ROWS:\n";
            // print_B();
            if (!check_covers_comp_rows() || check_empty()) {
                delete [] states.top();
                states.pop();
                return false;
            }
            Element candidate = find_the_least();
            Q.push_back(candidate);
            if (!eliminate_incompatible(candidate)) {
                return false;
            }
            // std::cout << "ELIM INCOMPAT:" << candidate.first << ' ' << candidate.second << "\n";
            // print_B();
        }
        for (coord i = 0; i < n; i++) {
            if (!(states.top()[i] & ST_IS_COV)) {
                cout << "NA";
                return false;
            }
        }
        return true;
    }
};


void AO2Zero_exp2(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2ZeroTrajectory_exp2 traj(n, m, R);
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

int main(int argc, char *argv[]) {
    coord HEIGHT = 20;
    coord WIDTH = 20;
    double SPARSITY = 0.5;

    double elapsed = 0;
    uint64_t n_cov1 = 0, n_cov2 = 0, n_cov3 = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;

    srand(time(NULL));
    int i = 0;
    while (true) {
        cout << i++ << '\n';
        generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
        ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

        if (has_zero_rows(R, HEIGHT, WIDTH)) {
            std::cout << "Matrix contains zero rows \n" << std::endl;
            continue;
        }

        clock_t start = clock();

        n_cov1 = n_extra = n_steps = 0;
        AO2Zero_exp(HEIGHT, WIDTH, R, n_cov1, n_extra, n_steps);

        clock_t stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        std::cout << elapsed << " & ";
        std::cout << uint64_t(n_cov1) << " & ";
        std::cout << uint64_t(n_extra) << " & ";
        std::cout << uint64_t(n_steps) << " \\\\ \n";

        start = clock();

        n_cov2 = n_extra = n_steps = 0;
        AO2Zero_exp2(HEIGHT, WIDTH, R, n_cov2, n_extra, n_steps);

        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        std::cout << elapsed << " & ";
        std::cout << uint64_t(n_cov2) << " & ";
        std::cout << uint64_t(n_extra) << " & ";
        std::cout << uint64_t(n_steps) << " \\\\ \n";

        start = clock();

        n_cov3 = n_extra = n_steps = 0;
        AO2Zero(HEIGHT, WIDTH, R, n_cov3, n_extra, n_steps);

        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        std::cout << elapsed << " & ";
        std::cout << uint64_t(n_cov3) << " & ";
        std::cout << uint64_t(n_extra) << " & ";
        std::cout << uint64_t(n_steps) << " \\\\ \n";
        for (coord i = 0; i < HEIGHT; i++) {
            delete [] R[i];
        }
        delete [] R;
        if ((n_cov1 == n_cov2) && (n_cov2 == n_cov3)) continue;
        break;
    }
}
