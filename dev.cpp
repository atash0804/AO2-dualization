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

class AO2Trajectory_lightest_row: public AO2Trajectory {
public:
    AO2Trajectory_lightest_row(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {}

    void complete_trajectory() {
        while (!check_empty()) {
            std::cout << "INITIAL:\n";
            print_B();
            eliminate_dominating_rows();
            std::cout << "ELIM DM ROWS:\n";
            print_B();
            if (check_empty()) {
                // std::cout << "EMPTY AFTER DOM ROWS" << '\n';
                return;
            }
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);

            std::cout << "ELIM INCOMPAT:" << candidate.first << ' ' << candidate.second << "\n";
            print_B();
            // if (!check_coverage()) {
                // std::cout << "NOT A COV AFTER ELIM INCOMPAT" << '\n';
                // find_neighbour();
            // }
        }
        std::cout << "NOT COV ROWS" << std::bitset<32>(not_cov_rows.top()[0]) << '\n';

    }

    //! Find non-zero element of B in the lightest row
    Coord find_the_least() {
        c_int least_d1 = -1;
        c_int least_d2 = -1;
        c_int least_weight = m+1;
        c_int cur_weight;
        ull tmp;

        ull *available = new ull[m]();
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (c_int col = 0; col < m; col++) {
                available[col] += ull(1) & (B[i][col/CH_SIZE] >> (CH_SIZE_1 - col % CH_SIZE));
            }
        }

        // print_B();
        // for (c_int col = 0; col < m; col++) {
        //     std::cout << available[col] << " ";
        // }
        // std::cout << "+++++++++++++++\n";

        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            cur_weight = 0;
            for (c_int j = 0; j < col_chunks; j++) {
                tmp = B[i][j];
                while (tmp) {
                    tmp &= (tmp - 1);
                    cur_weight++;
                }
            }
            if (cur_weight < least_weight) {
                for (c_int col = 0; col < m; col++) {
                    if (B[i][col/CH_SIZE] & (ull(1) << (CH_SIZE_1 - col % CH_SIZE)) && (available[col] > 1)) {
                        least_weight = cur_weight;
                        least_d1 = i;
                        least_d2 = col;
                    }
                }
            }
        }
        if (least_d1 != c_int(-1) && least_d2 != c_int(-1)) {
            return Coord(least_d1, least_d2);
        }
        else {
            least_d1 = -1;
            least_d2 = -1;
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
    }
};

void AO2_lightest_row (c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_lightest_row traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            for (auto q: traj.get_coverage()) std::cout << q << ' ';
            std::cout << '\n';
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    n_cov += coverages.size();
}

class AO2Trajectory_stop_not_upper: public AO2Trajectory {
public:
    AO2Trajectory_stop_not_upper(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {}

    bool check_covers_comp_rows() {
        // std::cout << "INTO COMPROWS\n" << std::flush;
        ull *mask = new ull[col_chunks]();
        for (Coord item: Q) {
            mask[item.second/CH_SIZE] |= ull(1) << (CH_SIZE_1 - item.second % CH_SIZE);
        }

        bool is_empty, valid;
        for (c_int i = 0; i < n; i++) {
            // if row is not covered
            if (not_cov_rows.top()[i/CH_SIZE] & (ull(1) << (CH_SIZE_1-i%CH_SIZE))) {
                // std::cout << "1st BR\n" << std::flush;
                is_empty = true;
                for (c_int i1 = 0; i1 < n; i1++) {
                    if (!B[i1]) continue;
                    for (c_int j = 0; j < col_chunks; j++) {
                        if (B[i1][j] & M[i][j]) {
                            is_empty = false;
                            break;
                        }
                    }
                    if (!is_empty) break;
                }

                if (is_empty) {
                    // std::cout << "OUT COMPROWS FAIL1\n" << std::flush;
                    return false;
                }
            } else {
                // std::cout << "2nd BR\n" << std::flush;
                if (deleted_by_domination.top()[i/CH_SIZE] & ull(1) << (CH_SIZE_1 - i % CH_SIZE)) continue;
                // check competing rows amongst those not deleted by domination
                for (Coord item: Q) {
                    if (i >= item.first) continue;
                    valid = true;
                    for (c_int j = 0; j < col_chunks; j++) {
                        if ((M[i][j] & mask[j]) != (M[item.first][j] & mask[j])) {
                            valid = false;
                            break;
                        }
                    }
                    // valid means that row i competes with row item.first
                    if (valid) {
                        is_empty = true;
                        for (c_int i1 = 0; i1 < n; i1++) {
                            if (!B[i1]) continue;
                            for (c_int j = 0; j < col_chunks; j++) {
                                if (B[i1][j] & M[i][j]) {
                                    is_empty = false;
                                    break;
                                }
                            }
                            if (!is_empty) break;
                        }
                        if (is_empty) {
                            // std::cout << "OUT COMPROWS FAIL2\n" << std::flush;
                            return false;
                        }
                        break;
                    }
                }
            }
        }
        // std::cout << "OUT COMPROWS OK\n" << std::flush;
        return true;
    }

    void complete_trajectory() {
        // std::cout << "INITIAL:\n";
        // print_B();
        while (!check_empty() && check_covers_comp_rows()) {
            eliminate_dominating_rows();
            // std::cout << "ELIM DM ROWS:\n";
            // print_B();
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);
            // std::cout << "ELIM INCOMPAT:" << candidate.first << ' ' << candidate.second << "\n";
            // print_B();
        }
    }
};

void AO2_stop_not_upper(c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_stop_not_upper traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    n_cov += coverages.size();
}

int main(int argc, char *argv[]) {
    c_int HEIGHT = 20;
    c_int WIDTH = 20;
    double SPARSITY = 0.5;

    double elapsed = 0;
    uint64_t n_cov = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;

    srand(time(NULL));
    generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
    ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

    if (has_zero_rows(R, HEIGHT, WIDTH)) {
        std::cout << "Matrix contains zero rows \n" << std::endl;
        return 1;
    }

    clock_t start = clock();

    AO2_stop_not_upper(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);

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
