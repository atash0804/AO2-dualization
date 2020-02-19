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

//! Class to represent changing trajectories
/*!
  This class implements the AO2 trajectory with an improvement:
    similar columns are deleted in constructor.
*/
class AO2Trajectory_delete_similar_cols: public AO2Trajectory {
protected:
    std::vector<ull*> similar_cols;
    std::vector<c_int> similar_cols_counter;
public:
    AO2Trajectory_delete_similar_cols(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {
        ull* mask;
        ull aij;
        c_int counter, min_pos;
        for (c_int col = 0; col < m; col++) {
            mask = new ull[col_chunks]();
            for (c_int i = 0; i < n; i++) {
                aij = Matrix[i][col / CH_SIZE] & (ull(1) << (CH_SIZE_1 - col % CH_SIZE)) ? ull(-1) : ull(0);
                for (c_int j = 0; j < col_chunks; j++) {
                    mask[j] |= aij ^ Matrix[i][j];
                }
            }
            mask[row_chunks - 1] |= ull(-1) >> (m % CH_SIZE);
            mask[col / CH_SIZE] |= ull(1) << (CH_SIZE_1 - col % CH_SIZE);

            // std::cout << std::bitset<32>(~mask[0]) << std::endl;
            min_pos = -1;
            bool has_similar = false;
            for (c_int j = 0; j < col_chunks; j++) {
                if (~mask[j]) {
                    min_pos = (j + 1)* CH_SIZE - c_int(log2(~mask[j])) - 1;
                    has_similar = true;
                    break;
                }
            }
            counter = 1;
            if (has_similar) {
                if (min_pos > col) {
                    for (c_int i = 0; i < n; i++) {
                        for (c_int j = 0; j < col_chunks; j++) {
                            B[i][j] &= mask[j];
                            M[i][j] &= mask[j];
                        }
                    }
                }
                for (c_int j = 0; j < col_chunks; j++) {
                    aij = ~mask[j];
                    while (aij) {
                        aij &= (aij - 1);
                        counter++;
                    }
                }
                // if (counter > 1) std::cout << "DELETED COLS" << counter << "\n";
            }
            similar_cols.push_back(mask);
            similar_cols_counter.push_back(counter);
        }
    }

    ~AO2Trajectory_delete_similar_cols() {
        for (c_int col = 0; col < m; col++) {
            delete [] similar_cols[col];
        }
    }

    uint64_t get_n_coverages() {
        uint64_t result = 1;
        for (Coord item: Q) {
            result *= similar_cols_counter[item.second];
            // if (similar_cols_counter[item.second] > 1) std::cout << result << ' ' << item.second << ' ' << std::bitset<32>(similar_cols[item.second][0]) << std::bitset<32>(similar_cols[item.second][1]) << std::endl;
        }
        return result;
    }
};

void AO2_delete_similar_cols (c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_delete_similar_cols traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
            n_cov += traj.get_n_coverages();
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
}

class AO2Trajectory_lightest_row: public AO2Trajectory {
public:
    AO2Trajectory_lightest_row(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {}

    void complete_trajectory() {
        while (!check_empty()) {
            // std::cout << "INITIAL:\n";
            // print_B();
            eliminate_dominating_rows();
            // std::cout << "ELIM DM ROWS:\n";
            // print_B();
            if (check_empty()) {
                // std::cout << "EMPTY AFTER DOM ROWS" << '\n';
                return;
            }
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);

            // std::cout << "ELIM INCOMPAT:" << candidate.first << ' ' << candidate.second << "\n";
            // print_B();
            // if (!check_coverage()) {
            //     std::cout << "NOT A COV AFTER ELIM INCOMPAT" << '\n';
            //     find_neighbour();
            // }
        }
        // std::cout << "NOT COV ROWS" << std::bitset<32>(not_cov_rows.top()[0]) << '\n';

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


class AO2Trajectory_heaviest_col: public AO2Trajectory {
public:
    AO2Trajectory_heaviest_col(c_int d1, c_int d2, ull** Matrix) :
        AO2Trajectory(d1, d2, Matrix) {}

    //! Find non-zero element of B in the heaviest col
    Coord find_the_least() {
        c_int least_d1 = -1;
        c_int least_d2 = -1;
        c_int best_weight = 0;

        ull *weights = new ull[m]();
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (c_int col = 0; col < m; col++) {
                weights[col] += ull(1) & (B[i][col/CH_SIZE] >> (CH_SIZE_1 - col % CH_SIZE));
            }
        }

        for (c_int col = 0; col < m; col++) {
            if (weights[col] > best_weight) {
                best_weight = weights[col];
                for (c_int i = 0; i < n; i++) {
                    if (!B[i]) continue;
                    if (B[i][col/CH_SIZE] & (ull(1) << (CH_SIZE_1 - col % CH_SIZE))) {
                        least_d1 = i;
                        break;
                    }
                }
                least_d2 = col;
            }
        }
        // print_B();
        // std::cout << "EL: " << least_d1 << ' ' << least_d2 << '\n' << std::flush;
        return Coord(least_d1, least_d2);
    }
};

void AO2_heaviest_col (c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_heaviest_col traj(n, m, R);
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

int main(int argc, char *argv[]) {
    c_int HEIGHT = 10;
    c_int WIDTH = 10;
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

    AO2_heaviest_col(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);

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
