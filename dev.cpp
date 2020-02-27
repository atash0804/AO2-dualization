#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

#include "matrix_utils.h"

/************************************************************
*   This file contains experimental implementation of AO2
*   algorithm
************************************************************/

//! Class to represent changing trajectories
/*!
  This class implements the trajectory - sequence of changing pairs (B, Q),
  which ends with a candidate for shortest coverage and all operations on
  trajectories which are needed for AO2.

*/
class AO2ZeroTrajectory {
protected:
    coord n; //!< First dimension of matrix
    coord m; //!< Second dimension of matrix
    coord col_chunks;
    coord row_chunks;

    //! Stores initial state of matrix
    ull** M;

    //! Representation of set B (elements that can be used to complete the trajectory)
    /*! B is represented by binary matrix of elements that can still be used
    to complete the trajectory*/
    ull** B;

    //! Stack of changes to B
    /*! Each change is represented as binary matrix of size mxn*/
    BMatrixStack changes;

    //! Representation of set Q (elements that from the trajectory)
    /*! Vector of pairs (i_r, j_r) which form the Q set*/
    QVector Q;

    //! Stack to represent rows not covered by Q at each point
    /*! Each element of this stack is set of row numbers, which correspond to
    rows not covered by Q*/
    RowSetStack not_cov_rows;

    RowSetStack deleted_by_domination;

    //! Checks if B matrix is empty (Stopping criterion)
    bool check_empty() {
        for (coord i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (coord j = 0; j < col_chunks; j++) {
                if (B[i][j]) return false;
            }
        }
        return true;
    }

    bool check_covers_comp_rows() {
        // std::cout << "INTO COMPROWS\n" << std::flush;
        ull *mask = new ull[col_chunks]();
        for (Element item: Q) {
            mask[item.second/CH_SIZE] |= ull(1) << (CH_SIZE_1 - item.second % CH_SIZE);
        }

        bool is_empty, valid;
        for (coord i = 0; i < n; i++) {
            // if row is not covered
            if (not_cov_rows.top()[i/CH_SIZE] & (ull(1) << (CH_SIZE_1-i%CH_SIZE))) {
                // std::cout << "1st BR\n" << std::flush;
                is_empty = true;
                for (coord i1 = 0; i1 < n; i1++) {
                    if (!B[i1]) continue;
                    for (coord j = 0; j < col_chunks; j++) {
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
                for (Element item: Q) {
                    if (i >= item.first) continue;
                    valid = true;
                    for (coord j = 0; j < col_chunks; j++) {
                        if ((M[i][j] & mask[j]) != (M[item.first][j] & mask[j])) {
                            valid = false;
                            break;
                        }
                    }
                    // valid means that row i competes with row item.first
                    if (valid) {
                        is_empty = true;
                        for (coord i1 = 0; i1 < n; i1++) {
                            if (!B[i1]) continue;
                            for (coord j = 0; j < col_chunks; j++) {
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


    //! Create delta_star set
    /*! Includes 1s from dominating rows which correspond to 0s in dominated row.*/
    void eliminate_dominating_rows() {
        ull **updated_B = new ull*[n];
        bool i1_dom_i2, i2_dom_i1;
        for (coord i1 = 0; i1 < n; i1++) {
            if (!B[i1]) {
                updated_B[i1] = NULL;
                continue;
            } else {
                updated_B[i1] = new ull[col_chunks]();
                for (coord j = 0; j < col_chunks; j++) updated_B[i1][j] = B[i1][j];
            }
        }
        ull *H_B = new ull[col_chunks]();
        for (coord i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (coord j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }

        ull* new_del = new ull[row_chunks];
        for (coord i = 0; i < row_chunks; i++) {
            new_del[i] = deleted_by_domination.top()[i];
        }

        for (coord i1 = 0; i1 < n; i1++) {
            if (!B[i1] || (!(not_cov_rows.top()[i1/CH_SIZE] & (ull(1) << (CH_SIZE_1-i1%CH_SIZE))))) continue;
            for (coord i2 = i1+1; i2 < n; i2++) {
                if (!B[i2] || (!(not_cov_rows.top()[i2/CH_SIZE] & (ull(1) << (CH_SIZE_1-i2%CH_SIZE))))) continue;
                i1_dom_i2 = i2_dom_i1 = true;
                for (coord j = 0; j < col_chunks; j++) {
                    if ((M[i2][j] & M[i1][j] & H_B[j]) != (M[i1][j] & H_B[j])) i2_dom_i1 = false;
                    if ((M[i1][j] & M[i2][j] & H_B[j]) != (M[i2][j] & H_B[j])) i1_dom_i2 = false;
                    if (!(i2_dom_i1 || i1_dom_i2)) break;
                }
                if (!(i2_dom_i1 || i1_dom_i2)) continue;
                if (i2_dom_i1) {
                    new_del[i2/CH_SIZE] ^= (~new_del[i2/CH_SIZE]) & ( ull(1) << (CH_SIZE_1 - i2 % CH_SIZE));
                    delete [] updated_B[i2];
                    updated_B[i2] = NULL;
                } else {
                    if (i1_dom_i2) {
                        new_del[i1/CH_SIZE] ^= (~new_del[i1/CH_SIZE]) & ( ull(1) << (CH_SIZE_1 - i1 % CH_SIZE));
                        delete [] updated_B[i1];
                        updated_B[i1] = NULL;
                    }
                }
            }
        }
        delete [] H_B;
        for (coord i = 0; i < n; i++) delete [] B[i];
        delete [] B;
        deleted_by_domination.push(new_del);
        B = updated_B;
        return;
    }

    //! Find non-zero element of B with the least index
    virtual Element find_the_least() {
        coord least_d1 = -1;
        coord least_d2 = -1;
        for (coord j = 0; j < col_chunks; j++) {
            for (coord i = 0; i < n; i++) {
                if (!B[i] || !B[i][j]) continue;
                if (B[i][j] && (least_d2 > (j+1)*CH_SIZE - coord(log2(B[i][j])) - 1)) {
                    least_d1 = i;
                    least_d2 = (j+1)*CH_SIZE - coord(log2(B[i][j])) - 1;
                }
            }
            if (least_d1 != coord(-1)) return Element(least_d1, least_d2);
        }
        return Element(least_d1, least_d2);
    }

    //! Create delta_ij set
    /*! Includes elements of B incompatible with element. Suppose the element
    is B_ij. Incompatible elements are all B_kl such that either B_kj or B_il
    equals 1.
    To identify incompatible elements, we cycle throw rows k=1..n and look at
    B_kj. If it equals 1, all 1s of the row are incomatible. Otherwise, if
    B_kj = 0, the incomatible elements are all 1s from B_k & B_i (elementwise
    conjunction of two rows)*/
    void eliminate_incompatible(Element& el) {
        ull **updated_B = new ull*[n];
        ull* new_cov = new ull[row_chunks];
        for (coord i = 0; i < row_chunks; i++) {
            new_cov[i] = not_cov_rows.top()[i];
        }

        for (coord k = 0; k < n; k++) {
            if (!B[k]) {
                if (M[k][el.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - el.second % CH_SIZE))) {
                    new_cov[k/CH_SIZE] ^= new_cov[k/CH_SIZE] & ( ull(1) << (CH_SIZE_1 - k % CH_SIZE));
                }
                updated_B[k] = NULL;
                continue;
            } else {
                if (M[k][el.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - el.second % CH_SIZE))) {
                    new_cov[k/CH_SIZE] ^= new_cov[k/CH_SIZE] & (ull(1) << (CH_SIZE_1 - k % CH_SIZE));
                    updated_B[k] = NULL;
                } else {
                    updated_B[k] = new ull[col_chunks]();
                    for (coord j = 0; j < col_chunks; j++)
                        updated_B[k][j] = B[k][j] ^ (B[k][j] & M[el.first][j]);
                }
            }
        }

        ull* tmp = new ull[row_chunks];
        for (coord i = 0; i < row_chunks; i++) {
            tmp[i] = deleted_by_domination.top()[i];
        }
        deleted_by_domination.push(tmp);
        not_cov_rows.push(new_cov);
        changes.push(B);
        B = updated_B;
        return;
    }

    //! Update B stack and B itself while finding a neighbouring trajectory
    /*! Converts stack from
        [...||delta_star_(t-1)||delta_ij_(t-1)||delta_star_(t)] to
        [...||delta_star_(t-1)||delta_ij_(t-1) + delta_star_(t) + latest_element|...]

        As delta_star_(t) is already applied to B, we should apply only latest_element*/
    void update_stack(Element& el) {
        if (changes.size() > 0) {
            ull** latest_state = changes.top();
            changes.pop();

            delete [] deleted_by_domination.top();
            deleted_by_domination.pop();
            ull* del_in_B = deleted_by_domination.top();
            deleted_by_domination.pop();
            delete [] deleted_by_domination.top();
            deleted_by_domination.pop();
            deleted_by_domination.push(del_in_B);

            for (coord p = 0; p < n; p++) {
                delete [] B[p];
            }
            delete [] B;
            B = latest_state;
            B[el.first][el.second / CH_SIZE] ^= ull(1) << (CH_SIZE_1 - el.second % CH_SIZE);
        } else {
            B[el.first][el.second / CH_SIZE] ^= ull(1) << (CH_SIZE_1 - el.second % CH_SIZE);
        }
        bool is_zero = true;
        for (coord j = 0; j < col_chunks; j++) {
            if (B[el.first][j]) is_zero = false;
        }
        if (is_zero) {
            delete [] B[el.first];
            B[el.first] = NULL;
        }
        return;
    }
public:
    AO2ZeroTrajectory(coord d1, coord d2, ull** Matrix) {
        n = d1;
        m = d2;
        col_chunks = (m-1) / CH_SIZE + 1;
        row_chunks = (n-1) / CH_SIZE + 1;
        B = new ull*[n];
        M = new ull*[n];
        for (coord i = 0; i < n; i++) {
            B[i] = new ull[col_chunks];
            M[i] = new ull[col_chunks];
            for (coord j = 0; j < col_chunks; j++) {
                B[i][j] = M[i][j] = Matrix[i][j];
            }
        }
        ull* tmp = new ull[row_chunks];
        for (coord i = 0; i <  row_chunks; i++) {
            tmp[i] = ull(-1);
        }
        not_cov_rows.push(tmp);
        tmp = new ull[row_chunks];
        for (coord i = 0; i <  row_chunks; i++) {
            tmp[i] = ull(0);
        }
        deleted_by_domination.push(tmp);
    }

    ~AO2ZeroTrajectory() {
        for (coord i = 0; i < n; i++) {
            if (B[i]) delete [] B[i];
            if (M[i]) delete [] M[i];
        }
        delete [] B;
        delete [] M;

        for (; !changes.empty(); changes.pop()) {
            for (coord i = 0; i < n; i++) {
                delete [] changes.top()[i];
            }
            delete [] changes.top();
        }

        for (; !not_cov_rows.empty(); not_cov_rows.pop()) {
            delete [] not_cov_rows.top();
        }

        for (; !deleted_by_domination.empty(); deleted_by_domination.pop()) {
            delete [] deleted_by_domination.top();
        }
    }

    //! Checks that all rows not covered by Q can be covered by H(B)
    bool check_coverage() {
        bool is_empty;
        ull *H_B = new ull[col_chunks]();
        for (coord i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (coord j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }
        for (coord i = 0; i < n; i++) {
            if (!(not_cov_rows.top()[i/CH_SIZE] & (ull(1) << (CH_SIZE_1 - i % CH_SIZE)))) continue;// row i is covered
            is_empty = true;
            for (coord j = 0; j < col_chunks; j++) {
                if (M[i][j] & H_B[j]) {
                    is_empty = false;
                    break;
                }
            }
            if (is_empty) {
                delete [] H_B;
                return false;
            }
        }
        delete [] H_B;
        return true;
    }

    //! Find a neighbour for constructed trajectory
    /*! Neighbour is the longest possible prefix of trajectory that can be used
    to construct a coverage. If the prefix ends with (B_t, Q_t), (B_t, Q_(t-1))
    is appended to the neighbour.*/
    bool find_neighbour() {
        while (Q.size() > 0) {
            delete [] not_cov_rows.top();
            not_cov_rows.pop(); // discard rows not covered by latest added element
            Element latest_element = Q.back(); // find latest added element
            Q.pop_back(); // discard latest added element
            update_stack(latest_element);
            if (check_covers_comp_rows()) return true;
        }
        return false;
    }

    //! Complete the trajectory
    /*! Builds the trajectory up to a point when it corresponds to coverage and
    is therefore complete. */
    void complete_trajectory() {
        // std::cout << "INITIAL:\n";
        // print_B();
        while (!check_empty() && check_covers_comp_rows()) {
            eliminate_dominating_rows();
            // std::cout << "ELIM DM ROWS:\n";
            // print_B();
            Element candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);
            // std::cout << "ELIM INCOMPAT:" << candidate.first << ' ' << candidate.second << "\n";
            // print_B();
        }
    }

    //! Checks if Q that was created is upper covering set
    /*! Only upper covering sets add coverages to results. That is done to
    avoid multiple copies of one coverage. */
    virtual bool check_upper() {
        ull *mask = new ull[col_chunks]();
        for (Element item: Q) {
            mask[item.second/CH_SIZE] |= ull(1) << (CH_SIZE_1 - item.second % CH_SIZE);
        }
        bool valid;
        for (Element item: Q) {
            for (coord i = 0; i < item.first; i++) {
                if (deleted_by_domination.top()[i/CH_SIZE] & ull(1) << (CH_SIZE_1 - i % CH_SIZE)) continue;
                valid = true;
                for (coord j = 0; j < col_chunks; j++) {
                    if ((M[i][j] & mask[j]) != (M[item.first][j] & mask[j])) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    return false;
                }
            }
        }
        return true;
    }

    //! Returns vector of columns which form Q
    std::vector<coord> get_coverage() {
        std::vector<coord> result;
        for (Element item: Q) {
            result.push_back(item.second);
        }
        return result;
    }

    void print_B() {
        std::cout << "+++++++++++++++++++++" << '\n';
        for (coord i = 0; i < n; i++) {
            if (!B[i]) {
                std::cout << "NULL" << '\n';
                continue;
            }
            for (coord j = 0; j < col_chunks-1; j++) {
                std::cout << std::bitset<CH_SIZE>(B[i][j]) << ' ';
            }
            size_t bits_left = m-(col_chunks-1)*CH_SIZE;
            for (size_t j = 0; j < bits_left; j++) {
                std::cout << ((B[i][col_chunks-1] >> (CH_SIZE_1 - j)) & ull(1));
            }
            std::cout << std::endl;
        }
        std::cout << "+++++++++++++++++++++" << '\n' << std::flush;
    }

    uint64_t get_changes_size() {
        return changes.size();
    }
};

//! Implementation of AO2 algorithm with zero extra steps
void AO2Zero(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
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
    coord HEIGHT = 50;
    coord WIDTH = 50;
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

    AO2Zero(HEIGHT, WIDTH, R, n_cov, n_extra, n_steps);

    clock_t stop = clock();
    elapsed += (double) (stop - start) / CLOCKS_PER_SEC;

    for (coord i = 0; i < HEIGHT; i++) {
        delete [] R[i];
    }
    delete [] R;

    std::cout << elapsed << " & ";
    std::cout << uint64_t(n_cov) << " & ";
    std::cout << uint64_t(n_extra) << " & ";
    std::cout << uint64_t(n_steps) << " \\\\ \n";
}
