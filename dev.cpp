#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

#include "matrix_utils.h"

typedef std::stack<uint64_t**> BMatrixStack;
typedef std::pair<uint32_t, uint32_t> Coord;
typedef std::vector<Coord> QVector;
typedef std::vector<std::vector<uint32_t>> CovCollector;
typedef std::stack<uint64_t*> RowSetStack;

//! Class to represent changing trajectories
/*!
  This class implements the trajectory - sequence of changing pairs (B, Q),
  which ends with a candidate for shortest coverage and all operations on
  trajectories which are needed for AO2.

*/
class Trajectory {
public:
    uint32_t n; //!< First dimension of matrix
    uint32_t m; //!< Second dimension of matrix
    uint32_t col_chunks;
    uint32_t row_chunks;

    //! Stores initial state of matrix
    uint64_t** M;

    //! Representation of set B (elements that can be used to complete the trajectory)
    /*! B is represented by binary matrix of elements that can still be used
    to complete the trajectory*/
    uint64_t** B;

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

    //! Checks if B matrix is empty (Stopping criterion)
    bool check_empty() {
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (uint32_t j = 0; j < col_chunks; j++) {
                if (B[i][j]) return false;
            }
        }
        return true;
    }

    //! Create delta_star set
    /*! Includes 1s from dominating rows which correspond to 0s in dominated row.*/
    void eliminate_dominating_rows() {
        uint64_t **updated_B = new uint64_t*[n];
        bool i1_dom_i2, i2_dom_i1;
        for (uint32_t i1 = 0; i1 < n; i1++) {
            if (!B[i1]) {
                updated_B[i1] = NULL;
                continue;
            } else {
                // () in the end sets memory to zeros
                updated_B[i1] = new uint64_t[col_chunks]();
                for (uint32_t j = 0; j < col_chunks; j++) updated_B[i1][j] = B[i1][j];
                // for (uint32_t j = 0; j < col_chunks; j++) delta_star[i1][j] = 0;
            }
        }
        uint64_t *H_B = new uint64_t[col_chunks]();
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (uint32_t j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }

        bool is_zero;
        for (uint32_t i1 = 0; i1 < n; i1++) {
            if (!B[i1] || (!(not_cov_rows.top()[i1/64] & (uint64_t(1) << (63-i1%64))))) continue;
            for (uint32_t i2 = i1+1; i2 < n; i2++) {
                if (!B[i2] || (!(not_cov_rows.top()[i2/64] & (uint64_t(1) << (63-i2%64))))) continue;
                i1_dom_i2 = i2_dom_i1 = true;
                for (uint32_t j = 0; j < col_chunks; j++) {
                    if ((M[i2][j] & M[i1][j] & H_B[j]) != (M[i1][j] & H_B[j])) i2_dom_i1 = false;
                    if ((M[i1][j] & M[i2][j] & H_B[j]) != (M[i2][j] & H_B[j])) i1_dom_i2 = false;
                    if (!(i2_dom_i1 || i1_dom_i2)) break;
                }
                // std::cout << "i1 = " << i1 << "; i2 = " << i2 << " i1 dom i2 " << bool(i1_dom_i2) << " i2 dom i1 " << bool(i2_dom_i1) << '\n';
                if (!(i2_dom_i1 || i1_dom_i2)) continue;
                if (i2_dom_i1 && i1_dom_i2) continue;
                if (i2_dom_i1) {
                    if (!updated_B[i2]) continue;
                    is_zero = true;
                    for (uint32_t j = 0; j < col_chunks; j++) {
                        updated_B[i2][j] ^= updated_B[i2][j] & (M[i1][j] ^ M[i2][j]);
                        if (updated_B[i2][j]) is_zero = false;
                    }
                    if (is_zero) {
                        delete [] updated_B[i2];
                        updated_B[i2] = NULL;
                    }
                }
                if (i1_dom_i2) {
                    is_zero = true;
                    if (!updated_B[i1]) continue;
                    for (uint32_t j = 0; j < col_chunks; j++) {
                        updated_B[i1][j] ^= updated_B[i1][j] & (M[i1][j] ^ M[i2][j]);
                        if (updated_B[i1][j]) is_zero = false;
                    }
                    if (is_zero) {
                        delete [] updated_B[i1];
                        updated_B[i1] = NULL;
                    }
                }
            }
        }
        changes.push(B);
        B = updated_B;
        return;
    }

    //! Find non-zero element of B with the least index
    Coord find_the_least() {
        uint32_t least_d1 = -1;
        uint32_t least_d2 = -1;
        for (uint32_t j = 0; j < col_chunks; j++) {
            for (uint32_t i = 0; i < n; i++) {
                if (!B[i]) continue;
                if (B[i][j] && (least_d2 > (j+1)*64 - int(log2(B[i][j])) - 1)) {
                    least_d1 = i;
                    least_d2 = (j+1)*64 - int(log2(B[i][j])) - 1;
                }
            }
            if (least_d1 != uint64_t(-1)) return Coord(least_d1, least_d2);
        }
        return Coord(least_d1, least_d2);
    }

    //! Create delta_ij set
    /*! Includes elements of B incompatible with element. Suppose the element
    is B_ij. Incompatible elements are all B_kl such that either B_kj or B_il
    equals 1.
    To identify incompatible elements, we cycle throw rows k=1..n and look at
    B_kj. If it equals 1, all 1s of the row are incomatible. Otherwise, if
    B_kj = 0, the incomatible elements are all 1s from B_k & B_i (elementwise
    conjunction of two rows)*/
    void eliminate_incompatible(Coord element) {
        // std::cout << "UNDER QUESTION\n" << std::flush;
        // print_B();
        uint64_t **updated_B = new uint64_t*[n];
        uint64_t* new_cov = new uint64_t[row_chunks];
        for (uint32_t i = 0; i < row_chunks; i++) {
            new_cov[i] = not_cov_rows.top()[i];
        }

        for (uint32_t k = 0; k < n; k++) {
            if (!B[k]) {
                if (M[k][element.second / 64] & (uint64_t(1) << (63 - element.second % 64))) {
                    new_cov[k/64] ^= new_cov[k/64] & ( uint64_t(1) << (63 - k % 64));
                }
                updated_B[k] = NULL;
                continue;
            } else {
                if (M[k][element.second / 64] & (uint64_t(1) << (63 - element.second % 64))) {
                    // std::cout << "COLUMN " << element.second << " COVERS " << k << '\n';
                    new_cov[k/64] ^= new_cov[k/64] & (uint64_t(1) << (63 - k % 64));
                    updated_B[k] = NULL;
                } else {
                    updated_B[k] = new uint64_t[col_chunks]();
                    for (uint32_t j = 0; j < col_chunks; j++)
                        updated_B[k][j] = B[k][j] ^ (B[k][j] & M[element.first][j]);
                }
            }
        }

        not_cov_rows.push(new_cov);
        changes.push(B);
        B = updated_B;
        // std::cout << "AFTER QUESTION\n" << std::flush;
        // print_B();
        return;
    }

public:
    Trajectory(uint32_t d1, uint32_t d2, uint64_t** Matrix) {
        n = d1;
        m = d2;
        col_chunks = (m-1) / 64 + 1;
        row_chunks = (n-1) / 64 + 1;
        B = new uint64_t*[n];
        M = new uint64_t*[n];
        for (uint32_t i = 0; i < n; i++) {
            B[i] = new uint64_t[col_chunks];
            M[i] = new uint64_t[col_chunks];
            for (uint32_t j = 0; j < col_chunks; j++) {
                B[i][j] = M[i][j] = Matrix[i][j];
            }
        }
        uint64_t* tmp = new uint64_t[row_chunks];
        for (uint32_t i = 0; i <  row_chunks; i++) {
            tmp[i] = uint64_t(-1);
        }
        not_cov_rows.push(tmp);
    }

    //! Recover changes in B matrix which were made earlier
    void recover_changes() {
        std::cout << "STACK SIZE:\n" << changes.size() << std::endl << std::flush;
        if (changes.empty()) {
            std::cout << "CHANGE STACK EMPTY:\n" << std::flush;
            throw "Empty change stack. Nothing to recover";
        }
        // if (changes.size() == 1) {
        //     uint64_t** latest_state = changes.top();
        //     print_matrix(std::cout, latest_state, n, m);
        //     // std::cout << "LATEST STATE:" << latest_state << std::endl << std::flush;
        //     for (uint32_t i = 0; i < n; i++) {
        //         // std::cout << "i:" << i << latest_state[i] << std::endl << std::flush;
        //
        //         if (!latest_state[i]) {
        //             B[i] = NULL;
        //             continue;
        //         }
        //         if (!B[i]) B[i] = new uint64_t[col_chunks];
        //         for (uint32_t j = 0; j < col_chunks; j++) {
        //             B[i][j] = latest_state[i][j];
        //         }
        //     }
        //     return;
        // }
        uint64_t** latest_state = changes.top();
        changes.pop();
        if (!latest_state) return;
        for (uint32_t i = 0; i < n; i++) {
            delete [] B[i];
        }
        delete [] B;
        B = latest_state;
    }

    //! Update B stack and B itself while finding a neighbouring trajectory
    /*! Converts stack from
        [...||delta_star_(t-1)||delta_ij_(t-1)||delta_star_(t)] to
        [...||delta_star_(t-1)||delta_ij_(t-1) + delta_star_(t) + latest_element|...]

        As delta_star_(t) is already applied to B, we should apply only latest_element*/
    void update_stack(Coord element) {
        if (changes.size() > 1) {
            uint64_t** latest_state = changes.top();
            changes.pop();
            uint64_t** previous_state = changes.top();
            changes.pop();
            for (uint32_t p = 0; p < n; p++) {
                delete [] previous_state[p];
                delete [] B[p];
            }
            delete [] previous_state;
            delete [] B;

            B = latest_state;
            B[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
        } else {
            B[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
        }
        bool is_zero = true;
        for (uint32_t j = 0; j < col_chunks; j++) {
            if (B[element.first][j]) is_zero = false;
        }
        if (is_zero) {
            delete [] B[element.first];
            B[element.first] = NULL;
        }
        // std::cout << "AFTER UPDATE STACK:\n" << std::flush;
        // print_B();
        return;
    }

    //! Checks that all rows not covered by Q can be covered by H(B)
    bool check_coverage() {
        uint64_t* not_cov = new uint64_t[row_chunks];
        for (uint32_t i = 0; i <  row_chunks; i++) {
            not_cov[i] = not_cov_rows.top()[i];
        }
        bool is_empty;
        uint64_t *H_B = new uint64_t[col_chunks]();
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (uint32_t j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }
        for (uint32_t i = 0; i < n; i++) {
            if (!(not_cov[i/64] & (uint64_t(1) << (63 - i % 64)))) continue;// row i is covered
            is_empty = true;
            for (uint32_t j = 0; j < col_chunks; j++) {
                if (M[i][j] & H_B[j]) {
                    is_empty = false;
                    break;
                }
            }
            // if (is_empty) std::cout << "CAN NOT COVER\n";
            if (is_empty) return false;
        }
        // std::cout << "CAN COVER\n";
        return true;
    }

    //! Find a neighbour for constructed trajectory
    /*! Neighbour is the longest possible prefix of trajectory that can be used
    to construct a coverage. If the prefix ends with (B_t, Q_t), (B_t, Q_(t-1))
    is appended to the neighbour.*/
    bool find_neighbour() {
        // std::cout << "FINDING NEIGHBOUR" << std::endl;
        while (Q.size() > 0) {
            // std::cout << "BEFORE REC CHANGES:\n" << std::flush;
            // recover_changes(); // recover from latest delta_ij
            // std::cout << "RECOVERED:\n" << std::flush;
            // print_B();
            delete [] not_cov_rows.top();
            not_cov_rows.pop(); // discard rows not covered by latest added element
            Coord latest_element = Q.back(); // find latest added element
            Q.pop_back(); // discard latest added element
            update_stack(latest_element);
            // std::cout << "BEFORE CHECK COV:\n" << std::flush;
            if (check_coverage()) return true;
            // std::cout << "AFTER CHECK COV:\n" << std::flush;
        }
        return false;
    }

    //! Complete the trajectory
    /*! Builds the trajectory up to a point when it corresponds to coverage and
    is therefore complete. */
    void complete_trajectory() {
        while (!check_empty()) {
            // std::cout << "INITIAL:\n";
            // print_B();
            eliminate_dominating_rows();
            // std::cout << "AFTER DELTA STAR:\n";
            // print_B();
            if (check_empty()) break;
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);
            // std::cout << "AFTER DELTA IJ:\n";
            // std::cout << "CANDIDATE " << candidate.first << ' ' << candidate.second << std::endl;
            // print_B();
        }
        // std::cout << "Trajectory completed" << '\n';
    }

    //! Checks if Q that was created is upper covering set
    /*! Only upper covering sets add coverages to results. That is done to
    avoid multiple copies of one coverage. */
    bool check_upper() {
        uint64_t *mask = new uint64_t[col_chunks]();
        for (Coord item: Q) {
            mask[item.second/64] |= uint64_t(1) << (63 - item.second % 64);
        }

        // std::cout << "MASK: " << std::bitset<10>(mask[0] >> (64-m)) << std::endl;

        bool valid;
        for (Coord item: Q) {
            for (uint32_t i = 0; i < item.first; i++) {
                valid = true;
                for (uint32_t j = 0; j < col_chunks; j++) {
                    if ((M[i][j] & mask[j]) != (M[item.first][j] & mask[j])) {
                        valid = false;
                        // std::cout << "NOT EQUAL" << '\n';
                        break;
                    }
                }
                if (valid) {
                    // std::cout << "NOT UPPER" << '\n';
                    return false;
                }
            }
        }
        return true;
    }

    //! Returns vector of columns which form Q
    std::vector<uint32_t> get_coverage() {
        std::vector<uint32_t> result;
        for (Coord item: Q) {
            result.push_back(item.second);
        }
        return result;
    }

    void print_B() {
        std::cout << "+++++++++++++++++++++" << '\n';
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) {
                std::cout << "NULL" << '\n';
                continue;
            }
            for (uint32_t j = 0; j < col_chunks-1; j++) {
                std::cout << std::bitset<64>(B[i][j]) << ' ';
            }
            size_t bits_left = m-(col_chunks-1)*64;
            for (size_t j = 0; j < bits_left; j++) {
                std::cout << ((B[i][col_chunks-1] >> (63 - j)) & uint64_t(1));
            }
            std::cout << std::endl;
        }
        std::cout << "+++++++++++++++++++++" << '\n' << std::flush;
    }
};

int main(int argc, char *argv[]) {
    uint32_t n = 30;
    uint32_t m = n;
    uint32_t col_chunks = m / 64 + 1 - (m % 64 == 0);

    generate_matrix(n, m, "matrix.txt", 0.5, 31);
    uint64_t** R = read_matrix("matrix.txt", n, m);

    CovCollector coverages;
    Trajectory traj(n, m, R);
    // std::cout << traj.find_the_least().first << '\n';
    // std::cout << traj.find_the_least().second << '\n';
    clock_t start = clock();
    do {
        traj.complete_trajectory();
        // std::cout << "COMPLETE STACK COND:";
        // for (auto dump = traj.changes; !dump.empty(); dump.pop())
        //     print_matrix(std::cout, dump.top(), n, m);
        // std::cout << "DONE\n";
        // std::cout << "||||||||||||||||||||||||||||TRAJ: \n";
        // for (auto item: traj.Q) {
        //     std::cout << '(' << item.first << ' ' << item.second << ") ";
        // }
        // std::cout << '\n';
        // std::cout << "RES COV: ";
        // for (uint32_t col: traj.get_coverage()) {
        //     std::cout << col << ' ';
        // }
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
            // std::cout << "SUCCESS";
        };
        // std::cout << '\n';
    } while (traj.find_neighbour());
    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    std::cout << "AO2 dualization time: " << elapsed << '\n';
    std::cout << "N coverages: " << coverages.size() << '\n';

    // std::cout << "FOUND COVERAGES:" << '\n';
    // for (auto cov: coverages) {
    //     for (uint32_t col: cov) {
    //         std::cout << col << ' ';
    //     }
    //     std::cout << '\n';
    // }
    // std::cout << "PERFORMING CHECKS" << coverages.size() << '\n';
    // for (auto cov: coverages) {
    //     uint64_t *H_cov = new uint64_t[col_chunks]();
    //     for (uint32_t col: cov) {
    //         H_cov[col/64] |= uint64_t(1) << (63 - col % 64);
    //     }
    //     bool is_cov;
    //     for (uint32_t i = 0; i < n; i++) {
    //         is_cov = false;
    //         for (uint32_t j = 0; j < col_chunks; j++) {
    //             if (R[i][j] & H_cov[j]) {is_cov = true; break;}
    //         }
    //         if (!is_cov) {
    //             std::cout << "ENTITY FOUND IS NOT COVERAGE\n";
    //         };
    //     }
    // }

}
