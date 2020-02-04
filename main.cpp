#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>
#include <cmath>
#include <set>

#include "matrix_utils.h"

typedef std::stack<uint64_t**> BMatrixStack;
typedef std::pair<uint32_t, uint32_t> Coord;
typedef std::vector<Coord> QVector;
typedef std::vector<std::vector<uint32_t>> CovCollector;
typedef std::stack<std::set<uint32_t>> RowSetStack;
// typedef std::stack<bool> Log_stack;

//! Class to represent changing trajectories
/*!
  This class implements the trajectory - sequence of changing pairs (B, Q),
  which ends with a candidate for shortest coverage.

  B is represented as

  Q is represented as
*/
class Trajectory {
public:
    uint32_t n; //!< First dimension of matrix
    uint32_t m; //!< Second dimension of matrix
    uint32_t nchunks;

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

    // //! Log of taken actions. Either Q or B is updated on each step
    // /*! More detailed description. */
    // Log_stack log;

    //! Checks if B matrix is empty (Stopping criterion)
    bool check_empty() {
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (uint32_t j = 0; j < nchunks; j++) {
                if (B[i][j]) return false;
            }
        }
        return true;
    }

    //! Create delta_star set
    /*! Includes 1s from dominating rows which correspond to 0s in dominated row.*/
    uint64_t** check_dominating_rows() {
        uint64_t **delta_star = new uint64_t*[n];

        bool i1_dom_i2, i2_dom_i1;
        for (uint32_t i1 = 0; i1 < n; i1++) {
            if (!B[i1]) {
                delta_star[i1] = NULL;
                continue;
            } else {
                // () in the end sets memory to zeros
                delta_star[i1] = new uint64_t[nchunks]();
                // for (uint32_t j = 0; j < nchunks; j++) delta_star[i1][j] = 0;
            }
        }
        uint64_t *H_B = new uint64_t[nchunks]();
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (uint32_t j = 0; j < nchunks; j++) {
                H_B[j] |= B[i][j];
            }
        }
        for (uint32_t i1: not_cov_rows.top()) {
            if (!B[i1]) continue;
            for (uint32_t i2: not_cov_rows.top()) {
                if (!B[i2]) continue;
                i1_dom_i2 = i2_dom_i1 = true;
                for (uint32_t j = 0; j < nchunks; j++) {
                    if ((M[i2][j] & M[i1][j] & H_B[j]) != (M[i1][j] & H_B[j])) i2_dom_i1 = false;
                    if ((M[i1][j] & M[i2][j] & H_B[j]) != (M[i2][j] & H_B[j])) i1_dom_i2 = false;
                    if (!(i2_dom_i1 || i1_dom_i2)) break;
                }
                // std::cout << "i1 = " << i1 << "; i2 = " << i2 << " i1 dom i2 " << bool(i1_dom_i2) << " i2 dom i1 " << bool(i2_dom_i1) << '\n';
                if (!(i2_dom_i1 || i1_dom_i2)) continue;
                if (i2_dom_i1 && i1_dom_i2) continue;
                if (i2_dom_i1) {
                    for (uint32_t j = 0; j < nchunks; j++) {
                        delta_star[i2][j] |= (M[i1][j] ^ M[i2][j]) & B[i2][j];
                    }
                }
                if (i1_dom_i2) {
                    for (uint32_t j = 0; j < nchunks; j++) {
                        delta_star[i1][j] |= (M[i1][j] ^ M[i2][j]) & B[i1][j];
                    }
                }
            }
        }
        // std::cout << "DELTA STAR" << std::endl;
        // print_matrix(std::cout, delta_star, n, m);
        return delta_star;
    }

    //! Find non-zero element of B with the least index
    Coord find_the_least() {
        uint32_t least_d1 = -1;
        uint32_t least_d2 = -1;
        for (uint32_t j = 0; j < nchunks; j++) {
            for (uint32_t i = 0; i < n; i++) {
                if (!B[i]) continue;
                if (B[i][j] && (least_d2 > (j+1)*64 - int(log2(B[i][j])) - 1)) {
                    least_d1 = i;
                    least_d2 = (j+1)*64 - int(log2(B[i][j])) - 1;
                }
            }
            if (least_d1 != -1) return Coord(least_d1, least_d2);
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
    uint64_t** eliminate_incompatible(Coord element) {
        uint64_t **delta_ij = new uint64_t*[n];
        std::set<uint32_t> new_cov = not_cov_rows.top();

        for (uint32_t k = 0; k < n; k++) {
            if (!B[k]) {
                delta_ij[k] = NULL;
                continue;
            } else {
                delta_ij[k] = new uint64_t[nchunks];
                if (B[k][element.second / 64] & (uint64_t(1) << (63 - element.second % 64))) {
                    std::cout << "COLUMN " << element.second << " COVERS " << k << '\n';
                    new_cov.erase(k);
                    for (uint32_t j = 0; j < nchunks; j++) {
                        delta_ij[k][j] = B[k][j];
                    }
                } else {
                    for (uint32_t j = 0; j < nchunks; j++)
                        delta_ij[k][j] = B[k][j] & B[element.first][j];
                }
            }
        }
        std::set<uint32_t> should_delete;
        for (uint32_t row: new_cov) {
            if (M[row][element.second / 64] & (uint64_t(1) << (63 - element.second % 64))) {
                should_delete.insert(row);
                std::cout << "COLUMN " << element.second << " IN FACT COVERS " << row << '\n';
            }
        }
        for (uint32_t row: should_delete) {
            new_cov.erase(row);
        }

        not_cov_rows.push(new_cov);
        return delta_ij;
    }

public:
    Trajectory(uint32_t d1, uint32_t d2, uint64_t** Matrix) {
        n = d1;
        m = d2;
        nchunks = (m-1) / 64 + 1;
        B = new uint64_t*[n];
        for (uint32_t i = 0; i < n; i++) {
            B[i] = new uint64_t[nchunks];
            for (uint32_t j = 0; j < nchunks; j++) {
                B[i][j] = Matrix[i][j];
            }
        }
        M = new uint64_t*[n];
        for (uint32_t i = 0; i < n; i++) {
            M[i] = new uint64_t[nchunks];
            for (uint32_t j = 0; j < nchunks; j++) {
                M[i][j] = Matrix[i][j];
            }
        }
        std::set<uint32_t> tmp;
        for (uint32_t i = 0; i < n; i++) {
            tmp.insert(i);
        }
        not_cov_rows.push(tmp);
    }

    //! Recover changes in B matrix which were made earlier
    void recover_changes() {
        if (changes.empty()) {
            throw "Empty change stack. Nothing to recover";
        }
        uint64_t** latest_change = changes.top();
        changes.pop();
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) {
                B[i] = latest_change[i];
            } else {
                for (uint32_t j = 0; j < nchunks; j++) {
                    B[i][j] ^= latest_change[i][j];
                }
                delete [] latest_change[i];
            }
        }
        delete [] latest_change;
    }

    //! Save changes made in B matrix. Typically sets more B elements to 0
    void apply_changes(uint64_t** recent_changes) {
        changes.push(recent_changes);
        bool zero_b_i;
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) continue;
            zero_b_i = true;
            for (uint32_t j = 0; j < nchunks; j++) {
                B[i][j] ^= recent_changes[i][j];
                if (B[i][j]) zero_b_i = false;
            }
            if (zero_b_i) {
                delete [] B[i];
                B[i] = NULL;
            }
        }
    }

    //! Update B stack and B itself while finding a neighbouring trajectory
    /*! Converts stack from
        [...||delta_star_(t-1)||delta_ij_(t-1)||delta_star_(t)] to
        [...||delta_star_(t-1)||delta_ij_(t-1) + delta_star_(t) + latest_element|...]

        As delta_star_(t) is already applied to B, we should apply only latest_element*/
    void update_stack(Coord element) {
        if (changes.size() > 1) {
            uint64_t** delta_star_t = changes.top();
            changes.pop();
            uint64_t** delta_ij_old = changes.top();

            // std::cout << "DELTA STAR T" << std::endl;
            // print_matrix(std::cout, delta_star_t, n, m);
            // std::cout << "DELTA IJ OLD" << std::endl;
            // print_matrix(std::cout, delta_ij_old, n, m);
            for (uint32_t p = 0; p < n; p++) {
                if (!delta_ij_old[p] || !delta_star_t[p]) continue;
                for (uint32_t q = 0; q < nchunks; q++) {
                    delta_ij_old[p][q] ^= delta_star_t[p][q];
                }
                delete [] delta_star_t[p];
            }
            delete [] delta_star_t;
            delta_ij_old[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
            // std::cout << "DELTA IJ NEW" << std::endl;
            // print_matrix(std::cout, delta_ij_old, n, m);

            B[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
        } else {
            uint64_t** delta_star_t = changes.top();
            delta_star_t[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
            B[element.first][element.second / 64] ^= uint64_t(1) << (63 - element.second % 64);
        }
        std::cout << "AFTER UPDATE STACK:\n";
        print_B();
        return;
    }

    //! Checks that all rows not covered by Q can be covered by H(B)
    bool check_coverage() {
        std::set<uint32_t> not_cov = not_cov_rows.top();
        uint64_t disjunct = 0;
        bool is_empty;
        for (uint32_t row: not_cov_rows.top()) {
            is_empty = true;
            if (!B[row]) return false;
            for (uint32_t j = 0; j < nchunks; j++) {
                if (B[row][j]) {
                    is_empty = false;
                    break;
                }
            }
            if (is_empty) std::cout << "CAN NOT COVER\n";
            if (is_empty) return false;
        }
        if (is_empty) std::cout << "CAN COVER\n";
        return true;
    }

    //! Find a neighbour for constructed trajectory
    /*! Neighbour is the longest possible prefix of trajectory that can be used
    to construct a coverage. If the prefix ends with (B_t, Q_t), (B_t, Q_(t-1))
    is appended to the neighbour.*/
    bool find_neighbour() {
        std::cout << "FINDING NEIGHBOUR" << std::endl;
        while (Q.size() > 0) {
            recover_changes(); // recover from latest delta_ij
            std::cout << "RECOVERED:\n";
            print_B();
            not_cov_rows.pop(); // discard rows not covered by latest added element
            Coord latest_element = Q.back(); // find latest added element
            Q.pop_back(); // discard latest added element
            update_stack(latest_element);
            if (check_coverage()) return true;
        }
        return false;
    }

    //! Complete the trajectory
    /*! Builds the trajectory up to a point when it corresponds to coverage and
    is therefore complete. */
    void complete_trajectory() {
        while (!check_empty()) {
            std::cout << "INITIAL:\n";
            print_B();
            uint64_t** delta_star = check_dominating_rows();
            apply_changes(delta_star);
            std::cout << "AFTER DELTA STAR:\n";
            print_B();
            if (check_empty()) break;
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            uint64_t** delta_ij = eliminate_incompatible(candidate);
            apply_changes(delta_ij);
            std::cout << "AFTER DELTA IJ:\n";
            print_B();
        }
        std::cout << "Trajectory completed" << '\n';
    }

    //! Checks if Q that was created is upper covering set
    /*! Only upper covering sets add coverages to results. That is done to
    avoid multiple copies of one coverage. */
    bool check_upper() {
        uint64_t *mask = new uint64_t[nchunks]();
        for (Coord item: Q) {
            mask[item.second/64] |= uint64_t(1) << (63 - item.second % 64);
        }

        // std::cout << "MASK: " << std::bitset<10>(mask[0] >> (64-m)) << std::endl;

        bool valid;
        for (Coord item: Q) {
            for (int i = 0; i < item.first; i++) {
                valid = true;
                for (uint32_t j = 0; j < nchunks; j++) {
                    // std::cout << "ELEM: " << std::bitset<10>(M[i][j] >> (64-m)) << ' ' << std::bitset<10>(M[item.first][j] >> (64-m))<< std::endl;
                    // std::cout << "ELEM: " << std::bitset<10>((M[i][j] & mask[j]) >> (64-m)) << ' ' << std::bitset<10>((M[item.first][j] & mask[j]) >> (64-m)) << std::endl;

                    if ((M[i][j] & mask[j]) != (M[item.first][j] & mask[j])) {
                        valid = false;
                        // std::cout << "NOT EQUAL" << '\n';
                        break;
                    }
                }
                if (valid) {
                    std::cout << "NOT UPPER" << '\n';
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
        //std::cout << "NCHUNKS:" << nchunks << '\n';
        std::cout << "+++++++++++++++++++++" << '\n';
        for (uint32_t i = 0; i < n; i++) {
            if (!B[i]) {
                std::cout << "NULL" << '\n';
                continue;
            }
            for (uint32_t j = 0; j < nchunks-1; j++) {
                std::cout << std::bitset<64>(B[i][j]) << ' ';
            }
            size_t bits_left = m-(nchunks-1)*64;
            for (size_t j = 0; j < bits_left; j++) {
                std::cout << ((B[i][nchunks-1] >> (63 - j)) & uint64_t(1));
            }
            std::cout << std::endl;
        }
        std::cout << "+++++++++++++++++++++" << '\n';
    }
};

int main(int argc, char *argv[]) {
    uint32_t n = 5;
    uint32_t m = 5;
    uint32_t nchunks = m / 64 + 1 - (m % 64 == 0);

    generate_matrix(n, m, "matrix.txt", 0.5, 36);
    uint64_t** R = read_matrix("matrix.txt", n, m);

    CovCollector coverages;
    Trajectory traj(n, m, R);
    // std::cout << traj.find_the_least().first << '\n';
    // std::cout << traj.find_the_least().second << '\n';

    int counter = 0;
    do {
        traj.complete_trajectory();
        std::cout << "TRAJ: ";
        for (auto item: traj.Q) {
            std::cout << '(' << item.first << ' ' << item.second << ") ";
        }
        std::cout << '\n';
        std::cout << "RES COV: ";
        for (uint32_t col: traj.get_coverage()) {
            std::cout << col << ' ';
        }
        if (traj.check_upper()) {
            coverages.push_back(traj.get_coverage());
            std::cout << "SUCCESS";
        };
        std::cout << '\n';
        counter ++;
    } while (traj.find_neighbour() && (counter < 3));

    std::cout << "FOUND COVERAGES:" << std::endl;
    for (auto cov: coverages) {
        for (uint32_t col: cov) {
            std::cout << col << ' ';
        }
        std::cout << std::endl;
    }
}
