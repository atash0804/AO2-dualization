#ifndef DUALIZATION_AO2_TRAJECTORY
#define DUALIZATION_AO2_TRAJECTORY

#include <iostream>
#include <stack>
#include <utility>      // std::pair
#include <bitset>
#include <vector>

#include "matrix_utils.h"

typedef uint32_t ull; // typedef for chunks in which matrix is stored
typedef uint32_t c_int; // typedef for matrix shape and coordinates

typedef std::stack<ull**> BMatrixStack;
typedef std::pair<c_int, c_int> Coord;
typedef std::vector<Coord> QVector;
typedef std::vector<std::vector<c_int>> CovCollector;
typedef std::stack<ull*> RowSetStack;

enum {
    CH_SIZE = sizeof(ull) * 8,
    CH_SIZE_1 = sizeof(ull) * 8 - 1
};

//! Class to represent changing trajectories
/*!
  This class implements the trajectory - sequence of changing pairs (B, Q),
  which ends with a candidate for shortest coverage and all operations on
  trajectories which are needed for AO2.

*/
class AO2Trajectory {
protected:
    c_int n; //!< First dimension of matrix
    c_int m; //!< Second dimension of matrix
    c_int col_chunks;
    c_int row_chunks;

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
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (c_int j = 0; j < col_chunks; j++) {
                if (B[i][j]) return false;
            }
        }
        return true;
    }

    //! Create delta_star set
    /*! Includes 1s from dominating rows which correspond to 0s in dominated row.*/
    void eliminate_dominating_rows() {
        ull **updated_B = new ull*[n];
        bool i1_dom_i2, i2_dom_i1;
        for (c_int i1 = 0; i1 < n; i1++) {
            if (!B[i1]) {
                updated_B[i1] = NULL;
                continue;
            } else {
                updated_B[i1] = new ull[col_chunks]();
                for (c_int j = 0; j < col_chunks; j++) updated_B[i1][j] = B[i1][j];
            }
        }
        ull *H_B = new ull[col_chunks]();
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (c_int j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }

        ull* new_del = new ull[row_chunks];
        for (c_int i = 0; i < row_chunks; i++) {
            new_del[i] = deleted_by_domination.top()[i];
        }

        for (c_int i1 = 0; i1 < n; i1++) {
            if (!B[i1] || (!(not_cov_rows.top()[i1/CH_SIZE] & (ull(1) << (CH_SIZE_1-i1%CH_SIZE))))) continue;
            for (c_int i2 = i1+1; i2 < n; i2++) {
                if (!B[i2] || (!(not_cov_rows.top()[i2/CH_SIZE] & (ull(1) << (CH_SIZE_1-i2%CH_SIZE))))) continue;
                i1_dom_i2 = i2_dom_i1 = true;
                for (c_int j = 0; j < col_chunks; j++) {
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
        for (c_int i = 0; i < n; i++) delete [] B[i];
        delete [] B;
        deleted_by_domination.push(new_del);
        B = updated_B;
        return;
    }

    //! Find non-zero element of B with the least index
    virtual Coord find_the_least() {
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

    //! Create delta_ij set
    /*! Includes elements of B incompatible with element. Suppose the element
    is B_ij. Incompatible elements are all B_kl such that either B_kj or B_il
    equals 1.
    To identify incompatible elements, we cycle throw rows k=1..n and look at
    B_kj. If it equals 1, all 1s of the row are incomatible. Otherwise, if
    B_kj = 0, the incomatible elements are all 1s from B_k & B_i (elementwise
    conjunction of two rows)*/
    void eliminate_incompatible(Coord& element) {
        ull **updated_B = new ull*[n];
        ull* new_cov = new ull[row_chunks];
        for (c_int i = 0; i < row_chunks; i++) {
            new_cov[i] = not_cov_rows.top()[i];
        }

        for (c_int k = 0; k < n; k++) {
            if (!B[k]) {
                if (M[k][element.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - element.second % CH_SIZE))) {
                    new_cov[k/CH_SIZE] ^= new_cov[k/CH_SIZE] & ( ull(1) << (CH_SIZE_1 - k % CH_SIZE));
                }
                updated_B[k] = NULL;
                continue;
            } else {
                if (M[k][element.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - element.second % CH_SIZE))) {
                    new_cov[k/CH_SIZE] ^= new_cov[k/CH_SIZE] & (ull(1) << (CH_SIZE_1 - k % CH_SIZE));
                    updated_B[k] = NULL;
                } else {
                    updated_B[k] = new ull[col_chunks]();
                    for (c_int j = 0; j < col_chunks; j++)
                        updated_B[k][j] = B[k][j] ^ (B[k][j] & M[element.first][j]);
                }
            }
        }

        ull* tmp = new ull[row_chunks];
        for (c_int i = 0; i < row_chunks; i++) {
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
    void update_stack(Coord& element) {
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

            for (c_int p = 0; p < n; p++) {
                delete [] B[p];
            }
            delete [] B;
            B = latest_state;
            B[element.first][element.second / CH_SIZE] ^= ull(1) << (CH_SIZE_1 - element.second % CH_SIZE);
        } else {
            B[element.first][element.second / CH_SIZE] ^= ull(1) << (CH_SIZE_1 - element.second % CH_SIZE);
        }
        bool is_zero = true;
        for (c_int j = 0; j < col_chunks; j++) {
            if (B[element.first][j]) is_zero = false;
        }
        if (is_zero) {
            delete [] B[element.first];
            B[element.first] = NULL;
        }
        return;
    }
public:
    AO2Trajectory(c_int d1, c_int d2, ull** Matrix) {
        n = d1;
        m = d2;
        col_chunks = (m-1) / CH_SIZE + 1;
        row_chunks = (n-1) / CH_SIZE + 1;
        B = new ull*[n];
        M = new ull*[n];
        for (c_int i = 0; i < n; i++) {
            B[i] = new ull[col_chunks];
            M[i] = new ull[col_chunks];
            for (c_int j = 0; j < col_chunks; j++) {
                B[i][j] = M[i][j] = Matrix[i][j];
            }
        }
        ull* tmp = new ull[row_chunks];
        for (c_int i = 0; i <  row_chunks; i++) {
            tmp[i] = ull(-1);
        }
        not_cov_rows.push(tmp);
        tmp = new ull[row_chunks];
        for (c_int i = 0; i <  row_chunks; i++) {
            tmp[i] = ull(0);
        }
        deleted_by_domination.push(tmp);
    }

    ~AO2Trajectory() {
        for (c_int i = 0; i < n; i++) {
            if (B[i]) delete [] B[i];
            if (M[i]) delete [] M[i];
        }
        delete [] B;
        delete [] M;

        for (; !changes.empty(); changes.pop()) {
            for (c_int i = 0; i < n; i++) {
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
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) continue;
            for (c_int j = 0; j < col_chunks; j++) {
                H_B[j] |= B[i][j];
            }
        }
        for (c_int i = 0; i < n; i++) {
            if (!(not_cov_rows.top()[i/CH_SIZE] & (ull(1) << (CH_SIZE_1 - i % CH_SIZE)))) continue;// row i is covered
            is_empty = true;
            for (c_int j = 0; j < col_chunks; j++) {
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
    virtual void complete_trajectory() {
        while (!check_empty()) {
            eliminate_dominating_rows();
            Coord candidate = find_the_least();
            Q.push_back(candidate);
            eliminate_incompatible(candidate);
        }
    }

    //! Checks if Q that was created is upper covering set
    /*! Only upper covering sets add coverages to results. That is done to
    avoid multiple copies of one coverage. */
    virtual bool check_upper() {
        ull *mask = new ull[col_chunks]();
        for (Coord item: Q) {
            mask[item.second/CH_SIZE] |= ull(1) << (CH_SIZE_1 - item.second % CH_SIZE);
        }
        bool valid;
        for (Coord item: Q) {
            for (c_int i = 0; i < item.first; i++) {
                if (deleted_by_domination.top()[i/CH_SIZE] & ull(1) << (CH_SIZE_1 - i % CH_SIZE)) continue;
                valid = true;
                for (c_int j = 0; j < col_chunks; j++) {
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
    std::vector<c_int> get_coverage() {
        std::vector<c_int> result;
        for (Coord item: Q) {
            result.push_back(item.second);
        }
        return result;
    }

    void print_B() {
        std::cout << "+++++++++++++++++++++" << '\n';
        for (c_int i = 0; i < n; i++) {
            if (!B[i]) {
                std::cout << "NULL" << '\n';
                continue;
            }
            for (c_int j = 0; j < col_chunks-1; j++) {
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

#endif
