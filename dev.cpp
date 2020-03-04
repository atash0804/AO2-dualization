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

//! Class to represent changing trajectories
/*!
  This class implements the trajectory - sequence of changing pairs (B, Q),
  which ends with a candidate for shortest coverage and all operations on
  trajectories which are needed for AO2.

*/
class AO2MoptimizedTrajectory {
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

    //! Stack of states of B rows
    /*! Each state is represented as *st */
    RowStatesStack states;

    //! Representation of set Q (elements that form the trajectory)
    /*! Vector of pairs (i_r, j_r) which form the Q set*/
    QVector Q;

    //! Element which was previously added to Q
    Element latest_element;

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
        bool is_empty, valid;
        // ull *H_B = new ull[col_chunks]();
        // for (coord i = 0; i < n; i++) {
        //     if (!B[i]) continue;
        //     for (coord j = 0; j < col_chunks; j++) {
        //         H_B[j] |= B[i][j];
        //     }
        // }

        // coord count_rows = 0;
        for (coord i = 0; i < n; i++) {
            // if row is not covered or is competing
            if ((!(states.top()[i] & ST_IS_COV)) || (states.top()[i] == ST_COMP)) {
                // count_rows++;
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
            }
        }
        return true;
    }

    //! Create delta_star set
    /*! Includes 1s from dominating rows which correspond to 0s in dominated row.*/
    bool eliminate_dominating_rows() {
        bool has_changed = false;
        ull **updated_B = new ull*[n];
        st* updated_state = new st[n];
        for (coord i = 0; i < n; i++) updated_state[i] = states.top()[i];
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

        for (coord i1 = 0; i1 < n; i1++) {
            if (!B[i1]) continue;
            for (coord i2 = i1+1; i2 < n; i2++) {
                if (!B[i2]) continue;
                i1_dom_i2 = i2_dom_i1 = true;
                for (coord j = 0; j < col_chunks; j++) {
                    if ((M[i2][j] & M[i1][j] & H_B[j]) != (M[i1][j] & H_B[j])) i2_dom_i1 = false;
                    if ((M[i1][j] & M[i2][j] & H_B[j]) != (M[i2][j] & H_B[j])) i1_dom_i2 = false;
                    if (!(i2_dom_i1 || i1_dom_i2)) break;
                }
                if (i2_dom_i1) {
                    has_changed = true;
                    delete [] updated_B[i2];
                    updated_state[i2] = (updated_state[i2] & 3) | ST_IS_DOM;
                    updated_B[i2] = NULL;
                } else {
                    if (i1_dom_i2) {
                        has_changed = true;
                        delete [] updated_B[i1];
                        updated_state[i1] = (updated_state[i1] & 3) | ST_IS_DOM;
                        updated_B[i1] = NULL;
                    }
                }
            }
        }
        delete [] H_B;
        for (coord i = 0; i < n; i++) delete [] B[i];
        delete [] B;
        states.push(updated_state);
        B = updated_B;
        return has_changed;
    }

    //! Find non-zero element of B with the least index
    virtual Element find_the_least() {
        bool use_cached = true;
        if (latest_element.first == coord(-1)) use_cached = false;
        if (use_cached) {
            for (coord i = latest_element.first + 1; i < n; i++) {
                if (!B[i]) continue;
                if (B[i][latest_element.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - latest_element.second % CH_SIZE))) {
                    latest_element = Element(-1, latest_element.second);
                    return Element(i, latest_element.second);
                }
            }
            use_cached = false;
        }
        if (!use_cached) {
            coord least_d1 = -1;
            coord least_d2 = -1;
            coord least_Er = n*m+1, curr_Er;
            coord *columns = new coord[m]();

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
    }

    //! Create delta_ij set
    /*! Includes elements of B incompatible with element. Suppose the element
    is B_ij. Incompatible elements are all B_kl such that either B_kj or B_il
    equals 1.
    To identify incompatible elements, we cycle throw rows k=1..n and look at
    B_kj. If it equals 1, all 1s of the row are incomatible. Otherwise, if
    B_kj = 0, the incomatible elements are all 1s from B_k & B_i (elementwise
    conjunction of two rows)*/
    bool eliminate_incompatible(Element& el) {
        ull **updated_B = new ull*[n];

        st* updated_state = new st[n];
        for (coord i = 0; i < n; i++) updated_state[i] = states.top()[i];
        bool has_comp = false, has_uncov = false;
        for (coord k = 0; k < n; k++) {
            if (M[k][el.second / CH_SIZE] & (ull(1) << (CH_SIZE_1 - el.second % CH_SIZE))) {
                // selected column covers row k
                switch (states.top()[k]) {
                case ST_EMPTY:
                    updated_state[k] = (k < el.first) ? ST_COMP : ST_COV;
                    break;
                case ST_DOM_NOT_COV:
                    updated_state[k] = ST_DOM_COV;
                    break;
                case ST_COMP:
                    updated_state[k] = ST_COV;
                    break;
                }
                updated_B[k] = NULL;
            } else {
                if (!B[k]) {
                    updated_B[k] = NULL;
                } else {
                    updated_B[k] = new ull[col_chunks]();
                    for (coord j = 0; j < col_chunks; j++)
                        updated_B[k][j] = B[k][j] ^ (B[k][j] & M[el.first][j]);
                }
            }

            if (updated_state[k] == ST_COMP) has_comp = true;
            if (!(updated_state[k] & ST_IS_COV)) {
                has_uncov = true;
            }
        }
        states.push(updated_state);
        changes.push(B);
        B = updated_B;
        if ((!has_uncov) && has_comp) {
            return false;
        }
        return true;
    }

    //! Update B stack and B itself while finding a neighbouring trajectory
    /*! Converts stack from
        [...||delta_star_(t-1)||delta_ij_(t-1)||delta_star_(t)] to
        [...||delta_star_(t-1)||delta_ij_(t-1) + delta_star_(t) + latest_element|...]

        As delta_star_(t) is already applied to B, we should apply only latest_element*/
    void update_stack(Element& el) {
        ull** latest_state = changes.top();
        changes.pop();

        delete [] states.top();
        states.pop();
        st* B_state = states.top();
        states.pop();
        delete [] states.top();
        states.pop();
        states.push(B_state);

        for (coord p = 0; p < n; p++) {
            delete [] B[p];
        }
        delete [] B;
        B = latest_state;
        B[el.first][el.second / CH_SIZE] ^= ull(1) << (CH_SIZE_1 - el.second % CH_SIZE);

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
    AO2MoptimizedTrajectory(coord d1, coord d2, ull** Matrix) {
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

        st* tmp = new st[n]();
        states.push(tmp);

        latest_element = Element(-1, -1);
    }

    ~AO2MoptimizedTrajectory() {
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

        for (; !states.empty(); states.pop()) {
            delete [] states.top();
        }
    }

    //! Find a neighbour for constructed trajectory
    /*! Neighbour is the longest possible prefix of trajectory that can be used
    to construct a coverage. If the prefix ends with (B_t, Q_t), (B_t, Q_(t-1))
    is appended to the neighbour.*/
    bool find_neighbour() {
        while (Q.size() > 0) {
            latest_element = Q.back(); // find latest added element
            Q.pop_back(); // discard latest added element
            update_stack(latest_element);
            if (check_covers_comp_rows()) return true;
        }
        return false;
    }

    //! Complete the trajectory
    /*! Builds the trajectory up to a point when it corresponds to coverage and
    is therefore complete. */
    virtual bool complete_trajectory() {
        while (!check_empty()) {
            if (!check_covers_comp_rows()) return false;
            if (eliminate_dominating_rows() && !check_covers_comp_rows()) {
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
                std::cout << "NULL";
            } else {
                for (coord j = 0; j < col_chunks-1; j++) {
                    std::cout << std::bitset<CH_SIZE>(B[i][j]) << ' ';
                }
                size_t bits_left = m-(col_chunks-1)*CH_SIZE;
                for (size_t j = 0; j < bits_left; j++) {
                    std::cout << ((B[i][col_chunks-1] >> (CH_SIZE_1 - j)) & ull(1));
                }
            }
            switch (states.top()[i]) {
            case ST_EMPTY:
                std::cout << "(NOT COV)\n";
                break;
            case ST_DOM_NOT_COV:
                std::cout << "(NOT COV, DEL AS DOM)\n";
                break;
            case ST_DOM_COV:
                std::cout << "(COV, DEL AS DOM)\n";
                break;
            case ST_COV:
                std::cout << "(COVERED)\n";
                break;
            case ST_COMP:
                std::cout << "(COMPETES)\n";
                break;
            }

        }
        std::cout << "+++++++++++++++++++++" << '\n' << std::flush;
    }

    uint64_t get_changes_size() {
        return changes.size();
    }
};

void AO2Moptimized(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2MoptimizedTrajectory traj(n, m, R);
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
    coord HEIGHT = 50;
    coord WIDTH = 50;
    double SPARSITY = 0.5;

    double elapsed = 0;
    uint64_t n_cov1 = 0, n_cov2 = 0, n_cov3 = 0, n_cov4 = 0;
    uint64_t n_extra = 0;
    uint64_t n_steps = 0;

    srand(time(NULL));
    int i = 0;
    while (true) {
        if (i > 0) break;
        cout << i++ << '\n';
        generate_matrix(HEIGHT, WIDTH, "matrix.txt", SPARSITY);
        ull** R = read_matrix("matrix.txt", HEIGHT, WIDTH);

        if (has_zero_rows(R, HEIGHT, WIDTH)) {
            std::cout << "Matrix contains zero rows \n" << std::endl;
            continue;
        }

        clock_t start = clock();

        n_cov1 = n_extra = n_steps = 0;
        AO2Mplus(HEIGHT, WIDTH, R, n_cov1, n_extra, n_steps);

        clock_t stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        std::cout << elapsed << " & ";
        std::cout << uint64_t(n_cov1) << " & ";
        std::cout << uint64_t(n_extra) << " & ";
        std::cout << uint64_t(n_steps) << " \\\\ \n";

        start = clock();

        n_cov2 = n_extra = n_steps = 0;
        AO2Moptimized(HEIGHT, WIDTH, R, n_cov2, n_extra, n_steps);

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
