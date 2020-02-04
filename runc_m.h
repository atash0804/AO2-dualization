#ifndef DUALIZATION_RUNC_M
#define DUALIZATION_RUNC_M

#include <map>
#include <set>
#include <ctime>
#include <cerrno>
#include <string>
#include <random>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

typedef uint64_t ull;
size_t cov_count = 0;

#ifndef DUALIZATION_CHUNK_SIZE
#define DUALIZATION_CHUNK_SIZE
enum {
    CHUNK_SIZE = sizeof(ull) * 8,
};
#endif

class customset {
    ull *data;
public:
    size_t sz;
    size_t chunks;

    explicit customset(size_t size) : sz(size) {
        chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull *) calloc(chunks, sizeof(ull));
    }

    customset(set<size_t> source, size_t size) : sz(size) {
        chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull *) calloc(chunks, sizeof(ull));
        for (auto entry: source) {
            this->sett(entry);
        }
    }

    customset(const customset &source) : sz(source.sz), chunks(source.chunks) {
        data = source.data;
    }

    customset &operator=(const customset &source) {
        data = source.data;
        sz = source.sz;
        return *this;
    }

//    ~customset() {
//        if (!used)
//            free(data);
//    }

    bool operator<(const customset &b) const {
        for (int i = 0; i < chunks; i++) {
            if (data[i] != b.data[i]) return (data[i] < b.data[i]);
        }
        return false;
    }

    void sett(ull k) {
        data[k / CHUNK_SIZE] |= ((ull) 1) << (k % CHUNK_SIZE);
    }

    void clear(ull k) {
        data[k / CHUNK_SIZE] &= ~((ull) 1 << (k % CHUNK_SIZE));
    }

    bool in(ull k) const {
        return bool((ull) 1 & (data[k / CHUNK_SIZE] >> (k % CHUNK_SIZE)));
    }

    friend bool check_intersection(customset s1, customset s2);
};

bool check_intersection(customset s1, customset s2) {
    if (s1.sz != s2.sz) {
        cerr << "Sets should be of equal size" << endl;
    }
    size_t chunks = (s1.sz + CHUNK_SIZE - 1) / CHUNK_SIZE;
    for (int i = 0; i < chunks; i++) {
        if (s1.data[i] & s2.data[i]) return false;
    }
    return true;
}

set<customset> default_coverage;
set<pair<set<size_t>, set<size_t>>> default_found_coverages;

/**
 * A simple class for binary matrix stored as bit sets
 */
class BitMatrix {
    /*! The matrix where all data is stored*/
    vector<vector<ull>> matrix;
    size_t height;
    size_t width;
    size_t chunks;

public:
    BitMatrix() : height(0), width(0), chunks(0) {}

    BitMatrix(istream &in, size_t n, size_t m) : height(n), width(m) {
        ull bit;
        chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

        for (int i = 0; i < height; i++) {
            vector<ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (int k = 1;
                     (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: "
                             << bit << endl;
                        throw out_of_range("");
                    }
                    chunk |= (bit << (CHUNK_SIZE - k));
                }
                row.push_back(chunk);
            }
            matrix.push_back(row);
        }
    }

    BitMatrix(const string &filename, size_t n, size_t m) : height(n),
                                                            width(m) {
        ifstream in;
        in.open(filename);
        if (!in.is_open()) {
            cerr << "Error: " << strerror(errno) << endl;
            cerr << "Failed to open input file" << endl;
            throw bad_exception();
        }
        ull bit;
        chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

        for (int i = 0; i < height; i++) {
            vector<ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (size_t k = 1;
                     (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: "
                             << bit << endl;
                        throw out_of_range("");
                    }
                    chunk |= (bit << (CHUNK_SIZE - k));
                }
                row.push_back(chunk);
            }
            matrix.push_back(row);
        }
        in.close();
    }

    ~BitMatrix() {
        for (auto a: matrix) a.clear();
        matrix.clear();
    }

    size_t getHeight() const {
        return height;
    }

    size_t getWidth() const {
        return width;
    }

    size_t getChunks() const {
        return chunks;
    }

    const vector<vector<ull>> &getMatrix() const {
        return matrix;
    }

    int at(size_t i, size_t j) const {
        size_t shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
        return int((matrix[i][j / CHUNK_SIZE] & (1ULL << shift)) >> shift);
    }

    friend ostream &operator<<(ostream &os, const BitMatrix &bm);
};

ostream &operator<<(ostream &os, const BitMatrix &bm) {
    cout << "Bit matrix with shape (" << bm.height << ", " << bm.width << ")"
         << endl << endl;
    for (size_t i = 0; i < bm.height; i++) {
        for (size_t j = 0; j < bm.width; j++) {
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << endl;
    return os;
}

class PartialBitMatrix : public BitMatrix {
    set<size_t> available_rows;
    set<size_t> available_cols;
    set<size_t> selected_cols;
    size_t cur_height;
    size_t cur_width;

    void delete_zero_columns() {
        vector<ull> rows_disjunction;

        for (size_t j = 0; j < this->getChunks(); j++)
            rows_disjunction.push_back(0);

        for (size_t i: available_rows) {
            for (size_t j = 0; j < this->getChunks(); j++) {
                rows_disjunction[j] |= this->getMatrix()[i][j];
            }
        }
        size_t shift;

        set<size_t> zero_cols;
        for (size_t j: available_cols) {
            shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
            if ((rows_disjunction[j / CHUNK_SIZE] & (1ULL << shift)) >> shift ==
                0) {
                zero_cols.insert(j);
            }
        }
        for (size_t j: zero_cols) {
            // cout << "Should delete zero col: " << j << endl << flush;
            delete_column(j, false);
        }
    }

    void delete_wider_rows() {
        bool fl; //fl is true if line1 >= line2
        bool eq; //eq is true if line1 == line2
        size_t shift;
        set<size_t> wider_rows;
        for (size_t i: available_rows) {
            for (size_t j: available_rows) {
                if (i == j) continue;
                fl = true;
                eq = true;
                for (size_t k: available_cols) {
                    shift = CHUNK_SIZE - 1 - k % CHUNK_SIZE;
                    if ((this->getMatrix()[i][k / CHUNK_SIZE] & (1ULL << shift))
                        != (this->getMatrix()[j][k / CHUNK_SIZE] &
                           (1ULL << shift))) {
                        eq = false;
                    }
                    if ((this->getMatrix()[i][k / CHUNK_SIZE] & (1ULL << shift))
                        < (this->getMatrix()[j][k / CHUNK_SIZE] &
                           (1ULL << shift))) {
                        fl = false;
                        break;
                    }
                }
                if (fl) {
                    if (eq && (i < j)) {
                        wider_rows.insert(i);
                    } else {
                        wider_rows.insert(i);
                    }
                    break;
                }
            }
        }
        for (size_t i: wider_rows) {
            // cout << "Should delete wider row: " << i << endl << flush;
            delete_row(i, false);
        }
    }

public:
    PartialBitMatrix(istream &in, size_t n, size_t m) : BitMatrix(in, n, m) {
        for (size_t i = 0; i < n; i++) available_rows.insert(i);
        for (size_t i = 0; i < m; i++) available_cols.insert(i);
        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix(const string &filename, size_t n, size_t m) : BitMatrix(
            filename, n, m) {
        for (size_t i = 0; i < n; i++) available_rows.insert(i);
        for (size_t i = 0; i < m; i++) available_cols.insert(i);

        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix() : cur_width(0), cur_height(0) {}

    void delete_column(size_t col, bool outside = true) {
        // cout << "BEFORE DELETING COLUMN " << col << endl << *this << endl;
        if (col >= this->getWidth()) {
            cerr << "Invalid column number" << endl;
            throw out_of_range("");
        }
        if (cur_width <= 0)
            cur_height = 0;
        else if (available_cols.find(col) != available_cols.end()) cur_width--;

        available_cols.erase(col);
        if (outside) selected_cols.insert(col);
    }

/**
 * @param row
 * @return
 */
    void delete_row(size_t row, bool outside = true) {
        // cout << "BEFORE DELETING ROW " << row << endl << *this << endl;
        if (row >= this->getHeight()) {
            cerr << "Invalid row number" << endl;
            throw out_of_range("");
        }
        if (cur_height <= 0)
            cur_height = 0;
        else if (available_rows.find(row) != available_rows.end()) cur_height--;

        available_rows.erase(row);
    }

    void update_matrix() {
        delete_zero_columns();
        delete_wider_rows();
    }

    size_t getCur_height() const {
        return cur_height;
    }

    size_t getCur_width() const {
        return cur_width;
    }

    const set<size_t> &getAvailable_rows() const {
        return available_rows;
    }

    const set<size_t> &getAvailable_cols() const {
        return available_cols;
    }

    const set<size_t> &getSelected_cols() const {
        return selected_cols;
    }

    pair<size_t, size_t> getLightestRow() const {
        ull tmp;
        size_t count, min = getHeight() + 1, min_n = 0;
        for (auto row: available_rows) {
            count = 0;
            for (int i = 0; i < getChunks(); i++) {
                tmp = getMatrix()[row][i];
                while (tmp != 0) {
                    tmp = tmp & (tmp - 1);
                    count++;
                }
            }
            if (count < min) {
                min = count;
                min_n = row;
            }
        }
        return {min_n, min};
    }

    friend ostream &operator<<(ostream &os, const PartialBitMatrix &bm);
};

ostream &operator<<(ostream &os, const PartialBitMatrix &bm) {
    cout << "Partial bit matrix with initial shape (" << bm.getHeight() << ", "
         << bm.getWidth() << ")" << endl;
    cout << "Current shape (" << bm.cur_height << ", " << bm.cur_width << ")"
         << endl << endl;
    for (auto i: bm.available_rows) {
        for (auto j: bm.available_cols) {
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << "Available cols: [";
    for (auto entry: bm.getAvailable_cols()) cout << entry << ",";
    cout << ']' << endl;
    cout << "Available rows: [";
    for (auto entry: bm.getAvailable_rows()) cout << entry << ",";
    cout << ']' << endl;
    cout << "Selected cols: [";
    for (auto entry: bm.getSelected_cols()) cout << entry << ",";
    cout << ']' << endl;
    cout << endl;
    return os;
}

bool check_support_rows(PartialBitMatrix &L,
                        map<size_t, set<size_t>> &supporting_rows, size_t col) {
    set<size_t> one_rows;
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.insert(i);
    }
    for (auto &entry: supporting_rows) {
        if (includes(one_rows.begin(), one_rows.end(), entry.second.begin(),
                     entry.second.end()))
            return false;
    }
    return true;
}

map<size_t, set<size_t>> update_support_rows(PartialBitMatrix &L,
                                             map<size_t, set<size_t>> &supporting_rows,
                                             size_t col) {
    set<size_t> one_rows;
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.insert(i);
    }
    map<size_t, set<size_t>> copy = supporting_rows;
    for (auto &entry: copy) {
        set<size_t> should_erase;
        for (auto row: entry.second) {
            if (one_rows.find(row) != one_rows.end()) {
                should_erase.insert(row);
                one_rows.erase(row);
            }
        }
        for (size_t row: should_erase) {
            entry.second.erase(row);
        }
    }
    copy[col] = one_rows;
    return copy;
}

void print_results(set<pair<set<size_t>, set<size_t>>> &found_coverages) {
    for (auto &cov_pair: found_coverages) {
        printf("{");
        for (auto entry: cov_pair.first) {
            printf("%ld ", entry);
        }
        printf("}  {");
        for (auto entry: cov_pair.second) {
            printf("%ld ", entry);
        }
        printf("}\n");
    }
}

void
dualization(PartialBitMatrix &L1, map<size_t, set<size_t>> &supporting_rows1,
            bool weights = false, \
    bool save = false, set<customset> &coverages = default_coverage) {
    // cout << "MATRIX:" << L1 << endl << endl << flush;
    PartialBitMatrix L_new(L1);
    map<size_t, set<size_t>> supporting_rows_new;
    bool L1_empty;
    size_t row_number;

    L1_empty = L1.getCur_height() == 0;
    if (L1_empty) {
        if (save) {
            coverages.insert(customset(L1.getSelected_cols(), L1.getWidth()));
        } else {
            printf("{");
            for (auto entry: L1.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}\n");
        }
        return;
    }
    if (weights) {
        row_number = L1.getLightestRow().first;
    } else {
        row_number = *L1.getAvailable_rows().begin();
    }
    for (size_t col: L1.getAvailable_cols()) {
        if (L1.at(row_number, col)) {
            if (check_support_rows(L1, supporting_rows1, col)) {
                supporting_rows_new = update_support_rows(L1, supporting_rows1,
                                                          col);
                L_new = L1;
                for (size_t i = 0; i < L1.getHeight(); i++) {
                    if (L1.at(i, col)) L_new.delete_row(i);
                }
                L_new.delete_column(col);
                L_new.update_matrix();
                dualization(L_new, supporting_rows_new, weights, save,
                            coverages);
            }
        }
    }
}
#endif
