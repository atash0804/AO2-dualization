#ifndef DUALIZATION_AO2
#define DUALIZATION_AO2

#include "AO2_trajectory.h"

//! Implementation of basic AO2 algorithm
/*! Cycles through all possible trajectories and collects coverages.*/
void AO2(c_int n, c_int m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory traj(n, m, R);
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

//! Implementation of AO2 algorithm which deletes similar cols
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

#endif
