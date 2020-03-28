#ifndef DUALIZATION_AO2
#define DUALIZATION_AO2

#include "AO2_trajectory.h"

//! Implementation of basic AO2 algorithm
/*! Cycles through all possible trajectories and collects coverages.*/
void AO2(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

//! Implementation of AO2 algorithm which deletes similar cols
void AO2_delete_similar_cols (coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_delete_similar_cols traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            // coverages.push_back(traj.get_coverage());
            n_cov += traj.get_n_coverages();
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
}

//! Implementation of AO2 algorithm which checks competing rows at each step
void AO2_stop_not_upper(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2Trajectory_stop_not_upper traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        traj.complete_trajectory();
        n_steps += traj.get_changes_size() - len_last;
        if (traj.check_upper()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        };
    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}


//! Implementation of AO2 algorithm with zero extra steps
void AO2Zero(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2ZeroTrajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

//! Implementation of AO2M algorithm with checking criterion on backtracing and no check_upper
void AO2Mplus(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2MplusTrajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

void AO2Mlightest(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2MlightestTrajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

void AO2LL(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2TrajectoryLL traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

void AO2Moptimized(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2MoptimizedTrajectory traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

void AO2Best(coord n, coord m, ull** R, uint64_t & n_cov, uint64_t& n_extra, uint64_t& n_steps) {
    uint64_t len_last = 0;
    CovCollector coverages;
    AO2TrajectoryBest traj(n, m, R);
    do {
        len_last = traj.get_changes_size();
        if (traj.complete_trajectory()) {
            n_cov++;
            // coverages.push_back(traj.get_coverage());
            // std::cout << "COVERAGE" << '\n';
            // for (auto q: traj.get_coverage()) std::cout << q << ' ';
            // std::cout << '\n';
        } else {
            n_extra++;
        }
        n_steps += traj.get_changes_size() - len_last;

    } while (traj.find_neighbour());
    // n_cov += coverages.size();
}

#endif
