#include "GreedyWorker.h"
#include "Parallel.h"

GreedyWorker::GreedyWorker(const ConfigParser &_parser) : parser(&_parser),
                                                          data(*parser),
                                                          scratch_dir(parser->getString("SCRATCH_DIR")),
                                                          world_rank(Parallel::get_world_rank()),
                                                          tag(0),
                                                          start(0),
                                                          stop(0),
                                                          min_obj(0.0),
                                                          sol_pool(parser->getSizeT("SOL_POOL_SIZE")),
                                                          end_(false) {}

GreedyWorker::~GreedyWorker() {}

void GreedyWorker::receive_problem() {
  MPI_Status status;
  // Check if signal to end was received
  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  if (status.MPI_TAG == Parallel::CONVERGE_TAG) {
    tag = Parallel::CONVERGE_TAG;

    char signal;
    MPI_Recv(&signal, 1, MPI_CHAR, 0, Parallel::CONVERGE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    end_ = true;

  } else if (status.MPI_TAG == Parallel::PS1_TAG) {
    tag = Parallel::PS1_TAG;
    
    MPI_Recv(&start, 1, CUSTOM_SIZE_T, 0, Parallel::PS1_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&stop, 1, CUSTOM_SIZE_T, 0, Parallel::PS1_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&min_obj, 1, MPI_DOUBLE, 0, Parallel::PS1_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  } else if (status.MPI_TAG == Parallel::PS2_TAG) {
    tag = Parallel::PS2_TAG;

    MPI_Recv(&start, 1, CUSTOM_SIZE_T, 0, Parallel::PS2_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&min_obj, 1, MPI_DOUBLE, 0, Parallel::PS2_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  } else if (status.MPI_TAG == Parallel::GREEDY_TAG) {
    tag = Parallel::GREEDY_TAG;

    // Receive number of markers in solution
    std::size_t sol_size;
    MPI_Recv(&sol_size, 1, CUSTOM_SIZE_T, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    sol.resize(sol_size);

    // Receive solution
    MPI_Recv(&sol[0], sol_size, CUSTOM_SIZE_T, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Receive lb
    MPI_Recv(&min_obj, 1, MPI_DOUBLE, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  } else {
    fprintf(stderr, "Unknown tag\n");
    exit(EXIT_FAILURE);
  }
}

void GreedyWorker::send_back_solution() {
  // Send number of solutions in pool
  const std::size_t num_sols = sol_pool.size();
  MPI_Send(&num_sols, 1, CUSTOM_SIZE_T, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD);

  // Send each solution
  for (std::size_t i = 0; i < num_sols; ++i) {
    auto pair = sol_pool.get_obj_sol_pair(i);
    
    MPI_Send(&pair.second[0], pair.second.size(), CUSTOM_SIZE_T, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD);
    MPI_Send(&pair.first, 1, MPI_DOUBLE, 0, Parallel::GREEDY_TAG, MPI_COMM_WORLD);
  }
}

void GreedyWorker::calc_ps1() {
  std::vector<std::size_t> sol(1);

  for (std::size_t i = start; i <= stop; ++i) {
    sol[0] = i;
    double f1 = data.get_grp1_freq(sol);

    if (f1 >= min_obj) {
      double obj = f1 - data.get_grp2_freq(sol);

      if (obj >= min_obj) {
        sol_pool.add_solution(obj, sol);
      
        if (sol_pool.size() == sol_pool.get_max_size() &&
            sol_pool.get_min_obj() > min_obj) {
          min_obj =  sol_pool.get_min_obj();
        }
      }
    }
  }
}

void GreedyWorker::calc_ps2() {
  std::string pairs_file = scratch_dir + + "markerPairs_part" + std::to_string(world_rank) + ".csv";
  FILE *pair_count_stream = open_file(pairs_file);

  std::vector<std::size_t> count(data.get_num_bins() - 1 - start, 0);
  std::vector<std::size_t> red_g1, red_g2;

  // Initialize solution to 'start' marker
  std::vector<std::size_t> sol = {start};
  
  // Reduce the individuals from group1 to only those that contain the 'start' marker
  for (std::size_t j = data.get_grp1_start(); j <= data.get_grp1_stop(); ++j) {
    if (data.indiv_has_pat(j, sol)) {
      red_g1.push_back(j);
    }
  }
  // Reduce the individuals from group1 to only those that contain the 'start' marker
  for (std::size_t j = data.get_grp2_start(); j <= data.get_grp2_stop(); ++j) {
    if (data.indiv_has_pat(j, sol)) {
      red_g2.push_back(j);
    }
  }

  // Add dummy marker to solution
  sol.push_back(0);

  // Loop through all markers after 'start'
  for (std::size_t i = start + 1; i < data.get_num_bins(); ++i) {
    count[i - start - 1] = 0; // Initialze count to zero
    
    // Count the number of individuals 'reduced group 1' that contain the i-th marker
    for (auto j : red_g1) {
      if (data.get_bin(i,j)) {
        ++count[i - start - 1];
      }
    }

    double f1 = static_cast<double>(count[i - start - 1])  / data.get_num_grp1();

    if (f1 >= min_obj) {
      double f2 = 0.0;
      for (auto j : red_g2) {
        if (data.get_bin(i,j)) {
          ++f2;
        }
      }
      f2 /= data.get_num_grp2();

      double obj = f1 - f2;
      sol[1] = i;
      sol_pool.add_solution(obj, sol);

      if (sol_pool.size() == sol_pool.get_max_size() &&
        sol_pool.get_min_obj() > min_obj) {
        min_obj =  sol_pool.get_min_obj();
      }
    }    
  }
  record_pair_count(pair_count_stream, count);
  close_file(pair_count_stream);
}

void GreedyWorker::calc() {
  std::vector<std::size_t> red_g1, red_g2;
  // Reduce the individuals from group1 to only those that contain the pattern
  for (std::size_t j = data.get_grp1_start(); j <= data.get_grp1_stop(); ++j) {
    if (data.indiv_has_pat(j, sol)) {
      red_g1.push_back(j);
    }
  }

  // // Check if f1 values are being recorded or if the f1 value is >= min_obj
  if (static_cast<double>(red_g1.size()) / data.get_num_grp1() >= min_obj) {
    // Reduce the individuals from group2 to only those that contain the pattern
    for (std::size_t j = data.get_grp2_start(); j <= data.get_grp2_stop(); ++j) {
      if (data.indiv_has_pat(j, sol)) {
        red_g2.push_back(j);
      }
    }

    auto new_sol = sol;
    new_sol.push_back(0);

    // Loop through all bins
    for (std::size_t i = 0; i < data.get_num_bins(); ++i) {
      // check that i is not in the solution
      if (std::find(sol.begin(), sol.end(), i) == sol.end()) {
        double f1 = 0.0, f2 = 0.0;

        for (auto j : red_g1) {
          if (data.get_bin(i,j)) {
            ++f1;
          }
        }

        f1 /= data.get_num_grp1();

        if (f1 >= min_obj) {
          for (auto j : red_g2) {
            if (data.get_bin(i,j)) {
              ++f2;
            }
          }
          f2 /= data.get_num_grp2();

          double obj = f1 - f2;

          if (obj >= min_obj) {
            new_sol[new_sol.size()-1] = i;
            sol_pool.add_solution(obj, new_sol);

            if (sol_pool.size() == sol_pool.get_max_size() &&
                sol_pool.get_min_obj() > min_obj) {
              min_obj =  sol_pool.get_min_obj();
            }
          }
        }
      }
    }
  }
}

FILE* GreedyWorker::open_file(const std::string &file_name) const {
  FILE *stream;
  if ((stream = fopen(file_name.c_str(), "a+")) == nullptr) {
    fprintf(stderr, "ERROR - Could not open file %s\n", file_name.c_str());
    exit(1);
  }
  return stream;
}

void GreedyWorker::close_file(FILE *stream) const {
  fclose(stream);
}

void GreedyWorker::record_pair_count(FILE *stream, const std::vector<std::size_t> &count) const {
  for (std::size_t i = 0; i < count.size()-1; ++i) {
    fprintf(stream, "%lu,", count[i]);
  }
  fprintf(stream, "%lu\n", count[count.size()-1]);
}

bool GreedyWorker::end() const {
  return end_;
}

void GreedyWorker::work() {
  receive_problem();
  if (end_) {
    return;
  }

  sol_pool.clear();

  if (tag == Parallel::PS1_TAG) {
    calc_ps1();
  } else if (tag == Parallel::PS2_TAG) {
    calc_ps2();
  } else {
    calc();
  }

  send_back_solution();
}
