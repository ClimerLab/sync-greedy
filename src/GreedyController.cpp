#include "GreedyController.h"

#include "Parallel.h"
#include <experimental/filesystem>

GreedyController::GreedyController(const ConfigParser &_parser) : parser(&_parser),
                                                                  data(*parser),
                                                                  world_size(Parallel::get_world_size()),
                                                                  scratch_dir(parser->getString("SCRATCH_DIR")),
                                                                  min_obj(parser->getDouble("MIN_OBJ")),
                                                                  ps(0),
                                                                  pool1(parser->getSizeT("SOL_POOL_SIZE")),
                                                                  pool2(parser->getSizeT("SOL_POOL_SIZE")),
                                                                  cur_pool(&pool1),
                                                                  old_pool(&pool2),
                                                                  lb(min_obj) {
  for (std::size_t i = 1; i < world_size; ++i) {
    available_workers.push(i);
  }
}

GreedyController::~GreedyController() {}

void GreedyController::send_ps1_problem(const std::size_t start, const std::size_t stop) {  
  // Get next worker
  const int worker = available_workers.top();

  // Send start and stop marker values
  MPI_Send(&start, 1, MPI_INT, worker, Parallel::PS1_TAG, MPI_COMM_WORLD);
  MPI_Send(&stop, 1, MPI_INT, worker, Parallel::PS1_TAG, MPI_COMM_WORLD);

  MPI_Send(&min_obj, 1, MPI_DOUBLE, worker, Parallel::PS1_TAG, MPI_COMM_WORLD);

  // Make worker unavailable
  available_workers.pop();
  unavailable_workers.insert(worker);
}

void GreedyController::send_ps2_problem(const std::size_t marker) {
  while (available_workers.empty()) {
    receive_completion();
  }

  // Get next worker
  const int worker = available_workers.top();

  // Send marker
  MPI_Send(&marker, 1, MPI_INT, worker, Parallel::PS2_TAG, MPI_COMM_WORLD);

  // Send min_obj
  MPI_Send(&lb, 1, MPI_DOUBLE, worker, Parallel::PS2_TAG, MPI_COMM_WORLD);

  // Make worker unavailable
  available_workers.pop();
  unavailable_workers.insert(worker);
}

void GreedyController::send_problem(const std::vector<std::size_t> &sol) {
  while (available_workers.empty()) {
    receive_completion();
  }

  // Get next worker
  const int worker = available_workers.top();

  // Send number of markers in solution
  const std::size_t sol_size = sol.size();
  MPI_Send(&sol_size, 1, CUSTOM_SIZE_T, worker, Parallel::GREEDY_TAG, MPI_COMM_WORLD);

  // Send solution
  MPI_Send(&sol[0], sol_size, CUSTOM_SIZE_T, worker, Parallel::GREEDY_TAG, MPI_COMM_WORLD);

  // Send lower bound
  MPI_Send(&lb, 1, MPI_DOUBLE, worker, Parallel::GREEDY_TAG, MPI_COMM_WORLD);

  // Make worker unavailable
  available_workers.pop();
  unavailable_workers.insert(worker);
}

void GreedyController::receive_completion() {
  assert(available_workers.size() < world_size - 1); // Cannot receive problem when no workers are working

  MPI_Status status;

  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  std::size_t num_sols;
  std::vector<std::size_t> sol(ps);
  double obj_value;

  // Receive number of solutions in pool
  MPI_Recv(&num_sols, 1, CUSTOM_SIZE_T, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

  for (std::size_t i = 0; i < num_sols; ++i) {
    // Receive markes in solution
    MPI_Recv(&sol[0], ps, CUSTOM_SIZE_T, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Receive solution objective value
    MPI_Recv(&obj_value, 1, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    cur_pool->add_solution(obj_value, sol);

    if (cur_pool->size() == cur_pool->get_max_size() && cur_pool->get_min_obj() > lb) {
      lb = cur_pool->get_min_obj();
    }
  }

  available_workers.push(status.MPI_SOURCE);
  unavailable_workers.erase(status.MPI_SOURCE);
}

void GreedyController::combine_marker_pair_files() {
  std::string file_name = scratch_dir + "markerPairs.csv";
  std::ofstream output(file_name, std::ios_base::binary);

  for (std::size_t p = 1; p < world_size; ++p) {
    std::string input_file = scratch_dir + "markerPairs_part" + std::to_string(p) + ".csv";
    std::ifstream input(input_file.c_str(), std::ios_base::binary);
    output << input.rdbuf();

    remove(input_file.c_str());
  }
}

void GreedyController::set_ps(const std::size_t _ps) {
  ps = _ps;
}

void GreedyController::solve_ps1() {
  set_ps(1);

  std::size_t delta = data.get_num_bins() / (world_size-1);

  for (std::size_t i = 1; i < world_size; ++i) {
    std::size_t start = delta * (i-1);
    std::size_t stop = i * delta - 1;
    if (data.get_num_bins() - 1 - stop < delta) {
      stop = data.get_num_bins() - 1;
    }

    send_ps1_problem(start, stop);
  }

  while (!unavailable_workers.empty()) {
    receive_completion();
  }

  std::string file_name = scratch_dir + "ps1.solPool";
  cur_pool->write_to_file(file_name, data);
}

void GreedyController::solve_ps2() {
  set_ps(2);

  cur_pool = &pool2;
  old_pool = &pool1;

  lb = min_obj;

  for (std::size_t i = 0; i < data.get_num_bins()-1; ++i) {
    send_ps2_problem(i);
  }

  while (!unavailable_workers.empty()) {
    receive_completion();
  }

  std::string file_name = scratch_dir + "ps2.solPool";
  cur_pool->write_to_file(file_name, data);

  fprintf(stderr, "Greedy max obj for PS=%lu: %lf\n", ps, cur_pool->get_max_obj());
  fprintf(stderr, "Greedy min obj for PS=%lu: %lf\n", ps, cur_pool->get_min_obj());

  combine_marker_pair_files();
}

void GreedyController::solve() {
  SolPool *tmp = cur_pool;
  cur_pool = old_pool;
  old_pool = tmp;
  cur_pool->clear();

  lb = min_obj;
  for (std::size_t i = 0; i < old_pool->size(); ++i) {
    send_problem(old_pool->get_sol(i));
  }

  while (!unavailable_workers.empty()) {
    receive_completion();
  }

  std::string file_name = scratch_dir + "ps" + std::to_string(ps) + ".solPool";
  cur_pool->write_to_file(file_name, data);

  fprintf(stderr, "Greedy max obj for PS=%lu: %lf\n", ps, cur_pool->get_max_obj());
  fprintf(stderr, "Greedy min obj for PS=%lu: %lf\n", ps, cur_pool->get_min_obj());
}

void GreedyController::signal_workers_to_end() const {
  char signal = 0;
  for (std::size_t i = 1; i < world_size; ++i) {
    MPI_Send(&signal, 1, MPI_CHAR, i, Parallel::CONVERGE_TAG, MPI_COMM_WORLD);
  }
}