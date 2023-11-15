#ifndef GREEDY_CONTROLLER_H
#define GREEDY_CONTROLLER_H

#include <string>
#include <stack>
#include <set>
#include "ConfigParser.h"
#include "ExprsData.h"
#include "SolPool.h"

class GreedyController {
  private:
    const ConfigParser *parser;
    const ExprsData data;
    const std::size_t world_size;
    const std::string scratch_dir;
    const double min_obj;

    std::stack<int> available_workers;
    std::set<int> unavailable_workers;

    std::size_t ps;

    SolPool pool1;
    SolPool pool2;
    SolPool *cur_pool;
    SolPool *old_pool;
    double lb;
  
    void send_ps1_problem(const std::size_t start, const std::size_t stop);
    void send_ps2_problem(const std::size_t marker);
    void send_problem(const std::vector<std::size_t> &sol);

    void receive_completion();

    void combine_marker_pair_files();

  public:
    GreedyController(const ConfigParser &_parser);
    ~GreedyController();

    void set_ps(const std::size_t _ps);
    void solve_ps1();
    void solve_ps2();
    void solve();

    void signal_workers_to_end() const;
};

#endif