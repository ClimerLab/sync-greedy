#ifndef GREEDY_WORKER_H
#define GREEDY_WORKER_H

#include <vector>
#include <string>
#include "ConfigParser.h"
#include "ExprsData.h"
#include "SolPool.h"

class GreedyWorker {
  private:
    const ConfigParser *parser;
    const ExprsData data;
    const std::string scratch_dir;
    const std::size_t world_rank;
    int tag;

    std::size_t start;
    std::size_t stop;
    double min_obj;
    std::vector<std::size_t> sol;

    SolPool sol_pool;
  
    bool end_;

    void receive_problem();
    void send_back_solution();

    void calc_ps1();
    void calc_ps2();
    void calc();

    FILE* open_file(const std::string &file_name) const;
    void close_file(FILE *stream) const;
    void record_pair_count(FILE *stream, const std::vector<std::size_t> &count) const;

  public:
    GreedyWorker(const ConfigParser &_parser);
    ~GreedyWorker();

    bool end() const;

    void work();
};

#endif