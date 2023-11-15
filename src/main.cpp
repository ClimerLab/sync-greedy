#include <cstdio>
#include <cstdlib>
#include "ConfigParser.h"
#include "Parallel.h"
#include "GreedyController.h"
#include "GreedyWorker.h"
#include "Timer.h"

int main(int argc, char *argv[]) {
  // MPI init
  MPI_Init(NULL, NULL);
  const int world_rank = Parallel::get_world_rank();
  const int world_size = Parallel::get_world_size();

  try {
    if (world_rank == 0) {
      // Check user inputs
      if (argc != 2) {
        fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
      } else if (world_size < 2) {
        fprintf(stderr, "world_size must be greater than 1.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    ConfigParser parser(argv[1]);

    switch (world_rank) {
      case 0: {
        Timer timer;        
        GreedyController controller(parser);
        fprintf(stderr, "Starting PS1\n");
        timer.start();
        controller.solve_ps1();
        timer.stop();
        fprintf(stderr, "PS1 took %lf\n", timer.elapsed_cpu_time());

        fprintf(stderr, "Starting PS2\n");
        timer.restart();
        controller.solve_ps2();
        timer.stop();
        fprintf(stderr, "PS2 took %lf\n", timer.elapsed_cpu_time());

        for (std::size_t PS=3; PS <= parser.getSizeT("MAX_PS"); ++PS) {
          fprintf(stderr, "Starting PS%lu\n", PS);
          controller.set_ps(PS);
          timer.restart();
          controller.solve();
          timer.stop();

          fprintf(stderr, "PS%lu took %lf\n\n", PS, timer.elapsed_cpu_time());
        }

        controller.signal_workers_to_end();
        break;
      }
      
      default: {
        GreedyWorker worker(parser);
        while (!worker.end()) {
          worker.work();
        }
        break;
      }
    }
  } catch (std::exception &e) {
    fprintf(stderr, "  *** Fatal error reported by rank_%d: %s *** \n", world_rank, e.what());
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Finalize();

  return 0;
}