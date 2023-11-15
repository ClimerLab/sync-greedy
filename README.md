# sync-greedy
Greedy implementation of the Sync feature selection algorithm.

## To Use
Configure the Makefile with the locaion of open mpi libraries amd binary

Compile with the Makefile by navigating to the root directory and entering: make

Update configuration file

Run the program. For an example enter: mpirun -np 4 ./sync <cfg_file>

## Configuration File
DATA_FILE - Tab seperated file where the first NUM_CASES columns are cases and the next NUM_CTRLS columns are controls. The row indicate features.
SCRATCH_DIR - Directory where ouputs are recorded

NUM_CASES - Number of cases in DATA_FILE
NUM_CTRLS - Numbr of controls in DATA_FILE
NUM_EXPRS - Number of features in DATA_FILE
NUM_HEAD_ROWS - Number of header rows in DATA_FILE
NUM_HEAD_COLS - Number of header columns in DATA_FILE

RISK - Boolean that indicates if risk patterns (true) or protective patterns should be found.
MAX_PS - Maximum pattern size to search.

MIN_OBJ - Minimum Younden J value that will be kept for PS=1 patterns
USE_SOL_POOL - Boolean that indicates if a solution pool will be kept.
SOL_POOL_SIZE - Size of solution pool to keep

MISSING_SYMBOL - String used to indicate missing data in DATA_FILE
SET_NA_TRUE - Boolean used to indicate if missing data is treated as both high and low
HIGH_VALUE - Value in DATA_FILE that indicates high expression
NORM_VALUE - Value in DATA_FILE that indicates normal expression
LOW_VALUE - Value in DATA_FILE that indicates low expression

## Outputs
PS#.solPool - Files containing a collection of patterns of size #
markerPairs.csv - Contains a count of individuals from G_1 that contain each pair of markers

## Notes
Requires Open MPI

DATA_FILE should be tab seperate, the columns represent individuals and the rows represent features