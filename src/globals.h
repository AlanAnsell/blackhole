#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <map>

#define N_CELLS 36
#define N_ROWS 8
#define RED 0
#define BLUE 1
#define NO_RESULT 400
#define MAX_DEGREE 6
#define N_STONES 15
#define N_BLOCKED_CELLS 5
#define MIN_RESULT -75
#define MAX_RESULT 75
#define NO_CELL 400
#define VALUE_INF 1000
#define TIME_ELAPSED 2

#define MASK(x) (1LL << (U64)(x))

typedef int Value;
typedef unsigned char U8;
typedef unsigned short int U16;
typedef unsigned int U32;
typedef unsigned long long U64;
typedef double Real;

extern const char * ENGINE_NAME;
extern const char * VERSION_NUMBER;

extern int PARITY[2];
extern Value OFFSET[2];

extern long long time_limit;
extern long long time_left;
extern long long time_started, time_ended;

long long get_time();

void send_move(char * move_str);

void get_move(char * move_str);

extern U32 ROW[N_CELLS];
extern U32 NUM[N_CELLS];
extern U32 ADJ[N_CELLS][MAX_DEGREE];
extern U32 N_ADJ[N_CELLS];
extern U64 ADJ_MASK[N_CELLS];

extern Value STONE_POWER[2][1 << N_STONES][MAX_DEGREE+1];

extern const U64 debruijn;
extern const U32 index64[64];
#define LSB(x) ((x) & -(x))
#define INDEX(x) (index64[((x) * debruijn) >> 58])


U32 row_and_num_to_id(U32 row, U32 num);

void cell_id_to_name(U32 id, char * name);

U32 cell_name_to_id(const char * name);

void init();


#define N_MCT_NODES 100000

//extern const Real UCB_C;
extern const Real TARGET_INCREMENT_THRESH[2];
extern const Real TARGET_DECREMENT_THRESH[2];
extern const Real CHOOSE_TARGET_THRESH[2];

#define MAX_SIMS 60000
extern Real SQRT[MAX_SIMS];
extern Real SQRT_LOG[MAX_SIMS];

extern bool ANALYSE_WITH_SOLVER;
extern bool ANALYSE;
extern bool PLAYOUT;
extern Value INIT_ALPHA;

#endif
