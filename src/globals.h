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

typedef int Value;
typedef size_t CellID;
typedef double Real;

extern const char * ENGINE_NAME;
extern const char * VERSION_NUMBER;

extern int parity[2];

extern long long time_limit;
extern long long time_left;
extern long long time_started, time_ended;

long long get_time();

void send_move(char * move_str);

void get_move(char * move_str);

extern size_t ROW[N_CELLS];
extern size_t NUM[N_CELLS];
extern CellID ADJ[N_CELLS][MAX_DEGREE];
extern size_t N_ADJ[N_CELLS];


CellID row_and_num_to_id(size_t row, size_t num);

void cell_id_to_name(CellID id, char * name);

CellID cell_name_to_id(const char * name);

void init();

#endif
