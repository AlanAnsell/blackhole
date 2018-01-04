#include "globals.h"

const char * ENGINE_NAME = "EventHorizon";
const char * VERSION_NUMBER = "2.4.1";

int parity[2] = {1, -1};
Value OFFSET[2] = {0, 1};

#ifdef FAST_
long long time_limit = 2000000;
long long time_left = 1950000;
#else
long long time_limit = 4900000;
long long time_left = time_limit;
#endif
long long time_started, time_ended;

long long get_time() {
	timeval current_time;
	gettimeofday(& current_time, NULL);
	return (long long)current_time.tv_sec * 1000000LL + (long long)current_time.tv_usec;
}

void send_move(char * move_str) {
	time_ended = get_time();
	time_left -= (time_ended - time_started);
	if (time_left < 0)
		time_left = 50000;
	fprintf(stderr, "%s\n", move_str);
	printf("%s\n", move_str);
	fflush(stdout);
}

void get_move(char * move_str) {
	assert(scanf("%s", move_str));
	time_started = get_time();
	if (! strcmp(move_str, "Quit\n"))
		exit(0);
}


U32 ROW[N_CELLS];
U32 NUM[N_CELLS];
U32 ADJ[N_CELLS][MAX_DEGREE];
U32 N_ADJ[N_CELLS];
U64 ADJ_MASK[N_CELLS];

Value STONE_POWER[2][1 << N_STONES][MAX_DEGREE+1];

U32 row_and_num_to_id(U32 row, U32 num) {
    return num + row * (row + 1) / 2;
}


void cell_id_to_name(U32 id, char * name) {
    name[0] = 'A' + (ROW[id] - NUM[id]);
    name[1] = '1' + NUM[id];
    name[2] = 0;
}


U32 cell_name_to_id(const char * name) {
    U32 num = name[1] - '1';
    U32 row = num + (name[0] - 'A');
    return row_and_num_to_id(row, num);
}


const U64 debruijn = 0x03f79d71b4cb0a89LL;

const U32 index64[64] = {
    0,  1, 48,  2, 57, 49, 28,  3,
   61, 58, 50, 42, 38, 29, 17,  4,
   62, 55, 59, 36, 53, 51, 43, 22,
   45, 39, 33, 30, 24, 18, 12,  5,
   63, 47, 56, 27, 60, 41, 37, 16,
   54, 35, 52, 21, 44, 32, 23, 11,
   46, 26, 40, 15, 34, 20, 31, 10,
   25, 14, 19,  9, 13,  8,  7,  6
};


void init() {
	fprintf(stderr, "R %s %s\n", ENGINE_NAME, VERSION_NUMBER);
#ifndef DEBUG_
    fflush(stderr);
    srand(time(NULL));
#endif

    U32 id = 0;
    U32 row, num;
    for (row = 0; row < N_ROWS; row++) {
        for (num = 0; num <= row; num++) {
            ROW[id] = row;
            NUM[id] = num;
            id++;
        }
    }

    for (id = 0; id < N_CELLS; id++) {
        row = ROW[id];
        num = NUM[id];
        if (num < row) {
            U32 next = row_and_num_to_id(row, num + 1);
            ADJ[id][N_ADJ[id]++] = next;
            ADJ[next][N_ADJ[next]++] = id;
        }
        if (row < N_ROWS - 1) {
            U32 left = row_and_num_to_id(row + 1, num);
            U32 right = row_and_num_to_id(row + 1, num + 1);
            ADJ[id][N_ADJ[id]++] = left;
            ADJ[id][N_ADJ[id]++] = right;
            ADJ[left][N_ADJ[left]++] = id;
            ADJ[right][N_ADJ[right]++] = id;
        }
    }
    
    for (id = 0; id < N_CELLS; id++) {
        ADJ_MASK[id] = 0;
        for (U32 j = 0; j < N_ADJ[id]; j++)
            ADJ_MASK[id] |= (1LL << (U64)ADJ[id][j]);
    }

    for (U32 p = 0; p < 2; p++) {
        Value m = parity[p];
        for (U32 mask = 0; mask < (1 << N_STONES); mask++) {
            Value * power = STONE_POWER[p][mask];
            power[0] = 0;
            U32 i = 1;
            for (int j = N_STONES - 1; j >= 0; j--) {
                if (mask & (1 << j)) {
                    power[i] = power[i - 1] + m * (j + 1);
                    i++;
                }
            }
        }
    }

}


const Real UCB_C = 1.2;

Real sigmoid(Real x) {
    return 1.0 / (1.0 + exp(-x));
}
