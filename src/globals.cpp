#include "globals.h"

const char * ENGINE_NAME = "Blackhole";
const char * VERSION_NUMBER = "0.0.0";

int parity[2] = {1, -1};
Value OFFSET[2] = {0, 1};

long long time_limit = 4900000;
long long time_left = time_limit;
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


size_t ROW[N_CELLS];
size_t NUM[N_CELLS];
CellID ADJ[N_CELLS][MAX_DEGREE];
size_t N_ADJ[N_CELLS];


CellID row_and_num_to_id(size_t row, size_t num) {
    return num + row * (row + 1) / 2;
}


void cell_id_to_name(CellID id, char * name) {
    name[0] = 'A' + (ROW[id] - NUM[id]);
    name[1] = '1' + NUM[id];
    name[2] = 0;
}


CellID cell_name_to_id(const char * name) {
    size_t num = name[1] - '1';
    size_t row = num + (name[0] - 'A');
    return row_and_num_to_id(row, num);
}


void init() {
	fprintf(stderr, "R %s %s\n", ENGINE_NAME, VERSION_NUMBER);
#ifndef DEBUG_
    srand(time(NULL));
#endif

    CellID id = 0;
    size_t row, num;
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
            CellID next = row_and_num_to_id(row, num + 1);
            ADJ[id][N_ADJ[id]++] = next;
            ADJ[next][N_ADJ[next]++] = id;
        }
        if (row < N_ROWS - 1) {
            CellID left = row_and_num_to_id(row + 1, num);
            CellID right = row_and_num_to_id(row + 1, num + 1);
            ADJ[id][N_ADJ[id]++] = left;
            ADJ[id][N_ADJ[id]++] = right;
            ADJ[left][N_ADJ[left]++] = id;
            ADJ[right][N_ADJ[right]++] = id;
        }
    }

}


const Real UCB_C = 1.6;

Real sigmoid(Real x) {
    return 1.0 / (1.0 + exp(-x));
}
