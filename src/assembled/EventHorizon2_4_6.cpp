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
#define TIME_ELAPSED 2

#define MASK(x) (1LL << (U64)(x))

typedef int Value;
typedef unsigned int U32;
typedef unsigned long long U64;
typedef double Real;

#define LSB(x) ((x) & -(x))
#define INDEX(x) (index64[((x) * debruijn) >> 58])

#define N_MCT_NODES 120000

const char * ENGINE_NAME = "EventHorizon";
const char * VERSION_NUMBER = "2.4.6";

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


U32 ROW[N_CELLS];
U32 NUM[N_CELLS];
U32 ADJ[N_CELLS][MAX_DEGREE];
U32 N_ADJ[N_CELLS];
U64 ADJ_MASK[N_CELLS];

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
#ifdef DEBUG_
    fprintf(stderr, "DEBUG mode\n");
#endif
#ifndef DEBUG_
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
}


const Real UCB_C = 1.2;

Real sigmoid(Real x) {
    return 1.0 / (1.0 + exp(-x));
}

struct AMAFRecord {
    U32 n_wins_;
    U32 n_playouts_;

    AMAFRecord(): n_wins_(0), n_playouts_(0) {}

};


class AMAFTable {
public:
    U32 n_stones_;
    std::vector<std::vector<AMAFRecord>> amaf_;

    void init(U32 n_open, U32 n_stones) {
        n_stones_ = n_stones;
        amaf_ = std::vector<std::vector<AMAFRecord>>(n_open, std::vector<AMAFRecord>());
    }

    U32 get_n_playouts(U32 cell_index, U32 stone_index) {
        if (amaf_[cell_index].empty())
            return 0;
        return amaf_[cell_index][stone_index].n_playouts_;
    }

    Real get_value(U32 cell_index, U32 stone_index) {
#ifdef DEBUG_
        assert(! amaf_[cell_index].empty());
#endif
        AMAFRecord& record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
        assert(record.n_playouts_ != 0);
#endif
        return (Real)record.n_wins_ / (Real)record.n_playouts_;
    }

    void update(U32 cell_index, U32 stone_index, bool win) {
        std::vector<AMAFRecord>& records = amaf_[cell_index];
        if (records.empty())
            records = std::vector<AMAFRecord>(n_stones_, AMAFRecord());
        AMAFRecord& record = records[stone_index];
        record.n_wins_ += win;
        record.n_playouts_++;
    }

};

class List {
public:
    U32 * val_;
    U32 * loc_;
    U32 len_;

    List(): val_(NULL), loc_(NULL), len_(0) {}

    void init(U32 max_len, U32 range) {
        len_ = 0;
        val_ = new U32[max_len];
        loc_ = new U32[range];
    }

    ~List() {
        delete [] val_;
        delete [] loc_;
    }

    inline void add(U32 val) {
        loc_[val] = len_;
        val_[len_] = val;
        len_++;
    }

    inline void remove(U32 val) {
        len_--;
        if (len_) {
            U32 end_val = val_[len_];
            U32 val_loc = loc_[val];
            loc_[end_val] = val_loc;
            val_[val_loc] = end_val;
        }
    }

    inline void readd(U32 val) {
        U32 old_loc = loc_[val];
        U32 new_val = val_[old_loc];
        loc_[new_val] = len_;
        val_[len_] = new_val;
        val_[old_loc] = val;
        len_++;
    }

};


class Move {
public:
    U32 cell_;
    Value stone_value_;

    Move(): cell_(0), stone_value_(0) {}

    Move(U32 cell, Value stone_value): cell_(cell), stone_value_(stone_value) {}

    Move(char * move_str, U32 turn);

    void to_str(char * str);

    bool operator < (const Move& other) const;
};


struct Snapshot {
    U32 n_open_;
    U64 open_mask_;
    U32 open_[N_CELLS];
    U32 adj_[N_CELLS][MAX_DEGREE];
    U32 n_adj_[N_CELLS];
    Value value_[N_CELLS];

    Value stones_[2][N_STONES];
    U32 n_stones_[2];
    U32 turn_;

    U64 controls_;
    U32 n_controls_[2];

    U64 dead_[2];
    U32 n_dead_[2];

    U64 stale_[2];
    U32 n_stale_[2];
};

class HashInfo;

#define HASH_COEFFICIENT 6787142722019916807ULL
#define HASH_SIZE_BASE 15
#define HASH_L (64 - HASH_SIZE_BASE)
#define HASH_SIZE (1 << HASH_SIZE_BASE)
#define BUCKET_SIZE 10
#define NO_HASH_ENTRY 2
// High Speed Hashing for Integers and Strings
// Mikkel Thorup July 15, 2014
#define HASH_BITMASKS(x, y) ((HASH_COEFFICIENT * (x) * (y)) >> HASH_L)


class HashRecord {
public:
    U64 red_;
    U64 blue_;
    U32 result_;
};


class HashTable {
public:
    U32 n_records_[HASH_SIZE];
    HashRecord table_[HASH_SIZE][BUCKET_SIZE];
    
    void init() {
        for (U32 i = 0; i < HASH_SIZE; i++)
            n_records_[i] = 0;
    }

    void add(U64 red, U64 blue, U32 result) {
        U64 hash = HASH_BITMASKS(red, blue);
        //fprintf(stderr, "Hash: %llu\n", hash);
        U32& n_records = n_records_[hash];
        if (n_records < BUCKET_SIZE - 1) {
            HashRecord& hash_record = table_[hash][n_records];
            hash_record.red_ = red;
            hash_record.blue_ = blue;
            hash_record.result_ = result;
            n_records++;
        }
#ifdef DEBUG_
        else {
            //fprintf(stderr, "Hash bucket %llu overflowed\n", hash);
        }
#endif
    }

    U32 find(U64 red, U64 blue) {
        U64 hash = HASH_BITMASKS(red, blue);
        U32 n_records = n_records_[hash];
        HashRecord * records = table_[hash];
        for (U32 i = 0; i < n_records; i++) {
            HashRecord& record = records[i];
            if (record.red_ == red && record.blue_ == blue)
                return record.result_;
        }
        return NO_HASH_ENTRY;
    }

};

class Position {
public:
    List open_;
    U64 open_mask_;
    List adj_[N_CELLS];
    Value value_[N_CELLS];

	Value stones_[2][N_STONES];
	U32 stone_loc_[2][N_STONES+1];
    U32 n_stones_[2];
	U32 turn_;
    int m_;

    Value alpha_;
    U64 controls_;
    U32 n_controls_[2];

    U64 dead_[2];
    U32 n_dead_[2];

    U64 stale_[2];
    U32 n_stale_[2];

    Snapshot snapshot_;


	Position(U32 blocked[N_BLOCKED_CELLS], Value alpha);
   
    void fill(U32 cell_id, Value stone_value);

    void unfill(U32 cell_id, Value stone_value);

    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    void get_stone_power(U32 p, Value * power);

    Move get_random_move();

    Move get_default_policy_move();

    Move get_expectation_maximising_move();

	Real calculate_win_prob(Value alpha) const;
  
    std::pair<Real, Move> get_best_alpha_move(Value alpha);

    std::pair<Real, Move> get_best_move();

    Real calculate_expectation() const;		

    std::pair<Real, Move> search_expectation(U32 depth, Real a, Real b);

    Real get_control_heuristic() const;

    void get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves);

    void get_all_moves(std::vector<Move>& moves);

    void get_all_reasonable_moves(std::vector<Move>& moves);

    U32 solve(long long end_time, U32& counter, U32& hash_hits,
            HashTable& table, HashInfo& hash_info, U64 stone_masks[2]);

    std::pair<U32, Move> get_optimal_move(
            long long end_time, U32& counter, U32& hash_hits, bool break_ties, HashTable& table);

    bool is_dead(U32 cell_id, U32 p, Value * other_power);

    bool all_adj_dead(U32 cell_id);

    void find_stale_cells();

    bool is_winning(U32 p) const;

    void set_alpha(Value alpha);

    void take_snapshot();

    void restore_snapshot();

    std::pair<U32, U32> get_cell_and_stone_indices(const Move& move);

    void print(FILE * f);    
};

class HashInfo {
public:
    U32 cell_index_[N_CELLS];
    U32 stone_shift_[2][N_STONES+1];

    HashInfo(const Position& pos) {
        U32 i;
        for (i = 0; i < pos.open_.len_; i++)
            cell_index_[pos.open_.val_[i]] = i + 1;
        for (U32 p = 0; p < 2; p++) {
            U32 n_stones = pos.n_stones_[p];
            const Value * stones = pos.stones_[p];
            for (U32 i = 0; i < n_stones; i++)
                stone_shift_[p][stones[i]] = 5 * i;
        }
    }

};

Move::Move(char * move_str, U32 turn) {
    U32 stone_number;
    char cell_str[5];
    sscanf(move_str, "%2s=%u", cell_str, &stone_number);
    cell_ = cell_name_to_id(cell_str);
    stone_value_ = parity[turn] * stone_number;
}


void Move::to_str(char * str) {
    char cell_str[5];
    cell_id_to_name(cell_, cell_str);
    sprintf(str, "%s=%d", cell_str, abs(stone_value_));
}


bool Move::operator < (const Move& other) const {
    if (cell_ != other.cell_)
        return cell_ < other.cell_;
    return stone_value_ < other.stone_value_;
}



Position::Position(U32 blocked[N_BLOCKED_CELLS], Value alpha): turn_(RED), m_(parity[RED]) {
    U32 i, j;
    open_.init(N_CELLS, N_CELLS);
    open_mask_ = 0;
    for (i = 0; i < N_CELLS; i++) {
        value_[i] = 0;
        U32 degree = N_ADJ[i];
        List& adj = adj_[i];
        adj.init(MAX_DEGREE, N_CELLS);
        for (j = 0; j < degree; j++)
            adj.add(ADJ[i][j]);
    }
    bool is_blocked[N_CELLS] = {false};
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        is_blocked[blocked[i]] = true;
    for (i = 0; i < N_CELLS; i++) {
        if (! is_blocked[i]) {
            open_mask_ |= MASK(i);
            open_.add(i);
        }
    }
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        fill(blocked[i], 0);

    for (U32 p = 0; p < 2; p++) {
        n_stones_[p] = N_STONES;
        for (i = 0; i < N_STONES; i++) {
            stones_[p][i] = i + 1;
            stone_loc_[p][i+1] = i;
        }
    }
    set_alpha(alpha);
}


void Position::fill(U32 cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (U32 i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        value_[adj_id] += stone_value;
        adj_[adj_id].remove(cell_id);
    }
}


void Position::unfill(U32 cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (U32 i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        value_[adj_id] -= stone_value;
        adj_[adj_id].readd(cell_id);
    }
}


void Position::make_move(const Move& move) {
    U32 cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    U32 * stone_loc = stone_loc_[turn_];
    U32& n_stones = n_stones_[turn_];
    n_stones--;
    U32 i = 0;
    while (stones[i] != stone_number)
        i++;
    for (; i < n_stones; i++) {
        stones[i] = stones[i+1];
        stone_loc[stones[i]] = i;
    }
    turn_ ^= 1;
    m_ *= -1;

    U64 cell_mask = MASK(cell_id);
    open_mask_ ^= cell_mask;
    open_.remove(cell_id);
    
    if (controls_ & cell_mask)
        n_controls_[BLUE]--;
    else
        n_controls_[RED]--;
    for (U32 p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]--;
            if (stale_[p] & cell_mask)
                n_stale_[p]--;
        }
    }
    fill(cell_id, move.stone_value_);
    
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = (1LL << (U64)adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        U32 open_id = open_.val_[i];
        U64 open_mask = (1LL << (U64)open_id);
        if (! (dead_[RED] & open_mask) && ! (dead_[BLUE] & open_mask)) {
            if (is_dead(open_id, RED, power[BLUE])) {
                n_dead_[RED]++;
                dead_[RED] |= open_mask;
            } else if (is_dead(open_id, BLUE, power[RED])) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= open_mask;
            }
        }
    }
    find_stale_cells();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    if (n_stale_[RED] > n_dead_[RED]) {
        print(stderr);
        fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[RED], n_dead_[RED]);
        fflush(stderr);
    }
    if (n_stale_[BLUE] > n_dead_[BLUE]) {
        print(stderr);
        fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[BLUE], n_dead_[BLUE]);
        fflush(stderr);
    }
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::unmake_move(const Move& move) {
    turn_ ^= 1;
    m_ *= -1;
    
    U32 cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    U32 * stone_loc = stone_loc_[turn_];
    U32& n_stones = n_stones_[turn_];
    U32 i;
    for (i = n_stones; i && stone_number < stones[i-1]; i--) {
        stones[i] = stones[i-1];
        stone_loc[stones[i]] = i;
    }
    stones[i] = stone_number;
    n_stones++;
   
    U64 cell_mask = (1LL << (U64)cell_id);
    open_mask_ |= cell_mask;
    open_.readd(cell_id); 
    
    unfill(cell_id, move.stone_value_);
    if (controls_ & cell_mask)
        n_controls_[BLUE]++;
    else
        n_controls_[RED]++;
    for (U32 p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]++;
            if (stale_[p] & cell_mask)
                n_stale_[p]++;
        }
    }
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = (1LL << (U64)adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        U32 open_id = open_.val_[i];
        U64 open_mask = (1LL << (U64)open_id);
        if ((dead_[RED] & open_mask)) {
            if (! is_dead(open_id, RED, power[BLUE])) {
                n_dead_[RED]--;
                dead_[RED] ^= open_mask;
            }
        } else if ((dead_[BLUE] & open_mask)) {
            if (! is_dead(open_id, BLUE, power[RED])) {
                n_dead_[BLUE]--;
                dead_[BLUE] ^= open_mask;
            }
        }
    }

    U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
    U64 stale = open_mask_ & (stale_[RED] | stale_[BLUE]);
    U64 stale_it = stale;
    while (stale_it) {
        U64 lsb = LSB(stale_it);
        stale_it ^= lsb;
        if (lsb & stale) {
            bool still_stale = true;
            if (! (dead & lsb))
                still_stale = false;
            else {
                U64 adj_mask = ADJ_MASK[INDEX(lsb)];
                if ((dead & adj_mask) != (open_mask_ & adj_mask))
                    still_stale = false;
            }
            if (! still_stale) {
                if (stale_[RED] & lsb) {
                    stale_[RED] ^= lsb;
                    n_stale_[RED]--;
                } else {
                    stale_[BLUE] ^= lsb;
                    n_stale_[BLUE]--;
                }
            }
        }
    }
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::get_stone_power(U32 p, Value * power) {
    U32 op = 1 - p;
    U32 n_stones = n_stones_[p];
    Value * stones = stones_[p];
    power[0] = 0;
    U32 i;
    for (i = 1; i <= MAX_DEGREE && i <= n_stones; i++)
        power[i] = power[i-1] + stones[n_stones-i] * parity[p];
    U32 max_our_stones = i;
    for (i = 0; i < n_stones_[op] && i + max_our_stones <= MAX_DEGREE; i++)
        power[i + max_our_stones] = power[i + max_our_stones - 1] + stones_[op][i] * parity[op];
}


Move Position::get_random_move() {
    return Move(open_.val_[rand() % open_.len_],
                m_ * stones_[turn_][rand() % n_stones_[turn_]]);
}


U32 lottery(U32 * tickets, U32 n_tickets) {
    int winner = rand() % n_tickets;
    U32 i;
    for (i = 0; winner >= 0; i++)
        winner -= tickets[i];
    return i - 1;
}


Move Position::get_default_policy_move() {
    U32 tickets[N_CELLS];
    U32 total_tickets = 0;
    U32 i;
    U32 best_stale = NO_CELL;
    U32 cell_id;
    Value best_stale_val = 1000;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        cell_id = open_.val_[i];
        U64 mask = (1LL << (U64)cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            Value stale_val = m_ * value_[cell_id];
            if (stale_val < best_stale_val) {
                best_stale = i;
                best_stale_val = stale_val;
            }
        } else {
            if ((dead_[1-turn_] & mask) || (n_stones_[turn_] > n_dead_[1-turn_])) {
                if ((bool)(controls_ & (1LL << U64(cell_id))) == (bool)turn_)
                    tickets[i] = 1;
                else
                    tickets[i] = 3;
                total_tickets += tickets[i];
            }
        }
    }
    if (best_stale != NO_CELL) {
        tickets[best_stale] = 1;
        total_tickets++;
    }
    
    U32 cell_index = lottery(tickets, total_tickets);
    cell_id = open_.val_[cell_index];

    Value reqs[MAX_DEGREE];
    U32 n_uncontrolled = 0;
    Value req_val = parity[turn_] * (alpha_ - OFFSET[turn_]);
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        Value adj_val = m_ * value_[adj.val_[i]];
        if (adj_val < req_val) {
            reqs[n_uncontrolled] = req_val - adj_val;
            n_uncontrolled++;
        }
    }
    std::sort(reqs, reqs + n_uncontrolled);
    
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    U32 stone_index;
    if (cell_index == best_stale)
        stone_index = 0;
    else {
        total_tickets = 0;
        for (i = n_stale_[1 - turn_]; i < n_stones; i++) {
            U32 changes;
            for (changes = 0; changes < n_uncontrolled && stones[i] <= reqs[changes]; changes++) {}
            tickets[i] = 20 * (changes + 1) - stones[i];
            total_tickets += tickets[i];
        }
#ifdef DEBUG_
        if (total_tickets == 0) {
            print(stderr);
            fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[RED], n_dead_[RED]);
            fprintf(stderr, "Blue: %u stale, %u dead\n", n_stale_[BLUE], n_dead_[BLUE]);
            fflush(stderr);
        }
        assert(total_tickets > 0);
#endif
        stone_index = lottery(tickets, total_tickets);
    }
    
    Value stone_value = m_ * stones[stone_index];
    return Move(cell_id, stone_value);
}


Move Position::get_expectation_maximising_move() {
    Move best_move;
    Real best_expectation = -1000.0;
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real  expectation = calculate_expectation();
            unmake_move(move);
            expectation *= m_;
            if (expectation > best_expectation) {
                best_move = move;
                best_expectation = expectation;
            }
        }
    }
    return best_move;
}


U32 ways[MAX_DEGREE + 1][2 * N_STONES + 1][2 * MAX_RESULT + 1];
U32 combos[MAX_DEGREE + 1];
Value stone_list[2 * N_STONES];

Real Position::calculate_win_prob(Value alpha) const {
    U32 i, j, n_stones = 0;
    Value k;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_list[n_stones++] = stones_[i][j] * parity[i];
    //for (i = 0; i < n_stones; i++)
    //    fprintf(stderr, "%d\n", stone_list[i]);
    memset(ways, 0, sizeof(ways));
    for (i = 0; i <= n_stones; i++)
        ways[0][i][MAX_RESULT] = 1;
    for (i = 1; i <= MAX_DEGREE && i <= n_stones; i++) {
        for (j = 1; j <= n_stones; j++) {
            Value val = stone_list[j-1];
            for (k = 0; k <= 2 * MAX_RESULT; k++) {
                U32& ans = ways[i][j][k];
                ans = ways[i][j-1][k];
                if (k - val >= 0 && k - val <= 2 * MAX_RESULT)
                    ans += ways[i-1][j-1][k-val];
            }
        }
    }
    //fprintf(stderr, "n_stones = %u\n", n_stones);
    combos[0] = 1;
    for (i = 1; i <= MAX_DEGREE && i <= n_stones; i++) {
        combos[i] = combos[i-1] * (n_stones + 1 - i) / i;
#ifdef DEBUG_
        U32 n_ways = 0;
        for (j = 0; j <= 2 * MAX_RESULT; j++) {
            //fprintf(stderr, "%u: %u\n", j, ways[i][2][j]);
            n_ways += ways[i][n_stones][j];
        }
        if (n_ways != combos[i]) {
            fprintf(stderr, "n_ways = %u, n_combos = %u\n", n_ways, combos[i]);
            fflush(stderr);
        }
        assert(n_ways == combos[i]);
#endif
    }
    Real prob_sum = 0.0;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        U32 n_adj = adj_[cell_id].len_;
        Value val = value_[cell_id];
        Value diff = alpha - val;
        U32 n_acceptable_outcomes = 0;
        for (k = diff + MAX_RESULT; k <= 2 * MAX_RESULT; k++)
            n_acceptable_outcomes += ways[n_adj][n_stones][k];
        prob_sum += (Real)n_acceptable_outcomes / (Real)combos[n_adj];
    }
    return prob_sum / open_.len_;
}


std::pair<Real, Move> Position::get_best_alpha_move(Value alpha) {
    Move best_move;
    Real best_prob = -1000.0;
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real prob = calculate_win_prob(alpha);
            unmake_move(move);
            if (turn_ == BLUE)
                prob = 1.0 - prob;
            if (prob > best_prob) {
                best_move = move;
                best_prob = prob;
            }
        }
    }
    return std::pair<Real, Move>(best_prob, best_move);
}


std::pair<Real, Move> Position::get_best_move() {
    Value alpha = 0;
    Value best_alpha = 400;
    Move best_move;
    Real best_prob = 0.0;
    while (alpha != best_alpha && alpha >= MIN_RESULT && alpha <= MAX_RESULT) {
        std::pair<Real, Move> result = get_best_alpha_move(alpha);
        fprintf(stderr, "Alpha = %d: %.5lf\n", alpha, result.first);
        if (result.first >= 0.5) {
            best_move = result.second;
            best_prob = result.first;
            best_alpha = alpha;
            alpha += m_;
        } else {
            alpha -= m_;
        }
    }
    return std::pair<Real, Move>(best_prob, best_move);
}


Real Position::calculate_expectation() const {
    if (open_.len_ == 1)
        return (Real)(value_[open_.val_[0]]);
    Real sum = 0.0;
    Real stone_sum = 0.0;
    U32 i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += parity[i] * stones_[i][j];
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        sum += value_[cell_id] + adj_[cell_id].len_ * stone_sum / (open_.len_ - 1);
    }
    return sum / open_.len_;
}


std::pair<Real, Move> Position::search_expectation(U32 depth, Real a, Real b) {
    if (depth == 0)
        return std::pair<Real, Move>(calculate_expectation(), Move());
    Real best_val = -1000.0;
    Move best_move;
    U32 n_stones = n_stones_[turn_];
    Value * stones = stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        for (U32 j = 0; j < n_stones; j++) {
            Move move(cell_id, stones[j] * m_);
            make_move(move);
            std::pair<Real, Move> result = search_expectation(depth - 1, b, a);
            unmake_move(move);
            Real val = result.first * m_;
            if (val > best_val) {
                best_move = move;
                best_val = val;
                if (val > a) {
                    a = val;
                    if (a >= -b)
                        return std::pair<Real, Move>(m_ * best_val, best_move);
                }
            }
        }
    }
    return std::pair<Real, Move>(m_ * best_val, best_move);
}

        
Real Position::get_control_heuristic() const {
    return (Real)n_controls_[RED] / (Real)open_.len_;
}


void Position::get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real expectation = calculate_expectation();
            unmake_move(move);
            moves.push_back(std::make_pair(m_ * expectation, move));
        }
    }
}


void Position::get_all_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            moves.push_back(move);
        }
    }
}


void Position::get_all_reasonable_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    Value stone_value;
    Value best_stale_val = 1000;
    U32 best_stale = NO_CELL;
    U32 i, j;
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        Value cell_value = value_[cell_id] * m_;
        U64 mask = MASK(cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            if (cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            if ((dead_[1-turn_] & mask) || (n_stones > n_dead_[1-turn_])) {
                bool done = false;
                if (adj_[cell_id].len_ == 1) {
                    U32 adj_id = adj_[cell_id].val_[0];
                    if (adj_[adj_id].len_ == 1) {
                        done = true;
                        Value adj_value = value_[adj_id] * m_;
                        if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                            stone_value = m_ * stones[0];
                            moves.push_back(Move(cell_id, stone_value));
                            Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                            if (adj_value + stones[0] < control_req) {
                                Value extra = control_req - adj_value;
                                for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                                if (j < n_stones) {
                                    stone_value = m_ * stones[j];
                                    moves.push_back(Move(cell_id, stone_value));
                                }
                            }       
                        }
                    } else {
                        if (we_control & mask)
                            done = true;
                    }
                }
                if (! done) {
                    for (j = n_stale_[1 - turn_]; j < n_stones; j++) {
                        stone_value = m_ * stones[j];
                        moves.push_back(Move(cell_id, stone_value));
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * m_;
        moves.push_back(Move(best_stale, stone_value));
    }
}


struct MoveInfo {
    Move move_;
    int swing_;

    MoveInfo(const Move& move, int swing): move_(move), swing_(swing) {}

    bool operator < (const MoveInfo& other) const {
        //if (stale_ != other.stale_)
        //    return stale_;
        //if (dead_ != other.dead_)
        //    return ! dead_;
        if (move_.stone_value_ != other.move_.stone_value_)
            return abs(move_.stone_value_) > abs(other.move_.stone_value_);
        if (swing_ != other.swing_)
            return swing_ > other.swing_;
        return move_.cell_ < other.move_.cell_;
    }
};


U32 Position::solve(long long end_time, U32& counter, U32& hash_hits,
        HashTable& table, HashInfo& hash_info, U64 stone_masks[2]) {
    counter++;
    if (counter % 1000 == 0 && get_time() > end_time)
        return TIME_ELAPSED;
    if (is_winning(RED))
        return RED;
    if (is_winning(BLUE))
        return BLUE;
    U32 hash_result = table.find(stone_masks[0], stone_masks[1]);
    if (hash_result != NO_HASH_ENTRY) {
        hash_hits++;
        return hash_result;
    }
    U64 child_stone_masks[2];
    child_stone_masks[0] = stone_masks[0];
    child_stone_masks[1] = stone_masks[1];
    std::vector<Move> moves;
    get_all_reasonable_moves(moves);
    std::vector<MoveInfo> move_infos;
    U32 i, j;
    U32 op = 1 - turn_;
    U64 opp_controls = controls_;
    if (turn_ == BLUE)
        opp_controls = ~opp_controls;
    for (i = 0; i < moves.size(); i++) {
        Move move = moves[i];
        U32 cell_id = move.cell_;
        U64 mask = MASK(cell_id);
        int swing = 0;
        if (dead_[turn_] & mask)
            swing = -1;
        else if (dead_[op] & mask)
            swing = 1;
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            U32 adj_id = adj.val_[j];
            U64 adj_mask = MASK(adj_id);
            if ((opp_controls & adj_mask) && ! (dead_[op] & adj_mask))
                swing++;
        }
        move_infos.push_back(MoveInfo(move, swing));
    }
    std::sort(move_infos.begin(), move_infos.end());
    U32 optimal_result = op;
    for (i = 0; i < move_infos.size(); i++) {
        Move move = move_infos[i].move_;
        Value stone_number = move.stone_value_ * m_;
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[move.cell_];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move(move);
        if (result == TIME_ELAPSED)
            return TIME_ELAPSED;
        if (result == turn_) {
            optimal_result = turn_;
            break;
        }
    }
    table.add(stone_masks[0], stone_masks[1], optimal_result);
    return optimal_result;
}


std::pair<U32, Move> Position::get_optimal_move(long long end_time, U32& counter, U32& hash_hits,
        bool break_ties, HashTable& table) {
    table.init();
    HashInfo hash_info(*this);
    U64 stone_masks[2] = {0};
    U64 child_stone_masks[2];
    child_stone_masks[0] = stone_masks[0];
    child_stone_masks[1] = stone_masks[1];
    std::vector<Move> moves;
    get_all_reasonable_moves(moves);
    
    std::vector<MoveInfo> move_infos;
    U32 i, j;
    U32 op = 1 - turn_;
    U64 opp_controls = controls_;
    if (turn_ == BLUE)
        opp_controls = ~opp_controls;
    Move move;
    for (i = 0; i < moves.size(); i++) {
        move = moves[i];
        U32 cell_id = move.cell_;
        U64 mask = MASK(cell_id);
        int swing = 0;
        if (dead_[turn_] & mask)
            swing = -1;
        else if (dead_[op] & mask)
            swing = 1;
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            U32 adj_id = adj.val_[j];
            U64 adj_mask = MASK(adj_id);
            if ((opp_controls & adj_mask) && ! (dead_[op] & adj_mask))
                swing++;
        }
        move_infos.push_back(MoveInfo(move, swing));
    }
    std::sort(move_infos.begin(), move_infos.end());
   
#ifdef DEBUG_
    fprintf(stderr, "Solving: alpha = %d\n", alpha_);
#endif
    std::vector<Move> winning_moves;
    for (i = 0; i < move_infos.size(); i++) {
        move = move_infos[i].move_;
#ifdef DEBUG_
        char move_str[10];
        move.to_str(move_str);
        fprintf(stderr, "%s\n", move_str);
#endif
        Value stone_number = move.stone_value_ * m_;
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[move.cell_];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move(move);
        if (result == TIME_ELAPSED)
            return std::pair<U32, Move>(TIME_ELAPSED, Move());
        if (result == turn_) {
            if (break_ties)
                winning_moves.push_back(move);
            else
                return std::pair<U32, Move>(turn_, move);
        }
    }
    
    if (winning_moves.size()) {
        fprintf(stderr, "%u optimal moves\n", winning_moves.size());
        Move best_move;
        Real best_val = -1000.0;
        for (i = 0; i < winning_moves.size(); i++) {
            move = winning_moves[i];
            make_move(move);
            Real val = calculate_expectation();
            unmake_move(move);
            val *= m_;
            if (val > best_val) {
                best_move = move;
                best_val = val;
            }
        }
        return std::pair<U32, Move>(turn_, best_move);
    }
    return std::pair<U32, Move>(1 - turn_, moves[0]);
}


bool Position::is_dead(U32 cell_id, U32 p, Value * other_power) {
    Value value = value_[cell_id];
    U32 n_adj = adj_[cell_id].len_;
    Value worst_val = parity[p] * (value + other_power[n_adj]);
    return (worst_val >= parity[p] * (alpha_ - OFFSET[p]));
}


bool Position::all_adj_dead(U32 cell_id) {
    List& adj = adj_[cell_id];
    U32 n_adj = adj.len_;
    for (U32 i = 0; i < n_adj; i++) {
        U64 adj_mask = (1LL << (U64)adj.val_[i]);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
            return false;
    }
    return true;
}


void Position::find_stale_cells() {
    U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
    U64 dead_it = dead;
    U64 stale = open_mask_ & (stale_[RED] | stale_[BLUE]);
    while (dead_it) {
        U64 lsb = LSB(dead_it);
        dead_it ^= lsb;
        if (! (stale & lsb)) {
            U64 adj_mask = ADJ_MASK[INDEX(lsb)];
            if ((dead & adj_mask) == (open_mask_ & adj_mask)) {
                if (dead_[RED] & lsb) {
                    stale_[RED] |= lsb;
                    n_stale_[RED]++;
                } else {
                    stale_[BLUE] |= lsb;
                    n_stale_[BLUE]++;
                }
            }
        }
    }
}


bool Position::is_winning(U32 p) const {
    return n_dead_[p] > n_stones_[1 - p];
}


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    Value power[2][MAX_DEGREE+1];
    dead_[0] = dead_[1] = 0;
    stale_[RED] = stale_[BLUE] = 0;
    n_stale_[RED] = n_stale_[BLUE] = 0;
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    U32 i, cell_id, n_adj;
    U64 mask;
    for (i = 0; i < open_.len_; i++) {
        cell_id = open_.val_[i];
        Value value = value_[cell_id];
        n_adj = adj_[cell_id].len_;
        mask = (1LL << (U64)cell_id);
        if (value >= alpha_) {
            n_controls_[RED]++;
            if (value + power[BLUE][n_adj] >= alpha_) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            }
        } else {
            controls_ |= mask;
            n_controls_[BLUE]++;
            if (value + power[RED][n_adj] < alpha_) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= mask;
            }
        }
    }
    find_stale_cells();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::take_snapshot() {
    U32 i, j, p;
    snapshot_.turn_ = turn_;
    snapshot_.open_mask_ = open_mask_;
    snapshot_.n_open_ = open_.len_;
    snapshot_.controls_ = controls_;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        snapshot_.open_[i] = cell_id;
        snapshot_.value_[cell_id] = value_[cell_id];
        List& adj = adj_[cell_id];
        U32 n_adj = adj.len_;
        snapshot_.n_adj_[cell_id] = n_adj;
        for (j = 0; j < n_adj; j++)
            snapshot_.adj_[cell_id][j] = adj.val_[j];
    }

    for (p = 0; p < 2; p++) {
        U32 n_stones = n_stones_[p];
        snapshot_.n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        Value * snapshot_stones = snapshot_.stones_[p];
        for (i = 0; i < n_stones; i++)
            snapshot_stones[i] = stones[i];

        snapshot_.n_controls_[p] = n_controls_[p];
        snapshot_.dead_[p] = dead_[p];
        snapshot_.n_dead_[p] = n_dead_[p];
        snapshot_.stale_[p] = stale_[p];
        snapshot_.n_stale_[p] = n_stale_[p];
    }

}


void Position::restore_snapshot() {
    U32 i, j, p;
    turn_ = snapshot_.turn_;
    m_ = parity[turn_];
    open_.len_ = snapshot_.n_open_;
    open_mask_ = snapshot_.open_mask_;
    controls_ = snapshot_.controls_;
    for (i = 0; i < snapshot_.n_open_; i++) {
        U32 cell_id = snapshot_.open_[i];
        value_[cell_id] = snapshot_.value_[cell_id];
        open_.val_[i] = cell_id;
        open_.loc_[cell_id] = i;
        U32 n_adj = snapshot_.n_adj_[cell_id];
        List& adj = adj_[cell_id];
        adj.len_ = n_adj; 
        for (j = 0; j < n_adj; j++) {
            U32 adj_id = snapshot_.adj_[cell_id][j];
            adj.val_[j] = adj_id;
            adj.loc_[adj_id] = j;
        }
    }
    
    for (p = 0; p < 2; p++) {
        U32 n_stones = snapshot_.n_stones_[p];
        n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        U32 * stone_loc = stone_loc_[p];
        Value * snapshot_stones = snapshot_.stones_[p];
        for (i = 0; i < n_stones; i++) {
            stones[i] = snapshot_stones[i];
            stone_loc[stones[i]] = i;
        }

        n_controls_[p] = snapshot_.n_controls_[p];
        dead_[p] = snapshot_.dead_[p];
        n_dead_[p] = snapshot_.n_dead_[p];
        stale_[p] = snapshot_.stale_[p];
        n_stale_[p] = snapshot_.n_stale_[p];
    }
}


std::pair<U32, U32> Position::get_cell_and_stone_indices(const Move& move) {
    return std::pair<U32, U32>(open_.loc_[move.cell_],
                                     stone_loc_[turn_][m_ * move.stone_value_]);
}


void Position::print(FILE * f) {
    U32 i, j;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        char cell_str[10];
        cell_id_to_name(cell_id, cell_str);
        fprintf(f, "%s (%d):", cell_str, value_[cell_id]);
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            cell_id_to_name(adj.val_[j], cell_str);
            fprintf(f, " %s", cell_str);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "Red:");
    for (i = 0; i < n_stones_[RED]; i++)
        fprintf(f, " %d", stones_[RED][i] * parity[RED]);
    fprintf(stderr, "\nBlue:");
    for (i = 0; i < n_stones_[BLUE]; i++)
        fprintf(stderr, " %d", stones_[BLUE][i] * parity[BLUE]);
    fprintf(stderr, "\n");
}

class MCTNode {
public:
    Move move_;
    MCTNode * parent_;

    Value alpha_;
    Real val_;
    U32 n_red_wins_;
    U32 n_playouts_;

    bool fully_explored_;
    bool expanded_;

    U32 n_children_fully_explored_;
    U32 n_children_;
    U32 n_child_moves_;
    std::vector<Move> child_moves_;
    std::vector<MCTNode*> children_;
    
    bool solved_;
    bool solve_attempted_;
    Move solution_;
    U32 solver_positions_;
    U32 solver_hash_hits_;
    long long solver_time_;

    AMAFTable amaf_[2];


    void init(MCTNode * parent, const Position& pos, const Move& move, Value alpha);

    void ucb(Position& pos);
    
    MCTNode * select(Position& pos, AMAFTable& amaf);

    MCTNode *  add_child(Position& pos, const Move& move);

    void get_children(Position& pos);

    MCTNode * expand(Position& pos, AMAFTable& amaf);

    bool light_playout(Position& pos, U32& move_count);

    bool attempt_solve(Position& pos, HashTable& table, long long allowed_time, bool break_ties);

    Move get_most_played_move();

    Move get_highest_value_move(Position& pos);

    void dispose();

    void print(FILE * f);

};

MCTNode * get_free();

void put_free(MCTNode * node);

void init_free_list();


class MCTSearch {
public:

    MCTNode * roots_[MAX_RESULT * 2 + 1];
    Value current_alpha_;
    HashTable table_;

    MCTSearch(const Position& pos, Value alpha);

    ~MCTSearch();

    MCTNode * select_alpha(const Position& pos);

    MCTNode * get_or_make_root(Value alpha, const Position& pos);

    Move get_best_move(Position& pos);

};


void MCTNode::dispose() {
    for (U32 i = 0; i < children_.size(); i++) {
        children_[i]->dispose();
        put_free(children_[i]);
    }
}


void MCTNode::init(MCTNode * parent, const Position& pos, const Move& move, Value alpha) {
    move_= move;
    parent_ = parent;
    alpha_ = alpha;
    if (pos.is_winning(RED)) {
        fully_explored_ = true;
        val_ = 1.0;
    } else if (pos.is_winning(BLUE)) {
        fully_explored_ = true;
        val_ = 0.0;
    } else {
        fully_explored_ = false;
        val_ = 0.0; 
    }
    n_red_wins_ = 0;
    n_playouts_ = 0;
    child_moves_.clear();
    children_.clear();
    expanded_ = false;
    n_children_ = 0;
    n_child_moves_ = 0;
    n_children_fully_explored_ = 0;
    solved_ = false;
    solve_attempted_ = false;
    solver_positions_ = 0;
    solver_hash_hits_ = 0;
    solver_time_ = 0;
    if (parent == NULL) {
        amaf_[0].init(pos.open_.len_, pos.n_stones_[RED]);
        amaf_[1].init(pos.open_.len_, pos.n_stones_[BLUE]);
    }
}


MCTNode * MCTNode::select(Position& pos, AMAFTable& amaf) {
    MCTNode * best_node = NULL;
    Real best_node_value = -1000.0, child_val;
    U32 i;
    for (i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            child_val = child->val_;
            if (pos.turn_ == BLUE)
                child_val = 1.0 - child_val;
            Real parent_playouts = (Real)n_playouts_;
            Real child_playouts = (Real)child->n_playouts_;
            Real ucb_value = child_val + UCB_C * sqrt(log(parent_playouts) / child_playouts);
            if (ucb_value > best_node_value) {
                best_node = child;
                best_node_value = ucb_value;
            }
        }
    }
    
    int best_move_index = -1;
    Real best_move_value = -1000.0;
    for (i = 0; i < n_child_moves_; i++) {
        Move& move = child_moves_[i];
        std::pair<U32, U32> indices = pos.get_cell_and_stone_indices(move);
        U32 n_amaf_playouts = amaf.get_n_playouts(indices.first, indices.second); 
        if (n_amaf_playouts > 0)
            child_val = amaf.get_value(indices.first, indices.second);
        else
            child_val = 0.5;
        Real parent_playouts = std::max(1.0, (Real)n_playouts_);
        child_val += UCB_C * sqrt(log(parent_playouts));
        if (child_val > best_move_value) {
            best_move_index = i;
            best_move_value = child_val;
        }
    }
    if (best_node == NULL || best_move_value > best_node_value) {
#ifdef DEBUG_
        assert(best_move_index != -1);
#endif
        n_child_moves_--;
        Move best_move = child_moves_[best_move_index];
        child_moves_[best_move_index] = child_moves_[n_child_moves_];
        return add_child(pos, best_move);
    }
#ifdef DEBUG_
    assert(best_node != NULL);
#endif
    return best_node;
}


MCTNode *  MCTNode::add_child(Position& pos, const Move& move) {
    pos.make_move(move);
    MCTNode * child = get_free();
    child->init(this, pos, move, alpha_);
    children_.push_back(child);
    pos.unmake_move(move);
    return child;
}


void MCTNode::get_children(Position& pos) {
    pos.get_all_reasonable_moves(child_moves_);
    n_child_moves_ = n_children_ = child_moves_.size();
}


MCTNode * MCTNode::expand(Position& pos, AMAFTable& amaf) {
    expanded_ = true;
    get_children(pos);
    return select(pos, amaf);
}


std::pair<U32, U32> simulation[N_CELLS];

bool MCTNode::light_playout(Position& pos, U32& move_count) {
    bool result = true;
    bool result_found = false;
    pos.take_snapshot();
    while (pos.open_.len_ > 1) {
        if (pos.is_winning(RED)) {
            result = true;
            result_found = true;
            break;
        }
        if (pos.is_winning(BLUE)) {
            result = false;
            result_found = true;
            break;
        } 
        Move move = pos.get_default_policy_move();
        simulation[move_count] = pos.get_cell_and_stone_indices(move);
        pos.make_move(move);
        move_count++;
    }
    if (! result_found)
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    pos.restore_snapshot();
    return result;
}


void MCTNode::ucb(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    MCTNode * search_root = this;
    U32 orig_turn = pos.turn_;
    U32 move_count = 0, i;
    while (search_root->expanded_) {
        search_root = search_root->select(pos, amaf_[pos.turn_]);
#ifdef DEBUG_
        assert(search_root != NULL);
#endif
        simulation[move_count] = pos.get_cell_and_stone_indices(search_root->move_);
        pos.make_move(search_root->move_);
        move_count++;
    }
    
    if (! search_root->fully_explored_) { 
        search_root = search_root->expand(pos, amaf_[pos.turn_]);
        simulation[move_count] = pos.get_cell_and_stone_indices(search_root->move_);
        move_count++;
        pos.make_move(search_root->move_);
    }
    U32 result;
    if (search_root->fully_explored_)
        result = (search_root->val_ == 1.0);
    else
        result = search_root->light_playout(pos, move_count);

    for (i = 0; i < move_count; i++) {
        U32 turn = (orig_turn + i) % 2;
        amaf_[turn].update(simulation[i].first, simulation[i].second, result != turn);
    }

    while (search_root != NULL) {
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (i = 0; i < search_root->children_.size(); i++)
                if (search_root->children_[i]->val_ * pos.m_ > search_root->val_)
                    search_root->val_ = search_root->children_[i]->val_ * pos.m_;
            search_root->val_ *= pos.m_;
        }
        if (! search_root->fully_explored_) {
            search_root->n_red_wins_ += result;
            search_root->n_playouts_++;
            search_root->val_ = Real(search_root->n_red_wins_) / Real(search_root->n_playouts_);
        }
        if (parent != NULL) {
            pos.unmake_move(search_root->move_);
            if (search_root->fully_explored_) {
                parent->n_children_fully_explored_++;
                if ((parent->n_children_fully_explored_ == parent->n_children_) ||
                        (pos.turn_ == RED && search_root->val_ == 1.0) ||
                        (pos.turn_ == BLUE && search_root->val_ == 0.0)) {
                    parent->fully_explored_ = true;
                }
            }
        }
        search_root = parent;
    }
}


Move MCTNode::get_highest_value_move(Position& pos) {
    if (solved_)
        return solution_;
    Move best_move;
    Real best_val = -1000.0;
    U32 n_playouts = 0;
    for (U32 i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if ((! fully_explored_ || child->fully_explored_) &&
                child->val_ * pos.m_ > best_val) {
            best_move = child->move_;
            best_val = child->val_ * pos.m_;
            n_playouts = child->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %u playouts\n", n_playouts);
    return best_move;
}


Move MCTNode::get_most_played_move() {
    Move best_move;
    U32 most_visits = 0;
    for (U32 i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (child->n_playouts_ > most_visits) {
            best_move = child->move_;
            most_visits = child->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %d playouts\n", most_visits);
    return best_move;
}


bool MCTNode::attempt_solve(Position& pos, HashTable& table, long long allowed_time, bool break_ties) {
    solve_attempted_ = true;
    long long start_time = get_time();
    long long end_time = start_time + allowed_time;
    std::pair<U32, Move> result = pos.get_optimal_move(end_time, solver_positions_, solver_hash_hits_, break_ties, table);
    solver_time_ = get_time() - start_time;
    if (result.first != TIME_ELAPSED) {
        if (result.first == RED)
            val_ = 1.0;
        else
            val_ = 0.0;
        solved_ = true;
        fully_explored_ = true;
        solution_ = result.second;
        return true;
    }
    return false;
}


void MCTNode::print(FILE * f) {
    char playout_status[20];
    char solve_status[100];
    if (solved_)
        sprintf(playout_status, "*");
    else
        sprintf(playout_status, "%u", n_playouts_);
    if (solve_attempted_) {
        Real solver_time = (Real)solver_time_ / 1.0e6;
        if (solved_)
            sprintf(solve_status, "(solver SUCCEEDED: %u nodes, %u hash hits, %.4lfs)",
                    solver_positions_, solver_hash_hits_, solver_time);
        else
            sprintf(solve_status, "(solver FAILED: %u nodes, %u hash hits, %.4lfs)",
                    solver_positions_, solver_hash_hits_, solver_time);
    } else {
        solve_status[0] = 0;
    }
    fprintf(f, "Alpha = %d (%s): %.4f %s\n", alpha_, playout_status, val_, solve_status);
}



MCTNode node_store[N_MCT_NODES];
MCTNode * free_list[N_MCT_NODES];
U32 n_free;


MCTNode * get_free() {
#ifdef DEBUG_
    assert(n_free > 0);
#endif
    return free_list[--n_free];
}


void put_free(MCTNode * node) {
    free_list[n_free++] = node;
}


void init_free_list() {
    for (U32 i = 0; i < N_MCT_NODES; i++)
        free_list[i] = &node_store[i];
    n_free = N_MCT_NODES;
}


MCTSearch::MCTSearch(const Position& pos, Value alpha) {
    current_alpha_ = alpha;
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++)
        roots_[i] = NULL;
}


MCTSearch::~MCTSearch() {
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++) {
        if (roots_[i] != NULL) {
            roots_[i]->dispose();
            put_free(roots_[i]);
        }
    }
}


MCTNode * MCTSearch::get_or_make_root(Value alpha, const Position& pos) {
    Value index = alpha + MAX_RESULT;
#ifdef DEBUG_
    assert(index >= 0 && index <= 2 * MAX_RESULT);
#endif
    if (roots_[index] == NULL) {
        roots_[index] = get_free();
        roots_[index]->init(NULL, pos, Move(), alpha);
        roots_[index]->fully_explored_ = false;
    }
    return roots_[index];
}

MCTNode * MCTSearch::select_alpha(const Position& pos) {
    MCTNode * test_root;
    MCTNode * root = NULL;
    Value alpha;
    if (pos.turn_ == RED) {
        for (alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL) {
                test_root->print(stderr);
                if (test_root->val_ >= 0.5 || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    } else {
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL) {
                test_root->print(stderr);
                if (test_root->val_ < 0.5 || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    }
    return root;
}


Move MCTSearch::get_best_move(Position& pos) {
    if (pos.open_.len_ >= 26) {
        std::pair<Real, Move> search_result = pos.get_best_move();
        //fprintf(stderr, "Eval = %.5lf\n", search_result.first);
        return search_result.second;
    }
    
    Real time_left_r = (Real) time_left;
    Real use_proportion;
    int moves_left = ((int)pos.open_.len_ - 16) / 2;
    if (moves_left < 1)
        use_proportion = 0.5;
    else {
        fprintf(stderr, "Expecting to make %d more moves before position solved\n", moves_left);
        use_proportion = 1.0 - pow(5e5 / time_left_r, 1.0 / (Real)moves_left);
    }
    long long move_time = (long long)(use_proportion * time_left_r);
    long long start_micros = get_time();
    long long time_now = start_micros;
    long long alpha_time = move_time / 3LL;
    long long solver_time = alpha_time / 2LL;
    fprintf(stderr, "Move time = %lld\n", move_time);
    fprintf(stderr, "Alpha time = %lld\n", alpha_time);
    U32 n_playouts = 0;
    U32 i;
    MCTNode * root;
    U32 solver_start = 19;
    
    while (time_now - start_micros < alpha_time) {
        pos.set_alpha(current_alpha_);
        root = get_or_make_root(current_alpha_, pos);
        if (pos.open_.len_ <= solver_start && ! root->solve_attempted_ && time_now - start_micros < solver_time)
            root->attempt_solve(pos, table_, solver_time - time_now + start_micros, false);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->ucb(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
                if (get_or_make_root(current_alpha_, pos)->fully_explored_)
                    break;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                if (get_or_make_root(current_alpha_, pos)->fully_explored_)
                    break;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
        } else {
            if (root->val_ >= 0.5) {
                if (current_alpha_ < MAX_RESULT)
                    current_alpha_++;
            } else {
                if (current_alpha_ > MIN_RESULT)
                    current_alpha_--;
            }
        }
        time_now = get_time();
    }

    root = select_alpha(pos);
    
    Value original_alpha = current_alpha_;
    fprintf(stderr, "Alpha = %d\n", current_alpha_);
    pos.set_alpha(current_alpha_);
    bool done = false;
    while (time_now - start_micros < move_time) {
        //if (pos.open_.len_ <= solver_start && ! root->solve_attempted_)
        //    root->attempt_solve(pos, table_);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->ucb(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
                if (root->fully_explored_) {
                    done = true;
                    break;
                }
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
                if (root->fully_explored_) {
                    done = true;
                    break;
                }
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
        }
        time_now = get_time();
    }

    current_alpha_ = original_alpha;
    if (done)
        root = select_alpha(pos);
    pos.set_alpha(current_alpha_);

    if (root->solved_) {
        long long time_remaining = move_time - (time_now - start_micros);
        if (time_remaining > 0) {
            if (root->attempt_solve(pos, table_, time_remaining, true))
                fprintf(stderr, "Chose most challenging optimal move\n");
        }
    }

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(time_now - start_micros) / 1e6);
        
    if (root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return root->get_highest_value_move(pos);
    }

    return root->get_most_played_move();
}


int main(int argc, char** argv) {
    time_started = get_time();
	init();
	init_free_list();
	
    char cell_str[5];
    U32 blocked[N_BLOCKED_CELLS];
    for (U32 i = 0; i < N_BLOCKED_CELLS; i++) {
        assert(scanf("%s", cell_str));
        blocked[i] = cell_name_to_id(cell_str);
    }
    
    Position pos(blocked, 0);
    char move_str[100];
    get_move(move_str);
    if (strcmp(move_str, "Start"))
        pos.make_move(Move(move_str, pos.turn_));
	Value alpha = 0;
    while (pos.open_.len_ > 1) {
        fprintf(stderr, "Time left: %.2f seconds\n", (float)time_left / 1e6);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES);
#endif
        MCTSearch search(pos, alpha);
        Move move = search.get_best_move(pos);
        fprintf(stderr, "Used %u nodes\n", N_MCT_NODES - n_free);
        alpha = search.current_alpha_;
		pos.make_move(move);
        move.to_str(move_str);
		send_move(move_str);
		get_move(move_str);
        pos.make_move(Move(move_str, pos.turn_));
	}
	return 0;
}




// EventHorizon0_0_0: plays move which maximises expected pure Monto-Carlo value
// EventHorizon1_0_0: untuned UCT implementation, scores -8 vs 0_0_0 (100 games)
// EventHorizon1_0_1: 0.3s/move, scores -3.66 vs 0_0_0 (100 games)
// EventHorizon1_0_2: remove move ordering, scores -2.82 vs 0_0_0 (100 games)
// EventHorizon1_0_3: use bitmask move table, scores -2.55  vs 0_0_0 (100 games)
// EventHorizon1_1_0: fixed target of 0, wins a lot of games by 1 point, loses horribly when it loses, scores -8.86  vs 0_0_0 (100 games)
// EventHorizon1_1_2: moving target, scores -6.02  vs 0_0_0 (100 games)
// EventHorizon1_1_3: don't play senseless moves in isolated cells, -4.74  vs 0_0_0 (100 games)
// EventHorizon1_1_4: maximise expectation for first 5 moves, create default policy, -0.39  vs 0_0_0 (100 games)
// EventHorizon1_1_5: disallow senseless moves with doubles, +0.2  vs 0_0_0 (100 games)
// EventHorizon1_1_6: default policy with lottery, fix propagation for fully explored positions, +0.86 vs 0_0_0, +0.46 vs player3
// EventHorizon1_1_7: add dead tiles to Position +0.62 vs player3
// EventHorizon1_2_0: remove Cell class, add stale cell tracking
// EventHorizon1_2_1: don't generate moves with stones which can be placed in stale cells, +0.95 vs player3
// EventHorizon2_0_0: use AMAF for move ordering, +0.37 vs 1_2_2 (400 games)
// EventHorizon2_1_0: use bit operations to find stale cells
// EventHorizon2_2_0: use win prob with random play at alpha level for first 3 moves, +0.55 vs 2_1_0 (400 gmes)
// EventHorizon2_3_0: implement endgame solver
// EventHorizon2_4_3: add hash table to endgame solver, adjust time management
// EventHorizon2_4_6: adjust UCB coefficient



