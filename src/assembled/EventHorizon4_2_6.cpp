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
typedef unsigned char U8;
typedef unsigned short int U16;
typedef unsigned int U32;
typedef unsigned long long U64;
typedef double Real;

#define LSB(x) ((x) & -(x))
#define INDEX(x) (index64[((x) * debruijn) >> 58])
#define N_MCT_NODES 100000
#define RAVE_C 1e-4
#define UCB_C 0.003125
#define AMAF_HORIZON 3
#define UNIFORM_TM

const char * ENGINE_NAME = "EventHorizon";
const char * VERSION_NUMBER = "4.2.6";

int parity[2] = {1, -1};
Value OFFSET[2] = {0, 1};

#ifdef FAST_
long long time_limit = 2000000;
long long time_left = 1900000;
#else
#ifdef SLOW_
long long time_limit = 10000000;
long long time_left =   9800000;
#else
long long time_limit = 5000000;
long long time_left =  4900000;
#endif
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
	if (! strcmp(move_str, "Quit"))
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
#ifndef SOLVER_IMPL_
    fprintf(stderr, "R %s %s\n", ENGINE_NAME, VERSION_NUMBER);
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

    for (U32 p = 0; p < 2; p++) {
        Value m = parity[p];
        for (U32 mask = 0; mask < (1 << N_STONES); mask++) {
            Value * power = STONE_POWER[p][mask];
            power[0] = 0;
            U32 i = 1;
            for (int j = N_STONES - 1; j >= 0 && i <= MAX_DEGREE; j--) {
                if (mask & (1 << j)) {
                    power[i] = power[i - 1] + m * (j + 1);
                    i++;
                }
            }
        }
    }

}


const Real TARGET_INCREMENT_THRESH[2] = {0.45, 0.65};
const Real TARGET_DECREMENT_THRESH[2] = {0.35, 0.55};
const Real CHOOSE_TARGET_THRESH[2] = {0.4, 0.6};


#define UNTRIED 0
#define TRIED 1


struct AMAFRecord {
    U16 n_wins_;
    U16 n_playouts_;
    U16 index_;
    U8 cell_index_;
    U8 stone_index_;
    bool played_;
    
    inline void init_generated() {
        n_playouts_ = 0;
        played_ = false;
    }
     
    void init(U32 cell_index, U32 stone_index, bool win) {
        cell_index_ = (U8)cell_index;
        stone_index_ = (U8)stone_index;
        n_wins_= (U16)win;
        n_playouts_ = 1;
    }

    inline bool operator < (const AMAFRecord& other) const {
        return (U64)n_wins_ * (U64)other.n_playouts_ < (U64)other.n_wins_ * (U64)n_playouts_;
    }

};

#define N_AMAF_RECORDS 5000000

void init_amaf_free_list();

AMAFRecord * amaf_alloc(U32 n);

AMAFRecord * get_amaf_record();

class AMAFTable {
public:
    U32 n_stones_;
    U32 n_open_;
    std::vector<AMAFRecord*> heap_;
    std::vector<AMAFRecord*> amaf_;
    
    void init(U32 n_open, U32 n_stones);
   
    AMAFRecord * cell_init(U32 cell_index);
    
    inline void check_initialised() {
        if (amaf_.size() == 0)
            amaf_ = std::vector<AMAFRecord*>(n_open_, NULL);
    }
     
    void play(U32 cell_index, U32 stone_index);
    
    U32 get_n_playouts(U32 cell_index, U32 stone_index);
    
    Real get_value(U32 cell_index, U32 stone_index);

    void update(U32 cell_index, U32 stone_index, bool win); 

    void get_best(U32& cell_index, U32& stone_index);

    bool is_played(U32 cell_index, U32 stone_index);

private:
    void heap_insert(AMAFRecord * record);

    void pop_heap();

    void heap_update(AMAFRecord * record, bool win);
};


AMAFRecord amaf_store[N_AMAF_RECORDS];
U32 amaf_pointer;


void init_amaf_free_list() {
    amaf_pointer = 0;
}


AMAFRecord * amaf_alloc(U32 n) {
#ifdef DEBUG_
    assert(amaf_pointer + n <= N_AMAF_RECORDS);
#endif
    amaf_pointer += n;
    return &amaf_store[amaf_pointer - n];
}


void AMAFTable::init(U32 n_open, U32 n_stones) {
    n_stones_ = n_stones;
    n_open_ = n_open;
    amaf_.clear(); 
    heap_ = std::vector<AMAFRecord*>();
}


AMAFRecord * AMAFTable::cell_init(U32 cell_index) {
#ifdef DEBUG_
    assert(amaf_.size() > 0);
    assert(amaf_[cell_index] == NULL);
#endif
    AMAFRecord * records = amaf_alloc(n_stones_);
    amaf_[cell_index] = records;
    for (U32 i = 0; i < n_stones_; i++)
        records[i].init_generated();
    return records;
}


void AMAFTable::play(U32 cell_index, U32 stone_index) {
    check_initialised();
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        records = cell_init(cell_index);
    AMAFRecord& record = records[stone_index];
#ifdef DEBUG_
    assert(! record.played_);
#endif
    if (record.n_playouts_ != 0)
        pop_heap();
    record.played_ = true;
}


U32 AMAFTable::get_n_playouts(U32 cell_index, U32 stone_index) {
    if (amaf_.empty())
        return 0;
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        return 0;    
    return records[stone_index].n_playouts_;
}

Real AMAFTable::get_value(U32 cell_index, U32 stone_index) {
#ifdef DEBUG_
    assert(amaf_.size() > 0 && amaf_[cell_index] != NULL);
#endif
    AMAFRecord& record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
    assert(record.n_playouts_ != 0);
#endif
    return (Real)record.n_wins_ / (Real)record.n_playouts_;
}

void AMAFTable::update(U32 cell_index, U32 stone_index, bool win) {
    check_initialised();
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        records = cell_init(cell_index);
    AMAFRecord& record = records[stone_index];
    if (record.n_playouts_ == 0) {
        record.init(cell_index, stone_index, win);
        if (! record.played_)
            heap_insert(&record);
    } else {
        record.n_wins_ += win;
        record.n_playouts_++;
#ifdef DEBUG_
        assert(record.n_playouts_ < 60000);
#endif
        if (! record.played_)
            heap_update(&record, win);
    }
}


void AMAFTable::get_best(U32& cell_index, U32& stone_index) {
    if (heap_.empty()) {
        cell_index = NO_CELL;
        return;
    }
    AMAFRecord * record = heap_[0];
    cell_index = record->cell_index_;
    stone_index = record->stone_index_;
    //record->tried_ = true;
    //pop_heap();
}


bool AMAFTable::is_played(U32 cell_index, U32 stone_index) {
    if (amaf_.size() == 0)
        return false;
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        return false; 
    return records[stone_index].played_;
}


void AMAFTable::heap_insert(AMAFRecord * record) {
    U32 index = heap_.size();
    heap_.push_back(NULL);
    while (index) {
        U32 parent = (index - 1) >> 1;
        AMAFRecord * rep = heap_[parent];
        if (*record < *rep)
            break;
        rep->index_ = index;
        heap_[index] = rep;
        index = parent;
    }
    record->index_ = index;
    heap_[index] = record;
}


void AMAFTable::pop_heap() {
#ifdef DEBUG_
    assert(heap_.size() > 0);
#endif
    //printf("In pop_heap\n");
    //fflush(stdout);
    AMAFRecord * record = heap_.back();
    heap_.pop_back();
    U32 size = heap_.size();
    U32 first_leaf = (size >> 1);
    U32 index = 0;
    while (index < first_leaf) {
        U32 left = (index << 1) + 1;
        U32 right = left + 1;
        AMAFRecord * child;
        if (right < size && *heap_[left] < *heap_[right]) {
            child = heap_[right];
            if (*record < *child) {
                child->index_ = index;
                heap_[index] = child;
                index = right;
            } else {
                break;
            }
        } else {
            child = heap_[left];
            if (*record < *child) {
                child->index_ = index;
                heap_[index] = child;
                index = left;
            } else {
                break;
            }
        }
    }
    record->index_ = index;
    heap_[index] = record;
}


void AMAFTable::heap_update(AMAFRecord * record, bool win) {
    //printf("In heap_update\n");
    //fflush(stdout);
    U32 index = record->index_;
    if (win) {
        //printf("Win\n");
        //fflush(stdout);
        while (index) {
            U32 parent = (index - 1) >> 1;
            AMAFRecord * rep = heap_[parent];
            if (*record < *rep)
                break;
            rep->index_ = index;
            heap_[index] = rep;
            index = parent;
        }
        record->index_ = index;
        heap_[index] = record;
    } else {
        //printf("Loss\n");
        //fflush(stdout);
        U32 size = heap_.size();
        U32 first_leaf = (size >> 1);
        while (index < first_leaf) {
            U32 left = (index << 1) + 1;
            U32 right = left + 1;
            AMAFRecord * child;
            if (right < size && *heap_[left] < *heap_[right]) {
                child = heap_[right];
                if (*record < *child) {
                    child->index_ = index;
                    heap_[index] = child;
                    index = right;
                } else {
                    break;
                }
            } else {
                child = heap_[left];
                if (*record < *child) {
                    child->index_ = index;
                    heap_[index] = child;
                    index = left;
                } else {
                    break;
                }
            }
        }
        record->index_ = index;
        heap_[index] = record;
    }
}


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

    void to_str(char * str) const;

    bool operator < (const Move& other) const;
};


struct Snapshot {
    U32 n_open_;
    U64 open_mask_;
    U32 open_[N_CELLS];
    U32 adj_[N_CELLS][MAX_DEGREE];
    U32 n_adj_[N_CELLS];
    U32 effective_adj_[N_CELLS];
    Value value_[N_CELLS];

    Value stones_[2][N_STONES];
	U32 stone_masks_[2];
    U32 n_stones_[2];
    Value * power_[2];
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

class MoveInfo;

class Position {
public:
    List open_;
    U64 open_mask_;
    List adj_[N_CELLS];
    U32 effective_adj_[N_CELLS];
    Value value_[N_CELLS];

	Value stones_[2][N_STONES];
	U32 stone_masks_[2];
    U32 stone_loc_[2][N_STONES+1];
    U32 n_stones_[2];
	Value * power_[2];
    
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
   
    inline void fill(U32 cell_id, Value stone_value) {
        List& adj = adj_[cell_id];
        for (U32 i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            value_[adj_id] += stone_value;
            adj_[adj_id].remove(cell_id);
        }
    }


    inline void unfill(U32 cell_id, Value stone_value) {
        List& adj = adj_[cell_id];
        for (U32 i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            value_[adj_id] -= stone_value;
            adj_[adj_id].readd(cell_id);
        }
    }

    bool is_legal(const Move& move);

    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    Move get_random_move();

    Move get_default_policy_move();

    Move get_expectation_maximising_move();

    Real calculate_expectation() const;		

    Move get_best_winning_move();

    bool is_reasonable_move(U32 cell_index, U32 stone_index) const;

    void get_all_reasonable_moves(std::vector<Move>& moves);

    void get_all_reasonable_moves_with_stone(U32 stone_index, std::vector<Move>& moves);

    bool add_solver_move(const Move& move, std::vector<MoveInfo>& moves);

    bool get_solver_moves(std::vector<MoveInfo>& moves, bool ignore_killer);

    void get_untried_move(U32& cell_index, U32& stone_index, AMAFTable& amaf);

    U32 solve(long long end_time, U32& counter, U32& hash_hits,
            HashTable& table, HashInfo& hash_info, U64 stone_masks[2]);

    std::pair<U32, Move> get_optimal_move(
            long long end_time, U32& counter, U32& hash_hits, bool break_ties, HashTable& table);

    bool all_adj_dead(U32 cell_id);

    void set_alpha(Value alpha);

    void take_snapshot();

    void restore_snapshot();

    void print(FILE * f);
        
    inline std::pair<U32, U32> get_cell_and_stone_indices(const Move& move) {
        return std::pair<U32, U32>(open_.loc_[move.cell_],
                                         stone_loc_[turn_][m_ * move.stone_value_]);
    }

    inline void find_stale_cells() {
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

    inline void find_effective_adj() {
        for (U32 i = 0; i < open_.len_; i++) {
            U32 cell_id = open_.val_[i]; 
            U32& effective_adj = effective_adj_[cell_id];
            effective_adj = 0;
            U64 non_stale_adj_mask = ADJ_MASK[cell_id] & open_mask_ &
                    (~stale_[RED]) & (~stale_[BLUE]);
            while (non_stale_adj_mask) {
                non_stale_adj_mask ^= LSB(non_stale_adj_mask);
                effective_adj++;
            }
        }
    }

    inline bool is_dead(U32 cell_id, U32 p) const {
        U32 op = 1 - p;
        Value * other_power = power_[op];
        Value value = value_[cell_id];
        U32 n_adj = std::min(n_stones_[op], adj_[cell_id].len_);
        Value worst_val = parity[p] * (value + other_power[n_adj]);
        return (worst_val >= parity[p] * (alpha_ - OFFSET[p]));
    }

    inline bool is_winning(U32 p) const {
        return n_dead_[p] > n_stones_[1 - p];
    }

    inline U32 get_non_stale_adj(const List& adj, U64 stale, U32 n) const {
        U32 count = 0;
        for (U32 i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            if (! (stale & MASK(adj_id))) {
                if (count == n)
                    return adj_id;
                else
                    count++;
            }
        }
#ifdef DEBUG_
        assert(false);
#endif
        return NO_CELL;
    }

    inline bool is_adj(const List& adj, U32 id2) const {
        for (U32 i = 0; i < adj.len_; i++)
            if (adj.val_[i] == id2)
                return true;
        return false;
    }

    inline bool is_valid(U32 cell_id, bool& duo) const {
        if (effective_adj_[cell_id] > 2)
            return true;
        const List& adj = adj_[cell_id];
        Value val = m_ * value_[cell_id];
        U32 adj_id;
        Value adj_val;
        U64 stale = stale_[RED] | stale_[BLUE];
        if (effective_adj_[cell_id] == 1) {
            adj_id = get_non_stale_adj(adj, stale, 0);
            adj_val = m_ * value_[adj_id];
            if (effective_adj_[adj_id] == 1) {
                duo = true;
                return (val < adj_val) || (val == adj_val && cell_id < adj_id);
            }
            const List& adj_adj = adj_[adj_id];
            U64 we_control;
            if (turn_ == RED)
                we_control = ~controls_;
            else
                we_control = controls_;
            if (we_control & MASK(cell_id))
                return false;
            if (effective_adj_[adj_id] == 2) {
                U32 adj_adj_id = get_non_stale_adj(adj_adj, stale, 0);
                if (adj_adj_id == cell_id)
                    adj_adj_id = get_non_stale_adj(adj_adj, stale, 1);
                if (effective_adj_[adj_adj_id] == 1) {
                    Value adj_adj_val = value_[adj_adj_id];
                    return (val < adj_adj_val) || (val == adj_adj_val && cell_id < adj_adj_id);
                }
            }
            return true;
        }
        if (effective_adj_[cell_id] == 2) {
            U32 adj1_id = get_non_stale_adj(adj, stale, 0);
            if (effective_adj_[adj1_id] == 2) {
                U32 adj2_id = get_non_stale_adj(adj, stale, 1);
                const List& adj2 = adj_[adj2_id];
                if (effective_adj_[adj2_id] == 2 && is_adj(adj2, adj1_id)) {
                    Value other_val = m_ * value_[adj1_id];
                    if (other_val < val || (other_val == val && adj1_id < cell_id))
                        return false;
                    other_val = m_ * value_[adj2_id];
                    if (other_val < val || (other_val == val && adj2_id < cell_id))
                        return false;
                }
            }
        }
        return true;
    }

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


void Move::to_str(char * str) const {
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
        stone_masks_[p] = (1 << N_STONES) - 1;
        power_[p] = STONE_POWER[p][stone_masks_[p]];
        for (i = 0; i < N_STONES; i++) {
            stones_[p][i] = i + 1;
            stone_loc_[p][i+1] = i;
        }
    }
    set_alpha(alpha);
}


bool Position::is_legal(const Move& move) {
    return move.cell_ >= 0 && move.cell_ < N_CELLS &&
            (open_mask_ & MASK(move.cell_)) &&
            abs(move.stone_value_) >= 1 && abs(move.stone_value_) <= N_STONES &&
            (stone_masks_[turn_] & (1 << (abs(move.stone_value_) - 1)));
}


void Position::make_move(const Move& move) {
#ifdef DEBUG_
    if (! is_legal(move)) { 
        print(stderr);
        char move_str[20];
        move.to_str(move_str);
        fprintf(stderr, "%s\n", move_str);
        fflush(stderr);
        assert(false);
    }
#endif
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
    stone_masks_[turn_] ^= (1 << (stone_number - 1));
    power_[turn_] = STONE_POWER[turn_][stone_masks_[turn_]];
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
        U64 adj_mask = MASK(adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    
    for (i = 0; i < open_.len_; i++) {
        U32 open_id = open_.val_[i];
        U64 open_mask = MASK(open_id);
        if (! (dead_[RED] & open_mask) && ! (dead_[BLUE] & open_mask)) {
            if (is_dead(open_id, RED)) {
                n_dead_[RED]++;
                dead_[RED] |= open_mask;
            } else if (is_dead(open_id, BLUE)) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= open_mask;
            }
        }
    }
    find_stale_cells();
    find_effective_adj();
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
    stone_masks_[turn_] |= (1 << (stone_number - 1));
    power_[turn_] = STONE_POWER[turn_][stone_masks_[turn_]];
    n_stones++;
   
    U64 cell_mask = MASK(cell_id);
    open_mask_ |= cell_mask;
    open_.readd(cell_id); 
    
    unfill(cell_id, move.stone_value_);
    if (value_[cell_id] >= alpha_) {
        controls_ &= ~cell_mask;
        n_controls_[RED]++;
    } else {
        controls_ |= cell_mask;
        n_controls_[BLUE]++;
    }
    //if (controls_ & cell_mask)
    //    n_controls_[BLUE]++;
    //else
    //    n_controls_[RED]++;
    U32 p;
    for (p = 0; p < 2; p++) {
        if (is_dead(cell_id, p)) {
            n_dead_[p]++;
            dead_[p] |= cell_mask;
            //if (stale_[p] & cell_mask)
            //    n_stale_[p]++;
        } else {
            dead_[p] &= ~cell_mask;
        }
    }
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = MASK(adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    
    for (i = 0; i < open_.len_; i++) {
        U32 open_id = open_.val_[i];
        if (open_id != cell_id) {
            U64 open_mask = MASK(open_id);
            if ((dead_[RED] & open_mask)) {
                if (! is_dead(open_id, RED)) {
                    n_dead_[RED]--;
                    dead_[RED] ^= open_mask;
                }
            } else if ((dead_[BLUE] & open_mask)) {
                if (! is_dead(open_id, BLUE)) {
                    n_dead_[BLUE]--;
                    dead_[BLUE] ^= open_mask;
                }
            }
        }
    }

    U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
    U64 stale = open_mask_ & (stale_[RED] | stale_[BLUE]) & (~cell_mask);
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
    for (p = 0; p < 2; p++) {
        bool is_stale = true;
        if (dead_[p] & cell_mask) {
            U64 adj_mask = ADJ_MASK[cell_id] & open_mask_;
            if ((adj_mask & dead) != adj_mask)
                is_stale = false;
        } else {
            is_stale = false;
        }
        if (is_stale) {
            stale_[p] |= cell_mask;
            n_stale_[p]++;
        } else {
            stale_[p] &= ~cell_mask;
        }
    }

    find_effective_adj();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


/*Move Position::get_random_move() {
    return Move(open_.val_[rand() % open_.len_],
                m_ * stones_[turn_][rand() % n_stones_[turn_]]);
}*/


inline U32 lottery(U32 * tickets, U32 n_tickets) {
    int winner = rand() % n_tickets;
    U32 i;
    for (i = 0; winner >= 0; i++)
        winner -= tickets[i];
    return i - 1;
}


Move Position::get_default_policy_move() {
    U32 tickets[N_CELLS];
    U32 total_tickets = 0;
    U32 i, j;
    U32 best_stale = NO_CELL;
    U32 cell_id;
    Value best_stale_val = 1000;
    U32 op = 1 - turn_;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        cell_id = open_.val_[i];
        U64 mask = MASK(cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            Value stale_val = m_ * value_[cell_id];
            if (stale_val < best_stale_val) {
                best_stale = i;
                best_stale_val = stale_val;
            }
        } else {
            bool duo = false;
            if (((dead_[op] & mask) || (n_stones_[turn_] > n_dead_[op])) && is_valid(cell_id, duo)) {
                int swing = 0;
                if (dead_[turn_] & mask)
                    swing = -1;
                else if (dead_[op] & mask)
                    swing = 1;
                List& adj = adj_[cell_id];
                for (j = 0; j < adj.len_; j++) {
                    U32 adj_id = adj.val_[j];
                    U64 adj_mask = MASK(adj_id);
                    if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
                        swing++;
                }
                tickets[i] = swing + 2;
                total_tickets += tickets[i];
            }
        }
    }
    if (best_stale != NO_CELL) {
        tickets[best_stale] = 1;
        total_tickets++;
    }
#ifdef DEBUG_
    if (total_tickets == 0)
        print(stderr);
    assert(total_tickets > 0);
#endif
     
    U32 cell_index = lottery(tickets, total_tickets);
    cell_id = open_.val_[cell_index];

    /*Value reqs[MAX_DEGREE];
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
    std::sort(reqs, reqs + n_uncontrolled);*/
    
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    U32 stone_index;
    if (cell_index == best_stale)
        stone_index = 0;
    else {
        total_tickets = 0;
        for (i = n_stale_[op]; i < n_stones; i++) {
            //U32 changes;
            //for (changes = 0; changes < n_uncontrolled && stones[i] <= reqs[changes]; changes++) {}
            tickets[i] = /*20 * (changes + 1) -*/ stones[i];
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


Move Position::get_best_winning_move() {
    U32 winner = RED;
    if (is_winning(RED))
        winner = RED;
    else if (is_winning(BLUE))
        winner = BLUE;
#ifdef DEBUG_
    else
        assert(false);
#endif
    U32 n_stones = n_stones_[turn_];
    Value * stones = stones_[turn_];
    Move best_move;
    Real best_value = -1000.0;
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = stones[j] * m_;
            Move move(cell_id, stone_value);
            make_move(move);
            bool still_winning = is_winning(winner);
            Real value = calculate_expectation();
            unmake_move(move);
            if (still_winning) {
                value *= m_;
                if (value > best_value) {
                    best_move = move;
                    best_value = value;
                }
            }
        }
    }
#ifdef DEBUG_
    assert(best_value > -1000.0);
#endif
    return best_move;
}


bool Position::is_reasonable_move(U32 cell_index, U32 stone_index) const {
    U32 cell_id = open_.val_[cell_index];
    Value stone_number = stones_[turn_][stone_index];
    //Value stone_value = m_ * stone_number;
   
    U32 op = 1 - turn_; 
    U64 mask = MASK(cell_id);
    if (stale_[op] & mask)
        //TODO: allow only the worst stale cell to be filled
        return stone_index == 0;
    if (stale_[turn_] & mask)
        return stone_index == 0 && dead_[turn_] == open_mask_;

    if (stone_index < n_stale_[op])
        return false;
    
    U32 n_stones = n_stones_[turn_];
    if (n_stones == n_dead_[op] && ! (dead_[op] & mask))
        return false;

    bool duo = false;
    if (! is_valid(cell_id, duo))
        return false;

    if (duo) {
        if (stone_index == 0)
            return true;
        U32 adj_id = get_non_stale_adj(adj_[cell_id], stale_[RED] | stale_[BLUE], 0);
        Value adj_value = value_[adj_id] * m_;
        Value prev_stone_number = stones_[turn_][stone_index - 1];
        Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
        return stone_number >= extra_req && prev_stone_number < extra_req;
    }

    return true;
}
 

void Position::get_all_reasonable_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    Value stone_value;
    Value best_stale_val = 1000;
    U32 best_stale = NO_CELL;
    U32 op = 1 - turn_;
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
            if ((dead_[op] & mask) || (n_stones > n_dead_[op])) {
                bool duo = false;
                if (is_valid(cell_id, duo)) {
                    if (duo) {
                        U32 adj_id = get_non_stale_adj(adj_[cell_id], stale_[RED] | stale_[BLUE], 0);
                        Value adj_value = value_[adj_id] * m_;
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
                    } else {
                        for (j = n_stale_[op]; j < n_stones; j++) {
                            stone_value = m_ * stones[j];
                            moves.push_back(Move(cell_id, stone_value));
                        }
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


void Position::get_all_reasonable_moves_with_stone(U32 stone_index, std::vector<Move>& moves) {
    Value stone_number = stones_[turn_][stone_index];
    Value prev_stone_number;
    if (stone_index == 0)
        prev_stone_number = 0;
    else
        prev_stone_number = stones_[turn_][stone_index - 1];
    Value stone_value = m_ * stone_number;
    Value best_stale_val = 1000;
    U32 best_stale = NO_CELL;
    U32 n_stones = n_stones_[turn_];
    U32 op = 1 - turn_;
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        Value cell_value = value_[cell_id] * m_;
        U64 mask = MASK(cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            if (stone_index == 0 && cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            if ((dead_[op] & mask) || (n_stones > n_dead_[op])) {
                bool duo = false;
                if (is_valid(cell_id, duo)) {
                    if (duo) {
                        U32 adj_id = get_non_stale_adj(adj_[cell_id], stale_[RED] | stale_[BLUE], 0);
                        Value adj_value = value_[adj_id] * m_;
                        if (stone_index == 0)
                            moves.push_back(Move(cell_id, stone_value));
                        else {
                            Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
                            if (stone_number >= extra_req && prev_stone_number < extra_req)
                                moves.push_back(Move(cell_id, stone_value));
                        }       
                    } else {
                        if (stone_index >= n_stale_[op])
                            moves.push_back(Move(cell_id, stone_value));
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL)
        moves.push_back(Move(best_stale, stone_value));
}


struct MoveInfo {
    Move move_;
    int kill_;
    int adj_;

    MoveInfo() {}

    MoveInfo(const Move& move, int kill, int adj):
            move_(move), kill_(kill), adj_(adj) {}

    inline bool operator < (const MoveInfo& other) const {
        if (kill_ != other.kill_)
            return kill_ > other.kill_;
        if (move_.stone_value_ != other.move_.stone_value_)
            return abs(move_.stone_value_) > abs(other.move_.stone_value_);
        if (adj_ != other.adj_)
            return adj_ > other.adj_;
        return move_.cell_ < other.move_.cell_;
    }

};


bool Position::add_solver_move(const Move& move, std::vector<MoveInfo>& moves) {
    U32 cell_id = move.cell_;
    U64 mask = MASK(cell_id);
    Value stone_number = m_ * move.stone_value_;
    U32 op = 1 - turn_;
    int swing = adj_[cell_id].len_;
    int kill = 0;
    if (dead_[turn_] & mask)
        kill = -1;
    //else if (dead_[op] & mask)
    //    kill = 1;
    List& adj = adj_[cell_id];
    bool all_adj_isolated = true;
    Value largest_val_req = 0;
    for (U32 j = 0; j < adj.len_; j++) {
        U32 adj_id = adj.val_[j];
        U64 adj_mask = MASK(adj_id);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask)) {
            U32 n_adj = adj_[adj_id].len_;
            if (n_adj != 1)
                all_adj_isolated = false;
            Value adj_value = (power_[op][std::min(n_adj-1, n_stones_[op])] + value_[adj_id]) * m_;
            Value control_req = m_ * (alpha_ - OFFSET[turn_]);
            if (adj_value + stone_number >= control_req) {
                kill++;
                if (control_req - adj_value > largest_val_req)
                    largest_val_req = control_req - adj_value;
            }
        }
    }
    if (kill + n_dead_[turn_] > n_stones_[op]) {
        moves.push_back(MoveInfo(move, kill, swing));
        return true;
    }
    if (all_adj_isolated) {
        U32 stone_index = stone_loc_[turn_][stone_number];
        if (stone_index > n_stale_[op] && stones_[turn_][stone_index-1] >= largest_val_req)
            return false;
    }
    moves.push_back(MoveInfo(move, kill, swing));
    return false;
}


bool Position::get_solver_moves(std::vector<MoveInfo>& moves, bool ignore_killer) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    Value stone_value;
    U32 op = 1 - turn_;
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
            if ((dead_[op] & mask) || (n_stones > n_dead_[op])) {
                bool duo = false;
                if (is_valid(cell_id, duo)) {
                    if (duo) {
                        U32 adj_id = get_non_stale_adj(adj_[cell_id], stale_[RED] | stale_[BLUE], 0);
                        Value adj_value = value_[adj_id] * m_;
                        Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                        if (adj_value + stones[0] < control_req) {
                            Value extra = control_req - adj_value;
                            for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                            if (j < n_stones) {
                                stone_value = m_ * stones[j];
                                if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                                    return true;
                            }
                        }
                        stone_value = m_ * stones[0];
                        if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                            return true;
                    } else {
                        for (j = n_stale_[op]; j < n_stones; j++) {
                            stone_value = m_ * stones[j];
                            if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                                return true;
                        }
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * m_;
        if (add_solver_move(Move(best_stale, stone_value), moves) && ! ignore_killer)
            return true;
    }
    return false;
}



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
    std::vector<MoveInfo> move_infos;
    if (get_solver_moves(move_infos, false))
        return turn_;
    U32 i;
    U32 op = 1 - turn_;
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
    
    std::vector<MoveInfo> move_infos;
    get_solver_moves(move_infos, true);
    U32 i;
    std::sort(move_infos.begin(), move_infos.end());
  
    Move move; 
    std::vector<Move> winning_moves;
    for (i = 0; i < move_infos.size(); i++) {
        move = move_infos[i].move_;
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
        fprintf(stderr, "%u optimal moves\n", (U32)winning_moves.size());
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
    return std::pair<U32, Move>(1 - turn_, move_infos[0].move_);
}


void Position::get_untried_move(U32& cell_index, U32& stone_index, AMAFTable& amaf) {
    U32 n_stones = n_stones_[turn_];
    stone_index = n_stones;
    U32 i;
    while (stone_index > 0) {
        stone_index--;
        std::vector<Move> moves;
        get_all_reasonable_moves_with_stone(stone_index, moves);
        std::vector<MoveInfo> move_infos;
        for (i = 0; i < moves.size(); i++)
            add_solver_move(moves[i], move_infos);
        std::sort(move_infos.begin(), move_infos.end());
        for (i = 0; i < move_infos.size(); i++) {
            Move move = move_infos[i].move_;
            cell_index = open_.loc_[move.cell_];
            if (! amaf.is_played(cell_index, stone_index))
                return;
        }
    }
    cell_index = NO_CELL;
}


bool Position::all_adj_dead(U32 cell_id) {
    List& adj = adj_[cell_id];
    U32 n_adj = adj.len_;
    for (U32 i = 0; i < n_adj; i++) {
        U64 adj_mask = MASK(adj.val_[i]);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
            return false;
    }
    return true;
}


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    dead_[0] = dead_[1] = 0;
    stale_[RED] = stale_[BLUE] = 0;
    n_stale_[RED] = n_stale_[BLUE] = 0;
    U32 i, cell_id, n_adj;
    U64 mask;
    for (i = 0; i < open_.len_; i++) {
        cell_id = open_.val_[i];
        Value value = value_[cell_id];
        n_adj = adj_[cell_id].len_;
        mask = MASK(cell_id);
        if (value >= alpha_) {
            n_controls_[RED]++;
            if (value + power_[BLUE][std::min(n_adj, n_stones_[BLUE])] >= alpha_) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            }
        } else {
            controls_ |= mask;
            n_controls_[BLUE]++;
            if (value + power_[RED][std::min(n_adj, n_stones_[RED])] < alpha_) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= mask;
            }
        }
    }
    find_stale_cells();
    find_effective_adj();
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
        snapshot_.effective_adj_[cell_id] = effective_adj_[cell_id];
        for (j = 0; j < n_adj; j++)
            snapshot_.adj_[cell_id][j] = adj.val_[j];
    }

    for (p = 0; p < 2; p++) {
        U32 n_stones = n_stones_[p];
        snapshot_.n_stones_[p] = n_stones;
        snapshot_.stone_masks_[p] = stone_masks_[p];
        snapshot_.power_[p] = power_[p];
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
        effective_adj_[cell_id] = snapshot_.effective_adj_[cell_id];
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
        stone_masks_[p] = snapshot_.stone_masks_[p];
        power_[p] = snapshot_.power_[p];
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
    fprintf(f, "\nBlue:");
    for (i = 0; i < n_stones_[BLUE]; i++)
        fprintf(f, " %d", stones_[BLUE][i] * parity[BLUE]);
    fprintf(f, "\n");
    fprintf(f, "Alpha = %d\n", alpha_);
    fprintf(f, "dead_[RED] = 0x%llu\n", dead_[RED]);
    fprintf(f, "dead_[BLUE] = 0x%llu\n", dead_[BLUE]);
    fprintf(f, "stale_[RED] = 0x%llu\n", stale_[RED]);
    fprintf(f, "stale_[BLUE] = 0x%llu\n", stale_[BLUE]);
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
    bool all_children_generated_;

    U32 n_children_fully_explored_;
    std::vector<MCTNode*> children_;
          
    bool solved_;
    bool solve_attempted_;
    Move solution_;
    U32 solver_positions_;
    U32 solver_hash_hits_;
    long long solver_time_;

    AMAFTable amaf_;


    void init(MCTNode * parent, const Position& pos, const Move& move, Value alpha);

    void simulate(Position& pos);
    
    MCTNode * select(Position& pos);

    MCTNode * add_child(Position& pos, const Move& move);

    bool is_now_fully_explored();

    bool is_playable();

    void generate_move(U32& cell_index, U32& stone_index, Position& pos);

    //MCTNode * expand(Position& pos);

    bool playout(Position& pos, U32& move_count);

    bool attempt_solve(Position& pos, HashTable& table, long long allowed_time, bool break_ties);

    Move get_most_played_move();

    Move get_highest_value_move(Position& pos);

    void dispose();

    void print(FILE * f);

    void print_pv(FILE * f);

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

    bool no_playable_move();

    MCTNode * select_alpha(const Position& pos);

    void display_roots();

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
    
    expanded_ = false;
    all_children_generated_ = false;
    
    n_children_fully_explored_ = 0;
    children_.clear();
    
    solved_ = false;
    solve_attempted_ = false;
    solver_positions_ = 0;
    solver_hash_hits_ = 0;
    solver_time_ = 0;
    amaf_.init(pos.open_.len_, pos.n_stones_[pos.turn_]);
}


MCTNode * MCTNode::select(Position& pos) {
    MCTNode * best_node = NULL;
    Real best_node_value = -1000.0, child_val;
    U32 i;
    Real log_pp = log((Real)n_playouts_);
    for (i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            std::pair<U32, U32> indices = pos.get_cell_and_stone_indices(child->move_);
            U32 n_amaf = amaf_.get_n_playouts(indices.first, indices.second);
#ifdef DEBUG_
            assert(n_amaf > 0);
#endif
            Real beta = (Real)n_amaf / (n_amaf + child->n_playouts_ + (RAVE_C * n_amaf) * child->n_playouts_);
            Real mc_val = child->val_;
            if (pos.turn_ == BLUE)
                mc_val = 1.0 - val_;
            child_val = (1.0 - beta) * mc_val + beta * amaf_.get_value(indices.first, indices.second) + 
                    UCB_C * sqrt(log_pp / child->n_playouts_);
            //Real child_playouts = (Real)child->n_playouts_;
            //Real ucb_value = child_val + UCB_C * sqrt(log_pp / child_playouts);
            if (child_val > best_node_value) {
                best_node = child;
                best_node_value = child_val;
            }
        }
    }
    
    if (! all_children_generated_) {
        U32 cell_index, stone_index;
        generate_move(cell_index, stone_index, pos);
        // since we don't know in advance whether all children have been generated,
        // we might discover during selection that this node is already fully explored.
        // This shouldn't happen during expansion however
        if (cell_index != NO_CELL) {
            U32 n_amaf_playouts = amaf_.get_n_playouts(cell_index, stone_index); 
            if (n_amaf_playouts > 0)
                child_val = amaf_.get_value(cell_index, stone_index);
            else
                child_val = 0.5;
            child_val += UCB_C * sqrt(log_pp);
            if (child_val > best_node_value) {
                Move move(pos.open_.val_[cell_index],
                        pos.m_ * pos.stones_[pos.turn_][stone_index]);
#ifdef DEBUG_
                if (! pos.is_legal(move)) { 
                    pos.print(stderr);
                    char move_str[20];
                    move.to_str(move_str);
                    fprintf(stderr, "%s\n", move_str);
                    fflush(stderr);
                    assert(false);
                }
#endif
                best_node = add_child(pos, move);
                amaf_.play(cell_index, stone_index);
            }
        }
    }

    if (best_node == NULL) {
#ifdef DEBUG_
        //fprintf(stderr, "select returning null node\n");
        //fflush(stderr);
        assert(children_.size() > 0);
        assert(is_now_fully_explored());
#endif
        fully_explored_ = true;
        return this;
    }
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


bool MCTNode::is_now_fully_explored() {
    return all_children_generated_ &&
            n_children_fully_explored_ == children_.size();
}


bool MCTNode::is_playable() {
    return fully_explored_ || children_.size() > 0;
}


void MCTNode::generate_move(U32& cell_index, U32& stone_index, Position& pos) {
    cell_index = NO_CELL;
    while (true) {
        amaf_.get_best(cell_index, stone_index);
        if (cell_index == NO_CELL || pos.is_reasonable_move(cell_index, stone_index))
            break;
        amaf_.play(cell_index, stone_index);
    }
    
    if (cell_index == NO_CELL) {
        pos.get_untried_move(cell_index, stone_index, amaf_);
        if (cell_index == NO_CELL) {
            all_children_generated_ = true;
            return;
        }
    }
}
    

/*MCTNode * MCTNode::expand(Position& pos) {
    expanded_ = true;
    return select(pos);
}*/


Move simulation[N_CELLS];

bool MCTNode::playout(Position& pos, U32& move_count) {
#ifdef DEBUG_
    if (pos.is_winning(RED) || pos.is_winning(BLUE)) {
        pos.print(stderr);
        for (U32 m = 0; m < move_count; m++) {
            char move_str[20];
            simulation[m].to_str(move_str);
            fprintf(stderr, "%s\n", move_str);
        }
        fprintf(stderr, "alpha = %d\n", alpha_);
        fprintf(stderr, "n_dead_[RED] = %u, mask = %llx\n", pos.n_dead_[RED], pos.dead_[RED]);
        fprintf(stderr, "n_dead_[BLUE] = %u, mask = %llx\n", pos.n_dead_[BLUE], pos.dead_[BLUE]);
        fflush(stderr);
    }
    assert(! pos.is_winning(RED) && ! pos.is_winning(BLUE));
#endif
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
        simulation[move_count] = pos.get_default_policy_move();
        pos.make_move(simulation[move_count]);
        move_count++;
    }
    if (! result_found) {
#ifdef DEBUG_
        //fprintf(stderr, "Got to end of game\n");
#endif
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    }
    pos.restore_snapshot();
    return result;
}


void MCTNode::simulate(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    MCTNode * search_root = this;
    MCTNode * new_search_root;
    U32 move_count = 0, depth = 0, i;
    while (search_root->expanded_) {
        new_search_root = search_root->select(pos);
#ifdef DEBUG_
        assert(new_search_root != NULL);
#endif
        if (new_search_root == search_root)
            break;
        search_root = new_search_root;
#ifdef DEBUG_
        if (! pos.is_legal(search_root->move_)) { 
            pos.print(stderr);
            char move_str[20];
            search_root->move_.to_str(move_str);
            fprintf(stderr, "%s\n", move_str);
            fflush(stderr);
            assert(false);
        }
#endif
        simulation[move_count] = search_root->move_;
        pos.make_move(search_root->move_);
        move_count++;
        depth++;
    }
    
    U32 result = NO_RESULT;
    if (! search_root->fully_explored_) {
        result = search_root->playout(pos, move_count);
#ifdef DEBUG_
        if (! pos.is_legal(simulation[depth])) { 
            pos.print(stderr);
            char move_str[20];
            simulation[depth].to_str(move_str);
            fprintf(stderr, "%s\n", move_str);
            fflush(stderr);
            assert(false);
        }
#endif
        search_root->expanded_ = true;
        search_root = search_root->add_child(pos, simulation[depth]);
        pos.make_move(simulation[depth]);
        depth++;
    } else {
        result = (search_root->val_ == 1.0);
    }

    while (search_root != NULL) {
        U32 last_move = std::min(move_count, depth + AMAF_HORIZON);
        for (i = depth; i < last_move; i += 2) {
            std::pair<int, int> indices = pos.get_cell_and_stone_indices(simulation[i]);
            search_root->amaf_.update(indices.first, indices.second, result != pos.turn_);
        }
        
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (i = 0; i < search_root->children_.size(); i++)
                if (search_root->children_[i]->val_ * pos.m_ > search_root->val_)
                    search_root->val_ = search_root->children_[i]->val_ * pos.m_;
            search_root->val_ *= pos.m_;
            result = (search_root->val_ == 1.0);
        }
        if (! search_root->fully_explored_) {
#ifdef DEBUG_
            assert(result != NO_RESULT);
#endif
            search_root->n_red_wins_ += result;
            search_root->n_playouts_++;
            search_root->val_ = Real(search_root->n_red_wins_) / Real(search_root->n_playouts_);
        }
        if (parent != NULL) {
            pos.unmake_move(search_root->move_);
            if (search_root->fully_explored_) {
                parent->n_children_fully_explored_++;
                if ((parent->is_now_fully_explored()) ||
                        (pos.turn_ == RED && search_root->val_ == 1.0) ||
                        (pos.turn_ == BLUE && search_root->val_ == 0.0)) {
                    parent->fully_explored_ = true;
                }
            }
        }
        search_root = parent;
        depth--;
    }
}


Move MCTNode::get_highest_value_move(Position& pos) {
    if (solved_)
        return solution_;
    if (fully_explored_ && children_.size() == 0)
        return pos.get_best_winning_move();
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


void MCTNode::print_pv(FILE * f) {
    fprintf(f, "PV:");
    MCTNode * curr = this;
    char move_str[20];
    while (curr->children_.size() > 0) {
        U32 most_visits = 0;
        MCTNode * best_child = NULL;
        for (U32 i = 0; i < curr->children_.size(); i++) {
            MCTNode * child = curr->children_[i];
            if (best_child == NULL || child->n_playouts_ > most_visits) {
                best_child = child;
                most_visits = child->n_playouts_;
            }
        }
        if (best_child == NULL) {
#ifdef DEBUG_
            fprintf(stderr, "PV child was NULL\n");
#endif
            break;
        }
        best_child->move_.to_str(move_str);
        fprintf(f, " %s", move_str);
        curr = best_child;
    }
    fprintf(f, "\n");
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


bool MCTSearch::no_playable_move() {
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++)
        if (roots_[i] != NULL && roots_[i]->is_playable())
            return false;
    return true;
}


MCTNode * MCTSearch::get_or_make_root(Value alpha, const Position& pos) {
    Value index = alpha + MAX_RESULT;
#ifdef DEBUG_
    assert(index >= 0 && index <= 2 * MAX_RESULT);
#endif
    if (roots_[index] == NULL) {
        roots_[index] = get_free();
        roots_[index]->init(NULL, pos, Move(), alpha);
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
            if (test_root != NULL && test_root->is_playable()) {
                if (test_root->val_ >= CHOOSE_TARGET_THRESH[RED] || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    } else {
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL && test_root->is_playable()) {
                if (test_root->val_ < CHOOSE_TARGET_THRESH[BLUE] || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    }
    return root;
}


void MCTSearch::display_roots() {
    for (Value alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
        MCTNode * root = roots_[alpha + MAX_RESULT];
        if (root != NULL)
            root->print(stderr);
    }
}


Move MCTSearch::get_best_move(Position& pos) {
    Real time_left_r = (Real) time_left;
    Real use_proportion;
    int moves_left = ((int)pos.open_.len_ - 18) / 2;
    if (moves_left < 1)
        use_proportion = 0.5;
    else {
        fprintf(stderr, "Expecting to make %d more moves before position solved\n", moves_left);
#ifdef UNIFORM_TM
        double use = (time_left_r - time_limit / 8.0) / (Real)moves_left;
        use_proportion = use / time_left_r;
#else
        use_proportion = 1.0 - pow(((Real)time_limit / 8.0)  / time_left_r, 1.0 / (Real)moves_left);
#endif
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
    
    init_amaf_free_list();
        
    while (time_now - start_micros < alpha_time || no_playable_move()) {
        pos.set_alpha(current_alpha_);
        root = get_or_make_root(current_alpha_, pos);
        if (pos.open_.len_ <= solver_start && ! root->solve_attempted_ && time_now - start_micros < solver_time)
            root->attempt_solve(pos, table_, solver_time - time_now + start_micros, false);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
            if (roots_[current_alpha_ + MAX_RESULT] != NULL &&
                    roots_[current_alpha_ + MAX_RESULT]->fully_explored_)
                break;
        } else {
            if (root->val_ >= TARGET_INCREMENT_THRESH[pos.turn_] && current_alpha_ < MAX_RESULT)
                current_alpha_++;
            if (root->val_ <= TARGET_DECREMENT_THRESH[pos.turn_] && current_alpha_ > MIN_RESULT)
                current_alpha_--;
        }
        time_now = get_time();
    }

    root = select_alpha(pos);
    display_roots();
     
    Value original_alpha = current_alpha_;
    fprintf(stderr, "Alpha = %d\n", current_alpha_);
    pos.set_alpha(current_alpha_);
    bool done = false;
    while (time_now - start_micros < move_time) {
        root = get_or_make_root(current_alpha_, pos);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
            pos.set_alpha(current_alpha_);
            if (roots_[current_alpha_ + MAX_RESULT] != NULL &&
                    roots_[current_alpha_ + MAX_RESULT]->fully_explored_) {
                done = true;
                break;
            }
        }
        time_now = get_time();
    }

    current_alpha_ = original_alpha;
    root = roots_[current_alpha_ + MAX_RESULT];
    if (done)
        root = select_alpha(pos);
    display_roots();
    pos.set_alpha(current_alpha_);

    if (root->solved_) {
        long long time_remaining = move_time - (time_now - start_micros);
        if (time_remaining > 0) {
            if (root->attempt_solve(pos, table_, time_remaining, true))
                fprintf(stderr, "Chose most challenging optimal move\n");
        }
    }

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(get_time() - time_started) / 1e6);
        
    if (root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return root->get_highest_value_move(pos);
    }
    
    root->print_pv(stderr);
    return root->get_most_played_move();
}


int main(int argc, char** argv) {
    time_started = get_time();
    if (argc > 1) {
        int millis;
        if (sscanf(argv[1], "-t=%d", &millis) == 1) {
            time_limit = 1000LL * (long long)millis;
            time_left = (49 * time_limit) / 50;
        }
    }
	init();
	init_free_list();
    init_amaf_free_list();
    	
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
        fprintf(stderr, "Time left: %.2f seconds\n", (float)(time_left - (get_time() - time_started)) / 1e6);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES);
#endif
        MCTSearch search(pos, alpha);
        Move move = search.get_best_move(pos);
        alpha = search.current_alpha_;
        fprintf(stderr, "Used %u nodes\n", N_MCT_NODES - n_free);
        fprintf(stderr, "Used %u AMAF records\n", amaf_pointer);
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
// EventHorizon3_0_0: generate only selected moves during expansion phase. Be more optimistic in target selection. 
// EventHorizon3_0_1: time usage
// EventHorizon3_0_2: tune UCB_C
// EventHorizon3_0_3: make opening heuristic calculation more efficient and more optimistic
// EventHorizon4_0_0: implement RAVE, modify default policy, alter solver move ordering 
// EventHorizon4_0_1: modify memory allocation to be more time and space efficient 
// EventHorizon4_0_2: fix bugs related to accelerated win detection
// EventHorizon4_1_1: improve solver and (hopefully) fix some crashes. 
// EventHorizon4_1_7: fix bugs relating to accessing invalid power array entries and root selection 
// EventHorizon4_2_2: more intelligently detect unreasonable moves 
// EventHorizon4_2_3: apply more intelligent unreasonable move detection to solver 
// EventHorizon4_2_4: tidy up a bit and inline a few functions 
