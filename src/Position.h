#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"


class List {
public:
    size_t * val_;
    size_t * loc_;
    size_t len_;

    List(): val_(NULL), loc_(NULL), len_(0) {}

    void init(size_t max_len, size_t range) {
        len_ = 0;
        val_ = new size_t[max_len];
        loc_ = new size_t[range];
    }

    ~List() {
        delete [] val_;
        delete [] loc_;
    }

    inline void add(size_t val) {
        loc_[val] = len_;
        val_[len_] = val;
        len_++;
    }

    inline void remove(size_t val) {
        len_--;
        if (len_) {
            size_t end_val = val_[len_];
            size_t val_loc = loc_[val];
            loc_[end_val] = val_loc;
            val_[val_loc] = end_val;
        }
    }

    inline void readd(size_t val) {
        size_t old_loc = loc_[val];
        size_t new_val = val_[old_loc];
        loc_[new_val] = len_;
        val_[len_] = new_val;
        val_[old_loc] = val;
        len_++;
    }

};


class Move {
public:
    CellID cell_;
    Value stone_value_;

    Move(): cell_(0), stone_value_(0) {}

    Move(CellID cell, Value stone_value): cell_(cell), stone_value_(stone_value) {}

    Move(char * move_str, size_t turn);

    void to_str(char * str);

    bool operator < (const Move& other) const;
};


struct Snapshot {
    size_t n_open_;
    int64 open_mask_;
    size_t open_[N_CELLS];
    size_t adj_[N_CELLS][MAX_DEGREE];
    size_t n_adj_[N_CELLS];
    Value value_[N_CELLS];

    Value stones_[2][N_STONES];
    size_t n_stones_[2];
    size_t turn_;

    int64 controls_;
    size_t n_controls_[2];

    int64 dead_[2];
    size_t n_dead_[2];

    int64 stale_[2];
    size_t n_stale_[2];
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
    int64 red_;
    int64 blue_;
    size_t result_;
};


class HashTable {
public:
    size_t n_records_[HASH_SIZE];
    HashRecord table_[HASH_SIZE][BUCKET_SIZE];
    
    void init() {
        for (size_t i = 0; i < HASH_SIZE; i++)
            n_records_[i] = 0;
    }

    void add(int64 red, int64 blue, size_t result) {
        int64 hash = HASH_BITMASKS(red, blue);
        //fprintf(stderr, "Hash: %llu\n", hash);
        size_t& n_records = n_records_[hash];
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

    size_t find(int64 red, int64 blue) {
        int64 hash = HASH_BITMASKS(red, blue);
        size_t n_records = n_records_[hash];
        HashRecord * records = table_[hash];
        for (size_t i = 0; i < n_records; i++) {
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
    int64 open_mask_;
    List adj_[N_CELLS];
    Value value_[N_CELLS];

	Value stones_[2][N_STONES];
	size_t stone_loc_[2][N_STONES+1];
    size_t n_stones_[2];
	size_t turn_;
    int m_;

    Value alpha_;
    int64 controls_;
    size_t n_controls_[2];

    int64 dead_[2];
    size_t n_dead_[2];

    int64 stale_[2];
    size_t n_stale_[2];

    Snapshot snapshot_;


	Position(CellID blocked[N_BLOCKED_CELLS], Value alpha);
   
    void fill(size_t cell_id, Value stone_value);

    void unfill(size_t cell_id, Value stone_value);

    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    void get_stone_power(size_t p, Value * power);

    Move get_random_move();

    Move get_default_policy_move();

    Move get_expectation_maximising_move();

	Real calculate_win_prob(Value alpha) const;
  
    std::pair<Real, Move> get_best_alpha_move(Value alpha);

    std::pair<Real, Move> get_best_move();

    Real calculate_expectation() const;		

    std::pair<Real, Move> search_expectation(size_t depth, Real a, Real b);

    Real get_control_heuristic() const;

    void get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves);

    void get_all_moves(std::vector<Move>& moves);

    void get_all_reasonable_moves(std::vector<Move>& moves);

    size_t solve(long long end_time, size_t& counter, size_t& hash_hits,
            HashTable& table, HashInfo& hash_info, int64 stone_masks[2]);

    std::pair<size_t, Move> get_optimal_move(
            long long end_time, size_t& counter, size_t& hash_hits, bool break_ties, HashTable& table);

    bool is_dead(size_t cell_id, size_t p, Value * other_power);

    bool all_adj_dead(size_t cell_id);

    void find_stale_cells();

    bool is_winning(size_t p) const;

    void set_alpha(Value alpha);

    void take_snapshot();

    void restore_snapshot();

    std::pair<size_t, size_t> get_cell_and_stone_indices(const Move& move);

    void print(FILE * f);    
};

class HashInfo {
public:
    size_t cell_index_[N_CELLS];
    size_t stone_shift_[2][N_STONES+1];

    HashInfo(const Position& pos) {
        size_t i;
        for (i = 0; i < pos.open_.len_; i++)
            cell_index_[pos.open_.val_[i]] = i + 1;
        for (size_t p = 0; p < 2; p++) {
            size_t n_stones = pos.n_stones_[p];
            const Value * stones = pos.stones_[p];
            for (size_t i = 0; i < n_stones; i++)
                stone_shift_[p][stones[i]] = 5 * i;
        }
    }

};


#endif
