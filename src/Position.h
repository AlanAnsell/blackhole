#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"


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


#endif
