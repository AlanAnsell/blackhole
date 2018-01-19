#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"
#include "AMAF.h"

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
   
    void fill(U32 cell_id, Value stone_value);

    void unfill(U32 cell_id, Value stone_value);

    bool is_legal(const Move& move);

    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    Move get_random_move();

    Move get_default_policy_move();

    Move get_expectation_maximising_move();

    void find_n_ways(U32 stone_index) const;

	Real calculate_win_prob(Value alpha, U32 stone_index) const;

    Move get_best_winning_move();
  
    std::pair<Real, Move> get_best_alpha_move(Value alpha);

    std::pair<Real, Move> get_best_move();

    Real calculate_expectation() const;		

    std::pair<Real, Move> search_expectation(U32 depth, Real a, Real b);

    Real get_control_heuristic() const;

    void get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves);

    void get_all_moves(std::vector<Move>& moves);

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

    void find_stale_cells();

    void find_effective_adj();

    void set_alpha(Value alpha);

    void take_snapshot();

    void restore_snapshot();

    std::pair<U32, U32> get_cell_and_stone_indices(const Move& move);

    void print(FILE * f);
        
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


#endif
