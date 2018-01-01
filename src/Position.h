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

#endif
