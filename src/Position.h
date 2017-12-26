#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"


class List {
public:
    size_t * val_;
    size_t * loc_;
    size_t len_;

    List(size_t max_len, size_t range): len_(0) {
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

class Cell {
public:
	CellID id_;
	Cell * board_;
	List adj_;
    Value value_;

    Cell(): adj_(MAX_DEGREE, N_CELLS) {}

    void init(CellID id, Cell * cells);

	void remove_neighbour(CellID neighbour_id, Value stone_value);

	void add_neighbour(CellID neighbour_id, Value stone_value);

	void fill(Value stone_value);
	
	void unfill(Value stone_value);
	
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


class Position {
public:
	Cell cells_[N_CELLS];
    List open_;

	Value stones_[2][N_STONES];
	size_t n_stones_[2];
	size_t turn_;

    Value alpha_;
    int64 controls_;
    size_t n_controls_[2];

	Position(CellID blocked[N_BLOCKED_CELLS], Value alpha);
    
    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    void get_stone_power(size_t p, Value alpha, Value * power);

    size_t dead_endgame_solve(Value alpha);
    
    Value dead_endgame_value();
    
    Move get_random_move();

    Move get_expectation_maximising_move_with_endgame_solve();

	Real search_expectation();

	Real calculate_expectation() const;		

    void get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves);

    void get_all_moves(std::vector<Move>& moves);

    void set_alpha(Value alpha);

    void print(FILE * f);    
};

#endif
