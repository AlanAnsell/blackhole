#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"


class Cell {
public:
	CellID id_;
	Cell * adj_[MAX_DEGREE];
	size_t loc_[N_CELLS];
	size_t degree_;
	Value value_;


    void init(CellID id, Cell * cells);

	void remove_neighbour(CellID neighbour_id, Value stone_value);

	void add_neighbour(CellID neighbour_id, Cell * neighbour, Value stone_value);

	void fill(Value stone_value);
	
	void unfill(Value stone_value);
	
};


class Move {
public:
    CellID cell_;
    Value stone_value_;

    Move(CellID cell, Value stone_value): cell_(cell), stone_value_(stone_value) {}

    Move(char * move_str, size_t turn);

    void to_str(char * str);
};


class Position {
public:
	Cell cells_[N_CELLS];
	CellID open_cells_[N_CELLS];
	size_t open_loc_[N_CELLS];
	size_t n_open_;
	
	Value stones_[2][N_STONES];
	size_t n_stones_[2];
	size_t turn_;
	
	Position(CellID blocked[N_BLOCKED_CELLS]);
    
    void make_move(const Move& move);
    
	void unmake_move(const Move& move);

    Move get_random_move();

    Move get_expectation_maximising_move();

	Real search_expectation();

	Real calculate_expectation();		

    void print();    
};

#endif
