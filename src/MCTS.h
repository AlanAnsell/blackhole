#include "globals.h"
#include "Position.h"


class MCTNode {
public:
    Move move_;
    MCTNode * parent_;

    Real val_;
    Real result_sum_;
    size_t n_playouts_;

    bool fully_explored_;
    bool expanded_;

    size_t n_children_fully_explored_;
    std::vector<MCTNode*> children_;
    Bitmask table_[N_CELLS];
    
    //MCTNode ** table_[N_CELLS];
    //std::vector<std::pair<Real, Move>> untried_moves_;
    //std::vector<Move> untried_moves_;
    unsigned char n_cell_moves_[N_CELLS];
    size_t n_untried_moves_;

    void init(MCTNode * parent, const Position& pos, const Move& move);

    void ucb(Position& pos);
    
    MCTNode * select(Position& pos);

    MCTNode * expand(Position& pos);

    Value light_playout(Position& pos);

    Move get_best_move(Position& pos);

    Move get_most_played_move();

    Move get_highest_value_move(Position& pos);

    void dispose();

};


extern MCTNode node_store[N_MCT_NODES];
extern MCTNode * free_list[N_MCT_NODES];
extern size_t n_free;

MCTNode * get_free();

void put_free(MCTNode * node);

void init_free_list();
