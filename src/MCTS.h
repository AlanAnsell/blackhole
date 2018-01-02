#ifndef MCTS_H_
#define MCTS_H_

#include "globals.h"
#include "Position.h"
#include "AMAF.h"


class MCTNode {
public:
    Move move_;
    MCTNode * parent_;

    Value alpha_;
    Real val_;
    size_t n_red_wins_;
    size_t n_playouts_;

    bool fully_explored_;
    bool expanded_;

    size_t n_children_fully_explored_;
    size_t n_children_;
    size_t n_child_moves_;
    std::vector<Move> child_moves_;
    std::vector<MCTNode*> children_;
    
    bool solved_;
    bool solve_attempted_;
    Move solution_;
    size_t solver_positions_;
    size_t solver_hash_hits_;
    long long solver_time_;

    AMAFTable amaf_[2];


    void init(MCTNode * parent, const Position& pos, const Move& move, Value alpha);

    void ucb(Position& pos);
    
    MCTNode * select(Position& pos, AMAFTable& amaf);

    MCTNode *  add_child(Position& pos, const Move& move);

    void get_children(Position& pos);

    MCTNode * expand(Position& pos, AMAFTable& amaf);

    bool light_playout(Position& pos, size_t& move_count);

    void attempt_solve(Position& pos, HashTable& table);

    Move get_most_played_move();

    Move get_highest_value_move(Position& pos);

    void dispose();

    void print(FILE * f);

};


extern MCTNode node_store[N_MCT_NODES];
extern MCTNode * free_list[N_MCT_NODES];
extern size_t n_free;

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

#endif
