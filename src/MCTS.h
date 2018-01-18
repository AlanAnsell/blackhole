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


extern MCTNode node_store[N_MCT_NODES];
extern MCTNode * free_list[N_MCT_NODES];
extern U32 n_free;

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

#endif
