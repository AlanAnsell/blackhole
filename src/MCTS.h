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
    std::vector<Move> untried_children_;
    //bool cached_child_;
    //U32 cached_cell_index_;
    //U32 cached_stone_index_;
        
    //bool get_untried_failed_;
    int untried_stone_index_;
    U64 valid_;
    U64 duo_;

    U32 n_children_fully_explored_;
    std::vector<MCTNode*> children_;
          
    bool solved_;
    bool solve_attempted_;
    Move solution_;
    U32 solver_positions_;
    U32 solver_hash_hits_;
    long long solver_time_;

    // tried_ is indexed first by cell ID, then by stone index
    std::vector<U16> tried_;
    AMAFTable amaf_[2];
    AMAFTable * my_amaf_;
    AMAFTable * par_amaf_;

    void init(MCTNode * parent, const Position& pos, Move move, Value alpha, AMAFTable * my_amaf);

    void get_best_from_amaf(U32& cell_id, U32& stone_index, Position& pos);

    void simulate(Position& pos);
    
    MCTNode * select(Position& pos);

    MCTNode * add_child(Position& pos, Move move);

    bool is_now_fully_explored();

    bool is_playable();

    void generate_move(U32& cell_id, U32& stone_index, Position& pos);

    //MCTNode * expand(Position& pos);

    bool playout(Position& pos, U32& move_count);

    bool attempt_solve(Position& pos, HashTable& table, long long allowed_time, bool break_ties);

    Move get_most_played_move();

    Move get_highest_value_move(Position& pos);

    void dispose();

    void print_rave_calc(char * S, Position& pos, MCTNode * child) const;

    void print_most_played_moves(Position& pos, U32 n) const;

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

    MCTSearch(Value alpha);

    ~MCTSearch();

    bool no_playable_move();

    MCTNode * select_alpha(const Position& pos, Real choose_target_thresh);

    void display_roots();

    MCTNode * get_or_make_root(Value alpha, const Position& pos);

    void analyse(Position& pos);

    Move get_best_move(Position& pos);

};

#endif
