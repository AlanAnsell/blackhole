#include "MCTS.h"


void MCTNode::dispose() {
    for (size_t i = 0; i < children_.size(); i++) {
        children_[i]->dispose();
        put_free(children_[i]);
    }
}


void MCTNode::init(MCTNode * parent, const Move& move, bool terminal, Real val) {
    move_= move;
    parent_ = parent;
    fully_explored_ = terminal;
    if (terminal)
        val_ = val;
    else
        val_ = 0.0;
    result_sum_ = 0.0;
    n_playouts_ = 0;
    children_.clear();
    untried_moves_.clear();
    expanded_ = false;
    n_children_fully_explored_ = 0;
    //for (size_t i = 0; i < N_CELLS; i++)
    //    table_[i] = NULL;
}


MCTNode * MCTNode::select(Position& pos) {
    // Should only be called if all moves have been tried
    MCTNode * best_node = NULL;
    Real best_value = -1000.0;
    for (size_t i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            Real ucb_value = child->val_ * parity[pos.turn_] + UCB_C * sqrt(log(child->n_playouts_) / n_playouts_);
            if (ucb_value > best_value) {
                best_node = child;
                best_value = ucb_value;
            }
        }
    }
    return best_node;
}


MCTNode * MCTNode::expand(Position& pos) {
    if (! expanded_) {
        pos.get_moves_with_heuristic(untried_moves_);
        std::sort(untried_moves_.begin(), untried_moves_.end());
        expanded_ = true;
    }
    Real val = untried_moves_.back().first * parity[pos.turn_];
    Move move = untried_moves_.back().second;
    untried_moves_.pop_back();
    MCTNode * new_node = get_free();
    new_node->init(this, move, pos.open_.len_ == 2, val);
    children_.push_back(new_node);
    /*if (new_node->fully_explored_) {
        fprintf(stderr, "Added new terminal node to tree\n");
        fflush(stdout);
    }*/
    return new_node;
}


Move playout_moves[N_CELLS];

Value MCTNode::light_playout(Position& pos) {
    int move_count;
    //char move_str[10];
    for (move_count = 0; pos.open_.len_ > 1; move_count++) {
        playout_moves[move_count] = pos.get_random_move();
        //playout_moves[move_count].to_str(move_str);
        //fprintf(stderr, "Making move %s\n", move_str);
        //fflush(stderr);
        pos.make_move(playout_moves[move_count]);
    }
    //fprintf(stderr, "Playout reached end of game\n");
    //fflush(stderr);
    Value result = pos.cells_[pos.open_.val_[0]].value_;
    //fprintf(stderr, "Playout result = %d\n", result);
    //fflush(stderr);
    for (move_count--; move_count >= 0; move_count--) {
        //playout_moves[move_count].to_str(move_str);
        //fprintf(stderr, "Unmaking move %s\n", move_str);
        //fflush(stderr);
        pos.unmake_move(playout_moves[move_count]);
    }
    //fprintf(stderr, "Finished unmaking moves\n");
    //fflush(stderr);
    return result;
}


void MCTNode::ucb(Position& pos) {
    MCTNode * search_root = this;
    while (search_root->expanded_ && search_root->untried_moves_.size() == 0) {
        //fprintf(stderr, "selecting\n");
        //fflush(stderr);
        search_root = search_root->select(pos);
#ifdef DEBUG_
        assert(search_root != NULL);
#endif
        pos.make_move(search_root->move_);
    }
    //fprintf(stderr, "expanding\n");
    //fflush(stderr);
    search_root = search_root->expand(pos);
    //char move_str[10];
    //search_root->move_.to_str(move_str);
    //fprintf(stderr, "Expanded move: %s\n", move_str);
    //fflush(stderr);
    /*if (search_root->fully_explored_) {
        fprintf(stderr, "About to move into final position\n");
        pos.print(stderr);
        char move_str[10];
        search_root->move_.to_str(move_str);
        fprintf(stderr, "Move: %s\n", move_str);
        fflush(stderr);
    }*/
    pos.make_move(search_root->move_);
    /*if (search_root->fully_explored_) {
        fprintf(stderr, "Should now be in final position\n");
        fflush(stderr);
    }*/
    Value result = 0;
    if (! search_root->fully_explored_)
        result = search_root->light_playout(pos);
    //fprintf(stderr, "Playout completed with result %d\n", result);
    //fflush(stderr);
    while (search_root != NULL) {
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (size_t i = 0; i < search_root->children_.size(); i++)
                if (search_root->children_[i]->val_ * parity[pos.turn_] > search_root->val_)
                    search_root->val_ = search_root->children_[i]->val_ * parity[pos.turn_];
            search_root->val_ *= parity[pos.turn_];
        }
        if (! search_root->fully_explored_) {
            search_root->result_sum_ += result;
            search_root->n_playouts_++;
            search_root->val_ = search_root->result_sum_ / search_root->n_playouts_;
        }
        if (parent !=  NULL) {
            pos.unmake_move(search_root->move_);
            if (search_root->fully_explored_) {
                parent->n_children_fully_explored_++;
                //fprintf(stderr, "Checking if search root fully explored\n");
                //fflush(stderr);
                if ((parent->untried_moves_.size() == 0) &&
                        (parent->n_children_fully_explored_ == parent->children_.size())) {
                    parent->fully_explored_ = true;
                    //fprintf(stderr, "Search root fully explored\n");
                    //fflush(stderr);
                }
            }
        }
        search_root = parent;
    }
}


Move MCTNode::get_highest_value_move(Position& pos) {
    Move best_move;
    Real best_val = -1000.0;
    fprintf(stderr, "Choosing highest value move\n");
    pos.print(stderr);
    for (size_t i = 0; i < children_.size(); i++) {
        char move_str[10];
        children_[i]->move_.to_str(move_str);
        fprintf(stderr, "%s: %lf\n", move_str, children_[i]->val_);
        if (children_[i]->val_ * parity[pos.turn_] > best_val) {
            best_move = children_[i]->move_;
            best_val = children_[i]->val_ * parity[pos.turn_];
        }
    }
    return best_move;
}


Move MCTNode::get_most_played_move() {
    Move best_move;
    size_t most_visits = 0;
    for (size_t i = 0; i < children_.size(); i++) {
        if (children_[i]->n_playouts_ > most_visits) {
            best_move = children_[i]->move_;
            most_visits = children_[i]->n_playouts_;
        }
    }
    return best_move;
}


Move MCTNode::get_best_move(Position& pos) {
    for (size_t i = 0; i < 5000 && ! fully_explored_; i++)
        ucb(pos);

    if (fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return get_highest_value_move(pos);
    }

    return get_most_played_move();
}


MCTNode node_store[N_MCT_NODES];
MCTNode * free_list[N_MCT_NODES];
size_t n_free;


MCTNode * get_free() {
#ifdef DEBUG_
    assert(n_free > 0);
#endif
    return free_list[--n_free];
}


void put_free(MCTNode * node) {
    free_list[n_free++] = node;
}


void init_free_list() {
    for (size_t i = 0; i < N_MCT_NODES; i++)
        free_list[i] = &node_store[i];
    n_free = N_MCT_NODES;
}






