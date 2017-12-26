#include "MCTS.h"


void MCTNode::dispose() {
    for (size_t i = 0; i < children_.size(); i++) {
        children_[i]->dispose();
        put_free(children_[i]);
    }
}


void MCTNode::init(MCTNode * parent, const Position& pos, const Move& move, Value alpha) {
    move_= move;
    parent_ = parent;
    alpha_ = alpha;
    fully_explored_ = pos.open_.len_ == 1;
    if (fully_explored_)
        val_ = (pos.cells_[pos.open_.val_[0]].value_ >= alpha_) ? 1.0 : 0.0;
    else
        val_ = 0.5;
    n_red_wins_ = 0;
    n_playouts_ = 0;
    children_.clear();
    //untried_moves_ = pos.open_.len_ * pos.n_stones_[pos.turn_];
    expanded_ = false;
    n_children_fully_explored_ = 0;
}


MCTNode * MCTNode::select(Position& pos) {
    MCTNode * best_node = NULL;
    Real best_value = -1000.0;
    bool fe = false;
    for (size_t i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            fe = true;
            Real child_val = child->val_;
            if (pos.turn_ == BLUE)
                child_val = 1.0 - child_val;
            Real parent_playouts = std::max(1.0, (Real)n_playouts_);
            Real child_playouts = std::max(1.0, (Real)child->n_playouts_);
            Real ucb_value = child_val + UCB_C * sqrt(log(parent_playouts) / child_playouts);
            if (ucb_value > best_value) {
                best_node = child;
                best_value = ucb_value;
            }
        }
    }
#ifdef DEBUG_
    assert(fe);
#endif
    return best_node;
}


// void MCTNode::fill_move_table(Position& pos) {
//     size_t p = pos.turn_;
//     size_t n_stones = pos.n_stones_[p];
//     Value * stones = pos.stones_[p];
//     //Value power[MAX_DEGREE+1];
// 
//     Bitmask mask = 0;
//     size_t i;
//     for (i = 0; i < n_stones; i++)
//         mask |= (1 << (stones[i] - 1));
// 
//     //pos.get_stone_power(1 - p, alpha_, power);
//     
//     //bool seen_dead = false;
//     n_untried_moves_ = 0;    
//     for (i = 0; i < pos.open_.len_; i++) {
//         CellID cell_id = pos.open_.val_[i];
//         Cell& cell = pos.cells_[cell_id];
//         /*Value worst_val = cell.value_ + power[cell.adj_.len_];
//         if (worst_val * parity[p] >= (alpha_ - OFFSET[p]) * parity[p]) {
//             if (seen_dead) {
//                 table_[i] = 0;
//                 n_cell_moves_[i] = 0;
//             } else {
//                 table_[i] = (1 << (pos.stones_[p][0] - 1));
//                 n_cell_moves_[i] = 1;
//                 n_untried_moves_++;
//                 seen_dead = true;
//             }*/
//         if (cell.adj_.len_ == 0) {
//             table_[i] = (1 << (pos.stones_[p][0] - 1));
//             n_cell_moves_[i] = 1;
//             n_untried_moves_++;
//         } else {
//             table_[i] = mask;
//             n_cell_moves_[i] = (unsigned char)n_stones;
//             n_untried_moves_ += n_stones;
//         }
//     }
// }


void MCTNode::get_children(Position& pos) {
    size_t p = pos.turn_;
    Value * stones = pos.stones_[p];
    size_t n_stones = pos.n_stones_[p];
    for (size_t i = 0; i < pos.open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[p] * stones[j];
            Move move(pos.open_.val_[i], stone_value);
            pos.make_move(move);
            MCTNode * child = get_free();
            child->init(this, pos, move, alpha_);
            children_.push_back(child);
            pos.unmake_move(move);
        }
    }
}


MCTNode * MCTNode::expand(Position& pos) {
    get_children(pos);
    std::random_shuffle(children_.begin(), children_.end());
    /*MCTNode * best_child = NULL;
    Real best_val = -1.0;
    for (size_t i = 0; i < children_.size(); i++) {
        Real child_val = children_[i]->val_;
        if (pos.turn_ == BLUE)
            child_val = 1 - child_val;
        if (child_val > best_val) {
            best_child = children_[i];
            best_val = child_val;
        }
    }*/
    expanded_ = true;
    
    return children_.back();
}


Move playout_moves[N_CELLS];

bool MCTNode::light_playout(Position& pos) {
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
    bool result = (pos.cells_[pos.open_.val_[0]].value_ >= alpha_);
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
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    MCTNode * search_root = this;
    while (search_root->expanded_) {
        //fprintf(stderr, "selecting\n");
        //fflush(stderr);
        search_root = search_root->select(pos);
#ifdef DEBUG_
        assert(search_root != NULL);
#endif
        pos.make_move(search_root->move_);
    }
    
    search_root = search_root->expand(pos);
    pos.make_move(search_root->move_);
    
    size_t result = 0;
    if (! search_root->fully_explored_)
        result = search_root->light_playout(pos);
    
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
            search_root->n_red_wins_ += result;
            search_root->n_playouts_++;
            search_root->val_ = Real(search_root->n_red_wins_) / Real(search_root->n_playouts_);
        }
        if (parent != NULL) {
            pos.unmake_move(search_root->move_);
            if (search_root->fully_explored_) {
                if (search_root->expanded_) {
                    parent->n_children_fully_explored_++;
                    if (parent->n_children_fully_explored_ == parent->children_.size()) {
                        parent->fully_explored_ = true;
                    }
                } else {
                    parent->fully_explored_ = true;
                }
            }
        }
        search_root = parent;
    }
}


Move MCTNode::get_highest_value_move(Position& pos) {
    Move best_move;
    Real best_val = -1000.0;
    size_t n_playouts = 0;
    //fprintf(stderr, "Choosing highest value move\n");
    //pos.print(stderr);
    for (size_t i = 0; i < children_.size(); i++) {
        //char move_str[10];
        //children_[i]->move_.to_str(move_str);
        //fprintf(stderr, "%s: %lf\n", move_str, children_[i]->val_);
        if (children_[i]->val_ * parity[pos.turn_] > best_val) {
            best_move = children_[i]->move_;
            best_val = children_[i]->val_ * parity[pos.turn_];
            n_playouts = children_[i]->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %u playouts\n", n_playouts);
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
    fprintf(stderr, "Selected move has %d playouts\n", most_visits);
    return best_move;
}


Move MCTNode::get_best_move(Position& pos) {
    long long start_micros = get_time();
    long long time_now = start_micros;
    long long move_time = 300000LL;
    while (! fully_explored_ && time_now - start_micros < move_time) {
        for (size_t i = 0; i < 100 && ! fully_explored_; i++) {
            //fprintf(stderr, "Running UCB\n");
            //fflush(stderr);
            ucb(pos);
            //fprintf(stderr, "Finished running UCB\n");
            //fflush(stderr);
        }
        time_now = get_time();
    }
    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts_, (double)(time_now - start_micros) / 1e6);

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


MCTSearch::MCTSearch(const Position& pos, Value alpha) {
    current_alpha_ = alpha;
    for (size_t i = 0; i < 2 * MAX_RESULT + 1; i++)
        roots_[i] = NULL;
}


MCTSearch::~MCTSearch() {
    for (size_t i = 0; i < 2 * MAX_RESULT + 1; i++) {
        if (roots_[i] != NULL) {
            roots_[i]->dispose();
            put_free(roots_[i]);
        }
    }
}


MCTNode * MCTSearch::get_or_make_root(Value alpha, const Position& pos) {
    Value index = alpha + MAX_RESULT;
#ifdef DEBUG_
    assert(index >= 0 && index <= 2 * MAX_RESULT);
#endif
    if (roots_[index] == NULL) {
        roots_[index] = get_free();
        roots_[index]->init(NULL, pos, Move(), alpha);
    }
    return roots_[index];
}


Move MCTSearch::get_best_move(Position& pos) {
    long long start_micros = get_time();
    long long time_now = start_micros;
    long long alpha_time = 100000LL;
    long long move_time = 300000LL;
    size_t n_playouts = 0;
    size_t i;
    MCTNode * root;
    while (time_now - start_micros < alpha_time) {
        //fprintf(stderr, "Iterating\n");
        //fflush(stderr);
        root = get_or_make_root(current_alpha_, pos);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            //fprintf(stderr, "Anout to UCB\n");
            //fflush(stderr);
            root->ucb(pos);
            n_playouts++;
        } 
        //fprintf(stderr, "Finished UCBs\n");
        //fflush(stderr);
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT)
                current_alpha_--;
            else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT)
                current_alpha_++;
            else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
        } else {
            if (root->val_ >= 0.5) {
                if (current_alpha_ < MAX_RESULT)
                    current_alpha_++;
            } else {
                if (current_alpha_ > MIN_RESULT)
                    current_alpha_--;
            }
        }
        time_now = get_time();
    }

    MCTNode * search_root = NULL;
    Value alpha;
    if (pos.turn_ == RED) {
        for (alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
            root = roots_[alpha + MAX_RESULT];
            if (root != NULL) {
                fprintf(stderr, "Alpha = %d (%u): %.4f\n", alpha, root->n_playouts_, root->val_);
                if (root->val_ >= 0.5 || search_root == NULL) {
                    current_alpha_ = alpha;
                    search_root = root;
                }
            }
        }
    } else {
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            root = roots_[alpha + MAX_RESULT];
            if (root != NULL) {
                fprintf(stderr, "Alpha = %d (%u): %.4f\n", alpha, root->n_playouts_, root->val_);
                if (root->val_ < 0.5 || search_root == NULL) {
                    current_alpha_ = alpha;
                    search_root = root;
                }
            }
        }
    }

    fprintf(stderr, "Alpha = %d\n", current_alpha_);

    while (time_now - start_micros < move_time) {
        for (i = 0; i < 100 && ! search_root->fully_explored_; i++) {
            //fprintf(stderr, "Anout to UCB\n");
            //fflush(stderr);
            search_root->ucb(pos);
            n_playouts++;
        } 
        time_now = get_time();
    }

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(time_now - start_micros) / 1e6);

    if (search_root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return search_root->get_highest_value_move(pos);
    }

    return search_root->get_most_played_move();
}






