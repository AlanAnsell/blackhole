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
        val_ = pos.get_control_heuristic(); 
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


void MCTNode::add_child(Position& pos, const Move& move) {
    pos.make_move(move);
    MCTNode * child = get_free();
    child->init(this, pos, move, alpha_);
    children_.push_back(child);
    pos.unmake_move(move);
}


void MCTNode::get_children(Position& pos) {
    size_t p = pos.turn_;
    Value * stones = pos.stones_[p];
    size_t n_stones = pos.n_stones_[p];
    Value stone_value;
    Value best_single_val = 1000;
    CellID best_single = 400;
    for (size_t i = 0; i < pos.open_.len_; i++) {
        CellID id = pos.open_.val_[i];
        Cell& cell = pos.cells_[id];
        Value cell_value = cell.value_ * parity[p];
        if (cell.adj_.len_ == 0) {
            if (cell_value < best_single_val) {
                best_single = id;
                best_single_val = cell_value;
            }
        } else {
            size_t j;
            bool pair = false;
            if (cell.adj_.len_ == 1) {
                CellID adj_id = cell.adj_.val_[0];
                Cell& adj = pos.cells_[adj_id];
                if (adj.adj_.len_ == 1) {
                    pair = true;
                    Value adj_value = adj.value_ * parity[p];
                    if (cell_value < adj_value || (cell_value == adj_value && id < adj_id)) {
                        stone_value = parity[p] * stones[0];
                        add_child(pos, Move(id, stone_value));
                        Value control_req = parity[p] * (alpha_ - OFFSET[p]);
                        if (adj_value + stones[0] < control_req) {
                            Value extra = control_req - adj_value;
                            for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                            if (j < n_stones) {
                                stone_value = parity[p] * stones[j];
                                add_child(pos, Move(id, stone_value));
                            }
                        }       
                    }
                }
            }
            if (! pair) {
                for (j = 0; j < n_stones; j++) {
                    stone_value = parity[p] * stones[j];
                    add_child(pos, Move(id, stone_value));
                }
            }
        }
    }
    if (best_single != 400) {
        stone_value = stones[0] * parity[p];
        add_child(pos, Move(best_single, stone_value));
    }
}


MCTNode * MCTNode::expand(Position& pos) {
    get_children(pos);
    MCTNode * best_child = NULL;
    Real best_val = -1.0;
    for (size_t i = 0; i < children_.size(); i++) {
        Real child_val = children_[i]->val_;
        if (pos.turn_ == BLUE)
            child_val = 1 - child_val;
        if (child_val > best_val) {
            best_child = children_[i];
            best_val = child_val;
        }
    }
    expanded_ = true;
    
    return best_child ;
}


Move playout_moves[N_CELLS];

bool MCTNode::light_playout(Position& pos) {
    int move_count;
    for (move_count = 0; pos.open_.len_ > 1; move_count++) {
        playout_moves[move_count] = pos.get_default_policy_move();
        //playout_moves[move_count] = pos.get_random_move();
        pos.make_move(playout_moves[move_count]);
    }
    bool result = (pos.cells_[pos.open_.val_[0]].value_ >= alpha_);
    for (move_count--; move_count >= 0; move_count--)
        pos.unmake_move(playout_moves[move_count]);
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
                    if ((parent->n_children_fully_explored_ == parent->children_.size()) ||
                            (pos.turn_ == RED && search_root->val_ == 1.0) ||
                            (pos.turn_ == BLUE && search_root->val_ == 0.0)) {
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
        if ((! fully_explored_ || children_[i]->fully_explored_) &&
                children_[i]->val_ * parity[pos.turn_] > best_val) {
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

MCTNode * MCTSearch::select_alpha(const Position& pos) {
    MCTNode * test_root;
    MCTNode * root = NULL;
    Value alpha;
    if (pos.turn_ == RED) {
        for (alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL) {
                fprintf(stderr, "Alpha = %d (%u): %.4f\n", alpha, test_root->n_playouts_, test_root->val_);
                if (test_root->val_ >= 0.5 || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    } else {
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL) {
                fprintf(stderr, "Alpha = %d (%u): %.4f\n", alpha, test_root->n_playouts_, test_root->val_);
                if (test_root->val_ < 0.5 || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    }
    return root;
}


Move MCTSearch::get_best_move(Position& pos) {
    long long start_micros = get_time();
    long long time_now = start_micros;
    long long alpha_time = 100000LL;
    long long move_time = 500000LL;
    size_t n_playouts = 0;
    size_t i;
    MCTNode * root;
    
    if (pos.open_.len_ > 20)
        return pos.get_expectation_maximising_move();
    
    while (time_now - start_micros < alpha_time) {
        root = get_or_make_root(current_alpha_, pos);
        pos.set_alpha(current_alpha_);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->ucb(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
                if (get_or_make_root(current_alpha_, pos)->fully_explored_)
                    break;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                if (get_or_make_root(current_alpha_, pos)->fully_explored_)
                    break;
            } else {
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

    root = select_alpha(pos);
    
    Value original_alpha = current_alpha_;
    fprintf(stderr, "Alpha = %d\n", current_alpha_);
    pos.set_alpha(current_alpha_);
    bool solved = false;
    while (time_now - start_micros < move_time) {
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->ucb(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
                root = get_or_make_root(current_alpha_, pos);
                pos.set_alpha(current_alpha_);
                if (root->fully_explored_) {
                    solved = true;
                    break;
                }
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                root = get_or_make_root(current_alpha_, pos);
                pos.set_alpha(current_alpha_);
                if (root->fully_explored_) {
                    solved = true;
                    break;
                }
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
        }
        time_now = get_time();
    }

    current_alpha_ = original_alpha;
    if (solved)
        root = select_alpha(pos);
    pos.set_alpha(current_alpha_);

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(time_now - start_micros) / 1e6);

    if (root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return root->get_highest_value_move(pos);
    }

    return root->get_most_played_move();
}






