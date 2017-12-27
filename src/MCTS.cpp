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
    if (pos.is_winning(RED)) {
        fully_explored_ = true;
        val_ = 1.0;
    } else if (pos.is_winning(BLUE)) {
        fully_explored_ = true;
        val_ = 0.0;
    } else {
        fully_explored_ = false;
        val_ = pos.get_control_heuristic(); 
    }
    n_red_wins_ = 0;
    n_playouts_ = 0;
    children_.clear();
    expanded_ = false;
    n_children_fully_explored_ = 0;
}


MCTNode * MCTNode::select(Position& pos) {
    MCTNode * best_node = NULL;
    Real best_value = -1000.0;
#ifdef DEBUG_
    bool fe = false;
#endif
    for (size_t i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
#ifdef DEBUG_
            fe = true;
#endif
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
    Value best_stale_val = 1000;
    size_t best_stale = NO_CELL;
    for (size_t i = 0; i < pos.open_.len_; i++) {
        size_t cell_id = pos.open_.val_[i];
        Value cell_value = pos.value_[cell_id] * parity[p];
        if (pos.stale_ & (1LL << (int64)cell_id)) {
            if (cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            size_t j;
            bool pair = false;
            if (pos.adj_[cell_id].len_ == 1) {
                size_t adj_id = pos.adj_[cell_id].val_[0];
                if (pos.adj_[adj_id].len_ == 1) {
                    pair = true;
                    Value adj_value = pos.value_[adj_id] * parity[p];
                    if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                        stone_value = parity[p] * stones[0];
                        add_child(pos, Move(cell_id, stone_value));
                        Value control_req = parity[p] * (alpha_ - OFFSET[p]);
                        if (adj_value + stones[0] < control_req) {
                            Value extra = control_req - adj_value;
                            for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                            if (j < n_stones) {
                                stone_value = parity[p] * stones[j];
                                add_child(pos, Move(cell_id, stone_value));
                            }
                        }       
                    }
                }
            }
            if (! pair) {
                for (j = 0; j < n_stones; j++) {
                    stone_value = parity[p] * stones[j];
                    add_child(pos, Move(cell_id, stone_value));
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * parity[p];
        add_child(pos, Move(best_stale, stone_value));
    }
}


MCTNode * MCTNode::expand(Position& pos) {
    expanded_ = true;
    get_children(pos);
    MCTNode * best_child = NULL;
    Real best_val = -1.0;
    size_t i;
    Real child_val;
    for (i = 0; i < children_.size(); i++) {
        if (children_[i]->fully_explored_)
            n_children_fully_explored_++;
        child_val = children_[i]->val_;
        if (pos.turn_ == BLUE)
            child_val = 1 - child_val;
        if (child_val > best_val) {
            best_child = children_[i];
            best_val = child_val;
        }
    }

    if (best_child->fully_explored_) {
        if ((n_children_fully_explored_ == children_.size()) ||
                (pos.turn_ == RED && best_child->val_ == 1.0) ||
                (pos.turn_ == BLUE && best_child->val_ == 0.0)) {
            fully_explored_ = true;
            return this;
        }
        best_val = -1.0;
        for (i = 0; i < children_.size(); i++) {
            if (! children_[i]->fully_explored_) {
                child_val = children_[i]->val_;
                if (pos.turn_ == BLUE)
                    child_val = 1 - child_val;
                if (child_val > best_val) {
                    best_child = children_[i];
                    best_val = child_val;
                }
            }
        }
    }

    return best_child ;
}


Move playout_moves[N_CELLS];

bool MCTNode::light_playout(Position& pos) {
    int move_count;
    bool result = true;
    bool result_found = false;
    //pos.take_snapshot();
    for (move_count = 0; pos.open_.len_ > 1; move_count++) {
        if (pos.is_winning(RED)) {
            result = true;
            result_found = true;
            break;
        }
        if (pos.is_winning(BLUE)) {
            result = false;
            result_found = true;
            break;
        } 
        //playout_moves[move_count] = pos.get_random_move();
        playout_moves[move_count] = pos.get_default_policy_move();
        pos.make_move(playout_moves[move_count]);
        //pos.make_move(pos.get_default_policy_move());
    }
    if (! result_found)
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    for (move_count--; move_count >= 0; move_count--)
        pos.unmake_move(playout_moves[move_count]);
    //pos.restore_snapshot();
    return result;
}


void MCTNode::ucb(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    MCTNode * search_root = this;
    while (search_root->expanded_) {
        search_root = search_root->select(pos);
#ifdef DEBUG_
        assert(search_root != NULL);
#endif
        pos.make_move(search_root->move_);
    }
    
    
    search_root = search_root->expand(pos);
    
    size_t result = 0;
    if (! search_root->fully_explored_) {
        pos.make_move(search_root->move_);
        result = search_root->light_playout(pos);
    }

    while (search_root != NULL) {
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (size_t i = 0; i < search_root->children_.size(); i++)
                if (search_root->children_[i]->val_ * pos.m_ > search_root->val_)
                    search_root->val_ = search_root->children_[i]->val_ * pos.m_;
            search_root->val_ *= pos.m_;
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
    for (size_t i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if ((! fully_explored_ || child->fully_explored_) &&
                child->val_ * pos.m_ > best_val) {
            best_move = child->move_;
            best_val = child->val_ * pos.m_;
            n_playouts = child->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %u playouts\n", n_playouts);
    return best_move;
}


Move MCTNode::get_most_played_move() {
    Move best_move;
    size_t most_visits = 0;
    for (size_t i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (child->n_playouts_ > most_visits) {
            best_move = child->move_;
            most_visits = child->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %d playouts\n", most_visits);
    return best_move;
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
        roots_[index]->fully_explored_ = false;
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
    long long move_time = std::min(time_left / 2LL, 1000000LL);
    long long alpha_time = std::min(move_time / 3LL, 330000LL);
    fprintf(stderr, "Move time = %lld\n", move_time);
    fprintf(stderr, "Alpha time = %lld\n", alpha_time);
    size_t n_playouts = 0;
    size_t i;
    MCTNode * root;
    
    if (pos.open_.len_ >= 22)
        return pos.get_expectation_maximising_move();
    
    while (time_now - start_micros < alpha_time) {
        pos.set_alpha(current_alpha_);
        root = get_or_make_root(current_alpha_, pos);
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
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
                if (root->fully_explored_) {
                    solved = true;
                    break;
                }
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
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






