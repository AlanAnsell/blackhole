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
        val_ = 0.0; 
    }
    n_red_wins_ = 0;
    n_playouts_ = 0;
    child_moves_.clear();
    children_.clear();
    expanded_ = false;
    n_children_ = 0;
    n_child_moves_ = 0;
    n_children_fully_explored_ = 0;
    if (parent == NULL) {
        amaf_[0].init(pos.open_.len_, pos.n_stones_[RED]);
        amaf_[1].init(pos.open_.len_, pos.n_stones_[BLUE]);
    }
}


MCTNode * MCTNode::select(Position& pos, AMAFTable& amaf) {
    MCTNode * best_node = NULL;
    Real best_node_value = -1000.0, child_val;
    size_t i;
    for (i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            child_val = child->val_;
            if (pos.turn_ == BLUE)
                child_val = 1.0 - child_val;
            Real parent_playouts = (Real)n_playouts_;
            Real child_playouts = (Real)child->n_playouts_;
            Real ucb_value = child_val + UCB_C * sqrt(log(parent_playouts) / child_playouts);
            if (ucb_value > best_node_value) {
                best_node = child;
                best_node_value = ucb_value;
            }
        }
    }
    
    int best_move_index = -1;
    Real best_move_value = -1000.0;
    for (i = 0; i < n_child_moves_; i++) {
        Move& move = child_moves_[i];
        std::pair<size_t, size_t> indices = pos.get_cell_and_stone_indices(move);
        size_t n_amaf_playouts = amaf.get_n_playouts(indices.first, indices.second); 
        if (n_amaf_playouts > 0)
            child_val = amaf.get_value(indices.first, indices.second);
        else
            child_val = 0.5;
        Real parent_playouts = std::max(1.0, (Real)n_playouts_);
        child_val += UCB_C * sqrt(log(parent_playouts));
        if (child_val > best_move_value) {
            best_move_index = i;
            best_move_value = child_val;
        }
    }
    if (best_node == NULL || best_move_value > best_node_value) {
#ifdef DEBUG_
        assert(best_move_index != -1);
#endif
        n_child_moves_--;
        Move best_move = child_moves_[best_move_index];
        child_moves_[best_move_index] = child_moves_[n_child_moves_];
        return add_child(pos, best_move);
    }
#ifdef DEBUG_
    assert(best_node != NULL);
#endif
    return best_node;
}


MCTNode *  MCTNode::add_child(Position& pos, const Move& move) {
    pos.make_move(move);
    MCTNode * child = get_free();
    child->init(this, pos, move, alpha_);
    children_.push_back(child);
    pos.unmake_move(move);
    return child;
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
        int64 mask = (1LL << (int64)cell_id);
        if ((pos.stale_[RED] & mask) || (pos.stale_[BLUE] & mask)) {
            if (cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            if ((pos.dead_[1-p] & mask) || (n_stones > pos.n_dead_[1-p])) {
                size_t j;
                bool pair = false;
                if (pos.adj_[cell_id].len_ == 1) {
                    size_t adj_id = pos.adj_[cell_id].val_[0];
                    if (pos.adj_[adj_id].len_ == 1) {
                        pair = true;
                        Value adj_value = pos.value_[adj_id] * parity[p];
                        if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                            stone_value = parity[p] * stones[0];
                            child_moves_.push_back(Move(cell_id, stone_value));
                            Value control_req = parity[p] * (alpha_ - OFFSET[p]);
                            if (adj_value + stones[0] < control_req) {
                                Value extra = control_req - adj_value;
                                for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                                if (j < n_stones) {
                                    stone_value = parity[p] * stones[j];
                                    child_moves_.push_back(Move(cell_id, stone_value));
                                }
                            }       
                        }
                    }
                }
                if (! pair) {
                    for (j = pos.n_stale_[1 - p]; j < n_stones; j++) {
                        stone_value = parity[p] * stones[j];
                        child_moves_.push_back(Move(cell_id, stone_value));
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * parity[p];
        child_moves_.push_back(Move(best_stale, stone_value));
    }
    n_child_moves_ = n_children_ = child_moves_.size();
}


MCTNode * MCTNode::expand(Position& pos, AMAFTable& amaf) {
    expanded_ = true;
    get_children(pos);
    return select(pos, amaf);
}


std::pair<size_t, size_t> simulation[N_CELLS];

bool MCTNode::light_playout(Position& pos, size_t& move_count) {
    bool result = true;
    bool result_found = false;
    pos.take_snapshot();
    while (pos.open_.len_ > 1) {
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
        Move move = pos.get_default_policy_move();
        simulation[move_count] = pos.get_cell_and_stone_indices(move);
        pos.make_move(move);
        move_count++;
    }
    if (! result_found)
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    pos.restore_snapshot();
    return result;
}


void MCTNode::ucb(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    MCTNode * search_root = this;
    size_t orig_turn = pos.turn_;
    size_t move_count = 0, i;
    while (search_root->expanded_) {
        search_root = search_root->select(pos, amaf_[pos.turn_]);
#ifdef DEBUG_
        assert(search_root != NULL);
#endif
        simulation[move_count] = pos.get_cell_and_stone_indices(search_root->move_);
        pos.make_move(search_root->move_);
        move_count++;
    }
    
    if (! search_root->fully_explored_) { 
        search_root = search_root->expand(pos, amaf_[pos.turn_]);
        simulation[move_count] = pos.get_cell_and_stone_indices(search_root->move_);
        move_count++;
        pos.make_move(search_root->move_);
    }
    size_t result;
    if (search_root->fully_explored_)
        result = (search_root->val_ == 1.0);
    else
        result = search_root->light_playout(pos, move_count);

    for (i = 0; i < move_count; i++) {
        size_t turn = (orig_turn + i) % 2;
        amaf_[turn].update(simulation[i].first, simulation[i].second, result != turn);
    }

    while (search_root != NULL) {
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (i = 0; i < search_root->children_.size(); i++)
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
                parent->n_children_fully_explored_++;
                if ((parent->n_children_fully_explored_ == parent->n_children_) ||
                        (pos.turn_ == RED && search_root->val_ == 1.0) ||
                        (pos.turn_ == BLUE && search_root->val_ == 0.0)) {
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
    
    if (pos.open_.len_ >= 24)
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






