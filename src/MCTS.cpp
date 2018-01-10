#include "MCTS.h"


void MCTNode::dispose() {
    amaf_.dispose();
    for (U32 i = 0; i < children_.size(); i++) {
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
    
    expanded_ = false;
    all_children_generated_ = false;
    
    n_children_fully_explored_ = 0;
    children_.clear();
    
    solved_ = false;
    solve_attempted_ = false;
    solver_positions_ = 0;
    solver_hash_hits_ = 0;
    solver_time_ = 0;
    amaf_.init(pos.open_.len_, pos.n_stones_[pos.turn_]);
}


MCTNode * MCTNode::select(Position& pos) {
    MCTNode * best_node = NULL;
    Real best_node_value = -1000.0, child_val;
    U32 i;
    Real log_pp = log((Real)n_playouts_);
    for (i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            std::pair<U32, U32> indices = pos.get_cell_and_stone_indices(child->move_);
            U32 n_amaf = amaf_.get_n_playouts(indices.first, indices.second);
            Real beta = (Real)n_amaf / (n_amaf + child->n_playouts_ + (RAVE_C * n_amaf) * child->n_playouts_);
            Real mc_val = child->val_;
            if (pos.turn_ == BLUE)
                mc_val = 1.0 - val_;
            child_val = (1.0 - beta) * mc_val + beta * amaf_.get_value(indices.first, indices.second) + 
                    UCB_C * sqrt(log_pp / child->n_playouts_);
            //Real child_playouts = (Real)child->n_playouts_;
            //Real ucb_value = child_val + UCB_C * sqrt(log_pp / child_playouts);
            if (child_val > best_node_value) {
                best_node = child;
                best_node_value = child_val;
            }
        }
    }
    
    if (! all_children_generated_) {
        U32 cell_index, stone_index;
        generate_move(cell_index, stone_index, pos);
        // since we don't know in advance whether all children have been generated,
        // we might discover during selection that this node is already fully explored.
        // This shouldn't happen during expansion however
        if (cell_index != NO_CELL) {
            U32 n_amaf_playouts = amaf_.get_n_playouts(cell_index, stone_index); 
            if (n_amaf_playouts > 0)
                child_val = amaf_.get_value(cell_index, stone_index);
            else
                child_val = 0.5;
            child_val += UCB_C * sqrt(log_pp);
            if (child_val > best_node_value) {
                Move move(pos.open_.val_[cell_index],
                        pos.m_ * pos.stones_[pos.turn_][stone_index]);
                best_node = add_child(pos, move);
            }
        }
    }

    if (best_node == NULL) {
/*#ifdef DEBUG_
        fprintf(stderr, "select returning this node\n");
        fflush(stderr);
        assert(is_now_fully_explored());
#endif*/
        fully_explored_ = true;
        return this;
    }
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


bool MCTNode::is_now_fully_explored() {
    return all_children_generated_ &&
            n_children_fully_explored_ == children_.size();
}


void MCTNode::generate_move(U32& cell_index, U32& stone_index, Position& pos) {
    cell_index = NO_CELL;
    do {
        amaf_.get_best(cell_index, stone_index);
    } while (cell_index != NO_CELL && ! pos.is_reasonable_move(cell_index, stone_index));
    
    if (cell_index == NO_CELL) {
        pos.get_untried_move(cell_index, stone_index, amaf_);
        if (cell_index == NO_CELL) {
            all_children_generated_ = true;
            return;
        }
    }
}
    

MCTNode * MCTNode::expand(Position& pos) {
    expanded_ = true;
    return select(pos);
}


Move simulation[N_CELLS];

bool MCTNode::playout(Position& pos, U32& move_count) {
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
        simulation[move_count] = pos.get_default_policy_move();
        pos.make_move(simulation[move_count]);
        move_count++;
    }
    if (! result_found)
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    pos.restore_snapshot();
    return result;
}


void MCTNode::simulate(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
    //fprintf(stderr, "Simulating\n");
    //fflush(stderr);
    MCTNode * search_root = this;
    MCTNode * new_search_root;
    U32 move_count = 0, depth = 0, i;
    while (search_root->expanded_) {
        new_search_root = search_root->select(pos);
#ifdef DEBUG_
        assert(new_search_root != NULL);
#endif
        if (new_search_root == search_root)
            break;
        search_root = new_search_root;
        simulation[move_count] = search_root->move_;
        pos.make_move(search_root->move_);
        move_count++;
        depth++;
    }
    
    /*if (! search_root->fully_explored_) { 
        new_search_root = search_root->expand(pos);
#ifdef DEBUG_
        assert(new_search_root != NULL && new_search_root != search_root);
#endif
        search_root = new_search_root;
        simulation[move_count] = search_root->move_;
        move_count++;
        depth++;
        pos.make_move(search_root->move_);
    }*/
    U32 result;
    if (search_root->fully_explored_)
        result = (search_root->val_ == 1.0);
    else {
        result = search_root->playout(pos, move_count);
        search_root->expanded_ = true;
        search_root = search_root->add_child(pos, simulation[depth]);
        pos.make_move(simulation[depth]);
        depth++;
    }

    while (search_root != NULL) {
        U32 last_move = std::min(move_count, depth + 11);
        for (i = depth; i < last_move; i += 2) {
            std::pair<int, int> indices = pos.get_cell_and_stone_indices(simulation[i]);
            search_root->amaf_.update(indices.first, indices.second, result != pos.turn_);
        }
        
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
                if ((parent->is_now_fully_explored()) ||
                        (pos.turn_ == RED && search_root->val_ == 1.0) ||
                        (pos.turn_ == BLUE && search_root->val_ == 0.0)) {
                    parent->fully_explored_ = true;
                }
            }
        }
        search_root = parent;
        depth--;
    }
}


Move MCTNode::get_highest_value_move(Position& pos) {
    if (solved_)
        return solution_;
    Move best_move;
    Real best_val = -1000.0;
    U32 n_playouts = 0;
    for (U32 i = 0; i < children_.size(); i++) {
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
    U32 most_visits = 0;
    for (U32 i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (child->n_playouts_ > most_visits) {
            best_move = child->move_;
            most_visits = child->n_playouts_;
        }
    }
    fprintf(stderr, "Selected move has %d playouts\n", most_visits);
    return best_move;
}


bool MCTNode::attempt_solve(Position& pos, HashTable& table, long long allowed_time, bool break_ties) {
    solve_attempted_ = true;
    long long start_time = get_time();
    long long end_time = start_time + allowed_time;
    std::pair<U32, Move> result = pos.get_optimal_move(end_time, solver_positions_, solver_hash_hits_, break_ties, table);
    solver_time_ = get_time() - start_time;
    if (result.first != TIME_ELAPSED) {
        if (result.first == RED)
            val_ = 1.0;
        else
            val_ = 0.0;
        solved_ = true;
        fully_explored_ = true;
        solution_ = result.second;
        return true;
    }
    return false;
}


void MCTNode::print(FILE * f) {
    char playout_status[20];
    char solve_status[100];
    if (solved_)
        sprintf(playout_status, "*");
    else
        sprintf(playout_status, "%u", n_playouts_);
    if (solve_attempted_) {
        Real solver_time = (Real)solver_time_ / 1.0e6;
        if (solved_)
            sprintf(solve_status, "(solver SUCCEEDED: %u nodes, %u hash hits, %.4lfs)",
                    solver_positions_, solver_hash_hits_, solver_time);
        else
            sprintf(solve_status, "(solver FAILED: %u nodes, %u hash hits, %.4lfs)",
                    solver_positions_, solver_hash_hits_, solver_time);
    } else {
        solve_status[0] = 0;
    }
    fprintf(f, "Alpha = %d (%s): %.4f %s\n", alpha_, playout_status, val_, solve_status);
}



MCTNode node_store[N_MCT_NODES];
MCTNode * free_list[N_MCT_NODES];
U32 n_free;


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
    for (U32 i = 0; i < N_MCT_NODES; i++)
        free_list[i] = &node_store[i];
    n_free = N_MCT_NODES;
}


MCTSearch::MCTSearch(const Position& pos, Value alpha) {
    current_alpha_ = alpha;
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++)
        roots_[i] = NULL;
}


MCTSearch::~MCTSearch() {
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++) {
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
                test_root->print(stderr);
                if (test_root->val_ >= CHOOSE_TARGET_THRESH[RED] || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    } else {
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL) {
                test_root->print(stderr);
                if (test_root->val_ < CHOOSE_TARGET_THRESH[BLUE] || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    }
    return root;
}


Move MCTSearch::get_best_move(Position& pos) {
    if (pos.open_.len_ >= 32) {
        std::pair<Real, Move> search_result = pos.get_best_move();
        //fprintf(stderr, "Eval = %.5lf\n", search_result.first);
        return search_result.second;
    }
    
    Real time_left_r = (Real) time_left;
    Real use_proportion;
    int moves_left = ((int)pos.open_.len_ - 16) / 2;
    if (moves_left < 1)
        use_proportion = 0.5;
    else {
        fprintf(stderr, "Expecting to make %d more moves before position solved\n", moves_left);
        use_proportion = 1.0 - pow(((Real)time_limit / 10.0)  / time_left_r, 1.0 / (Real)moves_left);
    }
    long long move_time = (long long)(use_proportion * time_left_r);
    long long start_micros = get_time();
    long long time_now = start_micros;
    long long alpha_time = move_time / 3LL;
    long long solver_time = alpha_time / 2LL;
    fprintf(stderr, "Move time = %lld\n", move_time);
    fprintf(stderr, "Alpha time = %lld\n", alpha_time);
    U32 n_playouts = 0;
    U32 i;
    MCTNode * root;
    U32 solver_start = 19;
    
    while (time_now - start_micros < alpha_time) {
        pos.set_alpha(current_alpha_);
        root = get_or_make_root(current_alpha_, pos);
        if (pos.open_.len_ <= solver_start && ! root->solve_attempted_ && time_now - start_micros < solver_time)
            root->attempt_solve(pos, table_, solver_time - time_now + start_micros, false);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
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
            if (root->val_ >= TARGET_INCREMENT_THRESH[pos.turn_] && current_alpha_ < MAX_RESULT)
                current_alpha_++;
            if (root->val_ <= TARGET_DECREMENT_THRESH[pos.turn_] && current_alpha_ > MIN_RESULT)
                current_alpha_--;
        }
        time_now = get_time();
    }

    root = select_alpha(pos);
    
    Value original_alpha = current_alpha_;
    fprintf(stderr, "Alpha = %d\n", current_alpha_);
    pos.set_alpha(current_alpha_);
    bool done = false;
    while (time_now - start_micros < move_time) {
        //if (pos.open_.len_ <= solver_start && ! root->solve_attempted_)
        //    root->attempt_solve(pos, table_);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
                if (root->fully_explored_) {
                    done = true;
                    break;
                }
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
                pos.set_alpha(current_alpha_);
                root = get_or_make_root(current_alpha_, pos);
                if (root->fully_explored_) {
                    done = true;
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
    if (done)
        root = select_alpha(pos);
    pos.set_alpha(current_alpha_);

    if (root->solved_) {
        long long time_remaining = move_time - (time_now - start_micros);
        if (time_remaining > 0) {
            if (root->attempt_solve(pos, table_, time_remaining, true))
                fprintf(stderr, "Chose most challenging optimal move\n");
        }
    }

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(time_now - start_micros) / 1e6);
        
    if (root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return root->get_highest_value_move(pos);
    }

    return root->get_most_played_move();
}






