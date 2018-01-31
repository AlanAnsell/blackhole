#include "MCTS.h"


void MCTNode::dispose() {
    for (U32 i = 0; i < children_.size(); i++) {
        children_[i]->dispose();
        put_free(children_[i]);
    }
}


void MCTNode::init(MCTNode * parent, const Position& pos, Move move, Value alpha, AMAFTable * my_amaf) {
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
    untried_children_.clear();

    //get_untried_failed_ = false;
    untried_stone_index_ = pos.n_stones_[pos.turn_] - 1;
    pos.get_validity_mask(valid_, duo_);
    
    n_children_fully_explored_ = 0;
    children_.clear();
    
    solved_ = false;
    solve_attempted_ = false;
    solver_positions_ = 0;
    solver_hash_hits_ = 0;
    solver_time_ = 0;

    tried_ = std::vector<U32>(pos.n_stones_[pos.turn_], 0);
    if (my_amaf) {
        amaf_node_ = false;
        my_amaf_ = my_amaf;
    } else {
        amaf_node_ = true;
        amaf_[RED].init(pos.open_.len_, pos.n_stones_[RED]);
        amaf_[BLUE].init(pos.open_.len_, pos.n_stones_[BLUE]);
        my_amaf_ = & amaf_[pos.turn_];
    }
}


MCTNode * MCTNode::select(Position& pos) {
    MCTNode * best_node = NULL;
    Real best_node_value = -1000.0, child_val;
    U32 i;
    Real sqrt_log_pp = SQRT_LOG[n_playouts_];
    for (i = 0; i < children_.size(); i++) {
        MCTNode * child = children_[i];
        if (! child->fully_explored_) {
            std::pair<U32, U32> indices = pos.get_cell_and_stone_indices(child->move_);
            U32 n_amaf = my_amaf_->get_n_playouts(indices.first, indices.second);
#ifdef DEBUG_
            assert(n_amaf > 0);
#endif
            Real beta = (Real)n_amaf / (n_amaf + child->n_playouts_ + (RAVE_C * n_amaf) * child->n_playouts_);
            Real mc_val = child->val_;
            if (pos.turn_ == BLUE)
                mc_val = 1.0 - val_;
            child_val = (1.0 - beta) * mc_val + beta * my_amaf_->get_value(indices.first, indices.second) + 
                    UCB_C * sqrt_log_pp / SQRT[child->n_playouts_];
            //Real child_playouts = (Real)child->n_playouts_;
            //Real ucb_value = child_val + UCB_C * sqrt(log_pp / child_playouts);
            if (child_val > best_node_value) {
                best_node = child;
                best_node_value = child_val;
            }
        }
    }
    
    if (! all_children_generated_) {
        U32 cell_index = NO_CELL, stone_index;
        while (true) {
            my_amaf_->get_best(cell_index, stone_index);
            if (cell_index == NO_CELL || pos.is_reasonable_move(cell_index, stone_index, valid_, duo_))
                break;
            my_amaf_->play(cell_index, stone_index);
        }
   
        if (cell_index == NO_CELL) 
            generate_move(cell_index, stone_index, pos);
        // since we don't know in advance whether all children have been generated,
        // we might discover during selection that this node is already fully explored.
        // This shouldn't happen during expansion however
        if (cell_index != NO_CELL) {
            U32 n_amaf_playouts = my_amaf_->get_n_playouts(cell_index, stone_index); 
            if (n_amaf_playouts > 0)
                child_val = my_amaf_->get_value(cell_index, stone_index);
            else {
                Real node_val;
                if (pos.turn_ == RED)
                    node_val = val_;
                else
                    node_val = 1.0 - val_;
                child_val = (2.0 * node_val + 0.5) / 3.0;
            }
            //if (pos.turn_ == BLUE)
            //    child_val = 1.0 - child_val;
            child_val += UCB_C * sqrt_log_pp;
            if (child_val > best_node_value) {
                Move move = CREATE_MOVE(pos.open_.val_[cell_index],
                                        pos.stones_[pos.turn_][stone_index]);
#ifdef DEBUG_
                if (! pos.is_legal(move)) { 
                    pos.print(stderr);
                    char move_str[20];
                    move_to_str(move, move_str);
                    fprintf(stderr, "%s\n", move_str);
                    fflush(stderr);
                    assert(false);
                }
#endif
                best_node = add_child(pos, move);
                my_amaf_->play(cell_index, stone_index);
            }
        }
    }

    if (best_node == NULL) {
#ifdef DEBUG_
        //fprintf(stderr, "select returning null node\n");
        //fflush(stderr);
        assert(children_.size() > 0);
        assert(is_now_fully_explored());
#endif
        fully_explored_ = true;
        return this;
    }
    return best_node;
}


MCTNode *  MCTNode::add_child(Position& pos, Move move) {
    pos.make_move(move, true);
    MCTNode * child = get_free();
    AMAFTable * child_amaf;
    if (amaf_node_)
        child_amaf = & amaf_[pos.turn_];
    else
        child_amaf = NULL; 
    child->init(this, pos, move, alpha_, child_amaf);
    children_.push_back(child);
    pos.unmake_move();
    return child;
}


bool MCTNode::is_now_fully_explored() {
    return all_children_generated_ &&
            n_children_fully_explored_ == children_.size();
}


bool MCTNode::is_playable() {
    return fully_explored_ || children_.size() > 0;
}


void MCTNode::generate_move(U32& cell_index, U32& stone_index, Position& pos) {
    while (true) {
        if (untried_children_.size() == 0) {
            pos.generate_untried_moves(untried_stone_index_, valid_, duo_, untried_children_);
            if (untried_children_.size() == 0) {
                cell_index = NO_CELL;
                all_children_generated_ = true;
                return;
            }
        }
        Move move = untried_children_.back();
        untried_children_.pop_back();
        std::pair<U32, U32> indices = pos.get_cell_and_stone_indices(move);
        cell_index = indices.first;
        stone_index = indices.second;
        if (! my_amaf_->is_played(cell_index, stone_index))
            return;
    }
}
    

Move simulation[N_CELLS];

bool MCTNode::playout(Position& pos, U32& move_count) {
#ifdef DEBUG_
    if (pos.is_winning(RED) || pos.is_winning(BLUE)) {
        pos.print(stderr);
        for (U32 m = 0; m < move_count; m++) {
            char move_str[20];
            move_to_str(simulation[m], move_str);
            fprintf(stderr, "%s\n", move_str);
        }
        fprintf(stderr, "alpha = %d\n", alpha_);
        fprintf(stderr, "n_dead_[RED] = %u, mask = %llx\n", pos.n_dead_[RED], pos.dead_[RED]);
        fprintf(stderr, "n_dead_[BLUE] = %u, mask = %llx\n", pos.n_dead_[BLUE], pos.dead_[BLUE]);
        fflush(stderr);
    }
    assert(! pos.is_winning(RED) && ! pos.is_winning(BLUE));
#endif
    bool result = true;
    bool result_found = false;
    U32 n_moves_made = pos.n_moves_made_;
    pos.save_history();
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
        U64 valid, duo;
        pos.get_validity_mask(valid, duo);
        simulation[move_count] = pos.get_default_policy_move(valid, duo);
        pos.make_move(simulation[move_count], false);
        move_count++;
    }
    if (! result_found) {
#ifdef DEBUG_
        //fprintf(stderr, "Got to end of game\n");
#endif
        result = (pos.value_[pos.open_.val_[0]] >= alpha_);
    }
    pos.rewind_to(n_moves_made);
    return result;
}


void MCTNode::simulate(Position& pos) {
#ifdef DEBUG_
    assert(! fully_explored_);
#endif
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
#ifdef DEBUG_
        if (! pos.is_legal(search_root->move_)) { 
            pos.print(stderr);
            char move_str[20];
            move_to_str(search_root->move_, move_str);
            fprintf(stderr, "%s\n", move_str);
            fflush(stderr);
            assert(false);
        }
#endif
        simulation[move_count] = search_root->move_;
        pos.make_move(search_root->move_, true);
        move_count++;
        depth++;
    }
    
    U32 result = NO_RESULT;
    if (! search_root->fully_explored_) {
        result = search_root->playout(pos, move_count);
#ifdef DEBUG_
        if (! pos.is_legal(simulation[depth])) { 
            pos.print(stderr);
            char move_str[20];
            move_to_str(simulation[depth], move_str);
            fprintf(stderr, "%s\n", move_str);
            fflush(stderr);
            assert(false);
        }
#endif
        search_root->expanded_ = true;
        search_root = search_root->add_child(pos, simulation[depth]);
        pos.make_move(simulation[depth], true);
        depth++;
    } else {
        result = (search_root->val_ == 1.0);
    }

    while (search_root != NULL) {
        U32 last_move = std::min(move_count, depth + AMAF_HORIZON);
        for (i = depth; i < last_move; i += 2) {
            std::pair<int, int> indices = pos.get_cell_and_stone_indices(simulation[i]);
            search_root->my_amaf_->update(indices.first, indices.second, result != pos.turn_);
        }
        
        MCTNode * parent = search_root->parent_;
        if (search_root->fully_explored_ && search_root->expanded_) {
            search_root->val_ = -1000.0;
            for (i = 0; i < search_root->children_.size(); i++)
                if (search_root->children_[i]->val_ * pos.m_ > search_root->val_)
                    search_root->val_ = search_root->children_[i]->val_ * pos.m_;
            search_root->val_ *= pos.m_;
            result = (search_root->val_ == 1.0);
        }
        if (! search_root->fully_explored_) {
#ifdef DEBUG_
            assert(result != NO_RESULT);
#endif
            search_root->n_red_wins_ += result;
            search_root->n_playouts_++;
            search_root->val_ = Real(search_root->n_red_wins_) / Real(search_root->n_playouts_);
        }
        if (parent != NULL) {
            pos.unmake_move();
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
    if (fully_explored_ && children_.size() == 0)
        return pos.get_best_winning_move();
    Move best_move = 0;
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
    Move best_move = 0;
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


void MCTNode::print_pv(FILE * f) {
    fprintf(f, "PV: ");
    MCTNode * curr = this;
    char move_str[20];
    while (curr->children_.size() > 0) {
        U32 most_visits = 0;
        MCTNode * best_child = NULL;
        for (U32 i = 0; i < curr->children_.size(); i++) {
            MCTNode * child = curr->children_[i];
            if (best_child == NULL || child->n_playouts_ > most_visits) {
                best_child = child;
                most_visits = child->n_playouts_;
            }
        }
        if (best_child == NULL) {
#ifdef DEBUG_
            fprintf(stderr, "PV child was NULL\n");
#endif
            break;
        }
        if (curr == this)
            fprintf(stderr, "(%.4lf)", best_child->val_);
        move_to_str(best_child->move_, move_str);
        fprintf(f, " %s (%u)", move_str, best_child->n_playouts_);
        curr = best_child;
    }
    fprintf(f, "\n");
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


void MCTNode::print_most_played_moves(U32 n) const {
    std::vector<std::pair<U32, U32> > moves;
    U32 i;
    MCTNode * child;
    for (i = 0; i < children_.size(); i++) {
        child = children_[i];
        moves.push_back(std::make_pair(child->n_playouts_, i));
    }
    sort(moves.begin(), moves.end(), std::greater<std::pair<U32, U32>>());
    fprintf(stderr, "Alpha = %d:\n", alpha_);
    for (i = 0; i < moves.size() && i < n; i++) {
        child = children_[moves[i].second];
        char move_str[10];
        move_to_str(child->move_, move_str);
        fprintf(stderr, "%s (%.4lf): %u\n", move_str, child->val_, child->n_playouts_);
    }
    if (children_.size() > 0) {
        fprintf(stderr, "Printing best child...\n");
        children_[moves[0].second]->print_most_played_moves(n);
    }
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


MCTSearch::MCTSearch(Value alpha) {
    fprintf(stderr, "Creating MCTSearch\n");
    fflush(stderr);
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


bool MCTSearch::no_playable_move() {
    for (U32 i = 0; i < 2 * MAX_RESULT + 1; i++)
        if (roots_[i] != NULL && roots_[i]->is_playable())
            return false;
    return true;
}


MCTNode * MCTSearch::get_or_make_root(Value alpha, const Position& pos) {
    Value index = alpha + MAX_RESULT;
#ifdef DEBUG_
    assert(index >= 0 && index <= 2 * MAX_RESULT);
#endif
    if (roots_[index] == NULL) {
        roots_[index] = get_free();
        roots_[index]->init(NULL, pos, Move(), alpha, NULL);
    }
    return roots_[index];
}

MCTNode * MCTSearch::select_alpha(const Position& pos, Real red_choose_target_thresh) {
    MCTNode * test_root;
    MCTNode * root = NULL;
    Value alpha;
    Real choose_target_thresh;
    if (pos.turn_ == RED) {
        choose_target_thresh = red_choose_target_thresh;
        for (alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL && test_root->is_playable()) {
                if (test_root->val_ >= choose_target_thresh || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    } else {
        choose_target_thresh = 1.0 - red_choose_target_thresh;
        for (alpha = MAX_RESULT; alpha >= MIN_RESULT; alpha--) {
            test_root = roots_[alpha + MAX_RESULT];
            if (test_root != NULL && test_root->is_playable()) {
                if (test_root->val_ < choose_target_thresh || root == NULL) {
                    current_alpha_ = alpha;
                    root = test_root;
                }
            }
        }
    }
    return root;
}


void MCTSearch::display_roots() {
    for (Value alpha = MIN_RESULT; alpha <= MAX_RESULT; alpha++) {
        MCTNode * root = roots_[alpha + MAX_RESULT];
        if (root != NULL)
            root->print(stderr);
    }
}


Move MCTSearch::get_best_move(Position& pos) {
    Real time_left_r = (Real) time_left;
    Real use_proportion;
    int moves_left = ((int)pos.open_.len_ - 16) / 2;
    if (moves_left < 1)
        use_proportion = 0.5;
    else {
        fprintf(stderr, "Expecting to make %d more moves before position solved\n", moves_left);
#ifdef UNIFORM_TM
        double use = (time_left_r - time_limit / 10.0) / (Real)moves_left;
        use_proportion = use / time_left_r;
#else
        use_proportion = 1.0 - pow(((Real)time_limit / 10.0)  / time_left_r, 1.0 / (Real)moves_left);
#endif
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
    
    init_amaf_free_list();
        
    while (time_now - start_micros < alpha_time || no_playable_move()) {
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
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
            if (roots_[current_alpha_ + MAX_RESULT] != NULL &&
                    roots_[current_alpha_ + MAX_RESULT]->fully_explored_)
                break;
        } else {
            Real target_increment_thresh, target_decrement_thresh;
            //if (pos.open_.len_ >= 22) {
                target_increment_thresh = TARGET_INCREMENT_THRESH[pos.turn_];
                target_decrement_thresh = TARGET_DECREMENT_THRESH[pos.turn_];
            /*} else {
                target_increment_thresh = 0.5;
                target_decrement_thresh = 0.5;
            }*/
            if (root->val_ >= target_increment_thresh && current_alpha_ < MAX_RESULT)
                current_alpha_++;
            if (root->val_ <= target_decrement_thresh && current_alpha_ > MIN_RESULT)
                current_alpha_--;
        }
        time_now = get_time();
    }

    root = select_alpha(pos, CHOOSE_TARGET_THRESH[RED]);
    display_roots();
     
    Value original_alpha = current_alpha_;
    fprintf(stderr, "Alpha = %d\n", current_alpha_);
    pos.set_alpha(current_alpha_);
    bool done = false;
    while (time_now - start_micros < move_time) {
        root = get_or_make_root(current_alpha_, pos);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
            pos.set_alpha(current_alpha_);
            if (roots_[current_alpha_ + MAX_RESULT] != NULL &&
                    roots_[current_alpha_ + MAX_RESULT]->fully_explored_) {
                done = true;
                break;
            }
        }
        time_now = get_time();
    }

    current_alpha_ = original_alpha;
    root = roots_[current_alpha_ + MAX_RESULT];
    if (done)
        root = select_alpha(pos, CHOOSE_TARGET_THRESH[RED]);
    display_roots();
    pos.set_alpha(current_alpha_);

    if (root->solved_) {
        long long time_remaining = move_time - (time_now - start_micros);
        if (time_remaining > 0) {
            if (root->attempt_solve(pos, table_, time_remaining, true))
                fprintf(stderr, "Chose most challenging optimal move\n");
        }
    }

    fprintf(stderr, "%u playouts in %.3lfs\n", n_playouts, (double)(get_time() - time_started) / 1e6);
        
    if (root->fully_explored_) {
        fprintf(stderr, "Search tree fully explored\n");
        return root->get_highest_value_move(pos);
    }
    
    root->print_pv(stderr);
    return root->get_most_played_move();
}


void MCTSearch::analyse(Position& pos) {
    U32 n_playouts = 0;
    U32 prev_print = 0;
    U32 i;
    MCTNode * root;
    init_amaf_free_list();
        
    while (true) {
        if (n_playouts - prev_print >= 10000) {
            fprintf(stderr, "*******************************************************\n");
            display_roots();
            fflush(stderr);
            root = select_alpha(pos, 0.5);
            root->print_most_played_moves(1000);
            root->print_pv(stderr);
            prev_print = n_playouts;
        }
        pos.set_alpha(current_alpha_);
        root = get_or_make_root(current_alpha_, pos);
        for (i = 0; i < 100 && ! root->fully_explored_; i++) {
            root->simulate(pos);
            n_playouts++;
        } 
        if (root->fully_explored_) {
            if (root->val_ == 0.0 && current_alpha_ > MIN_RESULT) {
                current_alpha_--;
            } else if (root->val_ == 1.0 && current_alpha_ < MAX_RESULT) {
                current_alpha_++;
            } else {
                fprintf(stderr, "MIN/MAX RESULT achieved\n");
                break;
            }
            if (roots_[current_alpha_ + MAX_RESULT] != NULL &&
                    roots_[current_alpha_ + MAX_RESULT]->fully_explored_)
                break;
        } else {
            if (root->val_ >= 0.5 && current_alpha_ < MAX_RESULT)
                current_alpha_++;
            if (root->val_ < 0.5 && current_alpha_ > MIN_RESULT)
                current_alpha_--;
        }
    }

    display_roots();
    root = select_alpha(pos, 0.5);
    root->print_pv(stderr);
}

