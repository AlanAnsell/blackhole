struct MoveInfo {
    Move move_;
    int kill_;
    int adj_;

    MoveInfo() {}

    MoveInfo(const Move& move, int kill, int adj):
            move_(move), kill_(kill), adj_(adj) {}

    bool operator < (const MoveInfo& other) const {
        if (kill_ != other.kill_)
            return kill_ > other.kill_;
        if (move_.stone_value_ != other.move_.stone_value_)
            return abs(move_.stone_value_) > abs(other.move_.stone_value_);
        if (adj_ != other.adj_)
            return adj_ > other.adj_;
        return move_.cell_ < other.move_.cell_;
    }

};


bool Position::add_solver_move(const Move& move, std::vector<MoveInfo>& moves) {
    U32 cell_id = move.cell_;
    U64 mask = MASK(cell_id);
    Value stone_number = m_ * move.stone_value_;
    U32 op = 1 - turn_;
    int swing = adj_[cell_id].len_;
    int kill = 0;
    if (dead_[turn_] & mask)
        kill = -1;
    //else if (dead_[op] & mask)
    //    kill = 1;
    List& adj = adj_[cell_id];
    bool all_adj_isolated = true;
    Value largest_val_req = 0;
    for (U32 j = 0; j < adj.len_; j++) {
        U32 adj_id = adj.val_[j];
        U64 adj_mask = MASK(adj_id);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask)) {
            U32 n_adj = adj_[adj_id].len_;
            if (n_adj != 1)
                all_adj_isolated = false;
            Value adj_value = (power_[op][std::min(n_adj-1, n_stones_[op])] + value_[adj_id]) * m_;
            Value control_req = m_ * (alpha_ - OFFSET[turn_]);
            if (adj_value + stone_number >= control_req) {
                kill++;
                if (control_req - adj_value > largest_val_req)
                    largest_val_req = control_req - adj_value;
            }
        }
    }
    if (kill + n_dead_[turn_] > n_stones_[op]) {
        moves.push_back(MoveInfo(move, kill, swing));
        return true;
    }
    if (all_adj_isolated) {
        U32 stone_index = stone_loc_[turn_][stone_number];
        if (stone_index > n_stale_[op] && stones_[turn_][stone_index-1] >= largest_val_req)
            return false;
    }
    moves.push_back(MoveInfo(move, kill, swing));
    return false;
}


bool Position::get_solver_moves(std::vector<MoveInfo>& moves, bool ignore_killer) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    Value stone_value;
    U32 op = 1 - turn_;
    Value best_stale_val = 1000;
    U32 best_stale = NO_CELL;
    U32 i, j;
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        Value cell_value = value_[cell_id] * m_;
        U64 mask = MASK(cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            if (cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            if ((dead_[1-turn_] & mask) || (n_stones > n_dead_[1-turn_])) {
                bool done = false;
                if (adj_[cell_id].len_ == 1) {
                    U32 adj_id = adj_[cell_id].val_[0];
                    if (adj_[adj_id].len_ == 1) {
                        done = true;
                        Value adj_value = value_[adj_id] * m_;
                        if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                            Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                            if (adj_value + stones[0] < control_req) {
                                Value extra = control_req - adj_value;
                                for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                                if (j < n_stones) {
                                    stone_value = m_ * stones[j];
                                    if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                                        return true;
                                }
                            }
                            stone_value = m_ * stones[0];
                            if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                                return true;
                        }
                    } else {
                        if (we_control & mask)
                            done = true;
                    }
                }
                if (! done) {
                    for (j = n_stale_[op]; j < n_stones; j++) {
                        stone_value = m_ * stones[j];
                        if (add_solver_move(Move(cell_id, stone_value), moves) && ! ignore_killer)
                            return true;
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * m_;
        if (add_solver_move(Move(best_stale, stone_value), moves) && ! ignore_killer)
            return true;
    }
    return false;
}



U32 Position::solve(long long end_time, U32& counter, U32& hash_hits,
        HashTable& table, HashInfo& hash_info, U64 stone_masks[2]) {
    counter++;
    if (counter % 1000 == 0 && get_time() > end_time)
        return TIME_ELAPSED;
    if (is_winning(RED))
        return RED;
    if (is_winning(BLUE))
        return BLUE;
    U32 hash_result = table.find(stone_masks[0], stone_masks[1]);
    if (hash_result != NO_HASH_ENTRY) {
        hash_hits++;
        return hash_result;
    }
    U64 child_stone_masks[2];
    child_stone_masks[0] = stone_masks[0];
    child_stone_masks[1] = stone_masks[1];
    std::vector<MoveInfo> move_infos;
    if (get_solver_moves(move_infos, false))
        return turn_;
    U32 i;
    U32 op = 1 - turn_;
    std::sort(move_infos.begin(), move_infos.end());
    U32 optimal_result = op;
    for (i = 0; i < move_infos.size(); i++) {
        Move move = move_infos[i].move_;
        Value stone_number = move.stone_value_ * m_;
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[move.cell_];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move(move);
        if (result == TIME_ELAPSED)
            return TIME_ELAPSED;
        if (result == turn_) {
            optimal_result = turn_;
            break;
        }
    }
    table.add(stone_masks[0], stone_masks[1], optimal_result);
    return optimal_result;
}


std::pair<U32, Move> Position::get_optimal_move(long long end_time, U32& counter, U32& hash_hits,
        bool break_ties, HashTable& table) {
    table.init();
    HashInfo hash_info(*this);
    U64 stone_masks[2] = {0};
    U64 child_stone_masks[2];
    child_stone_masks[0] = stone_masks[0];
    child_stone_masks[1] = stone_masks[1];
    
    std::vector<MoveInfo> move_infos;
    get_solver_moves(move_infos, true);
    U32 i;
    std::sort(move_infos.begin(), move_infos.end());
  
    Move move; 
    std::vector<Move> winning_moves;
    for (i = 0; i < move_infos.size(); i++) {
        move = move_infos[i].move_;
        Value stone_number = move.stone_value_ * m_;
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[move.cell_];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move(move);
        if (result == TIME_ELAPSED)
            return std::pair<U32, Move>(TIME_ELAPSED, Move());
        if (result == turn_) {
            if (break_ties)
                winning_moves.push_back(move);
            else
                return std::pair<U32, Move>(turn_, move);
        }
    }
    
    if (winning_moves.size()) {
        fprintf(stderr, "%u optimal moves\n", (U32)winning_moves.size());
        Move best_move;
        Real best_val = -1000.0;
        for (i = 0; i < winning_moves.size(); i++) {
            move = winning_moves[i];
            make_move(move);
            Real val = calculate_expectation();
            unmake_move(move);
            val *= m_;
            if (val > best_val) {
                best_move = move;
                best_val = val;
            }
        }
        return std::pair<U32, Move>(turn_, best_move);
    }
    return std::pair<U32, Move>(1 - turn_, move_infos[0].move_);
}

