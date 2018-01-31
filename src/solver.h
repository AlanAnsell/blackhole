struct MoveInfo {
    Move move_;
    int kill_;
    int adj_;

    MoveInfo() {}

    MoveInfo(const Move& move, int kill, int adj):
            move_(move), kill_(kill), adj_(adj) {}

    inline bool operator < (const MoveInfo& other) const {
        if (kill_ != other.kill_)
            return kill_ > other.kill_;
        if (GET_STONE_NUMBER(move_) != GET_STONE_NUMBER(other.move_))
            return GET_STONE_NUMBER(move_) > GET_STONE_NUMBER(other.move_);
        if (adj_ != other.adj_)
            return adj_ > other.adj_;
        return GET_CELL(move_) < GET_CELL(other.move_);
    }

};


bool Position::add_solver_move(Move move, MoveInfo * moves, U32& N) {
    U32 cell_id = GET_CELL(move);
    U64 mask = MASK(cell_id);
    Value stone_number = GET_STONE_NUMBER(move);
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
        moves[N++] = MoveInfo(move, kill, swing);
        return true;
    }
    if (all_adj_isolated) {
        U32 stone_index = stone_loc_[turn_][stone_number];
        if (stone_index > n_stale_[op] && stones_[turn_][stone_index-1] >= largest_val_req)
            return false;
    }
    moves[N++] = MoveInfo(move, kill, swing);
    return false;
}


bool Position::get_solver_moves(MoveInfo * moves, U32& N, bool ignore_killer) {
    N = 0;
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    U32 op = 1 - turn_;
    U32 i, j;
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        //Value cell_value = value_[cell_id] * m_;
        //U64 mask = MASK(cell_id);
        bool duo = false;
        if (is_valid(cell_id, duo)) {
            if (cell_id == worst_stale_)
                add_solver_move(CREATE_MOVE(cell_id, stones[0]), moves, N);
            else if (duo) {
                U32 adj_id = get_non_stale_adj(adj_[cell_id], 0);
                Value adj_value = value_[adj_id] * m_;
                Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                if (adj_value + stones[0] < control_req) {
                    Value extra = control_req - adj_value;
                    for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                    if (j < n_stones) {
                        if (add_solver_move(CREATE_MOVE(cell_id, stones[j]), moves, N) && ! ignore_killer)
                            return true;
                    }
                }
                if (add_solver_move(CREATE_MOVE(cell_id, stones[0]), moves, N) && ! ignore_killer)
                    return true;
            } else {
                for (j = n_stale_[op]; j < n_stones; j++) {
                    if (add_solver_move(CREATE_MOVE(cell_id, stones[j]), moves, N) && ! ignore_killer)
                        return true;
                }
            }
        }
    }
    return false;
}


MoveInfo MOVE_INFO_STORE[N_CELLS][N_CELLS * N_STONES];

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
    
    MoveInfo * move_infos = MOVE_INFO_STORE[n_moves_made_];
    U32 N = 0;
    if (get_solver_moves(move_infos, N, false))
        return turn_;
    std::sort(move_infos, move_infos + N);
    
    U32 i;
    U32 op = 1 - turn_;
    U32 optimal_result = op;
    for (i = 0; i < N; i++) {
        Move move = move_infos[i].move_;
        Value stone_number = GET_STONE_NUMBER(move);
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[GET_CELL(move)];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move, true);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move();
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
    
    MoveInfo * move_infos = MOVE_INFO_STORE[n_moves_made_];
    U32 N = 0;
    get_solver_moves(move_infos, N, true);
    std::sort(move_infos, move_infos + N);
  
    Move move; 
    std::vector<Move> winning_moves;
    U32 i;
    for (i = 0; i < N; i++) {
        move = move_infos[i].move_;
        Value stone_number = GET_STONE_NUMBER(move);
        U32 shift = hash_info.stone_shift_[turn_][stone_number];
        U32 cell_index = hash_info.cell_index_[GET_CELL(move)];
        U64 move_mask = ((U64)cell_index << (U64)shift);
        child_stone_masks[turn_] = stone_masks[turn_] | move_mask;
        make_move(move, true);
        U32 result = solve(end_time, counter, hash_hits, table, hash_info, child_stone_masks);
        unmake_move();
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
        Move best_move = 0;
        Real best_val = -1000.0;
        for (i = 0; i < winning_moves.size(); i++) {
            move = winning_moves[i];
            make_move(move, true);
            Real val = calculate_expectation();
            unmake_move();
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

