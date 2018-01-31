#include "Position.h"
#include "AMAF.h"


Move move_from_str(char * move_str) {
    U32 stone_number;
    char cell_str[5];
    sscanf(move_str, "%2s=%u", cell_str, &stone_number);
    return CREATE_MOVE(cell_name_to_id(cell_str), stone_number);
}


void move_to_str(Move move, char * str) {
    char cell_str[5];
    cell_id_to_name(GET_CELL(move), cell_str);
    sprintf(str, "%s=%d", cell_str, GET_STONE_NUMBER(move));
}


Position::Position(U32 blocked[N_BLOCKED_CELLS], Value alpha):
        turn_(RED), n_moves_made_(0), m_(PARITY[RED]) {
    U32 i, j;
    open_.init(N_CELLS, N_CELLS);
    open_mask_ = 0;
    for (i = 0; i < N_CELLS; i++) {
        value_[i] = 0;
        U32 degree = N_ADJ[i];
        List& adj = adj_[i];
        adj.init(MAX_DEGREE, N_CELLS);
        for (j = 0; j < degree; j++)
            adj.add(ADJ[i][j]);
    }
    bool is_blocked[N_CELLS] = {false};
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        is_blocked[blocked[i]] = true;
    for (i = 0; i < N_CELLS; i++) {
        if (! is_blocked[i]) {
            open_mask_ |= MASK(i);
            open_.add(i);
        }
    }
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        fill(blocked[i], 0);

    for (U32 p = 0; p < 2; p++) {
        n_stones_[p] = N_STONES;
        stone_masks_[p] = (1 << N_STONES) - 1;
        power_[p] = STONE_POWER[p][stone_masks_[p]];
        for (i = 0; i < N_STONES; i++) {
            stones_[p][i] = i + 1;
            stone_loc_[p][i+1] = i;
        }
    }
    set_alpha(alpha);
}


bool Position::is_legal(Move move) {
    U32 cell_id = GET_CELL(move);
    U32 stone_number = GET_STONE_NUMBER(move);
    return cell_id >= 0 && cell_id < N_CELLS &&
            (open_mask_ & MASK(cell_id)) &&
            stone_number >= 1 && stone_number <= N_STONES &&
            (stone_masks_[turn_] & (1 << (stone_number - 1)));
}


void Position::save_history() {
    PositionHistory& hist = history_[n_moves_made_];
    
    hist.open_mask_ = open_mask_;

    hist.power_[0] = power_[0];
    hist.power_[1] = power_[1];

    hist.controls_ = controls_;
    hist.n_controls_[0] = n_controls_[0];
    hist.n_controls_[1] = n_controls_[1];

    hist.dead_[0] = dead_[0];
    hist.dead_[1] = dead_[1];
    hist.n_dead_[0] = n_dead_[0];
    hist.n_dead_[1] = n_dead_[1];

    hist.stale_[0] = stale_[0];
    hist.stale_[1] = stale_[1];
    hist.n_stale_[0] = n_stale_[0];
    hist.n_stale_[1] = n_stale_[1];
    hist.either_stale_ = either_stale_;
    hist.worst_stale_ = worst_stale_;

    hist.stone_masks_[0] = stone_masks_[0];
    hist.stone_masks_[1] = stone_masks_[1];

    U32 i;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        hist.effective_adj_[cell_id] = effective_adj_[cell_id];
        hist.value_[i] = value_[i];
    }

    for (U32 p = 0; p < 2; p++) {
        U32 n_stones = n_stones_[p];
        hist.n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        Value * hist_stones = hist.stones_[p];
        for (i = 0; i < n_stones; i++)
            hist_stones[i] = stones[i];
    }

}


void Position::restore_history() {
    // regress_move should already have been called to restore
    // n_moves_made_ and turn_ to the correct value
    PositionHistory& hist = history_[n_moves_made_];
    
    open_mask_ = hist.open_mask_;
    controls_ = hist.controls_;

    power_[0] = hist.power_[0];
    power_[1] = hist.power_[1];

    n_controls_[0] = hist.n_controls_[0];
    n_controls_[1] = hist.n_controls_[1];

    dead_[0] = hist.dead_[0];
    dead_[1] = hist.dead_[1];
    n_dead_[0] = hist.n_dead_[0];
    n_dead_[1] = hist.n_dead_[1];

    stale_[0] = hist.stale_[0];
    stale_[1] = hist.stale_[1];
    n_stale_[0] = hist.n_stale_[0];
    n_stale_[1] = hist.n_stale_[1];
    either_stale_ = hist.either_stale_;
    worst_stale_ = hist.worst_stale_;

    stone_masks_[0] = hist.stone_masks_[0];
    stone_masks_[1] = hist.stone_masks_[1];

    U32 i;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        effective_adj_[cell_id] = hist.effective_adj_[cell_id];
        value_[i] = hist.value_[i];
    }

    for (U32 p = 0; p < 2; p++) {
        U32 n_stones = hist.n_stones_[p];
        n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        U32 * stone_loc = stone_loc_[p];
        Value * hist_stones = hist.stones_[p];
        for (i = 0; i < n_stones; i++) {
            stones[i] = hist_stones[i];
            stone_loc[hist_stones[i]] = i;
        }
    }
    

}


void Position::make_move(Move move, bool save) {
#ifdef DEBUG_
    if (! is_legal(move)) {
        char move_str[20];
        fprintf(stderr, "Making move\n");
        print(stderr);
        move_to_str(move, move_str);
        fprintf(stderr, "%s\n", move_str);
        fflush(stderr);
    }
    assert(is_legal(move));
#endif
    if (save)
        save_history();
    moves_made_[n_moves_made_] = move;
    U32 cell_id = GET_CELL(move);
    Value stone_number = (Value)GET_STONE_NUMBER(move);
    Value stone_value = m_ * stone_number;
    Value * stones = stones_[turn_];
    U32 * stone_loc = stone_loc_[turn_];
    U32& n_stones = n_stones_[turn_];
    n_stones--;
    U32 i = 0;
    while (stones[i] != stone_number)
        i++;
    for (; i < n_stones; i++) {
        stones[i] = stones[i+1];
        stone_loc[stones[i]] = i;
    }
    stone_masks_[turn_] ^= (1 << (stone_number - 1));
    power_[turn_] = STONE_POWER[turn_][stone_masks_[turn_]];
    turn_ ^= 1;
    m_ *= -1;
    n_moves_made_++;

    U64 cell_mask = MASK(cell_id);
    open_mask_ ^= cell_mask;
    open_.remove(cell_id);
    
    if (controls_ & cell_mask)
        n_controls_[BLUE]--;
    else
        n_controls_[RED]--;
    for (U32 p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]--;
            if (stale_[p] & cell_mask)
                n_stale_[p]--;
        }
    }
    fill(cell_id, stone_value);
    
    List& adj = adj_[cell_id];
    //bool is_stale = (cell_mask & stale_[RED]) || (cell_mask & stale_[BLUE]);
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = MASK(adj_id);
        //if (! is_stale)
        effective_adj_[adj_id]--;
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    
    for (i = 0; i < open_.len_; i++) {
        U32 open_id = open_.val_[i];
        U64 open_mask = MASK(open_id);
        if (! (dead_[RED] & open_mask) && ! (dead_[BLUE] & open_mask)) {
            if (is_dead(open_id, RED)) {
                n_dead_[RED]++;
                dead_[RED] |= open_mask;
            } else if (is_dead(open_id, BLUE)) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= open_mask;
            }
        }
    }
    find_stale_cells();
    find_worst_stale();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    if (n_stale_[RED] > n_dead_[RED]) {
        print(stderr);
        fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[RED], n_dead_[RED]);
        fflush(stderr);
    }
    if (n_stale_[BLUE] > n_dead_[BLUE]) {
        print(stderr);
        fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[BLUE], n_dead_[BLUE]);
        fflush(stderr);
    }
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::regress_board() {
    n_moves_made_--;
    turn_ ^= 1;
    m_ *= -1;
        
    Move move = moves_made_[n_moves_made_]; 
    U32 cell_id = GET_CELL(move);
    Value stone_number = GET_STONE_NUMBER(move);
    open_.readd(cell_id); 
    unfill(cell_id, m_ * stone_number);
}


void Position::rewind_to(U32 n_moves_made) {
    while (n_moves_made_ > n_moves_made)
        regress_board();
    restore_history();
}


inline U32 lottery(U32 * tickets, U32 n_tickets) {
    int winner = rand() % n_tickets;
    U32 i;
    for (i = 0; winner >= 0; i++)
        winner -= tickets[i];
    return i - 1;
}


Move Position::get_default_policy_move(U64 valid, U64 duo) {
    U32 tickets[N_CELLS];
    U32 total_tickets = 0;
    U32 i, j;
    U32 cell_id;
    U32 op = 1 - turn_;
    U64 mask;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        cell_id = open_.val_[i];
        mask = MASK(cell_id);
        if ((valid & mask) && (cell_id != worst_stale_)) {
            int swing = 0;
            if (dead_[turn_] & mask)
                swing = -1;
            else if (dead_[op] & mask)
                swing = 1;
            List& adj = adj_[cell_id];
            for (j = 0; j < adj.len_; j++) {
                U32 adj_id = adj.val_[j];
                U64 adj_mask = MASK(adj_id);
                if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
                    swing++;
            }
            tickets[i] = swing + 2;
            total_tickets += tickets[i];
        }
    }
    if (worst_stale_ != NO_CELL) {
        tickets[open_.loc_[worst_stale_]] = 1;
        total_tickets++;
    }
#ifdef DEBUG_
    if (total_tickets == 0)
        print(stderr);
    assert(total_tickets > 0);
#endif
     
    U32 cell_index = lottery(tickets, total_tickets);
    cell_id = open_.val_[cell_index];
    mask = MASK(cell_id);
    
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    U32 stone_index;
    if (cell_id == worst_stale_)
        stone_index = 0;
    else {
        if (duo & mask) {
            tickets[0] = 1;
            total_tickets = 1;
            U32 adj_id = get_non_stale_adj(adj_[cell_id], 0);
            Value adj_value = value_[adj_id] * m_;
            Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
            for (i = 1; i < n_stones && stones[i] < extra_req; i++) {}
            if (i < n_stones) {
                tickets[i] = 3;
                total_tickets += 3;
            }
        } else {
            total_tickets = 0;
            for (i = n_stale_[op]; i < n_stones; i++) {
                //U32 changes;
                //for (changes = 0; changes < n_uncontrolled && stones[i] <= reqs[changes]; changes++) {}
                tickets[i] = /*20 * (changes + 1) -*/ stones[i];
                total_tickets += tickets[i];
            }
        }
#ifdef DEBUG_
        if (total_tickets == 0) {
            char cell_str[5];
            cell_id_to_name(cell_id, cell_str);
            fprintf(stderr, "Cell chosen: %s\n", cell_str);
            print(stderr);
            fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[RED], n_dead_[RED]);
            fprintf(stderr, "Blue: %u stale, %u dead\n", n_stale_[BLUE], n_dead_[BLUE]);
            fflush(stderr);
        }
        assert(total_tickets > 0);
#endif
        stone_index = lottery(tickets, total_tickets);
    }
    
    return CREATE_MOVE(cell_id, stones[stone_index]);
}


Move Position::get_expectation_maximising_move() {
    Move best_move = 0;
    Real best_expectation = -1000.0;
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Move move = CREATE_MOVE(open_.val_[i], stones[j]);
            make_move(move, true);
            Real expectation = calculate_expectation();
            unmake_move();
            expectation *= m_;
            if (expectation > best_expectation) {
                best_move = move;
                best_expectation = expectation;
            }
        }
    }
    return best_move;
}



Real Position::calculate_expectation() const {
    if (open_.len_ == 1)
        return (Real)(value_[open_.val_[0]]);
    Real sum = 0.0;
    Real stone_sum = 0.0;
    U32 i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += PARITY[i] * stones_[i][j];
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        sum += value_[cell_id] + adj_[cell_id].len_ * stone_sum / (open_.len_ - 1);
    }
    return sum / open_.len_;
}


Move Position::get_best_winning_move() {
    U32 winner = RED;
    if (is_winning(RED))
        winner = RED;
    else if (is_winning(BLUE))
        winner = BLUE;
#ifdef DEBUG_
    else
        assert(false);
#endif
    U32 n_stones = n_stones_[turn_];
    Value * stones = stones_[turn_];
    Move best_move = 0;
    Real best_value = -1000.0;
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        for (U32 j = 0; j < n_stones; j++) {
            Move move = CREATE_MOVE(cell_id, stones[j]);
            make_move(move, true);
            bool still_winning = is_winning(winner);
            Real value = calculate_expectation();
            unmake_move();
            if (still_winning) {
                value *= m_;
                if (value > best_value) {
                    best_move = move;
                    best_value = value;
                }
            }
        }
    }
#ifdef DEBUG_
    assert(best_value > -1000.0);
#endif
    return best_move;
}


void Position::get_all_reasonable_moves_with_stone(U32 stone_index, U64 valid, U64 duo, Move * moves, U32& N) {
    N = 0;
    Value stone_number = stones_[turn_][stone_index];
    Value prev_stone_number;
    if (stone_index == 0)
        prev_stone_number = 0;
    else
        prev_stone_number = stones_[turn_][stone_index - 1];

    U32 op = 1 - turn_;

    while (valid) {
        U64 lsb = LSB(valid);
        valid ^= lsb;
        U32 cell_id = INDEX(lsb);
        if (cell_id == worst_stale_) {
            if (stone_index == 0)
                moves[N++] = CREATE_MOVE(cell_id, stone_number);
        } else {
            if (duo & lsb) {
                U32 adj_id = get_non_stale_adj(adj_[cell_id], 0);
                Value adj_value = value_[adj_id] * m_;
                if (stone_index == 0)
                    moves[N++] = CREATE_MOVE(cell_id, stone_number);
                else {
                    Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
                    if (stone_number >= extra_req && prev_stone_number < extra_req)
                        moves[N++] = CREATE_MOVE(cell_id, stone_number);
                }
            } else {
                if (stone_index >= n_stale_[op])
                    moves[N++] = CREATE_MOVE(cell_id, stone_number);
            }
        }
    }
}


#ifdef SOLVER_IMPL_
#include SOLVER_IMPL_
#else
#include "solver.h"
#endif

void Position::generate_untried_moves(int& stone_index, U64 valid, U64 duo,
        std::vector<Move>& untried_children) {
    // untried_children should be empty
#ifdef DEBUG_
    assert(untried_children.empty());
#endif
    //if (valid == 0)
    //    get_validity_mask(valid, duo);
    U32 i;
    U32 op = 1 - turn_;
    while (untried_children.size() == 0 && stone_index >= 0) {
        U32 N;
        Move moves[N_CELLS];
        MoveInfo move_infos[N_CELLS];
        get_all_reasonable_moves_with_stone(stone_index, valid, duo, moves, N);
        Value stone_number = stones_[turn_][stone_index];
        for (i = 0; i < N; i++) {
            U32 cell_id = GET_CELL(moves[i]);
            U64 mask = MASK(cell_id);
            int swing = adj_[cell_id].len_;
            int kill = 0;
            if (dead_[turn_] & mask)
                kill = -1;
            else if (dead_[op] & mask)
                kill = 1;
            List& adj = adj_[cell_id];
            for (U32 j = 0; j < adj.len_; j++) {
                U32 adj_id = adj.val_[j];
                U64 adj_mask = MASK(adj_id);
                if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask)) {
                    U32 n_adj = adj_[adj_id].len_;
                    Value adj_value = (power_[op][std::min(n_adj-1, n_stones_[op])] + value_[adj_id]) * m_;
                    Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                    if (adj_value + stone_number >= control_req)
                        kill++;
                }
            }
            move_infos[i] = MoveInfo(moves[i], kill, swing);
        }
        std::sort(move_infos, move_infos + N);
        for (i = N; i > 0; i--)
            untried_children.push_back(move_infos[i-1].move_);
        stone_index--;
    }
}


/*bool Position::all_adj_dead(U32 cell_id) {
    List& adj = adj_[cell_id];
    U32 n_adj = adj.len_;
    for (U32 i = 0; i < n_adj; i++) {
        U64 adj_mask = MASK(adj.val_[i]);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
            return false;
    }
    return true;
}*/


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    dead_[0] = dead_[1] = 0;
    U32 i, cell_id, n_adj;
    U64 mask;
    for (i = 0; i < open_.len_; i++) {
        cell_id = open_.val_[i];
        Value value = value_[cell_id];
        n_adj = adj_[cell_id].len_;
        mask = MASK(cell_id);
        if (value >= alpha_) {
            n_controls_[RED]++;
            if (value + power_[BLUE][std::min(n_adj, n_stones_[BLUE])] >= alpha_) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            }
        } else {
            controls_ |= mask;
            n_controls_[BLUE]++;
            if (value + power_[RED][std::min(n_adj, n_stones_[RED])] < alpha_) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= mask;
            }
        }
    }
    stale_[RED] = stale_[BLUE] = 0;
    n_stale_[RED] = n_stale_[BLUE] = 0;
    either_stale_ = 0;
    find_stale_cells();
    find_worst_stale();
    find_effective_adj();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::print(FILE * f) const {
    U32 i, j;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        char cell_str[10];
        cell_id_to_name(cell_id, cell_str);
        fprintf(f, "%s (%d; %u):", cell_str, value_[cell_id], effective_adj_[cell_id]);
        const List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            cell_id_to_name(adj.val_[j], cell_str);
            fprintf(f, " %s", cell_str);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "Red:");
    for (i = 0; i < n_stones_[RED]; i++)
        fprintf(f, " %d", stones_[RED][i] * PARITY[RED]);
    fprintf(f, "\nBlue:");
    for (i = 0; i < n_stones_[BLUE]; i++)
        fprintf(f, " %d", stones_[BLUE][i] * PARITY[BLUE]);
    fprintf(f, "\n");
    fprintf(f, "Alpha = %d\n", alpha_);
    fprintf(f, "open_mask_ = 0x%llx\n", open_mask_);
    fprintf(f, "controls_ = 0x%llx\n", controls_);
    fprintf(f, "dead_[RED] = 0x%llx\n", dead_[RED]);
    fprintf(f, "dead_[BLUE] = 0x%llx\n", dead_[BLUE]);
    fprintf(f, "stale_[RED] = 0x%llx\n", stale_[RED]);
    fprintf(f, "stale_[BLUE] = 0x%llx\n", stale_[BLUE]);
    fprintf(f, "stone_masks_[RED] = 0x%x\n", stone_masks_[RED]);
    fprintf(f, "stone_masks_[BLUE] = 0x%x\n", stone_masks_[BLUE]);
    fprintf(f, "Worst stale cell: ");
    if (worst_stale_ == NO_CELL)
        fprintf(f, "none\n");
    else {
        char cell_str[5];
        cell_id_to_name(worst_stale_, cell_str);
        fprintf(f, "%s\n", cell_str);
    }
}

