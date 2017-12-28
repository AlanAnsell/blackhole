#include "Position.h"


Move::Move(char * move_str, size_t turn) {
    size_t stone_number;
    char cell_str[5];
    sscanf(move_str, "%2s=%u", cell_str, &stone_number);
    cell_ = cell_name_to_id(cell_str);
    stone_value_ = parity[turn] * stone_number;
}


void Move::to_str(char * str) {
    char cell_str[5];
    cell_id_to_name(cell_, cell_str);
    sprintf(str, "%s=%d", cell_str, abs(stone_value_));
}


bool Move::operator < (const Move& other) const {
    if (cell_ != other.cell_)
        return cell_ < other.cell_;
    return stone_value_ < other.stone_value_;
}


Position::Position(size_t blocked[N_BLOCKED_CELLS], Value alpha): turn_(RED), m_(parity[RED]) {
    size_t i, j;
    open_.init(N_CELLS, N_CELLS);
    for (i = 0; i < N_CELLS; i++) {
        value_[i] = 0;
        size_t degree = N_ADJ[i];
        List& adj = adj_[i];
        adj.init(MAX_DEGREE, N_CELLS);
        for (j = 0; j < degree; j++)
            adj.add(ADJ[i][j]);
    }
    bool is_blocked[N_CELLS] = {false};
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        is_blocked[blocked[i]] = true;
    for (i = 0; i < N_CELLS; i++)
        if (! is_blocked[i])
            open_.add(i);
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        fill(blocked[i], 0);

    for (size_t p = 0; p < 2; p++) {
        n_stones_[p] = N_STONES;
        for (i = 0; i < N_STONES; i++) {
            stones_[p][i] = i + 1;
            stone_loc_[p][i+1] = i;
        }
    }
    set_alpha(alpha);
}


void Position::fill(size_t cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (size_t i = 0; i < adj.len_; i++) {
        size_t adj_id = adj.val_[i];
        value_[adj_id] += stone_value;
        adj_[adj_id].remove(cell_id);
    }
}


void Position::unfill(size_t cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (size_t i = 0; i < adj.len_; i++) {
        size_t adj_id = adj.val_[i];
        value_[adj_id] -= stone_value;
        adj_[adj_id].readd(cell_id);
    }
}


void Position::make_move(const Move& move) {
    size_t cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    size_t * stone_loc = stone_loc_[turn_];
    size_t& n_stones = n_stones_[turn_];
    n_stones--;
    size_t i = 0;
    while (stones[i] != stone_number)
        i++;
    for (; i < n_stones; i++) {
        stones[i] = stones[i+1];
        stone_loc[stones[i]] = i;
    }
    turn_ ^= 1;
    m_ *= -1;

    open_.remove(cell_id);
    
    int64 cell_mask = (1LL << (int64)cell_id);
    if (controls_ & cell_mask)
        n_controls_[BLUE]--;
    else
        n_controls_[RED]--;
    for (size_t p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]--;
            if (stale_[p] & cell_mask)
                n_stale_[p]--;
        }
    }
    fill(cell_id, move.stone_value_);
    
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        size_t adj_id = adj.val_[i];
        int64 adj_mask = (1LL << (int64)adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        size_t open_id = open_.val_[i];
        int64 open_mask = (1LL << (int64)open_id);
        if (! (dead_[RED] & open_mask) && ! (dead_[BLUE] & open_mask)) {
            if (is_dead(open_id, RED, power[BLUE])) {
                n_dead_[RED]++;
                dead_[RED] |= open_mask;
            } else if (is_dead(open_id, BLUE, power[RED])) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= open_mask;
            }
        }
    }
    find_stale_cells();
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


void Position::unmake_move(const Move& move) {
    turn_ ^= 1;
    m_ *= -1;
    
    size_t cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    size_t * stone_loc = stone_loc_[turn_];
    size_t& n_stones = n_stones_[turn_];
    size_t i;
    for (i = n_stones; i && stone_number < stones[i-1]; i--) {
        stones[i] = stones[i-1];
        stone_loc[stones[i]] = i;
    }
    stones[i] = stone_number;
    n_stones++;
   
    open_.readd(cell_id); 
    
    int64 cell_mask = (1LL << (int64)cell_id);
    unfill(cell_id, move.stone_value_);
    if (controls_ & cell_mask)
        n_controls_[BLUE]++;
    else
        n_controls_[RED]++;
    for (size_t p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]++;
            if (stale_[p] & cell_mask)
                n_stale_[p]++;
        }
    }
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        size_t adj_id = adj.val_[i];
        int64 adj_mask = (1LL << (int64)adj_id);
        bool new_control = (value_[adj_id] < alpha_);
        bool old_control = (controls_ & adj_mask);
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= adj_mask;
        }
    }
    
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        size_t open_id = open_.val_[i];
        int64 open_mask = (1LL << (int64)open_id);
        if ((dead_[RED] & open_mask)) {
            if (! is_dead(open_id, RED, power[BLUE])) {
                n_dead_[RED]--;
                dead_[RED] ^= open_mask;
            }
        } else if ((dead_[BLUE] & open_mask)) {
            if (! is_dead(open_id, BLUE, power[RED])) {
                n_dead_[BLUE]--;
                dead_[BLUE] ^= open_mask;
            }
        }
    }
    for (i = 0; i < open_.len_; i++) {
        size_t id = open_.val_[i];
        int64 mask = (1LL << (int64)id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            if ((! (dead_[RED] & mask) && ! (dead_[BLUE] & mask)) || ! all_adj_dead(id)) {
                if (stale_[RED] & mask) {
                    stale_[RED] ^= mask;
                    n_stale_[RED]--;
                } else {
                    stale_[BLUE] ^= mask;
                    n_stale_[BLUE]--;
                }
            }
        }
    }
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::get_stone_power(size_t p, Value * power) {
    size_t op = 1 - p;
    size_t n_stones = n_stones_[p];
    Value * stones = stones_[p];
    power[0] = 0;
    size_t i;
    for (i = 1; i <= MAX_DEGREE && i <= n_stones; i++)
        power[i] = power[i-1] + stones[n_stones-i] * parity[p];
    size_t max_our_stones = i;
    for (i = 0; i < n_stones_[op] && i + max_our_stones <= MAX_DEGREE; i++)
        power[i + max_our_stones] = power[i + max_our_stones - 1] + stones_[op][i] * parity[op];
}


Move Position::get_random_move() {
    return Move(open_.val_[rand() % open_.len_],
                m_ * stones_[turn_][rand() % n_stones_[turn_]]);
}


size_t lottery(size_t * tickets, size_t n_tickets) {
    int winner = rand() % n_tickets;
    size_t i;
    for (i = 0; winner >= 0; i++)
        winner -= tickets[i];
    return i - 1;
}


Move Position::get_default_policy_move() {
    size_t tickets[N_CELLS];
    size_t total_tickets = 0;
    size_t i;
    size_t best_stale = NO_CELL;
    size_t cell_id;
    Value best_stale_val = 1000;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        cell_id = open_.val_[i];
        int64 mask = (1LL << (int64)cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            Value stale_val = m_ * value_[cell_id];
            if (stale_val < best_stale_val) {
                best_stale = i;
                best_stale_val = stale_val;
            }
        } else {
            if ((dead_[1-turn_] & mask) || (n_stones_[turn_] > n_dead_[1-turn_])) {
                if ((bool)(controls_ & (1LL << int64(cell_id))) == (bool)turn_)
                    tickets[i] = 1;
                else
                    tickets[i] = 3;
                total_tickets += tickets[i];
            }
        }
    }
    if (best_stale != NO_CELL) {
        tickets[best_stale] = 1;
        total_tickets++;
    }
    
    size_t cell_index = lottery(tickets, total_tickets);
    cell_id = open_.val_[cell_index];

    Value reqs[MAX_DEGREE];
    size_t n_uncontrolled = 0;
    Value req_val = parity[turn_] * (alpha_ - OFFSET[turn_]);
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        Value adj_val = m_ * value_[adj.val_[i]];
        if (adj_val < req_val) {
            reqs[n_uncontrolled] = req_val - adj_val;
            n_uncontrolled++;
        }
    }
    std::sort(reqs, reqs + n_uncontrolled);
    
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    size_t stone_index;
    if (cell_index == best_stale)
        stone_index = 0;
    else {
        total_tickets = 0;
        for (i = n_stale_[1 - turn_]; i < n_stones; i++) {
            size_t changes;
            for (changes = 0; changes < n_uncontrolled && stones[i] <= reqs[changes]; changes++) {}
            tickets[i] = 20 * (changes + 1) - stones[i];
            total_tickets += tickets[i];
        }
#ifdef DEBUG_
        if (total_tickets == 0) {
            print(stderr);
            fprintf(stderr, "Red: %u stale, %u dead\n", n_stale_[RED], n_dead_[RED]);
            fprintf(stderr, "Blue: %u stale, %u dead\n", n_stale_[BLUE], n_dead_[BLUE]);
            fflush(stderr);
        }
        assert(total_tickets > 0);
#endif
        stone_index = lottery(tickets, total_tickets);
    }
    
    Value stone_value = m_ * stones[stone_index];
    return Move(cell_id, stone_value);
}


Move Position::get_expectation_maximising_move() {
    Move best_move;
    Real best_expectation = -1000.0;
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real  expectation = calculate_expectation();
            unmake_move(move);
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
    size_t i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += parity[i] * stones_[i][j];
    for (i = 0; i < open_.len_; i++) {
        size_t cell_id = open_.val_[i];
        sum += value_[cell_id] + adj_[cell_id].len_ * stone_sum / (open_.len_ - 1);
    }
    return sum / open_.len_;
}


Real Position::get_control_heuristic() const {
    return (Real)n_controls_[RED] / (Real)open_.len_;
}


void Position::get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves) {
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real expectation = calculate_expectation();
            unmake_move(move);
            moves.push_back(std::make_pair(m_ * expectation, move));
        }
    }
}


void Position::get_all_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            moves.push_back(move);
        }
    }
}


bool Position::is_dead(size_t cell_id, size_t p, Value * other_power) {
    Value value = value_[cell_id];
    size_t n_adj = adj_[cell_id].len_;
    Value worst_val = parity[p] * (value + other_power[n_adj]);
    return (worst_val >= parity[p] * (alpha_ - OFFSET[p]));
}


bool Position::all_adj_dead(size_t cell_id) {
    List& adj = adj_[cell_id];
    size_t n_adj = adj.len_;
    for (size_t i = 0; i < n_adj; i++) {
        int64 adj_mask = (1LL << (int64)adj.val_[i]);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
            return false;
    }
    return true;
}


void Position::find_stale_cells() {
    for (size_t i = 0; i < open_.len_; i++) {
        size_t cell_id = open_.val_[i];
        int64 mask = (1LL << (int64)cell_id);
        if (! (stale_[RED] & mask) && ! (stale_[BLUE] & mask) &&
                ((dead_[RED] & mask) || (dead_[BLUE] & mask))) {
            if (all_adj_dead(cell_id)) {
                if (dead_[RED] & mask) {
                    stale_[RED] |= mask;
                    n_stale_[RED]++;
                } else {
                    stale_[BLUE] |= mask;
                    n_stale_[BLUE]++;
                }
            }
        }
    }
}


bool Position::is_winning(size_t p) const {
    return n_dead_[p] > n_stones_[1 - p];
}


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    Value power[2][MAX_DEGREE+1];
    dead_[0] = dead_[1] = 0;
    stale_[RED] = stale_[BLUE] = 0;
    n_stale_[RED] = n_stale_[BLUE] = 0;
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    size_t i, cell_id, n_adj;
    int64 mask;
    for (i = 0; i < open_.len_; i++) {
        cell_id = open_.val_[i];
        Value value = value_[cell_id];
        n_adj = adj_[cell_id].len_;
        mask = (1LL << (int64)cell_id);
        if (value >= alpha_) {
            n_controls_[RED]++;
            if (value + power[BLUE][n_adj] >= alpha_) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            }
        } else {
            controls_ |= mask;
            n_controls_[BLUE]++;
            if (value + power[RED][n_adj] < alpha_) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= mask;
            }
        }
    }
    find_stale_cells();
#ifdef DEBUG_
    assert((dead_[RED] & dead_[BLUE]) == 0LL);
    assert((stale_[RED] & stale_[BLUE]) == 0LL);
    assert(n_stale_[RED] <= n_dead_[RED]);
    assert(n_stale_[BLUE] <= n_dead_[BLUE]);
#endif
}


void Position::take_snapshot() {
    size_t i, j, p;
    snapshot_.turn_ = turn_;
    snapshot_.n_open_ = open_.len_;
    snapshot_.controls_ = controls_;
    for (i = 0; i < open_.len_; i++) {
        size_t cell_id = open_.val_[i];
        snapshot_.open_[i] = cell_id;
        snapshot_.value_[cell_id] = value_[cell_id];
        List& adj = adj_[cell_id];
        size_t n_adj = adj.len_;
        snapshot_.n_adj_[cell_id] = n_adj;
        for (j = 0; j < n_adj; j++)
            snapshot_.adj_[cell_id][j] = adj.val_[j];
    }

    for (p = 0; p < 2; p++) {
        size_t n_stones = n_stones_[p];
        snapshot_.n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        Value * snapshot_stones = snapshot_.stones_[p];
        for (i = 0; i < n_stones; i++)
            snapshot_stones[i] = stones[i];

        snapshot_.n_controls_[p] = n_controls_[p];
        snapshot_.dead_[p] = dead_[p];
        snapshot_.n_dead_[p] = n_dead_[p];
        snapshot_.stale_[p] = stale_[p];
        snapshot_.n_stale_[p] = n_stale_[p];
    }

}


void Position::restore_snapshot() {
    size_t i, j, p;
    turn_ = snapshot_.turn_;
    m_ = parity[turn_];
    open_.len_ = snapshot_.n_open_;
    controls_ = snapshot_.controls_;
    for (i = 0; i < snapshot_.n_open_; i++) {
        size_t cell_id = snapshot_.open_[i];
        value_[cell_id] = snapshot_.value_[cell_id];
        open_.val_[i] = cell_id;
        open_.loc_[cell_id] = i;
        size_t n_adj = snapshot_.n_adj_[cell_id];
        List& adj = adj_[cell_id];
        adj.len_ = n_adj; 
        for (j = 0; j < n_adj; j++) {
            size_t adj_id = snapshot_.adj_[cell_id][j];
            adj.val_[j] = adj_id;
            adj.loc_[adj_id] = j;
        }
    }
    
    for (p = 0; p < 2; p++) {
        size_t n_stones = snapshot_.n_stones_[p];
        n_stones_[p] = n_stones;
        Value * stones = stones_[p];
        size_t * stone_loc = stone_loc_[p];
        Value * snapshot_stones = snapshot_.stones_[p];
        for (i = 0; i < n_stones; i++) {
            stones[i] = snapshot_stones[i];
            stone_loc[stones[i]] = i;
        }

        n_controls_[p] = snapshot_.n_controls_[p];
        dead_[p] = snapshot_.dead_[p];
        n_dead_[p] = snapshot_.n_dead_[p];
        stale_[p] = snapshot_.stale_[p];
        n_stale_[p] = snapshot_.n_stale_[p];
    }
}


std::pair<size_t, size_t> Position::get_cell_and_stone_indices(const Move& move) {
    return std::pair<size_t, size_t>(open_.loc_[move.cell_],
                                     stone_loc_[turn_][m_ * move.stone_value_]);
}


void Position::print(FILE * f) {
    size_t i, j;
    for (i = 0; i < open_.len_; i++) {
        size_t cell_id = open_.val_[i];
        char cell_str[10];
        cell_id_to_name(cell_id, cell_str);
        fprintf(f, "%s (%d):", cell_str, value_[cell_id]);
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            cell_id_to_name(adj.val_[j], cell_str);
            fprintf(f, " %s", cell_str);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "Red:");
    for (i = 0; i < n_stones_[RED]; i++)
        fprintf(f, " %d", stones_[RED][i] * parity[RED]);
    fprintf(stderr, "\nBlue:");
    for (i = 0; i < n_stones_[BLUE]; i++)
        fprintf(stderr, " %d", stones_[BLUE][i] * parity[BLUE]);
    fprintf(stderr, "\n");
}

