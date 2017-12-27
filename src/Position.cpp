#include "Position.h"


void Cell::init(CellID id, Cell * cells) {
    id_ =  id;
    value_ = 0;
    board_ = cells;
    size_t degree = N_ADJ[id_];
    for (size_t i = 0; i < degree; i++)
        adj_.add(ADJ[id_][i]);
}


void Cell::remove_neighbour(CellID neighbour_id, Value stone_value) {
    value_ += stone_value;
    adj_.remove(neighbour_id);
}
	

void Cell::add_neighbour(CellID neighbour_id, Value stone_value) {
    value_ -= stone_value;
    adj_.readd(neighbour_id);
}


void Cell::fill(Value stone_value) {
    for (size_t i = 0; i < adj_.len_; i++)
        board_[adj_.val_[i]].remove_neighbour(id_, stone_value);
}


void Cell::unfill(Value stone_value) {
    for (size_t i = 0; i < adj_.len_; i++)
        board_[adj_.val_[i]].add_neighbour(id_, stone_value);
}


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


Position::Position(CellID blocked[N_BLOCKED_CELLS], Value alpha): open_(N_CELLS, N_CELLS), turn_(RED) {
    size_t i;
    bool is_blocked[N_CELLS] = {false};
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        is_blocked[blocked[i]] = true;
    for (i = 0; i < N_CELLS; i++) {
        cells_[i].init(i, cells_);
        if (! is_blocked[i])
            open_.add(i);
    }
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        cells_[blocked[i]].fill(0);

    for (size_t p = 0; p < 2; p++) {
        n_stones_[p] = N_STONES;
        for (i = 0; i < N_STONES; i++)
            stones_[p][i] = i + 1;
    }
    set_alpha(alpha);
}

void Position::make_move(const Move& move) {
    Value stone_number = parity[turn_] * move.stone_value_;
    Value * stones = stones_[turn_];
    size_t& n_stones = n_stones_[turn_];
    n_stones--;
    size_t i = 0;
    while (stones[i] != stone_number)
        i++;
    for (; i < n_stones; i++)
        stones[i] = stones[i+1];
    turn_ ^= 1;
    
    open_.remove(move.cell_);
    
    if (controls_ & (1LL << (int64)move.cell_))
        n_controls_[BLUE]--;
    else
        n_controls_[RED]--;
    for (size_t p = 0; p < 2; p++)
        if (dead_[p] & (1LL << (int64)move.cell_))
            n_dead_[p]--;
    Cell& cell = cells_[move.cell_];
    cell.fill(move.stone_value_);
    
    for (i = 0; i < cell.adj_.len_; i++) {
        CellID id = cell.adj_.val_[i];
        Cell& adj = cells_[id];
        bool new_control = (adj.value_ < alpha_);
        bool old_control = (controls_ & (1LL << (int64)id));
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= (1LL << (int64)id);
        }
    }
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        CellID id = open_.val_[i];
        int64 mask = (1LL << (int64)id);
        if (! (dead_[RED] & mask) && ! (dead_[BLUE] & mask)) {
            Cell& cell = cells_[id];
            if (is_dead(cell, RED, power[BLUE])) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            } else if (is_dead(cell, BLUE, power[RED])) {
                n_dead_[BLUE]++;
                dead_[BLUE] |= mask;
            }
        }
    }
}


void Position::unmake_move(const Move& move) {
    turn_ ^= 1;
    
    Value stone_number = parity[turn_] * move.stone_value_;
    Value * stones = stones_[turn_];
    size_t& n_stones = n_stones_[turn_];
    size_t i;
    for (i = n_stones; i && stone_number < stones[i-1]; i--)
        stones[i] = stones[i-1];
    stones[i] = stone_number;
    n_stones++;
   
    open_.readd(move.cell_); 
    
    Cell& cell = cells_[move.cell_];
    cell.unfill(move.stone_value_);
    if (controls_ & (1LL << (int64)move.cell_))
        n_controls_[BLUE]++;
    else
        n_controls_[RED]++;
    for (size_t p = 0; p < 2; p++)
        if (dead_[p] & (1LL << (int64)move.cell_))
            n_dead_[p]++;
    for (i = 0; i < cell.adj_.len_; i++) {
        CellID id = cell.adj_.val_[i];
        Cell& adj = cells_[id];
        bool new_control = (adj.value_ < alpha_);
        bool old_control = (controls_ & (1LL << (int64)id));
        if (new_control != old_control) {
            n_controls_[old_control]--;
            n_controls_[new_control]++;
            controls_ ^= (1LL << (int64)id);
        }
    }
    
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (i = 0; i < open_.len_; i++) {
        CellID id = open_.val_[i];
        Cell& cell = cells_[id];
        int64 mask = (1LL << (int64)id);
        if ((dead_[RED] & mask)) {
            if (! is_dead(cell, RED, power[BLUE])) {
                n_dead_[RED]--;
                dead_[RED] ^= mask;
            }
        } else if ((dead_[BLUE] & mask)) {
            if (! is_dead(cell, BLUE, power[RED])) {
                n_dead_[BLUE]--;
                dead_[BLUE] ^= mask;
            }
        }
    }
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


size_t Position::dead_endgame_solve(Value alpha) {
    size_t controls[2] = {0};
    for (size_t p = 0; p < 2; p++) {
        size_t op = 1 - p;
        Value power[MAX_DEGREE+1];
        get_stone_power(p, power);
        for (size_t i = 0; i < open_.len_; i++) {
            Cell& cell = cells_[open_.val_[i]];
            Value best_val = cell.value_ + power[cell.adj_.len_];
            if (best_val * parity[p] < (alpha - OFFSET[p]) * parity[p])
                controls[op]++;
        }
    }
    if (controls[RED] > n_stones_[BLUE])
        return RED;
    if (controls[BLUE] > n_stones_[RED])
        return BLUE;
    return NO_RESULT;
}


Value Position::dead_endgame_value() {
    Value lb = MIN_RESULT;
    Value ub = MAX_RESULT;
    while (lb < ub) {
        Value mid = (lb + ub + 1) / 2;
        size_t result = dead_endgame_solve(mid);
        if (result == NO_RESULT)
            return NO_RESULT;
        if (result == RED)
            lb = mid;
        else
            ub = mid - 1;
    }
    return lb;
}


Move Position::get_random_move() {
    return Move(open_.val_[rand() % open_.len_],
                parity[turn_] * stones_[turn_][rand() % n_stones_[turn_]]);
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
    size_t best_single = 400;
    Value best_single_val = 1000;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        CellID id = open_.val_[i];
        Cell& cell = cells_[id];
        if (cell.adj_.len_ == 0) {
            Value single_val = parity[turn_] * cell.value_;
            if (single_val < best_single_val) {
                best_single = i;
                best_single_val = single_val;
            }
        } else {
            if ((bool)(controls_ & (1LL << int64(id))) == (bool)turn_)
                tickets[i] = 3;
            else
                tickets[i] = 9;
            total_tickets += tickets[i];
        }
    }
    if (best_single != 400) {
        tickets[best_single] = 1;
        total_tickets++;
    }
    
    size_t cell_index = lottery(tickets, total_tickets);
    CellID pid = open_.val_[cell_index];
    Cell& pcell = cells_[pid];

    Value reqs[MAX_DEGREE];
    size_t n_uncontrolled = 0;
    Value req_val = parity[turn_] * (alpha_ - OFFSET[turn_]);
    for (i = 0; i < pcell.adj_.len_; i++) {
        Value adj_val = parity[turn_] * cells_[pcell.adj_.val_[i]].value_;
        if (adj_val < req_val) {
            reqs[n_uncontrolled] = req_val - adj_val;
            n_uncontrolled++;
        }
    }
    std::sort(reqs, reqs + n_uncontrolled);
    
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    total_tickets = 0;
    for (i = 0; i < n_stones; i++) {
        size_t changes;
        for (changes = 0; changes < n_uncontrolled && stones[i] <= reqs[changes]; changes++) {}
        tickets[i] = 20 * (changes + 1) - stones[i];
        total_tickets += tickets[i];
    }

    size_t stone_index = lottery(tickets, total_tickets);
    Value stone_value = stones[stone_index] * parity[turn_];
    return Move(pid, stone_value);
}


Move Position::get_expectation_maximising_move() {
    Move best_move;
    Real best_expectation = -1000.0;
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real  expectation = calculate_expectation();
            unmake_move(move);
            expectation *= parity[turn_];
            if (expectation > best_expectation) {
                best_move = move;
                best_expectation = expectation;
            }
        }
    }
    return best_move;
}


Real Position::search_expectation() {
    //printf("Entering search_expectation\n");
    //fflush(stdout);
    if (open_.len_ == 1) {
        //print();
        return (Real)(cells_[open_.val_[0]].value_);
    }
    Real sum = 0.0;
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            sum += search_expectation();
            unmake_move(move);
        }
    }
    //printf("Exiting search_expectation\n");
    //fflush(stdout);
    return sum / (open_.len_ * n_stones);
}


Real Position::calculate_expectation() const {
    if (open_.len_ == 1)
        return (Real)(cells_[open_.val_[0]].value_);
    Real sum = 0.0;
    Real stone_sum = 0.0;
    size_t i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += parity[i] * stones_[i][j];
    for (i = 0; i < open_.len_; i++) {
        const Cell& cell = cells_[open_.val_[i]];
        sum += cell.value_ + cell.adj_.len_ * stone_sum / (open_.len_ - 1);
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
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real expectation = calculate_expectation();
            unmake_move(move);
            moves.push_back(std::make_pair(parity[turn_] * expectation, move));
        }
    }
}


void Position::get_all_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < open_.len_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_.val_[i], stone_value);
            moves.push_back(move);
        }
    }
}


bool Position::is_dead(Cell& cell, size_t p, Value * other_power) {
    Value value = cell.value_;
    size_t n_adj = cell.adj_.len_;
    Value worst_val = parity[p] * (value + other_power[n_adj]);
    return (worst_val >= parity[p] * (alpha_ - OFFSET[p]));
}


bool Position::is_winning(size_t p) const {
    size_t op = 1 - p;
    return n_dead_[p] > n_stones_[op];
}


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    dead_[0] = dead_[1] = 0;
    Value power[2][MAX_DEGREE+1];
    get_stone_power(RED, power[0]);
    get_stone_power(BLUE, power[1]);
    for (size_t i = 0; i < open_.len_; i++) {
        CellID id = open_.val_[i];
        Cell& cell = cells_[id];
        Value value = cell.value_;
        size_t n_adj = cell.adj_.len_;
        int64 mask = (1LL << (int64)id);
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
}


void Position::print(FILE * f) {
    size_t i, j;
    for (i = 0; i < open_.len_; i++) {
        char cell_str[10];
        cell_id_to_name(open_.val_[i], cell_str);
        Cell& cell = cells_[open_.val_[i]];
        fprintf(f, "%s (%d):", cell_str, cell.value_);
        for (j = 0; j < cell.adj_.len_; j++) {
            cell_id_to_name(cell.adj_.val_[j], cell_str);
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


/*size_t n_dead(size_t player, Value alpha) {
    Value * stones = stones_[player];
    size_t n_stones = n_stones_[player];
    Value cuml[N_STONES+1];
    cuml[0] = 0;
    size_t i;
    for (i = 0; i < n_stones; i++)
        cuml[i+1] = cuml[i] + stones[i];*/

