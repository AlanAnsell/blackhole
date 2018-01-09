#include "Position.h"
#include "AMAF.h"


Move::Move(char * move_str, U32 turn) {
    U32 stone_number;
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


Position::Position(U32 blocked[N_BLOCKED_CELLS], Value alpha): turn_(RED), m_(parity[RED]) {
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


void Position::fill(U32 cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (U32 i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        value_[adj_id] += stone_value;
        adj_[adj_id].remove(cell_id);
    }
}


void Position::unfill(U32 cell_id, Value stone_value) {
    List& adj = adj_[cell_id];
    for (U32 i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        value_[adj_id] -= stone_value;
        adj_[adj_id].readd(cell_id);
    }
}


void Position::make_move(const Move& move) {
    U32 cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
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
    fill(cell_id, move.stone_value_);
    
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = MASK(adj_id);
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
    
    U32 cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    U32 * stone_loc = stone_loc_[turn_];
    U32& n_stones = n_stones_[turn_];
    U32 i;
    for (i = n_stones; i && stone_number < stones[i-1]; i--) {
        stones[i] = stones[i-1];
        stone_loc[stones[i]] = i;
    }
    stones[i] = stone_number;
    stone_masks_[turn_] |= (1 << (stone_number - 1));
    power_[turn_] = STONE_POWER[turn_][stone_masks_[turn_]];
    n_stones++;
   
    U64 cell_mask = MASK(cell_id);
    open_mask_ |= cell_mask;
    open_.readd(cell_id); 
    
    unfill(cell_id, move.stone_value_);
    if (controls_ & cell_mask)
        n_controls_[BLUE]++;
    else
        n_controls_[RED]++;
    for (U32 p = 0; p < 2; p++) {
        if (dead_[p] & cell_mask) {
            n_dead_[p]++;
            if (stale_[p] & cell_mask)
                n_stale_[p]++;
        }
    }
    List& adj = adj_[cell_id];
    for (i = 0; i < adj.len_; i++) {
        U32 adj_id = adj.val_[i];
        U64 adj_mask = MASK(adj_id);
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
        if ((dead_[RED] & open_mask)) {
            if (! is_dead(open_id, RED)) {
                n_dead_[RED]--;
                dead_[RED] ^= open_mask;
            }
        } else if ((dead_[BLUE] & open_mask)) {
            if (! is_dead(open_id, BLUE)) {
                n_dead_[BLUE]--;
                dead_[BLUE] ^= open_mask;
            }
        }
    }

    U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
    U64 stale = open_mask_ & (stale_[RED] | stale_[BLUE]);
    U64 stale_it = stale;
    while (stale_it) {
        U64 lsb = LSB(stale_it);
        stale_it ^= lsb;
        if (lsb & stale) {
            bool still_stale = true;
            if (! (dead & lsb))
                still_stale = false;
            else {
                U64 adj_mask = ADJ_MASK[INDEX(lsb)];
                if ((dead & adj_mask) != (open_mask_ & adj_mask))
                    still_stale = false;
            }
            if (! still_stale) {
                if (stale_[RED] & lsb) {
                    stale_[RED] ^= lsb;
                    n_stale_[RED]--;
                } else {
                    stale_[BLUE] ^= lsb;
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


/*void Position::light_unmake_move(const Move& move) {
    turn_ ^= 1;
    m_ *= -1;
    
    U32 cell_id = move.cell_;
    Value stone_number = m_ * move.stone_value_;
    Value * stones = stones_[turn_];
    U32 * stone_loc = stone_loc_[turn_];
    U32& n_stones = n_stones_[turn_];
    U32 i;
    for (i = n_stones; i && stone_number < stones[i-1]; i--) {
        stones[i] = stones[i-1];
    }
    stones[i] = stone_number;
    n_stones++;
   
    open_.readd(cell_id); 
}*/


Move Position::get_random_move() {
    return Move(open_.val_[rand() % open_.len_],
                m_ * stones_[turn_][rand() % n_stones_[turn_]]);
}


U32 lottery(U32 * tickets, U32 n_tickets) {
    int winner = rand() % n_tickets;
    U32 i;
    for (i = 0; winner >= 0; i++)
        winner -= tickets[i];
    return i - 1;
}


Move Position::get_default_policy_move() {
    U32 tickets[N_CELLS];
    U32 total_tickets = 0;
    U32 i;
    U32 best_stale = NO_CELL;
    U32 cell_id;
    Value best_stale_val = 1000;
    for (i = 0; i < open_.len_; i++) {
        tickets[i] = 0;
        cell_id = open_.val_[i];
        U64 mask = (1LL << (U64)cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            Value stale_val = m_ * value_[cell_id];
            if (stale_val < best_stale_val) {
                best_stale = i;
                best_stale_val = stale_val;
            }
        } else {
            if ((dead_[1-turn_] & mask) || (n_stones_[turn_] > n_dead_[1-turn_])) {
                if ((bool)(controls_ & (1LL << U64(cell_id))) == (bool)turn_)
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
    
    U32 cell_index = lottery(tickets, total_tickets);
    cell_id = open_.val_[cell_index];

    Value reqs[MAX_DEGREE];
    U32 n_uncontrolled = 0;
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
    U32 n_stones = n_stones_[turn_];
    U32 stone_index;
    if (cell_index == best_stale)
        stone_index = 0;
    else {
        total_tickets = 0;
        for (i = n_stale_[1 - turn_]; i < n_stones; i++) {
            U32 changes;
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
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
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


U32 ways[N_STONES][MAX_DEGREE + 1][2 * N_STONES + 1][2 * MAX_RESULT + 1];
U32 combos[MAX_DEGREE + 1];
Value stone_list[2 * N_STONES];


void Position::find_n_ways(U32 stone_index) const {
    U32 i, j, n_stones = 0;
    Value k;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            if (! (i == turn_ && j == stone_index))
                stone_list[n_stones++] = stones_[i][j] * parity[i];
    U32 (*dp)[2*N_STONES+1][2*MAX_RESULT+1] = ways[stone_index];
    for (i = 0; i <= n_stones; i++)
        dp[0][i][MAX_RESULT] = 1;
    for (i = 1; i <= MAX_DEGREE && i <= n_stones; i++) {
        for (j = 1; j <= n_stones; j++) {
            Value val = stone_list[j-1];
            for (k = 0; k <= 2 * MAX_RESULT; k++) {
                U32& ans = dp[i][j][k];
                ans = dp[i][j-1][k];
                if (k - val >= 0 && k - val <= 2 * MAX_RESULT)
                    ans += dp[i-1][j-1][k-val];
            }
        }
    }
}


Real Position::calculate_win_prob(Value alpha, U32 stone_index) const {
    U32 total_stones = n_stones_[RED] + n_stones_[BLUE];
    Real prob_sum = 0.0;
    U32 (*dp)[2*N_STONES+1][2*MAX_RESULT+1] = ways[stone_index];
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        U32 n_adj = adj_[cell_id].len_;
        Value val = value_[cell_id];
        Value diff = alpha - val;
        U32 n_acceptable_outcomes = 0;
        for (Value k = diff + MAX_RESULT; k <= 2 * MAX_RESULT; k++)
            n_acceptable_outcomes += dp[n_adj][total_stones][k];
        prob_sum += (Real)n_acceptable_outcomes / (Real)combos[n_adj];
    }
    return prob_sum / open_.len_;
}


std::pair<Real, Move> Position::get_best_alpha_move(Value alpha) {
    U32 i, j;
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    U32 total_stones = n_stones_[RED] + n_stones_[BLUE] - 1;
    
    memset(ways, 0, sizeof(ways));
    combos[0] = 1;
    for (i = 1; i <= MAX_DEGREE && i <= total_stones; i++)
        combos[i] = combos[i-1] * (total_stones + 1 - i) / i;
    for (i = 0; i < n_stones; i++)
        find_n_ways(i);
    
    Move best_move;
    Real best_prob = -1000.0;
    for (i = 0; i < open_.len_; i++) {
        for (j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            make_move(move);
            Real prob = calculate_win_prob(alpha, j);
            unmake_move(move);
            if (m_ * prob > best_prob) {
                best_move = move;
                best_prob = m_ * prob;
            }
        }
    }
    return std::pair<Real, Move>(m_ * best_prob, best_move);
}


std::pair<Real, Move> Position::get_best_move() {
    Value alpha = 0;
    Value best_alpha = 400;
    Value direction = 0;
    Move best_move;
    Real best_prob = 0.0;
    while (alpha >= MIN_RESULT && alpha <= MAX_RESULT) {
        std::pair<Real, Move> result = get_best_alpha_move(alpha);
        Real red_win_prob = result.first;
        fprintf(stderr, "Alpha = %d: %.5lf\n", alpha, red_win_prob);
        if (best_alpha == 400 ||
                (turn_ == RED && red_win_prob >= CHOOSE_TARGET_THRESH[RED]) ||
                (turn_ == BLUE && red_win_prob < CHOOSE_TARGET_THRESH[BLUE])) {
            best_move = result.second;
            best_prob = red_win_prob;
            best_alpha = alpha;
        }
        Value change;
        if (red_win_prob >= CHOOSE_TARGET_THRESH[turn_])
            change = 1;
        else 
            change = -1;
        if (direction == 0)
            direction = change;
        else if (change != direction)
            break;
        alpha += change;
    }
    return std::pair<Real, Move>(best_prob, best_move);
}


Real Position::calculate_expectation() const {
    if (open_.len_ == 1)
        return (Real)(value_[open_.val_[0]]);
    Real sum = 0.0;
    Real stone_sum = 0.0;
    U32 i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += parity[i] * stones_[i][j];
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        sum += value_[cell_id] + adj_[cell_id].len_ * stone_sum / (open_.len_ - 1);
    }
    return sum / open_.len_;
}


std::pair<Real, Move> Position::search_expectation(U32 depth, Real a, Real b) {
    if (depth == 0)
        return std::pair<Real, Move>(calculate_expectation(), Move());
    Real best_val = -1000.0;
    Move best_move;
    U32 n_stones = n_stones_[turn_];
    Value * stones = stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        for (U32 j = 0; j < n_stones; j++) {
            Move move(cell_id, stones[j] * m_);
            make_move(move);
            std::pair<Real, Move> result = search_expectation(depth - 1, b, a);
            unmake_move(move);
            Real val = result.first * m_;
            if (val > best_val) {
                best_move = move;
                best_val = val;
                if (val > a) {
                    a = val;
                    if (a >= -b)
                        return std::pair<Real, Move>(m_ * best_val, best_move);
                }
            }
        }
    }
    return std::pair<Real, Move>(m_ * best_val, best_move);
}

        
Real Position::get_control_heuristic() const {
    return (Real)n_controls_[RED] / (Real)open_.len_;
}


void Position::get_moves_with_heuristic(std::vector<std::pair<Real, Move>>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
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
    U32 n_stones = n_stones_[turn_];
    for (U32 i = 0; i < open_.len_; i++) {
        for (U32 j = 0; j < n_stones; j++) {
            Value stone_value = m_ * stones[j];
            Move move(open_.val_[i], stone_value);
            moves.push_back(move);
        }
    }
}


bool Position::is_reasonable_move(U32 cell_index, U32 stone_index) {
    U32 cell_id = open_.val_[cell_index];
    Value stone_number = stones_[turn_][stone_index];
    //Value stone_value = m_ * stone_number;
    
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    
    U64 mask = MASK(cell_id);
    if (stale_[1 - turn_] & mask)
        //TODO: allow only the worst stale cell to be filled
        return stone_index == 0;
    if (stale_[turn_] & mask)
        return stone_index == 0 && dead_[turn_] == open_mask_;

    U32 op = 1 - turn_;
    if (stone_index < n_stale_[op])
        return false;
    
    U32 n_stones = n_stones_[turn_];
    if (n_stones == n_dead_[1 - turn_] && ! (dead_[1 - turn_] & mask))
        return false;

    if (adj_[cell_id].len_ == 1) {
        Value cell_value = value_[cell_id] * m_;
        U32 adj_id = adj_[cell_id].val_[0];
        if (adj_[adj_id].len_ == 1) {
            Value adj_value = value_[adj_id] * m_;
            if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                if (stone_index == 0)
                    return true;
                else {
                    Value prev_stone_number;
                    if (stone_index == 0)
                        prev_stone_number = 0;
                    else
                        prev_stone_number = stones_[turn_][stone_index - 1];
                    Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
                    return stone_number >= extra_req && prev_stone_number < extra_req;
                }       
            }
        } else {
            return ! (we_control & mask);
        }
    }
    return true;
}
 

void Position::get_all_reasonable_moves(std::vector<Move>& moves) {
    Value * stones = stones_[turn_];
    U32 n_stones = n_stones_[turn_];
    Value stone_value;
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
                            stone_value = m_ * stones[0];
                            moves.push_back(Move(cell_id, stone_value));
                            Value control_req = m_ * (alpha_ - OFFSET[turn_]);
                            if (adj_value + stones[0] < control_req) {
                                Value extra = control_req - adj_value;
                                for (j = 1; j < n_stones && stones[j] < extra; j++) {}
                                if (j < n_stones) {
                                    stone_value = m_ * stones[j];
                                    moves.push_back(Move(cell_id, stone_value));
                                }
                            }       
                        }
                    } else {
                        if (we_control & mask)
                            done = true;
                    }
                }
                if (! done) {
                    for (j = n_stale_[1 - turn_]; j < n_stones; j++) {
                        stone_value = m_ * stones[j];
                        moves.push_back(Move(cell_id, stone_value));
                    }
                }
            }
        }
    }
    if (best_stale != NO_CELL) {
        stone_value = stones[0] * m_;
        moves.push_back(Move(best_stale, stone_value));
    }
}


void Position::get_all_reasonable_moves_with_stone(U32 stone_index, std::vector<Move>& moves) {
    Value stone_number = stones_[turn_][stone_index];
    Value prev_stone_number;
    if (stone_index == 0)
        prev_stone_number = 0;
    else
        prev_stone_number = stones_[turn_][stone_index - 1];
    Value stone_value = m_ * stone_number;
    Value best_stale_val = 1000;
    U32 best_stale = NO_CELL;
    U32 n_stones = n_stones_[turn_];
    U32 op = 1 - turn_;
    U64 we_control = controls_;
    if (turn_ == RED)
        we_control = ~we_control;
    for (U32 i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        Value cell_value = value_[cell_id] * m_;
        U64 mask = MASK(cell_id);
        if ((stale_[RED] & mask) || (stale_[BLUE] & mask)) {
            if (stone_index == 0 && cell_value < best_stale_val) {
                best_stale = cell_id;
                best_stale_val = cell_value;
            }
        } else {
            if ((dead_[op] & mask) || (n_stones > n_dead_[op])) {
                bool done = false;
                if (adj_[cell_id].len_ == 1) {
                    U32 adj_id = adj_[cell_id].val_[0];
                    if (adj_[adj_id].len_ == 1) {
                        done = true;
                        Value adj_value = value_[adj_id] * m_;
                        if (cell_value < adj_value || (cell_value == adj_value && cell_id < adj_id)) {
                            if (stone_index == 0)
                                moves.push_back(Move(cell_id, stone_value));
                            else {
                                Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
                                if (stone_number >= extra_req && prev_stone_number < extra_req)
                                    moves.push_back(Move(cell_id, stone_value));
                            }       
                        }
                    } else {
                        if (we_control & mask)
                            done = true;
                    }
                }
                if (! done && stone_index >= n_stale_[op])
                    moves.push_back(Move(cell_id, stone_value));
            }
        }
    }
    if (best_stale != NO_CELL)
        moves.push_back(Move(best_stale, stone_value));
}


void Position::get_untried_move(U32& cell_index, U32& stone_index, AMAFTable& amaf) {
    U32 n_stones = n_stones_[turn_];
    stone_index = n_stones;
    while (stone_index > 0) {
        stone_index--;
        std::vector<Move> moves;
        get_all_reasonable_moves_with_stone(stone_index, moves);
        for (U32 i = 0; i < moves.size(); i++) {
            cell_index = open_.loc_[moves[i].cell_];
            if (! amaf.is_tried(cell_index, stone_index)) {
                amaf.set_tried(cell_index, stone_index);
                return;
            }
        }
    }
    cell_index = NO_CELL;
}


struct MoveInfo {
    Move move_;
    int swing_;

    MoveInfo(const Move& move, int swing): move_(move), swing_(swing) {}

    bool operator < (const MoveInfo& other) const {
        if (move_.stone_value_ != other.move_.stone_value_)
            return abs(move_.stone_value_) > abs(other.move_.stone_value_);
        if (swing_ != other.swing_)
            return swing_ > other.swing_;
        return move_.cell_ < other.move_.cell_;
    }
};


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
    std::vector<Move> moves;
    get_all_reasonable_moves(moves);
    std::vector<MoveInfo> move_infos;
    U32 i, j;
    U32 op = 1 - turn_;
    //U64 opp_controls = controls_;
    //if (turn_ == BLUE)
    //    opp_controls = ~opp_controls;
    for (i = 0; i < moves.size(); i++) {
        Move move = moves[i];
        U32 cell_id = move.cell_;
        int swing = adj_[cell_id].len_;
        //U64 mask = MASK(cell_id);
        //int swing = 0;
        /*if (dead_[turn_] & mask)
            swing = -1;
        else if (dead_[op] & mask)
            swing = 1;
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            U32 adj_id = adj.val_[j];
            U64 adj_mask = MASK(adj_id);
            if ((opp_controls & adj_mask) && ! (dead_[op] & adj_mask))
                swing++;
        }*/
        move_infos.push_back(MoveInfo(move, swing));
    }
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
    std::vector<Move> moves;
    get_all_reasonable_moves(moves);
    
    std::vector<MoveInfo> move_infos;
    U32 i, j;
    U32 op = 1 - turn_;
    U64 opp_controls = controls_;
    if (turn_ == BLUE)
        opp_controls = ~opp_controls;
    Move move;
    for (i = 0; i < moves.size(); i++) {
        move = moves[i];
        U32 cell_id = move.cell_;
        U64 mask = MASK(cell_id);
        int swing = 0;
        if (dead_[turn_] & mask)
            swing = -1;
        else if (dead_[op] & mask)
            swing = 1;
        List& adj = adj_[cell_id];
        for (j = 0; j < adj.len_; j++) {
            U32 adj_id = adj.val_[j];
            U64 adj_mask = MASK(adj_id);
            if ((opp_controls & adj_mask) && ! (dead_[op] & adj_mask))
                swing++;
        }
        move_infos.push_back(MoveInfo(move, swing));
    }
    std::sort(move_infos.begin(), move_infos.end());
   
#ifdef DEBUG_
    fprintf(stderr, "Solving: alpha = %d\n", alpha_);
#endif
    std::vector<Move> winning_moves;
    for (i = 0; i < move_infos.size(); i++) {
        move = move_infos[i].move_;
#ifdef DEBUG_
        char move_str[10];
        move.to_str(move_str);
        fprintf(stderr, "%s\n", move_str);
#endif
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
    return std::pair<U32, Move>(1 - turn_, moves[0]);
}


bool Position::is_dead(U32 cell_id, U32 p) {
    Value * other_power = power_[1 - p];
    Value value = value_[cell_id];
    U32 n_adj = adj_[cell_id].len_;
    Value worst_val = parity[p] * (value + other_power[n_adj]);
    return (worst_val >= parity[p] * (alpha_ - OFFSET[p]));
}


bool Position::all_adj_dead(U32 cell_id) {
    List& adj = adj_[cell_id];
    U32 n_adj = adj.len_;
    for (U32 i = 0; i < n_adj; i++) {
        U64 adj_mask = MASK(adj.val_[i]);
        if (! (dead_[RED] & adj_mask) && ! (dead_[BLUE] & adj_mask))
            return false;
    }
    return true;
}


void Position::find_stale_cells() {
    U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
    U64 dead_it = dead;
    U64 stale = open_mask_ & (stale_[RED] | stale_[BLUE]);
    while (dead_it) {
        U64 lsb = LSB(dead_it);
        dead_it ^= lsb;
        if (! (stale & lsb)) {
            U64 adj_mask = ADJ_MASK[INDEX(lsb)];
            if ((dead & adj_mask) == (open_mask_ & adj_mask)) {
                if (dead_[RED] & lsb) {
                    stale_[RED] |= lsb;
                    n_stale_[RED]++;
                } else {
                    stale_[BLUE] |= lsb;
                    n_stale_[BLUE]++;
                }
            }
        }
    }
}


bool Position::is_winning(U32 p) const {
    return n_dead_[p] > n_stones_[1 - p];
}


void Position::set_alpha(Value alpha) {
    alpha_ = alpha;
    n_controls_[0] = n_controls_[1] = 0;
    n_dead_[0] = n_dead_[1] = 0;
    controls_ = 0;
    dead_[0] = dead_[1] = 0;
    stale_[RED] = stale_[BLUE] = 0;
    n_stale_[RED] = n_stale_[BLUE] = 0;
    U32 i, cell_id, n_adj;
    U64 mask;
    for (i = 0; i < open_.len_; i++) {
        cell_id = open_.val_[i];
        Value value = value_[cell_id];
        n_adj = adj_[cell_id].len_;
        mask = MASK(cell_id);
        if (value >= alpha_) {
            n_controls_[RED]++;
            if (value + power_[BLUE][n_adj] >= alpha_) {
                n_dead_[RED]++;
                dead_[RED] |= mask;
            }
        } else {
            controls_ |= mask;
            n_controls_[BLUE]++;
            if (value + power_[RED][n_adj] < alpha_) {
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
    U32 i, j, p;
    snapshot_.turn_ = turn_;
    snapshot_.open_mask_ = open_mask_;
    snapshot_.n_open_ = open_.len_;
    snapshot_.controls_ = controls_;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
        snapshot_.open_[i] = cell_id;
        snapshot_.value_[cell_id] = value_[cell_id];
        List& adj = adj_[cell_id];
        U32 n_adj = adj.len_;
        snapshot_.n_adj_[cell_id] = n_adj;
        for (j = 0; j < n_adj; j++)
            snapshot_.adj_[cell_id][j] = adj.val_[j];
    }

    for (p = 0; p < 2; p++) {
        U32 n_stones = n_stones_[p];
        snapshot_.n_stones_[p] = n_stones;
        snapshot_.stone_masks_[p] = stone_masks_[p];
        snapshot_.power_[p] = power_[p];
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
    U32 i, j, p;
    turn_ = snapshot_.turn_;
    m_ = parity[turn_];
    open_.len_ = snapshot_.n_open_;
    open_mask_ = snapshot_.open_mask_;
    controls_ = snapshot_.controls_;
    for (i = 0; i < snapshot_.n_open_; i++) {
        U32 cell_id = snapshot_.open_[i];
        value_[cell_id] = snapshot_.value_[cell_id];
        open_.val_[i] = cell_id;
        open_.loc_[cell_id] = i;
        U32 n_adj = snapshot_.n_adj_[cell_id];
        List& adj = adj_[cell_id];
        adj.len_ = n_adj; 
        for (j = 0; j < n_adj; j++) {
            U32 adj_id = snapshot_.adj_[cell_id][j];
            adj.val_[j] = adj_id;
            adj.loc_[adj_id] = j;
        }
    }
    
    for (p = 0; p < 2; p++) {
        U32 n_stones = snapshot_.n_stones_[p];
        n_stones_[p] = n_stones;
        stone_masks_[p] = snapshot_.stone_masks_[p];
        power_[p] = snapshot_.power_[p];
        Value * stones = stones_[p];
        U32 * stone_loc = stone_loc_[p];
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


std::pair<U32, U32> Position::get_cell_and_stone_indices(const Move& move) {
    return std::pair<U32, U32>(open_.loc_[move.cell_],
                                     stone_loc_[turn_][m_ * move.stone_value_]);
}


void Position::print(FILE * f) {
    U32 i, j;
    for (i = 0; i < open_.len_; i++) {
        U32 cell_id = open_.val_[i];
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

