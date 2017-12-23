#include "Position.h"


void Cell::init(CellID id, Cell * cells) {
    id_ =  id;
    degree_ = N_ADJ[id_];
    value_ = 0;
    for (size_t i = 0; i < degree_; i++) {
        adj_[i] = cells + ADJ[id_][i];
        loc_[ADJ[id_][i]] = i;
    }
}


void Cell::remove_neighbour(CellID neighbour_id, Value stone_value) {
    value_ += stone_value;
    degree_--;
    if (degree_) {
        size_t neighbour_loc = loc_[neighbour_id];
        Cell * last_adj = adj_[degree_];
        loc_[last_adj->id_] = neighbour_loc;
        adj_[neighbour_loc] = last_adj;
    }
}
	

void Cell::add_neighbour(CellID neighbour_id, Cell * neighbour, Value stone_value) {
    value_ -= stone_value;
    size_t old_loc = loc_[neighbour_id];
    Cell * new_neighbour = adj_[old_loc];
    loc_[new_neighbour->id_] = degree_;
    adj_[degree_] = new_neighbour;
    adj_[old_loc] = neighbour;
    degree_++;
}


void Cell::fill(Value stone_value) {
    for (size_t i = 0; i < degree_; i++)
        adj_[i]->remove_neighbour(id_, stone_value);
}


void Cell::unfill(Value stone_value) {
    for (size_t i = 0; i < degree_; i++)
        adj_[i]->add_neighbour(id_, this, stone_value);
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


Position::Position(CellID blocked[N_BLOCKED_CELLS]): n_open_(0), turn_(RED) {
    size_t i;
    bool is_blocked[N_CELLS] = {false};
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        is_blocked[blocked[i]] = true;
    for (i = 0; i < N_CELLS; i++) {
        cells_[i].init(i, cells_);
        if (! is_blocked[i]) {
            open_cells_[n_open_] = i;
            open_loc_[i] = n_open_;
            n_open_++;
        }
    }
    for (i = 0; i < N_BLOCKED_CELLS; i++)
        cells_[blocked[i]].fill(0);

    for (size_t p = 0; p < 2; p++) {
        n_stones_[p] = N_STONES;
        for (i = 0; i < N_STONES; i++)
            stones_[p][i] = i + 1;
    }
}

void Position::make_move(const Move& move) {
    Value stone_number = parity[turn_] * move.stone_value_;
    Value * stones = stones_[turn_];
    size_t& n_stones = n_stones_[turn_];
    n_stones--;
    size_t i = 0;
    //printf("Deleting stone\n");
    //fflush(stdout);
    while (stones[i] != stone_number)
        i++;
    for (; i < n_stones; i++)
        stones[i] = stones[i+1];
    turn_ ^= 1;
    
    //printf("Deleting cell\n");
    //fflush(stdout);
    n_open_--;
    size_t move_loc = open_loc_[move.cell_];
    CellID last_cell = open_cells_[n_open_];
    open_cells_[move_loc] = last_cell;
    open_loc_[last_cell] = move_loc;
    
    //printf("Filling cell\n");
    //fflush(stdout);
    cells_[move.cell_].fill(move.stone_value_);
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
    
    size_t old_loc = open_loc_[move.cell_];
    size_t new_cell = open_cells_[old_loc];
    open_loc_[new_cell] = n_open_;
    open_cells_[n_open_] = new_cell;
    open_cells_[old_loc] = move.cell_;
    n_open_++;
    
    cells_[move.cell_].unfill(move.stone_value_);
}


Move Position::get_random_move() {
    return Move(open_cells_[rand() % n_open_],
                parity[turn_] * stones_[turn_][rand() % n_stones_[turn_]]);
}


Move Position::get_expectation_maximising_move() {
    Move best_move(0, 0);
    Real best_expectation = -1000.0;
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < n_open_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_cells_[i], stone_value);
            make_move(move);
            Real expectation = calculate_expectation();
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
    if (n_open_ == 1) {
        //print();
        return (Real)(cells_[open_cells_[0]].value_);
    }
    Real sum = 0.0;
    Value * stones = stones_[turn_];
    size_t n_stones = n_stones_[turn_];
    for (size_t i = 0; i < n_open_; i++) {
        for (size_t j = 0; j < n_stones; j++) {
            Value stone_value = parity[turn_] * stones[j];
            Move move(open_cells_[i], stone_value);
            make_move(move);
            sum += search_expectation();
            unmake_move(move);
        }
    }
    //printf("Exiting search_expectation\n");
    //fflush(stdout);
    return sum / (n_open_ * n_stones);
}


Real Position::calculate_expectation() {
    if (n_open_ == 1)
        return (Real)(cells_[open_cells_[0]].value_);
    Real sum = 0.0;
    Real stone_sum = 0.0;
    size_t i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < n_stones_[i]; j++)
            stone_sum += parity[i] * stones_[i][j];
    for (i = 0; i < n_open_; i++) {
        Cell& cell = cells_[open_cells_[i]];
        sum += cell.value_ + cell.degree_ * stone_sum / (n_open_ - 1);
    }
    return sum / n_open_;
}
    

void Position::print() {
    size_t i, j;
    for (i = 0; i < n_open_; i++) {
        char cell_str[10];
        cell_id_to_name(open_cells_[i], cell_str);
        printf("%s (%d):", cell_str, cells_[open_cells_[i]].value_);
        for (j = 0; j < cells_[open_cells_[i]].degree_; j++) {
            cell_id_to_name(cells_[open_cells_[i]].adj_[j]->id_, cell_str);
            printf(" %s", cell_str);
        }
        printf("\n");
    }
    printf("Red:");
    for (i = 0; i < n_stones_[RED]; i++)
        printf(" %d", stones_[RED][i] * parity[RED]);
    printf("\nBlue:");
    for (i = 0; i < n_stones_[BLUE]; i++)
        printf(" %d", stones_[BLUE][i] * parity[BLUE]);
    printf("\n");
}


/*size_t n_dead(size_t player, Value alpha) {
    Value * stones = stones_[player];
    size_t n_stones = n_stones_[player];
    Value cuml[N_STONES+1];
    cuml[0] = 0;
    size_t i;
    for (i = 0; i < n_stones; i++)
        cuml[i+1] = cuml[i] + stones[i];*/

