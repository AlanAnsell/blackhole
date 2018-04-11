#ifndef POSITION_H_
#define POSITION_H_

#include "globals.h"
#include "AMAF.h"

class List {
    /*
     * Stores a subset of the values 0...max_len-1 which supports constant time insertion and deletion
     * and maintains its order when an item is deleted and re-added assuming that items
     * are re-added in the reverse order in which they were deleted.
     *
     * It's possible that a bitmask would work better in hindsight.
     */
public:
    U32 * val_; // list of cell indices
    U32 * loc_; // loc_[i] = j <=> val_[j] = i
    U32 len_; // size of list

    List(): val_(NULL), loc_(NULL), len_(0) {}

    void init(U32 max_len, U32 range) {
        len_ = 0;
        val_ = new U32[max_len];
        loc_ = new U32[range];
    }

    ~List() {
        delete [] val_;
        delete [] loc_;
    }

    inline void add(U32 val) {
        // can only be called before remove or readd
        loc_[val] = len_;
        val_[len_] = val;
        len_++;
    }

    inline void remove(U32 val) {
        len_--;
        if (len_) {
            U32 end_val = val_[len_];
            U32 val_loc = loc_[val];
            loc_[end_val] = val_loc;
            val_[val_loc] = end_val;
        }
    }

    inline void readd(U32 val) {
        U32 old_loc = loc_[val];
        U32 new_val = val_[old_loc];
        loc_[new_val] = len_;
        val_[len_] = new_val;
        val_[old_loc] = val;
        len_++;
    }

};


typedef U32 Move;

#define CREATE_MOVE(cell, stone) ((cell) | ((stone) << 6))
#define GET_CELL(move) ((move) & 0x3F)
#define GET_STONE_NUMBER(move) ((move) >> 6)

Move move_from_str(char * move_str);

void move_to_str(Move move, char * str);

struct PositionHistory {
    /*
     * All the data stored in a Position which is easier to cache and restore than to
     * recalculate when unmaking a move.
     */
    
    // These variables are changed by make_move and restore_history
    U64 open_mask_;
    Value value_[N_CELLS];
	Value stones_[2][N_STONES];
	U32 stone_masks_[2];
    U32 stone_loc_[2][N_STONES+1];
    U32 n_stones_[2];
	Value * power_[2];
    
    // These variables are changed by make_move, restore_history and set_alpha
    U64 controls_;
    U32 n_controls_[2];
    U64 dead_[2];
    U32 n_dead_[2];
    U64 stale_[2];
    U32 n_stale_[2];
    U64 either_stale_;
    U32 effective_adj_[N_CELLS];
    U32 worst_stale_;
};

class HashInfo;

#define HASH_COEFFICIENT 6787142722019916807ULL
#define HASH_SIZE_BASE 14
#define HASH_L (64 - HASH_SIZE_BASE)
#define HASH_SIZE (1 << HASH_SIZE_BASE)
#define BUCKET_SIZE 10
#define NO_HASH_ENTRY 2
// High Speed Hashing for Integers and Strings
// Mikkel Thorup July 15, 2014
#define HASH_BITMASKS(x, y) ((HASH_COEFFICIENT * (x) * (y)) >> HASH_L)


class HashRecord {
public:
    U64 red_;
    U64 blue_;
    U32 result_;
};


class HashTable {
public:
    U32 n_records_[HASH_SIZE];
    HashRecord table_[HASH_SIZE][BUCKET_SIZE];
    
    void init() {
        for (U32 i = 0; i < HASH_SIZE; i++)
            n_records_[i] = 0;
    }

    void add(U64 red, U64 blue, U32 result) {
        U64 hash = HASH_BITMASKS(red, blue);
        //fprintf(stderr, "Hash: %llu\n", hash);
        U32& n_records = n_records_[hash];
        if (n_records < BUCKET_SIZE - 1) {
            HashRecord& hash_record = table_[hash][n_records];
            hash_record.red_ = red;
            hash_record.blue_ = blue;
            hash_record.result_ = result;
            n_records++;
        }
#ifdef DEBUG_
        else {
            //fprintf(stderr, "Hash bucket %llu overflowed\n", hash);
        }
#endif
    }

    U32 find(U64 red, U64 blue) {
        U64 hash = HASH_BITMASKS(red, blue);
        U32 n_records = n_records_[hash];
        HashRecord * records = table_[hash];
        for (U32 i = 0; i < n_records; i++) {
            HashRecord& record = records[i];
            if (record.red_ == red && record.blue_ == blue)
                return record.result_;
        }
        return NO_HASH_ENTRY;
    }

};

class MoveInfo;

class Position {
public:
    //***These variables are changed by make_move and regress_move***
    // list of currently empty cells
    List open_;
    // list of empty neighbours for each cell
    List adj_[N_CELLS];
    U32 turn_;
    U32 n_moves_made_;
    // -1 if blue to move, 1 if red
    int m_;
    
    //***These variables are changed by make_move and restore_history***
    // bitmask of empty cells
    U64 open_mask_;
    // sum of values of stones adjacent to each cell
    Value value_[N_CELLS];
    // stones not yet played by each player
	Value stones_[2][N_STONES];
    // bitmask of stones not yet played by each player
	U32 stone_masks_[2];
    // maps stone values to their index in stones_
    U32 stone_loc_[2][N_STONES+1];
    // number of stones left for each player
    U32 n_stones_[2];
    // power_[p][i] = sum of player p's most valuable i stones
	Value * power_[2];
    
    
    //***These variables are changed by make_move, restore_history and set_alpha***
    // bitmasks for control, death and staleness
    U64 controls_;
    U32 n_controls_[2];
    U64 dead_[2];
    U32 n_dead_[2];
    U64 either_stale_;
    U64 stale_[2];
    U32 n_stale_[2];
    // number of non-stale neighbours for each cell
    U32 effective_adj_[N_CELLS];
    // index of worst stale cell (i.e. the only one the current player could reasonably play in)
    U32 worst_stale_;
    
    //***This variable is changed only by set_alpha***
    Value alpha_;

    // History variables
    Move moves_made_[N_CELLS];
    PositionHistory history_[N_CELLS];


	Position(U32 blocked[N_BLOCKED_CELLS], Value alpha);
   
    inline void fill(U32 cell_id, Value stone_value) {
        List& adj = adj_[cell_id];
        for (U32 i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            value_[adj_id] += stone_value;
            adj_[adj_id].remove(cell_id);
        }
    }


    inline void unfill(U32 cell_id, Value stone_value) {
        List& adj = adj_[cell_id];
        for (U32 i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            value_[adj_id] -= stone_value;
            adj_[adj_id].readd(cell_id);
        }
    }

    bool is_legal(Move move);

    void save_history();

    void restore_history();

    void make_move(Move move, bool save);
    
	void regress_board();

    void rewind_to(U32 n_moves_made);

    inline void unmake_move() {
        rewind_to(n_moves_made_ - 1);
    }

    Move get_default_policy_move(U64 valid, U64 duo);

    Move get_expectation_maximising_move();

    Real calculate_expectation() const;		

    Move get_best_winning_move();

    void get_all_reasonable_moves(Move * moves, U32& N);

    void get_all_reasonable_moves_with_stone(U32 stone_index, U64 valid, U64 duo, Move * moves, U32& N);

    bool add_solver_move(Move move, MoveInfo * moves, U32& N);

    bool get_solver_moves(MoveInfo * moves, U32& N, bool ignore_killer);

    void generate_untried_moves(int& stone_index, U64 valid, U64 duo, std::vector<Move>& untried_moves);

    U32 solve(long long end_time, U32& counter, U32& hash_hits,
            HashTable& table, HashInfo& hash_info, U64 stone_masks[2]);

    std::pair<U32, Move> get_optimal_move(
            long long end_time, U32& counter, U32& hash_hits, bool break_ties, HashTable& table);

    bool all_adj_dead(U32 cell_id);

    void set_alpha(Value alpha);

    void print(FILE * f) const;
        
    // updates variables relating to staleness
    inline void find_stale_cells() {
        U64 dead = open_mask_ & (dead_[RED] | dead_[BLUE]);
        U64 dead_it = dead;
        // a cell can only be stale if it is dead
        while (dead_it) {
            U64 lsb = LSB(dead_it);
            dead_it ^= lsb;
            // check the cell has not already been identified as stale
            if (! (either_stale_ & lsb)) {
                U64 adj_mask = ADJ_MASK[INDEX(lsb)] & open_mask_;
                // cell is stale if all neighbours are dead
                if ((dead & adj_mask) == adj_mask) {
                    if (dead_[RED] & lsb) {
                        stale_[RED] |= lsb;
                        n_stale_[RED]++;
                    } else {
                        stale_[BLUE] |= lsb;
                        n_stale_[BLUE]++;
                    }
                    // decrement effective_adj_ for all neighbours
                    while (adj_mask) {
                        lsb = LSB(adj_mask);
                        adj_mask ^= lsb;
                        effective_adj_[INDEX(lsb)]--;
                    }
                }
            }
        }
        either_stale_ = stale_[RED] | stale_[BLUE];
    }

    inline void find_worst_stale() {
        worst_stale_ = NO_CELL;
        Value worst_val = VALUE_INF;
        U64 stale = either_stale_ & open_mask_;
        while (stale) {
            U64 lsb = LSB(stale);
            stale ^= lsb;
            U32 cell_id = INDEX(lsb);
            Value val = m_ * value_[cell_id];
            if (val < worst_val) {
                worst_stale_ = cell_id;
                worst_val = val;
            }
        }
    }

    inline void find_effective_adj() {
        for (U32 i = 0; i < open_.len_; i++) {
            U32 cell_id = open_.val_[i]; 
            U32& effective_adj = effective_adj_[cell_id];
            effective_adj = 0;
            U64 non_stale_adj_mask = ADJ_MASK[cell_id] & open_mask_ & ~either_stale_;
            while (non_stale_adj_mask) {
                non_stale_adj_mask ^= LSB(non_stale_adj_mask);
                effective_adj++;
            }
        }
    }

    // checks whether a cell is dead if favour of player p, i.e. whether the cell's value
    // will definitely be worse than the target if it becomes the black hole.
    inline bool is_dead(U32 cell_id, U32 p) const {
        U32 op = 1 - p;
        Value * other_power = power_[op];
        Value value = value_[cell_id];
        // the maximum number of stones the opponent can place adjacent to this stone
        // in the future
        U32 n_adj = std::min(n_stones_[op], adj_[cell_id].len_);
        // the worst-case value of the cell at the end of the game
        Value worst_val = PARITY[p] * (value + other_power[n_adj]);
        // check if the worst-case value is still as good as or better than the target
        return (worst_val >= PARITY[p] * (alpha_ - OFFSET[p]));
    }

    inline bool is_winning(U32 p) const {
        return n_dead_[p] > n_stones_[1 - p];
    }

    // gets the nth non-stale cell in a list
    inline U32 get_non_stale_adj(const List& adj, U32 n) const {
        U32 count = 0;
        U32 i;
        for (i = 0; i < adj.len_; i++) {
            U32 adj_id = adj.val_[i];
            if (! (either_stale_ & MASK(adj_id))) {
                if (count == n)
                    return adj_id;
                else
                    count++;
            }
        }
#ifdef DEBUG_
        // n should always be less than the length of the list
        assert(false);
#endif
        return NO_CELL;
    }

    inline bool is_adj(const List& adj, U32 id2) const {
        for (U32 i = 0; i < adj.len_; i++)
            if (adj.val_[i] == id2)
                return true;
        return false;
    }

    // determines whether there is any reasonable move in a particular cell
    inline bool is_valid(U32 cell_id, bool& duo) const {
        duo = false;
        U64 mask = MASK(cell_id);
        
        // it's only reasonable to play in a stale cell if it's the least favourable one for you
        if (either_stale_ & mask)
            return cell_id == worst_stale_;
        
        // if all your remaining stones have to go in dead cells, then don't play in a non-dead
        // cell
        if (n_stones_[turn_] == n_dead_[1-turn_] &&
                ! (dead_[1-turn_] & mask))
            return false;

        // all other unreasonable moves we detect involve cells with two or fewer neighbours
        if (effective_adj_[cell_id] > 2)
            return true;

        const List& adj = adj_[cell_id];
        Value val = m_ * value_[cell_id];
        U32 adj_id;
        Value adj_val;
        if (effective_adj_[cell_id] == 1) {
            adj_id = get_non_stale_adj(adj, 0);
            adj_val = m_ * value_[adj_id];
            if (effective_adj_[adj_id] == 1) {
                // the cell is one of an isolated pair, so it is only reasonable to play in it
                // if it has a value not better than its neighbour's
                duo = true;
                return (val < adj_val) || (val == adj_val && cell_id < adj_id);
            }

            const List& adj_adj = adj_[adj_id];
            U64 we_control;
            if (turn_ == RED)
                we_control = ~controls_;
            else
                we_control = controls_;
            // there's never any point playing in a cell which we control which has a single neighbour
            // which itself has several neighbours because it is never worse to play in said
            // neighbour 
            if (we_control & mask)
                return false;
            
            if (effective_adj_[adj_id] == 2) {
                U32 adj_adj_id = get_non_stale_adj(adj_adj, 0);
                if (adj_adj_id == cell_id)
                    adj_adj_id = get_non_stale_adj(adj_adj, 1);
                if (effective_adj_[adj_adj_id] == 1) {
                    // if we have a line of 3 isolated cells, only play at the least favourable end
                    Value adj_adj_val = value_[adj_adj_id];
                    return (val < adj_adj_val) || (val == adj_adj_val && cell_id < adj_adj_id);
                }
            }
            return true;
        }
        if (effective_adj_[cell_id] == 2) {
            U32 adj1_id = get_non_stale_adj(adj, 0);
            if (effective_adj_[adj1_id] == 2) {
                U32 adj2_id = get_non_stale_adj(adj, 1);
                const List& adj2 = adj_[adj2_id];
                if (effective_adj_[adj2_id] == 2 && is_adj(adj2, adj1_id)) {
                    // if we have an isolated triangle shape, only play in the worst cell
                    // of the triangle
                    Value other_val = m_ * value_[adj1_id];
                    if (other_val < val || (other_val == val && adj1_id < cell_id))
                        return false;
                    other_val = m_ * value_[adj2_id];
                    if (other_val < val || (other_val == val && adj2_id < cell_id))
                        return false;
                }
            }
        }
        return true;
    }

    
    inline bool is_reasonable_move(U32 cell_id, U32 stone_index, U64 valid, U64 duo) const {
        U64 mask = MASK(cell_id);
        
        // first check that it is reasonable to play in the cell
        if (! (valid & mask))
            return false;
        
        // only put your least valuable stone in a stale cell
        if (cell_id == worst_stale_)
            return stone_index == 0;

        if (duo & mask) {
            // if you play in an isolated pair, only use your least valuable stone or
            // the least valuable stone needed to control the adjacent stone
            if (stone_index == 0)
                return true;
            Value stone_number = stones_[turn_][stone_index];
            U32 adj_id = get_non_stale_adj(adj_[cell_id], 0);
            Value adj_value = value_[adj_id] * m_;
            Value prev_stone_number = stones_[turn_][stone_index - 1];
            Value extra_req = m_ * (alpha_ - OFFSET[turn_]) - adj_value;
            return stone_number >= extra_req && prev_stone_number < extra_req;
        }

        return stone_index >= n_stale_[1 - turn_];
    }
 

    void get_validity_mask(U64& valid, U64& duo) const {
        // get all cells which are reasonable to play in
        valid = 0;
        duo = 0;
        for (U32 i = 0; i < open_.len_; i++) {
            U32 cell_id = open_.val_[i];
            bool is_duo = false;
            if (is_valid(cell_id, is_duo)) {
                valid |= MASK(cell_id);
                if (is_duo)
                    duo |= MASK(cell_id);
            }
        }
    }

};

class HashInfo {
public:
    U32 cell_index_[N_CELLS];
    U32 stone_shift_[2][N_STONES+1];

    HashInfo(const Position& pos) {
        U32 i;
        for (i = 0; i < pos.open_.len_; i++)
            cell_index_[pos.open_.val_[i]] = i + 1;
        for (U32 p = 0; p < 2; p++) {
            U32 n_stones = pos.n_stones_[p];
            const Value * stones = pos.stones_[p];
            for (U32 i = 0; i < n_stones; i++)
                stone_shift_[p][stones[i]] = 5 * i;
        }
    }

};


#endif
