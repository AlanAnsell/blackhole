#include "globals.h"
#include "Position.h"

int main(int argc, char** argv) {
	init();
	// Read in board info
    time_started = get_time();
	
    char cell_str[5];
    CellID blocked[N_BLOCKED_CELLS];
    for (size_t i = 0; i < N_BLOCKED_CELLS; i++) {
        assert(scanf("%s", cell_str));
        blocked[i] = cell_name_to_id(cell_str);
    }
    
    Position pos(blocked);
    char move_str[100];
    get_move(move_str);
	Value endgame_value;
    if (strcmp(move_str, "Start"))
        pos.make_move(Move(move_str, pos.turn_));
	while (true) {
		fprintf(stderr, "Time left: %.2f seconds\n", (float)time_left / 1e6);
		endgame_value = pos.dead_endgame_value();
        if (endgame_value != NO_RESULT)
            fprintf(stderr, "Solved endgame: %d\n", endgame_value);
        Move move = pos.get_expectation_maximising_move_with_endgame_solve();
		pos.make_move(move);
		endgame_value = pos.dead_endgame_value();
        if (endgame_value != NO_RESULT)
            fprintf(stderr, "Solved endgame: %d\n", endgame_value);
        move.to_str(move_str);
		send_move(move_str);
		get_move(move_str);
        pos.make_move(Move(move_str, pos.turn_));
	}
	return 0;
}





