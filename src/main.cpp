#include "globals.h"
#include "Position.h"
#include "MCTS.h"


int main(int argc, char** argv) {
    time_started = get_time();
	init();
	init_free_list();
	
    char cell_str[5];
    CellID blocked[N_BLOCKED_CELLS];
    for (size_t i = 0; i < N_BLOCKED_CELLS; i++) {
        assert(scanf("%s", cell_str));
        blocked[i] = cell_name_to_id(cell_str);
    }
    
    Position pos(blocked);
    char move_str[100];
    get_move(move_str);
	//Value endgame_value;
    if (strcmp(move_str, "Start"))
        pos.make_move(Move(move_str, pos.turn_));
	MCTNode * root = get_free();
    while (true) {
		fprintf(stderr, "Time left: %.2f seconds\n", (float)time_left / 1e6);
		//endgame_value = pos.dead_endgame_value();
        //if (endgame_value != NO_RESULT)
        //    fprintf(stderr, "Solved endgame: %d\n", endgame_value);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES - 1);
#endif
        //fprintf(stderr, "Creating MCTS tree\n");
        //fflush(stderr);
        root->init(NULL, Move(), false, 0.0);
        //fprintf(stderr, "Performing MCTS search\n");
        //fflush(stderr);
        Move move = root->get_best_move(pos);
        //fprintf(stderr, "Disposing MCTS tree\n");
        //fflush(stderr);
        root->dispose();
        //fprintf(stderr, "Done\n");
        //fflush(stderr);
		pos.make_move(move);
		//endgame_value = pos.dead_endgame_value();
        //if (endgame_value != NO_RESULT)
        //    fprintf(stderr, "Solved endgame: %d\n", endgame_value);
        move.to_str(move_str);
		send_move(move_str);
		get_move(move_str);
        pos.make_move(Move(move_str, pos.turn_));
	}
	return 0;
}




// EventHorizon0_0_0: plays move which maximises expected pure Monto-Carlo value
// EventHorizon1_0_0: untuned UCT implementation