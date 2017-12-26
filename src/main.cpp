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
    
    Position pos(blocked, 0);
    char move_str[100];
    get_move(move_str);
	//Value endgame_value;
    if (strcmp(move_str, "Start"))
        pos.make_move(Move(move_str, pos.turn_));
	Value alpha = 0;
    while (true) {
        fprintf(stderr, "Time left: %.2f seconds\n", (float)time_left / 1e6);
		//endgame_value = pos.dead_endgame_value();
        //if (endgame_value != NO_RESULT)
        //    fprintf(stderr, "Solved endgame: %d\n", endgame_value);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES);
#endif
        //fprintf(stderr, "Creating MCTS tree\n");
        //fflush(stderr);
        pos.set_alpha(0);
        fprintf(stderr, "Control heuristic: %.3f\n", pos.get_control_heuristic());
        MCTSearch search(pos, alpha);
        //fprintf(stderr, "Performing MCTS search\n");
        //fflush(stderr);
        Move move = search.get_best_move(pos);
        alpha = search.current_alpha_;
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
// EventHorizon1_0_0: untuned UCT implementation, scores -8 vs 0_0_0 (100 games)
// EventHorizon1_0_1: 0.3s/move, scores -3.66 vs 0_0_0 (100 games)
// EventHorizon1_0_2: remove move ordering, scores -2.82 vs 0_0_0 (100 games)
// EventHorizon1_0_3: use bitmask move table, scores -2.55  vs 0_0_0 (100 games)
// EventHorizon1_1_0: fixed target of 0, wins a lot of games by 1 point, loses horribly when it loses, scores -8.86  vs 0_0_0 (100 games)
// EventHorizon1_1_2: moving target, scores -6.02  vs 0_0_0 (100 games)
// EventHorizon1_1_3: don't play senseless moves in isolated cells, -4.74  vs 0_0_0 (100 games)
// EventHorizon1_1_4: maximise expectation for first 5 moves, create default policy, -0.39  vs 0_0_0 (100 games)
// EventHorizon1_1_5: disallow senseless moves with doubles, +0.2  vs 0_0_0 (100 games)
