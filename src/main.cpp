#include "globals.h"
#include "Position.h"
#include "MCTS.h"


int main(int argc, char** argv) {
    time_started = get_time();
    if (argc > 1) {
        int millis;
        if (sscanf(argv[1], "-t=%d", &millis) == 1) {
            time_limit = 1000LL * (long long)millis;
            time_left = (49 * time_limit) / 50;
        }
    }
	init();
	init_free_list();
	
    char cell_str[5];
    U32 blocked[N_BLOCKED_CELLS];
    for (U32 i = 0; i < N_BLOCKED_CELLS; i++) {
        assert(scanf("%s", cell_str));
        blocked[i] = cell_name_to_id(cell_str);
    }
    
    Position pos(blocked, 0);
    char move_str[100];
    get_move(move_str);
    if (strcmp(move_str, "Start"))
        pos.make_move(Move(move_str, pos.turn_));
	Value alpha = 0;
    while (pos.open_.len_ > 1) {
        fprintf(stderr, "Time left: %.2f seconds\n", (float)time_left / 1e6);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES);
#endif
        MCTSearch search(pos, alpha);
        Move move = search.get_best_move(pos);
        alpha = search.current_alpha_;
        fprintf(stderr, "Used %u nodes\n", N_MCT_NODES - n_free);
		pos.make_move(move);
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
// EventHorizon1_1_6: default policy with lottery, fix propagation for fully explored positions, +0.86 vs 0_0_0, +0.46 vs player3
// EventHorizon1_1_7: add dead tiles to Position +0.62 vs player3
// EventHorizon1_2_0: remove Cell class, add stale cell tracking
// EventHorizon1_2_1: don't generate moves with stones which can be placed in stale cells, +0.95 vs player3
// EventHorizon2_0_0: use AMAF for move ordering, +0.37 vs 1_2_2 (400 games)
// EventHorizon2_1_0: use bit operations to find stale cells
// EventHorizon2_2_0: use win prob with random play at alpha level for first 3 moves, +0.55 vs 2_1_0 (400 gmes)
// EventHorizon2_3_0: implement endgame solver
// EventHorizon2_4_3: add hash table to endgame solver, adjust time management
// EventHorizon2_4_6: adjust UCB coefficient
// EventHorizon3_0_0: generate only selected moves during expansion phase. Be more optimistic in target selection. 
// EventHorizon3_0_1: time usage
// EventHorizon3_0_2: tune UCB_C

