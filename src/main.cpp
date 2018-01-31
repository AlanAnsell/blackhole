#include "globals.h"
#include "Position.h"
#include "MCTS.h"


void analyse_pos(Position& pos) {
    char move_str[10];
    while (scanf("%s", move_str) == 1) {
        Move move = move_from_str(move_str);
        pos.make_move(move, false);
    }
    fprintf(stderr, "Analysing...\n");
    fflush(stderr);
    pos.print(stderr);
    MCTSearch search(0);
    search.analyse(pos);
}


int main(int argc, char** argv) {
    time_started = get_time();
    if (argc > 1) {
        int millis;
        if (sscanf(argv[1], "-t=%d", &millis) == 1) {
            time_limit = 1000LL * (long long)millis;
            time_left = (49 * time_limit) / 50;
        }
    }
	fprintf(stderr, "Before init\n");
    fflush(stderr);
    init();
	fprintf(stderr, "After init\n");
    fflush(stderr);
	init_free_list();
    init_amaf_free_list();
    	
    char cell_str[5];
    U32 blocked[N_BLOCKED_CELLS];
    for (U32 i = 0; i < N_BLOCKED_CELLS; i++) {
        assert(scanf("%s", cell_str));
        blocked[i] = cell_name_to_id(cell_str);
    }
    
	fprintf(stderr, "Making position\n");
    fflush(stderr);
    Position pos(blocked, 0);
	fprintf(stderr, "Made position\n");
    fflush(stderr);
    if (argc > 1 && ! strcmp(argv[1], "-a")) {
        fprintf(stderr, "About to analyse\n");
        fflush(stderr);
        analyse_pos(pos);
        return 0;
    }
    fprintf(stderr, "Didn't analyse\n");
    fflush(stderr);

    char move_str[100];
    get_move(move_str);
    if (strcmp(move_str, "Start"))
        pos.make_move(move_from_str(move_str), false);
	Value alpha = 0;
    while (pos.open_.len_ > 1) {
        fprintf(stderr, "Time left: %.2f seconds\n", (float)(time_left - (get_time() - time_started)) / 1e6);
#ifdef DEBUG_
        assert(n_free == N_MCT_NODES);
#endif
        MCTSearch search(alpha);
        Move move = search.get_best_move(pos);
        alpha = search.current_alpha_;
        fprintf(stderr, "Used %u nodes\n", N_MCT_NODES - n_free);
        fprintf(stderr, "Used %u AMAF records\n", amaf_pointer);
		pos.make_move(move, false);
        move_to_str(move, move_str);
		send_move(move_str);
		get_move(move_str);
        pos.make_move(move_from_str(move_str), false);
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
// EventHorizon3_0_3: make opening heuristic calculation more efficient and more optimistic
// EventHorizon4_0_0: implement RAVE, modify default policy, alter solver move ordering 
// EventHorizon4_0_1: modify memory allocation to be more time and space efficient 
// EventHorizon4_0_2: fix bugs related to accelerated win detection
// EventHorizon4_1_1: improve solver and (hopefully) fix some crashes. 
// EventHorizon4_1_7: fix bugs relating to accessing invalid power array entries and root selection 
// EventHorizon4_2_2: more intelligently detect unreasonable moves 
// EventHorizon4_2_3: apply more intelligent unreasonable move detection to solver 
// EventHorizon4_2_4: tidy up a bit and inline a few functions 


