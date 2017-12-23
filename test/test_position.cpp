#include "../src/globals.h"
#include "../src/Position.h"

void test_expectation_calculation() {
    printf("test_expectation_calculation\n");
    std::vector<CellID> V;
    size_t i;
    for (i = 0; i < N_CELLS; i++)
        V.push_back(i);
    std::random_shuffle(V.begin(), V.end());
    CellID blocked[N_BLOCKED_CELLS];
    //printf("Blocked cells:");
    for (i = 0; i < N_BLOCKED_CELLS; i++) {
        //char blocked_str[5];
        blocked[i] = V[i];
        //cell_id_to_name(blocked[i], blocked_str);
        //printf(" %s", blocked_str);
    }
    //printf("\nCreating position\n");
    //fflush(stdout);
    Position position(blocked);
    //position.print();
    //printf("Created position\n");
    //fflush(stdout);
    for (i = 0; i < N_CELLS - 12; i++) {
        //printf("Getting random move\n");
        //fflush(stdout);
        Move move = position.get_random_move();
        char move_str[10];
        //move.to_str(move_str);
        //printf("Got %s\n", move_str);
        //fflush(stdout);
        position.make_move(move);
        //printf("Made move\n");
        //fflush(stdout);
        //position.print();
    }
    //position.print();
    Real calculated_expectation = position.calculate_expectation();
    Real searched_expectation = position.search_expectation();
    //position.print();
    printf("searched_expectation = %.5lf\n", searched_expectation);
    printf("calculated_expectation = %.5lf\n", calculated_expectation);
    assert(fabs(searched_expectation - calculated_expectation) < 1e-6);
}



int main() {
    init();
    test_expectation_calculation();
    return 0;
}
