#include "../src/globals.h"
#include "../src/AMAF.h"
#include <set>

typedef std::pair<U32, U32> ii;

int O = 1000;

int main() {
    init();
    init_amaf_free_list();
    srand(42);
    for (int t = 0; t < 100; t++) {
        printf("Running test case %d\n", t);
        fflush(stdout);
        int n_open = rand() % 30 + 2;
        int n_stones = n_open / 2;
        printf("%d open, %d stones\n", n_open, n_stones);
        std::map<ii, int> val;
        int i, j;
        for (i = 0; i < n_open; i++)
            for (j = 0; j < n_stones; j++)
                val[ii(i, j)] = (Real)(rand() % (O + 1));
        AMAFTable amaf;
        amaf.init(n_open, n_stones);
        std::map<ii, ii> suc;
        std::set<std::pair<Real, ii>> best;
        std::set<ii> played;
        for (i = 0; i < n_open * n_stones; i++) {
            ii move;
            //printf("Doing AMAF updates\n");
            //fflush(stdout);
            for (j = 0; j < 10; j++) {
                move.first = rand() % n_open;
                move.second = rand() % n_stones;
                bool win = (rand() % O) < val[move];
                //printf("Move: %d %d (%s)\n", move.first, move.second, win ? "win" : "loss");
                //fflush(stdout);
                amaf.update(move.first, move.second, win);
                //printf("Did an update\n");
                //fflush(stdout);
                if (! played.count(move)) {
                    ii new_suc;
                    if (suc.count(move)) {
                        ii old_suc = suc[move];
                        Real old_suc_val = (Real)old_suc.first / (Real)old_suc.second;
                        assert(best.count(std::make_pair(old_suc_val, move)));
                        best.erase(std::make_pair(old_suc_val, move));
                        new_suc = old_suc;
                        new_suc.first += win;
                        new_suc.second++;
                    } else {
                        new_suc = ii((U32)win, 1);
                    }
                    suc[move] = new_suc;
                    Real new_suc_val = (Real)new_suc.first / (Real)new_suc.second;
                    best.insert(std::make_pair(new_suc_val, move));
                }
            }
            //printf("Done AMAF updates\n");
            //fflush(stdout);
            amaf.get_best(move.first, move.second);
            //printf("Selected move %d %d\n", move.first, move.second);
            //fflush(stdout);
            assert(move.first != NO_CELL);
            assert(! played.count(move));
            played.insert(move);
            auto it = best.end();
            it--;
            assert(fabs(it->first - amaf.get_value(move.first, move.second)) < 1e-6);
            assert(suc.count(move));
            Real suc_val = (Real)suc[move].first / (Real)suc[move].second;
            //for (it = best.begin(); it != best.end(); ++it)
            //    printf("%d %d (%.4f)\n", it->second.first, it->second.second, it->first);
            assert(best.count(std::make_pair(suc_val, move)));
            best.erase(std::make_pair(suc_val, move));
        }
        printf("all good\n");
        fflush(stdout);
    }
}




