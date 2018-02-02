g++ -std=c++11 -D UNIFORM_TM -D AMAF_HORIZON=5 -D UCB_C=0.03 -D RAVE_C=1e-3 -D NO_FAST_ -D NO_DEBUG_ -Wall -O2 -o ../../bin/EventHorizon5_0_0 main.cpp Position.cpp globals.cpp MCTS.cpp AMAF.cpp
