g++ -std=c++11 -D UNIFORM_TM -D AMAF_HORIZON=3 -D UCB_C=0.0044 -D RAVE_C=2e-4 -D NO_FAST_ -D NO_DEBUG_ -Wall -O2 -o ../../bin/EventHorizon4_2_5 main.cpp Position.cpp globals.cpp MCTS.cpp AMAF.cpp
