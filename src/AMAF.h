#ifndef AMAF_H_
#define AMAF_H_

#include "globals.h"

#define UNTRIED 0
#define TRIED 1


struct AMAFRecord {
    U16 n_wins_;
    U16 n_playouts_;
    U8 cell_id_;
    U8 stone_index_;

    inline void init_generated() {
        n_playouts_ = 0;
    }
     
    void init(U32 cell_id, U32 stone_index, bool win) {
        cell_id_ = (U8)cell_id;
        stone_index_ = (U8)stone_index;
        n_wins_= (U16)win;
        n_playouts_ = 1;
    }
    
    inline bool operator < (const AMAFRecord& other) const {
        return (U64)n_wins_ * (U64)other.n_playouts_ < (U64)other.n_wins_ * (U64)n_playouts_;
    }

};

#define N_AMAF_RECORDS 5000000

extern AMAFRecord amaf_store[N_AMAF_RECORDS];
extern U32 amaf_pointer;

void init_amaf_free_list();

AMAFRecord * amaf_alloc(U32 n);

AMAFRecord * get_amaf_record();

class AMAFTable {
public:
    U32 n_stones_;
    std::vector<AMAFRecord*> amaf_;
    std::vector<AMAFRecord*> list_;
     
    void init(U32 n_stones);
   
    AMAFRecord * cell_init(U32 cell_id);
    
    inline void check_initialised() {
        if (amaf_.size() == 0)
            amaf_ = std::vector<AMAFRecord*>(N_CELLS, NULL);
    }
     
    U32 get_n_playouts(U32 cell_id, U32 stone_index);
    
    Real get_value(U32 cell_id, U32 stone_index);

    void update(U32 cell_id, U32 stone_index, bool win); 

    void get_best(U32& cell_id, U32& stone_index, const std::vector<U16>& tried) const;

    void print(FILE * f) const;
};

#endif
