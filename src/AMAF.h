#ifndef AMAF_H_
#define AMAF_H_

#include "globals.h"

#define UNTRIED 0
#define TRIED 1


struct AMAFRecord {
    U16 n_wins_;
    U16 n_playouts_;
    U16 index_;
    U8 cell_index_;
    U8 stone_index_;
    bool played_;
    
    inline void init_generated() {
        n_playouts_ = 0;
        played_ = false;
    }
     
    void init(U32 cell_index, U32 stone_index, bool win) {
        cell_index_ = (U8)cell_index;
        stone_index_ = (U8)stone_index;
        n_wins_= (U16)win;
        n_playouts_ = 1;
    }

    inline bool operator < (const AMAFRecord& other) const {
        return (U64)n_wins_ * (U64)other.n_playouts_ < (U64)other.n_wins_ * (U64)n_playouts_;
    }

};

#define N_AMAF_RECORDS 10000000

extern AMAFRecord amaf_store[N_AMAF_RECORDS];
extern U32 amaf_pointer;

void init_amaf_free_list();

AMAFRecord * amaf_alloc(U32 n);

AMAFRecord * get_amaf_record();

class AMAFTable {
public:
    U32 n_stones_;
    U32 n_open_;
    std::vector<AMAFRecord*> heap_;
    std::vector<AMAFRecord*> amaf_;
    
    void init(U32 n_open, U32 n_stones);
   
    AMAFRecord * cell_init(U32 cell_index);
    
    inline void check_initialised() {
        if (amaf_.size() == 0)
            amaf_ = std::vector<AMAFRecord*>(n_open_, NULL);
    }
     
    void play(U32 cell_index, U32 stone_index);
    
    U32 get_n_playouts(U32 cell_index, U32 stone_index);
    
    Real get_value(U32 cell_index, U32 stone_index);

    void update(U32 cell_index, U32 stone_index, bool win); 

    void get_best(U32& cell_index, U32& stone_index);

    bool is_played(U32 cell_index, U32 stone_index);

private:
    void heap_insert(AMAFRecord * record);

    void pop_heap();

    void heap_update(AMAFRecord * record, bool win);
};

#endif
