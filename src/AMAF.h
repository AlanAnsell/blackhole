#ifndef AMAF_H_
#define AMAF_H_

#include "globals.h"

struct AMAFRecord {
    U32 n_wins_;
    U32 n_playouts_;
    U32 index_;
    U32 cell_index_;
    U32 stone_index_;
    bool tried_;

    void init(U32 cell_index, U32 stone_index) {
        cell_index_ = cell_index;
        stone_index_ = stone_index;
        n_wins_= 0;
        n_playouts_ = 0;
        tried_ = false;
    }

    inline bool operator < (const AMAFRecord& other) const {
        return (U64)n_wins_ * (U64)other.n_playouts_ < (U64)other.n_wins_ * (U64)n_playouts_;
    }

};

#define N_AMAF_RECORDS 2000000

extern AMAFRecord amaf_store[N_AMAF_RECORDS];
extern AMAFRecord * amaf_free_list[N_AMAF_RECORDS];
extern U32 n_amaf_free;

void init_amaf_free_list();

AMAFRecord * get_amaf_record();

void put_amaf_record(AMAFRecord * record);


class AMAFTable {
public:
    U32 n_stones_;
    std::vector<AMAFRecord*> heap_;
    std::vector<std::vector<AMAFRecord*>> amaf_;
    
    void init(U32 n_open, U32 n_stones);
     
    void dispose();
    
    U32 get_n_playouts(U32 cell_index, U32 stone_index);
    
    Real get_value(U32 cell_index, U32 stone_index);

    void update(U32 cell_index, U32 stone_index, bool win); 

    void get_best(U32& cell_index, U32& stone_index);

    void set_tried(U32 cell_index, U32 stone_index);

    bool is_tried(U32 cell_index, U32 stone_index);

private:
    void heap_insert(AMAFRecord * record);

    void pop_heap();

    void heap_update(AMAFRecord * record, bool win);
};

#endif
