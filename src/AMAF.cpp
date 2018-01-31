#include "AMAF.h"

AMAFRecord amaf_store[N_AMAF_RECORDS];
U32 amaf_pointer;


void init_amaf_free_list() {
    amaf_pointer = 0;
}


AMAFRecord * amaf_alloc(U32 n) {
#ifdef DEBUG_
    assert(amaf_pointer + n <= N_AMAF_RECORDS);
#endif
    amaf_pointer += n;
    return &amaf_store[amaf_pointer - n];
}


void AMAFTable::init(U32 n_open, U32 n_stones) {
    n_stones_ = n_stones;
    n_open_ = n_open;
    amaf_.clear(); 
    list_.clear();
}


AMAFRecord * AMAFTable::cell_init(U32 cell_index) {
#ifdef DEBUG_
    assert(amaf_.size() > 0);
    assert(amaf_[cell_index] == NULL);
#endif
    AMAFRecord * records = amaf_alloc(n_stones_);
    amaf_[cell_index] = records;
    for (U32 i = 0; i < n_stones_; i++)
        records[i].init_generated();
    return records;
}


U32 AMAFTable::get_n_playouts(U32 cell_index, U32 stone_index) {
    if (amaf_.empty())
        return 0;
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        return 0;    
    return records[stone_index].n_playouts_;
}

Real AMAFTable::get_value(U32 cell_index, U32 stone_index) {
#ifdef DEBUG_
    assert(amaf_.size() > 0 && amaf_[cell_index] != NULL);
#endif
    AMAFRecord& record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
    assert(record.n_playouts_ != 0);
#endif
    return (Real)record.n_wins_ / (Real)record.n_playouts_;
}

void AMAFTable::update(U32 cell_index, U32 stone_index, bool win) {
    check_initialised();
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        records = cell_init(cell_index);
    AMAFRecord& record = records[stone_index];
    if (record.n_playouts_ == 0) {
        record.init(cell_index, stone_index, win);
        list_.push_back(record);
    } else {
        record.n_wins_ += win;
        record.n_playouts_++;
#ifdef DEBUG_
        assert(record.n_playouts_ < 60000);
#endif
    }
}


/*void AMAFTable::get_best(U32& cell_index, U32& stone_index) {
    if (heap_.empty()) {
        cell_index = NO_CELL;
        return;
    }
    AMAFRecord * record = heap_[0];
    cell_index = record->cell_index_;
    stone_index = record->stone_index_;
    //record->tried_ = true;
    //pop_heap();
}*/


