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
    heap_ = std::vector<AMAFRecord*>();
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


void AMAFTable::play(U32 cell_index, U32 stone_index) {
    check_initialised();
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        records = cell_init(cell_index);
    AMAFRecord& record = records[stone_index];
#ifdef DEBUG_
    assert(! record.played_);
#endif
    if (record.n_playouts_ != 0)
        pop_heap();
    record.played_ = true;
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
        if (! record.played_)
            heap_insert(&record);
    } else {
        record.n_wins_ += win;
        record.n_playouts_++;
#ifdef DEBUG_
        assert(record.n_playouts_ < 60000);
#endif
        if (! record.played_)
            heap_update(&record, win);
    }
}


void AMAFTable::get_best(U32& cell_index, U32& stone_index) {
    if (heap_.empty()) {
        cell_index = NO_CELL;
        return;
    }
    AMAFRecord * record = heap_[0];
    cell_index = record->cell_index_;
    stone_index = record->stone_index_;
    //record->tried_ = true;
    //pop_heap();
}


bool AMAFTable::is_played(U32 cell_index, U32 stone_index) {
    if (amaf_.size() == 0)
        return false;
    AMAFRecord * records = amaf_[cell_index];
    if (records == NULL)
        return false; 
    return records[stone_index].played_;
}


void AMAFTable::heap_insert(AMAFRecord * record) {
    U32 index = heap_.size();
    heap_.push_back(NULL);
    while (index) {
        U32 parent = (index - 1) >> 1;
        AMAFRecord * rep = heap_[parent];
        if (*record < *rep)
            break;
        rep->index_ = index;
        heap_[index] = rep;
        index = parent;
    }
    record->index_ = index;
    heap_[index] = record;
}


void AMAFTable::pop_heap() {
#ifdef DEBUG_
    assert(heap_.size() > 0);
#endif
    //printf("In pop_heap\n");
    //fflush(stdout);
    AMAFRecord * record = heap_.back();
    heap_.pop_back();
    U32 size = heap_.size();
    U32 first_leaf = (size >> 1);
    U32 index = 0;
    while (index < first_leaf) {
        U32 left = (index << 1) + 1;
        U32 right = left + 1;
        AMAFRecord * child;
        if (right < size && *heap_[left] < *heap_[right]) {
            child = heap_[right];
            if (*record < *child) {
                child->index_ = index;
                heap_[index] = child;
                index = right;
            } else {
                break;
            }
        } else {
            child = heap_[left];
            if (*record < *child) {
                child->index_ = index;
                heap_[index] = child;
                index = left;
            } else {
                break;
            }
        }
    }
    record->index_ = index;
    heap_[index] = record;
}


void AMAFTable::heap_update(AMAFRecord * record, bool win) {
    //printf("In heap_update\n");
    //fflush(stdout);
    U32 index = record->index_;
    if (win) {
        //printf("Win\n");
        //fflush(stdout);
        while (index) {
            U32 parent = (index - 1) >> 1;
            AMAFRecord * rep = heap_[parent];
            if (*record < *rep)
                break;
            rep->index_ = index;
            heap_[index] = rep;
            index = parent;
        }
        record->index_ = index;
        heap_[index] = record;
    } else {
        //printf("Loss\n");
        //fflush(stdout);
        U32 size = heap_.size();
        U32 first_leaf = (size >> 1);
        while (index < first_leaf) {
            U32 left = (index << 1) + 1;
            U32 right = left + 1;
            AMAFRecord * child;
            if (right < size && *heap_[left] < *heap_[right]) {
                child = heap_[right];
                if (*record < *child) {
                    child->index_ = index;
                    heap_[index] = child;
                    index = right;
                } else {
                    break;
                }
            } else {
                child = heap_[left];
                if (*record < *child) {
                    child->index_ = index;
                    heap_[index] = child;
                    index = left;
                } else {
                    break;
                }
            }
        }
        record->index_ = index;
        heap_[index] = record;
    }
}

