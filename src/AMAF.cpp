#include "AMAF.h"

AMAFRecord amaf_store[N_AMAF_RECORDS];
AMAFRecord * amaf_free_list[N_AMAF_RECORDS];
U32 n_amaf_free;

//const AMAFRecord * UNREASONABLE - new AMAFRecord;

void init_amaf_free_list() {
    for (n_amaf_free = 0; n_amaf_free < N_AMAF_RECORDS; n_amaf_free++)
        amaf_free_list[n_amaf_free] = &amaf_store[n_amaf_free];
}


AMAFRecord * get_amaf_record() {
#ifdef DEBUG_
    assert(n_amaf_free > 0);
#endif
    n_amaf_free--;
    return amaf_free_list[n_amaf_free];
}


void put_amaf_record(AMAFRecord * record) {
    amaf_free_list[n_amaf_free++] = record;
}


void AMAFTable::init(U32 n_open, U32 n_stones) {
    n_stones_ = n_stones;
    amaf_ = std::vector<std::vector<AMAFRecord*>>(n_open, std::vector<AMAFRecord*>());
    heap_ = std::vector<AMAFRecord*>();
}


void AMAFTable::dispose() {
    for (U32 i = 0; i < amaf_.size(); i++) {
        std::vector<AMAFRecord*> records = amaf_[i];
        for (U32 j = 0; j < records.size(); j++)
            if (records[j] != NULL)
                put_amaf_record(records[j]);
    }
    amaf_.clear();
    heap_.clear();
}


U32 AMAFTable::get_n_playouts(U32 cell_index, U32 stone_index) {
    std::vector<AMAFRecord*>& records = amaf_[cell_index];
    if (records.empty())
        return 0;
    AMAFRecord * record = records[stone_index];
    if (record == NULL /*|| record == UNREASONABLE*/)
        return 0;
    return record->n_playouts_;
}

Real AMAFTable::get_value(U32 cell_index, U32 stone_index) {
#ifdef DEBUG_
    assert(! amaf_[cell_index].empty());
#endif
    AMAFRecord * record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
    assert(record != NULL /*&& record != UNREASONABLE*/ && record->n_playouts_ != 0);
#endif
    return (Real)record->n_wins_ / (Real)record->n_playouts_;
}

void AMAFTable::update(U32 cell_index, U32 stone_index, bool win) {
    //fprintf(stderr, "update called\n");
    //fflush(stderr);
    std::vector<AMAFRecord*>& records = amaf_[cell_index];
    if (records.empty())
        records = std::vector<AMAFRecord*>(n_stones_, NULL);
    AMAFRecord * record;
    if (records[stone_index] == NULL) {
        //printf("Getting new record\n");
        //fflush(stdout);
        record = get_amaf_record();
        record->init(cell_index, stone_index);
        record->n_playouts_ = 1;
        if (win)
            record->n_wins_ = 1;
        records[stone_index] = record;
        //printf("About to heap_insert\n");
        //fflush(stdout);
        heap_insert(record);
        //printf("Did heap_insert\n");
        //fflush(stdout);    
    } else {
        record = records[stone_index];
/*#ifdef DEBUG_
        assert(record != UNREASONABLE);
#endif*/
        record->n_wins_ += win;
        record->n_playouts_++;
        if (! record->tried_)
            heap_update(record, win);
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
    record->tried_ = true;
    pop_heap();
}


void AMAFTable::set_tried(U32 cell_index, U32 stone_index) {
    //fprintf(stderr, "set_tried called\n");
    //fflush(stderr);
    std::vector<AMAFRecord*>& records = amaf_[cell_index];
    if (records.empty())
        records = std::vector<AMAFRecord*>(n_stones_, NULL);
#ifdef DEBUG_
    assert(records[stone_index] == NULL);
#endif
    AMAFRecord * record = get_amaf_record();
    records[stone_index] = record;
    record->init(cell_index, stone_index);
    record->tried_ = true;
}


bool AMAFTable::is_tried(U32 cell_index, U32 stone_index) {
    std::vector<AMAFRecord*>& records = amaf_[cell_index];
    if (records.empty())
        return false; 
    AMAFRecord * record = records[stone_index];
    if (record == NULL)
        return false;
    return record->tried_;
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

