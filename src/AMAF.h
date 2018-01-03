#ifndef AMAF_H_
#define AMAF_H_

#include "globals.h"


struct AMAFRecord {
    U32 n_wins_;
    U32 n_playouts_;

    AMAFRecord(): n_wins_(0), n_playouts_(0) {}

};


class AMAFTable {
public:
    U32 n_stones_;
    std::vector<std::vector<AMAFRecord>> amaf_;

    void init(U32 n_open, U32 n_stones) {
        n_stones_ = n_stones;
        amaf_ = std::vector<std::vector<AMAFRecord>>(n_open, std::vector<AMAFRecord>());
    }

    U32 get_n_playouts(U32 cell_index, U32 stone_index) {
        if (amaf_[cell_index].empty())
            return 0;
        return amaf_[cell_index][stone_index].n_playouts_;
    }

    Real get_value(U32 cell_index, U32 stone_index) {
#ifdef DEBUG_
        assert(! amaf_[cell_index].empty());
#endif
        AMAFRecord& record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
        assert(record.n_playouts_ != 0);
#endif
        return (Real)record.n_wins_ / (Real)record.n_playouts_;
    }

    void update(U32 cell_index, U32 stone_index, bool win) {
        std::vector<AMAFRecord>& records = amaf_[cell_index];
        if (records.empty())
            records = std::vector<AMAFRecord>(n_stones_, AMAFRecord());
        AMAFRecord& record = records[stone_index];
        record.n_wins_ += win;
        record.n_playouts_++;
    }

};

#endif
