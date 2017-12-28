#ifndef AMAF_H_
#define AMAF_H_

#include "globals.h"


struct AMAFRecord {
    size_t n_wins_;
    size_t n_playouts_;

    AMAFRecord(): n_wins_(0), n_playouts_(0) {}

};


class AMAFTable {
public:
    size_t n_stones_;
    std::vector<std::vector<AMAFRecord>> amaf_;

    void init(size_t n_open, size_t n_stones) {
        n_stones_ = n_stones;
        amaf_ = std::vector<std::vector<AMAFRecord>>(n_open, std::vector<AMAFRecord>());
    }

    size_t get_n_playouts(size_t cell_index, size_t stone_index) {
        if (amaf_[cell_index].empty())
            return 0;
        return amaf_[cell_index][stone_index].n_playouts_;
    }

    Real get_value(size_t cell_index, size_t stone_index) {
#ifdef DEBUG_
        assert(! amaf_[cell_index].empty());
#endif
        AMAFRecord& record = amaf_[cell_index][stone_index];
#ifdef DEBUG_
        assert(record.n_playouts_ != 0);
#endif
        return (Real)record.n_wins_ / (Real)record.n_playouts_;
    }

    void update(size_t cell_index, size_t stone_index, bool win) {
        std::vector<AMAFRecord>& records = amaf_[cell_index];
        if (records.empty())
            records = std::vector<AMAFRecord>(n_stones_, AMAFRecord());
        AMAFRecord& record = records[stone_index];
        record.n_wins_ += win;
        record.n_playouts_++;
    }

};

#endif
