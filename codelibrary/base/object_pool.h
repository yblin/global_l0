//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_OBJECT_POOL_H_
#define BASE_OBJECT_POOL_H_

#include <algorithm>
#include <cassert>
#include <climits>

namespace cl {

/**
 * A simple but efficient object memory pool.
 *
 * ObjectPool is used for dynamic management of objects of the same size, it
 * can allocate and deallocate memory for objects quickly.
 *
 * ObjectPool implements a lazy lifecycle strategy. In this strategy objects are
 * default-constructed the first time they are allocated and destroyed when the
 * pool itself is destroyed.
 *
 * The time comparison between ObjectPool and new/delete is given below:
 *                        ObjectPool       new/delete
 * 10,000,000 allocate     0.093(s)         1.609(s)
 * 10,000,000 deallocate   0.032(s)         1.625(s)
 */
template <typename T>
class ObjectPool {
    // The chunk for ObjectPool.
    struct Chunk {
        explicit Chunk(int n)
            : size(n) {
            data = new (std::nothrow) T[n];
            assert(data && "Memory is not enough");

            ptr_data = new (std::nothrow) T*[n];
            assert(ptr_data && "Memory is not enough");

            for (int i = 0; i < size; ++i) {
                ptr_data[i] = &data[i];
            }
        }

        ~Chunk() {
            delete[] ptr_data;
            delete[] data;
        }

        int size = 0;           // The size of chunk.
        Chunk* prev = nullptr;  // The previous chunk of chunk list.
        Chunk* next = nullptr;  // The next chunk of chunk list.
        int used_size = 0;      // The size of used memory.
        T* data = nullptr;      // The actual data pointer.
        T** ptr_data = nullptr; // The pointer to the actual data.
    };

public:
    /**
     * ObjectPool is perform as a free chunk list.
     * The default size of first chunk is 1024, if user do not give it.
     */
    explicit ObjectPool(int first_chunk_size = 1024)
        : first_chunk_size_(first_chunk_size) {
        assert(first_chunk_size_ > 0);
        first_chunk_ = cur_chunk_ = new Chunk(first_chunk_size_);
    }

    ~ObjectPool() {
        ClearChunks();
    }

    /**
     * Allocate an object from pool.
     */
    T* Allocate() {
        assert(n_available_ < INT_MAX);

        if (cur_chunk_->used_size == cur_chunk_->size) {
            if (cur_chunk_->next == nullptr) {
                int n = cur_chunk_->size + cur_chunk_->size;
                if (n < 0) n = INT_MAX;

                // If not enough memory, we create a new chunk, this chunk's
                // size is double of the previous one.
                auto t = new Chunk(n);
                cur_chunk_->next = t;
                t->prev = cur_chunk_;
            }
            cur_chunk_ = cur_chunk_->next;
        }
        ++n_available_;
        n_allocated_ = std::max(n_allocated_, n_available_);

        return cur_chunk_->ptr_data[cur_chunk_->used_size++];
    }

    /**
     * Recycle a object, by putting it back to the pool.
     * Will assert() if object is null.
     */
    void Deallocate(T* object) {
        assert(object);
        assert(n_available_ > 0);

        if (cur_chunk_->used_size == 0) {
            cur_chunk_ = cur_chunk_->prev;
        }
        --n_available_;
        cur_chunk_->ptr_data[--cur_chunk_->used_size] = object;
    }

    /**
     * Clear the object pool.
     */
    void clear() {
        ClearChunks();
        n_available_ = 0;
        n_allocated_ = 0;
        cur_chunk_ = first_chunk_ = new Chunk(first_chunk_size_);
    }

    /**
     * Current available objects.
     */
    int n_available() const {
        return n_available_;
    }

    /**
     * The total number of allocated objects.
     */
    int n_allocated() const {
        return n_allocated_;
    }

private:
    /**
     * Clear the chunks.
     */
    void ClearChunks() {
        while (first_chunk_) {
            Chunk* p = first_chunk_->next;
            delete first_chunk_;
            first_chunk_ = p;
        }
    }

    int n_available_ = 0;          // The number of available objects.
    int n_allocated_ = 0;          // The number of allocated objects.
    Chunk* cur_chunk_ = nullptr;   // The current chunk.
    Chunk* first_chunk_ = nullptr; // The first chunk.
    int first_chunk_size_ = 0;     // The size of first chunk.
};

} // namespace cl

#endif // BASE_OBJECT_POOL_H_
