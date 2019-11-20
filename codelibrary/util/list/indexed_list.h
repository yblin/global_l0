//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_LIST_INDEXED_LIST_H_
#define UTIL_LIST_INDEXED_LIST_H_

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "codelibrary/base/array.h"
#include "codelibrary/base/object_pool.h"

namespace cl {

/**
 * Indexed list handles nodes with pluggable properties.
 *
 * Indexed list has the following features:
 *   1. Allocate a node and generate a unique ID for this node in O(1) time.
 *   2. Deallocate a node in O(1) time.
 *   3. Access a node by its ID in O(1) time.
 *   4. Dynamic add or delete one or more properties in O(1) time.
 *   5. Access the property of a node in O(1) time.
 *
 * Indexed list is useful for the data structures that need to dynamic access
 * properties of each node.
 */
template <typename BaseNode>
class IndexedList {
    static_assert(std::is_class<BaseNode>::value,
                  "'BaseNode' should be a class.");

    struct BaseProperty {
        BaseProperty() = default;
        virtual ~BaseProperty() = default;
        virtual void Add() = 0;
        virtual void Reset(int) = 0;
        virtual void Resize(int) = 0;
    };

    template <typename T>
    class PropertyArray : public BaseProperty {
        friend class IndexedList;

    public:
        PropertyArray(int size, const T& initial_value)
            : BaseProperty(),
              data_(size, initial_value),
              initial_value_(initial_value) {}

        T& operator[] (int index) {
            assert(index < data_.size());

            return data_[index];
        }

        const T& operator[] (int index) const {
            assert(index < data_.size());

            return data_[index];
        }

        /**
         * Add a new element.
         */
        virtual void Add() {
            data_.push_back(initial_value_);
        }

        /**
         * Reset the index-th element.
         */
        virtual void Reset(int index) {
            data_[index] = initial_value_;
        }

        /**
         * Resize the property data.
         */
        virtual void Resize(int size) {
            data_.resize(size);
        }

    private:
        Array<T> data_;
        T initial_value_;
    };

public:
    /**
     * Node equipped with a unique ID.
     */
    class Node : public BaseNode {
        friend class IndexedList;
        friend class ObjectPool<Node>;

    public:
        Node() = default;

        int id() const { return id_; }

    private:
        int id_ = -1;
    };

    /**
     * Property for IndexedList.
     *
     * Sample usage:
     *   auto a = o.AddProperty("ta", 0);
     */
    template <typename T>
    class Property {
        friend class IndexedList;

        explicit Property(std::shared_ptr<PropertyArray<T> > property_array)
            : property_array_(property_array) {}

    public:
        Property() = default;

        T& operator[] (const Node* o) {
            return property_array_->operator [](o->id());
        }

        const T& operator[] (const Node* o) const {
            return property_array_->operator [](o->id());
        }

    private:
        std::shared_ptr<PropertyArray<T> > property_array_;
    };

    /**
     * Iterator for available nodes.
     */
    class Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = Node*;
        using difference_type   = int;
        using pointer           = const value_type*;
        using reference         = const value_type&;

        Iterator() = default;

        Iterator(const Array<Node*>* nodes, int index)
            : nodes_(nodes), index_(index) {}

        bool operator ==(const Iterator& rhs) const {
            return nodes_ == rhs.nodes_ && index_ == rhs.index_;
        }

        bool operator !=(const Iterator& rhs) const {
            return !(*this == rhs);
        }

        Node* operator*()  const { return (*nodes_)[index_]; }
        Node* operator->() const { return (*nodes_)[index_]; }

        Iterator& operator++() {
            ++index_;
            return *this;
        }

        Iterator& operator--() {
            --index_;
            return *this;
        }

    protected:
        const Array<Node*>* nodes_ = nullptr;
        int index_ = -1;
    };

    /**
     * Const iterator for available nodes.
     */
    class ConstIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = const Node*;
        using difference_type   = int;
        using pointer           = const value_type*;
        using reference         = const value_type&;

        ConstIterator() = default;

        ConstIterator(const Array<Node*>* nodes, int index)
            : nodes_(nodes), index_(index) {}

        bool operator ==(const ConstIterator& rhs) const {
            return nodes_ == rhs.nodes_ && index_ == rhs.index_;
        }

        bool operator !=(const ConstIterator& rhs) const {
            return !(*this == rhs);
        }

        const Node* operator*()  const { return (*nodes_)[index_]; }
        const Node* operator->() const { return (*nodes_)[index_]; }

        ConstIterator& operator++() {
            ++index_;
            return *this;
        }

        ConstIterator& operator--() {
            --index_;
            return *this;
        }

    protected:
        const Array<Node*>* nodes_ = nullptr;
        int index_ = -1;
    };

    IndexedList() = default;

    virtual ~IndexedList() = default;

    IndexedList(const IndexedList& list) {
        list.Clone(this);
    }

    IndexedList& operator =(const IndexedList& list) {
        list.Clone(this);
        return *this;
    }

    /**
     * Allocate a node.
     */
    Node* Allocate() {
        Node* o = pool_.Allocate();
        nodes_.push_back(o);
        if (++n_available_ > n_allocated_) {
            o->id_ = n_allocated_;
            id_to_index_.push_back(n_allocated_);
            ++n_allocated_;

            for (auto& i : property_map_) {
                i.second->Add();
            }
        }

        return o;
    }

    /**
     * Deallocate a node.
     */
    void Deallocate(Node* node) {
        assert(node);
        assert(node->id_ < n_allocated_);

        int pos = id_to_index_[node->id_];
        assert(pos < n_available_ && "Invalid object");

        Node* back = nodes_.back();
        nodes_[pos] = back;
        id_to_index_[back->id_] = pos;
        id_to_index_[node->id_] = --n_available_;
        nodes_.pop_back();

        pool_.Deallocate(node);

        for (auto& i : property_map_) {
            i.second->Reset(node->id_);
        }
    }

    /**
     * Clone indexed list.
     */
    void Clone(IndexedList* list) const {
        assert(list);

        if (this == list) return;

        list->clear();
        list->n_allocated_ = n_allocated_;
        list->n_available_ = n_available_;
        list->id_to_index_ = id_to_index_;
        list->nodes_.resize(n_available_);

        Array<Node*> tmp(n_allocated_);
        for (int i = 0; i < n_allocated_; ++i) {
            Node* o = list->pool_.Allocate();
            o->id_ = i;
            tmp[i] = o;
        }

        for (int i = 0; i < n_available_; ++i) {
            *tmp[nodes_[i]->id_] = *nodes_[i];
        }

        Array<bool> used(n_allocated_, false);
        for (int i = 0; i < n_allocated_; ++i) {
            if (id_to_index_[i] < n_available_) {
                used[i] = true;
                list->nodes_[id_to_index_[i]] = tmp[i];
            }
        }

        for (int i = 0; i < n_allocated_; ++i) {
            if (!used[i]) list->pool_.Deallocate(tmp[i]);
        }

        // We do not copy the properties, but just resize them.
        for (auto i : property_map_) {
            i.second->Resize(n_allocated_);
        }
    }

    /**
     * Add a new property.
     */
    template <typename T>
    Property<T> AddProperty(const std::string& name, const T& initial_value) {
        if (property_map_.find(name) != property_map_.end()) {
            assert(false && "The property is already exist.");
        }

        std::shared_ptr<PropertyArray<T> > property_array(
                new PropertyArray<T>(n_allocated_, initial_value));
        property_map_[name] = property_array;
        property_type_[name] = typeid(PropertyArray<T>).name();
        return Property<T>(property_array);
    }
    template <typename T>
    Property<T> AddProperty(const std::string& name) {
        return AddProperty<T>(name, T());
    }

    /**
     * Add a property for const reference of list.
     */
    template <typename T>
    Property<T> AddProperty(const T& initial_value = T()) const {
        std::shared_ptr<PropertyArray<T> > property_array(
                new PropertyArray<T>(n_allocated_, initial_value));
        return Property<T>(property_array);
    }

    /**
     * Get a property by name.
     */
    template <typename T>
    Property<T> GetProperty(const std::string& name) {
        auto iter = property_map_.find(name);
        if (iter == property_map_.end()) {
            assert(false && "The property is not exist.");
        }

        assert(std::string(typeid(PropertyArray<T>).name()) ==
               property_type_[name] && 
               "The value type of property is not matched");

        auto res = std::dynamic_pointer_cast<PropertyArray<T> >(iter->second);
        return  Property<T>(res);
    }

    /**
     * Remove a property by its name.
     */
    void EraseProperty(const std::string& name) {
        property_map_.erase(name);
        property_type_.erase(name);
    }

    /**
     * Return the total number of allocated nodes.
     */
    int n_allocated() const {
        return n_allocated_;
    }

    /**
     * Return the number of available nodes.
     */
    int n_available() const {
        return n_available_;
    }

    /**
     * Return available nodes.
     */
    const Array<Node*>& nodes() const {
        return nodes_;
    }

    Iterator begin()            { return Iterator(&nodes_, 0);      }
    ConstIterator begin() const { return ConstIterator(&nodes_, 0); }
    Iterator end()            { return Iterator(&nodes_, n_available_);      }
    ConstIterator end() const { return ConstIterator(&nodes_, n_available_); }

    /**
     * Return true if no object is allocated.
     */
    bool empty() const {
        return n_allocated_ == 0;
    }

    /**
     * Clear all nodes. Note that, we do not clear the properties.
     */
    void clear() {
        n_available_ = 0;
        n_allocated_ = 0;
        id_to_index_.clear();
        nodes_.clear();
        pool_.clear();
    }

    Node* operator[] (int id) {
        assert(id >= 0 && id < n_allocated_);
        assert(id_to_index_[id] < n_available_);

        return nodes_[id_to_index_[id]];
    }

    const Node* operator[] (int id) const {
        assert(id >= 0 && id < n_allocated_);
        assert(id_to_index_[id] < n_available_);

        return nodes_[id_to_index_[id]];
    }

    /**
     * Return true if the pointer is available.
     */
    bool IsAvailable(const Node* o) const {
        assert(o);

        return o->id_ >= 0 && o->id_ < n_allocated_ &&
               id_to_index_[o->id_] < n_available_;
    }

private:
    // The number of available nodes.
    int n_available_ = 0;

    // The number of allocated nodes.
    int n_allocated_ = 0;

    // Convert the unique ID of each object to the index of the data array.
    Array<int> id_to_index_;

    // Store the pointer to each object. It can be used to fast access the
    // available nodes.
    Array<Node*> nodes_;

    // ObjectPool to handle the memory.
    ObjectPool<Node> pool_;

    // Properties for each object.
    std::unordered_map<std::string, std::shared_ptr<BaseProperty> >
        property_map_;

    // The type of each property.
    std::unordered_map<std::string, std::string> property_type_;
};

} // namespace cl

#endif // UTIL_LIST_INDEXED_LIST_H_
