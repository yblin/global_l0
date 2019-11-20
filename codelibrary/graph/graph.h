//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GRAPH_GRAPH_H_
#define GRAPH_GRAPH_H_

#include <cassert>
#include <memory>
#include <string>

#include "codelibrary/base/array.h"
#include "codelibrary/util/list/indexed_list.h"

namespace cl {

/**
 * Graph, represented by adjacency list.
 */ 
class Graph {
public:
    class BaseEdge;
    using Edge = typename IndexedList<BaseEdge>::Node;

    class BaseEdge {
        friend class Graph;
        friend class EdgeList;

    public:
        int source()       const { return source_; }
        int target()       const { return target_; }
        Edge* twin()             { return twin_;   }
        const Edge* twin() const { return twin_;   }

    private:
        int source_; // Source vertex of this edge.
        int target_; // Target vertex of this edge.
        Edge* twin_; // Twin edge of this edge.
        Edge* prev_; // For EdgeList.
        Edge* next_; // For EdgeList.
    };

    /**
     * Adjacency list of an incident vertex.
     */ 
    class EdgeList {
        friend class Graph;

    public:
        /**
         * Iterator for EdgeList.
         */ 
        class Iterator {
            using iterator_category = std::forward_iterator_tag;
            using value_type        = Edge*;
            using difference_type   = int;
            using pointer           = const value_type*;
            using reference         = const value_type&;

        public:
            Iterator() = default;

            explicit Iterator(Edge* edge) : edge_(edge) {}

            bool operator == (const Iterator& rhs) const {
                return edge_ == rhs.edge_;
            }

            bool operator != (const Iterator& rhs) const {
                return !(*this == rhs);
            }

            Edge* operator*()  const { return edge_; }
            Edge* operator->() const { return edge_; }

            Iterator& operator++() {
                edge_ = edge_->next_;
                return *this;
            }

            Iterator& operator--() {
                edge_ = edge_->prev_;
                return *this;
            }

        protected:
            Edge* edge_ = nullptr;
        };

        /**
         * Const iterator for EdgeList.
         */ 
        class ConstIterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type        = const Edge*;
            using difference_type   = int;
            using pointer           = const value_type*;
            using reference         = const value_type&;

            ConstIterator() = default;
            explicit ConstIterator(Iterator i) : edge_(*i) {}
            explicit ConstIterator(const Edge* e) : edge_(e) {}

            bool operator == (const ConstIterator& rhs) const {
                return edge_ == rhs.edge_;
            }

            bool operator != (const ConstIterator& rhs) const {
                return !(*this == rhs);
            }

            const Edge* operator * () const { return edge_; }
            const Edge* operator ->() const { return edge_; }

            ConstIterator& operator++() {
                edge_ = edge_->next_;
                return *this;
            }

            ConstIterator& operator--() {
                edge_ = edge_->prev_;
                return *this;
            }

        protected:
            const Edge* edge_ = nullptr;
        };

        EdgeList() {
            head_.prev_ = &head_;
            head_.next_ = &head_;
        }

        Iterator begin() {
            return Iterator(head_.next_);
        }

        Iterator end() {
            return Iterator(&head_);
        }

        ConstIterator begin() const {
            return ConstIterator(head_.next_);
        }

        ConstIterator end() const {
            return ConstIterator(&head_);
        }

        /**
         * Return the number of edges in the list.
         */
        int size() const {
            return size_;
        }

        /**
         * Return if the list is empty.
         */
        bool empty() const {
            return size_ == 0;
        }

    private:
        /**
         * Push back a node into list.
         */
        void PushBack(Edge* e) {
            e->prev_ = head_.prev_;
            e->next_ = &head_;
            head_.prev_->next_ = e;
            head_.prev_ = e;
            ++size_;
        }

        /**
         * Erase the node from list.
         */
        void Erase(Edge* e) {
            e->next_->prev_ = e->prev_;
            e->prev_->next_ = e->next_;
            --size_;
        }

        Edge head_;    // The head of node list.
        int size_ = 0; // Number of edges in the list.
    };

    template <typename T>
    using EdgeProperty = typename IndexedList<BaseEdge>::template Property<T>;

    explicit Graph(int n_vertices = 0)
        : n_vertices_(n_vertices),
          edges_from_(n_vertices) {
        assert(n_vertices_ >= 0);
    }

    Graph(const Graph& graph) {
        graph.Clone(this);
    }

    Graph& operator=(const Graph& graph) {
        graph.Clone(this);
        return *this;
    }

    /**
     * Return the number of vertices.
     */
    int n_vertices() const {
        return n_vertices_;
    }

    /**
     * Return the number of edges.
     */
    int n_edges() const {
        return edges_.n_available();
    }

    /**
     * Return the number of allocated edges.
     */
    int n_allocated_edges() const {
        return edges_.n_allocated();
    }

    /**
     * Return true if the graph is empty.
     */
    bool empty() const {
        return n_vertices_ == 0;
    }

    /**
     * Clear the data of graph.
     */
    void clear() {
        n_vertices_ = 0;
        n_two_way_edges_ = 0;
        edges_from_.clear();
        edges_.clear();
    }

    /**
     * Resize the vertices of graph and clear the edges.
     */
    void Resize(int n_vertices) {
        assert(n_vertices >= 0);

        clear();
        n_vertices_ = n_vertices;
        Array<EdgeList> t(n_vertices);
        edges_from_.swap(t);
    }

    /**
     * Clone this graph.
     *
     * Note that, we do NOT clone the properties.
     */
    void Clone(Graph* graph) const {
        assert(graph);

        if (this == graph) return;

        graph->Resize(n_vertices_);
        graph->edges_ = edges_;
        graph->n_two_way_edges_ = n_two_way_edges_;

        using Iterator = typename IndexedList<BaseEdge>::Iterator;
        using ConstIterator = typename IndexedList<BaseEdge>::ConstIterator;

        ConstIterator e1 = edges_.begin();
        Iterator e2 = graph->edges_.begin();
        for (; e1 != edges_.end(); ++e1, ++e2) {
            if (e1->twin_)
                e2->twin_ = graph->edges_[e1->twin()->id()];
        }

        for (int i = 0; i < n_vertices_; ++i) {
            for (const Edge* e : edges_from_[i]) {
                graph->edges_from_.at(i).PushBack(graph->edges_[e->id()]);
            }
        }
    }

    /**
     * Insert a one-way edge, return the inserted edge.
     */
    Edge* InsertOneWayEdge(int source, int target) {
        assert(source >= 0 && source < n_vertices_);
        assert(target >= 0 && target < n_vertices_);

        Edge* e = edges_.Allocate();
        e->source_ = source;
        e->target_ = target;
        e->twin_   = nullptr;

        edges_from_[source].PushBack(e);
        return e;
    }

    /**
     * Insert a one-way edge, return the inserted edge.
     */
    Edge* InsertTwoWayEdge(int source, int target) {
        Edge* e1 = InsertOneWayEdge(source, target);
        Edge* e2 = InsertOneWayEdge(target, source);
        e1->twin_ = e2;
        e2->twin_ = e1;
        ++n_two_way_edges_;
        return e1;
    }

    /**
     * Find the first edge equal to (source, target).
     */
    Edge* FindEdge(int source, int target) {
        assert(source >= 0 && source < n_vertices_);
        assert(target >= 0 && target < n_vertices_);

        for (Edge* e : edges_from_[source]) {
            if (e->target_ == target) return e;
        }
        return nullptr;
    }

    /**
     * Find the first edge equal to (source, target).
     */
    const Edge* FindEdge(int source, int target) const {
        assert(source >= 0 && source < n_vertices_);
        assert(target >= 0 && target < n_vertices_);

        for (const Edge* e : edges_from_[source]) {
            if (e->target_ == target) return e;
        }
        return nullptr;
    }

    /**
     * Erase a one-way edge.
     */
    void EraseOneWayEdge(Edge* e) {
        if (e->twin_ != nullptr) {
            e->twin_->twin_ = nullptr;
        }
        edges_from_[e->source_].Erase(e);
        edges_.Deallocate(e);
    }

    /**
     * Erase a two-way edge.
     */
    void EraseTwoWayEdge(Edge* e) {
        EraseOneWayEdge(e);
        if (e->twin_ != nullptr) {
            EraseOneWayEdge(e->twin_);
            --n_two_way_edges_;
        }
    }

    /**
     * Erase all of edges associated with vertex i.
     */
    void EraseEdgesOfVertex(int i) {
        Array<Edge*> edges;
        for (Edge* e : edges_from(i)) {
            edges.push_back(e);
        }
        for (Edge* e : edges) {
            EraseTwoWayEdge(e);
        }
    }

    /**
     * Return the edges whose source vertex equal to v.
     */
    EdgeList& edges_from(int v) {
        assert(v >= 0 && v < n_vertices_);

        return edges_from_[v];
    }

    /**
     * Return the edges whose source vertex equal to v.
     */
    const EdgeList& edges_from(int v) const {
        assert(v >= 0 && v < n_vertices_);

        return edges_from_[v];
    }

    /**
     * Add a new edge property with a given name.
     */
    template <typename T>
    EdgeProperty<T> AddEdgeProperty(const std::string& name,
                                    const T& initial_value) {
        return edges_.AddProperty(name, initial_value);
    }
    template <typename T>
    EdgeProperty<T> AddEdgeProperty(const std::string& name) {
        return edges_.AddProperty(name, T());
    }

    /**
     * Add a const edge property.
     */
    template <typename T>
    EdgeProperty<T> AddEdgeProperty(const T& initial_value = T()) const {
        return edges_.AddProperty(initial_value);
    }

    /**
     * Erase a edge property.
     */
    template <typename T>
    void EraseEdgeProperty(const std::string& name) {
        edges_.EraseProperty(name);
    }

    /**
     * Check if the graph is bidirectional.
     */
    bool IsBidirectional() const {
        return n_two_way_edges_ * 2 == n_edges();
    }

private:
    // Number of vertices.
    int n_vertices_ = 0;

    // Number of two way edges.
    int n_two_way_edges_ = 0;

    // edge_lists_[v] stores the edges with the source vertex equal to v.
    Array<EdgeList> edges_from_;

    // Edge pool.
    IndexedList<BaseEdge> edges_;
};

} // namespace cl

#endif // GRAPH_GRAPH_H_
