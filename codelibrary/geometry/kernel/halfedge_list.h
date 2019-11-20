//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_HALFEDGE_LIST_H_
#define GEOMETRY_KERNEL_HALFEDGE_LIST_H_

#include <cassert>
#include <string>

#include "codelibrary/base/array.h"
#include "codelibrary/util/list/indexed_list.h"

namespace cl {
namespace geometry {

/**
 * Halfedge list data structure.
 *
 * A halfedge list data structure also known as the doubly connected edge list
 * (DCEL), it is an edge-centered data structure capable of maintaining
 * incidence information of vertices, edges and faces, for example for planar
 * maps, polyhedral, or other orientable, two-dimensional surfaces embedded in
 * arbitrary dimension.
 *
 * Halfedge contains a record for halfedges, vertices, and faces of the
 * subdivision.
 *
 * Each halfedge bounds a single face and stores following informations.
 * 1. A pointer to the source vertex.
 * 2. A pointer to its twin edge.
 * 3. A pointer to the next edge on the boundary of the incident face.
 * 4. A pointer to the previous edge on the boundary of the incident face.
 * 5. A pointer to the incident face (Note that, the halfedges along the
 *    boundary of a hole are called border halfedges and have no incident face).
 */
template <typename Point>
class HalfedgeList {
public:
    class BaseVertex;
    class BaseHalfedge;
    class BaseFace;

    using VertexList = IndexedList<BaseVertex>;
    using EdgeList   = IndexedList<BaseHalfedge>;
    using FaceList   = IndexedList<BaseFace>;

    using Vertex   = typename VertexList::Node;
    using Halfedge = typename EdgeList::Node;
    using Face     = typename FaceList::Node;

    // Base Vertex of HalfedgeList.
    class BaseVertex {
        friend class HalfedgeList;

    public:
        BaseVertex() = default;

        bool is_isolated()         const { return halfedge_ == nullptr; }
        const Point& point()       const { return point_;               }
        Halfedge* halfedge()             { return halfedge_;            }
        const Halfedge* halfedge() const { return halfedge_;            }

    private:
        Point point_;                  // The position of this vertex.
        Halfedge* halfedge_ = nullptr; // The pointer to incident halfedge.
    };

    // Base Halfedge of HalfedgeList.
    class BaseHalfedge {
        friend class HalfedgeList;

    public:
        BaseHalfedge() = default;

        Vertex* source()                  { return vertex_;                }
        Vertex* target()                  { return twin_->vertex_;         }
        Halfedge* twin()                  { return twin_;                  }
        const Halfedge* twin()      const { return twin_;                  }
        Halfedge* next()                  { return next_;                  }
        const Halfedge* next()      const { return next_;                  }
        Halfedge* prev()                  { return prev_;                  }
        const Halfedge* prev()      const { return prev_;                  }
        Face* face()                      { return face_;                  }
        const Face* face()          const { return face_;                  }
        const Vertex* source()      const { return vertex_;                }
        const Vertex* target()      const { return twin_->vertex_;         }
        const Point& source_point() const { return vertex_->point_;        }
        const Point& target_point() const { return twin_->vertex_->point_; }
        bool is_boundary()          const { return face_ == nullptr;       }

    private:
        Vertex* vertex_ = nullptr; // The pointer to source vertex.
        Halfedge* twin_ = nullptr; // The pointer to twin halfedge.
        Halfedge* next_ = nullptr; // The pointer to next halfedge.
        Halfedge* prev_ = nullptr; // The pointer to previous halfedge.
        Face* face_     = nullptr; // The pointer to the incident face.
    };

    // Base Face of HalfedgeList.
    class BaseFace {
        friend class HalfedgeList;

    public:
        BaseFace() = default;

        Halfedge* halfedge()             { return halfedge_; }
        const Halfedge* halfedge() const { return halfedge_; }

    private:
        Halfedge* halfedge_ = nullptr; // The pointer to incident halfedge.
    };

    // Iterators.
    using VertexIterator        = typename VertexList::Iterator;
    using VertexConstIterator   = typename VertexList::ConstIterator;
    using HalfedgeIterator      = typename EdgeList::Iterator;
    using HalfedgeConstIterator = typename EdgeList::ConstIterator;
    using FaceIterator          = typename FaceList::Iterator;
    using FaceConstIterator     = typename FaceList::ConstIterator;

    // Properties.
    template<class T>
    using VertexProperty   = typename VertexList::template Property<T>;
    template<class T>
    using HalfedgeProperty = typename EdgeList::template Property<T>;
    template<class T>
    using FaceProperty     = typename FaceList::template Property<T>;

public:
    HalfedgeList() = default;

    HalfedgeList(const HalfedgeList& list) {
        list.Clone(this);
    }

    HalfedgeList& operator =(const HalfedgeList& list) {
        list.Clone(this);
        return *this;
    }

    /**
     * Clear the data.
     */
    void clear() {
        vertices_.clear();
        halfedges_.clear();
        faces_.clear();
    }

    /**
     * Check if the halfedge list is empty.
     */
    bool empty() const {
        return n_vertices() == 0;
    }

    /**
     * Return the number of vertices.
     */
    int n_vertices() const {
        return vertices_.n_available();
    }

    /**
     * Return the number of halfedges.
     */
    int n_halfedges() const {
        return halfedges_.n_available();
    }

    /**
     * Return the number of faces.
     */
    int n_faces() const {
        return faces_.n_available();
    }

    /**
     * Return the number of allocated vertices.
     */
    int n_allocated_vertices() const {
        return vertices_.n_allocated();
    }

    /**
     * Return the number of allocated halfedges.
     */
    int n_allocated_halfedges() const {
        return halfedges_.n_allocated();
    }

    /**
     * Return the number of allocated faces.
     */
    int n_allocated_faces() const {
        return faces_.n_allocated();
    }

    /**
     * Add a new vertex with position p.
     *
     * Return the pointer to the new vertex.
     */
    Vertex* AddVertex(const Point& p) {
        Vertex* v = CreateVertex();
        v->point_ = p;
        return v;
    }

    /**
     * Add a new pair of unattached halfedges.
     * This function only insert a pair of halfedges into HalfedgeList, but
     * don't join these halfedges to the other existed halfedges.
     *
     * Parameters:
     *  source - the pointer to the source vertex.
     *  target - the pointer to the target vertex.
     *
     * Return the halfedge from source to target.
     */
    Halfedge* AddEdge(Vertex* source, Vertex* target) {
        Halfedge* e1 = CreateEdge();
        Halfedge* e2 = e1->twin_;

        e1->vertex_ = source;
        e2->vertex_ = target;

        if (source->is_isolated()) {
            source->halfedge_ = e1;
        }
        if (target->is_isolated()) {
            target->halfedge_ = e2;
        }

        return e1;
    }

    /**
     * Add a new face by its incident halfedge 'e'.
     */
    Face* AddFace(Halfedge* e) {
        assert(e->is_boundary());

        Face* f = CreateFace();
        f->halfedge_ = e;

        Halfedge* e1 = e;
        do {
            e1->face_ = f;
            e1 = e1->next_;
        } while (e1 != e);

        return f;
    }

    /**
     * Add a triangle by three halfedges (e1->e2->e3).
     */
    void AddTriangle(Halfedge* e1, Halfedge* e2, Halfedge* e3) {
        set_next(e1, e2);
        set_next(e2, e3);
        set_next(e3, e1);
    }

    /**
     * Find the halfedge from v1 to v2.
     */
    Halfedge* FindHalfedge(Vertex* v1, Vertex* v2) {
        assert(v1 != v2);

        Halfedge* e = v1->halfedge_;
        if (e == nullptr) return nullptr;

        Halfedge* tmp = e;
        do {
            if (tmp->target() == v2) {
                return tmp;
            }
            tmp = tmp->twin_->next_;
        } while (tmp != e);

        return nullptr;
    }

    /**
     * Find the halfedge with two given vertices v1 and v2.
     *
     * Return the const pointer to halfedge from v1 to v2.
     */
    const Halfedge* FindHalfedge(const Vertex* v1, const Vertex* v2) const {
        assert(v1 != v2);

        const Halfedge* e = v1->halfedge_;
        if (e == nullptr) return nullptr;

        const Halfedge* tmp = e;
        do {
            if (tmp->target() == v2) {
                return tmp;
            }
            tmp = tmp->twin_->next_;
        } while (tmp != e);

        return nullptr;
    }

    /**
     * Create a new pair of halfedges to join the a's target to b's source.
     *
     * Note that, the halfedges a and b have to belong to the same face.
     *
     * Parameters:
     *  a is the halfedge which target point will be joined.
     *  b is the halfedge which source point will be joined.
     *
     * Return the halfedge that from a's target to b's source.
     */
    Halfedge* JoinEdge(Halfedge* a, Halfedge* b) {
        assert(a->target() != b->source());
        assert(a->face() == b->face());

        Halfedge *e1, *e2;
        e1 = AddEdge(a->target(), b->source());
        e2 = e1->twin_;

        set_next(e2, a->next_);
        set_next(a, e1);
        set_next(b->prev_, e2);
        set_next(e1, b);

        if (a->face_) {
            e1->face_ = a->face_;
            a->face_->halfedge_ = a;

            AddFace(e2);
        }

        return e1;
    }

    /**
     * Erase an existing pair of halfedge.
     */
    void EraseEdge(Halfedge* e) {
        bool need_add_face = false;

        Halfedge* e1 = e->twin_;
        if (e1->face_ && !e->face_) {
            EraseFace(e1->face_);
        } else if (e->face_ && !e1->face_) {
            EraseFace(e->face_);
        } else if (e->face_ && e1->face_) {
            Face* f = e->face_;
            Face* f1 = e1->face_;
            if (f != f1) {
                // Erase two faces.
                EraseFace(f);
                EraseFace(f1);
                need_add_face = true;
            } else {
                if (f->halfedge_ == e)
                    f->halfedge_ = e->next_ == e1 ? e1->next_ : e->next_;
                if (f->halfedge_ == e1)
                    f->halfedge_ = e1->next_ == e ? e->next_ : e1->next_;
            }
        }

        Halfedge* e_prev = e->prev_;
        Halfedge* e1_next = e1->next_;
        set_next(e_prev, e1_next);

        Halfedge* e1_prev = e1->prev_;
        Halfedge* e_next = e->next_;
        set_next(e1_prev, e_next);

        if (e->source()->halfedge_ == e) {
            e->source()->halfedge_ = (e1_next == e ? nullptr : e1_next);
        }
        if (e1->source()->halfedge_ == e1) {
            e1->source()->halfedge_ = (e_next == e1 ? nullptr : e_next);
        }

        if (need_add_face) AddFace(e_next);

        halfedges_.Deallocate(e);
        halfedges_.Deallocate(e1);
    }

    /**
     * Erase face.
     */
    void EraseFace(Face* f) {
        // Reset the incident halfedges of face.
        Halfedge* e = f->halfedge_;
        assert(e);

        Halfedge* e1 = e;
        do {
            e1->face_ = nullptr;
            e1 = e1->next_;
        } while (e1 != e);

        faces_.Deallocate(f);
    }

    /**
     * Clear all faces.
     */
    void ClearFaces() {
        for (Face* f : faces_) {
            // Reset the incident halfedges of face.
            Halfedge* e = f->halfedge_;
            assert(e);

            Halfedge* e1 = e;
            do {
                e1->face_ = nullptr;
                e1 = e1->next_;
            } while (e1 != e);
        }

        faces_.clear();
    }

    /**
     * Split the edge 'e' by vertex 'v', required 'v' is isolated.
     *
     * Return one of the split halfedge (from e->source to v).
     */
    Halfedge* SplitEdge(Halfedge* e, Vertex* v) {
        assert(v->is_isolated());

        Halfedge* e_next = e->next_;
        Halfedge* e_prev = e->prev_;
        Halfedge* e_twin = e->twin_;
        Halfedge* e_twin_prev = e_twin->prev_;
        Halfedge* e_twin_next = e_twin->next_;

        Halfedge *e1, *e2, *e3, *e4;
        e1 = AddEdge(e->source(), v);
        e2 = e1->twin_;
        e3 = AddEdge(v, e->target());
        e4 = e3->twin_;

        set_next(e_prev, e1);
        set_next(e1, e3);
        set_next(e3, e_next);

        set_next(e_twin_prev, e4);
        set_next(e4, e2);
        set_next(e2, e_twin_next);

        e1->face_ = e->face_;
        e2->face_ = e_twin->face_;
        e3->face_ = e->face_;
        e4->face_ = e_twin->face_;
        if (e->face_ && e->face_->halfedge_ == e)
            e->face_->halfedge_ = e1;
        if (e_twin->face_ && e_twin->face_->halfedge_ == e_twin)
            e_twin->face_->halfedge_ = e2;

        if (e->source()->halfedge_ == e)
            e->source()->halfedge_ = e1;
        if (e_twin->source()->halfedge_ == e_twin) {
            e_twin->source()->halfedge_ = e4;
        }

        halfedges_.Deallocate(e);
        halfedges_.Deallocate(e_twin);

        return e1;
    }

    /**
     * Clone this HalfedgeList.
     *
     * Note that we do NOT clone the properties.
     */
    void Clone(HalfedgeList* list) const {
        assert(list);

        if (this == list) return;

        list->clear();
        list->vertices_  = vertices_;
        list->halfedges_ = halfedges_;
        list->faces_     = faces_;

        VertexConstIterator v1 = vertices_.begin();
        VertexIterator v2 = list->vertices_.begin();
        for (; v1 != vertices_.end(); ++v1, ++v2) {
            if (v1->halfedge_) {
                v2->halfedge_ = list->halfedges_[v1->halfedge_->id()];
            }
        }

        HalfedgeConstIterator e1 = halfedges_.begin();
        HalfedgeIterator e2 = list->halfedges_.begin();
        for (; e1 != halfedges_.end(); ++e1, ++e2) {
            if (e1->face_) {
                e2->face_ = list->faces_[e1->face_->id()];
            }
            e2->vertex_ = list->vertices_[e1->vertex_->id()];
            e2->next_   = list->halfedges_[e1->next_->id()];
            e2->prev_   = list->halfedges_[e1->prev_->id()];
            e2->twin_   = list->halfedges_[e1->twin_->id()];
        }

        FaceConstIterator f1 = faces_.begin();
        FaceIterator f2 = list->faces_.begin();
        for (; f1 != faces_.end(); ++f1, ++f2) {
            f2->halfedge_ = list->halfedges_[f1->halfedge_->id()];
        }
    }

    // Access functions.
    const VertexList& vertices() const { return vertices_;  }
    const EdgeList& halfedges()  const { return halfedges_; }
    const FaceList& faces()      const { return faces_;     }

          Vertex* vertex(int id)           { return vertices_[id];  }
    const Vertex* vertex(int id)     const { return vertices_[id];  }
          Halfedge* halfedge(int id)       { return halfedges_[id]; }
    const Halfedge* halfedge(int id) const { return halfedges_[id]; }
          Face* face(int id)               { return faces_[id];     }
    const Face* face(int id)         const { return faces_[id];     }

    HalfedgeIterator      begin()       { return halfedges_.begin(); }
    HalfedgeConstIterator begin() const { return halfedges_.begin(); }
    HalfedgeIterator      end()         { return halfedges_.end();   }
    HalfedgeConstIterator end()   const { return halfedges_.end();   }

    VertexIterator      vertex_begin()       { return vertices_.begin(); }
    VertexConstIterator vertex_begin() const { return vertices_.begin(); }
    VertexIterator      vertex_end()         { return vertices_.end();   }
    VertexConstIterator vertex_end()   const { return vertices_.end();   }

    HalfedgeIterator      halfedge_begin()       { return halfedges_.begin(); }
    HalfedgeConstIterator halfedge_begin() const { return halfedges_.begin(); }
    HalfedgeIterator      halfedge_end()         { return halfedges_.end();   }
    HalfedgeConstIterator halfedge_end()   const { return halfedges_.end();   }

    FaceIterator      face_begin()       { return faces_.begin(); }
    FaceConstIterator face_begin() const { return faces_.begin(); }
    FaceIterator      face_end()         { return faces_.end();   }
    FaceConstIterator face_end()   const { return faces_.end();   }

    /**
     * Add a vertex property.
     */
    template <typename T>
    VertexProperty<T> AddVertexProperty(const std::string& name,
                                        const T& initial_value) {
        return vertices_.AddProperty(name, initial_value);
    }

    /**
     * Add a const vertex property.
     */
    template <typename T>
    VertexProperty<T> AddConstVertexProperty(const T& initial_value) const {
        return vertices_.AddProperty(initial_value);
    }

    /**
     * Add a halfedge property.
     */
    template <typename T>
    HalfedgeProperty<T> AddHalfedgeProperty(const std::string& name,
                                            const T& initial_value = T()) {
        return halfedges_.AddProperty(name, initial_value);
    }

    /**
     * Add a const halfedge property.
     */
    template <typename T>
    HalfedgeProperty<T>
    AddHalfedgeProperty(const T& initial_value = T()) const {
        return halfedges_.AddProperty(initial_value);
    }

    /**
     * Add a face property.
     */
    template <typename T>
    FaceProperty<T> AddFaceProperty(const std::string& name,
                                    const T& initial_value = T()) {
        return faces_.AddProperty(name, initial_value);
    }

    /**
     * Add a const face property.
     */
    template <typename T>
    FaceProperty<T> AddFaceProperty(const T& initial_value = T()) const {
        return faces_.AddProperty(initial_value);
    }

    /**
     * Erase a vertex property with the given name.
     */
    template <typename T1>
    void EraseVertexProperty(const std::string& name) {
        vertices_.EraseProperty(name);
    }

    /**
     * Erase a halfedge property with the given name.
     */
    void EraseHalfedgeProperty(const std::string& name) {
        halfedges_.EraseProperty(name);
    }

    /**
     * Erase a face property with the given name.
     */
    void EraseFaceProperty(const std::string& name) {
        faces_.EraseProperty(name);
    }

    /**
     * Return true if the given face is available.
     */
    bool IsAvailable(const Vertex* v) const {
        return vertices_.IsAvailable(v);
    }

    /**
     * Return true if the given face is available.
     */
    bool IsAvailable(const Halfedge* e) const {
        return halfedges_.IsAvailable(e);
    }

    /**
     * Return true if the given face is available.
     */
    bool IsAvailable(const Face* f) const {
        return faces_.IsAvailable(f);
    }

    /**
     * set e1->next = e2, and e2->prev = e1.
     */
    void set_next(Halfedge* e1, Halfedge* e2) {
        assert(e1->target() == e2->source());

        e1->next_ = e2;
        e2->prev_ = e1;
    }

    /**
     * set e1->next = e2, and e2->prev = e1.
     */
    void set_prev(Halfedge* e1, Halfedge* e2) {
        assert(e2->target() == e1->source());

        e1->prev_ = e2;
        e2->next_ = e1;
    }

protected:
    /**
     * Create a new vertex, but does not change the other information in the
     * HalfedgeList.
     */
    Vertex* CreateVertex() {
        Vertex* v = vertices_.Allocate();
        v->halfedge_ = nullptr;
        return v;
    }

    /**
     * Create a new pairs of halfedges, but does not change the other
     * information in the HalfedgeList.
     */
    Halfedge* CreateEdge() {
        Halfedge* e1 = halfedges_.Allocate();
        Halfedge* e2 = halfedges_.Allocate();

        e1->twin_ = e2;
        e2->twin_ = e1;
        set_next(e1, e2);
        set_prev(e1, e2);
        e1->face_ = nullptr;
        e2->face_ = nullptr;

        return e1;
    }

    /**
     * Create a new pairs of halfedges, but does not change the other
     * information in the HalfedgeList.
     */
    Face* CreateFace() {
        Face* f = faces_.Allocate();
        f->halfedge_ = nullptr;
        return f;
    }

    VertexList vertices_;
    EdgeList   halfedges_;
    FaceList   faces_;
};

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_HALFEDGE_LIST_H_
