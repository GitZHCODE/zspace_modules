// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TRIANGLE_TRIANGLE_ADJACENCY_H
#define IGL_TRIANGLE_TRIANGLE_ADJACENCY_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // Constructs the triangle-triangle adjacency matrix for a given
  // mesh (V,F).
  //
  // Inputs:
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  // Outputs:
  //   TT   #F by #3 adjacent matrix, the element i,j is the id of the triangle
  //        adjacent to the j edge of triangle i
  //   TTi  #F by #3 adjacent matrix, the element i,j is the id of edge of the
  //        triangle TT(i,j) that is adjacent with triangle i
  //
  // NOTE: the first edge of a triangle is [0,1] the second [1,2] and the third
  //       [2,3].  this convention is DIFFERENT from
  //       cotmatrix_entries.h/edge_lengths.h/etc. To fix this you could use:
  //           // Fix mis-match convention
  //           {
  //             Eigen::PermutationMatrix<3,3> perm(3);
  //             perm.indices() = Eigen::Vector3i(1,2,0);
  //             TT = (TT*perm).eval();
  //             TTi = (TTi*perm).eval();
  //             for(int i=0;i<TTi.rows();i++)
  //               for(int j=0;j<TTi.cols();j++)
  //                 TTi(i,j)=TTi(i,j)==-1?-1:(TTi(i,j)+3-1)%3;
  //           }
  template <typename DerivedF, typename DerivedTT, typename DerivedTTi>
  IGL_INLINE void triangle_triangle_adjacency(
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedTT>& TT,
    Eigen::PlainObjectBase<DerivedTTi>& TTi);
  template <typename DerivedF, typename DerivedTT>
  IGL_INLINE void triangle_triangle_adjacency(
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedTT>& TT);
  // Preprocessing
  template <typename DerivedF, typename TTT_type>
  IGL_INLINE void triangle_triangle_adjacency_preprocess(
    const Eigen::MatrixBase<DerivedF>& F,
    std::vector<std::vector<TTT_type> >& TTT);
  // Extract the face adjacencies
  template <typename DerivedF, typename TTT_type, typename DerivedTT>
  IGL_INLINE void triangle_triangle_adjacency_extractTT(
    const Eigen::MatrixBase<DerivedF>& F,
    std::vector<std::vector<TTT_type> >& TTT,
    Eigen::PlainObjectBase<DerivedTT>& TT);
  // Extract the face adjacencies indices (needed for fast traversal)
  template <typename DerivedF, typename TTT_type, typename DerivedTTi>
  IGL_INLINE void triangle_triangle_adjacency_extractTTi(
    const Eigen::MatrixBase<DerivedF>& F,
    std::vector<std::vector<TTT_type> >& TTT,
    Eigen::PlainObjectBase<DerivedTTi>& TTi);
  // Adjacency list version, which works with non-manifold meshes
  //
  // Inputs:
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   TT  #F by 3 list of lists so that TT[i][c] --> {j,k,...} means that
  //     faces j and k etc. are edge-neighbors of face i on face i's edge
  //     opposite corner c
  //   TTj  #F list of lists so that TTj[i][c] --> {j,k,...} means that face
  //     TT[i][c][0] is an edge-neighbor of face i incident on the edge of face
  //     TT[i][c][0] opposite corner j, and TT[i][c][1] " corner k, etc.
  template <
    typename DerivedF, 
    typename TTIndex, 
    typename TTiIndex>
    IGL_INLINE void triangle_triangle_adjacency(
      const Eigen::MatrixBase<DerivedF> & F,
      std::vector<std::vector<std::vector<TTIndex> > > & TT,
      std::vector<std::vector<std::vector<TTiIndex> > > & TTi);
  template < typename DerivedF, typename TTIndex>
    IGL_INLINE void triangle_triangle_adjacency(
      const Eigen::MatrixBase<DerivedF> & F,
      std::vector<std::vector<std::vector<TTIndex> > > & TT);
  // Wrapper with bool to choose whether to compute TTi (this prototype should
  // be "hidden").
  template <
    typename DerivedF, 
    typename TTIndex, 
    typename TTiIndex>
    IGL_INLINE void triangle_triangle_adjacency(
      const Eigen::MatrixBase<DerivedF> & F,
      const bool construct_TTi,
      std::vector<std::vector<std::vector<TTIndex> > > & TT,
      std::vector<std::vector<std::vector<TTiIndex> > > & TTi);
  // Inputs:
  //   E  #F*3 by 2 list of all of directed edges in order (see
  //     `oriented_facets`)
  //   EMAP #F*3 list of indices into uE, mapping each directed edge to unique
  //     undirected edge
  //   uE2E  #uE list of lists of indices into E of coexisting edges
  // See also: unique_edge_map, oriented_facets
  template <
    typename DerivedE, 
    typename DerivedEMAP,
    typename uE2EType,
    typename TTIndex, 
    typename TTiIndex>
    IGL_INLINE void triangle_triangle_adjacency(
      const Eigen::MatrixBase<DerivedE> & E,
      const Eigen::MatrixBase<DerivedEMAP> & EMAP,
      const std::vector<std::vector<uE2EType > > & uE2E,
      const bool construct_TTi,
      std::vector<std::vector<std::vector<TTIndex> > > & TT,
      std::vector<std::vector<std::vector<TTiIndex> > > & TTi);
  template <
    typename DerivedEMAP,
    typename DeriveduEC,
    typename DeriveduEE,
    typename TTIndex, 
    typename TTiIndex>
  IGL_INLINE void triangle_triangle_adjacency(
    const Eigen::MatrixBase<DerivedEMAP> & EMAP,
    const Eigen::MatrixBase<DeriveduEC> & uEC,
    const Eigen::MatrixBase<DeriveduEE> & uEE,
    const bool construct_TTi,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,
    std::vector<std::vector<std::vector<TTiIndex> > > & TTi);
}

#ifndef IGL_STATIC_LIBRARY
#  include "triangle_triangle_adjacency.cpp"
#endif

#endif
