/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include "Event/ZipUtils.h"
#include <limits>

/** @file  ChildRelation.h
 *  @brief Types and functions for associating physics objects each other.
 */

/** @class ChildRelation
 *  @brief Stores a link to a child object.
 */
struct ChildRelation {
  ChildRelation() = default;
  ChildRelation( std::size_t index, Zipping::ZipFamilyNumber family ) : m_index{index}, m_family{family} {}
  std::size_t              index() const { return m_index; }
  Zipping::ZipFamilyNumber zipIdentifier() const { return m_family; }
  void                     setIndex( std::size_t index ) { m_index = index; }
  void                     setZipIdentifier( Zipping::ZipFamilyNumber family ) { m_family = family; }

private:
  std::size_t              m_index{std::numeric_limits<std::size_t>::max()};
  Zipping::ZipFamilyNumber m_family{Zipping::generateZipIdentifier()};
};
