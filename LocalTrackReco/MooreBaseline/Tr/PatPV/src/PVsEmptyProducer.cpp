/*****************************************************************************\
* (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/PrimaryVertices.h"
#include "Event/RecVertex.h"
#include "LHCbAlgs/EmptyProducer.h"

/** @class PVsEmptyProducer
 * @brief dummy producer of an empty container of PVs
 */

DECLARE_COMPONENT_WITH_ID( EmptyProducer<LHCb::RecVertices>, "RecVertexEmptyProducer" )
DECLARE_COMPONENT_WITH_ID( EmptyProducer<LHCb::Event::PV::PrimaryVertexContainer>, "PVsEmptyProducer" )
