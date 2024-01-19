/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
/** @class VertexListRefiner VertexListRefiner.h
 *
 *  Make a subselection of a track list
 *
 *  @author Wouter Hulsbergen
 *  @date   05/01/2010
 */

#include "Event/RecVertex.h"
#include "LHCbAlgs/Transformer.h"
#include "TrackKernel/TrackPredicates.h"
#include <limits>
#include <string>

class VertexListRefiner
    : public LHCb::Algorithm::Transformer<LHCb::RecVertex::Selection( const LHCb::RecVertex::Range& )> {
public:
  VertexListRefiner( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, KeyValue{"InputLocation", {}}, KeyValue{"OutputLocation", {}} ) {}

  LHCb::RecVertex::Selection operator()( const LHCb::RecVertex::Range& verticesin ) const override {
    LHCb::RecVertex::Selection verticesselected;
    m_PVtotal += verticesin.size();

    std::copy_if( verticesin.begin(), verticesin.end(), std::back_inserter( verticesselected ),
                  [&]( const LHCb::RecVertex* vertex ) {
                    // check vertrex chi2
                    if ( m_maxChi2PerDoF > 0 && vertex->chi2PerDoF() > m_maxChi2PerDoF ) return false;

                    // check if vertex in fiducial region
                    bool inFiducialRegion = true;
                    if ( inFiducialRegion && m_minX < m_maxX )
                      inFiducialRegion = m_minX < vertex->position().x() && vertex->position().x() < m_maxX;
                    if ( inFiducialRegion && m_minY < m_maxY )
                      inFiducialRegion = m_minY < vertex->position().y() && vertex->position().y() < m_maxY;
                    if ( inFiducialRegion && m_minZ < m_maxZ )
                      inFiducialRegion = m_minZ < vertex->position().z() && vertex->position().z() < m_maxZ;
                    if ( inFiducialRegion == m_vetoFiducialRegion ) return false;

                    // unfortunately stl doesn't work with the smartrefs in
                    // vertex. furthermore, when reading a dst, track pointers can be
                    // zero.

                    std::vector<const LHCb::Track*> tracks;
                    tracks.reserve( vertex->tracks().size() );
                    std::copy_if( vertex->tracks().begin(), vertex->tracks().end(), std::back_inserter( tracks ),
                                  []( const LHCb::Track* t ) { return t != nullptr; } );
                    bool accept = int( tracks.size() ) >= m_minNumTracks && int( tracks.size() ) <= m_maxNumTracks;
                    accept      = accept && ( m_minNumLongTracks == 0 ||
                                         std::count_if( tracks.begin(), tracks.end(),
                                                        TrackPredicates::Type( LHCb::Track::Types::Long ) ) >=
                                             m_minNumLongTracks );

                    if ( accept && ( m_minNumBackwardTracks > 0 || m_minNumForwardTracks > 0 ) ) {
                      int numback    = std::count_if( tracks.begin(), tracks.end(),
                                                   TrackPredicates::Type( LHCb::Track::Types::VeloBackward ) );
                      int numforward = tracks.size() - numback;
                      accept         = numback >= m_minNumBackwardTracks && numforward >= m_minNumForwardTracks;
                    }

                    return accept;
                  } );
    m_PVpassed += verticesselected.size();
    return verticesselected;
  }

private:
  Gaudi::Property<int>                       m_minNumTracks{this, "MinNumTracks", 0};
  Gaudi::Property<int>                       m_maxNumTracks{this, "MaxNumTracks", std::numeric_limits<int>::max()};
  Gaudi::Property<int>                       m_minNumBackwardTracks{this, "MinNumBackwardTracks", 0};
  Gaudi::Property<int>                       m_minNumForwardTracks{this, "MinNumForwardTracks", 0};
  Gaudi::Property<int>                       m_minNumLongTracks{this, "MinNumLongTracks", 0};
  Gaudi::Property<double>                    m_maxChi2PerDoF{this, "MaxChi2PerDoF", -1};
  Gaudi::Property<double>                    m_minX{this, "MinX", 1};
  Gaudi::Property<double>                    m_maxX{this, "MaxX", -1};
  Gaudi::Property<double>                    m_minY{this, "MinY", 1};
  Gaudi::Property<double>                    m_maxY{this, "MaxY", -1};
  Gaudi::Property<double>                    m_minZ{this, "MinZ", 1};
  Gaudi::Property<double>                    m_maxZ{this, "MaxZ", -1};
  Gaudi::Property<bool>                      m_vetoFiducialRegion{this, "VetoFiducialRegion", false};
  mutable Gaudi::Accumulators::StatCounter<> m_PVtotal{this, "#total"};
  mutable Gaudi::Accumulators::StatCounter<> m_PVpassed{this, "#passed"};
};

DECLARE_COMPONENT( VertexListRefiner )
