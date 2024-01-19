/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include <Event/PrHits.h>
#include <Gaudi/Accumulators/Histogram.h>
#include <LHCbAlgs/Consumer.h>
class MonitorDetectorCorrelationsVeloSciFi
    : public LHCb::Algorithm::Consumer<void( LHCb::Pr::VP::Hits const&, LHCb::Pr::FT::Hits const& )> {
public:
  MonitorDetectorCorrelationsVeloSciFi( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator, {KeyValue{"VeloHits", ""}, KeyValue{"SciFiHits", ""}} ){};

  void operator()( LHCb::Pr::VP::Hits const& velo_hits, LHCb::Pr::FT::Hits const& scifi_hits ) const override;
  mutable Gaudi::Accumulators::Histogram<2> m_velo_scifi_hits_correlation{
      this,
      "VeloScifiHitsCorrelation",
      "VeloScifiHitsCorrelation",
      {{100, 0, 6000, "Velo Hits"}, {100, 0, 10000, "ScifiHits"}}};
};
DECLARE_COMPONENT( MonitorDetectorCorrelationsVeloSciFi )

void MonitorDetectorCorrelationsVeloSciFi::operator()( LHCb::Pr::VP::Hits const& velo_hits,
                                                       LHCb::Pr::FT::Hits const& scifi_hits ) const {
  ++m_velo_scifi_hits_correlation[{velo_hits.size(), scifi_hits.size()}];
}

class MonitorDetectorCorrelations
    : public LHCb::Algorithm::Consumer<void( LHCb::Pr::VP::Hits const&, LHCb::Pr::FT::Hits const&,
                                             MuonHitContainer const& )> {
public:
  MonitorDetectorCorrelations( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"VeloHits", ""}, KeyValue{"SciFiHits", ""}, KeyValue{"MuonHits", ""}} ){};

  void operator()( LHCb::Pr::VP::Hits const& velo_hits, LHCb::Pr::FT::Hits const& scifi_hits,
                   MuonHitContainer const& muon_hits ) const override;
  mutable Gaudi::Accumulators::Histogram<2> m_velo_scifi_hits_correlation{
      this,
      "VeloScifiHitsCorrelation",
      "VeloScifiHitsCorrelation",
      {{200, 0, 6000, "Velo Hits"}, {200, 0, 10000, "ScifiHits"}}};
  mutable Gaudi::Accumulators::Histogram<2> m_velo_muon_hits_correlation{
      this,
      "VeloMuonHitsCorrelation",
      "VeloMuonHitsCorrelation",
      {{200, 0, 6000, "Velo Hits"}, {200, 0, 600, "MuonHits"}}};
  mutable Gaudi::Accumulators::Histogram<2> m_scifi_muon_hits_correlation{
      this,
      "ScifiMuonHitsCorrelation",
      "ScifiMuonHitsCorrelation",
      {{200, 0, 10000, "SciFi Hits"}, {200, 0, 600, "MuonHits"}}};
};
DECLARE_COMPONENT( MonitorDetectorCorrelations )

void MonitorDetectorCorrelations::operator()( LHCb::Pr::VP::Hits const& velo_hits, LHCb::Pr::FT::Hits const& scifi_hits,
                                              MuonHitContainer const& muon_hits ) const {
  ++m_velo_scifi_hits_correlation[{velo_hits.size(), scifi_hits.size()}];
  auto n_muon_hits =
      muon_hits.hits( 0 ).size() + muon_hits.hits( 1 ).size() + muon_hits.hits( 2 ).size() + muon_hits.hits( 3 ).size();
  ++m_velo_muon_hits_correlation[{velo_hits.size(), n_muon_hits}];
  ++m_scifi_muon_hits_correlation[{scifi_hits.size(), n_muon_hits}];
}