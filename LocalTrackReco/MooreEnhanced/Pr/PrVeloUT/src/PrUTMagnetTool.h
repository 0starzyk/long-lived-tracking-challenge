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
#pragma once
#include <DetDesc/GenericConditionAccessorHolder.h>
#include <GaudiAlg/GaudiTool.h>
#include <Kernel/LUTForFunction.h>
#include <Magnet/DeMagnet.h>
#include <PrKernel/IPrUTMagnetTool.h>
#include <TrackInterfaces/ITrackExtrapolator.h>
#include <UTDet/DeUTDetector.h>
#include <limits>

/** @class PrUTMagnetTool PrUTMagnetTool.h
 *
 *  Tool used by the upstream pattern recognition to obtain integrals
 *  of the B-field as a function of parameters of a Velo track
 *
 *  @author Mariusz Witek
 *  @date   2006-09-25
 *  @update for A-Team framework 2007-08-20 SHM
 *
 */

namespace LHCb::Pr {

  class UTMagnetTool final : public LHCb::DetDesc::ConditionAccessorHolder<extends<GaudiTool, IPrUTMagnetTool>> {
  public:
    /// Standard constructor
    using ConditionAccessorHolder::ConditionAccessorHolder;

    StatusCode initialize() override;

    float bdlIntegral( float ySlopeVelo, float zOrigin, float zVelo ) const override {
      auto& c = cache();
      return c.noField ? s_bdlIntegral_NoB : c.lutBdl->getInterpolatedValue( {ySlopeVelo, zOrigin, zVelo} );
    }
    float zBdlMiddle( float ySlopeVelo, float zOrigin, float zVelo ) const override {
      auto& c = cache();
      return c.noField ? s_zBdlMiddle_NoB : c.lutZHalfBdl->getInterpolatedValue( {ySlopeVelo, zOrigin, zVelo} );
    }
    float dist2mom( float ySlope ) const override {
      auto& c = cache();
      return c.noField ? s_averageDist2mom_NoB : c.lutDxToMom->getValue( {ySlope} );
    }

    //=========================================================================
    // z middle of UT
    //=========================================================================
    float zMidUT() const override { return cache().zCenterUT; }
    //=========================================================================
    // z middle of B field betweenz=0 and zMidUT
    //=========================================================================
    float zMidField() const override { return cache().zMidField; }
    //=========================================================================
    // averageDist2mom
    //=========================================================================
    float averageDist2mom() const override { return cache().dist2mom; }

    const LUTForFunction<3, 30>&      DxLayTable() const override { return *cache().lutDxLay; }
    const LUTForFunction<30, 10, 10>& BdlTable() const override { return *cache().lutBdl; }

    const std::string& cacheLocation() const override { return m_cache.key(); }

  private:
    void  prepareBdlTables( const DeMagnet& magnet, const DeUTDetector& utdet, Cache& c ) const;
    void  prepareDeflectionTables( IGeometryInfo const& geometry, const DeMagnet& magnet, Cache& c ) const;
    Cache makeCache( const LHCb::Detector::DeLHCb& lhcb, const DeMagnet& magnet, const DeUTDetector& utdet ) const;

    ConditionAccessor<Cache> m_cache{this, "PrUTMagnetTool-Cache-" + name()};
    const Cache&             cache() const;

    // set parameters for no field run
    constexpr static float s_bdlIntegral_NoB     = 0.1;
    constexpr static float s_zBdlMiddle_NoB      = 1900.;
    constexpr static float s_zMidField_NoB       = 1900.;
    constexpr static float s_averageDist2mom_NoB = 0.00004;

    // Retrieve extrapolators
    PublicToolHandle<ITrackExtrapolator> m_linear{this, "LinearExtrapolator", "TrackLinearExtrapolator"};
    PublicToolHandle<ITrackExtrapolator> m_parabolic{this, "Extrapolator", "TrackRungeKuttaExtrapolator"};
  };
} // namespace LHCb::Pr
