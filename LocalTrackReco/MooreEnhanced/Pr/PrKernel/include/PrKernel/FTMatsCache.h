/*****************************************************************************\
* (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

#include "FTDAQ/FTInfo.h"
#include "FTDet/DeFTDetector.h"
#include <array>

namespace FTMatsCache {

  const std::string Location = "AlgorithmSpecific-FTMatsCache";

  struct MatsCache {
    /**
     * partial SoA cache for mats, reserve enough (here 4096 which is more than enough)
     * space for all mats ( all mats should be less than 2 * 8 mats * 12 modules * 12 layers)
     */
    std::array<float, LHCb::Detector::FT::maxNumberMats>                  dxdy{};
    std::array<float, LHCb::Detector::FT::maxNumberMats>                  dzdy{};
    std::array<float, LHCb::Detector::FT::maxNumberMats>                  globaldy{};
    std::array<ROOT::Math::XYZPointF, LHCb::Detector::FT::maxNumberMats>  mirrorPoint{};
    std::array<ROOT::Math::XYZVectorF, LHCb::Detector::FT::maxNumberMats> ddx{};

    float uBegin{};
    float halfChannelPitch{};
    float dieGap{};
    float sipmPitch{};

    std::array<std::vector<double>, LHCb::Detector::FT::maxNumberMats> matContractionParameterVector{};

    MatsCache(){}; // Needed for DD4HEP
    MatsCache( const DeFT& ftDet ) {
      const auto first_mat = ftDet.firstMat();
      // This parameters are constant accross all mats:
#ifdef USE_DD4HEP
      this->dieGap           = first_mat.dieGap();
      this->sipmPitch        = first_mat.sipmPitch();
      this->uBegin           = first_mat.uBegin();
      this->halfChannelPitch = first_mat.halfChannelPitch();
#else
      this->dieGap           = first_mat->dieGap();
      this->sipmPitch        = first_mat->sipmPitch();
      this->uBegin           = first_mat->uBegin();
      this->halfChannelPitch = first_mat->halfChannelPitch();
#endif
      auto func = [this]( const DeFTMat& mat ) {
        assert( this->dieGap == mat.dieGap() && "Unexpected difference in dieGap" );
        assert( this->sipmPitch == mat.sipmPitch() && "Unexpected difference in sipmPitch" );
        assert( this->uBegin == mat.uBegin() && "Unexpected difference in uBegin" );
        assert( this->halfChannelPitch == mat.halfChannelPitch() && "Unexpected difference in halfChannelPitch" );
        const auto index = mat.elementID().globalMatID();
        // FIXME
        this->mirrorPoint[index]                   = mat.mirrorPoint();
        this->ddx[index]                           = mat.ddx();
        this->dxdy[index]                          = mat.dxdy();
        this->dzdy[index]                          = mat.dzdy();
        this->globaldy[index]                      = mat.globaldy();
        this->matContractionParameterVector[index] = mat.getmatContractionParameterVector();
      };
      ftDet.applyToAllMats( func );
    };
  };
} // namespace FTMatsCache
