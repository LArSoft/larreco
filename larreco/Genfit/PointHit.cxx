// This Class' Header ------------------
#include "Genfit/PointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "Genfit/GeaneTrackRep2.h"
#include "Genfit/RKTrackRep.h"
#include "Genfit/SlTrackRep.h"
#include "Genfit/GFDetPlane.h"
#include "TRandom.h"
#include "math.h"

// Class Member definitions -----------

//ClassImp(PointHit)


genf::PointHit::~PointHit()
{}

genf::PointHit::PointHit()
  : SpacepointRecoHit(NparHitRep)
{}

genf::PointHit::PointHit(TVector3 point,double res)
  : SpacepointRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  fHitCov[2][2] = res*res;
  GFDetPlane d;

  fHitCoord[0][0] = point.X();
  fHitCoord[1][0] = point.Y();
  fHitCoord[2][0] = point.Z();

}


genf::PointHit::PointHit(TVector3 point, std::vector<double> &res)
  : SpacepointRecoHit(NparHitRep){

  //assert (res.size()==4);

  fHitCov[0][0] = res.at(0);
  fHitCov[1][1] = res.at(1);
  fHitCov[2][2] = res.at(3);
  fHitCov[1][2] = res.at(2); // yz cov element.
  fHitCov[2][1] = res.at(2); // yz cov element.

  GFDetPlane d;

  fHitCoord[0][0] = point.X();
  fHitCoord[1][0] = point.Y();
  fHitCoord[2][0] = point.Z();
}

genf::GFAbsRecoHit* genf::PointHit::clone(){
  return new PointHit(*this);
}


TMatrixT<Double_t>
genf::PointHit::getHMatrix(const GFAbsTrackRep* stateVector, const Double_t& betac, const Double_t& dist)
{
  if ((dynamic_cast<const genf::GeaneTrackRep2*>(stateVector) != NULL) ||
      (dynamic_cast<const genf::RKTrackRep*>(stateVector) != NULL)){
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);


    //   TMatrixT<Double_t> HMatrix(2,6);
    // WTF? Making this 2,5 not 2,6. EC, 3-Jan-2011.
    //    TMatrixT<Double_t> HMatrix(2,5);
    TMatrixT<Double_t> HMatrix(5,5);

    Double_t C = 0.0136/betac*sqrt(dist/14.0)*(1+0.038*log(dist/14.0)); // EC, 2-Jan-2012.

    HMatrix[0][0] = C; 
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;
    HMatrix[0][4] = 0.;

    HMatrix[1][0] = 0.; 
    HMatrix[1][1] = 1.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 0.;

    HMatrix[2][0] = 0.; 
    HMatrix[2][1] = 0.;
    HMatrix[2][2] = 1.;
    HMatrix[2][3] = 0.;
    HMatrix[2][4] = 0.;

    HMatrix[3][0] = 0.; 
    HMatrix[3][1] = 0.;
    HMatrix[3][2] = 0.;
    HMatrix[3][3] = 1.;
    HMatrix[3][4] = 0.;

    HMatrix[4][0] = 0.; 
    HMatrix[4][1] = 0.;
    HMatrix[4][2] = 0.;
    HMatrix[4][3] = 0.;
    HMatrix[4][4] = 1.;


    return HMatrix;
  }
  else if(dynamic_cast<const genf::SlTrackRep*>(stateVector)){
    TMatrixT<Double_t> HMatrix(2,4);
    
    HMatrix[0][0] = 1.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;

    
    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 1.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    return HMatrix;
  }

 else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
   throw;
 }
 
}

TMatrixT<Double_t>
genf::PointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const genf::GeaneTrackRep2*>(stateVector) != NULL) ||
      (dynamic_cast<const genf::RKTrackRep*>(stateVector) != NULL)){
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);


    //   TMatrixT<Double_t> HMatrix(2,6);
    // WTF? Making this 2,5 not 2,6. EC, 3-Jan-2011.
    TMatrixT<Double_t> HMatrix(2,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    //HMatrix[0][5] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;
    //HMatrix[1][5] = 0.;

    return HMatrix;
  }
  else if(dynamic_cast<const genf::SlTrackRep*>(stateVector)){
    TMatrixT<Double_t> HMatrix(2,4);
    
    HMatrix[0][0] = 1.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;

    
    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 1.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    return HMatrix;
  }

 else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
   throw;
 }
 
}


