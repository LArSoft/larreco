#include "Utilities/AssociationUtil.h"
#include "AnalysisBase/Calorimetry.h"
#include "RecoBase/Track.h"
#include "RecoBase/Seed.h"

#ifndef BEZIERTRACK_h
#define BEZIERTRACK_h

class TVector3;

namespace recob
{
  class Hit;
  class SpacePoint;
}

namespace calo
{
  class CalorimetryAlg;
}

namespace trkf {
  
  class BezierTrack: public recob::Track
  {
  public:
    BezierTrack();
    BezierTrack(const recob::Track& btb);
    BezierTrack(std::vector<TVector3> Pos, 
		std::vector<TVector3> Dir, 
		std::vector<std::vector<double> > dQdx,
		int const& id);
    BezierTrack(std::vector<recob::Seed> const& );
     
    ~BezierTrack();

    int NSegments()                                     const;

    double GetLength()                                  const;
    double GetRMSCurvature()                            const;

    double GetTotalCharge( unsigned int View )          const;
    double GetViewdQdx(    unsigned int View )          const; 

    TVector3 GetTrackPointV     (  double s )           const;
    TVector3 GetTrackDirectionV (  double s )           const;
    void   GetTrackPoint    (  double s, double* xyz )  const;
    void   GetTrackDirection(  double s, double* xyz )  const;  

    void   FillTrackVectors( std::vector<TVector3>& xyzVector,
			     std::vector<TVector3>& dirVector,
			     double const ds=0.1 ) const;
    
    double GetCurvature(double s)                       const;
    double GetdQdx(double s, unsigned int View)         const ;
    
    void   GetProjectedPointUVWX( double s, double* uvw, double * x,  int c, int t )    const;  
    void   GetProjectedPointUVWT( double s, double* uvw, double * ticks, int c, int t ) const;  

    void   GetClosestApproach( recob::Hit const&  hit,            double &s,  double& Distance) const;
    void   GetClosestApproach( art::Ptr<recob::Hit> const& hit,   double &s,  double& Distance) const;
    void   GetClosestApproach( recob::SpacePoint * sp,            double &s,  double& Distance) const;
    void   GetClosestApproach( recob::Seed const& seed,           double &s,  double& Distance) const;
    void   GetClosestApproach( TVector3 vec,                      double &s,  double& Distance) const;
    void   GetClosestApproach( uint32_t w, int p, int t, int c, float x, double& s,  double& Distance) const;

    
    double GetTrackPitch( geo::View_t view, double s, double WirePitch, unsigned int c=0, unsigned int t=0);

    void   GetClosestApproaches( art::PtrVector<recob::Hit> const& hits,     std::vector<double>& s, std::vector<double>& Distances) const;
    void   GetClosestApproaches( std::vector< art::Ptr<recob::Hit> > const& hits,     std::vector<double>& s, std::vector<double>& Distances) const;

    
    void   CalculatedQdx(art::PtrVector<recob::Hit>const &);   
    void   CalculatedQdx(art::PtrVector<recob::Hit>const &, std::vector<double> const & SValues);   
    
    anab::Calorimetry GetCalorimetryObject(std::vector<art::Ptr<recob::Hit> > const & Hits, geo::SigType_t sigtype, calo::CalorimetryAlg const& );
    anab::Calorimetry GetCalorimetryObject(art::PtrVector<recob::Hit > const & Hits, geo::SigType_t sigtype, calo::CalorimetryAlg const& );
    anab::Calorimetry GetCalorimetryObject(art::PtrVector<recob::Hit> const & Hits, geo::View_t view, calo::CalorimetryAlg const& calalg);


    BezierTrack Reverse();
    
    std::vector<recob::SpacePoint> GetSpacePointTrajectory(int N);
 
    BezierTrack GetPartialTrack(double LowS, double HighS);
    std::unique_ptr<recob::Track> GetBaseTrack();
    
    std::vector<recob::Seed> GetSeedVector();
    
    int WhichSegment(double S) const;
    
    
  private:
 
    void CalculateSegments();    
    
    void FillTrajectoryVectors();
    void FillSeedVector();
    
    std::vector<double>       fSegmentLength;
    std::vector<double>       fCumulativeLength;

    std::vector<recob::Seed>  fSeedCollection;
    
    double  fTrackLength;
    int     fBezierResolution;
    
  };

  class HitPtrVec
  {
  public:
    HitPtrVec() {}
    ~HitPtrVec() {}
    std::vector<art::Ptr<recob::Hit> > Hits;
  };

}

#endif
