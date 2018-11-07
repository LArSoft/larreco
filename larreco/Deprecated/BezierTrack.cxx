#include "larreco/Deprecated/BezierTrack.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierCurveHelper.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "cetlib_except/exception.h"

#include <utility> // std::forward()
#include <cassert>


namespace trkf {

  //----------------------------------------------------------------------
  // Constructor from Base object (for analysis)
  //
  BezierTrack::BezierTrack(int id, const recob::Trajectory& traj)
    : fTraj(traj)
    // Clang: private field 'fID' is not used
  //  , fID(id)
  {
    fBezierResolution=1000;
    CalculateSegments();
  }


  //----------------------------------------------------------------------
  BezierTrack::BezierTrack(const recob::Track& track)
    : BezierTrack(track.ID(), track.Trajectory().Trajectory()) {}
  
  //----------------------------------------------------------------------
  // Constructor from seed vector (for production)
  //
  BezierTrack::BezierTrack(std::vector<recob::Seed> const& SeedCol )
    : fTraj()
  {
    if(SeedCol.size()>0)
      {
        double Pt[3], Dir[3], PtErr[3], DirErr[3];
        double FirstPt[3];
        double FirstDir[3];
        SeedCol.at(0).GetPoint(     Pt,  PtErr  );
        SeedCol.at(0).GetDirection( Dir, DirErr );
        for(int i=0; i!=3; ++i)
          {
            double BigN=1000;
            FirstPt[i]=Pt[i]-((1.+1./BigN) * Dir[i]);
            FirstDir[i]=Dir[i]/BigN;
          }
        fSeedCollection.push_back(recob::Seed(FirstPt, FirstDir, PtErr, DirErr));
      }
    for(size_t i=0; i!=SeedCol.size(); ++i)
      {
        fSeedCollection.push_back(SeedCol.at(i));
      }
    if(SeedCol.size()>0)
      {
        double Pt[3], Dir[3], PtErr[3], DirErr[3];
        double LastPt[3], LastDir[3];
        SeedCol.at(SeedCol.size()-1).GetPoint(     Pt,  PtErr  );
        SeedCol.at(SeedCol.size()-1).GetDirection( Dir, DirErr );
        for(int i=0; i!=3; ++i)
          {
            double BigN=1000;
            LastPt[i]=Pt[i]+((1.+1./BigN) * Dir[i]);
            LastDir[i]=Dir[i]/BigN;
          }
        fSeedCollection.push_back(recob::Seed(LastPt, LastDir, PtErr, DirErr));
      }
   
    CalculateSegments();
    fBezierResolution=1000;
  }


  //----------------------------------------------------------------------
  // Constructor from track coordinates
  //
  BezierTrack::BezierTrack(std::vector<TVector3> const& Pos,
                           std::vector<TVector3> const& Dir,
                           std::vector<std::vector<double> > const& dQdx,
                           int id)
    : fTraj(geo::vect::convertCollToPoint(Pos), geo::vect::convertCollToVector(Dir), false)
    // Clang: private field 'fID' is not used
  //  , fID(id)
    , fdQdx(dQdx)
  {
    fBezierResolution=1000;
    CalculateSegments();
  }


  //----------------------------------------------------------------------
  //  Return the seed collection (with spurious end seeds cut off
  //
  std::vector<recob::Seed> BezierTrack::GetSeedVector()
  {
    std::vector<recob::Seed> ReturnVector=fSeedCollection;
    if(fSeedCollection.size()>2)
      {
        ReturnVector.erase(ReturnVector.begin());
        ReturnVector.pop_back();
        return ReturnVector;
      }
    else
      return std::vector<recob::Seed>();
    
  }



  //----------------------------------------------------------------------
  // Get the track pitch at some s for a particular view
  //

  double BezierTrack::GetTrackPitch(geo::View_t view, double s, double /* WirePitch */, unsigned int c, unsigned int t)
  {
    static std::map<geo::View_t, bool>          DoneCalc;
    static std::map<geo::View_t, TVector3>      PitchVecs;

    if(!(DoneCalc[view]))
      {
        unsigned int p;
        unsigned int pTop;

        art::ServiceHandle<geo::Geometry> geom;

        pTop=geom->Cryostat(c).TPC(t).Nplanes();
        for(p=0; p!=pTop; ++p)
          if(geom->Cryostat(c).TPC(t).Plane(p).View() == view) break;

        double WireEnd1[3], WireEnd2[3];

        geom->WireEndPoints(c, t, p, 0, WireEnd1, WireEnd2);

        double WirePitch = geom->WirePitch(view);

        TVector3 WireVec;
        for(size_t i=0; i!=3; ++i)
          WireVec[i]= WireEnd2[i] -WireEnd1[i];
        PitchVecs[view] = (1.0 / WirePitch) * WireVec.Cross(TVector3(1,0,0)).Unit();

        DoneCalc[view]=true;
      }

    TVector3 TrackDir = GetTrackDirectionV(s).Unit();

    return 1. / (PitchVecs[view].Dot(TrackDir));
  }


  //------------------------------------------------------
  // Build an anab::Calorimetry object for this track (wrapper for
  //  vector<art::Ptr>

  anab::Calorimetry BezierTrack::GetCalorimetryObject(std::vector<art::Ptr<recob::Hit> > const& Hits, geo::SigType_t sigtype, calo::CalorimetryAlg const& calalg)
  {
    art::PtrVector<recob::Hit> Hitv;
    for(size_t i=0; i!=Hitv.size(); ++i)
      {
        Hitv.push_back(Hits.at(i));
      }
    Hitv.clear();
    
    return GetCalorimetryObject(Hitv, sigtype, calalg);
  }

  //------------------------------------------------------
  // Build an anab::Calorimetry object for this track

  anab::Calorimetry BezierTrack::GetCalorimetryObject(art::PtrVector<recob::Hit > const & Hits, geo::SigType_t sigtype, calo::CalorimetryAlg const& calalg)
  {

    
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<float>  dEdx;
    std::vector<float>  dQdx;
    std::vector<float>  resRange;
    std::vector<float>  deadwire;
    float Range = GetLength();
    std::vector<float>  TrkPitch;



    double s, d;

    bool do_now=true;
    geo::PlaneID planeID;

    if(Hits.size()>0)
      {
        for(size_t i=0; i!=Hits.size(); ++i)
          {
            if(Hits.at(i)->SignalType()==sigtype)
              {
                geo::View_t view = Hits.at(i)->View();
                double WirePitch = geom->WirePitch(view);
                GetClosestApproach(Hits.at(i),s,d);
                double ThisPitch = GetTrackPitch( view, s, WirePitch );
                TrkPitch.push_back(ThisPitch);
                resRange.push_back(Range - s*Range);
                dQdx.push_back(Hits.at(i)->Integral()/ThisPitch);
                dEdx.push_back( calalg.dEdx_AMP(Hits.at(i), ThisPitch) );

                if(do_now) {
                  planeID = Hits.at(i)->WireID().planeID();
                  do_now=false;
                }
              }
            
          }
      }


    double KineticEnergy = 0.;

    for(size_t i=0; i!=dEdx.size(); ++i)
      KineticEnergy+=dEdx.at(i)*TrkPitch.at(i);

    return anab::Calorimetry(KineticEnergy, dEdx, dQdx, resRange, deadwire, Range, TrkPitch,planeID);


  }


  //------------------------------------------------------
  // Build an anab::Calorimetry object for this track

  anab::Calorimetry BezierTrack::GetCalorimetryObject(art::PtrVector<recob::Hit > const & Hits, geo::View_t view, calo::CalorimetryAlg const& calalg)
  {

   
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<float>  dEdx;
    std::vector<float>  dQdx;
    std::vector<float>  resRange;
    std::vector<float>  deadwire;
    float Range = GetLength();
    std::vector<float>  TrkPitch;

    double WirePitch = geom->WirePitch(view);

    double KineticEnergy = 0;                
 
    double s, d;

    bool do_now=true;
    geo::PlaneID planeID;

    if(Hits.size()>0)
      {
        for(size_t i=0; i!=Hits.size(); ++i)
          {
            if(Hits.at(i)->View()==view)
              {
                GetClosestApproach(Hits.at(i),s,d);
                double ThisPitch = GetTrackPitch( view, s, WirePitch );
                TrkPitch.push_back(ThisPitch);
                resRange.push_back( (1.-s) * Range);
                dQdx.push_back(Hits.at(i)->Integral()/ThisPitch);
                double ThisdEdx = calalg.dEdx_AMP(Hits.at(i), ThisPitch);
                dEdx.push_back( ThisdEdx);
                KineticEnergy+=ThisdEdx*ThisPitch;
                                
                if(do_now) {
                  planeID = Hits.at(i)->WireID().planeID();
                  do_now=false;
                }

              }
            
          }
      }


    return anab::Calorimetry(KineticEnergy, dEdx, dQdx, resRange, deadwire, Range, TrkPitch, planeID);

  }


  //----------------------------------------------------------------------
  // Given track points, fill the seed vector.  The seeds are then used
  // for geomtry calculations / bezier fitting
  //
  void BezierTrack::FillSeedVector()
  {
    int NSeg=NSegments();



    double Pt[3], Dir[3];
    for(int i=0; i!=NSeg; i++)
      {

        auto const& pos = fTraj.LocationAtPoint(i);
        Pt[0]=pos.X();
        Pt[1]=pos.Y();
        Pt[2]=pos.Z();

        auto const& dir = fTraj.DirectionAtPoint(i);
        Dir[0]=dir.X();
        Dir[1]=dir.Y();
        Dir[2]=dir.Z();

        fSeedCollection.push_back(recob::Seed(Pt,Dir));
      }

  }




  //----------------------------------------------------------------------
  void BezierTrack::FillTrajectoryVectors()
  {
    art::ServiceHandle<geo::Geometry> geom;
    int NPlanes = geom->Nplanes();
    fdQdx.resize(NPlanes) ;
    double Pt[3], Dir[3], ErrPt[3], ErrDir[3];
    recob::Trajectory::Positions_t positions;
    recob::Trajectory::Momenta_t directions;
    for(std::vector<recob::Seed>::const_iterator it=fSeedCollection.begin();
        it!=fSeedCollection.end(); ++it)
      {

        it->GetPoint(Pt,ErrPt);
        it->GetDirection(Dir, ErrDir);
        recob::Track::Point_t Point(Pt[0],Pt[1],Pt[2]);
        recob::Track::Vector_t Direction(Dir[0],Dir[1],Dir[2]);

        for(int i=0; i!=NPlanes; ++i)
          fdQdx.at(i).push_back(0);

        positions.push_back(Point);
        directions.push_back(Direction);
      }
    
    fTraj
      = recob::Trajectory(std::move(positions), std::move(directions), false);
    
  }


  //----------------------------------------------------------------------
  // Calculate the lengths of each bezier segment and fill useful
  // calculational data members
  //

  void BezierTrack::CalculateSegments()
  {

    if(fSeedCollection.size()==0)
      {
        if(NSegments()!=0) FillSeedVector();
        else
          {
            throw cet::exception("no points in track")
              <<"CalculateSegments method of Bezier track called with no"
              <<" track information loaded.  You must fill track with "
              <<"poisition and direction data before calling this method."
              <<"\n";
          }
      }
    if(NSegments()==0)
      {
        if(fSeedCollection.size()!=0)  FillTrajectoryVectors();
      }

    BezierCurveHelper bhlp(100);

    int Segments = NSegments();

    fTrackLength=0;
    bool FirstSeg=true;
    for(int i=0; i!=Segments; ++i)
      {
        if(!FirstSeg)
          {
            float SegmentLength=bhlp.GetSegmentLength(fSeedCollection.at(i-1),fSeedCollection.at(i));

            fTrackLength+=SegmentLength;
            fSegmentLength.push_back(SegmentLength);
          }
        FirstSeg=false;
        fCumulativeLength.push_back(fTrackLength);
      }
  }


  //----------------------------------------------------------------------
  // Find the point which is fraction s along the track and return
  // as a double[3]
  //
  void BezierTrack::GetTrackPoint(double s, double * xyz) const
  {
    if((s>=1.)||(s<=0.))
      {
        // catch these easy floating point errors
        if((s>0.9999)&&(s<1.00001))
          {
            auto End1 = fTraj.End();
            xyz[0]=End1.X();
            xyz[1]=End1.Y();
            xyz[2]=End1.Z();

            return;
          }
        else if((s<0.0001)&&(s>-0.00001))
          {
            auto End0 = fTraj.Start();
            xyz[0]=End0.X();
            xyz[1]=End0.Y();
            xyz[2]=End0.Z();

            return;
          }

        // otherwise complain about the mistake
        else
          throw cet::exception("track point out of range")<<" s = "<<s <<" out of range \n";
      }
    else
      {
        BezierCurveHelper bhlp;
        for(unsigned int i=1; i!=fCumulativeLength.size(); i++)
          {
            if(  (   (fCumulativeLength.at(i-1) / fTrackLength) <=  s)
                 &&( (fCumulativeLength.at(i)   / fTrackLength) > s))
              {
                double locals = (s * fTrackLength - fCumulativeLength[i-1])/fSegmentLength[i-1];


                bhlp.GetBezierPointXYZ(fSeedCollection.at(i-1),fSeedCollection.at(i),locals, xyz);
              }

          }
      }
  }





  //----------------------------------------------------------------------
  //  Find the point which is fraction s along the track and get its
  //   projected point in the wire view, in system uvwx
  //

  void BezierTrack::GetProjectedPointUVWX(double s, double *uvw, double*x, int t=0, int c=0) const
  {

    art::ServiceHandle<geo::Geometry>   geo;

    double  xyz[3];

    // Get point in 3D space
    GetTrackPoint(s , xyz);

    int NPlanes=geo->Cryostat(c).TPC(t).Nplanes();

    for(int p=0; p!=NPlanes; p++)
      {        
        uvw[p]= geo->NearestWire(xyz,p,t,c);
      }
    x[0]=xyz[0];
  }






  //----------------------------------------------------------------------
  //  Find the point which is fraction s along the track and get its
  //   projected point in the wire view, in system uvwt
  //

  void BezierTrack::GetProjectedPointUVWT(double s, double *uvw, double*ticks, int t=0, int c=0) const
  {
    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    art::ServiceHandle<geo::Geometry>            geo;

    double  xyz[3];

    // Get point in 3D space
    GetTrackPoint(s , xyz);
    
    int NPlanes=geo->Cryostat(c).TPC(t).Nplanes();

    for(int p=0; p!=NPlanes; p++)
      {
        uvw[p]= geo->NearestWire(xyz,p,t,c);
        ticks[p]=det->ConvertXToTicks(xyz[0],p,t,c);
      }
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   - this fast version is for ennumerating which members of a hit
  //     collection are within d, how far each is and where the closest
  //     approach occurs.  This version is optimized for speed for this
  //     application
  void BezierTrack::GetClosestApproaches( art::PtrVector<recob::Hit> const & hits ,     std::vector<double>& s,  std::vector<double>& Distances) const
  {

    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    art::ServiceHandle<geo::Geometry>            geo;

    s.clear();
    Distances.clear();


    // Pull all the relevant information out of hits and into simple vectors

    std::vector<TVector3> HitEnd1s, HitEnd2s;
    std::vector<double> WireLengths;

    HitEnd1s.resize(hits.size());
    HitEnd2s.resize(hits.size());
    WireLengths.resize(hits.size());

    //unsigned int c1, t1, p1, w1;
    double End1[3], End2[3];

    size_t NHits = hits.size();

    for(size_t i=0; i!=NHits; i++)
      {
        Distances.push_back(10000);
        s.push_back(-1);

        geo::WireID hitWireID = hits.at(i)->WireID();

        geo->WireEndPoints(hitWireID.Cryostat,hitWireID.TPC,hitWireID.Plane,hitWireID.Wire,End1,End2);

        HitEnd1s.at(i)[0]= HitEnd2s.at(i)[0]= det->ConvertTicksToX(hits.at(i)->PeakTime(),hitWireID.Plane,hitWireID.TPC,hitWireID.Cryostat);
        HitEnd1s.at(i)[1]= End1[1];
        HitEnd2s.at(i)[1]= End2[1];
        HitEnd1s.at(i)[2]= End1[2];
        HitEnd2s.at(i)[2]= End2[2];

        WireLengths.at(i)=((HitEnd1s.at(i)-HitEnd2s.at(i)).Mag());
      }


    double iS;

    for(int ipt=0; ipt<=fBezierResolution; ++ipt)
      {
        iS=float(ipt)/fBezierResolution;
        TVector3 trackpt = GetTrackPointV(iS);

        for(size_t ihit=0; ihit!=NHits; ++ihit)
          {
            float d = ((trackpt-HitEnd1s.at(ihit)).Cross(trackpt-HitEnd2s.at(ihit))).Mag()/WireLengths.at(ihit);

            if(d<Distances.at(ihit))
              {
                Distances.at(ihit)=d;
                s.at(ihit)=iS;
              }
          }
      }


  }

  //--------------------------------------------------------

  void BezierTrack::GetClosestApproaches( std::vector< art::Ptr<recob::Hit> > const & hits,     std::vector<double>& s,  std::vector<double>& Distances) const
  {
    art::PtrVector<recob::Hit> hitv;
    for(size_t h = 0; h < hits.size(); ++h) hitv.push_back(hits[h]);
    this->GetClosestApproaches(hitv, s, Distances);
    hitv.clear();
  }




  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   and also the point where this occurs
  //

  void BezierTrack::GetClosestApproach( recob::Hit const & hit,       double& s,  double& Distance) const
  {
    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    
    art::ServiceHandle<geo::Geometry>            geo;

    //unsigned int c1, t1, p1, w1;

    double xyzend1[3], xyzend2[3];

    geo::WireID hitWireID = hit.WireID();
    geo->WireEndPoints(hitWireID.Cryostat,hitWireID.TPC,hitWireID.Plane,hitWireID.Wire,xyzend1,xyzend2);
    
    xyzend1[0] = xyzend2[0] = det->ConvertTicksToX(hit.PeakTime(),hitWireID.Plane,hitWireID.TPC,hitWireID.Cryostat);
    
    double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i<=fBezierResolution; ++i)
      {
        iS=float(i)/fBezierResolution;
        GetTrackPoint(iS, xyz);
        // calculate line to point distance in 3D
        TVector3 end1(xyzend1[0],xyzend1[1],xyzend1[2]);
        TVector3 end2(xyzend2[0],xyzend2[1],xyzend2[2]);
        TVector3 trackpt(xyz[0],xyz[1],xyz[2]);

        float d = ((trackpt-end1).Cross(trackpt-end2)).Mag()/(end2-end1).Mag();

        if(d<MinDistanceToPoint)
          {
            MinDistanceToPoint=d;
            MinS=iS;
          }
      }

    s = MinS;
    Distance = MinDistanceToPoint;

  }




  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   and also the point where this occurs
  //

  void BezierTrack::GetClosestApproach( art::Ptr<recob::Hit> const& hit,       double& s,  double& Distance) const
  {
    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    art::ServiceHandle<geo::Geometry>            geo;

    //unsigned int c1, t1, p1, w1;

    double xyzend1[3], xyzend2[3];

    geo::WireID hitWireID = hit->WireID();
    geo->WireEndPoints(hitWireID.Cryostat,hitWireID.TPC,hitWireID.Plane,hitWireID.Wire,xyzend1,xyzend2);

    xyzend1[0] = xyzend2[0] = det->ConvertTicksToX(hit->PeakTime(),hitWireID.Plane,hitWireID.TPC,hitWireID.Cryostat);

    double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i<=fBezierResolution; ++i)
      {
        iS=float(i)/fBezierResolution;
        GetTrackPoint(iS, xyz);
        // calculate line to point distance in 3D
        TVector3 end1(xyzend1[0],xyzend1[1],xyzend1[2]);
        TVector3 end2(xyzend2[0],xyzend2[1],xyzend2[2]);
        TVector3 trackpt(xyz[0],xyz[1],xyz[2]);

        float d = ((trackpt-end1).Cross(trackpt-end2)).Mag()/(end2-end1).Mag();

        if(d<MinDistanceToPoint)
          {
            MinDistanceToPoint=d;
            MinS=iS;
          }
      }

    s = MinS;
    Distance = MinDistanceToPoint;

  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   and also the point where this occurs
  //

  void BezierTrack::GetClosestApproach( uint32_t w, int p, int t, int c, float x, double& s,  double& Distance) const
  {
    //    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    art::ServiceHandle<geo::Geometry>            geo;
    
    
    double xyzend1[3], xyzend2[3];

    geo->WireEndPoints(c,t,p,w,xyzend1,xyzend2);

    xyzend1[0] = xyzend2[0] = x;

    double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i<=fBezierResolution; ++i)
      {
        iS=float(i)/fBezierResolution;
        GetTrackPoint(iS, xyz);
        // calculate line to point distance in 3D
        TVector3 end1(xyzend1[0],xyzend1[1],xyzend1[2]);
        TVector3 end2(xyzend2[0],xyzend2[1],xyzend2[2]);
        TVector3 trackpt(xyz[0],xyz[1],xyz[2]);

        float d = ((trackpt-end1).Cross(trackpt-end2)).Mag()/(end2-end1).Mag();

        if(d<MinDistanceToPoint)
          {
            MinDistanceToPoint=d;
            MinS=iS;
          }
      }

    s = MinS;
    Distance = MinDistanceToPoint;

  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given spacepoint
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( recob::SpacePoint* sp, double& s,  double& Distance) const
  {
    const double* xyz = sp->XYZ();
    TVector3 Vec(xyz[0],xyz[1],xyz[2]);
    GetClosestApproach(Vec, s, Distance);
  }

  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given seed
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( recob::Seed const& seed, double& s,  double& Distance) const
  {
    double SeedPt[3], SeedErr[3];
    seed.GetPoint(SeedPt,SeedErr);
    TVector3 Vec(SeedPt[0],SeedPt[1],SeedPt[2]);
    GetClosestApproach(Vec, s, Distance);
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given position
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( TVector3 vec,          double& s,  double& Distance) const
  {
    //    const detinfo::DetectorProperties* det = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    art::ServiceHandle<geo::Geometry>            geo;

    double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i<=fBezierResolution; ++i)
      {
        iS=float(i)/fBezierResolution;
        GetTrackPoint(iS, xyz);
        TVector3 trackpt(xyz[0],xyz[1],xyz[2]);

        float d = (vec-trackpt).Mag();

        if(d<MinDistanceToPoint)
          {
            MinDistanceToPoint=d;
            MinS=iS;
          }
      }

    s = MinS;
    Distance = MinDistanceToPoint;


  }



  //----------------------------------------------------------------------
  // Calculate the direction of the track by finding the difference
  //  in position between two points. Dir is normalized to 1
  //

  void BezierTrack::GetTrackDirection(double s, double * xyz) const
  {

    if((s<0.5/fBezierResolution)||(s>(1.-0.5/fBezierResolution)))
      {
        if     (( s < 0.5 / fBezierResolution ) && ( s > -0.00001 ) )
          s = 0.5 / fBezierResolution;
        else if((s>(1.-0.5/fBezierResolution)) && ( s < 1.00001)) s = 1.-0.5/fBezierResolution;
        else
          throw cet::exception("BezierTrack error: s out of range")<<
            " cannot query gradient within "<< 0.5/fBezierResolution<<
            " of track end.  You asked for s = "<<s <<
            ", which is out of range \n";
      }
    double xyz1[3], xyz2[3];
    GetTrackPoint(s - 0.5/fBezierResolution, xyz1);
    GetTrackPoint(s + 0.5/fBezierResolution, xyz2);

    double dx = pow(pow(xyz1[0]-xyz2[0],2)+
                    pow(xyz1[1]-xyz2[1],2)+
                    pow(xyz1[2]-xyz2[2],2),0.5);

    for(int i=0; i!=3; ++i)
      {
        xyz[i] = (xyz2[i]-xyz1[i])/dx;
      }

  }


  //----------------------------------------------------------------------
  // Friendly versions returning TVector3's
  //


  //----------------------------------------------------------------------
  TVector3 BezierTrack::GetTrackPointV(double s) const
  {
    double xyz[3];
    GetTrackPoint(s,xyz);
    TVector3 ReturnVec(xyz[0],xyz[1],xyz[2]);
    return ReturnVec;
  }

  //----------------------------------------------------------------------
  TVector3 BezierTrack::GetTrackDirectionV(double s) const
  {
    double xyz[3];
    GetTrackDirection(s,xyz);
    TVector3 ReturnVec(xyz[0],xyz[1],xyz[2]);
    return ReturnVec;
  }



  //----------------------------------------------------------------------
  // Method for finding rate at which the track direction changes
  //  (output is local dtheta/dx)

  double BezierTrack::GetCurvature(double s) const
  {
    if((s<1./fBezierResolution)||(s>(1.-1./fBezierResolution)))
      {
        throw cet::exception("BezierTrack error: s out of range")<<
          " cannot query curvature within "<< 1./fBezierResolution<<
          " of track end.  You asked for s = "<<s <<
          ", which is out of range \n";
      }

    TVector3 Pos1 = GetTrackPointV(s - 0.5/fBezierResolution);
    TVector3 Pos2 = GetTrackPointV(s );
    TVector3 Pos3 = GetTrackPointV(s + 0.5/fBezierResolution);

    double dx     = (1./fBezierResolution)*fTrackLength;
    double dtheta = (Pos3-Pos2).Angle(Pos2-Pos1);

    return (dtheta/dx);
  }


  //----------------------------------------------------------------------
  // Find the RMS curvature along the entire track
  //  (possible measure of multiple scattering)

  double BezierTrack::GetRMSCurvature() const
  {
    double RMS =0.;
    for(int i=1; i!=(fBezierResolution-1); ++i)
      {
        RMS += pow(GetCurvature( float(i)/fBezierResolution),2);
      }
    return (pow(RMS/(fBezierResolution-2),0.5));
  }



  //----------------------------------------------------------------------
  //  return the track length (already calculated, so easy!)
  //

  double BezierTrack::GetLength() const
  {
    return fTrackLength;
  }

  //--------------------------------------------

  BezierTrack BezierTrack::GetPartialTrack(double LowS, double HighS)
  {
    if(!((LowS>=0)&&(HighS<=1.0)&&(LowS<HighS)))
      {
        mf::LogError("BezierTrack")<<"Error in partial track calc - S out of range";
      }
   
    int LowSegment  = WhichSegment(LowS);
    int HighSegment = WhichSegment(HighS);
  
    TVector3 HighEnd = GetTrackPointV(HighS);
    TVector3 LowEnd  = GetTrackPointV(LowS);
    
    if(LowSegment!=HighSegment)
      {
        
        std::vector<recob::Seed> NewSeedCollection;
        for(int seg=LowSegment+1; seg<=HighSegment; ++seg)
          {
            NewSeedCollection.push_back(fSeedCollection.at(seg));
          }
        
        double PtLow[3],  DirLow[3], Err[3];
        double PtHigh[3], DirHigh[3];
        double LengthLow, LengthHigh;
        
        NewSeedCollection.at(0).GetPoint(PtLow,Err);
        NewSeedCollection.at(0).GetDirection(DirLow,Err);
        LengthLow = NewSeedCollection.at(0).GetLength();
        
        NewSeedCollection.at(NewSeedCollection.size()-1).GetPoint(PtHigh,Err);
        NewSeedCollection.at(NewSeedCollection.size()-1).GetDirection(DirHigh,Err);
        LengthHigh = NewSeedCollection.at(NewSeedCollection.size()-1).GetLength();
        
        double LowScale = fabs((TVector3(PtLow[0],PtLow[1],PtLow[2]) - LowEnd).Dot(TVector3(DirLow[0],DirLow[1],DirLow[2]).Unit()));
        
        double HighScale = fabs((TVector3(PtHigh[0],PtHigh[1],PtHigh[2]) - HighEnd).Dot(TVector3(DirHigh[0],DirHigh[1],DirHigh[2]).Unit()));
        

        for(size_t n=0; n!=3; ++n)
          {
            DirLow[n]  *= LowScale  / LengthLow;
            DirHigh[n] *= HighScale / LengthHigh;
            Err[n]      = 0;
          } 
        
        NewSeedCollection.at(0).SetDirection(DirLow, Err);
        NewSeedCollection.at(NewSeedCollection.size()-1).SetDirection(DirHigh, Err);
        return trkf::BezierTrack(NewSeedCollection);
      }
    else
      {
        double Pt[3], Dir[3], Err[3];
        for(size_t n=0; n!=3; ++n)
          {
            Pt[n]  = 0.5*(HighEnd[n] + LowEnd[n]);
            Dir[n] = 0.5*(HighEnd[n] - LowEnd[n]);
            Err[n] = 0;
          }
        recob::Seed OneSeed(Pt, Dir, Err, Err);
        std::vector<recob::Seed> NewSeedCollection;
        NewSeedCollection.push_back(OneSeed);
        
        return trkf::BezierTrack(NewSeedCollection);
      }
   
  }

  //----------------------------------------------------------------------
  //  Fill the dQdx vector for the track, based on a set of hits
  //    provided

  void BezierTrack::CalculatedQdx(art::PtrVector<recob::Hit> const &Hits)
  {
    fdQdx.clear();

    std::map<int, std::map<int, double> > hitmap;
    //        ^              ^       ^
    //       view            seg    charge

    for(size_t i=0; i!=Hits.size(); ++i)
      {
        double Distance, S;
        GetClosestApproach(Hits.at(i), S, Distance);

        (hitmap[Hits.at(i)->View()])[WhichSegment(S)] += Hits.at(i)->Integral();
      }

    int NSeg = NSegments();

    for(std::map<int,std::map<int,double> >::const_iterator itview=hitmap.begin();
        itview!=hitmap.end(); ++itview)
      {
        std::vector<double> ThisViewdQdx;
        ThisViewdQdx.resize(NSeg);
        for(std::map<int,double>::const_iterator itseg=itview->second.begin();
            itseg!=itview->second.end(); ++itseg)
          {
            int seg = itseg->first;

            // need to nudge hits which fell outside the track

            if((seg>-1) && (seg<NSeg))
              ThisViewdQdx[seg] = (itseg->second / fSegmentLength[seg]);
          }
        fdQdx.push_back(ThisViewdQdx);
      }

  }



  //----------------------------------------------------------------------
  //  Fill the dQdx vector for the track, based on a set of hits
  //    provided - optimized version

  void BezierTrack::CalculatedQdx(art::PtrVector<recob::Hit> const& Hits, std::vector<double> const& SValues)
  {
    fdQdx.clear();

    std::map<int, std::map<int, double> > hitmap;
    //        ^              ^       ^
    //       view            seg    charge

    for(size_t i=0; i!=Hits.size(); ++i)
      {

        (hitmap[Hits.at(i)->View()])[WhichSegment(SValues.at(i))] += Hits.at(i)->Integral();
      }

    int NSeg = NSegments();

    for(std::map<int,std::map<int,double> >::const_iterator itview=hitmap.begin();
        itview!=hitmap.end(); ++itview)
      {
        std::vector<double> ThisViewdQdx;
        ThisViewdQdx.resize(NSeg);
        for(std::map<int,double>::const_iterator itseg=itview->second.begin();
            itseg!=itview->second.end(); ++itseg)
          {
          
            int seg = itseg->first;

            // need to nudge hits which fell outside the track

            if((seg>-1) && (seg<NSeg))
              ThisViewdQdx[seg] = (itseg->second / fSegmentLength[seg]);
          }
        fdQdx.push_back(ThisViewdQdx);
      }
   

  }


  //----------------------------------------------------------------------
  //  Get dQdx for a particular S value
  //

  double BezierTrack::GetdQdx(double s, unsigned int view) const
  {
    view++;
    if((s<0.)||(s>1.))
      {
        throw cet::exception("Bezier dQdx: S out of range")
          <<"Bezier track S value of " << s <<" is not in the range 0<S<1"
          <<"\n";
      }
    if( /* (view<0)|| */ view>(fdQdx.size()-1))
      {
        throw cet::exception("Bezier dQdx: view out of range")
          <<"Bezier track view value of " << view <<" is not in the range "
          <<"of stored views, 0 < view < " << fdQdx.size()
          <<"\n";
      }
    int Segment = WhichSegment(s);
    //    if((Segment>=fdQdx[view].size())||(Segment<0))
    //      {
    //  throw cet::exception("Bezier dQdx: segment of range")
    //    <<"Bezier track segment value of " << Segment <<" is not in the range "
    //    <<"of stored segments, 0 < seg < " << fdQdx[view].size()
    //    <<"\n";
    //     }
    if( ( Segment < (int)fdQdx[view].size() ) && (Segment>1))
      return fdQdx[view][Segment];
    else
      {
        mf::LogVerbatim("BezierTrack")<<"GetdQdx : Bad s  " <<s<<", " <<Segment<<std::endl;
        return 0;
      }
  }




  //----------------------------------------------------------------------
  //  Get RMS dQdx for a particular view
  //

  double BezierTrack::GetViewdQdx(unsigned int view) const
  {
    view++;
    if(/* (view<0)|| */ view>fdQdx.size()-1)
      {
        throw cet::exception("view out of range")
          <<"Bezier track view value of " << view <<" is not in the range "
          <<"of stored views, 0 < view < " << fdQdx.size()
          <<"\n";
      }

    size_t NSeg = NSegments();
    double TotaldQdx = 0.;
    for(size_t i=0; i!=NSeg; ++i)
      {
        TotaldQdx+=fdQdx.at(view).at(i)*fSegmentLength[i];
      }
    return TotaldQdx / fTrackLength;
  }



  //----------------------------------------------------------------------
  //  Get total charge for a particular view
  //

  double BezierTrack::GetTotalCharge(unsigned int View) const
  {
    View++;
    return GetViewdQdx(View) * fTrackLength;
  }



  //----------------------------------------------------------------------
  //  Find which track segment a particular S lies in
  //
  int BezierTrack::WhichSegment(double s) const
  {
    int ReturnVal=-1;
    for(size_t i=0; i!=fCumulativeLength.size()-1; i++)
      {
        if( (fCumulativeLength.at(i)/fTrackLength <= s)
            &&(fCumulativeLength.at(i+1)/fTrackLength >= s))
          ReturnVal=i;
      }
    if( ( s<(fCumulativeLength.at(0)/fTrackLength)) && (s>-0.0001)) ReturnVal=0;
    return ReturnVal;
  }


  //----------------------------------------------------------------------
  // Return the number of track segments
  //

  int BezierTrack::NSegments() const
  {
    return fTraj.NumberTrajectoryPoints();
  }


  //----------------------------------------------------------------------

  std::vector<recob::SpacePoint> BezierTrack::GetSpacePointTrajectory(int N)
  {
    std::vector<recob::SpacePoint> spts(N);
    for(int i=0; i!=N; ++i)
      {
        double xyz[3];
        double ErrXYZ[3]={0.1,0.1,0.1};
        GetTrackPoint(float(i)/N,xyz);
        recob::SpacePoint TheSP(xyz, ErrXYZ, util::kBogusD, i);

        spts[i] = TheSP;

      }

    return spts;
  }


  //-----------------------------------------
  trkf::BezierTrack BezierTrack::Reverse()
  {
    std::vector<recob::Seed> SeedCol = GetSeedVector();
    std::reverse(SeedCol.begin(), SeedCol.end());
    for(size_t i=0; i!=SeedCol.size(); ++i)
      {
        SeedCol.at(i)=SeedCol.at(i).Reverse();
      }
    return BezierTrack(SeedCol);
  }

  //-----------------------------------------
  void BezierTrack::FillTrackVectors(std::vector<TVector3>& xyzVector,
                                     std::vector<TVector3>& dirVector,
                                     double const ds) const
  {
    const double s = ds / GetLength();
    const size_t n_traj_pts = (size_t)(GetLength()/ds);

    //+2: evenly space points, plus one for start and one for end
    xyzVector.resize(n_traj_pts+2);
    dirVector.resize(n_traj_pts+2);

    for(size_t i_traj=0; i_traj<=n_traj_pts; i_traj++){
      xyzVector[i_traj] = this->GetTrackPointV(i_traj*s);
      dirVector[i_traj] = this->GetTrackDirectionV(i_traj*s);
    }

    xyzVector[n_traj_pts+1] = this->GetTrackPointV(1.);
    dirVector[n_traj_pts+1] = this->GetTrackDirectionV(1.);


  }
  
} 




