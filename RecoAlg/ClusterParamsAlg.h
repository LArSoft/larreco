////////////////////////////////////////////////////////////////////////
// ClusterParamsAlg.h
//
// ClusterParamsAlg class
//
// Andrzej Szelc (andrzej.szelc@yale.edu)
//
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERPARAMSALG_H
#define CLUSTERPARAMSALG_H

#include "TMath.h"
#include <vector>



#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/HoughBaseAlg.h"


namespace recob { 
  class Hit;
  class Cluster; 
}

class TH2F;
class TF1;
class TH1F;
class TPrincipal;

namespace cluster {
   
  class ClusterParamsAlg {

  public:
    
  ClusterParamsAlg(fhicl::ParameterSet const& pset,std::string mother);
  
  void reconfigure(fhicl::ParameterSet const& pset); 
  
  int Find2DAxisRough(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist); /**Calculate 2D angle histograms, provided vertex is know */ 
    
  int Find2DAxisRoughHighCharge(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist);  
    //  void CalculateAxisParameters(unsigned nClust, std::vector < art::Ptr < recob::Hit> >  hitlist,double wstart,double tstart,double wend,double tend);
    
  int  Find2DStartPointsBasic(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end);

  int Find2DStartPointsBasic(double lineslope,double lineintercept,std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end);

  int Find2DStartPointsHighCharge(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end);
  
  
  int FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,double &wstn,double &tstn,double &wendn,double &tendn,double lineslope,double lineintercept);
  
  int FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,double &wstn,double &tstn,double &wendn,double &tendn); //overloaded interface, will calculate the lineslope and basic points using Find2DStartPointsBasic.
  
  void FindDirectionWeights(double lineslope,double wstn,double tstn,double wendn,double tendn, std::vector < art::Ptr < recob::Hit> > hitlist,double &HiBin,double &LowBin,double &invHiBin,double &invLowBin, double *altWeight=0); 
  
  void FindSideWeights(double lineslope,double lineintercept,double wstn,double tstn,int direction, std::vector < art::Ptr < recob::Hit> > hitlist,double &sideWeightCharge); 
  
  void RefineStartPointsHough(std::vector< art::Ptr < recob::Hit> > hitlist, 								double & wire_start,
				  double & time_start, 
				  double & wire_end,
				  double & time_end, 
				  int &direction);
 
  int DecideClusterDirection(std::vector < art::Ptr<recob::Hit> > hitlist,
				double lineslope,double wstn,double tstn,double wendn,double tendn);
 
  int FindHitCountDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end, double lineslope);
  
  int FindMultiHitDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end, double lineslope);
  
  bool isShower(double lineslope,double wstn,double tstn,double wendn,double tendn, std::vector < art::Ptr < recob::Hit> > hitlist);
  
  int MultiHitWires(std::vector < art::Ptr < recob::Hit> > hitlist);
  
//     
//     void RefineStartPoints(unsigned int nClust, std::vector< art::Ptr < recob::Hit> >  hitlist, double  wire_start,double time_start);
//     
//     double Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist);
//     
//     double Get2DAngleForHit( unsigned int wire, double time,std::vector < art::Ptr < recob::Hit> > hitlist);
int FindPrincipalDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end,double lineslope);
void GetPrincipal(std::vector< art::Ptr < recob::Hit> > hitlist, TPrincipal * pc);
int Find2DAxisPrincipal(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist);
int Find2DAxisPrincipalHighCharge(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist);
double Get2DAngleForHit( unsigned int swire,double stime,std::vector < art::Ptr < recob::Hit> > hitlist);
    
  private:
   
   const std::string extname; 
   std::string funcname; 
   HoughBaseAlg fHBAlg;  
   double fWiretoCm,fTimetoCm,fWireTimetoCmCm;
   double fWirePitch ;   // wire pitch in cm
   double fTimeTick;
   TF1 * linefittest_cm;
   TH2F *  tgxtest;
   TH1F *  hithistinv, * hitinv2, * hitreinv2;
   TH1F *  fh_omega_single;
   int fMinHitListSize;
   double fOutlierRadius;
   util::GeometryUtilities gser;
   art::ServiceHandle<util::DetectorProperties> detp;
   art::ServiceHandle<util::LArProperties> larp;
   art::ServiceHandle<geo::Geometry> geo;
   std::vector<double> fChargeCutoffThreshold;
   double fSelectBoxSizePar;
   double fSelectBoxSizePerp;
  // double fSelectBoxDiff;
   bool fForceRightGoing;
     //helper functions
  void FindExtremeIntercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low);  

  art::Ptr<recob::Hit> FindClosestHit(std::vector < art::Ptr < recob::Hit > > hitlist,
			     double wire_online,
			     double time_online);
   
   int GetPlaneAndTPC(art::Ptr<recob::Hit> a,
			unsigned int &p,
			unsigned int &cs,
			unsigned int &t,
			unsigned int &w);
   
   void SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double radlimit);
   
   void SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double linearlimit,   double ortlimit, double lineslopetest);
   
   std::vector< art::Ptr<recob::Hit> > CreateHighHitlist(double threshold,std::vector< art::Ptr<recob::Hit> > hitlist);
   
   void FindMinMaxWiresTimes(double &MinWire, double &MaxWire,double &MinTime,double &Maxtime,double &MeanCharge,std::vector< art::Ptr < recob::Hit> > hitlist);
   
   double fHitDensityCutoff;
   double fMultiHitWireCutoff;
   double fOffAxisHitsCutoff;
   double fPrincipalComponentCutoff;
   bool fShowerSelisORorAND;
   
   
   
  }; //class ClusterParamsAlg
  
} //namespace cluster








#endif