#include "RecoAlg/ClusterParamsAlg.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH1F.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "TPrincipal.h"


cluster::ClusterParamsAlg::ClusterParamsAlg(fhicl::ParameterSet const& pset,std::string mothermodule)
 : extname(mothermodule),fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
{

/**Get TFileService and define output Histograms*/
art::ServiceHandle<art::TFileService> tfs;
tgxtest=tfs->make<TH2F>(Form("%s charge hi distrib",extname.c_str()),"charge hi distribution per wires",
 			  10,0,10,10,0,10);
funcname=Form("%s linefit",extname.c_str());   //new name for the function - module specific
linefittest_cm=tfs->make<TF1>(funcname.c_str(),"pol1", 0,10); 

std::string histname=Form("%s omega_single",extname.c_str()); 
fh_omega_single= tfs->make<TH1F>(histname.c_str(),"Theta distribution Hit",720,-180., 180.) ;

//fChargeCutoffThreshold	=pset.get<   std::vector<double> > ("ChargeCutoffThreshold");



}




void cluster::ClusterParamsAlg::reconfigure(fhicl::ParameterSet const& pset) 
{
  fChargeCutoffThreshold	=pset.get<   std::vector<double> > ("ChargeCutoffThreshold");    
  fSelectBoxSizePar    		=pset.get<   double > ("SelectBoxSizePar");
  fSelectBoxSizePerp		=pset.get<   double > ("SelectBoxSizePerp");
  //fSelectBoxDiff		=pset.get<   double > ("SelectBoxDiff");
  fMinHitListSize		=pset.get<   double > ("MinHitListSize"); //40;
  fOutlierRadius		=pset.get<   double > ("fOutlierRadius"); //5;
  fForceRightGoing 		=pset.get<   bool > ("ForceRightGoing");
  fHitDensityCutoff		=pset.get<   double > ("HitDensityCutoff"); // Showers is high, tracks low, 0 don't use cut
  fMultiHitWireCutoff		=pset.get<   double > ("MultiHitWireCutoff"); // Showers is high, tracks low, 0 don't use cut
  fOffAxisHitsCutoff		=pset.get<   double > ("OffAxisHitsCutoff"); // Showers is high, tracks low, 0 don't use cut
  fPrincipalComponentCutoff	=pset.get<   double > ("PrincipalComponentCutoff"); // Tracks is high, showers low, 0 don't use cut
  fShowerSelisORorAND		=pset.get<   double > ("ShowerSelisORorAND"); // 0 = OR, 1 = AND
  
  fHBAlg.reconfigure(pset.get< fhicl::ParameterSet >("HoughBaseAlg"));
  
  fWirePitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  fTimeTick=detp->SamplingRate()/1000.; 
 
  //define conversion constants
  double fDriftVelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature());
  fWiretoCm=fWirePitch;
  fTimetoCm=fTimeTick*fDriftVelocity;
  fWireTimetoCmCm=(fTimeTick*fDriftVelocity)/fWirePitch;
  
 }




//---------------------------------------------------------------------------------------------
//Calculate Minimum and Maximum wires and times of the cluster. And Mean charge
//--------------------------------------------------------------------------------------------

void cluster::ClusterParamsAlg::FindMinMaxWiresTimes(double &MinWire, double &MaxWire,double &MinTime,double &MaxTime,double &MeanCharge,std::vector< art::Ptr < recob::Hit> > hitlist){

  unsigned int wire,tpc, cstat;
  unsigned int plane;

  double mcharge=0;
  
  if(hitlist.size()==0)
    return;
  
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    double time = (*hitIter)->PeakTime();  
    //time_C -= (presamplings+10.1);
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
    
    if(time>MaxTime)
	MaxTime=time;

    if(time<MinTime)
	MinTime=time;
    
    if(wire>MaxWire)
	MaxWire=wire;

    if(wire<MinWire)
	MinWire=wire;

    mcharge+=(*hitIter)->Charge();
    
  }
  
  MeanCharge=mcharge/hitlist.size();
  
  //fTotalCharge=mcharge;
  
  //fClusterPlane=plane;  // time coordinate of last point for each plane
      
 // //////std::cout << " +++ maximums " << MinWire << " " << MaxWire << " " <<MinTime << " " << MaxTime << std::endl;
 // //////std::cout << " +++ meancharge " << MeanCharge << std::endl;
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster, using hits above the average charge. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::Find2DAxisRoughHighCharge(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double kMinWire=99999,kMinTime=999999,kMaxWire=0,kMaxTime=0,kMeanCharge=0;
  FindMinMaxWiresTimes(kMinWire,kMaxWire,kMinTime,kMaxTime,kMeanCharge,hitlist);
  std::vector< art::Ptr<recob::Hit> > highhitlist=CreateHighHitlist(kMeanCharge,hitlist);
 // //////std::cout << "$$$$$$$ highhitlist size: " << highhitlist.size() << " " << kMeanCharge <<  std::endl;
    
  if(Find2DAxisRough(lineslope,lineintercept,goodness,highhitlist)==-1)
    return -1;
  
  return 0;
}

//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster, using hits above the average charge. Using Principal components
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::Find2DAxisPrincipalHighCharge(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double kMinWire=99999,kMinTime=999999,kMaxWire=0,kMaxTime=0,kMeanCharge=0;
  FindMinMaxWiresTimes(kMinWire,kMaxWire,kMinTime,kMaxTime,kMeanCharge,hitlist);
  std::vector< art::Ptr<recob::Hit> > highhitlist=CreateHighHitlist(kMeanCharge,hitlist);
//   //////std::cout << "$$$$$$$ highhitlist size: " << highhitlist.size() << " " << kMeanCharge <<  std::endl;
  if(Find2DAxisRough(lineslope,lineintercept,goodness,highhitlist)==-1)
    return -1;
  else
    return 0;
  
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::Find2DAxisPrincipal(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  //double time;
  //unsigned int wire,tpc, cstat;
  //unsigned int plane;

  
     //TVector2 HitMean=GetHitCenter(hitlist);
  if(hitlist.size()<5)
  {
   //////std::cout << " hitlist too small," <<hitlist.size() << " bailing " << std::endl;
    return -1;
  }
 
 TPrincipal  pc(2, "D");
 GetPrincipal( hitlist,&pc);  //get Principal component object from hitlist.

 TVector2 PrincipalDirection;
// double PrincipalEigenvalue;
 //PrincipalEigenvalue = (*pc.GetEigenValues())[0];
 PrincipalDirection.Set( (*pc.GetEigenVectors())[0][0],(*pc.GetEigenVectors())[0][1]);
 
 double a =0,c=0;
 if(PrincipalDirection.X()!=0) 
    a=-PrincipalDirection.Y()/PrincipalDirection.X();
 else
    a = 99999;
 
 //calculate intercept:
  for(unsigned int i=0; i!=hitlist.size(); i++)
   {
    c+=  hitlist[i]->PeakTime()*fTimetoCm -  hitlist[i]->WireID().Wire*fWiretoCm*a; 
   }

   c/=hitlist.size();
   
  lineslope=a/fWireTimetoCmCm;
  lineintercept=c/fTimetoCm;  
  goodness=(*pc.GetEigenValues())[1]*10;
  
 return 0; 
}







//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::Find2DAxisRough(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double time;
  unsigned int wire,tpc, cstat;
  unsigned int plane;

 
  if(hitlist.size()<5)
  {
   //////std::cout << " hitlist too small," <<hitlist.size() << " bailing " << std::endl;
    return -1;
  }
  
  
  
 // padding of the selected TGraph in wires and time. 
  int wirepad=20;
  int timepad=wirepad*fWiretoCm/fTimetoCm+0.5;
 
  double kMinWire=99999,kMinTime=999999,kMaxWire=0,kMaxTime=0,kMeanCharge=0;
  
  FindMinMaxWiresTimes(kMinWire,kMaxWire,kMinTime,kMaxTime,kMeanCharge,hitlist);
  
  
  //>\todo check that min/max wires are ok
  int nbinsx= (kMaxWire-kMinWire+2*wirepad)*fWiretoCm;  // nbins to have 
  int nbinsy= (kMaxTime-kMinTime+2.*(double)timepad)*fTimetoCm;  // nbins to have 
 
 
  tgxtest->Reset();
 // linefittest_cm->Reset();
  tgxtest->SetBins(nbinsx,((double)kMinWire-(double)wirepad)*fWiretoCm,			((double)kMaxWire+(double)wirepad)*fWiretoCm,nbinsy,
			  (kMinTime-timepad)*fTimetoCm,(kMaxTime+timepad)*fTimetoCm);
    
  linefittest_cm->SetRange(((double)kMinWire-(double)wirepad)*fWiretoCm,										((double)kMaxWire+(double)wirepad)*fWiretoCm); 
//     linefit2_cm->SetRange(((double)kMinWire-(double)wirepad)*fWiretoCm,										((double)kMaxWire+(double)wirepad)*fWiretoCm);
//   }
  
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != 				hitlist.end();  hitIter++){
    
    time =  (*hitIter)->PeakTime();  
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
  
    tgxtest->Fill((double)wire*fWiretoCm,
			      time*fTimetoCm,(*hitIter)->Charge());
  }
  
 
  tgxtest->Fit(funcname.c_str(),"+QMRNCFrob=0.95");


  
 double kRMS_wire=tgxtest->GetRMS(1);
 double kRMS_time=tgxtest->GetRMS(2);
 double kChisq;
 
  if(linefittest_cm->GetNDF())
    kChisq=linefittest_cm->GetChisquare()/linefittest_cm->GetNDF();
  else
    kChisq=1;
    
  double kCorrelation=tgxtest->GetCorrelationFactor();
  double kCovariance=tgxtest->GetCovariance();
  
  
  lineslope=linefittest_cm->GetParameter(1)/fWireTimetoCmCm;
  lineintercept=linefittest_cm->GetParameter(0)/fTimetoCm;

  
    if(kRMS_time && (kMaxTime-kMinTime) ){
  
      if(linefittest_cm->GetNDF() && tgxtest->GetCorrelationFactor() && tgxtest->GetCovariance())
	goodness=kRMS_time/kRMS_wire*((kMaxTime-kMinTime)*fTimetoCm)/((kMaxWire-kMinWire)*fWiretoCm)*kChisq/100000*1/kCorrelation*1/kCovariance*10;
      else
	goodness=-1;
    }

 
 

    
   return 0; 
}


//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Interface for version without the slope yet
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::Find2DStartPointsBasic(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
 double lineslope, lineintercept,goodness;

 if(Find2DAxisRough(lineslope,lineintercept,goodness,hitlist)==-1)
   return -1;
 Find2DStartPointsBasic(lineslope,lineintercept,hitlist,wire_start,time_start,wire_end,time_end);
 
 return 0;
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Interface for version without the slope yet
//--------------------------------------------------------------------------------------------

int   cluster::ClusterParamsAlg::Find2DStartPointsHighCharge(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
 double lineslope, lineintercept,goodness;

 if(Find2DAxisRoughHighCharge(lineslope,lineintercept,goodness,hitlist)==-1 )
    return -1;
 ////////std::cout << " lineslope: " << lineslope << " " << lineintercept << std::endl;
 Find2DStartPointsBasic(lineslope,lineintercept,hitlist,wire_start,time_start,wire_end,time_end);
 
 return 0;   
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------


int   cluster::ClusterParamsAlg::Find2DStartPointsBasic(double lineslope,double lineintercept,std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
  //  double time;
  unsigned int wire,plane,tpc,cstat;
  double a,c;
 
//  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != 				hitlist.end();  hitIter++){
//     
//     double time =  (*hitIter)->PeakTime();  
//     int wire= (*hitIter)->WireID().Wire;
//   
//    //////std::cout << "in StartPoints Basic time " << time << " " << wire << std::endl; 
//   }
//   
  
  
  GetPlaneAndTPC((*hitlist.begin()),plane,cstat,tpc,wire);
   
  //get paramters of the straight line fit. (and rescale them to cm/cm)
 
   a=lineslope*fWireTimetoCmCm;
   c=lineintercept*fTimetoCm;

  //get slope of lines orthogonal to those found crossing the shower.
  double aprim=0;
  if(a){
    aprim=-1./a*fWireTimetoCmCm*fWireTimetoCmCm;
  }
 else
   aprim = -999999.;
  

 
 ////////std::cout << "in starting points basic, line slope, intercept, a,c: " << lineslope << " " << lineintercept << " " << a << " " << c << std::endl;
 // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high,extreme_intercept_low;
  
  FindExtremeIntercepts(hitlist,aprim,extreme_intercept_high,extreme_intercept_low);
  double extreme_intercept_start=-999999;
  double extreme_intercept_end=999999;
  
 // int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.
  
  if(a>=0) {  // for the time being assuming forward going showers
   //  multiplier=1;
     extreme_intercept_end=extreme_intercept_high;
     extreme_intercept_start=extreme_intercept_low;
  }
  else if(a<0){
   // multiplier=-1;
    extreme_intercept_start=extreme_intercept_high;
    extreme_intercept_end=extreme_intercept_low;
  }
  
 // //////std::cout << " extreme intercepts: " << extreme_intercept_start << " " << extreme_intercept_end << std::endl;

  
  double wire_online_end,wire_online_begin,time_online_end,time_online_begin;
  
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_end,wire_online_end,time_online_end);
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_start,wire_online_begin,time_online_begin);
  
// //////std::cout << " ---------:::::::: wire_online begin point " << wire_online_begin << " " << time_online_begin << std::endl;
// //////std::cout << " ----------:::::::: wire_online end point " << wire_online_end << " " << time_online_end << std::endl;
  //calculate the first and last cluster points on the through line:
  
  art::Ptr<recob::Hit> startHit=FindClosestHit(hitlist, wire_online_begin,time_online_begin);
  art::Ptr<recob::Hit> endHit=FindClosestHit(hitlist, wire_online_end,time_online_end);
  
  GetPlaneAndTPC(startHit,plane,cstat,tpc,wire);
  wire_start=wire;
  time_start=startHit->PeakTime();
  
  GetPlaneAndTPC(endHit,plane,cstat,tpc,wire);
  wire_end=wire;
  time_end=endHit->PeakTime();
  
 return 0;
//  CalculateAxisParameters(nClust,hitlist,wire_start,time_start,wire_end,time_end);
}



//---------------------------------------------------------------------------------------------
// Find the trunk of the cluster and use it to determine new rought start/end points.
// uses, rough strarting points found earlier and the slope of the line found earlier
// returns preliminary shower direction //interface that executes all the stuff before it
//--------------------------------------------------------------------------------------------
int cluster::ClusterParamsAlg::FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,double &wstn,double &tstn,double &wendn,double &tendn)
{
 double lineslope, lineintercept,goodness;
  
 
 if(Find2DAxisRough(lineslope,lineintercept,goodness,hitlist)==-1 )
  return -99;
 Find2DStartPointsBasic(lineslope,lineintercept, hitlist,wstn,tstn,wendn,tendn); 
 return FindTrunk(hitlist,wstn,tstn,wendn,tendn,lineslope,lineintercept); 
  
}



//---------------------------------------------------------------------------------------------
// Find the trunk of the cluster and use it to determine new rought start/end points.
// uses, rough strarting points found earlier and the slope of the line found earlier
// returns preliminary shower direction //interface that needs earlier input
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,
					 double &wstn,
					 double &tstn,
					 double &wendn,
					 double &tendn,
					 double lineslope,
					 double lineintercept)
  {
  int fDirection=0;
  unsigned int currplane=999;
  double fTotalCharge;
  
 
  
  if(hitlist.size()==0)
    return 0;

  double at=lineslope*fWireTimetoCmCm;
  double ct=lineintercept*fTimetoCm;

  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
  
 
  // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high,extreme_intercept_low;
    
  FindExtremeIntercepts(hitlist,aprimt,extreme_intercept_high,extreme_intercept_low);
    
  double extreme_intercept_start=-999999;
  double extreme_intercept_end=999999;
  
 
  if(at>=0) {  // for the time being assuming forward going showers
     extreme_intercept_end=extreme_intercept_high;
     extreme_intercept_start=extreme_intercept_low;
  }
  else if(at<0){
    extreme_intercept_start=extreme_intercept_high;
    extreme_intercept_end=extreme_intercept_low;
  }
  
//projections of the start and end points onto the rough cluster axis
  double wire_online_end,wire_online_begin,time_online_end,time_online_begin;

  gser.GetPointOnLineWSlopes(at,ct,extreme_intercept_end,wire_online_end,time_online_end);
  gser.GetPointOnLineWSlopes(at,ct,extreme_intercept_start,wire_online_begin,time_online_begin);
  

  double projectedlength=TMath::Sqrt( ((double)wire_online_end-(double)wire_online_begin)*((double)wire_online_end-(double)wire_online_begin)*fWiretoCm*fWiretoCm +  ((double)wire_online_end-(double)wire_online_begin)*((double)wire_online_end-(double)wire_online_begin)*at*at/fWireTimetoCmCm/fWireTimetoCmCm*fTimetoCm*fTimetoCm  );;

 //number of bins for a test histogram
  //int nbins = (hitlist.size() > 1000) ? 100 : hitlist.size()/10;

  hithistinv=new TH1F(Form("hithistinv_ev"),Form("hithistinv_pl"),(int)projectedlength,0.,projectedlength);  
  hithistinv->Clear();
  hithistinv->SetName(Form("hithistinv_pl"));
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
    if(currplane==999) currplane=plane;
    
    double wire_on_line,time_on_line;
    
    gser.GetPointOnLine(at/fWireTimetoCmCm,wire_online_begin,time_online_begin,wire,time,wire_on_line,time_on_line);
    
  //get linear distance from beginning and orthogonal distance from cluster axis
    double linedist=gser.Get2DDistance(wire_on_line,time_on_line,wire_online_begin,time_online_begin);
    double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
    
   ////////////////////////////////////////////////////////////////////// 
   //calculate the weight along the axis, to guess the Shower direction  
   /////////////////////////////////////////////////////////////////////// 
    double weight= (ortdist<1.) ? 10*theHit->Charge() : 5*theHit->Charge()/ortdist;
    hithistinv->Fill(linedist,weight);
    
    fTotalCharge+=(*hitIter)->Charge();  //calculate total charge while we're at it.
        
  }
 
  //Now we have a histogram with the projection of the charge along the axis.
  
    double start_off=0., end_off=((double)wire_online_end-(double)wire_online_begin); 
    
    double fMaxweight=hithistinv->GetMaximum();  
  
   // calculate integrals along the axis: total, from maximum to start and from maximum to end 
   // double fullinteg=hithistinv->Integral();
    double startinteg=hithistinv->Integral(0,hithistinv->GetMaximumBin());
    double endinteg=hithistinv->Integral(hithistinv->GetMaximumBin(),hithistinv->GetNbinsX());
  
    
    if(fMaxweight>fChargeCutoffThreshold[currplane]) //check whether this operation makes sense
    {
      for(int x= hithistinv->GetMaximumBin();x>0;x--)
	{
	  if( hithistinv->GetBinContent(x)<fChargeCutoffThreshold[currplane] || x==1)
	  {
	    start_off=x;
	    if(  x>1 && hithistinv->Integral(0,x)/startinteg <0.01  ){
	      break;
	  }
	}
      }
 
      for(int x=hithistinv->GetMaximumBin();x<hithistinv->GetNbinsX();x++)
	{
	if( hithistinv->GetBinContent(x)<fChargeCutoffThreshold[currplane] || x==hithistinv->GetNbinsX()-1)
	  {
	  end_off=x;
	  if(  x<hithistinv->GetNbinsX()-1  && hithistinv->Integral(x,hithistinv->GetNbinsX())/endinteg <0.01  ){
	  break;
	  }
	}
      }
    }	
    
    
    
 //now create new rough points on axis based on the findings of the histogram above:
  double w_on_end,w_on_begin,t_on_end,t_on_begin;

  //calculate intervals in wire and time
  double dws= start_off*fabs(((double)wire_online_begin-(double)wire_online_end))/projectedlength;
  double dwe= end_off*fabs(((double)wire_online_begin-(double)wire_online_end))/projectedlength;

  double dts=at*dws/fWireTimetoCmCm;
  double dte=at*dwe/fWireTimetoCmCm;

  w_on_begin=wire_online_begin+dws;
  w_on_end=wire_online_begin+dwe;
  t_on_begin=time_online_begin+dts;
  t_on_end=time_online_begin+dte;

  ////////std::cout << "inside findTrunk, start points: " << w_on_begin << " "<<t_on_begin << " | "<< w_on_end << "  " << t_on_end << std::endl;
  /////// create new, two bin histogram, to test where the front and back are:
   double lowbin=0,highbin=0;
   double relowbin=0,rehighbin=0;
   double side_weight_start_charge=0,side_weight_end_charge=0;
   
  FindDirectionWeights(lineslope,w_on_begin, t_on_begin,w_on_end,t_on_end, hitlist,highbin,lowbin,rehighbin,relowbin); 
  
  FindSideWeights(lineslope,lineintercept,w_on_begin,t_on_begin,1,hitlist,side_weight_start_charge);
  FindSideWeights(lineslope,lineintercept,w_on_end,t_on_end,-1,hitlist,side_weight_end_charge);

  

  
  ///////////////////////////////////////////////////////
  // At this point the start still means most leftwise, end rightwise.
  //////////////////////////////////////////////////////
  
 
//   //////std::cout << "==== low and highbin " << lowbin << " " << highbin << std::endl;
//   //////std::cout << "==== reverse low and highbin " << relowbin << " " << rehighbin << std::endl;
//   //////std::cout << "==== reverse low and highbin hist" << hitreinv2->GetBinContent(1)<< " " << hitreinv2->GetBinContent(2) << std::endl;
//   
//   //////std::cout << " histo coords: LowE " <<  hitreinv2->GetBinLowEdge(1) << " " <<  hitreinv2->GetBinLowEdge(2) << " " << hitreinv2->GetBinWidth(1) << " " <<  hitreinv2->GetBinWidth(2) << std::endl;
  
  /////////Refine start points to real hits:
  double wstartnorm=0,tstartnorm=0,wendnorm=0,tendnorm=0;
  double wstartwght=0,tstartwght=0,wendwght=0,tendwght=0;
  double mindiststart=999999,mindistend=999999,minwghtdstart=999999,minwghtdend=999999;
  double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
  double inv_intercept_end=t_on_end*fTimetoCm-aprimt*(double)w_on_end*fWiretoCm; 
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
    double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
    // we're basically only interested in the sign of these.
    double inter_distance_start = (inv_intercept - inv_intercept_start)*at;   // the a makes the sign work for negative.
    double inter_distance_end = (inv_intercept_end - inv_intercept)*at;   // the a makes the sign work for negative.
  
    if(inter_distance_start<0 || inter_distance_end<0)
      continue;
    
    
    double linedistst=gser.Get2DDistance(wire,time,w_on_begin,t_on_begin);
    double linedistend=gser.Get2DDistance(wire,time,w_on_end,t_on_end);
    double intercept=time*fTimetoCm-at*(double)wire*fWiretoCm;
    double locweight=(intercept-ct);
   
    //start points
    
    if(linedistst<mindiststart)
      {
      wstartnorm=wire; tstartnorm=time;
      mindiststart=linedistst;
      }
      // ony use the weight if it is big enough - otherwise you risk going the wrong way.
    if(linedistst<minwghtdstart && fabs(side_weight_start_charge/fTotalCharge) > 0.003 && side_weight_start_charge*locweight >=0 )
      {
      wstartwght=wire; tstartwght=time;
      minwghtdstart=linedistst;
      }
    
     //end points 
    
    if(linedistend<mindistend)
      {
      wendnorm=wire; tendnorm=time;
      mindistend=linedistend;
      }
     // ony use the weight if it is big enough - otherwise you risk going the wrong way.
   if(linedistend<minwghtdend &&  fabs(side_weight_start_charge/fTotalCharge) > 0.003 && side_weight_end_charge*locweight >=0 )
      {
      wendwght=wire; tendwght=time;
      minwghtdend=linedistend;
      }
 } // end refining of starting points
  
  if(wstartwght==0)
    {
    wstartwght=wstartnorm;
    tstartwght=tstartnorm;
    }
  if(wendwght==0)
    {
    wendwght=wendnorm;
    tendwght=tendnorm;
    }
  
 // //////std::cout << "**** starting points. standard: " << wstartnorm << "," << tstartnorm << " weighted: " << wstartwght <<","<<tstartwght << std::endl;
 // //////std::cout << "**** ending points. standard: " << wendnorm << "," << tendnorm << " weighted: " << wendwght <<","<<tendwght << std::endl;
  
  if(highbin>lowbin)  // cluster is leftward going
  {
    wstn=wendwght; 
    tstn=tendwght;
    wendn=wstartwght;
    tendn=tstartwght;
  
    fDirection=-1;
  }
  else  // cluster is rightward going 
  {
   wstn=wstartwght; 
   tstn=tstartwght;
   wendn=wendwght;
   tendn=tendwght; 
   fDirection=1;
  }
    
//  //////std::cout << "**** Final starting points. standard: " << wstn << "," << tstn << std::endl;
    

return fDirection;

}





void cluster::ClusterParamsAlg::FindDirectionWeights(double lineslope,double w_on_begin,double t_on_begin,double w_on_end,double t_on_end, std::vector < art::Ptr < recob::Hit> > hitlist,double &HiBin,double &LowBin,double &invHiBin,double &invLowBin, double *altWeight)
{
  
  if(altWeight!=0) *altWeight=0.;
  
  double at=lineslope*fWireTimetoCmCm;
  //double ct=lineintercept*fTimetoCm;
  double smallprojlength= TMath::Sqrt( ((double)w_on_end-(double)w_on_begin)*((double)w_on_end-(double)w_on_begin)*fWiretoCm*fWiretoCm +  ((double)w_on_end-(double)w_on_begin)*((double)w_on_end-(double)w_on_begin)*at*at/fWireTimetoCmCm/fWireTimetoCmCm*fTimetoCm*fTimetoCm  );
  
  ////////std::cout << " smallprojlength: "<< smallprojlength << std::endl;
  
  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
  
   //histogram for weight, i.e. hits along the axis:
  hitinv2=new TH1F(Form("hitinv2_ev"),Form("hitinv2_ev"),2,0.,smallprojlength);  
  double lowbin=0,highbin=0;

  //histogram for inverse weight:
  hitreinv2=new TH1F(Form("rehitinv2_ev"),Form("rehitinv2_ev"),2,0.,smallprojlength);  
  double relowbin=0,rehighbin=0;

  hitinv2->Clear();
  hitinv2->SetName(Form("hitinv2_ev"));
  hitreinv2->Clear();
  hitreinv2->SetName(Form("rehitinv2_"));
  
  ///////////////////////////////////////////////////////
  // At this point the start still means most leftwise, end rightwise.
  //////////////////////////////////////////////////////
  
  
   //// these are the intercepts of perpendicular lines going through the start and end point.
   //// I will use the perpendicular intercept as a measurement of distance for simplicity.
    
   double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
   double inv_intercept_end=t_on_end*fTimetoCm-aprimt*(double)w_on_end*fWiretoCm;
    
  // //////std::cout << "intercepts in weights: " << inv_intercept_start << " " << inv_intercept_end << " " << hitlist.size()  <<std::endl;
 
   ///////// Calculate weights
    for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = theHit->PeakTime() ;  
      unsigned int wire,cstat, tpc,plane;
      GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
      double wire_on_line,time_on_line;
    
      gser.GetPointOnLine(at/fWireTimetoCmCm,w_on_begin,t_on_begin,wire,time,wire_on_line,time_on_line);
      double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
      // we're basically only interested in the sign of these.
      double inter_distance_start = (inv_intercept - inv_intercept_start)*at;   // the a makes the sign work for negative.
      double inter_distance_end = (inv_intercept_end - inv_intercept)*at;   // the a makes the sign work for negative.
  
  
      // startEnd checking
  
      double linedist=gser.Get2DDistance(wire_on_line,time_on_line,w_on_begin,t_on_begin);
      double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
      double weight= (ortdist<1.) ? 10*theHit->Charge() : 5*theHit->Charge()/ortdist;
      double revweight= (ortdist<1.) ? 0.1 : ortdist;

      if(altWeight!=0 && ortdist>0.5) *altWeight+=ortdist; 
      
      hitreinv2->Fill(linedist,revweight);
      hitinv2->Fill(linedist,weight);
 
      if(inter_distance_start>0. && linedist<= smallprojlength*0.5 )
      {
      lowbin+=weight;
      relowbin+=revweight;
      }
      else if (linedist>smallprojlength*0.5 && inter_distance_end > 0. ) 
      { highbin+=weight;
	rehighbin+=revweight;
      }
     
 
       
  }  // end weight calculating loop
  
  
  HiBin=highbin;
  LowBin=lowbin;
  invHiBin=rehighbin;
  invLowBin=relowbin;
  
  
}

// ---------------------------------------------------------------------------------------
// function to decide whether cluster is shower like or track like
// ---------------------------------------------------------------------------------------

bool cluster::ClusterParamsAlg::isShower(double lineslope,double wstn,double tstn,double wendn,double tendn, std::vector < art::Ptr < recob::Hit> > hitlist)
{
  
//   fHitDensityCutoff	 // Showers is high, tracks low, 0 don't use cut
//   fMultiHitWireCutoff	 // Showers is high, tracks low, 0 don't use cut
//   fOffAxisHitsCutoff	 // Showers is high, tracks low, 0 don't use cut
//   fPrincipalComponentCutoff // Tracks is high, showers low, 0 don't use cut
//   fShowerSelisORorAND	  //// 0 = OR, 1 = AND
  
   double HiBin,LowBin,invHiBin,invLowBin,altWeight=0.;
   double PrincipalEigenvalue=1.,ModPV=-900;
   double multihit=0; 
   double MODHitDensity=0;
   
   // no hits = no shower
   if(hitlist.size()==0)
      return false;
   
   /////////////////////////////////////////////////////////////////////////
   // calculate cluster 2D length and 2D angle 
   double length=TMath::Sqrt( (wstn-wendn)*(wstn-wendn)*fWiretoCm +(tstn-tendn)*(tstn-tendn)*fTimetoCm );
   double xangle=Get2DAngleForHit( wstn,tstn, hitlist);
    if(xangle>90) xangle-=180;
    if(xangle<-90) xangle+=180;  
   
   
 //////////////////////////////////////////////////////////  
   // do we apply the off AxisHitsCut? 
   if(fOffAxisHitsCutoff)
   {
      FindDirectionWeights(lineslope,wstn,tstn,wendn,tendn,hitlist,HiBin,LowBin,invHiBin,invLowBin,&altWeight); 
      altWeight/=length;
   }
   else if(fShowerSelisORorAND==false) //using or
       altWeight=0;   // use small value so that the OR doesn't select gibberish.
   else if(fShowerSelisORorAND==true)  //using AND
       altWeight=1;   // use large value, so that AND doesn't balk.
  //////////////////////////////////////////////////////// 
    // do we apply the Principal Component cut? 
   if( fPrincipalComponentCutoff)
    {
       TPrincipal pc(2,"D");
       GetPrincipal(hitlist,&pc);
       PrincipalEigenvalue = (*pc.GetEigenValues())[0];
       ModPV=TMath::Log(1-PrincipalEigenvalue);
    }
    else
    {
    if(fShowerSelisORorAND==false) //using or
       { 
       PrincipalEigenvalue=1;   // use small value so that the OR doesn't select gibberish.
       ModPV=-999;
       }
    }

 
    //////////////////////////////////////////////////////// 
    // do we apply the MultiHit cut? 
   if( fMultiHitWireCutoff)
    {
	multihit=MultiHitWires(hitlist)/length;
    }
    else
    {
    if(fShowerSelisORorAND==false) //using or
       multihit=0;   // use small value so that the OR doesn't select gibberish.
    }
   
    //////////////////////////////////////////////////////// 
    // do we apply the HitDensity cut? 
    if( fHitDensityCutoff)
    {
	double HitDensity=hitlist.size()/length;
// 	MODHitDensity=HitDensity+std::abs(xangle)*0.0299;
	MODHitDensity=HitDensity-(1.82 * cosh(3.14/180*abs(xangle)-1.24)-1.56);
    }
    else
    {
    if(fShowerSelisORorAND==false) //using or
       MODHitDensity=0;   // use small value so that the OR doesn't select gibberish.
    }
   
   std::string ororand= (fShowerSelisORorAND) ? " && " : " || ";
         
   //////std::cout << "##########3 Cutting on: ModHitDensit:" << MODHitDensity << " > " << fHitDensityCutoff << ororand 
//   << " MultiHitWires " << multihit << " > " << fMultiHitWireCutoff << ororand 
//   << " OffAxisHits " << altWeight << " > " << fOffAxisHitsCutoff << ororand 
//   << " PrincipalComponent " << ModPV << " < " << fPrincipalComponentCutoff
//   << std::endl;
   
      //(fHitDensity[0]>1.9  || fOffAxisNorm[0] > 1.) &&(fHitDensity[1]>1.9  || fOffAxisNorm[1] > 1.)
    if( !fShowerSelisORorAND && 
      (MODHitDensity>fHitDensityCutoff  || 
      multihit > fMultiHitWireCutoff || 
      PrincipalEigenvalue < fPrincipalComponentCutoff || 
      altWeight > fOffAxisHitsCutoff ) )
	return true;
    else if ( fShowerSelisORorAND && 
      (MODHitDensity>fHitDensityCutoff  && 
      multihit > fMultiHitWireCutoff && 
      ModPV < fPrincipalComponentCutoff && ModPV >= -4 &&
      altWeight > fOffAxisHitsCutoff ) )
	return true;
    else
      return false;
  
}


int cluster::ClusterParamsAlg::MultiHitWires(std::vector < art::Ptr < recob::Hit> > hitlist)
{
  int multihit=0;
    for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      //art::Ptr<recob::Hit> theHit = (*hitIter);
      unsigned int wire=(*hitIter)->WireID().Wire;   //current hit wire/time
      for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitItersecond = hitIter; hitItersecond != hitlist.end();  hitItersecond++){
	unsigned int wiresecond=(*hitItersecond)->WireID().Wire;   //current hit wire/time
	if (wiresecond==wire)
	  multihit++;             // not breaking - want multicounting to blow up the showers for now.
      }
  
    } 
  
return multihit;   
  
}


double cluster::ClusterParamsAlg::Get2DAngleForHit( unsigned int swire,double stime,std::vector < art::Ptr < recob::Hit> > hitlist) {
  
  fh_omega_single->Reset();
  util::GeometryUtilities gser;
  unsigned int wire;
  // this should changed on the loop on the cluster of the shower
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = theHit->PeakTime();  
      wire=theHit->WireID().Wire; 
      double omx=gser.Get2Dangle((double)wire,(double)swire,time,stime);
      fh_omega_single->Fill(180*omx/TMath::Pi(), theHit->Charge());
     }
    
  double omega = fh_omega_single->GetBinCenter(fh_omega_single->GetMaximumBin());// Mean value of the fit
   
  return omega; // in degrees.
}




void cluster::ClusterParamsAlg::FindSideWeights(double lineslope,double lineintercept,double w_on_begin,double t_on_begin, int direction,std::vector < art::Ptr < recob::Hit> > hitlist,double &side_weight_start_charge)
{
 
 
  double at=lineslope*fWireTimetoCmCm;
  double ct=lineintercept*fTimetoCm;
   
  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
   

   //// these are the intercepts of perpendicular lines going through the start and end point.
   //// I will use the perpendicular intercept as a measurement of distance for simplicity.
    
   double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
    
   int startctr=0;

   ///////// Calculate weights
    for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = theHit->PeakTime() ;  
      unsigned int wire,cstat, tpc,plane;
      GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
      double wire_on_line,time_on_line;
    
      gser.GetPointOnLine(at/fWireTimetoCmCm,w_on_begin,t_on_begin,wire,time,wire_on_line,time_on_line);
      double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
      // we're basically only interested in the sign of these.
      double inter_distance_start = (inv_intercept - inv_intercept_start)*at*direction;   // the a makes the sign work for negative., the direction, differentiates between points at start and end.
      
      double linedist=gser.Get2DDistance(wire_on_line,time_on_line,w_on_begin,t_on_begin);
    
 
           
     //now calculate the side weight, i.e. on which side of the axis more charge is found (on the start side)
      if( linedist < 10. && inter_distance_start > 0.  )
	{
	// this is the intercept the line would have if it went through the point
	double intercept=time*fTimetoCm-at*(double)wire*fWiretoCm;
    
	side_weight_start_charge+=(intercept-ct)*theHit->Charge();
	startctr++;
     	}
    
     
     
       
  }  // end weight calculating loop
 //////////////////////////////////////////////////
   if(startctr>0)
    {
      side_weight_start_charge/=startctr;
    }
  
        
}





//---------------------------------------------------------------------------------------------
//Create new Hitlist consisting of hits higher than a Threshold and not further away than 
// radius: 
// and the list not smaller then 40: hits:
//--------------------------------------------------------------------------------------------


std::vector< art::Ptr<recob::Hit> > cluster::ClusterParamsAlg::CreateHighHitlist(double threshold,std::vector< art::Ptr<recob::Hit> > hitlist)
{
  
  std::vector< art::Ptr<recob::Hit> > hitlist_high;
    
  //first select a subset of points that is close to the selected start point:
  std::vector < art::Ptr<recob::Hit> > hitlistlocal;
  
  for(unsigned int ix = 0; ix<  hitlist.size();  ix++){
    	art::Ptr<recob::Hit> theHit = hitlist[ix];
	if(theHit->Charge()>threshold)
	  hitlist_high.push_back(theHit);
    }
  
  //subtracting outliers
    
  // this is the erasing loop
 for(unsigned int ix = 0; ix<  hitlist_high.size();  ix++){
   if(hitlist_high.size()<(unsigned int)fMinHitListSize)
	  break;
   art::Ptr<recob::Hit> theHit = hitlist_high[ix];
   double time = theHit->PeakTime() ;  
   unsigned int plane,cstat,tpc,wire;
   GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
   //ok. now I have a sub selection of hits that are close to the perceived start of the shower.
   SelectLocalHitlist(hitlist_high, hitlistlocal, (double)wire,time,fOutlierRadius);
   
  
   if(hitlistlocal.size()<10 && hitlist_high.size() > ix && hitlist_high.size()>0 ){
      hitlist_high.erase(hitlist_high.begin()+ix);
      ix--;
    //  //////std::cout << " erasing hit @ w,t" << wire << " "<< time << " (start) "  << std::endl;
    }
   hitlistlocal.clear();
	

 } 
  
return hitlist_high;  
}



int cluster::ClusterParamsAlg::FindPrincipalDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end, double lineslope)
{
  
  double linearlimit=fSelectBoxSizePar; 
  double ortlimit=fSelectBoxSizePerp;
  int direction_local=1;
  
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_start;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_end;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  
  ///////////////////////////// find trunk of shower:
  //double wstn=0.,tstn=0.,wendn=0.,tendn=0.;
  
  //FindTrunk(nClust,hitlist,wstn,tstn,wendn,tendn);
  
  /////
    
  //first select a subset of points that is close to the selected start point:
  
 SelectLocalHitlist(hitlist, hitlistlocal_start, wire_start,time_start,linearlimit,ortlimit,lineslope);
 SelectLocalHitlist(hitlist, hitlistlocal_end, wire_end,time_end,linearlimit,ortlimit,lineslope);
  
 TPrincipal  pc_start(2, "D");
 TPrincipal  pc_end(2, "D");
 GetPrincipal( hitlistlocal_start,&pc_start);  //get Principal component object from hitlist.
 GetPrincipal( hitlistlocal_end,&pc_end);  //get Principal component object from hitlist.
 
 double PrincipalEigenvalue_start = (*pc_start.GetEigenValues())[0];
 double PrincipalEigenvalue_end = (*pc_end.GetEigenValues())[0];
 
 //////std::cout << " EigenValue_start: " << PrincipalEigenvalue_start << " "<< PrincipalEigenvalue_end << std::endl;
 
 if(PrincipalEigenvalue_start/hitlistlocal_start.size() > PrincipalEigenvalue_end/hitlistlocal_end.size() )
   direction_local=1;
 else 
   direction_local=-1;
 
 if(wire_start>wire_end)
 { 
   direction_local*=-1;     //make the principal calculated direction absolute.
 }
 
  
  return direction_local;
}


int cluster::ClusterParamsAlg::FindMultiHitDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end, double lineslope)
{
  
  double linearlimit=fSelectBoxSizePar; 
  double ortlimit=fSelectBoxSizePerp;
  int direction_local=1;
  
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_start;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_end;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  

 SelectLocalHitlist(hitlist, hitlistlocal_start, wire_start,time_start,linearlimit,ortlimit,lineslope);
 SelectLocalHitlist(hitlist, hitlistlocal_end, wire_end,time_end,linearlimit,ortlimit,lineslope);
  
 int multihit_start=MultiHitWires(hitlistlocal_start);
 int multihit_end=MultiHitWires(hitlistlocal_start);
 
 //////std::cout << "-------- multihits; start: " << multihit_start << " end: " << multihit_end << std::endl;
 
 
 if( multihit_start <  multihit_end )
   direction_local=1;
 else 
   direction_local=-1;

 if(wire_start>wire_end)
 { 
   direction_local*=-1;     //make the principal calculated direction absolute.
 }
  //////std::cout << " -- multihit direction: " << direction_local << std::endl;
  
  return direction_local;
}


int cluster::ClusterParamsAlg::FindHitCountDirection(std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double  time_start, double  wire_end,double  time_end, double lineslope)
{
  
  double linearlimit=fSelectBoxSizePar; 
  double ortlimit=fSelectBoxSizePerp;
  int direction_local=1;
  
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_start;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_end;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  

 SelectLocalHitlist(hitlist, hitlistlocal_start, wire_start,time_start,linearlimit,ortlimit,lineslope);
 SelectLocalHitlist(hitlist, hitlistlocal_end, wire_end,time_end,linearlimit,ortlimit,lineslope);
  
int hitcount_start=hitlistlocal_start.size();
int hitcount_end=hitlistlocal_end.size(); 
 
 //////std::cout << "+++++++++ HitCount; start: " << hitcount_start << " end: " << hitcount_end << std::endl;
 
 
 if( hitcount_start <  hitcount_end )
   direction_local=1;
 else 
   direction_local=-1;

 if(wire_start>wire_end)
 { 
   direction_local*=-1;     //make the principal calculated direction absolute.
 }
 
   //////std::cout << " -- HitCount direction: " << direction_local << std::endl;
  return direction_local;
}



/////////////////////////////////////////////////////////////////////////
// Helper function to return Principal Components from a Hitlist
// cluster.
//////////////////////////////////////////////////////////////////////////

void cluster::ClusterParamsAlg::GetPrincipal(std::vector< art::Ptr < recob::Hit> > hitlist, TPrincipal * pc)
{
   
  
  for(unsigned int i=0; i!=hitlist.size(); i++)
    {
     double ThisStep[2];
//      ThisStep[0] = hitlist[i]->WireID().Wire*fWiretoCm - HitMean.X();
//      ThisStep[1] = hitlist[i]->PeakTime()*fTimetoCm - HitMean.Y();
     ThisStep[0] = hitlist[i]->WireID().Wire*fWiretoCm ;
     ThisStep[1] = hitlist[i]->PeakTime()*fTimetoCm ;
     pc->AddRow(ThisStep);
    }

 pc->MakePrincipals();

  
 return;
 }



/////////////////////////////////////////////////////////////////////////
// refine the provided start points by looking for Hough Lines at the start and end of 
// cluster.
//////////////////////////////////////////////////////////////////////////

void cluster::ClusterParamsAlg::RefineStartPointsHough(std::vector< art::Ptr < recob::Hit> > hitlist, double & wire_start,double & time_start, double & wire_end,double & time_end, int &direction)
{
  //parameters to be set later?
  double linearlimit=fSelectBoxSizePar; 
  double ortlimit=fSelectBoxSizePerp;
  int direction_local=1;
  
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_start;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_end;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  
  ///////////////////////////// find trunk of shower:
  //double wstn=0.,tstn=0.,wendn=0.,tendn=0.;
  
  //FindTrunk(nClust,hitlist,wstn,tstn,wendn,tendn);
  
  /////
  double lineslope=0; //either add as parameter or pick up from Find Rough
  
  //first select a subset of points that is close to the selected start point:
  
 SelectLocalHitlist(hitlist, hitlistlocal_start, wire_start,time_start,linearlimit,ortlimit,lineslope);
 SelectLocalHitlist(hitlist, hitlistlocal_end, wire_end,time_end,linearlimit,ortlimit,lineslope);
 ////////std::cout << "Hough hitlist_local start size " <<  hitlistlocal_start.size() << std::endl;
 ////////std::cout << "Hough hitlist_local end size " <<  hitlistlocal_end.size() << std::endl;
 std::vector < art::PtrVector<recob::Hit> > houghlines; 
 std::vector < art::PtrVector<recob::Hit> > houghlinesend; 
 
 //save sizes of hitlists found
 unsigned int lochitlistsize=hitlistlocal_start.size();
 unsigned int lochitlistendsize=hitlistlocal_end.size();
   
   
//  for(unsigned int ix=0;ix<hitlistlocal_start.size();ix++) {
//    unsigned int w;
//    //unsigned int c = (*hitlistlocal_start[ix]).Wire()->RawDigit()->Channel(); 
//    //geo->ChannelToWire(c,cs,t,p,w);
//    w  = (*hitlistlocal_start[ix]).WireID().Wire;
//   //  //////std::cout << " local hits for: " << " plane: " << p << " w,t: "  << w<< " "<< (*hitlistlocal_start[ix]).PeakTime() << std::endl;
//    }
  // //////std::cout<<std::endl<<std::endl;    
    
  unsigned int iplane,cstat,tpc,wire;
  GetPlaneAndTPC((*hitlist.begin()),iplane,cstat,tpc,wire);
// 
  size_t numclus =0; 
    if(hitlistlocal_start.size()>5) numclus = fHBAlg.FastTransform(hitlistlocal_start, houghlines);
    
  size_t numclusend =0;
  if(hitlistlocal_end.size()>5) numclusend=fHBAlg.FastTransform(hitlistlocal_end, houghlinesend);
 
 
 // //////std::cout << "found " << numclus << "lines with HoughBaseAlg ";
 // if(numclus>0)
 //  //////std::cout << houghlines[0].size();
 // //////std::cout << std::endl;
//  //////std::cout << "found " << numclusend << "end lines with HoughBaseAlg ";
 // if(numclusend>0)
 //  //////std::cout << houghlinesend[0].size();
//  //////std::cout << std::endl; 

  //int nhoughlines=numclus;
  //int nhoughlinesend=numclusend;
 
  int starthoughclus=-1,endhoughclus=-1;  //index of largest houghline
 // int numclusfunc;
 // std::vector < art::PtrVector<recob::Hit> > *houghutil; 
 // numclusfunc=numclus;
 // houghutil=&houghlines;
 
 // assume that direction inserted is correct (wire_start is start), unless found otherwise
 //if(numclusend>=1 && ( (fReLowBin > fReHighBin && fDirection==1 ) || (fReLowBin < fReHighBin && fDirection==-1 ) ) )
  if(numclusend>=1 )
  {
   unsigned int stsize=0; 
   unsigned int endsize=0;
   for(unsigned int i=0;i<numclus;i++)
   {
     if(houghlines[i].size()>stsize)
     {
	stsize= houghlines[i].size();
	starthoughclus=i; 
     }
   }
   /// end sizes:
    for(unsigned int i=0;i<numclusend;i++)
   {
     if(houghlinesend[i].size()>endsize)
     {

	endsize= houghlinesend[i].size();
	endhoughclus=i;
     }
   }
   
//    for(int i=0;i<numclusend;i++)
//    endsize= houghlinesend[i].size();
  
   
   
   //if(endsize>stsize )
   if(endsize/lochitlistendsize  > stsize/lochitlistsize  )
   {
     direction_local=-1; 
     //////std::cout << "!!! Hough direction is different! " << direction_local << std::endl;
   //  numclusfunc=numclusend;
   //  houghutil=&houghlinesend; 
   }
   
  }
   
 
 
 
 ////////std::cout << "new parameters: " << numclusfunc << " " << (*houghutil)[0].size() << std::endl;
 //////std::cout << " direction in:  " << direction << std::endl;
 double houghdirection=1;    //direction as it currently is
 if(wire_start>wire_end)
 { houghdirection=-1; 
   direction_local*=-1;     //make the hough calculated direction absolute.
 }
   
  double startwire=10000*(1+houghdirection),starttime=-1,endwire=10000*(1-houghdirection),endtime=-1;
 
  
  
  //use longest houghline to find new start
  if(endhoughclus!=-1)
  {
    for(unsigned int ix=0;ix<houghlinesend[endhoughclus].size();ix++) {
	unsigned int w;
	w=houghlinesend[endhoughclus][ix]->WireID().Wire; 
   
	if(w*houghdirection>endwire*houghdirection) //multiply by original direction to find most extreme point 
	  {endwire=w;
	  endtime=houghlinesend[endhoughclus][ix]->PeakTime();
	  }
      }
  wire_end=endwire;
  time_end=endtime;
  }
  
  if(starthoughclus!=-1)
  {
    for(unsigned int ix=0;ix<houghlines[starthoughclus].size();ix++) {
	unsigned int w;
	//  unsigned int c = houghlines[starthoughclus][ix]->Wire()->RawDigit()->Channel(); 
	//  geo->ChannelToWire(c,cs,t,p,w);
    
	w=houghlines[starthoughclus][ix]->WireID().Wire; 
   
	if(w*houghdirection<startwire*houghdirection) //multiply by original direction to find most extreme point 
	  {startwire=w;
	  starttime=houghlines[starthoughclus][ix]->PeakTime();
	  }
    }
    wire_start=startwire;
    time_start=starttime;
  }
  
  
 //////std::cout << " low, high wire,time: " << startwire <<","<<starttime << " | " << endwire<<","<<endtime<<std::endl;
 direction=direction_local;  // returning the newfound direction for a decision outside.
 
 

 
}




////////////////////// add override if forcing rightwise:

int cluster::ClusterParamsAlg::DecideClusterDirection(std::vector < art::Ptr<recob::Hit> > hitlist,
				double lineslope,double wstn,double tstn,double wendn,double tendn)  
{
  int fDirection=1;
  double HiBin,LowBin,invHiBin,invLowBin;
  double wire_start,time_start,wire_end,time_end;
  
  if(wstn < wendn)  // make sure we're looking at general direction values.
  {
    wire_start=wstn;
    wire_end=wendn;
    time_start=tstn;
    time_end=tendn;
  }
  else  // if we received a start point further right than the end point, switch locally, so that we can calculate an absolute value.
  {
    wire_start=wendn ;
    wire_end=wstn;
    time_start=tendn;
    time_end=tstn;
  }
  
  //doesn't change wire,time values - assumes it's already the right ones.
  FindDirectionWeights(lineslope,wire_start,time_start,wire_end,time_end,hitlist,HiBin,LowBin,invHiBin,invLowBin); 
  
  
 if(HiBin>LowBin)  // shower is leftward going (contrary to our previous assumptions)
  {
   fDirection=-1;
  }
  
  //this one does refine the start and end point, although it leaves them in the same order.
  // it will now be up to the algorithm here to decide whether we should override the decision  
  // based on FindDirectionWeights
  int houghdirection=fDirection;
  RefineStartPointsHough(hitlist,wire_start,time_start,wire_end,time_end,houghdirection);
  
//  int principaldirection=FindPrincipalDirection(hitlist,wire_start,time_start,wire_end,time_end,lineslope);
  
 //////std::cout << " direction based on principal components. " << principaldirection << std::endl;
  
  //check whether it makes sense to reverse direction based on the findings of RefineStartPointsHough
   if(fDirection!=houghdirection 
     && ( (invLowBin > invHiBin && fDirection==1 ) || (invLowBin < invHiBin && fDirection==-1 ) ))
   {
    //overriding basic with HoughLine Direction, flipping start and end 
   /*     wstn=wire_end;
	tstn=time_end;
	wendn=wire_start;
	tendn=time_start;
   */     fDirection = houghdirection; 
    }
  /*else   // not flipping
   {  // assigning to pick up the hough refined values
      wstn=wire_start;
      tstn=time_start;
      wendn=wire_end;
      tendn=time_end;
   }
  */  
  
  if(fForceRightGoing)
  {
   if(wstn > wendn  || fDirection==-1 )  //change only if decided it's backwards up to now. otherwise it's fine as is.
   {
     //////std::cout << "forcing right going shower" << std::endl;
     fDirection=1;
   
     wstn=wire_end;
     tstn=time_end;
      
     wendn=wire_start;
     tendn=time_start;
   }
   
  }
 

 
 return fDirection;
}




void cluster::ClusterParamsAlg::FindExtremeIntercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*fTimetoCm-perpslope*(double)wire*fWiretoCm;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}




void cluster::ClusterParamsAlg::SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double linearlimit,   double ortlimit, double lineslopetest)
{
  
  double locintercept=time_start-wire_start*lineslopetest;
  
  ////////std::cout << "$$$$$$$$$$$$$ locintercept: " << locintercept << " " << lineinterctest << std::endl;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
	double wonline=wire,tonline=time;
	//gser.GetPointOnLine(lineslopetest,lineinterctest,wire,time,wonline,tonline);
	gser.GetPointOnLine(lineslopetest,locintercept,wire,time,wonline,tonline);
	
	//calculate linear distance from start point and orthogonal distance from axis
	double lindist=gser.Get2DDistance(wonline,tonline,wire_start,time_start);
	double ortdist=gser.Get2DDistance(wire,time,wonline,tonline);
	
	////////std::cout << " w,t: " << wire << " " << time << " ws,ts " << wonline << " "<< tonline <<" "<< lindist << " " << ortdist << std::endl;
	
	if(lindist<linearlimit && ortdist<ortlimit)
	  hitlistlocal.push_back(theHit);
    
    
    }
    
}
///////////////////

void cluster::ClusterParamsAlg::SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double radlimit)
{
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
	//calculate linear distance from start point and orthogonal distance from axis
	double lindist=gser.Get2DDistance(wire,time,wire_start,time_start);
	
	if(lindist<radlimit )
	  hitlistlocal.push_back(theHit);
    
    
    }
    
}

/*
void cluster::ClusterParamsAlg::FindExtremeIntercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*fTimetoCm-perpslope*(double)wire*fWiretoCm;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}*/




art::Ptr<recob::Hit> cluster::ClusterParamsAlg::FindClosestHit(std::vector < art::Ptr < recob::Hit > > hitlist,
			     double wire_online,
			     double time_online)
{
  
  double min_length_from_start=99999;
  art::Ptr<recob::Hit> nearHit;
   
  unsigned int plane,tpc,wire,cstat;
   
   
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
  //  //////std::cout << " in find closest hit " << std::endl;
    double dist_mod=gser.Get2DDistance(wire_online,time_online,wire,time);
    //TMath::Sqrt( pow(((double)wire_online-(double)wire*fWirePitch),2)+pow((time_online-time*fDriftVelocity*fTimeTick),2) );	

    if(dist_mod<min_length_from_start){
	//wire_start[plane]=wire;
	//time_start[plane]=time;
	nearHit=(*hitIter);
	min_length_from_start=dist_mod;
	}	

  } 
  
return nearHit;    
}









// ******************************* //
int cluster::ClusterParamsAlg::GetPlaneAndTPC(art::Ptr<recob::Hit> a,
						unsigned int &p,
						unsigned int &cs,
						unsigned int &t,
						unsigned int &w)
{
   
    
     w=a->WireID().Wire; 
     p=a->WireID().Plane; 
     t=a->WireID().TPC; 
     cs=a->WireID().Cryostat; 
     
  return 0;
}




