////////////////////////////////////////////////////////////////////////
//
// HarrisVertexFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding. 
//  A weak vertex is a vertex that has been found using a dedicated vertex finding algorithm only. A 
//  strong vertex is a vertex that has been found using a dedicated vertex finding algorithm and matched 
//  to a crossing of two or more HoughLineFinder lines. The VertexMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey 
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors". 
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions
#ifndef HarrisVertexFinder_H
#define HarrisVertexFinder_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TMath.h"
#include "TH2.h"

///vertex reconstruction
namespace vertex {
   
 class HarrisVertexFinder :  public art::EDProducer {
    
  public:
    
   explicit HarrisVertexFinder(fhicl::ParameterSet const& pset); 
   virtual ~HarrisVertexFinder();        
   void beginJob();
   
   void produce(art::Event& evt);
   
  private:

   double Gaussian(int x, int y, double sigma);
   double GaussianDerivativeX(int x, int y);
   double GaussianDerivativeY(int x, int y);
   void   SaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy);

   std::string fDBScanModuleLabel;
   int         fTimeBins;
   int         fMaxCorners;
   double      fGsigma;
   int         fWindow;
   double      fThreshold;
   int         fSaveVertexMap;
   TH2F*       fNoVertices;
   
  };
    
}

//-----------------------------------------------------------------------------
vertex::HarrisVertexFinder::HarrisVertexFinder(fhicl::ParameterSet const& pset) 
  : fDBScanModuleLabel  (pset.get< std::string >("DBScanModuleLabel"))
  , fTimeBins        (pset.get< int         >("TimeBins")      )
  , fMaxCorners      (pset.get< int         >("MaxCorners")    )
  , fGsigma          (pset.get< double      >("Gsigma")        )
  , fWindow          (pset.get< int         >("Window")        )
  , fThreshold       (pset.get< double      >("Threshold")     )
  , fSaveVertexMap   (pset.get< int         >("SaveVertexMap") )
{
  produces< std::vector<recob::EndPoint2D> >();
  produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
}

//-----------------------------------------------------------------------------
vertex::HarrisVertexFinder::~HarrisVertexFinder()
{
}

void vertex::HarrisVertexFinder::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fNoVertices= tfs->make<TH2F>("fNoVertices", ";Event No; No of vertices", 100,0, 100, 30, 0, 30);

}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::Gaussian(int x, int y, double sigma)
{
  double Norm  = 1./std::sqrt(2*TMath::Pi()*pow(sigma,2));
  double value = Norm*exp(-(pow(x,2)+pow(y,2))/(2*pow(sigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::GaussianDerivativeX(int x,int y)
{
  double Norm  = 1./(std::sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value = Norm*(-x)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::GaussianDerivativeY(int x,int y)
{
  double Norm  = 1./(std::sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value = Norm*(-y)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
void vertex::HarrisVertexFinder::produce(art::Event& evt)
{

  art::ServiceHandle<geo::Geometry> geom;
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel, clusterListHandle);

  //Point to a collection of vertices to output and the associations with the hits
  std::unique_ptr<std::vector<recob::EndPoint2D> > vtxcol(new std::vector<recob::EndPoint2D>);
  std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assn(new art::Assns<recob::EndPoint2D, recob::Hit>);

  std::vector< art::Ptr<recob::Hit> > cHits;
  art::PtrVector<recob::Hit> hit;
   
  art::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){
    art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
    clusIn.push_back(cluster);
  }

  art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fDBScanModuleLabel);

  int flag   = 0;
  int n      = 0; //index of window cell. There are 49 cells in the 7X7 Gaussian and Gaussian derivative windows
  unsigned int numberwires = 0;
  double MatrixAAsum = 0.;
  double MatrixBBsum = 0.; 
  double MatrixCCsum = 0.;
  std::vector<double> Cornerness2;
  //gaussian window definitions. The cell weights are calculated here to help the algorithm's speed
  double  w[49] = {0.};
  double wx[49] = {0.};
  double wy[49] = {0.};
  int ctr = 0;
  for(int i = -3; i < 4; ++i){
    for(int j = 3; j > -4; --j){
      w[ctr] = Gaussian(i, j, fGsigma);
      wx[ctr] = GaussianDerivativeX(i,j);
      wy[ctr] = GaussianDerivativeY(i,j);
      ++ctr;
    }
  }
  
  const unsigned int numbertimesamples
    = lar::providerFrom<detinfo::DetectorPropertiesService>()->ReadOutWindowSize();
  
  const float BinsPerTick = fTimeBins / numbertimesamples;
  const float TicksPerBin = numbertimesamples / fTimeBins;
  
  for(auto pid : geom->PlaneIDs()){
    art::PtrVector<recob::Hit> vHits;
    geo::View_t view = geom->View(pid);
    hit.clear();
    cHits.clear();      
    for(size_t ci = 0; ci < clusIn.size(); ++ci) {
      cHits = fmh.at(ci);
      if(cHits.size() > 0)
	for(unsigned int i = 0; i < cHits.size(); ++i) hit.push_back(cHits[i]);
      
    } 
    if(hit.size() == 0) continue;
    
    numberwires       = geom->Cryostat(pid.Cryostat).TPC(pid.TPC).Plane(pid.Plane).Nwires();
    std::vector< std::vector<double> > MatrixAsum(numberwires);
    std::vector< std::vector<double> > MatrixBsum(numberwires);
    std::vector< std::vector<double> > hit_map(numberwires);//the map of hits 
    std::vector< std::vector<int>    > hit_loc(numberwires);//the index of the hit that corresponds to the potential corner
    std::vector< std::vector<double> > Cornerness(numberwires);//the "weight" of a corner
    
    for(unsigned int wi = 0; wi < numberwires; ++wi){
      MatrixAsum[wi].resize(fTimeBins);
      MatrixBsum[wi].resize(fTimeBins);
      hit_map[wi].resize(fTimeBins);
      hit_loc[wi].resize(fTimeBins);
      Cornerness[wi].resize(fTimeBins);
      for(int timebin = 0; timebin < fTimeBins; ++timebin){
	hit_map[wi][timebin] = 0.;
	hit_loc[wi][timebin] = -1;
	Cornerness[wi][timebin] = 0.;
	MatrixAsum[wi][timebin] = 0.;
	MatrixBsum[wi][timebin] = 0.;
      }
    }
    
    for(unsigned int i = 0; i < hit.size(); ++i){
      unsigned int wire = hit[i]->WireID().Wire;
      // pixelization using a Gaussian (3 sigmas around peak)
      const float center = hit[i]->PeakTime(), sigma = hit[i]->RMS();
      const int iFirstBin = int((center - 3*sigma) / TicksPerBin),
        iLastBin = int((center + 3*sigma) / TicksPerBin);
      for (int iBin = iFirstBin; iBin <= iLastBin; ++iBin) {
        const float bin_center = iBin * TicksPerBin;
        hit_map[wire][iBin] += Gaussian(bin_center, center, sigma);
      }
    }
    
    // Gaussian derivative convolution  
    for(unsigned int wire = 1;wire < numberwires-1; ++wire){
      for(int timebin = 1;timebin < fTimeBins-1; ++timebin){
	MatrixAsum[wire][timebin] = 0.;
	MatrixBsum[wire][timebin] = 0.;
	n = 0;
	for(int i = -3; i <= 3; ++i) {
	  unsigned int windex = ((int)wire < -i)? 0: wire + i;
	  if (windex >= numberwires) windex = numberwires - 1;
	  for(int j = -3; j <= 3; ++j) {
	    unsigned int tindex = (unsigned int) std::min(std::max(timebin +j, 0), fTimeBins - 1);
	    
	    MatrixAsum[wire][timebin] += wx[n]*hit_map[windex][tindex];
	    MatrixBsum[wire][timebin] += wy[n]*hit_map[windex][tindex];
	    ++n;
	  }// end loop over j
	}// end loop over i
      }// end loop over time bins
    }// end loop over wires
    
    //calculate the cornerness of each pixel while making sure not to fall off the hit map.
    for(unsigned int wire = 1; wire < numberwires-1; ++wire){
      for(int timebin = 1; timebin < fTimeBins-1; ++timebin){     
	MatrixAAsum = 0.;
	MatrixBBsum = 0.;
	MatrixCCsum = 0.;
	//Gaussian smoothing convolution
	n = 0;
	
	for(int i = -3; i <= 3; i++) {
	  unsigned int windex = ((int)wire < -i)? 0: wire + i;
	  if (windex >= numberwires) windex = numberwires - 1;
	  for(int j = -3; j <= 3; j++) {
	    unsigned int tindex = (unsigned int) std::min(std::max(timebin +j, 0), fTimeBins - 1);
	    MatrixAAsum += w[n]*pow(MatrixAsum[windex][tindex],2);
	    MatrixBBsum += w[n]*pow(MatrixBsum[windex][tindex],2);
	    MatrixCCsum += w[n]*MatrixAsum[windex][tindex]*MatrixBsum[windex][tindex]; 
	    ++n;
	  }// end loop over j
	}// end loop over i
	
	if((MatrixAAsum+MatrixBBsum) > 0)		
	  Cornerness[wire][timebin] = (MatrixAAsum*MatrixBBsum-pow(MatrixCCsum,2))/(MatrixAAsum+MatrixBBsum);
	else
	  Cornerness[wire][timebin] = 0;
	
	if(Cornerness[wire][timebin]>0){	  
	  for(unsigned int i = 0;i < hit.size(); ++i){
	    unsigned int wire2 = hit[i]->WireID().Wire;
	    // make sure the vertex candidate coincides with an actual hit
	    // (time within one RMS from the center of the hit)
	    const float time = timebin / TicksPerBin;
	    if(wire == wire2 && std::abs(hit[i]->TimeDistanceAsRMS(time)) < 1.){
	      //this index keeps track of the hit number
	      hit_loc[wire][timebin] = i;
	      Cornerness2.push_back(Cornerness[wire][timebin]);
	      break;
	    } 	        
	  }// end loop over hits	     
	}// end if cornerness > 0	    
      }// end loop over time bins
    }// end loop over wires
    
    std::sort(Cornerness2.rbegin(),Cornerness2.rend());
    
    for(int vertexnum = 0; vertexnum < fMaxCorners; ++vertexnum){
      flag = 0;
      for(unsigned int wire = 0;wire < numberwires && flag == 0; ++wire){
	for(int timebin = 0;timebin < fTimeBins && flag == 0; ++timebin){    
	  if(Cornerness2.size() > (unsigned int)vertexnum)
	    if(Cornerness[wire][timebin] == Cornerness2[vertexnum] 
	       && Cornerness[wire][timebin] > 0. 
	       && hit_loc[wire][timebin]>-1){
	      ++flag;      
	      //thresholding
	      if(Cornerness2.size())
		if(Cornerness[wire][timebin] < (fThreshold*Cornerness2[0]))
		  vertexnum = fMaxCorners;
	      
	      vHits.push_back(hit[hit_loc[wire][timebin]]);
	      // get the total charge from the associated hits
	      double totalQ = 0.;
	      for(size_t vh = 0; vh < vHits.size(); ++vh) totalQ += vHits[vh]->Integral();
	      
	      recob::EndPoint2D vertex(hit[hit_loc[wire][timebin]]->PeakTime(),
				       hit[hit_loc[wire][timebin]]->WireID(),  //for update to EndPoint2D ... WK 4/22/13
				       Cornerness[wire][timebin],
				       vtxcol->size(),
				       view,
				       totalQ);
	      vtxcol->push_back(vertex);
	      
	      util::CreateAssn(*this, evt, *(vtxcol.get()), vHits, *(assn.get()));
	      
	      vHits.clear();
	      // non-maximal suppression on a square window. The wire coordinate units are 
	      // converted to time ticks so that the window is truly square. 
	      // Note that there are 1/0.0743=13.46 time samples per 4.0 mm (wire pitch in ArgoNeuT), 
	      // assuming a 1.5 mm/us drift velocity for a 500 V/cm E-field 
	      
	      /// \todo this is completely specific to ArgoNeuT!!!!
	      
	      for(unsigned int wireout=wire-(int)((fWindow*BinsPerTick*.0743)+.5);
		  wireout <= wire+(int)((fWindow*BinsPerTick*.0743)+.5) ; wireout++)
		for(int timebinout=timebin-fWindow;timebinout <= timebin+fWindow; timebinout++)
		  if(std::sqrt(pow(wire-wireout,2)+pow(timebin-timebinout,2))<fWindow)//circular window 
		    Cornerness[wireout][timebinout]=0;	  
	    }
	}// end loop over time bins     
      }// end loop over wires
    }// end loop over vertices
    Cornerness2.clear();
    hit.clear();
    
    if(pid.Plane == (unsigned int)fSaveVertexMap){ 
      unsigned char *outPix = new unsigned char [fTimeBins*numberwires];
      //finds the maximum cell in the map for image scaling
      int cell, pix=0, maxCell=0;
      for (int y = 0; y < fTimeBins; ++y){
	for (unsigned int x = 0; x < numberwires; ++x){
	  cell = (int)(hit_map[x][y]*1000);
	  if (cell > maxCell){
	    maxCell = cell;
	  }
	}
      }
      
      for (unsigned int y = 0; y < (unsigned int)fTimeBins; ++y){
	for (unsigned int x = 0; x < numberwires; ++x){ 
	  //scales the pixel weights based on the maximum cell value     
	  if(maxCell>0)
	    pix = (int)((1500000*hit_map[x][y])/maxCell);
	  outPix[y*numberwires + x] = pix;
	}
      }
      SaveBMPFile("harrisvertexmap.bmp", outPix, numberwires, fTimeBins);
      delete [] outPix;
    }   
  }// end loop over planes

  mf::LogInfo("HarrisVertexFinder") << "Size of vtxcol= " << vtxcol->size();
  fNoVertices->Fill(evt.id().event(),vtxcol->size());
  
  evt.put(std::move(vtxcol));   
}

//-----------------------------------------------------------------------------
//this method saves a BMP image of the vertex map space, which can be viewed with gimp
void vertex::HarrisVertexFinder::SaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
{
  std::ofstream bmpFile(fileName, std::ios::binary);
  bmpFile.write("B", 1);
  bmpFile.write("M", 1);
  int bitsOffset = 54 +256*4; 
  int size = bitsOffset + dx*dy; //header plus 256 entry LUT plus pixels
  bmpFile.write((const char *)&size, 4);
  int reserved = 0;
  bmpFile.write((const char *)&reserved, 4);
  bmpFile.write((const char *)&bitsOffset, 4);
  int bmiSize = 40;
  bmpFile.write((const char *)&bmiSize, 4);
  bmpFile.write((const char *)&dx, 4);
  bmpFile.write((const char *)&dy, 4);
  short planes = 1;
  bmpFile.write((const char *)&planes, 2);
  short bitCount = 8;
  bmpFile.write((const char *)&bitCount, 2);
  int i, temp = 0;
  for (i=0; i<6; i++)
    bmpFile.write((const char *)&temp, 4);  // zero out optional color info
  // write a linear LUT
  char lutEntry[4]; // blue,green,red
  lutEntry[3] = 0;  // reserved part
  for (i=0; i<256; i++)
    {
      lutEntry[0] = lutEntry[1] = lutEntry[2] = i;
      bmpFile.write(lutEntry, sizeof lutEntry);
    }
  // write the actual pixels
  bmpFile.write((const char *)pix, dx*dy);
}

namespace vertex{

  DEFINE_ART_MODULE(HarrisVertexFinder)

} // end of vertex namespace

#endif // HarrisVertexFinder_H
