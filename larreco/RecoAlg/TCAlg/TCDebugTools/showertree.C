// ***************************************************
// trajcluster debug viewer
// 
// Rory Fitzpatrick
// 
// ***************************************************

#include "showertree.h"
#include <sys/stat.h>

using namespace std;

int getColor(int showerID) {
  int showerColor = 2 + showerID;
  if (showerColor == 10) showerColor = 28;
  else if (showerColor > 10) showerColor += 27; // start again with color 38
  return showerColor;
}

void showertree::Loop() {
  gROOT->SetBatch();
  gStyle->SetOptStat(0);
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "NUMBER OF ENTRIES IN TREE: " << nentries << endl;
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    fChain->LoadTree(jentry);
    fChain->GetEntry(jentry);
    
    const int nstages = nStages;
    const int nplanes = nPlanes;
    
    if (nstages == 0) {
      cout << "Error: There is no shower tree information saved in this file." << endl;
      return;
    }
    
    // define histograms
    //TH2F* TPCview[nplanes][nstages];
    std::vector<std::vector<TH2F*>> TPCview(nplanes);
    for (size_t i = 0; i<TPCview.size(); ++i){
      TPCview[i].reserve(nstages);
    }
    vector<vector<vector<TLine*>>> tjlines(nplanes);  
    for (int i = 0; i < nplanes; ++i) {
      for (int j = 0; j < nstages; ++j) {
	stringstream name, title;
	name << "pln" << i << "_" << j+1 << "_" << event;
	title << "Plane " << i << ", Stage " << j+1 << ", " << StageName->at(j) << ";wire;tick";
	TPCview[i][j] = new TH2F(name.str().c_str(), title.str().c_str(), 500, 0, 5000, 500, 0, 10000);
	vector<TLine*> tempvec;
	tjlines[i].push_back(tempvec);
      }
    } // end hist definitions
    
    // loop through all trajectories
    for (unsigned int i = 0; i < TjID->size(); ++i) {
      // shower trajectories aren't real. we don't want to draw them
      if (IsShowerTj->at(i) == 1) continue;
      
      int pl = (int)PlaneNum->at(i);
      int st = (int)StageNum->at(i) - 1;
      
      float x1 = BeginWir->at(i);
      float y1 = BeginTim->at(i);
      float x2 = EndWir->at(i);
      float y2 = EndTim->at(i);
      
      TPCview[pl][st]->Fill(x1, y1);
      TPCview[pl][st]->Fill(x2, y2);
      
      TLine* tj = new TLine(x1, y1, x2, y2);
      tj->SetLineColor(getColor(ShowerID->at(i)));
      tj->SetLineWidth(2);
      
      // check for parent and update draw settings accordingly
      if (IsShowerParent->at(i) == 1) {
	tj->SetLineStyle(3);
	tj->SetLineWidth(4);
      }
      
      tjlines[pl][st].push_back(tj);
      
    } // end tj loop
    
    // **** draw ****//
    TCanvas** bothViews = new TCanvas*[nstages];
    // vars for making envelopes
    float xx[5], yy[5];
    vector<TPolyLine*> senv;
    
    for (int i = 0; i < nstages; ++i) {
      stringstream name;
      name << "c" << i << "_" << event;
      bothViews[i] = new TCanvas(name.str().c_str(), name.str().c_str(), 800, 800);
      bothViews[i]->Divide(1,nplanes);
      
      for (int j = 0; j < nplanes; ++j) {	
	senv.clear();
	int vtxcnt = 0;

	// get envelope:
	for (int k = 0; k < (int)Envelope->size(); ++k) {
	  
	  if (EnvStage->at(k) != i+1) continue;
	  if (EnvPlane->at(k) != j) continue;
	  
	  if (k%2 == 0) xx[(int)floor(vtxcnt/2)] = Envelope->at(k);
	  else yy[(int)floor(vtxcnt/2)] = Envelope->at(k);
	  
	  vtxcnt++;
	  
	  if (vtxcnt == 8) {
	    xx[4] = xx[0];  yy[4] = yy[0];
	    TPolyLine* env = new TPolyLine(5, xx, yy);
	    env->SetLineColor(getColor(EnvShowerID->at(k)));
	    senv.push_back(env);
	    vtxcnt = 0;
	    continue;
	  } 
	} // end envelope loop
	
	bothViews[i]->cd(nplanes-j);
	// scale the canvase to the region of interest
	int xmin, xmax, ymin, ymax;
	xmin = TPCview[j][i]->GetXaxis()->GetBinCenter( TPCview[j][i]->FindFirstBinAbove(0,1) );
	xmax = TPCview[j][i]->GetXaxis()->GetBinCenter( TPCview[j][i]->FindLastBinAbove(0,1) );
	ymin = TPCview[j][i]->GetYaxis()->GetBinCenter( TPCview[j][i]->FindFirstBinAbove(0,2) );
	ymax = TPCview[j][i]->GetYaxis()->GetBinCenter( TPCview[j][i]->FindLastBinAbove(0,2) );
	
	TPCview[j][i]->GetXaxis()->SetRangeUser(0.9*xmin, 1.1*xmax);
	TPCview[j][i]->GetYaxis()->SetRangeUser(0.9*ymin, 1.1*ymax);

	TPCview[j][i]->SetMarkerColor(0); // we don't want to see the points in the TH2s
	
	TPCview[j][i]->Draw();
	// draw envelops
	for (int k = 0; k < (int)senv.size(); ++k) {
	  senv.at(k)->Draw("SAME");
	}
	// draw trajectories
	for (int k = 0; k < (int)tjlines[j][i].size(); ++k) { 
	  tjlines[j][i].at(k)->Draw("SAME");
	}
	
      } // nplanes
    } // nstages

    // save output
    struct stat sb;
    if (!(stat("rsf_evd_output", &sb) == 0 && S_ISDIR(sb.st_mode))) {
      cout << "Creating folder 'rsf_evd_output'" << endl;
      mkdir("rsf_evd_output", 0777);
    }
    bothViews[0]->SaveAs(Form("rsf_evd_output/tjs_%d.pdf[", event));
    for (int i = 0; i < nstages; ++i) {
      bothViews[i]->SaveAs(Form("rsf_evd_output/tjs_%d.pdf", event));
    }
    bothViews[nstages-1]->SaveAs(Form("rsf_evd_output/tjs_%d.pdf]", event));
    
  }  // end event loop
  
} // Loop

int main(int argc, char *argv[]) {
  
  if (!argv[1]) {
    cout << "Error: No input file was provided! Exiting." << endl;
    return 0;
  }
  
  string infile = argv[1];
  TTree* shtree;
  
  TFile* f = new TFile(infile.c_str());
  if (f->IsZombie()) {
    cout << "Error: Could not open input file. Exiting." << endl;
    return 0;
  }
  
  TDirectory * dir = (TDirectory*)f->Get(Form("%s:/trajcluster", infile.c_str()));
  dir->GetObject("showervarstree",shtree);
  
  showertree t(shtree);
  
  t.Loop();
} // main
