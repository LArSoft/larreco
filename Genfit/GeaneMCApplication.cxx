#include "Genfit/GeaneMCApplication.h"
#include "TVirtualMC.h"
#include <iostream>
#include "TSystem.h"
#include"assert.h"
#include"GFFieldManager.h"



genf::GeaneMCApplication::GeaneMCApplication()
{
}

void genf::GeaneMCApplication::ConstructGeometry(){
  gGeoManager->CloseGeometry();
  gMC->SetRootGeometry();  
}

void genf::GeaneMCApplication::InitMC(){
  gMC->Init();
  gMC->BuildPhysics();
}

void genf::GeaneMCApplication::Field(const Double_t* x, Double_t* b) const
{
  //assert(field!=NULL);
  TVector3 pos(x[0],x[1],x[2]);
  //  TVector3 B = FieldManager::getInstance()->getField()->get(pos);
  TVector3 B = GFFieldManager::getFieldVal(pos);
  //  std::cout << "GeaneMCApplication::Field()" << std::endl;
  //B.Print();
  b[0]=B.X();
  b[1]=B.Y();
  b[2]=B.Z();
}


//ClassImp(GeaneMCApplication)


