/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "GFBookkeeping.h"
#include "GFException.h"
#include "TBuffer.h"
#include "TString.h"
#include <sstream>

genf::GFBookkeeping::GFBookkeeping(const GFBookkeeping& bk) : TObject(bk)
{
  fNhits = bk.fNhits;
  fMatrices = bk.fMatrices; //implicit copy constructor call
  fNumbers = bk.fNumbers;
  fPlanes = bk.fPlanes;
  fFailedHits = bk.fFailedHits;

  //deep copy

  std::map<std::string, TMatrixT<Double_t>*>::const_iterator it;
  std::map<std::string, TMatrixT<Double_t>*>::iterator it_here;

  it = bk.fMatrices.begin();
  it_here = fMatrices.begin();
  while (it != bk.fMatrices.end()) {
    it_here->second = new TMatrixT<Double_t>[fNhits];
    for (int i = 0; i < fNhits; ++i) {
      (it_here->second)[i] = (it->second)[i];
    }
    it++;
    it_here++;
  }

  it = bk.fNumbers.begin();
  it_here = fNumbers.begin();
  while (it != bk.fNumbers.end()) {
    it_here->second = new TMatrixT<Double_t>[fNhits];
    for (int i = 0; i < fNhits; ++i) {
      (it_here->second)[i] = (it->second)[i];
    }
    it++;
    it_here++;
  }

  std::map<std::string, genf::GFDetPlane*>::const_iterator ip;
  std::map<std::string, genf::GFDetPlane*>::iterator ip_here;

  ip = bk.fPlanes.begin();
  ip_here = fPlanes.begin();
  while (ip != bk.fPlanes.end()) {
    ip_here->second = new genf::GFDetPlane[fNhits];
    for (int i = 0; i < fNhits; ++i) {
      (ip_here->second)[i] = ((ip->second)[i]);
    }
    ip++;
    ip_here++;
  }
}

void genf::GFBookkeeping::Streamer(TBuffer& R__b)
{

  // Stream an object of class GFBookkeeping.
  if (R__b.IsReading()) {

    //     Version_t R__v = R__b.ReadVersion();
    TObject::Streamer(R__b);

    clearAll();

    R__b >> fNhits;

    TString s;
    std::string key;
    unsigned int nkeys;
    TMatrixT<Double_t> mat;
    genf::GFDetPlane pl;

    { //reading matrices
      R__b >> nkeys;
      for (unsigned int i = 0; i < nkeys; ++i) {
        s.Streamer(R__b);
        key = s.Data();
        bookMatrices(key);
        for (int j = 0; j < fNhits; ++j) {
          mat.Streamer(R__b);
          setMatrix(key, j, mat);
        }
      }
    } //done reading matrices
    { //reading planes
      R__b >> nkeys;
      for (unsigned int i = 0; i < nkeys; ++i) {
        s.Streamer(R__b);
        key = s.Data();
        bookGFDetPlanes(key);
        for (int j = 0; j < fNhits; ++j) {
          pl.Streamer(R__b);
          setDetPlane(key, j, pl);
        }
      }
    } //done reading planes
    { //reading numbers
      R__b >> nkeys;
      for (unsigned int i = 0; i < nkeys; ++i) {
        s.Streamer(R__b);
        key = s.Data();
        bookNumbers(key);
        for (int j = 0; j < fNhits; ++j) {
          mat.Streamer(R__b);
          setNumber(key, j, mat[0][0]);
        }
      }
    } //done reading numbers
    { //read failed hits
      clearFailedHits();
      unsigned int nFailedHits;
      R__b >> nFailedHits;
      unsigned int aFailedHit;
      for (unsigned int i = 0; i < nFailedHits; ++i) {
        R__b >> aFailedHit;
        fFailedHits.push_back(aFailedHit);
      }
    } //done reading failed hits
  }
  else {
    //     R__b.WriteVersion(genf::GFBookkeeping::IsA());
    TObject::Streamer(R__b);

    //write number of hits
    R__b << fNhits;

    std::vector<std::string> keys;
    { //save matrices
      keys = getMatrixKeys();
      R__b << (unsigned int)(keys.size());
      for (unsigned int i = 0; i < keys.size(); ++i) {
        TString s(keys.at(i));
        s.Streamer(R__b);
        for (int j = 0; j < fNhits; ++j) {
          ((fMatrices[keys.at(i)])[j]).Streamer(R__b);
        }
      }
    }
    keys.clear();
    { //save GFDetPlanes
      keys = getGFDetPlaneKeys();
      R__b << (unsigned int)(keys.size());
      for (unsigned int i = 0; i < keys.size(); ++i) {
        TString s(keys.at(i));
        s.Streamer(R__b);
        for (int j = 0; j < fNhits; ++j) {
          ((fPlanes[keys.at(i)])[j]).Streamer(R__b);
        }
      }
    } //done saving GFDetPlanes
    keys.clear();
    { //save numbers
      keys = getNumberKeys();
      R__b << (unsigned int)(keys.size());
      for (unsigned int i = 0; i < keys.size(); ++i) {
        TString s(keys.at(i));
        s.Streamer(R__b);
        for (int j = 0; j < fNhits; ++j) {
          ((fNumbers[keys.at(i)])[j]).Streamer(R__b);
        }
      }
    } //done saving numbers
    { //save failedHits
      R__b << ((unsigned int)fFailedHits.size());
      for (unsigned int i = 0; i < fFailedHits.size(); ++i) {
        R__b << fFailedHits.at(i);
      }
    } //done saving failed Hits
  }
}

void genf::GFBookkeeping::bookMatrices(std::string key)
{
  if (fNhits < 0) {
    GFException exc("fNhits not defined", __LINE__, __FILE__);
    throw exc;
  }
  if (fMatrices[key] != NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is already occupied in genf::GFBookkeeping::bookMatrices()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  fMatrices[key] = new TMatrixT<Double_t>[fNhits];
}

void genf::GFBookkeeping::bookGFDetPlanes(std::string key)
{
  if (fNhits < 0) {
    GFException exc("fNhits not defined", __LINE__, __FILE__);
    throw exc;
  }
  if (fPlanes[key] != NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is already occupied in genf::GFBookkeeping::bookGFDetPlanes()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  fPlanes[key] = new genf::GFDetPlane[fNhits];
}

//val is default set to 0.
void genf::GFBookkeeping::bookNumbers(std::string key, double val)
{
  if (fNhits < 0) {
    GFException exc("fNhits not defined", __LINE__, __FILE__);
    throw exc;
  }
  if (fPlanes[key] != NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is already occupied in genf::GFBookkeeping::bookNumbers()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  fNumbers[key] = new TMatrixT<Double_t>[fNhits];
  for (int i = 0; i < fNhits; ++i) {
    ((fNumbers[key])[i]).ResizeTo(1, 1);
    ((fNumbers[key])[i])[0][0] = val;
  }
}

void genf::GFBookkeeping::setMatrix(std::string key,
                                    unsigned int index,
                                    const TMatrixT<Double_t>& mat)
{
  if (fMatrices[key] == NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::setMatrix()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::setMatrix()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  (fMatrices[key])[index].ResizeTo(mat);
  (fMatrices[key])[index] = mat;
}
void genf::GFBookkeeping::setDetPlane(std::string key,
                                      unsigned int index,
                                      const genf::GFDetPlane& pl)
{
  if (fPlanes[key] == NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  (fPlanes[key])[index] = pl;
}
void genf::GFBookkeeping::setNumber(std::string key, unsigned int index, const double& num)
{
  if (fNumbers[key] == NULL) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::setNumber()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::setNumber()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  ((fNumbers[key])[index])[0][0] = num;
}

bool genf::GFBookkeeping::getMatrix(std::string key,
                                    unsigned int index,
                                    TMatrixT<Double_t>& mat) const
{
  auto iMatrix = fMatrices.find(key);
  if (iMatrix == fMatrices.end()) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  mat.ResizeTo(iMatrix->second[index]);
  mat = iMatrix->second[index];
  return true;
}
bool genf::GFBookkeeping::getDetPlane(std::string key,
                                      unsigned int index,
                                      genf::GFDetPlane& pl) const
{

  auto iPlane = fPlanes.find(key);
  if (iPlane == fPlanes.end()) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  pl = iPlane->second[index];
  return true;
}
bool genf::GFBookkeeping::getNumber(std::string key, unsigned int index, double& num) const
{
  auto iNumber = fNumbers.find(key);
  if (iNumber == fNumbers.end()) {
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in genf::GFBookkeeping::getNumber()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  if (index >= (unsigned int)fNhits) {
    std::ostringstream ostr;
    ostr << "The index " << index << " is out of range in genf::GFBookkeeping::getNumber()";
    GFException exc(ostr.str(), __LINE__, __FILE__);
    throw exc;
  }
  num = iNumber->second[index][0][0];
  return true;
}

void genf::GFBookkeeping::addFailedHit(unsigned int id)
{
  fFailedHits.push_back(id);
}

unsigned int genf::GFBookkeeping::getNumFailed()
{
  return fFailedHits.size();
}

unsigned int genf::GFBookkeeping::hitFailed(unsigned int id)
{
  unsigned int retVal = 0;
  for (unsigned int i = 0; i < fFailedHits.size(); ++i) {
    if (fFailedHits.at(i) == id) { ++retVal; }
  }
  return retVal;
}

void genf::GFBookkeeping::clearFailedHits()
{
  fFailedHits.clear();
}

void genf::GFBookkeeping::reset()
{
  std::vector<std::string> matKeys = getMatrixKeys();
  std::vector<std::string> planeKeys = getGFDetPlaneKeys();
  std::vector<std::string> numKeys = getNumberKeys();

  clearAll();
  clearFailedHits();

  for (unsigned int i = 0; i < matKeys.size(); ++i) {
    bookMatrices(matKeys.at(i));
  }
  for (unsigned int i = 0; i < planeKeys.size(); ++i) {
    bookGFDetPlanes(planeKeys.at(i));
  }
  for (unsigned int i = 0; i < numKeys.size(); ++i) {
    bookNumbers(numKeys.at(i));
  }
}

void genf::GFBookkeeping::clearAll()
{
  std::map<std::string, TMatrixT<Double_t>*>::iterator itMat;
  for (itMat = fMatrices.begin(); itMat != fMatrices.end(); itMat++) {
    if (itMat->second != NULL) delete[] itMat->second;
  }
  std::map<std::string, genf::GFDetPlane*>::iterator itPl;
  for (itPl = fPlanes.begin(); itPl != fPlanes.end(); itPl++) {
    if (itPl->second != NULL) delete[] itPl->second;
  }
  std::map<std::string, TMatrixT<Double_t>*>::iterator itNum;
  for (itNum = fNumbers.begin(); itNum != fNumbers.end(); itNum++) {
    if (itNum->second != NULL) delete[] itNum->second;
  }
  fMatrices.clear();
  fPlanes.clear();
  fNumbers.clear();
}

std::vector<std::string> genf::GFBookkeeping::getMatrixKeys() const
{
  std::vector<std::string> keys;
  std::map<std::string, TMatrixT<Double_t>*>::const_iterator it;
  for (it = fMatrices.begin(); it != fMatrices.end(); it++) {
    if (it->second != NULL) keys.push_back(it->first);
  }
  return keys;
}
std::vector<std::string> genf::GFBookkeeping::getGFDetPlaneKeys() const
{
  std::vector<std::string> keys;
  std::map<std::string, genf::GFDetPlane*>::const_iterator it;
  for (it = fPlanes.begin(); it != fPlanes.end(); it++) {
    if (it->second != NULL) keys.push_back(it->first);
  }
  return keys;
}
std::vector<std::string> genf::GFBookkeeping::getNumberKeys() const
{
  std::vector<std::string> keys;
  std::map<std::string, TMatrixT<Double_t>*>::const_iterator it;
  for (it = fNumbers.begin(); it != fNumbers.end(); it++) {
    if (it->second != NULL) keys.push_back(it->first);
  }
  return keys;
}

void genf::GFBookkeeping::Print(std::ostream& out /* = std::cout */) const
{
  out << "=============genf::GFBookkeeping::print()==============" << std::endl;
  out << "-----printing all matrices:------" << std::endl;
  std::vector<std::string> keys = getMatrixKeys();
  for (unsigned int i = 0; i < keys.size(); ++i) {
    out << "key " << keys.at(i) << " has " << fNhits << " entries:" << std::endl;
    for (int j = 0; j < fNhits; ++j) {
      TMatrixT<Double_t> m;
      getMatrix(keys.at(i), j, m);
      m.Print(); // TODO print the matrix outself
    }
  }
  out << "-----printing all GFDetPlanes:------" << std::endl;
  keys = getGFDetPlaneKeys();
  for (unsigned int i = 0; i < keys.size(); ++i) {
    out << "key " << keys.at(i) << " has " << fNhits << " entries:" << std::endl;
    for (int j = 0; j < fNhits; ++j) {
      genf::GFDetPlane p;
      getDetPlane(keys.at(i), j, p);
      p.Print(out);
    }
  }
  out << "-----printing all numbers:------" << std::endl;
  keys = getNumberKeys();
  for (unsigned int i = 0; i < keys.size(); ++i) {
    out << "key " << keys.at(i) << " has " << fNhits << " entries:" << std::endl;
    for (int j = 0; j < fNhits; ++j) {
      double n(-1111.);
      getNumber(keys.at(i), j, n);
      out << n << std::endl;
    }
  }
  out << "-----failed hits:------" << std::endl;
  for (unsigned int i = 0; i < fFailedHits.size(); ++i) {
    out << fFailedHits.at(i) << " ";
  }
  out << std::endl;
}

//ClassImp(GFBookkeeping)
