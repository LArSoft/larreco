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
/** @addtogroup genfit
 * @{
 */

#ifndef GFTRACK_H
#define GFTRACK_H

#include <vector>

#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/GFException.h"

#include "TObjArray.h"

#include "larreco/Genfit/GFBookkeeping.h"
#include "larreco/Genfit/GFDetPlane.h"
#include "larreco/Genfit/GFTrackCand.h"

class TVirtualGeoTrack;

namespace genf {
  class GFAbsRecoHit;
}

/** @brief Track object for genfit. genfit algorithms work on these objects.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * Can be used as transient (only in memory) or
 * persistent (written to ROOT-file) object.
 *
 * A GFTrack contains a collection of RecoHits plus a collection of
 * track representations. The GFTrackCand member is a helper object to store the
 * indices of the hits in the GFTrack.
 *
 * For a GFTrack one so called "cardinal representation" can be defined. It is
 * that track representation that is used to access the fit results. Usually
 * one will after the fit choose the best fitting representation to be
 * the cardinal rep.
 *
 * The GFTrack takes ownership over the GFAbsRecoHit pointers it holds.
 */

namespace genf {

  class GFTrack : public TObject {
  private:
    /** @brief Collection of track representations
   *
   * this array is only to be added to in the addTrackRep method
   * because the synchronized construction of bookkeeping objects
   * and repAtHit array is ensured there. NEVER delete elements from
   * this array!
   * If this functionality will be need, it has to be done synchronized
   * with bookkeeping!!
   */
    TObjArray* fTrackReps; //->

    /** @brief Collection of RecoHits
   */
    std::vector<GFAbsRecoHit*> fHits;
    //!

    /** @brief Collection of Bookeeping objects for failed hits
   * in every trackrep
   */
    std::vector<GFBookkeeping*> fBookkeeping;

    /** @brief repAtHit keeps track of at which hit index which rep
   * is currently defined, to avoid null extrapolations
   */
    std::vector<int> fRepAtHit;

    /** @brief Helper to store the indices of the hits in the track.
   * See GFTrackCand for details.
   */
    GFTrackCand fCand; // list of hits
    // Cov and Updates at each step (on last iteration!) as calculated
    // in GFKalman.cxx.
    std::vector<TMatrixT<Double_t>> fHitMeasCov;
    std::vector<TMatrixT<Double_t>> fHitUpdate;
    std::vector<TMatrixT<Double_t>> fHitState;
    std::vector<Double_t> fHitChi2;
    std::vector<TMatrixT<Double_t>> fHitCov7x7;
    std::vector<TMatrixT<Double_t>> fHitCov;
    std::vector<TVector3> fHitPlaneXYZ;
    std::vector<TVector3> fHitPlaneUxUyUz;
    std::vector<TVector3> fHitPlaneU;
    std::vector<TVector3> fHitPlaneV;

    int fPDG;

    static const int fDefNumTrackReps = 10; //!
    unsigned int fCardinal_rep;             // THE selected rep, default=0;

    unsigned int fNextHitToFit;

  public:
    /** @brief Default constructor -- needed for compatibility with ROOT */
    GFTrack();

    /** @brief Copy constructor */
    GFTrack(const GFTrack&);

    /** @brief assignement operator */
    GFTrack& operator=(const GFTrack&);

    /** @brief Initializing constructor
   *
   * A track needs at least one track representation to be useable.
   * The given track representation has to contain starting values for the fit!
   */
    GFTrack(GFAbsTrackRep*);

    virtual ~GFTrack();

    // -----------------------
    // Accessors
    // -----------------------

    /** @brief Resets the GFTrack -- deletes RecoHits!
   */
    void reset(); // deletes the RecoHits!

    /** @brief return the number of failed Hits in track fit
   * repId == -1 will use cardinal rep
  */
    int getFailedHits(int repId = -1)
    {
      int theRep;
      if (repId == -1)
        theRep = fCardinal_rep;
      else
        theRep = repId;
      return fBookkeeping.at(theRep)->getNumFailed();
    }

    std::vector<GFAbsRecoHit*> getHits() { return fHits; }

    const GFTrackCand& getCand() const { return fCand; }

    GFAbsRecoHit* getHit(int id) const { return fHits.at(id); }

    unsigned int getNumHits() const { return fHits.size(); }

    /** @brief Merge two GFTracks. Only hits will be merged.
   *
   * All hits from trk will be merged into this GFTrack.
   * trk will be empty afterwards.
   *
   * Kalman::continueTrack can be used to include the newly added hits
   * in the fit.
   *
   * Note that the new hits are inserted at the end of the present track!
   */
    void mergeHits(GFTrack* trk);

    /** @brief Clear hit vector. Note that hits will not be deleted!
   *
   * Be carefull not to create memory leaks here.
   */
    void releaseHits() { fHits.clear(); }

    /** @brief Accessor for fNextHitToFit
   */
    unsigned int getNextHitToFit() const { return fNextHitToFit; }

    /** @brief Set next hit to be used in a fit
   */
    void setNextHitToFit(unsigned int i) { fNextHitToFit = i; }

    /** @brief Accessor for track representations
   */
    GFAbsTrackRep* getTrackRep(int id) const
    {
      return reinterpret_cast<GFAbsTrackRep*>(fTrackReps->At(id));
    }

    /** @brief Get number of track represenatations
   */
    unsigned int getNumReps() const { return fTrackReps->GetEntriesFast(); }

    /** @brief Get cardinal track representation
   *
   * The user has to choose which track rep should be considered the
   * best one after the fit. Usually the track representation giving the
   * smallest chi2 is choosen. By default the first in the list is returned.
   */
    GFAbsTrackRep* getCardinalRep() const
    {
      return ((GFAbsTrackRep*)fTrackReps->At(fCardinal_rep));
    }

    /** @brief Get momentum at the present position
   *
   * Cardinal representation is used.
   */
    TVector3 getMom() const { return getCardinalRep()->getMom(); }

    /** @brief Get momentum at GFDetPlane
   *
   * The track will be extrapolated to GFDetPlane to get the momentum there.
   * The track will not be modified. Cardinal representation is used.
   */
    TVector3 getMom(const GFDetPlane& pl) const { return getCardinalRep()->getMom(pl); }

    /** @brief Get present position
   *
   * Cardinal representation is used.
   */
    TVector3 getPos() const { return getCardinalRep()->getPos(); }

    /** @brief Get position at GFDetPlane
   *
   * The track will be extrapolated to GFDetPlane to get the position there.
   * The track will not be modified. Cardinal representation is used.
   */
    TVector3 getPos(const GFDetPlane& pl) const { return getCardinalRep()->getPos(pl); }

    /** @brief Get position, momentum, and 6x6 covariance at current position
   *
   * Cardinal representation is used.
   */
    void getPosMomCov(TVector3& pos, TVector3& mom, TMatrixT<Double_t>& cov)
    {
      getCardinalRep()->getPosMomCov(pos, mom, cov);
    }

    /** @brief Get position, momentum, and 6x6 covariance at GFDetPlane
   *
   * The track will be extrapolated to GFDetPlane to get everything there.
   * The track will not be modified. Cardinal representation is used.
   */
    void getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixT<Double_t>& cov)
    {
      getCardinalRep()->getPosMomCov(pl, pos, mom, cov);
    }

    /** @brief Get chi2
   *
   * Cardinal representation is used.
   */
    double getChiSqu() const { return getCardinalRep()->getChiSqu(); }

    /** @brief Get NDF
   *
   * Cardinal representation is used.
   */
    unsigned int getNDF() const { return getCardinalRep()->getNDF(); }

    /** @brief Get chi2/NDF
   *
   * Cardinal representation is used.
   */
    double getRedChiSqu() const { return getCardinalRep()->getRedChiSqu(); }

    /** @brief Get charge from fit
   *
   * Cardinal representation is used.
   */
    double getCharge() const { return getCardinalRep()->getCharge(); }

    /** @brief Fill TVirtualGeoTrack object Cardinal representation is used.
   */
    void fillGeoTrack(TVirtualGeoTrack* tr) const { fillGeoTrack(tr, fCardinal_rep); }

    /** @brief Fill TVirtualGeoTrack object with data from specified track rep
   */
    void fillGeoTrack(TVirtualGeoTrack* tr, unsigned int repid) const;

    // ---------------------
    // Modifiers
    // ---------------------

    void addFailedHit(unsigned int irep, unsigned int id)
    {
      fBookkeeping.at(irep)->addFailedHit(id);
    }

    /** @brief deprecated!
   */
    inline void addHit(GFAbsRecoHit* theHit) { fHits.push_back(theHit); }

    /** @brief Add single hit. Updates the GFTrackCand
   */
    void addHit(GFAbsRecoHit* theHit,
                unsigned int detId,
                unsigned int hitId,
                double rho = 0.,
                unsigned int planeId = 0)
    {
      fHits.push_back(theHit);
      fCand.addHit(detId, hitId, rho, planeId);
    }

    /** @brief Add collection of hits
   *
   * This is the standard way to fill the track with hit data
   */
    void addHitVector(std::vector<GFAbsRecoHit*> hits) { fHits = hits; }

    /** @brief Add track represenation
   *
   * The given track representation has to contain starting values for fit!
   */
    void addTrackRep(GFAbsTrackRep* theTrackRep)
    {
      if (fTrackReps == NULL) fTrackReps = new TObjArray(fDefNumTrackReps);
      fTrackReps->Add(theTrackRep);
      fBookkeeping.push_back(new GFBookkeeping());
      fRepAtHit.push_back(-1);
    }

    //! get GFBookKeeping object for particular track rep (default is cardinal rep)
    GFBookkeeping* getBK(int index = -1)
    {
      if (index == -1) return fBookkeeping.at(fCardinal_rep);
      if ((unsigned int)index >= getNumReps())
        throw GFException("genf::GFTrack::getBK(): index out of range", __LINE__, __FILE__)
          .setFatal();
      return fBookkeeping.at(index);
    }

    //! set track candidate
    void setCandidate(const GFTrackCand& cand, bool doreset = false);

    /** @brief Choose cardinal track represenatation
   *
   * @sa getCardinalRep
   */
    void setCardinalRep(unsigned int r)
    {
      if ((int)r < fTrackReps->GetEntriesFast()) fCardinal_rep = r;
    }

    void setHitMeasuredCov(TMatrixT<Double_t> mat) { fHitMeasCov.push_back(mat); }
    void setHitUpdate(TMatrixT<Double_t> mat) { fHitUpdate.push_back(mat); }
    void setHitChi2(Double_t mat) { fHitChi2.push_back(mat); }
    void setHitState(TMatrixT<Double_t> mat) { fHitState.push_back(mat); }
    void setHitCov(TMatrixT<Double_t> mat) { fHitCov.push_back(mat); }
    void setHitCov7x7(TMatrixT<Double_t> mat) { fHitCov7x7.push_back(mat); }
    void setHitPlaneXYZ(TVector3 pl) { fHitPlaneXYZ.push_back(pl); }
    void setHitPlaneUxUyUz(TVector3 pl) { fHitPlaneUxUyUz.push_back(pl); }
    void setHitPlaneU(TVector3 pl) { fHitPlaneU.push_back(pl); }
    void setHitPlaneV(TVector3 pl) { fHitPlaneV.push_back(pl); }
    void setPDG(int pdgt) { fPDG = pdgt; }
    std::vector<TMatrixT<Double_t>> getHitMeasuredCov() { return fHitMeasCov; }
    std::vector<TMatrixT<Double_t>> getHitUpdate() { return fHitUpdate; }
    std::vector<Double_t> getHitChi2() { return fHitChi2; }
    std::vector<TMatrixT<Double_t>> getHitState() { return fHitState; }
    std::vector<TMatrixT<Double_t>> getHitCov() { return fHitCov; }
    std::vector<TMatrixT<Double_t>> getHitCov7x7() { return fHitCov; }
    std::vector<TVector3> getHitPlaneXYZ() { return fHitPlaneXYZ; }
    std::vector<TVector3> getHitPlaneUxUyUz() { return fHitPlaneUxUyUz; }
    std::vector<TVector3> getHitPlaneU() { return fHitPlaneU; }
    std::vector<TVector3> getHitPlaneV() { return fHitPlaneV; }

    int getPDG() { return fPDG; }
    /** @brief Get residuals
   *
   * @param detId which detector?
   * @param dim = index of coordinate to choose from resiudal vector
   * @param rep which track representation?
   * @param result results are written to this vector
   */
    void getResiduals(unsigned int detId, // which detector?
                      unsigned int dim,   // which projection?
                      unsigned int rep,   // which trackrep ?
                      std::vector<double>& result);

    /** @brief set the hit index at which plane,state&cov of rep irep is defined
   */
    void setRepAtHit(unsigned int irep, int ihit)
    {
      if (irep >= getNumReps())
        throw GFException("genf::GFTrack::setRepAtHit(): index out of range", __LINE__, __FILE__)
          .setFatal();
      fRepAtHit.at(irep) = ihit;
    }

    /** @brief get the hit index at which plane,state&cov of rep irep is defined
   */
    int getRepAtHit(unsigned int irep)
    {
      if (irep >= getNumReps())
        throw GFException("genf::GFTrack::getRepAtHit(): index out of range", __LINE__, __FILE__)
          .setFatal();
      return fRepAtHit.at(irep);
    }

    /** @brief clear the hit indices at which plane,state&cov of reps are defined
   */
    void clearRepAtHit()
    {
      for (unsigned int i = 0; i < getNumReps(); ++i) {
        fRepAtHit.at(i) = -1;
      }
    }

    /** @brief print bookkeeping
   */
    void printBookkeeping(std::ostream& out = std::cout) const;

    void Print(std::ostream& out = std::cout) const;

    void clearBookkeeping()
    {
      for (unsigned int i = 0; i < getNumReps(); ++i) {
        fBookkeeping.at(i)->clearAll();
      }
    }

    void clearFailedHits()
    {
      for (unsigned int i = 0; i < getNumReps(); ++i) {
        fBookkeeping.at(i)->clearFailedHits();
      }
    }

    //! use planeId information of GFTrackCand and return by ref groups of hit ids which
    //! are in the same planes.
    void getHitsByPlane(std::vector<std::vector<int>*>& retVal);

  private:
    virtual void Print(Option_t*) const
    {
      throw std::logic_error(std::string(__func__) + "::Print(Option_t*) not available");
    }

    //public:
    //ClassDef(GFTrack,1)
  };
} // namespace genf

#endif

/** @} */
