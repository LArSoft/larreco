////////////////////////////////////////////////////////////////////////////
//
// \brief Definition of 3D cluster object for LArSoft
//
// \author usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

namespace reco {

  ClusterHit2D::ClusterHit2D()
    : m_statusBits(0)
    , m_docaToAxis(9999.)
    , m_arcLenToPoca(0.)
    , m_xPosition(0.)
    , m_timeTicks(0.)
    , m_wireID(geo::WireID())
    , m_hit(nullptr)
  {}

  ClusterHit2D::ClusterHit2D(unsigned statusBits,
                             float doca,
                             float poca,
                             float xPosition,
                             float timeTicks,
                             const geo::WireID& wireID,
                             const recob::Hit* hit)
    : m_statusBits(statusBits)
    , m_docaToAxis(doca)
    , m_arcLenToPoca(poca)
    , m_xPosition(xPosition)
    , m_timeTicks(timeTicks)
    , m_wireID(wireID)
    , m_hit(hit)
  {}

  ClusterHit2D::ClusterHit2D(const ClusterHit2D& toCopy)
  {
    m_statusBits = toCopy.m_statusBits;
    m_docaToAxis = toCopy.m_docaToAxis;
    m_arcLenToPoca = toCopy.m_arcLenToPoca;
    m_xPosition = toCopy.m_xPosition;
    m_timeTicks = toCopy.m_timeTicks;
    m_wireID = toCopy.m_wireID;
    m_hit = toCopy.m_hit;
  }

  std::ostream& operator<<(std::ostream& o, const ClusterHit2D& c)
  {
    o << c.getHit();

    return o;
  }

  bool operator<(const ClusterHit2D& a, const ClusterHit2D& b) { return a.getHit() < b.getHit(); }

  ClusterHit3D::ClusterHit3D()
    : fID(std::numeric_limits<size_t>::max())
    , fStatusBits(0)
    , fPosition(Eigen::Vector3f::Zero())
    , fTotalCharge(0.)
    , fAvePeakTime(-1.)
    , fDeltaPeakTime(0.)
    , fSigmaPeakTime(0.)
    , fHitChiSquare(0.)
    , fOverlapFraction(0.)
    , fChargeAsymmetry(0.)
    , fDocaToAxis(0.)
    , fArclenToPoca(0.)
  {
    fHitDelTSigVec.clear();
    fWireIDVector.clear();
    fHitVector.clear();
    fHitDelTSigVec.resize(3, 0.);
    fWireIDVector.resize(3, geo::WireID());
    fHitVector.resize(3, NULL);
  }

  ClusterHit3D::ClusterHit3D(size_t id,
                             unsigned int statusBits,
                             const Eigen::Vector3f& position,
                             float totalCharge,
                             float avePeakTime,
                             float deltaPeakTime,
                             float sigmaPeakTime,
                             float hitChiSquare,
                             float overlapFraction,
                             float chargeAsymmetry,
                             float docaToAxis,
                             float arclenToPoca,
                             const ClusterHit2DVec& hitVec,
                             const std::vector<float>& hitDelTSigVec,
                             const std::vector<geo::WireID>& wireIDs)
    : fID(id)
    , fStatusBits(statusBits)
    , fPosition(position)
    , fTotalCharge(totalCharge)
    , fAvePeakTime(avePeakTime)
    , fDeltaPeakTime(deltaPeakTime)
    , fSigmaPeakTime(sigmaPeakTime)
    , fHitChiSquare(hitChiSquare)
    , fOverlapFraction(overlapFraction)
    , fChargeAsymmetry(chargeAsymmetry)
    , fDocaToAxis(docaToAxis)
    , fArclenToPoca(arclenToPoca)
    , fHitDelTSigVec(hitDelTSigVec)
    , fWireIDVector(wireIDs)
  {
    fHitVector.resize(3, NULL);
    std::copy(hitVec.begin(), hitVec.end(), fHitVector.begin());
  }

  ClusterHit3D::ClusterHit3D(const ClusterHit3D& toCopy)
  {
    fID = toCopy.fID;
    fStatusBits = toCopy.fStatusBits;
    fPosition = toCopy.fPosition;
    fTotalCharge = toCopy.fTotalCharge;
    fAvePeakTime = toCopy.fAvePeakTime;
    fDeltaPeakTime = toCopy.fDeltaPeakTime;
    fSigmaPeakTime = toCopy.fSigmaPeakTime;
    fHitChiSquare = toCopy.fHitChiSquare;
    fOverlapFraction = toCopy.fOverlapFraction;
    fChargeAsymmetry = toCopy.fChargeAsymmetry;
    fDocaToAxis = toCopy.fDocaToAxis;
    fArclenToPoca = toCopy.fArclenToPoca;
    fHitVector = toCopy.fHitVector;
    fHitDelTSigVec = toCopy.fHitDelTSigVec;
    fWireIDVector = toCopy.fWireIDVector;
  }

  void ClusterHit3D::initialize(size_t id,
                                unsigned int statusBits,
                                const Eigen::Vector3f& position,
                                float totalCharge,
                                float avePeakTime,
                                float deltaPeakTime,
                                float sigmaPeakTime,
                                float hitChiSquare,
                                float overlapFraction,
                                float chargeAsymmetry,
                                float docaToAxis,
                                float arclenToPoca,
                                const ClusterHit2DVec& hitVec,
                                const std::vector<float>& hitDelTSigVec,
                                const std::vector<geo::WireID>& wireIDs)
  {
    fID = id;
    fStatusBits = statusBits;
    fPosition = position;
    fTotalCharge = totalCharge;
    fAvePeakTime = avePeakTime;
    fDeltaPeakTime = deltaPeakTime;
    fSigmaPeakTime = sigmaPeakTime;
    fHitChiSquare = hitChiSquare;
    fOverlapFraction = overlapFraction;
    fChargeAsymmetry = chargeAsymmetry;
    fDocaToAxis = docaToAxis;
    fArclenToPoca = arclenToPoca;
    fHitVector = hitVec;
    fHitDelTSigVec = hitDelTSigVec;
    fWireIDVector = wireIDs;

    return;
  }

  void ClusterHit3D::setWireID(const geo::WireID& wid) const { fWireIDVector[wid.Plane] = wid; }

  std::ostream& operator<<(std::ostream& o, const ClusterHit3D& c)
  {
    o << "ClusterHit3D has " << c.getHits().size() << " hits associated";

    return o;
  }

  //bool operator < (const ClusterHit3D & a, const ClusterHit3D & b)
  //{
  //    if (a.m_position[2] != b.m_position[2]) return a.m_position[2] < b.m_position[2];
  //    else return a.m_position[0] < b.m_position[0];
  //}

  PrincipalComponents::PrincipalComponents()
    : m_svdOK(false)
    , m_numHitsUsed(0)
    , m_eigenValues(EigenValues::Zero())
    , m_eigenVectors(EigenVectors::Zero())
    , m_avePosition(Eigen::Vector3f::Zero())
    , m_aveHitDoca(9999.)
  {}

  PrincipalComponents::PrincipalComponents(bool ok,
                                           int nHits,
                                           const EigenValues& eigenValues,
                                           const EigenVectors& eigenVecs,
                                           const Eigen::Vector3f& avePos,
                                           const float aveHitDoca)
    : m_svdOK(ok)
    , m_numHitsUsed(nHits)
    , m_eigenValues(eigenValues)
    , m_eigenVectors(eigenVecs)
    , m_avePosition(avePos)
    , m_aveHitDoca(aveHitDoca)
  {}

  void PrincipalComponents::flipAxis(size_t axisDir)
  {
    m_eigenVectors.row(axisDir) = -m_eigenVectors.row(axisDir);

    return;
  }

  std::ostream& operator<<(std::ostream& o, const PrincipalComponents& a)
  {
    if (a.m_svdOK) {
      o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
      o << " PCAxis ID run with " << a.m_numHitsUsed << " space points" << std::endl;
      o << "   - center position: " << std::setw(6) << a.m_avePosition(0) << ", "
        << a.m_avePosition(1) << ", " << a.m_avePosition(2) << std::endl;
      o << "   - eigen values: " << std::setw(8) << std::right << a.m_eigenValues(0) << ", "
        << a.m_eigenValues(1) << ", " << a.m_eigenValues(1) << std::endl;
      o << "   - average doca: " << a.m_aveHitDoca << std::endl;
      o << "   - Principle axis: " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors(0, 0)
        << ", " << a.m_eigenVectors(0, 1) << ", " << a.m_eigenVectors(0, 2) << std::endl;
      o << "   - second axis:    " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors(1, 0)
        << ", " << a.m_eigenVectors(1, 1) << ", " << a.m_eigenVectors(1, 2) << std::endl;
      o << "   - third axis:     " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors(2, 0)
        << ", " << a.m_eigenVectors(2, 1) << ", " << a.m_eigenVectors(2, 2) << std::endl;
    }
    else
      o << " Principal Components Axis is not valid" << std::endl;

    return o;
  }

  bool operator<(const PrincipalComponents& a, const PrincipalComponents& b)
  {
    if (a.m_svdOK && b.m_svdOK) return a.m_eigenValues(0) > b.m_eigenValues(0);

    return false; //They are equal
  }

  Cluster3D::Cluster3D()
    : m_statusBits(0)
    , m_pcaResults(PrincipalComponents())
    , m_totalCharge(0.)
    , m_startPosition{0., 0., 0.}
    , m_endPosition{0., 0., 0.}
    , m_clusterIdx(0)
  {}

  Cluster3D::Cluster3D(unsigned statusBits,
                       const PrincipalComponents& pcaResults,
                       float totalCharge,
                       const float* startPosition,
                       const float* endPosition,
                       int idx)
    : m_statusBits(statusBits)
    , m_pcaResults(pcaResults)
    , m_totalCharge(totalCharge)
    , m_startPosition{startPosition[0], startPosition[1], startPosition[2]}
    , m_endPosition{endPosition[0], endPosition[1], endPosition[2]}
    , m_clusterIdx(idx)
  {}

  //----------------------------------------------------------------------
  //  Addition operator.
  //
  Cluster3D Cluster3D::operator+(Cluster3D a)
  {
    /*
    // throw exception if the clusters are not from the same plane
    if( a.View() != this->View() )
      throw cet::exception("Cluster+operator") << "Attempting to sum clusters from "
                 << "different views is not allowed\n";

    // check the start and end positions - for now the
    // smallest wire number means start position, largest means end position
    std::vector<float> astart(a.StartPos());
    std::vector<float> aend  (a.EndPos()  );
    std::vector<float> start(StartPos());
    std::vector<float> end  (EndPos()  );
    std::vector<float> sigstart(SigmaStartPos());
    std::vector<float> sigend  (SigmaEndPos()  );

    if(astart[0] < fStartPos[0]){
      start = astart;
      sigstart = a.SigmaStartPos();
    }

    if(aend[0] > fEndPos[0]){
      end = aend;
      sigend = a.SigmaEndPos();
    }

    //take weighted mean in obtaining average slope and differential charge,
    //based on total charge each cluster
    float dtdw = ((this->Charge()*dTdW()) + (a.Charge()*a.dTdW()))/(this->Charge() + a.Charge());
    float dqdw = ((this->Charge()*dQdW()) + (a.Charge()*a.dQdW()))/(this->Charge() + a.Charge());

    //hits.sort();//sort the PtrVector to organize Hits of new Cluster
    float sigdtdw = TMath::Max(SigmadTdW(), a.SigmadTdW());
    float sigdqdw = TMath::Max(SigmadQdW(), a.SigmadQdW());

    Cluster sum(//hits,
    start[0], sigstart[0],
    start[1], sigstart[1],
    end[0],   sigend[0],
    end[1],   sigend[1],
    dtdw, sigdtdw,
    dqdw, sigdqdw,
    this->Charge() + a.Charge(),
    this->View(),
    ID());
*/
    //return sum;
    return a;
  }

  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator<<(std::ostream& o, const Cluster3D& c)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << "Cluster ID " << std::setw(5) << std::right << c.getClusterIdx();
    //      << " : View = "     << std::setw(3)  << std::right << c.View()
    //      << " StartWire = "  << std::setw(7)  << std::right << c.StartPos()[0]
    //      << " EndWire = "    << std::setw(7)  << std::right << c.EndPos()[0]
    //      << " StartTime = "  << std::setw(9)  << std::right << c.StartPos()[1]
    //      << " EndTime = "    << std::setw(9)  << std::right << c.EndPos()[1]
    //      << " dTdW = "       << std::setw(9)  << std::right << c.dTdW()
    //      << " dQdW = "       << std::setw(9)  << std::right << c.dQdW()
    //      << " Charge = "     << std::setw(10) << std::right << c.Charge();

    return o;
  }

  //----------------------------------------------------------------------
  // < operator.
  //
  bool operator<(const Cluster3D& a, const Cluster3D& b)
  {
    /*
    if(a.View() != b.View())
      return a.View() < b.View();
    if(a.ID() != b. ID())
      return a.ID() < b.ID();
    if(a.StartPos()[0] != b.StartPos()[0])
      return a.StartPos()[0] < b.StartPos()[0];
    if(a.EndPos()[0] != b.EndPos()[0])
      return a.EndPos()[0] < b.EndPos()[0];
*/
    if (a.getStartPosition()[2] < b.getStartPosition()[2]) return true;

    return false; //They are equal
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void RecobClusterParameters::UpdateParameters(const reco::ClusterHit2D* clusterHit)
  {
    /**
     *  @brief a utility routine for building 3D clusters to keep basic info up to date
     *         (a candidate for a better way to do this)
     */
    const recob::Hit* hit = clusterHit->getHit();

    // Need to keep track of stuff so we can form cluster
    if (clusterHit->WireID().Wire < m_startWire) {
      m_startWire = clusterHit->WireID().Wire;
      m_startTime = hit->PeakTimeMinusRMS();
      m_sigmaStartTime = hit->SigmaPeakTime();
    }

    if (clusterHit->WireID().Wire > m_endWire) {
      m_endWire = clusterHit->WireID().Wire;
      m_endTime = hit->PeakTimePlusRMS();
      m_sigmaEndTime = hit->SigmaPeakTime();
    }

    m_totalCharge += hit->Integral();
    m_plane = clusterHit->WireID().planeID();
    m_view = hit->View();

    m_hitVector.push_back(clusterHit);

    return;
  }

} // namespace
