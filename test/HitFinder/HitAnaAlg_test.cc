#define BOOST_TEST_MODULE (HitAnaAlg_test)
#include "boost/test/unit_test.hpp"

#include "larreco/HitFinder/HitAnaAlg.h"

namespace hit {

  class HitAnaAlgTest {
  public:
    WireROIInfo const& GetWireDataStruct() const { return alg.wireData; }
    std::vector<std::string> const& GetHitModuleLabels() const { return alg.HitModuleLabels; }
    std::vector<HitAnaAlg::HitAssocPair> const& GetHitProcessingQueue() const
    {
      return alg.HitProcessingQueue;
    }

    void LoadHitAssocPair(std::vector<recob::Hit> const& HitVector,
                          std::vector<std::vector<int>> const& AssocVector,
                          std::string const& HitModuleLabel)
    {
      alg.LoadHitAssocPair(HitVector, AssocVector, HitModuleLabel);
    }

    void InitWireData(unsigned int e, unsigned int r) { alg.InitWireData(e, r); }

    void ClearWireDataHitInfo() { alg.ClearWireDataHitInfo(); }

    void AddHitModuleLabel(std::string str) { alg.HitModuleLabels.push_back(str); }

  private:
    HitAnaAlg alg;
  };

}

struct HitAnaAlgFixture {

  HitAnaAlgFixture() : myHitAnaAlgTest() {}
  hit::HitAnaAlgTest myHitAnaAlgTest;
};

BOOST_FIXTURE_TEST_SUITE(HitAnaAlg_test, HitAnaAlgFixture)

BOOST_AUTO_TEST_CASE(checkConstructor)
{

  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHitModules == 0);
  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels().size() == 0U);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue().size() == 0U);
}

BOOST_AUTO_TEST_CASE(LoadHitAssocPair_FirstTime)
{
  size_t nHits = 10;
  size_t nWires = 10;

  std::vector<recob::Hit> HitVector(nHits);

  std::vector<std::vector<int>> AssocVector(nWires);
  for (size_t iter = 0; iter < nWires; iter++)
    AssocVector[iter].push_back(iter);

  std::string HitModuleLabel = "hit";

  myHitAnaAlgTest.LoadHitAssocPair(HitVector, AssocVector, HitModuleLabel);
  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels().size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue().size() == 1U);
}

BOOST_AUTO_TEST_CASE(LoadHitAssocPair_MultipleHitModules)
{
  size_t nHits1 = 10;
  size_t nHits2 = 10;
  size_t nWires = 10;

  std::vector<recob::Hit> HitVector1(nHits1);
  std::vector<recob::Hit> HitVector2(nHits2);

  std::vector<std::vector<int>> AssocVector1(nWires);
  std::vector<std::vector<int>> AssocVector2(nWires);
  for (size_t iter = 0; iter < nWires; iter++) {
    AssocVector1[iter].push_back(iter);
    AssocVector2[iter].push_back(nWires - iter - 1);
  }
  std::string HitModuleLabel1 = "hit1";
  std::string HitModuleLabel2 = "hit2";

  myHitAnaAlgTest.LoadHitAssocPair(HitVector1, AssocVector1, HitModuleLabel1);
  myHitAnaAlgTest.LoadHitAssocPair(HitVector2, AssocVector2, HitModuleLabel2);

  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels().size() == 2U);
  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels()[0].compare(HitModuleLabel1) == 0);
  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels()[1].compare(HitModuleLabel2) == 0);

  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue().size() == 2U);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[0].second.size() == nWires);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[0].second[0].size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[0].second[0][0] == 0);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[1].second.size() == nWires);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[1].second[0].size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetHitProcessingQueue()[1].second[0][0] == (int)(nWires - 1));
}

BOOST_AUTO_TEST_CASE(LoadHitAssocPair_VectorSizesOff)
{
  std::string str = "test";
  myHitAnaAlgTest.AddHitModuleLabel(str);
  BOOST_TEST(myHitAnaAlgTest.GetHitModuleLabels().size() == 1U);

  size_t nHits = 10;
  size_t nWires = 10;

  std::vector<recob::Hit> HitVector(nHits);

  std::vector<std::vector<int>> AssocVector(nWires);
  for (size_t iter = 0; iter < nWires; iter++)
    AssocVector[iter].push_back(iter);

  std::string HitModuleLabel = "hit";

  BOOST_CHECK_THROW(myHitAnaAlgTest.LoadHitAssocPair(HitVector, AssocVector, HitModuleLabel),
                    hit::HitAnaAlgException);
}

BOOST_AUTO_TEST_CASE(InitWireData_OneModule)
{
  size_t nHits = 10;
  size_t nWires = 10;

  std::vector<recob::Hit> HitVector(nHits);

  std::vector<std::vector<int>> AssocVector(nWires);
  for (size_t iter = 0; iter < nWires; iter++)
    AssocVector[iter].push_back(iter);

  std::string HitModuleLabel = "hit";

  myHitAnaAlgTest.LoadHitAssocPair(HitVector, AssocVector, HitModuleLabel);

  unsigned int event = 50;
  unsigned int run = 200;

  myHitAnaAlgTest.InitWireData(event, run);

  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().event == event);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().run == run);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHitModules == 1);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().HitModuleLabels.size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().HitModuleLabels[0] == HitModuleLabel);
}

BOOST_AUTO_TEST_CASE(InitWireData_NoModules)
{
  unsigned int event = 50;
  unsigned int run = 200;

  myHitAnaAlgTest.InitWireData(event, run);

  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().event == event);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().run == run);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHitModules == 0);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().HitModuleLabels.size() == 0U);
}

BOOST_AUTO_TEST_CASE(ClearWireDataHitInfo_NoModules)
{
  unsigned int event = 50;
  unsigned int run = 200;

  myHitAnaAlgTest.InitWireData(event, run);
  myHitAnaAlgTest.ClearWireDataHitInfo();

  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().event == event);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().run == run);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHitModules == 0);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHits.size() == 0U);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().Hits_IntegratedCharge.size() == 0U);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().Hits.size() == 0U);
}

BOOST_AUTO_TEST_CASE(ClearWireDataHitInfo_OneModule)
{
  size_t nHits = 10;
  size_t nWires = 10;

  std::vector<recob::Hit> HitVector(nHits);

  std::vector<std::vector<int>> AssocVector(nWires);
  for (size_t iter = 0; iter < nWires; iter++)
    AssocVector[iter].push_back(iter);

  std::string HitModuleLabel = "hit";

  myHitAnaAlgTest.LoadHitAssocPair(HitVector, AssocVector, HitModuleLabel);

  unsigned int event = 50;
  unsigned int run = 200;

  myHitAnaAlgTest.InitWireData(event, run);
  myHitAnaAlgTest.ClearWireDataHitInfo();

  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().event == event);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().run == run);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHitModules == 1);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHits.size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().NHits[0] == 0);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().Hits_IntegratedCharge.size() == 1U);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().Hits_IntegratedCharge[0] == 0);
  BOOST_TEST(myHitAnaAlgTest.GetWireDataStruct().Hits.size() == 1U);
}

BOOST_AUTO_TEST_SUITE_END()
