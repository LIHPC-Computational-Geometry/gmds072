/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
//#include "FacetedGeometryTest.h"
//#include "GeometryInsertionTest.h"
//#include "GeometryInsertionSplitterTest.h"
//#include "GeometryInsertionSheetTest.h"
//#include "GeometryInsertionMeshTrends.h"
//#include "GeometryInsertionMeshTrends2.h"
#include "GeometryInsertionSumUp.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

