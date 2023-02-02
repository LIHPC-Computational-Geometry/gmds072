/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "MathTest.h"
#include "AdjacencyTest.h"
#include "BooleanMarkTest.h"
#include "CellCreationTest.h"
#include "ModelChangeTest.h"
#include "VariableTest.h"
#include "MeshVariableTest.h"  
#include "WriterTest.h"
#include "ReaderTest.h"
#include "GroupTest.h"
#include "IGMeshDoctorTest.h"
#include "QualityTest.h"
#include "SheetOperatorTest.h"
#include "SerializationTest.h"
#include "Cross2DTest.h"
#include "CrossTest.h"
#include "ChartTest.h"
#include "QuaternionTest.h"
#include "GETMeTest.h"
#include "FacetedModelTest.h"
#include "PillowingTest.h"
#include "HexahedronTest.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

