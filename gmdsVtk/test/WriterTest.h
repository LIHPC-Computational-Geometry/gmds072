/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/IOVTK/VTKReader.h>
#include <GMDS/IOVTK/VTKWriter.h>
#include <GMDS/IO/MatrixMarketWriter.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class WriterTest: public ::testing::Test {

  protected:
	WriterTest(){;}
    virtual ~WriterTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest ,vtkWrite2Triangles3D)
{
	MeshModel mod = DIM3|N|F|F2N;
	IGMesh mesh(mod);

	Node t[5];
	Face f[2];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);

	t[3] = mesh.newNode(0.5,0.5,-1.);
	t[4] = mesh.newNode(0.5,0.5,1.);

	f[0] = mesh.newTriangle(t[0],t[1],t[4]);
	f[1] = mesh.newTriangle(t[2],t[3],t[4]);

	EXPECT_EQ(2,mesh.getNbFaces());

	VTKWriter<IGMesh> writer(mesh);
	writer.write("Samples/OUT/vtkWrite2Triangles3D",N|F);
}
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest ,vtkWrite3D)
{
	EXPECT_NO_THROW({
	MeshModel mod = DIM3|N|R|R2N;
	IGMesh mesh(mod);


	Node n0 = mesh.newNode(0,0,0);
	Node n1 = mesh.newNode(0,1,0);
	Node n2 = mesh.newNode(1,1,0);
	Node n3 = mesh.newNode(1,0,0);
	Node n4 = mesh.newNode(0,0,1);
	Node n5 = mesh.newNode(0,1,1);
	Node n6 = mesh.newNode(1,1,1);
	Node n7 = mesh.newNode(1,0,1);

	Node n8 = mesh.newNode(3,0,0);
	Node n9 = mesh.newNode(3,1,0);
	Node n10= mesh.newNode(4,1,0);
	Node n11= mesh.newNode(3,1,1);

	Region r1 = mesh.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
	Region r2 = mesh.newTet(n8,n9,n10,n11);

	EXPECT_EQ(2,mesh.getNbRegions());

	VTKWriter<IGMesh> writer(mesh);
	writer.write("Samples/OUT/vtkWrite3D",N|R);

	VTKReader<IGMesh> reader(mesh);
	reader.read("Samples/OUT/vtkWrite3D.vtu");

	EXPECT_EQ(2,mesh.getNbRegions());
	EXPECT_EQ(12,mesh.getNbNodes());
	writer.write("Samples/OUT/vtkWrite3D_2",N|R);
	}
	);
}

/*----------------------------------------------------------------------------*/
TEST_F(WriterTest ,vtkWriteGrid2D)
{
	EXPECT_NO_THROW({
		MeshModel mod = DIM2|N|F|F2N;
		IGMesh mesh(mod);
		const TInt nb=10;
		Node t[nb][nb];

		for(unsigned int i=0;i<nb;i++)
			for(unsigned int j=0;j<nb;j++)
				t[i][j]= mesh.newNode(i,j,0);

		for(unsigned int i=0;i<nb-1;i++)
			for(unsigned int j=0;j<nb-1;j++)
				mesh.newQuad(t[i][j],t[i+1][j],t[i+1][j+1],t[i][j+1]);

		VTKWriter<IGMesh> writer(mesh);
		writer.write("Samples/OUT/vtkWriteGrid2D",N|F);
	}
	);
}
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest ,vtkWriteVariables)
{
	EXPECT_NO_THROW({
		MeshModel mod = DIM2|N|F|F2N;
		IGMesh mesh(mod);
		Variable<TInt>*   v1 = mesh.newVariable<TInt>(GMDS_NODE,"variable entiere");
		Variable<double>* v2 = mesh.newVariable<double>(GMDS_NODE,"variable flottante");
		Variable<double>* v3 = mesh.newVariable<double>(GMDS_FACE,"variable flottante");
		Variable<math::Vector>* v4 = mesh.newVariable<math::Vector>(GMDS_NODE,"variable vectorielle");
		const TInt nb=10;
		Node t[nb][nb];

		for(unsigned int i=0;i<nb;i++)
			for(unsigned int j=0;j<nb;j++){
				t[i][j]= mesh.newNode(i,j,0);
				(*v1)[t[i][j].getID()]=i;
				(*v2)[t[i][j].getID()]=(double)(i+j)*2.22;
				if(i>j){
					(*v4)[t[i][j].getID()][0]=1;
					(*v4)[t[i][j].getID()][1]=0;
					(*v4)[t[i][j].getID()][2]=0;
				}
				else {
					(*v4)[t[i][j].getID()][0]=1;
					(*v4)[t[i][j].getID()][1]=1;
					(*v4)[t[i][j].getID()][2]=1;
				}
			}

		for(unsigned int i=0;i<nb-1;i++)
			for(unsigned int j=0;j<nb-1;j++){
				Face f= mesh.newQuad(t[i][j],t[i+1][j],t[i+1][j+1],t[i][j+1]);
				if(i>j)
					(*v3)[f.getID()]=(double)(i+0.5);
			}

		VTKWriter<IGMesh> writer(mesh);
		writer.write("Samples/OUT/vtkWrite_variables",N|F);
	}
	);
}
/*----------------------------------------------------------------------------*/
