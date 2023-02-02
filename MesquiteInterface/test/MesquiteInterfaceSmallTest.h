/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * MesquiteInterfaceSmallTest.h
 *
 *  Created on: 30 oct 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/MesquiteCaller.h>
#include <GMDS/IG/IG.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Utils/RandomGenerator.h>
/*----------------------------------------------------------------------------*/
class MesquiteInterfaceSmallTest: public ::testing::Test {

  protected:
	MesquiteInterfaceSmallTest(){;}
    virtual ~MesquiteInterfaceSmallTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(MesquiteInterfaceSmallTest,DISABLED_cas_quad2D) {

	gmds::MeshModel mod = gmds::DIM2|gmds::N|gmds::F|gmds::F2N;
        gmds::IGMesh mesh(mod);

        const int nx = 3;
        const int ny = 3;

        gmds::Node nodes[nx+1][ny+1];

        const double xmin = 0.;
        const double ymin = 0.;

        const double dx = 1.;
        const double dy = 1.;

        int markFixedNodes = mesh.getNewMark<gmds::Node>();

        gmds::RandomGenerator randGen;
        randGen.init();

        for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {

                                double xpos = xmin+i*dx;
                                double ypos = ymin+j*dy;

                                if(i!= 0 && j!=0 && i!=nx && j!=ny) {

                                        double xrand = (2.*randGen.value()) -1.;
                                        double yrand = (2.*randGen.value()) -1.;

                                        xpos += 3*xrand;
                                        ypos += 3*yrand;
                                }

                                nodes[i][j] = mesh.newNode(xpos,ypos);

                                if(i== 0 || j==0 || i==nx || j==ny) {
                                        mesh.mark(nodes[i][j],markFixedNodes);
                                }
                        }
        }


        gmds::Face faces[nx][ny];

	for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        gmds::Face f = mesh.newQuad(
                                                    nodes[i][j]    ,
                                                nodes[i+1][j]     ,
                                                nodes[i+1][j+1]     ,
                                                nodes[i][j+1]
                                       );
                        }
                }

        gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("MesquiteInterfaceSmallTest_quad2D_before",gmds::N|gmds::F);
	//double minScaledJacobian = getMe.computeMinScaledJacobian2D();
	double minScaledJacobian = 0;

	gmds::geom::FacetedGeomManager geomManager;

	gmds::MesquiteCaller mesquiteCaller(mesh,geomManager);
	mesquiteCaller.execSurf(markFixedNodes);

	writer.write("MesquiteInterfaceSmallTest_quad2D_after",gmds::N|gmds::F);

        //minScaledJacobian = getMe.computeMinScaledJacobian2D();
        std::cout<<"minScaledJacobian "<<minScaledJacobian<<std::endl;

        mesh.unmarkAll<gmds::Node>(markFixedNodes);
        mesh.freeMark<gmds::Node>(markFixedNodes);

        EXPECT_LT(0.8,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
