/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
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
 * Patching2D.cpp
 *
 *  Created on: Sept 10, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/Math/Line.h>
/*----------------------------------------------------------------------------*/
#include "Patch2D.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Patch2D::Patch2D(PatchComplex2D* AOwner,std::vector<PatchEdge2D*>& AEdges):
m_complex(AOwner),
m_X(10),m_Y(10),
m_has_no_grid(true),
m_edges(AEdges)
{;}
/*----------------------------------------------------------------------------*/
Patch2D::Patch2D()
:  m_complex(0),
m_X(10),m_Y(10),
m_has_no_grid(true)
{;}

/*----------------------------------------------------------------------------*/
Patch2D::Patch2D(const Patch2D& AP):
m_complex(AP.m_complex),
m_X(AP.m_X),
m_Y(AP.m_Y),
m_grid(AP.m_grid),
m_has_no_grid(AP.m_has_no_grid),
m_edges(AP.m_edges)
{;}
/*----------------------------------------------------------------------------*/
Patch2D::~Patch2D()
{;}
/*----------------------------------------------------------------------------*/
int Patch2D::X() const
{
    return m_X;
}
/*----------------------------------------------------------------------------*/
int Patch2D::Y() const
{
    return m_Y;
}
/*----------------------------------------------------------------------------*/
Patch2D::Grid& Patch2D::grid()
{
    return m_grid;
}

/*----------------------------------------------------------------------------*/
bool Patch2D::buildInnerGrid()
{
    
    double stepX = 1.0/(m_X-1);
    double stepY = 1.0/(m_Y-1);

    //============================================
    // Build inner nodes
    //============================================
    for(int i=1; i<m_X-1;i++){
        for(int j=1; j<m_Y-1; j++){
            double alphaX = i*stepX;
            double alphaY = j*stepY;
            math::Point pi = (1-alphaY)*m_grid.points[i][0]+alphaY*m_grid.points[i][m_Y-1];
            math::Point pj = (1-alphaX)*m_grid.points[0][j]+alphaX*m_grid.points[m_X-1][j];
            m_grid.points[i][j] = 0.5*(pi+pj);
        }
    }
    return true;
}
/*----------------------------------------------------------------------------*/
bool Patch2D::buildBoundaryGrid()
{
    
    //initialization of the grid structure
    m_grid =Grid(m_X,m_Y);
    
    double stepX = 1.0/(m_X-1);
    double stepY = 1.0/(m_Y-1);
    
//    math::Point c0 = m_corners[0];
//    math::Point c1 = m_corners[1];
//    math::Point c2 = m_corners[2];
//    math::Point c3 = m_corners[3];
//    
//    
//    //============================================
//    // Build the corner points
//    //============================================
//    m_grid.points[0    ][0    ] = c0;
//    m_grid.points[m_X-1][0    ] = c1;
//    m_grid.points[m_X-1][m_Y-1] = c2;
//    m_grid.points[0    ][m_Y-1] = c3;
//    
//    //============================================
//    // Build 2 extrema lines in constant X
//    //============================================
//    for(int i=1; i<m_X-1;i++){
//        double alpha = i*stepX;
//        m_grid.points[i][0    ] = (1-alpha)*m_grid.points[0][0    ] + alpha*m_grid.points[m_X-1][0    ];
//        m_grid.points[i][m_Y-1] = (1-alpha)*m_grid.points[0][m_Y-1] + alpha*m_grid.points[m_X-1][m_Y-1];
//    }
//    //============================================
//    // Build 2 extrema lines in constant Y
//    //============================================
//    for(int i=1; i<m_Y-1;i++){
//        double alpha = i*stepY;
//        m_grid.points[0    ][i] = (1-alpha)*m_grid.points[0    ][0] + alpha*m_grid.points[0    ][m_Y-1];
//        m_grid.points[m_X-1][i] = (1-alpha)*m_grid.points[m_X-1][0] + alpha*m_grid.points[m_X-1][m_Y-1];
//    }
//    
    m_has_no_grid=false;
    return true;
}
/*----------------------------------------------------------------------------*/
Patch2D::Grid::Grid():X(0),Y(0)
{;}
/*----------------------------------------------------------------------------*/
Patch2D::Grid::Grid(const int AX, const int AY){
    X=AX;
    Y=AY;
    points.resize(X);
    for(unsigned int i=0;i<X;i++)
        points[i].resize(Y);
}

/*----------------------------------------------------------------------------*/
Patch2D::Grid::Grid(const Grid& AG){
    X= AG.X;
    Y=AG.Y;
    points.resize(X);
    for(unsigned int i=0;i<X;i++){
        points[i].resize(Y);
    }
    for(unsigned int i=0;i<X;i++){
        for(unsigned int j=0;j<Y;j++){
            
            points[i][j] = AG.points[i][j];
        }
    }
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

