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
 * P1Element.cpp
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/P1Elements.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        double P1Element::m_tolerance = 1e-6;
        /*--------------------------------------------------------------------*/
        void P1Element::setTolerance(const double ATolerance)
        {
            m_tolerance=ATolerance;
        }
        /*--------------------------------------------------------------------*/
        double P1Element::getTolerance()
        {
            return m_tolerance;
        }
        /*--------------------------------------------------------------------*/
        P1Element::P1Element(const std::vector<Point>& APnts)
        {
            auto nb_points = APnts.size();
            m_x = new double[nb_points];
            m_y = new double[nb_points];
            m_z = new double[nb_points];
            for(auto i = 0; i < nb_points; i++){
                math::Point pi = APnts[i];
                m_x[i] = pi.X();
                m_y[i] = pi.Y();
                m_z[i] = pi.Z();
            }
   
        }
        /*--------------------------------------------------------------------*/
        P1Element::~P1Element()
        {
            if(m_x!=NULL)
                delete [] m_x;
            if(m_y!=NULL)
                delete [] m_y;
            if(m_z!=NULL)
                delete [] m_z;

        }
        /*--------------------------------------------------------------------*/
        Point P1Element::getPoint(int AIndex) const
        {
            Point p(0,0,0);
            int nb_points = 0;//getNbPoints();
            if(AIndex >=0 && AIndex <= nb_points) {
                p = math::Point(m_x[AIndex], m_y[AIndex], m_z[AIndex]);
            }
            
            return p;
        }

    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
