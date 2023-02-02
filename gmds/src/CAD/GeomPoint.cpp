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
/** \file    GeomPoint.cpp
 *  \author  N. LE GOFF
 *  \date    03/07/2013
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/GeomPoint.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
void GeomPoint::getOrdered(std::vector<GeomCurve*>& ACur) const
{
	// we empty the result vector
	ACur.empty();

	// we get the curves adjacent to this point
	std::vector<GeomCurve*> curves;
	this->get(curves);

	// this is a check to ensure that there are no two following curves
	// that have only one same surface;
	// in case there is we do not yet know what to do.
	{
		std::vector<unsigned int> nbOfSurfaces(curves.size(),0);
		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			std::vector<GeomSurface* > surfaces;
			curves[iCurve]->get(surfaces);

			nbOfSurfaces[iCurve] = surfaces.size();
		}

		// if two following curves have both only one surface, it is necessarily
		// the same surface hence it is not allowed at this moment
		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			if(nbOfSurfaces[iCurve] == nbOfSurfaces[(iCurve+1)%nbOfSurfaces.size()]) {
				GMDSException("GeomPoint::getOrdered this call is not allowed because "
						"two following curves both have only one same surface");
			}
		}
	}

	// we duplicate the loop curves in order to make them appear twice in the vector.
	std::vector<GeomCurve*> curves_tmp;
	for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
		curves_tmp.push_back(curves[iCurve]);
		if(curves[iCurve]->isALoop()) {
			curves_tmp.push_back(curves[iCurve]);
		}
	}
	curves = curves_tmp;

	std::vector<bool> curves_treated(curves.size(),false);

	// starting from the first curve, we jump from surface to surface
	GeomCurve* current_curve = curves[0];
	ACur.push_back(current_curve);
	curves_treated[0] = true;

	GeomSurface* current_surface = NULL;

	bool loopDone = false;

	while (!loopDone)
	{
		GeomCurve* next_curve = NULL;

		std::vector<GeomSurface*> surfaces;
		current_curve->get(surfaces);

		// this is the case for the first curve
		// we start from one of the surface(s), we do not care which one
		if (current_surface == NULL)
		{
			current_surface = surfaces[0];
		}
		else
		{
			// only one surface means that one of the end-points has only one curve
			// or the case of a periodic surface (see cylinder)
			if(1 == surfaces.size()) {
				current_surface = surfaces[0];
			} else {
				// loops and normal two end-points curves
				if(current_surface == surfaces[0])
				{
					current_surface = surfaces[1];
				}
				else
				{
					current_surface = surfaces[0];
				}
			}
		}

		std::vector<GeomCurve*> curves_pool;
		current_surface->get(curves_pool);

		// find the next curve
		for (size_t iCurve_pool=0; iCurve_pool<curves_pool.size(); iCurve_pool++)
		{
            // if there is a curve in the surface that is adjacent only to
			// this surface then it is selected in priority.
			// TODO : case when there are two or more such curves on the surface,
			// we need to determine which one will be taken first.
            for (size_t iCurve=1; iCurve<curves.size(); iCurve++)
            {
                if ((curves_pool[iCurve_pool] == curves[iCurve]) && (!curves_treated[iCurve]))
                {
                	std::vector<GeomSurface*> surfaces_pool;
                	curves[iCurve]->get(surfaces_pool);
                	if (1 == surfaces_pool.size()) {

                		curves_treated[iCurve] = true;
                		next_curve = curves[iCurve];
                		break;
                	}
                }
            }
		}
		if (next_curve == NULL) {
			for (size_t iCurve_pool=0; iCurve_pool<curves_pool.size(); iCurve_pool++)
			{
            	// we exclude the first curve that was the initial curve
            	for (size_t iCurve=1; iCurve<curves.size(); iCurve++)
            	{
            		if ((curves_pool[iCurve_pool] == curves[iCurve]) && (!curves_treated[iCurve]))
            		{
            			curves_treated[iCurve] = true;
            			next_curve = curves[iCurve];
            			break;
            		}
            	}
            }
        }
//			if (next_curve != NULL)
//			{
//				break;
//			}


		if (next_curve != NULL)
		{
			current_curve = next_curve;
			if(ACur.back() != current_curve) {
				ACur.push_back(current_curve);
			}
//			std::cout<<"selected curve"<<std::endl;
//			std::cout<<current_curve->getFirstPoint()->getPoint()<<" "<<current_curve->getSecondPoint()->getPoint()<<std::endl;
		}
		else
		{
			// in this case we have turned around and treated all the curves
			loopDone = true;
		}
	} // while (!loopDone)

	if((ACur.size() > 1) && (ACur.back() == ACur.front())) {
		ACur.pop_back();
	}

	// we check that no curve was overlooked
	for (size_t iCurve=0; iCurve<curves_treated.size(); iCurve++)
	{
		if (!curves_treated[iCurve])
		{
			throw GMDSException("GeomPoint::getOrdered get ordered curves failed");
		}
	}
}
/*----------------------------------------------------------------------------*/
void GeomPoint::getOrderedDirect(std::vector<GeomCurve*>& ACur) const
{
	// we empty the result vector
	ACur.clear();

	// we get the ordered curves
	std::vector<GeomCurve* > curvesOrdered;
	getOrdered(curvesOrdered);

	// we now determine whether the order is direct or indirect
	if(curvesOrdered.size() == 0)
		return;

	if(curvesOrdered.size() == 1)
	{
		ACur.push_back(curvesOrdered[0]);
		return;
	}

	if(curvesOrdered.size() == 2) {
		ACur.push_back(curvesOrdered[0]);
		ACur.push_back(curvesOrdered[1]);
		return;
	}

	// find a curve that has more than one surface
	// and take the next curve
	gmds::geom::GeomCurve* curve1 = NULL;
	gmds::geom::GeomCurve* curve2 = NULL;

	for(unsigned int iCurve=0; iCurve<curvesOrdered.size(); iCurve++) {
		std::vector<GeomSurface* > surfaces;
		curve1 = curvesOrdered[iCurve];
		curve1->get(surfaces);

		// if there is only one adjacent surface, this curve is not a good choice 
		// because the surface will both be on the left and right.
		if(surfaces.size()>1) {
			curve2 = curvesOrdered[(iCurve+1)%curvesOrdered.size()];

			// we here check that there is only one surface in common between curve1 and cuve2;
			// necessary for the time being because we will determine the orientation of the 
			// order based on the position (left or right) of the surface in common related to curve1;
			// with two surfaces in common, which one should we choose?
			std::vector<GeomSurface*> surfacesInCommon;
			int nbCommonSurfaces = curve1->commonSurfaces(curve2,surfacesInCommon);
			if(nbCommonSurfaces == 1) {
				break;
			} else {
				curve1 = NULL;
				curve2 = NULL;
			}	
		}
	}

	if(NULL == curve2) {
		throw GMDSException("GeomPoint::getOrderedDirect there is no curve "
				"that has more than one surface or all ordered curves pairs " 
				"had more than one surface in common.");
	}

	GeomPoint* p1_first = curve1->getFirstPoint();

	// we first the surface shared by these two curves
	// and we check whether it is on the left or right of curve1
	std::vector<GeomSurface*> surfacesInCommon;
	curve1->commonSurfaces(curve2,surfacesInCommon);
	gmds::geom::GeomSurface* commonSurface = surfacesInCommon[0];	

	bool isDirect;

	if(this == p1_first) {
		if(commonSurface == curve1->getLeftSurface()) {
			isDirect = true;
		} else {
			isDirect = false;
		}
	} else {
		if(commonSurface == curve1->getLeftSurface()) {
			isDirect = false;
		} else {
			isDirect = true;
		}
	}

	if(isDirect)
	{
		ACur.insert(ACur.begin(),curvesOrdered.begin(),curvesOrdered.end());
	}
	else
	{
		ACur.insert(ACur.begin(),curvesOrdered.rbegin(),curvesOrdered.rend());
	}

}
/*----------------------------------------------------------------------------*/
TInt GeomPoint::getNbCurvesOfSurface(GeomSurface& ASurf) const
{
	std::vector<GeomCurve*> curves;
	this->get(curves);

	std::vector<GeomCurve*> surfCurves;
	ASurf.get(surfCurves);

	unsigned int nbCurvesOfSurface = 0;

	for(unsigned int iSurfCurve=0; iSurfCurve<surfCurves.size(); iSurfCurve++) {
		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			if(surfCurves[iSurfCurve] == curves[iCurve]) {
				nbCurvesOfSurface++;
			}
		}
	}

	return nbCurvesOfSurface;
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
