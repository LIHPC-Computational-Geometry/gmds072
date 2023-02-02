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
/** \file    FacetedPoint.t.h
 *  \author  N. LE GOFF
 *  \date    03/04/2012
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedPoint.h>
#include <GMDS/CAD/FacetedCurve.h>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
int FacetedPoint::next_id_=1;
/*----------------------------------------------------------------------------*/
FacetedPoint::FacetedPoint()
: id_(next_id_++)
{

}
/*----------------------------------------------------------------------------*/
FacetedPoint::FacetedPoint(Node ANode, const std::string& AName)
: GeomPoint(AName),id_(next_id_++),pnt_(ANode)
{

}
/*----------------------------------------------------------------------------*/
FacetedPoint::~FacetedPoint()
{

}
/*----------------------------------------------------------------------------*/
void FacetedPoint::get(std::vector<GeomCurve*>& ACur) const
{
	ACur.clear();
	ACur.resize(curves_.size());
	for(unsigned int i=0;i<curves_.size();i++)
		ACur[i]=curves_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedPoint::get(std::vector<GeomSurface*>& ASurf) const{
	std::set<GeomSurface* > surf_set;
	for(unsigned int i=0;i<curves_.size();i++)
	{
		std::vector<GeomSurface* > surf_i;
		curves_[i]->get(surf_i);
		surf_set.insert(surf_i.begin(),surf_i.end());
	}
	ASurf.clear();
	ASurf.insert(ASurf.begin(),surf_set.begin(),surf_set.end());
}
/*----------------------------------------------------------------------------*/
void FacetedPoint::get(std::vector<GeomVolume*>& AVol) const{
	std::set<GeomVolume* > vol_set;
	for(unsigned int i=0;i<curves_.size();i++)
	{
		std::vector<GeomVolume* > vol_i;
		curves_[i]->get(vol_i);
		vol_set.insert(vol_i.begin(),vol_i.end());
	}
	AVol.clear();
	AVol.insert(AVol.begin(),vol_set.begin(),vol_set.end());
}
/*----------------------------------------------------------------------------*/
TInt
FacetedPoint::getNbCurves() const
{
	return curves_.size();
}
/*----------------------------------------------------------------------------*/
TCoord
FacetedPoint::X() const
{
	return pnt_.X();
}
/*----------------------------------------------------------------------------*/
TCoord
FacetedPoint::Y() const
{
        return pnt_.Y();
}
/*----------------------------------------------------------------------------*/
TCoord
FacetedPoint::Z() const
{
        return pnt_.Z();
}
/*----------------------------------------------------------------------------*/
void
FacetedPoint::XYZ(TCoord ACoordinates[3]) const
{
	ACoordinates[0] = pnt_.X();
        ACoordinates[1] = pnt_.Y();
        ACoordinates[2] = pnt_.Z();
}
/*----------------------------------------------------------------------------*/
TCoord
FacetedPoint::computeArea() const
{
	return 0;
}
/*----------------------------------------------------------------------------*/
void
FacetedPoint::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	minXYZ[0]=pnt_.X(); maxXYZ[0]=pnt_.X();
	minXYZ[1]=pnt_.Y(); maxXYZ[1]=pnt_.Y();
	minXYZ[2]=pnt_.Z(); maxXYZ[2]=pnt_.Z();
}
/*----------------------------------------------------------------------------*/
gmds::math::Point
FacetedPoint::getPoint() const
{
	return gmds::math::Point(pnt_.X(),pnt_.Y(),pnt_.Z());
}
/*----------------------------------------------------------------------------*/
void
FacetedPoint::set(Node ANode)
{
	pnt_ = ANode;
}
/*----------------------------------------------------------------------------*/
void
FacetedPoint::add(FacetedCurve* ACurve)
{
	curves_.push_back(ACurve);
}
/*----------------------------------------------------------------------------*/
void
FacetedPoint::add(GeomCurve*)
{
	throw GMDSException("FacetedPoint::add(GeomCurve*) should not be called");
}
/*----------------------------------------------------------------------------*/
Node
FacetedPoint::getNode() const
{
	return pnt_;
}
/*----------------------------------------------------------------------------*/
int
FacetedPoint::getId() const
{
	return id_;
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
