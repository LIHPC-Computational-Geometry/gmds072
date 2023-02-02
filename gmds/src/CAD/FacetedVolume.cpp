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
/** \file    FacetedVolume.t.h
 *  \author  F. LEDOUX
 *  \date    29/06/2011
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedVolume.h>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
int FacetedVolume::next_id_=1;
/*----------------------------------------------------------------------------*/
FacetedVolume::FacetedVolume()
{}
/*----------------------------------------------------------------------------*/
FacetedVolume::FacetedVolume(std::vector<FacetedSurface* >& AS,
		const std::string& AName)
:GeomVolume(AName),surfaces_(AS)
{}
/*----------------------------------------------------------------------------*/
FacetedVolume::~FacetedVolume()
{}
/*----------------------------------------------------------------------------*/
void FacetedVolume::get(std::vector<GeomSurface*>& ASurf) const
{
	ASurf.clear();
	ASurf.resize(surfaces_.size());
	for(unsigned int i=0;i<surfaces_.size();i++){
		ASurf[i]=surfaces_[i];
	}
}
/*----------------------------------------------------------------------------*/
void FacetedVolume::get(std::vector<GeomCurve*>& AC) const
{
	std::set<GeomCurve* > curves_set;
	for(unsigned int i=0;i<surfaces_.size();i++)
	{
		std::vector<GeomCurve* > c_i;
		surfaces_[i]->get(c_i);
		curves_set.insert(c_i.begin(),c_i.end());
	}
	AC.clear();
	AC.insert(AC.begin(),curves_set.begin(),curves_set.end());
}
/*----------------------------------------------------------------------------*/
void FacetedVolume::get(std::vector<GeomPoint*>& AC) const
{
	std::set<GeomPoint* > points_set;
	for(unsigned int i=0;i<surfaces_.size();i++)
	{
		std::vector<GeomPoint* > p_i;
		surfaces_[i]->get(p_i);
		points_set.insert(p_i.begin(),p_i.end());
	}
	AC.clear();
	AC.insert(AC.begin(),points_set.begin(),points_set.end());
}
/*----------------------------------------------------------------------------*/
TCoord FacetedVolume::computeArea() const
{
	throw GMDSException("Not yet implemented!");
}
/*----------------------------------------------------------------------------*/
void FacetedVolume::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	minXYZ[0] =  HUGE_VALF;
	minXYZ[1] =  HUGE_VALF;
	minXYZ[2] =  HUGE_VALF;
	maxXYZ[0] = -HUGE_VALF;
	maxXYZ[1] = -HUGE_VALF;
	maxXYZ[2] = -HUGE_VALF;
	
	for(unsigned int i=0;i<surfaces_.size();i++){
		TCoord minXYZbis[3];
		TCoord maxXYZbis[3];
		surfaces_[i]->computeBoundingBox(minXYZbis,maxXYZbis);
	
		if(minXYZbis[0] < minXYZ[0]) minXYZ[0] = minXYZbis[0];
		if(minXYZbis[1] < minXYZ[1]) minXYZ[1] = minXYZbis[1];
		if(minXYZbis[2] < minXYZ[2]) minXYZ[2] = minXYZbis[2];
		if(maxXYZbis[0] > maxXYZ[0]) maxXYZ[0] = maxXYZbis[0];
                if(maxXYZbis[1] > maxXYZ[1]) maxXYZ[1] = maxXYZbis[1];
                if(maxXYZbis[2] > maxXYZ[2]) maxXYZ[2] = maxXYZbis[2];
	}

}
/*----------------------------------------------------------------------------*/
//void FacetedVolume::getTriangulation(std::vector<math::Triangle<3,TCoord> >& ATri) const
//{
//	ATri.clear();
//	for(unsigned int i=0;i<surfaces_.size();i++){
//		std::vector<math::Triangle<3,TCoord> > tri;
//		surfaces_[i]->getTriangulation(tri);
//		ATri.insert(ATri.begin(),tri.begin(),tri.end());
//	}
//}
/*----------------------------------------------------------------------------*/
void FacetedVolume::getMeshFaces(std::vector<Face>& AFaces) const
{
	AFaces.clear();
	for(unsigned int i=0;i<surfaces_.size();i++){
		std::vector<Face> faces;
		surfaces_[i]->getMeshFaces(faces);
		AFaces.insert(AFaces.begin(),faces.begin(),faces.end());
	}
}
/*----------------------------------------------------------------------------*/
void 
FacetedVolume::reorient()
{
	for(int iSurf=0; iSurf<surfaces_.size(); iSurf++) {
		surfaces_[iSurf]->reorient(this);
	}
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
