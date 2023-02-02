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
 * IGMeshQualityEvaluation.cpp
 *
 *  Created on: 22 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMeshQualityEvaluation.h>
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Matrix.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
IGMeshQualityEvaluation::IGMeshQualityEvaluation()
{}
/*----------------------------------------------------------------------------*/
IGMeshQualityEvaluation::~IGMeshQualityEvaluation()
{}
/*----------------------------------------------------------------------------*/
TCoord IGMeshQualityEvaluation::minAngle(const Face& AF) const
{
	std::vector<Node> nodes;
	AF.get<Node>(nodes);
	std::vector<math::Point> p;
	p.resize(nodes.size());
	for(unsigned int i=0;i<nodes.size();i++)
		p[i]=nodes[i].getPoint();

	TCoord min = math::Vector(p[0],p[1]).angle(math::Vector(p[0],p.back()));
	for(unsigned int i=1;i<p.size();i++){

		math::Point current = p[i];
		math::Point prev    = p[i-1];
		math::Point next    = p[(i+1)%p.size()];

		TCoord current_angle = math::Vector(current,next).angle(math::Vector(current,prev));
		if(current_angle<min)
			min=current_angle;
	}
	return min;
}
/*----------------------------------------------------------------------------*/
TCoord IGMeshQualityEvaluation::maxAngle(const Face& AF) const
{
	std::vector<Node> nodes;
	AF.get<Node>(nodes);
	std::vector<math::Point> p;
	p.resize(nodes.size());
	for(unsigned int i=0;i<nodes.size();i++)
		p[i]=nodes[i].getPoint();

	TCoord max = math::Vector(p[0],p[1]).angle(math::Vector(p[0],p.back()));
	for(unsigned int i=1;i<p.size();i++){

		math::Point current = p[i];
		math::Point prev    = p[i-1];
		math::Point next    = p[(i+1)%p.size()];

		TCoord current_angle = math::Vector(current,next).angle(math::Vector(current,prev));
		if(current_angle>max)
			max=current_angle;
	}
	return max;
}
/*----------------------------------------------------------------------------*/
TCoord IGMeshQualityEvaluation::aspectRatio(const Face& AF) const
{
	std::vector<Node> nodes;
	AF.get<Node>(nodes);
	std::vector<math::Point> p;
	p.resize(nodes.size());
	for(unsigned int i=0;i<nodes.size();i++)
		p[i]=nodes[i].getPoint();


	TCoord min_length = math::Vector(p[0],p[1]).norm();
	TCoord max_length = math::Vector(p[0],p[1]).norm();
	for(unsigned int i=1;i<p.size();i++){

		math::Point current = p[i];
		math::Point next    = p[(i+1)%p.size()];
		TCoord cur_length = math::Vector(current,next).norm();
		if(cur_length<min_length)
			min_length=cur_length;
		if(cur_length>max_length)
			max_length=cur_length;

	}
	return max_length/min_length;
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::jacobian(const Region& AR) const
{
	if(AR.getType()==GMDS_TETRA)
		return jacobianTet(AR);
	else if(AR.getType()==GMDS_HEX)
		return jacobianHex(AR);
	else
		throw GMDSException("Jacobian can be only computed for tetrahedral and hexahedral elements");
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::scaledJacobian(const Region& AR) const
{
	if(AR.getType()==GMDS_TETRA)
		return scaledJacobianTet(AR);
	else if(AR.getType()==GMDS_HEX)
		return scaledJacobianHex(AR);
	else
		throw GMDSException("Scaled jacobian can be only computed for tetrahedral and hexahedral elements");}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::jacobianTet(const Region& AR) const
{
	std::vector<Node> n;
	AR.get<Node>(n);
	return jacobian(n[0],n[1],n[2],n[3]);
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::jacobianHex(const Region& AR) const
{
	std::vector<Node> n;
	AR.get<Node>(n);
	double j[8]= {jacobian(n[0],n[1],n[3],n[4]),
			jacobian(n[1],n[2],n[0],n[5]),
			jacobian(n[2],n[3],n[1],n[6]),
			jacobian(n[3],n[0],n[2],n[7]),
			jacobian(n[4],n[7],n[5],n[0]),
			jacobian(n[5],n[4],n[6],n[1]),
			jacobian(n[6],n[5],n[7],n[2]),
			jacobian(n[7],n[6],n[4],n[3])};
	double min = j[0];
	for(int i=1;i<8;i++){
		if(j[i]<min)
			min=j[i];
	}
	return min;
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::scaledJacobianTet(const Region& AR) const
{
	std::vector<Node> n;
	AR.get<Node>(n);
	TCoord norm1 = math::Vector(n[0].getPoint(),n[1].getPoint()).norm();
	TCoord norm2 = math::Vector(n[0].getPoint(),n[2].getPoint()).norm();
	TCoord norm3 = math::Vector(n[0].getPoint(),n[3].getPoint()).norm();
	if(norm1==0. || norm2==0. || norm3==0.) {
		return 0.;
	} else {
		return jacobian(n[0],n[1],n[2],n[3])/(norm1*norm2*norm3);
	}
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::scaledJacobianHex(const Region& AR) const
{
	std::vector<Node> n;
	AR.get<Node>(n);
	TCoord norm01 = math::Vector(n[0].getPoint(),n[1].getPoint()).norm();
	TCoord norm03 = math::Vector(n[0].getPoint(),n[3].getPoint()).norm();
	TCoord norm04 = math::Vector(n[0].getPoint(),n[4].getPoint()).norm();
	TCoord norm12 = math::Vector(n[1].getPoint(),n[2].getPoint()).norm();
	TCoord norm15 = math::Vector(n[1].getPoint(),n[5].getPoint()).norm();
	TCoord norm23 = math::Vector(n[2].getPoint(),n[3].getPoint()).norm();
	TCoord norm26 = math::Vector(n[2].getPoint(),n[6].getPoint()).norm();
	TCoord norm37 = math::Vector(n[3].getPoint(),n[7].getPoint()).norm();
	TCoord norm45 = math::Vector(n[4].getPoint(),n[5].getPoint()).norm();
	TCoord norm47 = math::Vector(n[4].getPoint(),n[7].getPoint()).norm();
	TCoord norm56 = math::Vector(n[5].getPoint(),n[6].getPoint()).norm();
	TCoord norm67 = math::Vector(n[6].getPoint(),n[7].getPoint()).norm();

	double j[8];
	j[0] = (norm12*norm01*norm15==0.) ? 0. : jacobian(n[0],n[1],n[3],n[4])/(norm01*norm03*norm04);
	j[1] = (norm12*norm01*norm15==0.) ? 0. : jacobian(n[1],n[2],n[0],n[5])/(norm12*norm01*norm15);
	j[2] = (norm23*norm12*norm26==0.) ? 0. : jacobian(n[2],n[3],n[1],n[6])/(norm23*norm12*norm26);
	j[3] = (norm03*norm23*norm37==0.) ? 0. : jacobian(n[3],n[0],n[2],n[7])/(norm03*norm23*norm37);
	j[4] = (norm47*norm45*norm04==0.) ? 0. : jacobian(n[4],n[7],n[5],n[0])/(norm47*norm45*norm04);
	j[5] = (norm45*norm56*norm15==0.) ? 0. : jacobian(n[5],n[4],n[6],n[1])/(norm45*norm56*norm15);
	j[6] = (norm56*norm67*norm26==0.) ? 0. : jacobian(n[6],n[5],n[7],n[2])/(norm56*norm67*norm26);
	j[7] = (norm67*norm47*norm37==0.) ? 0. : jacobian(n[7],n[6],n[4],n[3])/(norm67*norm47*norm37);
	
	double min = j[0];
	for(int i=1;i<8;i++){
		if(j[i]<min)
			min=j[i];
	}
	return min;
}
/*----------------------------------------------------------------------------*/
double IGMeshQualityEvaluation::
jacobian(const Node& AN0,const Node& AN1,const Node& AN2,const Node& AN3) const
{
	TCoord mat[3][3] = {
			{AN1.X()-AN0.X(), AN2.X()-AN0.X(), AN3.X()-AN0.X()},
			{AN1.Y()-AN0.Y(), AN2.Y()-AN0.Y(), AN3.Y()-AN0.Y()},
			{AN1.Z()-AN0.Z(), AN2.Z()-AN0.Z(), AN3.Z()-AN0.Z()} };

	math::Matrix<3,3,double> jacobianMatrix(mat);
	return jacobianMatrix.det();
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/



