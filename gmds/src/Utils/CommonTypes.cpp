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
 * CommonTypes.cpp
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <algorithm>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
bool isZero(const TCoord a){
	return fabs(a)<TCoord_Epsilon;
}
/*----------------------------------------------------------------------------*/
bool areEquals(const TCoord a,const TCoord b){
	return fabs(a-b)<TCoord_Epsilon;
}
/*----------------------------------------------------------------------------*/
TCoord min3(const TCoord a,const TCoord b, const TCoord c){
       	if(a<b) {
		if(a<c) {
			return a;
		} else {
			return c;	
		}
	} else {
		if(b<c) {
			return b;
		} else {
			return c;
		}
	}
}
/*----------------------------------------------------------------------------*/
TCoord max3(const TCoord a,const TCoord b, const TCoord c){
        if(a>b) {
                if(a>c) {
                        return a;
                } else {
                        return c;
                }
        } else {
                if(b>c) {
                        return b;
                } else {
                        return c;
                }
        }
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
getCommonBut(const std::vector<TCellID>& AS1,const std::vector<TCellID>& AS2, const TCellID ABut)
{
	std::vector<TCellID> result;
	unsigned int size1 = AS1.size();
	unsigned int size2 = AS2.size();
	for (unsigned int i1=0; i1<size1; i1++){
		TCellID id1= AS1[i1];
		if(id1!=ABut)
		{
			bool found = false;
			for (unsigned int i2=0; i2<size2 && !found; i2++){
				TCellID id2= AS2[i2];
				if(id1==id2){
					result.push_back(id1);
					found =true;
				}
			}
		}
	}
	return result;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
keepFilter(const std::vector<TCellID>& ASet, const TInt ANb)
{
	std::vector<TCellID> result;

	unsigned int s = ASet.size();

	for (unsigned int i=0; i<s; i++){
		TCellID current_id = ASet[i];
		int nb_occ=1;
		for (unsigned int j=i+1; j<s; j++){
			TCellID other_id = ASet[j];
			if(other_id==current_id)
				nb_occ++;
		}

		if(nb_occ>=ANb && std::find(result.begin(),result.end(),current_id)==result.end())
			result.push_back(current_id);
	}
	return result;
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, const MeshModel& model) {
	int m = model.getDef();
	stream << "|";
	if(m&N)
		stream<<"N|";
	if(m&E)
		stream<<"E|";
	if(m&F)
		stream<<"F|";
	if(m&R)
		stream<<"R|";
	if(m&N2E)
		stream<<"N2E|";
	if(m&N2F)
		stream<<"N2F|";
	if(m&N2R)
		stream<<"N2R|";
	if(m&E2N)
		stream<<"E2N|";
	if(m&E2F)
		stream<<"E2F|";
	if(m&E2R)
		stream<<"E2R|";
	if(m&F2N)
		stream<<"F2N|";
	if(m&F2E)
		stream<<"F2E|";
	if(m&F2F)
		stream<<"F2F|";
	if(m&F2R)
		stream<<"F2R|";
	if(m&R2N)
		stream<<"R2N|";
	if(m&R2E)
		stream<<"R2E|";
	if(m&R2F)
		stream<<"R2F|";
	if(m&R2R)
		stream<<"R2R|";

	return stream;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

