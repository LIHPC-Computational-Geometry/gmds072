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
 * VariableManager.cpp
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/VariableManager.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
VariableManager::VariableManager(){;}
/*----------------------------------------------------------------------------*/
VariableManager::~VariableManager()
{
	for(unsigned int k=0;k<m_variables.size();k++)
	{
		VariableItf *v =m_variables[k];
		delete v;
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::deleteVariable(const std::string& AName){
	for(unsigned int k=0;k<m_variables.size();k++)
	{
		if (m_variables[k]->getName()==AName){
			VariableItf *v =m_variables[k];
			if(k!=m_variables.size()-1)
				m_variables[k] = m_variables.back();
			m_variables.pop_back();
			delete v;
		}
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::deleteVariable(VariableItf* AVar){
	for(unsigned int k=0;k<m_variables.size();k++)
	{
		if (m_variables[k]==AVar){
			VariableItf *v =m_variables[k];
			if(k!=m_variables.size()-1)
				m_variables[k] = m_variables.back();
			m_variables.pop_back();
			delete v;
		}
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::addEntry(const int i){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->addEntry(i);
}
/*----------------------------------------------------------------------------*/
void VariableManager::removeEntry(const int i){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->removeEntry(i);
}
/*----------------------------------------------------------------------------*/
void VariableManager::setDomain(const int size){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->setDomain(size);
}
/*----------------------------------------------------------------------------*/
void VariableManager::clearVariables(){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->clear();
}
/*----------------------------------------------------------------------------*/
void VariableManager::compact(){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->compact();
}
/*----------------------------------------------------------------------------*/
void VariableManager::serialize(std::ostream& AStr) {
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void VariableManager::unserialize(std::istream& AStr){

	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->unserialize(AStr);
}
/*----------------------------------------------------------------------------*/
int VariableManager::getNbVariables() const{
	return m_variables.size();
}
/*----------------------------------------------------------------------------*/
VariableItf* VariableManager::getVariable(const TInt i){
	if(i<0 || i>=m_variables.size())
		throw GMDSException("VariableManager::getVariable - Out of the bounds");
	return m_variables[i];
}
/*----------------------------------------------------------------------------*/
bool VariableManager::doesVariableExist(const std::string& AName){
	unsigned int k=0;
	for(;k<m_variables.size();k++)
		if (m_variables[k]->getName()==AName)
			return true;

	return false;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
