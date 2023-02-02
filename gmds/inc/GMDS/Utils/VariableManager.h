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
/*
 * VariableManager.h
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VARIABLEMANAGER_H_
#define GMDS_VARIABLEMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/Variable.h>
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Utils/Exception.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
	/*----------------------------------------------------------------------------*/
	/** \class VariableManager
	 *  \brief Handle the creation and update of a collection of variables. A
	 *  	   variable is defined as a set of discrete values associated to a key.
	 *  	   Few holes are in the key numerotation.
	 */
	class EXPORT_GMDS VariableManager{

	public:
		/*------------------------------------------------------------------------*/
		/** \brief  Constructor.
		*/
		VariableManager();

		/*------------------------------------------------------------------------*/
		/** \brief  Destructor.
		*/
		~VariableManager();

		/*------------------------------------------------------------------------*/
		/** \brief  creation of a variable allocated in the stack. The domain of the
		 * 			variable is initialized to [0,initSize].
		 */
		template<typename T> Variable<T>* newVariable(const std::string& AName,
			const int initSize = 2, const std::vector<int>* ref = 0);

		/*------------------------------------------------------------------------*/
		/** \brief  Returns whether a variable exists.
		*/
		bool doesVariableExist(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  Access to a variable.
		*/
		template<typename T> EXPORT_GMDS Variable<T>* getVariable(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  suppression of a variable. the memory used in the stack is free.
		*/
		void deleteVariable(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  suppression of a variable. the memory used in the stack is free.
		*/
		void deleteVariable(VariableItf* AVar);

		/*------------------------------------------------------------------------*/
		/** \brief  Add a default entry to all the variables. Each variable has the
		 * 			responsability to initialize the corresponding value.
		 */
		void addEntry(const int i);

		/*------------------------------------------------------------------------*/
		/** \brief  Remove the i th entry of all the managed variables.
		 */
		void removeEntry(const int i);

		/*------------------------------------------------------------------------*/
		/** \brief  set the domain of all the variables to size.
		 */
		void setDomain(const int size);

		/*------------------------------------------------------------------------*/
		/** \brief  set the domain of all the variables to size.
		 */
		int getNbVariables() const;

		/*------------------------------------------------------------------------*/
		/** \brief  get access to the i^th variable in a abstract form
		 */
		VariableItf* getVariable(const TInt i);

		/*------------------------------------------------------------------------*/
		/** \brief indicates if variables are attached
		 */
		bool empty() const;

		/*------------------------------------------------------------------------*/
		/** \brief clear all the variable
		 */
		void clearVariables();

		/*------------------------------------------------------------------------*/
		/** \brief compact all the variables
		 */
		void compact();

		/*------------------------------------------------------------------------*/
		/** \brief serialize (*this) in AStr
		 *
		 * \param AStr an output streammap
		 */
		void serialize(std::ostream& AStr);

		/*------------------------------------------------------------------------*/
		/** \brief unserialize (*this) from AStr
		 *
		 * \param AStr an input stream
		 */
		void unserialize(std::istream& AStr);

		/*------------------------------------------------------------------------*/
		/** \brief  get the list of variables
		*/
		std::vector<VariableItf*> getAllVariables(){ return m_variables; }

	private:

		std::vector<VariableItf*>  m_variables;
	};
	/*----------------------------------------------------------------------------*/
	template<typename T>
	Variable<T>* VariableManager::newVariable(const std::string& AName,
		const int initSize, const std::vector<int>* ref){

		for (unsigned int k = 0; k < m_variables.size(); k++)
            if (m_variables[k]->getName() == AName){
                std::string mess= "Impossible to create a variable "+AName+": name already used";
			throw GMDSException(mess);
            }

		Variable<T>* v = new Variable<T>(AName);
		m_variables.push_back(v);

		if (ref == 0)
			v->setDomain(initSize + 1);
		else
			v->setDomainWithDefault(initSize + 1, *ref);

		return v;
	}
	/*----------------------------------------------------------------------------*/
	template<typename T>
	Variable<T>* VariableManager::getVariable(const std::string& AName){
		unsigned int k = 0;
		for (; k < m_variables.size(); k++)
		if (m_variables[k]->getName() == AName)
			return dynamic_cast<Variable<T>*>(m_variables[k]);
        std::string mess= "No variable named "+AName;
		throw GMDSException(mess);
	}
	/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_VARIABLEMANAGER_H_ */
/*----------------------------------------------------------------------------*/
