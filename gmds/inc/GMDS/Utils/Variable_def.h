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
/* Variable_def.h
 *
 *  Created on: 26 July 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
/** \class VariableItf
 *  \brief Defines the API of a mesh variable.
 */
class EXPORT_GMDS VariableItf{
public:

	enum export_type{
		var_int, 
		var_double, 
		var_double_vec,
		var_cross,
		var_cross_2D,
		var_quaternion,
		var_unknown
	};

	VariableItf(){;}
	virtual ~VariableItf(){;}


	/*------------------------------------------------------------------------*/
	/** \brief  Gives acces to the type of variable if it has been specified
	*/
	export_type getType();

	/*------------------------------------------------------------------------*/
	/** \brief  Accessor to the variable name
	 */
	virtual std::string getName() const=0;
	/*------------------------------------------------------------------------*/
	/** \brief  Add a new entry for this variable
	 */
	virtual void addEntry	(const int& i)=0;
	/*------------------------------------------------------------------------*/
	/** \brief  Remove an entry for this variable
	 */
	virtual void removeEntry(const int& i)=0;
	/*------------------------------------------------------------------------*/
	/** \brief Set the domain entry to [0,i]
	 */
	virtual void setDomain	(const int& i)=0;
	/*------------------------------------------------------------------------*/
	/** \brief Get the domain size
	 */
	virtual int getDomainSize() const=0;
	/*------------------------------------------------------------------------*/
	/** \brief Set the domain entry to [0,i] with some default values
	 */
	virtual void setDomainWithDefault(const int& i,
									  const std::vector<int>& def)=0;
	/*------------------------------------------------------------------------*/
	/** \brief Clear all the domain
	 */
	virtual void clear()=0;
	/*------------------------------------------------------------------------*/
	/** \brief Serialize the variable into stream str. Warning this method does
	 * 		   not support typename T where pointers would be present.
	 */
	virtual void serialize(std::ostream& stream)  =0;
	/*------------------------------------------------------------------------*/
	/** \brief Unserialize the variable from stream str. Warning this method
	 * 		   does not support typename T where pointers would be present.
	 */
	virtual void unserialize(std::istream& stream)=0;
	/*------------------------------------------------------------------------*/
	/** \brief compact the variable data
	 */
	virtual void compact()=0;

};
/*----------------------------------------------------------------------------*/
/** \class Variable
 *  \brief Defines what a mesh variable is
 *
 *  \param T  type of the variable
 */
template<typename T>
class EXPORT_GMDS Variable : public VariableItf{
public:

	/*------------------------------------------------------------------------*/
	/** \brief  Constructor. Every variable has a name
	 */
	Variable(const std::string& AName="no name"):m_name(AName){;}

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.
	 */
	virtual ~Variable(){;}

	/*------------------------------------------------------------------------*/
	/** \brief  overload operator[] const
	 */
	T const & operator[](const int& i) const {
		return m_data[i];
	}

	/*------------------------------------------------------------------------*/
	/** \brief  overload operator[]. It allows to modify a variable value
	 */
	T& operator[](const int& i)  {
		return m_data[i];
	}

	/*------------------------------------------------------------------------*/
	/** \brief Update all the elements to value AVal
	 */
	void setValuesTo(const T& AVal){
		typename SmartVector<T>::iterator it =m_data.begin();

		for(;!it.isDone();it.next()){
			T ptr = it.currentItem();
			(ptr)=AVal;
		}
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Accessor to the variable name
	 */
	std::string getName() const{
		return m_name;
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of values the variable is defined on
	 */
	int getNbValues() const{
		return m_data.size();
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Add a new entry for this variable
	 */
	virtual void addEntry(const int& i){
		if(m_data.isOutOfContainer(i))
			m_data.resize(2*i);
		//This choice is maybe too expensive
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Remove an entry for the variable
	 */
	virtual void removeEntry(const int& i){
		m_data.remove(i);
	}
	/*------------------------------------------------------------------------*/
	/** \brief Set entry i to value val
	 */
	void set(const int& i, const T& val)
	{
		m_data.assign(val,i);
	}

	/*------------------------------------------------------------------------*/
	/** \brief Set the domain entry to [0,i]
	 */
	void setDomain(const int& i){
		m_data.resize(i);
	}

	/*------------------------------------------------------------------------*/
	/** \brief Get the domain size
	 */
	virtual int getDomainSize() const{
		return m_data.capacity();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Set the domain entry to [0,i] with some default values
	 */
	void setDomainWithDefault(const int& i, const std::vector<int>& def){
		m_data.resize(i);

		if(def.size()>i)
			return;

		m_data.mark(def);

		// to make iterators consistent
		m_data.update();
	}
	/*------------------------------------------------------------------------*/
	/** \brief Clear all the domain
	 */
	void clear(){
		m_data.clear();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Serialize the variable into stream str. Warning this method does
	 * 		   not support typename T where pointers would be present.
	 */
	void serialize(std::ostream& stream);

	/*------------------------------------------------------------------------*/
	/** \brief Unserialize the variable from stream str. Warning this method
	 * 		   does not support typename T where pointers would be present.
	 */
	void unserialize(std::istream& stream);

	/*------------------------------------------------------------------------*/
	/** \brief compact the variable data
	 */
	void compact();

private:

	/* variable name*/
	std::string m_name;

	/* data*/
	SmartVector<T> m_data;
};
/*----------------------------------------------------------------------------*/


