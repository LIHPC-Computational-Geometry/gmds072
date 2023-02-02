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
 *  Parameters.h
 *
 *  Created on: April 12, 2016
 *  Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PARAMETERS_H_
#define GMDS_PARAMETERS_H_
/*----------------------------------------------------------------------------*/
// STL Headers
#include <string>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    class Parameters {
        
    public:
        /*--------------------------------------------------------------------*/
        /** \enum ETypeParam
         *  \brief enum type defining the different type of parameters we handle
         */
        enum ETypeParam{
            INT_P, DOUBLE_P, STRING_P, BOOL_P
        };
        
        /*--------------------------------------------------------------------*/
        /** \struct Entry
         *  \brief describe a parameters entry
         */
        struct Entry{
            std::string section;
            std::string name;
            ETypeParam type;
            
        };
        
        /*--------------------------------------------------------------------*/
        /** \brief Default constructor
         */
        Parameters();
        /*--------------------------------------------------------------------*/
        /** \brief Destructor
         */
        ~Parameters();
        
        /*--------------------------------------------------------------------*/
        /** \brief add a new param which must be in section \p ASection with
         *         name \p AName
         *
         * \param[in] ASection parameter section
         * \param[in] AName    parameter name
         * \param[in] AType    parameter type
         *
         * \return true if the insertion succeeds, false if an existing entry with
         *         the same section and name already exists
         */
        bool add(const std::string& ASection,
                 const std::string& AName,
                 const ETypeParam AType);
        /*--------------------------------------------------------------------*/
        /** \brief get the list of parameter entries, that is the description of
         *         the parameters (but not their value)
         *
         * \return the list of entries
         */
        std::vector<Entry> getEntries() const;
        
        
        /*--------------------------------------------------------------------*/
        /** \brief Returns the value of the parameter (\p ASection, \p AName) if
         *         it exists. Return value is in \p AOut
         *
         * \param[in]  ASection parameter section
         * \param[in]  AName    parameter name
         * \param[out] AOut     return value (int, double or string)
         *
         * \return true if the value is found, false otherwise
         */
        bool get(const std::string& ASection, const std::string& AName,
                 int& AOut);
        
        bool get(const std::string& ASection, const std::string& AName,
                 double& AOut);
        
        bool get(const std::string& ASection, const std::string& AName,
                 std::string& AOut);
        
        bool get(const std::string& ASection, const std::string& AName,
                 bool& AOut);
        
        /*--------------------------------------------------------------------*/
        /** \briefInitialize the parameters from the file \p AFileName
         *
         * \param[in] AFileName file to be read
         *
         * \return the list of not found entries, or found with a wrong type
         */
        std::vector<std::string> parseIni(const std::string& AFileName);
        
        
        private :
        
        /** value indicating that the parameters have been initialized */
        bool m_initialized;
        
        /** Entries that define the parameter we look for */
        std::vector<Entry> m_entries;
        /** each entry value is stored as a std::string */
        std::vector<std::string> m_values;
        /*--------------------------------------------------------------------*/
        /** \brief Find the entry (\p ASection, \p AName) and returns its
         *         position in m_entries
         *
         * \param[in]  ASection parameter section
         * \param[in]  AName    parameter name
         *
         * \return entry position and -1 if it is not found
         */
        int find(const std::string& ASection, const std::string& AName);
        
    };
    /*------------------------------------------------------------------------*/
} //namespace gmds
/*----------------------------------------------------------------------------*/
#endif // GMDS_PARAMETERS_H_
/*----------------------------------------------------------------------------*/
