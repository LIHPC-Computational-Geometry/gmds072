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
/** \file    MeditReader_def.h
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */

/*----------------------------------------------------------------------------*/
template<typename TMesh>
class MeditReader{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param mesh the mesh in which we want to copy the content of a Lima
     *  	   file.
     */
	MeditReader(TMesh& mesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~MeditReader();

    /*------------------------------------------------------------------------*/
    /** \brief  Read the content of the file named fileName and write it in
     *   		m_mesh.
     *
     *   \param fileName 	name of the file to read
     *   \param mask	 	definition of what we read in the file
     *   				 	(for instance, N|F = nodes + faces)
     *   \param useLabels	uses or not the labels of the medit file
     */
	void read(const std::string& fileName, int mask, bool useLabels=false);

protected:

	void readNodes(std::ifstream& input, bool useLabels);
	void readEdges(std::ifstream& input, bool useLabels);
	void readTriangles(std::ifstream& input, bool useLabels);
	void readQuadrilaterals(std::ifstream& input, bool useLabels);
	void readTetrahedra(std::ifstream& input, bool useLabels);
	void readHexahedra(std::ifstream& input, bool useLabels);

	void storeFaceLabel  (Face&   f, const int& lbl);
	void storeRegionLabel(Region& r, const int& lbl);

	/* a mesh */
	TMesh& m_mesh;

	/** mesh dimension */
	int m_mesh_dimension;

	/** labels of surfaces */
	std::map<int,std::string> m_2D_labels;

	/** labels of volumes */
	std::map<int,std::string> m_3D_labels;
};
/*----------------------------------------------------------------------------*/
