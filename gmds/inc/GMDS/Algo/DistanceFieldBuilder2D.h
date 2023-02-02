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
/** \file    DistanceFieldBuilder2D.h
 *  \author  F. LEDOUX
 *  \date    12/17/2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_DISTANCE_FIELD_BUILDER_2D_H_
#define GMDS_DISTANCE_FIELD_BUILDER_2D_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Cross2D.h>
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    /** \class DistanceFieldBuilder2D
     *  \brief Compute the distance field to a set of nodes for a mesh. The
     obtained distance is stored in a mesh variable. This operation can be used for
     2D Meshes. In the latter case, the algorithm only work on the mesh
     surface
     */
    /*----------------------------------------------------------------------------*/
    
    class EXPORT_GMDS DistanceFieldBuilder2D
    {
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         *  \param AMesh the mesh we work on
         */
        DistanceFieldBuilder2D(IGMesh* AMesh = 0);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.	*/
        virtual ~DistanceFieldBuilder2D();
        
        /*------------------------------------------------------------------------*/
        /** \brief  Check if we can apply the algorithm on the provided mesh
         *		(dimension, available connectivities)
         */
        bool isValid();
        
        /*------------------------------------------------------------------------*/
        /** \brief this operation constructs the distance for each node in AToCompute
         *	   and the nodes in AFrom. In thi
         *
         *	\param AFrom      a set of nodes we want the distance to
         *  \param AToCompute the nodes we want assigned a distance value
         *	\param AMarkFace  faces we can work with are marked
         */
        Variable<double>* computeDistance(std::vector<Node>& AFrom,
                                          std::vector<Node>& AToCompute,
                                          const int AMarkFace);
        
        
        /*------------------------------------------------------------------------*/
        /** \brief this operation constructs the distance for each node in AToCompute
         *	   and the nodes in AFrom. In thi
         *
         * \param AFrom      a set of nodes we want the distance to
         * \param AToCompute the nodes we want assigned a distance value
         * \param AMarkFace  faces we can work with are marked
         */
        Variable<double>*
        computeLInfDistanceFromCrossField(std::vector<Node>& AFrom,
                                          std::vector<Node>& AToCompute,
                                          Variable<math::Cross2D>* AField,
                                          const int AMarkFace);
        
        /*------------------------------------------------------------------------*/
        /** \brief this operation provides a vector giving the order of the node
         *		   insertion process during the distance computation process.
         */
        std::vector<TCellID> getInsertionOrder() const;
        
    protected:
        
        /*------------------------------------------------------------------------*/
        /** \brief this operation constructs the distance for each node in AToCompute
         *	   and the nodes in AFrom. In thi
         *
         *	\param AFrom      a set of nodes we want the distance to
         *  \param AToCompute the nodes we want assigned a distance value
         *	\param AMarkFace  faces we can work with are marked
         */
        void internalComputeDistance(std::vector<Node>& AFrom,
                                     std::vector<Node>& AToCompute,
                                     const int AMarkFace);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Computes the set of nodes sharing a face with ANode and being
         *          marked with AMark
         *
         *  \param  ANode the nodes we start from
         *  \param  AMarkNode resulting nodes must be mark with it
         \param	AMarkFace they also have to belong to faces marked AMarkFace
         *  \return the nodes both adjacent to ANode and marked with AMark
         */
        std::vector<Node> getAdjacentNodes(Node& ANode, const int AMarkNode,
                                           const int AMarkFace);
        
        /*------------------------------------------------------------------------*/
        /** \brief Computes the distance at AN using values defined in the other
         *		   nodes  AN1 and AN2
         *         The distance value in AN is returned
         */
        double extrapolateDistance(Node& AN, Node& AN1, Node& AN2,
                                   Variable<double>*& ADist,
                                   const bool AOnlyOne=false);
        
        /*------------------------------------------------------------------------*/
        /** \brief Computes the distance at AN using values defined in the other
         *		   nodes  AN1 and AN2
         *         The distance value in AN is returned
         */
        double extrapolateLInfDistance(Node& AN, Node& AN1, Node& AN2,
                                       Variable<double>*& ADist,
                                       const bool AOnlyOne=false);

        
        /*------------------------------------------------------------------------*/
        /** \brief Initializes the set of nodes that will make the first layer for
         *		   the advancing front algorithm
         *	\param AFrom the nodes we start from (distance 0 so)
         *	\param ANarrowBand the set of nodes that will be put in the narrow band,
         *		   i.e. the set of candidates or the active front
         *	\param AExtrapolateDistance distance computed for each node (id) in the
         *	       narrow band
         *	\param AMarkNodeToWorkOn the nodes we will consider only
         *	\param AMarkNodeAlive the nodes where a distance is already defined
         *	\param AMarkNodeNarrow the nodes that are candidates to be next alive node
         *	\param AMarkFace faces we can use nodes
         */
        void initNarrowBand(std::vector<Node>& AFrom,
                            std::multimap<double, Node > &ANarrowBand,
                            std::map<TCellID, double>& AExtrapolateDistance,
                            const int AMarkNodeToWorkOn,
                            const int AMarkNodeAlive,
                            const int AMarkNodeNarrow,
                            const int AMarkFace);
        
        /*------------------------------------------------------------------------*/
        /** \brief Advancing front algorithm process
         *	\param ANarrowBand the set of nodes that are candidates to advance, or
         the active front
         *	\param AExtrapolateDistance distance computed for each node (id) in the
         *	       narrow band
         *	\param AMarkNodeToWorkOn the nodes we will consider only
         *	\param AMarkNodeAlive the nodes where a distance is already defined
         *	\param AMarkNodeNarrow the nodes that are candidates to be next alive node
         *	\param AMarkFace faces we can use nodes
         */
        void advanceDistanceFront(std::multimap<double,Node>& ANarrowBand,
                                  std::map<TCellID, double>& AExtrapolateDistance,
                                  const int AMarkNodeToWorkOn,
                                  const int AMarkNodeAlive,
                                  const int AMarkNodeNarrow,
                                  const int AMarkFace);
        /*------------------------------------------------------------------------*/
        /** \brief Check if 3 nodes belongs to the same face
         *	\param AN1 a node
         *	\param AN2 a node
         *	\param AN3 a node
         *	\param AF  a face
         *	\return true if AN1, AN2 and AN3 belongs to AF, false otherwise
         */
        bool belongToTheSameFace(Node& AN1, Node& AN2, Node& AN3, Face& AF);
        
        /*------------------------------------------------------------------------*/
        /** \brief Among a set of nodes (AIN), we get the pair of nodes that share
         *		   a common face in AOUT
         */
        void getNodesInTheSameFace(const std::vector<Node>& AIN,
                                   std::vector<std::pair<TCellID, TCellID> >& AOUT);
    private:
        //the mesh we work on
        IGMesh* m_mesh;
        Variable<double>* m_distance;
        Variable<math::Cross2D>* m_cross_field;
        Variable<math::Vector>* m_adv_field[2];
        Variable<double>* m_adv_dist[2];
        Variable<int>* m_adv_order;
        std::vector<TCellID> m_insertion_order;
        
        
        typedef double (DistanceFieldBuilder2D::*DistFuncPtr)
        (Node& AN, Node& AN1, Node& AN2, Variable<double>*& ADist, const bool);
        
        DistFuncPtr m_extrapolate_func_ptr;
        
    };
    
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_DISTANCE_FIELD_BUILDER_2D_H_ */
/*----------------------------------------------------------------------------*/

