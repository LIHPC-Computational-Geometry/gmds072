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
/** \file    BoundaryOperator.h
 *  \author  F. LEDOUX/
 *  \date    12/18/2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_BOUNDARY_OPERATOR_H_
#define GMDS_BOUNDARY_OPERATOR_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{
    
    /*----------------------------------------------------------------------------*/
    /** \class  BoundaryOperator
     *  \brief  Class gathering operations used to get and mark the cells that
     belongs to the boundary of a 2D, 3D mesh
     */
    class EXPORT_GMDS BoundaryOperator
    {
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         *  \param AMesh the mesh where sheet operations are performed.
         */
        BoundaryOperator(IGMesh* AMesh);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.	*/
        virtual ~BoundaryOperator();
        
        
        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Mark all the boundary cells.
         *			Marks are created inside the method, but must be free by the
         *          caller. This function also initialize and fill in the variable
         *          BND_SURFACE_COLOR, which is assigned to faces.
         *
         * 	\param AMarkFOnSurf faces on boundary surfaces are marked with it
         * 	\param AMarkEOnSurf edges on boundary surfaces are marked with it
         * 	\param AMarkNOnSurf nodes on boundary surfaces are marked with it
         * 	\param AMarkEOnCurv edges on boundary curves are marked with it
         * 	\param AMarkNOnCurv nodes on boundary curves are marked with it
         * 	\param AMarkNOnPnt  nodes on geometric point are marked with it
         * 	\param AMarkAN isolated nodes (connected to nothing) are marked with it
         */
        
        void markCellOnGeometry(const int AMarkFOnSurf,
                                const int AMarkEOnSurf,
                                const int AMarkNOnSurf,
                                const int AMarkEOnCurve,
                                const int AMarkNOnCurve,
                                const int AMarkNOnPnt,
                                const int AMarkAN);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Mark cells on the boundary surface
         *
         * 	\param AMarkBF faces on the boundary surfs are marked with it
         * 	\param AMarkBE edges on the boundary surfs are marked with it
         * 	\param AMarkBE edges on the boundary surfs are marked with it
         */
        
        void markCellsOnSurfaces(const int AMarkBF, const int AMarkBE, const int AMarkBN);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Mark the boundary edges that fit geometric curves
         * \param  AMarkBF IN faces on geom surf
         * \param  AMarkBE IN edges on geom surf
         * \param  AMarkCE OUT edges on geom curves
         * \param  AMarkCN OUT nodes on geom curves
         */
        void markCellsOnCurves(const int AMarkBF , const int AMarkBE,
                               const int AMarkCE, const int AMarkCN);
        
        /*------------------------------------------------------------------------*/
        /** \brief Mark the boundary edges that fit geometric curves for surface
         *         mesh
         * \param  AMarkCE OUT edges on geom curves
         * \param  AMarkCN OUT nodes on geom curves
         */
        void markCellsOnCurves(const int AMarkCE, const int AMarkCN);
        
        
        /*------------------------------------------------------------------------*/
        /** \brief  Mark the boundary nodes that fit geometric vertices
         * \param  AMarkCE IN  edges on geom curves
         * \param  AMarkCN IN  nodes on geom curves
         * \param  AMarkPN OUT nodes on geom points
         */
        void markNodesOnPoint(const int AMarkCE, const int AMarkCN,
                              const int AMarkPN);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Mark the boundary nodes that do not have adjacent edges
         * \param  AMarkAlone IN mark of alone nodes
         */
        void markAloneNodes(const int AMarkAlone);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Color the boundary faces with different colors (using a mesh
         variable called "BND_SURFACE_COLOR") to get one color per geometric surface
         *
         * \param[in] AMarkFOnSurf all the boundary faces are marked with it
         * \param[in] AMarkEOnCurv all the edges classified on curves are marked with 
         *            it
         */
        void colorFaces(const int AMarkFOnSurf, const int AMarkEOnCurv);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Color the boundary edges with different colors (using a mesh
         variable called "BND_SURFACE_COLOR") to get one color per geometric surface
         *
         * \param[in] AMarkEOnCurv all the boundary edges are marked with it
         * \param[in] AMarkNOnPnt  all the nodes classified on points are marked
         *            with it
         */
        void colorEdges(const int AMarkEOnCurv, const int AMarkNOnPnt);
        
        /*------------------------------------------------------------------------*/
        /** \brief  Color the boundary boundary nodes (using a mesh
         variable called "BND_VERTEX_COLOR") to get one color per geometric vertex
         *
         * \param[in] AMarkNOnPnt  all the nodes classified on points are marked
         *            with it
         */
        void colorNodes(const int AMarkNOnPnt);
        
        /*------------------------------------------------------------------------*/
        /** \brief Return the boundary nodes
         *
         * \param  ANodeIDs the ids of boundary nodes
         */
        void getBoundaryNodes(std::vector<TCellID>& ANodeIDs);
        
        /*------------------------------------------------------------------------*/
        /** \brief Compute the normal to AFace going out of ARegion
         * \param  AFace a face
         * \param  ARegion a region adj. to AFace
         * \return the normal to AFace going out of ARegion
         */
        math::Vector getOutputNormal(Face& AFace, Region& ARegion);
        
        /*------------------------------------------------------------------------*/
        /** \brief Compute the normal to AFace going out of a mesh
         * \param  AFace a face
         * \return the normal to AFace going out
         */
        math::Vector getOutputNormalOfABoundaryFace(Face& AFace);
        /*------------------------------------------------------------------------*/
        /** \brief Compute the normal to ANode as the average of the normal to the
         *         adjacent boundary faces (weighted by the face area)
         * \param  ANode a face
         * \return the normal to AFace going out
         */
        math::Vector getOutputNormalOfABoundaryNode(Node& ANode);
        
    protected:
        /* a mesh */
        IGMesh* m_mesh;
        
        Variable<int>* m_var_color_surf;
        Variable<int>* m_var_color_curve;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_BOUNDARY_OPERATOR_H_ */
/*----------------------------------------------------------------------------*/
