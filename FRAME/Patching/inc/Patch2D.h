/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
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
 * Patch2D.h
 *
 *  Created on: sept. 10, 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCH_2D_H_
#define PATCH_2D_H_
/*----------------------------------------------------------------------------*/
// STL Headers
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
// FRAME Headers
#include "PatchEdge2D.h"
#include "PatchVertex2D.h"
/*----------------------------------------------------------------------------*/
// Forward declarations to solve compilation loops.
class PatchComplex2D;
/*----------------------------------------------------------------------------*/
/* \class This class describes what is a single 2D Patch
 */
class EXPORT_GMDS Patch2D
{
public:
    
    struct Grid{
        /*------------------------------------------------------------------*/
        /* \brief Constructor frome size grid. As a result an AXxAY grid is 
         *        initialized
         */
        Grid(const int AX, const int AY);
        /*------------------------------------------------------------------*/
        /* \brief Default constructor
         */
        Grid();
        /*------------------------------------------------------------------*/
        /* \brief Copy constructor
         */
        Grid(const Grid& AG);
        
        int X;
        int Y;
        std::vector<std::vector< gmds::math::Point> > points;
    };
    
    /*------------------------------------------------------------------------*/
    /* \brief Default constructor for STL container compliance
     */
    Patch2D();
    /*------------------------------------------------------------------------*/
    /* \brief Constructor from a list of edges
     *
     * \param AOwner      the complex the patch belongs to
     * \param A vector of edges forming a closed path
     */
    Patch2D(PatchComplex2D* AOwner,
            std::vector<PatchEdge2D*>& AEdges);
    /*------------------------------------------------------------------------*/
    /* \brief Copy constructor
     */
    Patch2D(const Patch2D& AP);
    
    /*------------------------------------------------------------------------*/
    /* \brief Destructor
     */
    virtual ~Patch2D();
    
    /*------------------------------------------------------------------------*/
    /* \brief build the inner grid of the patch
     *
     * \return true if the grid is built
     */
    bool buildInnerGrid();
    
    /*------------------------------------------------------------------------*/
    /* \brief build the boundary lines of the grid
     *
     * \return true if the grid is built
     */
    bool buildBoundaryGrid();
    
    /*------------------------------------------------------------------------*/
    /* \brief Get access to the patch grid
     *
     * \param[OUT] the corners points
     */
    Grid& grid();
    
    /*------------------------------------------------------------------------*/
    /* \brief Return the number of mesh points expected along X direction
     */
    int X() const;
    
    /*------------------------------------------------------------------------*/
    /* \brief Return the number of mesh points expected along Y direction
     */
    int Y() const;
private:
    
    /** complex in which (*this) is buit */
    PatchComplex2D* m_complex;

    /** nb mesh point in X direction */
    int m_X;
    /** nb mesh point in X direction */
    int m_Y;
    /** internal grid of points */
    Grid m_grid;
    /* boolean indicating if the grid is filled in */
    bool m_has_no_grid;

    /** boundary patch edges */
    std::vector<PatchEdge2D*> m_edges;
    
};
/*----------------------------------------------------------------------------*/
#endif /* defined(__GMDSSuite__Patch2D__) */
/*----------------------------------------------------------------------------*/
