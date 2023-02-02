/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    OrderedSmartLaplacianSmoothing.h
 *  \author  legoff
 *  \date    06/30/2015
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ORDEREDSMARTLAPLACIANSMOOTHING_H_
#define GMDS_ORDEREDSMARTLAPLACIANSMOOTHING_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// CaGe File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/CAD/GeomManager.h"
#include "GeomMeshIntersectionService.h"

#include "Smoothing.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class OrderedSmartLaplacianSmoothing
 *  \brief This is a dummy class.
 */
/*----------------------------------------------------------------------------*/
class OrderedSmartLaplacianSmoothing: public Smoothing {

public:


    /*------------------------------------------------------------------------*/
    /** \brief  Constructor
     */
	OrderedSmartLaplacianSmoothing(
			gmds::IGMesh& AMesh,
			gmds::geom::GeomManager& AManager,
			gmds::geom::GeomMeshIntersectionService& AService);

    /*------------------------------------------------------------------------*/
    /** \brief Copy constructor
     */
	OrderedSmartLaplacianSmoothing(const OrderedSmartLaplacianSmoothing&);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor
     */
	~OrderedSmartLaplacianSmoothing();

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator=
     */
	OrderedSmartLaplacianSmoothing& operator=(const OrderedSmartLaplacianSmoothing&);

    /*------------------------------------------------------------------------*/
    /** \brief  Execute the smoothing algorithm one time.
     */
	virtual void exec(int AMarkFixedNodes=-1);

    /*------------------------------------------------------------------------*/
    /** \brief  Execute the smoothing algorithm one time.
     */
	virtual void initNodesAdjacencies();

private:

	bool checkMeshModel() const;

	struct NodeAndQuality {
		TCellID cellID_;
		double quality_;
	};

	struct nodeQualityCompare{
		bool operator()(const NodeAndQuality& AFirstCell, const NodeAndQuality& ASecondCell) const {
			return AFirstCell.quality_ < ASecondCell.quality_;
		}
	};

	

	void initOrderedNodeList(
			std::map<gmds::TCellID,double>& ARegionsQuality,
			std::multiset<NodeAndQuality,nodeQualityCompare>& AOrderedNodes);

	bool isQualityImproved(gmds::Node& ANode, gmds::math::Point& ANewPosition);

	/* a mesh */
	IGMesh& mesh_;

	/* a geometric model */
	gmds::geom::GeomManager& manager_;

	/* service associated to the geometric model */
	gmds::geom::GeomMeshIntersectionService& service_;

	std::map<gmds::Node,std::vector<gmds::Node> > node2nodes_;

	gmds::Variable<double>* nodeQuality_;
	gmds::Variable<double>* regionQuality_;

};
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_ORDEREDSMARTLAPLACIANSMOOTHING_H_ */
