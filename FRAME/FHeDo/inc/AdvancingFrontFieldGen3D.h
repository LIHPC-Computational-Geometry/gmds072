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
#ifndef ADVANCING_FRONT_FIELD_GEN_3D_H_
#define ADVANCING_FRONT_FIELD_GEN_3D_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <map>
#include <list>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
#include <GMDS/Algo/DistanceFieldBuilder3D.h>
/*----------------------------------------------------------------------------*/
#include "AdvancingFrontFieldCommon.h"
/*----------------------------------------------------------------------------*/
/** \class  AdvancingFrontFieldGen3D
*  \brief  Computes the frame field onto a 3D mesh 
*/
class EXPORT_GMDS AdvancingFrontFieldGen3D : public AdvancingFrontFieldCommon
{
public:
	/*------------------------------------------------------------------------*/
	/** \brief Constructor.
	*
	*  \param AMesh the mesh where we work on
	*/
	AdvancingFrontFieldGen3D(gmds::IGMesh* AMesh);


   /*------------------------------------------------------------------------*/
   /** \brief Function to be called for executing a frame field generation 
    *		   algorithm. This function follows the template method DP. It
    *		   means that it drives the global execution algorithm, each 
    *		   specific part being delegated to child classes.
    */
   virtual void execute();
	/*------------------------------------------------------------------------*/
	/** \brief Two types of methods can be used to compute the distance field
	*		   The first consists in building a 2D distance field on the 
	*		   surface. Each surface node will have its distance to the closest
	*		   curve. And a 3D distance field, where each inner node will store
	*		   the distance to the closest surface, curve or point.
	*		   The second consists in building for all the nodes their distance
	*		   to the closest geometric curve.
	*
	*  \param AMesh the mesh where we work on
	*/
	enum EDistanceFieldMethod {
		EDistanceField_DisconnectedSurfAndVol,
		EDistanceField_VolFromEdge,
	} ;
	/*------------------------------------------------------------------------*/
	/** \brief Set the way we compute the distance field
	*
	*  \param AMesh the mesh where we work on
	*/
	void setDistanceFieldMethod(const EDistanceFieldMethod& AMethod);

protected:
    
    /*----------------------------------------------------------------------------*/
    /** \struct StabilityInfo
     *  \brief  Stores a relation between a node and the stability value, 
     *          which is computed for it. It also provices a relation order useful
     *          for some STL containers.
     */

    struct StabilityInfo
    {
        /// the node id we interest on
        gmds::TCellID  nodeID;
        /// stability value
        double stability;
        /*-----------------------------------------------------------------*/
        /** \brief Constructor.
         *  \param AN a node
         *  \paral AV the associated value
         */
        StabilityInfo(const gmds::Node& AN, const double AV=0)
        :nodeID(AN.getID()), stability(AV){}
        /*-----------------------------------------------------------------*/
        /** \brief Constructor.
         *  \param AI a node id
         *  \paral AV the associated value
         */
        StabilityInfo(const gmds::TCellID AI=gmds::NullID, const double AV=0)
        :nodeID(AI), stability(AV){}
        
        /*-----------------------------------------------------------------*/
        /** \brief Comparison operator < betwee two stability infos
         */
        friend bool
        operator<(const StabilityInfo& AS1, const StabilityInfo& AS2)
        {
            if(AS1.stability < AS2.stability)
                return true;
            if(AS1.stability == AS2.stability)
                if(AS1.nodeID < AS2.nodeID)
                    return true;
            return false;
        }
        
        /*-----------------------------------------------------------------*/
        /** \brief Comparison operator  == betwee two stability infos
         */
        friend bool
        operator==(const StabilityInfo& AS1, const StabilityInfo& AS2)
        {
            return(AS1.nodeID == AS2.nodeID);
        }
    };
    
    /*------------------------------------------------------------------------*/
    /** \brief extracts the most stable candidate from AL. In practice, it is the
     *         with the smalles value for stability.
     */
    StabilityInfo extractElectFrom(std::list<StabilityInfo>& AL);
    

	/*----------------------------------------------------------------------*/
	double distancePonderation(const double AVal, const double ADist);
	void initDistanceBallsOnSurf();
	virtual void initDistanceBallsInVol();
	void putInBalls(gmds::Node& AN1, gmds::Node& AN2);
	virtual void computeDistanceFields();
	virtual void initQuaternionsBySnapping();

	void buildFrameField();

	void  initNarrowBand(std::list<StabilityInfo>& ANarrowBand);

	void initNarrowBandSurf(std::list<StabilityInfo>& ANarrowBand);

	void propagateQuaternions(std::list<StabilityInfo>& ANarrowBand);

	void propagateQuaternionsSurf(std::list<StabilityInfo>& ANarrowBand);

	virtual void smoothAll();
    
    virtual void colorSimplices();

	std::vector<gmds::Node> getAdjacentNodesByRegion(gmds::Node& ANode);

	void getNodesInTheSameTetrahedron(
		const std::vector<gmds::Node>& AIN,
		std::vector<std::vector<gmds::Node> >& AOUT);
	
private:


	void surfaceSmoothing();
	void smoothNodeAndBall(gmds::Node& ANode);
	void internalSmooth(gmds::Node& ANode);

	/*------------------------------------------------------------------------*/
	/** \brief Indicates if we must update the stability info for the current
			   ball. 
	*
	*  \param ANarrowBandSize the actual size of the narrow band. 
	*
	*/
	bool mustUpdateStabilityInfo(const int ANarrowBandSize);
	/*------------------------------------------------------------------------*/
	/** \brief Returns the distance computed in ANode considering the distance
	*		   field strategy (surf+vol or global to curves)
	*
	*  \param ANode the node we want to compute the distance on
	*  \return a distance
	*
	*/
	double getDistance(const gmds::Node& Anode);

	void smooth(const int AMark);

	// Compute the frame in ANode using the frames of adjacent nodes, marked with AMark
	// ANode must be marked at the end of this step
	bool computeLocalFrame(gmds::Node& ANode, const int& AMark);


	/* Using the distance field, all the nodes that belongs to the "medial axis" are
	* computed and returned
	*/
	std::vector<gmds::Node> computeApproximateMedialAxis(const std::vector<gmds::Node>& ANodes);



	/*---------------------------------------------------------------------------*/
	bool belongToTheSameRegion(
		gmds::Node& ATo,
		gmds::Node& AFrom1,
		gmds::Node& AFrom2,
		gmds::Node& AFrom3,
		gmds::Region& ARegion);


	void extrapolateQuaternionInTet(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);


private:

	EDistanceFieldMethod m_distance_method;
	gmds::Variable<double>* m_distance_field_2D; //dist. field on the surface
	gmds::Variable<double>* m_distance_field_3D; //dist. field inside the volume

	double m_max_distance_2D;
	double m_max_distance_3D;
	double m_max_distance_2D_div2;
	double m_max_distance_3D_div2;
	//gmds::Variable<double>*  m_medial_axis_distance_field; //dist_MA

	gmds::DistanceFieldBuilder2D m_distanceBuilder2D;
	gmds::DistanceFieldBuilder3D m_distanceBuilder3D;

	std::vector<gmds::Node> m_inner_nodes; //all the nodes in the volume

	std::vector<gmds::TCellID> m_nodes_by_creation_order;
};
/*----------------------------------------------------------------------------*/
#endif //ADVANCING_FRONT_FIELD_GEN_3D_H_
/*----------------------------------------------------------------------------*/
