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
/** \file    GETMe.h
 *  \author  legoff
 *  \date    10/16/2015
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GETME_H_
#define GMDS_GETME_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IG.h"

#include "GMDS/Math/Hexahedron.h"
#include "GMDS/Math/Tetrahedron.h"
#include "GMDS/Math/Pyramid.h"
#include "GMDS/Math/Prism3.h"
#include "GMDS/Math/Quadrilateral.h"
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  GETMe
 *  \brief  Class that smoothes meshes according to the GETMe method
 */
/*----------------------------------------------------------------------------*/
class GETMe{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param AMesh the mesh where sheet operations are performed.
     */
	GETMe(IGMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     */
	virtual ~GETMe();

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.
 	 */
        void setVolumeTarget(bool ABool);

	/*------------------------------------------------------------------------*/
        /** \brief  Destructor.
         */
        void setVolumeTargetVariableName(std::string AVariableName);

    /*------------------------------------------------------------------------*/
    /** \brief  Smooth the mesh 
     *
     * \param AMarkFixedNodes a mark carried by nodes that should not move
     *
     */
	void exec(
			int ANbMaxIterSimult,
			int ANbMaxIterSeq,
			double AQualityThreshold,
			bool AAvoidTangledElements,
			int AMarkFixedNodes=-1);

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the minimum of the mean ratio of the mesh 
	 *
	 * \return the minimum of the mean ratio of the mesh
	 *
	 */
        double computeMinMeanRatio();

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the minimum of the scaled jacobian of the mesh 
         *
         * \return the minimum of the scaled jacobian of the mesh
         *
         */
        double computeMinNormalizedScaledJacobian();

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the minimum of the scaled jacobian of the mesh 
         *
         * \return the minimum of the scaled jacobian of the mesh
         *
         */
        double computeMinScaledJacobian2D();

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the max difference between cell volume and their expected 
	 *          volume.
         *
         * \return the max of the volume discrepancy
         *
         */
	double computeVolumeTargetDiscrepancy();

	/*------------------------------------------------------------------------*/
	/** \brief  Apply the simultaneous GETME smoothing
	*
	* \param AMarkFixedNodes a mark carried by nodes that should not move
	*
	*/
        void execSimult(
			int ANbMaxIter,
			double AQualityThreshold,
			double ARelaxFactor,
			double AQualityDependantFactorMin,
	                double AQualityDependantFactorMax,
        	        double AWeightExponent,
                	bool AIsGeomAssociated,
	                bool AIsGeomAssocForced,
			bool AAvoidTangledElements,
                        int AMarkFixedNodes=-1);	

       /*------------------------------------------------------------------------*/
       /** \brief  Apply the simultaneous GETME smoothing
	*
	* \param AMarkFixedNodes a mark carried by nodes that should not move
	*
	*/
        void execSimult2D(
			int ANbMaxIter,
			double AQualityThreshold,
			double ARelaxFactor,
			double AQualityDependantFactorMin,
	                double AQualityDependantFactorMax,
        	        double AWeightExponent,
                	bool AIsGeomAssociated,
	                bool AIsGeomAssocForced,
			bool AAvoidTangledElements,
                        int AMarkFixedNodes=-1);

	/*------------------------------------------------------------------------*/
        /** \brief  Apply the sequential GETME smoothing
        *
        * \param AMarkFixedNodes a mark carried by nodes that should not move
        *
        */
        void execSeq(
                        int ANbMaxIter,
                        double AQualityThreshold,
                        double ARelaxFactor,
                        double AQualityDependantFactor,
			double ADeltaInvalid,
			double ADeltaReselect,
			double ADeltaSuccess,
                        bool AIsGeomAssociated,
                        bool AIsGeomAssocForced,
                        bool AAvoidTangledElements,
                        int AMarkFixedNodes=-1);

protected:

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the new position of the vertices of a cell according to the GETMe smoothing scheme
         *
	 * \param ARegion the region on which we work on
         * \param AScalingFactor the scaling factor applied
	 * \param APoints the computed positions
         * 
         */
	void computeGETMePoint(
			const gmds::Cell* ACell, 
			double AScalingFactor, 
			std::vector<gmds::math::Point>& APoints) const;

	void computeGETMePoint(
                        const gmds::math::Hexahedron AHex,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const;

	void computeGETMePoint(
                        const gmds::math::Tetrahedron ATet,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const;

	void computeGETMePoint(
                        const gmds::math::Pyramid APyr,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const;

	void computeGETMePoint(
                        const gmds::math::Prism3 APrism,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const;

	void computeGETMePoint(
                        const gmds::math::Quadrilateral AQuad,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Apply a size modification factor to the new position of the vertices of a cell 
         *
	 * \param ACell the cell on which we work on
	 * \param APoints the computed positions
	 * \param AModifySizeFactor the scaling factor applied
         * 
         */
	void applyScaling(
			  const gmds::Cell* ACell,
			  std::vector<gmds::math::Point>& APoints,
			  double AModifySizeFactor=1.) const;


	/* a mesh */
	IGMesh& m_mesh;

	/* indicates whether a target volume should be aimed at for each cell */
	bool m_isVolumeTargetOn;

	/* name of the cell variable that carries the target volume */
	std::string m_volumeTargetVariableName;

	/* stopping threshold regarding volume target */
	double m_volumeTargetAccuracy;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GETME_H_ */
/*----------------------------------------------------------------------------*/
