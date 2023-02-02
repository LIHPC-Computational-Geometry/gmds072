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
/** \file    SubMapping.h
 *  \author  legoff
 *  \date    01/02/2016
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SUBMAPPING_SUBMAPPING_H_
#define GMDS_SUBMAPPING_SUBMAPPING_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/GeomCurve.h>
#include <GMDS/CAD/GeomPoint.h>
#include <GMDS/CAD/GeomSurface.h>
#include <GMDS/SubMapping/SubMappingCommon.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace submapping {
/*----------------------------------------------------------------------------*/
/** \class SubMapping
 *  \brief This is a dummy class.
 */
/*----------------------------------------------------------------------------*/
class SubMapping
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        SubMapping();

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor
         */
        SubMapping(const SubMapping&);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        ~SubMapping();

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator=
         */
        SubMapping& operator=(const SubMapping&);

        /*------------------------------------------------------------------------*/
        /** \brief  Init
         */
        void init(gmds::geom::GeomVolume* AVol);

        /*------------------------------------------------------------------------*/
        /** \brief  Init
         */
        void init(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Run the submapping algorithm on a Volume (its surfaces)
         */
        void exec(gmds::geom::GeomVolume* AVol);

        /*------------------------------------------------------------------------*/
        /** \brief  Run the submapping algorithm on a surface
         */
        void exec(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Print vertices classification for a volume
         */
        void printVerticesClassification(gmds::geom::GeomVolume* AVol);

        /*------------------------------------------------------------------------*/
        /** \brief  Print curves classification for a volume
         */
        void printCurvesClassification(gmds::geom::GeomVolume* AVol);

        /*------------------------------------------------------------------------*/
        /** \brief  Print vertices classification for a surface
         */
        void printVerticesClassification(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Print curves classification for a surface
         */
        void printCurvesClassification(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Print curves discretization for a surface
         */
        void printCurvesDiscretization(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Print curves discretization for a surface
         */
        void printSubmappingInfo(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Print curves discretization for a surface
         */
        void exportVTKSubmappingInfo(gmds::geom::GeomVolume* AVol, const std::string& AFile);

 private:
        /*------------------------------------------------------------------------*/
        /** \brief  Classify the vertices
         */
        void verticesClassification(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Correct the vertices classification
         */
        void correctVerticesClassification(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief
         */
        void computeRhoAndOmega(gmds::geom::GeomSurface* ASurf, double* rho, double* omega);

        /*------------------------------------------------------------------------*/
        /** \brief  Classify the curves
         */
        void curvesClassification(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Assign a target discretization to each curve (depending on size)
         */
        void intervalAssignment(gmds::geom::GeomVolume* AVol, double ASize = 1.);

        /*------------------------------------------------------------------------*/
        /** \brief  Assign a target discretization to each curve (depending on size)
         */
        void intervalAssignment(gmds::geom::GeomSurface* ASurf, double ASize = 1.);

        /*------------------------------------------------------------------------*/
        /** \brief  Solve the integer problem of curves discretization
         */
        void boundaryDiscretization(gmds::geom::GeomVolume* AVol);

        /*------------------------------------------------------------------------*/
        /** \brief  Solve the integer problem of curves discretization
         */
        void boundaryDiscretization(gmds::geom::GeomSurface* ASurf);

        /*------------------------------------------------------------------------*/
        /** \brief  Solve the integer problem of curves discretization
         */
        void boundaryDiscretization(std::vector<gmds::geom::GeomSurface*> ASurfaces,
                                    std::vector<gmds::geom::GeomCurve*> ACurves);

        /*------------------------------------------------------------------------*/
        /** \brief  Get a vertex classification (seen from a surface)
         */
        int getClassification(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint);

        /*------------------------------------------------------------------------*/
        /** \brief  Set a vertex classification (seen from a surface)
         */
        void setClassification(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint, const int AInt);

        /*------------------------------------------------------------------------*/
        /** \brief  Get a vertex angle (seen from a surface)
         */
        double getAngle(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint);

        std::map<gmds::geom::GeomSurface*, std::vector<gmds::geom::GeomPoint*> > m_orderedDirectVertices;
        std::map<gmds::geom::GeomSurface*, std::vector<gmds::geom::GeomCurve*> > m_orderedDirectCurves;

        std::map<gmds::geom::GeomSurface*, std::map<gmds::geom::GeomCurve*, bool> > m_isOrderedDirectCurves;

        std::map<gmds::geom::GeomSurface*, std::map<gmds::geom::GeomPoint*, double> > m_angleAtVertexFromSurface;

        std::map<gmds::geom::GeomSurface*, std::map<gmds::geom::GeomPoint*, int> > m_verticesClassification;
        std::map<gmds::geom::GeomSurface*, std::map<gmds::geom::GeomCurve*, gmds::submapping::CurveClassification> >
            m_curvesClassification;

        std::map<gmds::geom::GeomCurve*, int> m_targetDiscretization;
        std::map<gmds::geom::GeomCurve*, int> m_computedDiscretization;
};
/*----------------------------------------------------------------------------*/
}  // end namespace submapping
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_SUBMAPPING__SUBMAPPING_H_ */
