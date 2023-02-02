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
/** \file    SubMapping.cpp
 *  \author  legoff
 *  \date    01/02/2016
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <vector>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <glpk.h>
/*----------------------------------------------------------------------------*/
// Project File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/SubMapping/SubMapping.h>

#include <GMDS/CAD/FacetedCurve.h>
#include <GMDS/CAD/GeomCurve.h>
#include <GMDS/IG/IG.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
#include <GMDS/SubMapping/SubMappingCommon.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace submapping {
/*----------------------------------------------------------------------------*/
SubMapping::SubMapping()
{
}
/*----------------------------------------------------------------------------*/
SubMapping::SubMapping(const SubMapping& ASubMapping)
{
}
/*----------------------------------------------------------------------------*/
SubMapping::~SubMapping()
{
}
/*----------------------------------------------------------------------------*/
SubMapping&
SubMapping::operator=(const SubMapping& ASubMapping)
{
        return *this;
}
/*----------------------------------------------------------------------------*/
void
SubMapping::init(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface*> surfaces;
        AVol->get(surfaces);

        for (int iSurf = 0; iSurf < surfaces.size(); iSurf++) {
                init(surfaces[iSurf]);
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::init(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::init ----" << std::endl;

        std::vector<gmds::geom::GeomPoint*> vertices;
        std::vector<gmds::geom::GeomCurve*> curves;

        ASurf->get(vertices);
        ASurf->get(curves);

        // a few checks are performed in order to ensure that we are in a "classic" configuration
        if (curves.empty()) {
                throw GMDSException("SubMapping::init cannot work when there is no curve.");
        }

        if (vertices.size() == 3) {
                throw GMDSException("SubMapping::init cannot work when there are only 3 vertices on a surface.");
        }

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                std::vector<gmds::geom::GeomSurface*> surfaces;
                curves[iCurve]->get(surfaces);

                // we will first check that no curve has only one surface (which should be ASurf);
                // it can happen for example in the cylinder case and the truncated sphere case.
                // It may be possible to handle the case but it is a bit tricky
                if (surfaces.size() == 1) {
                        throw GMDSException(
                            "SubMapping::init a curve is adjacent to only one surface. SubMapping currently does not "
                            "handle that.");
                }

                // check that the curve is not a loop. In this case it could be a good idea to use a predifined template
                // introducing "virtual vertices"
                if (curves[iCurve]->isALoop()) {
                        throw GMDSException(
                            "SubMapping::init a curve is a loop. SubMapping currently does not handle that.");
                }
        }

        // now we order the curves and the vertices
        std::cout << "---- SubMapping::init ordering curves and vertices ----" << std::endl;
        std::vector<gmds::geom::GeomCurve*> orderedCurves;
        std::vector<gmds::geom::GeomPoint*> orderedVertices;
        std::map<gmds::geom::GeomCurve*, bool> isOrderedDirectCurves;
        {
                gmds::geom::GeomCurve* current_curve = curves[0];
                gmds::geom::GeomPoint* current_vertex = current_curve->getFirstPoint();

                do {
                        // std::cout<<current_vertex->getPoint()<<std::endl;
                        // std::cout<<current_curve<<std::endl;
                        orderedCurves.push_back(current_curve);
                        orderedVertices.push_back(current_vertex);

                        gmds::geom::GeomPoint* next_vertex;

                        if (current_vertex == current_curve->getFirstPoint()) {
                                next_vertex = current_curve->getSecondPoint();
                                isOrderedDirectCurves[current_curve] = true;
                        } else {
                                next_vertex = current_curve->getFirstPoint();
                                isOrderedDirectCurves[current_curve] = false;
                        }

                        // identify the next curve
                        bool next_curve_found = false;
                        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                                gmds::geom::GeomPoint* firstVertex = curves[iCurve]->getFirstPoint();
                                gmds::geom::GeomPoint* secondVertex = curves[iCurve]->getSecondPoint();

                                if (firstVertex == next_vertex && secondVertex != current_vertex) {
                                        current_curve = curves[iCurve];
                                        current_vertex = next_vertex;
                                        next_curve_found = true;
                                        break;
                                } else {
                                        if (secondVertex == next_vertex && firstVertex != current_vertex) {
                                                current_curve = curves[iCurve];
                                                current_vertex = next_vertex;
                                                next_curve_found = true;
                                                break;
                                        }
                                }
                        }
                        if (not next_curve_found) {
                                throw GMDSException("SubMapping::init next curve not found.");
                        }

                } while (current_curve != curves[0]);

                // check post ordering
                if (orderedCurves.size() != curves.size()) {
                        throw GMDSException(
                            "SubMapping::init orderedCurves.size() != curves.size(). Most likely there is a hole, "
                            "which is not currently handled.");
                }
                if (orderedVertices.size() != vertices.size()) {
                        throw GMDSException(
                            "SubMapping::init orderedVertices.size() != vertices.size(). Most likely there is a hole, "
                            "which is not currently handled.");
                }
        }

        // determine the direct order
        std::cout << "---- SubMapping::init determining the direct order ----" << std::endl;
        bool isDirect = true;

        if (orderedVertices[0] == orderedCurves[0]->getFirstPoint()) {
                if (ASurf == orderedCurves[0]->getLeftSurface()) {
                        // we do nothing, we are in direct order
                } else {
                        isDirect = false;
                }
        } else {
                throw GMDSException(
                    "SubMapping::init orderedVertices[0] != orderedCurves[0]->getFirstPoint() that should not happen.");
        }

        std::cout << "---- SubMapping::init determining the direct order ---- " << isDirect << std::endl;

        if (!isDirect) {
                std::reverse(++orderedVertices.begin(), orderedVertices.end());
                std::reverse(orderedCurves.begin(), orderedCurves.end());
                for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                        isOrderedDirectCurves[curves[iCurve]] = not(isOrderedDirectCurves[curves[iCurve]]);
                }
        }

        // offset in order to always start from the same leftest / bottomest / frontest vertex
        int firstVertexIndex = 0;
        gmds::geom::GeomPoint* firstVertex = orderedVertices[firstVertexIndex];

        for (int iVertex = 1; iVertex < orderedVertices.size(); iVertex++) {
                gmds::geom::GeomPoint* vertex = orderedVertices[iVertex];
                if (vertex->X() < firstVertex->X()) {
                        firstVertex = vertex;
                        firstVertexIndex = iVertex;
                } else if (vertex->X() == firstVertex->X()) {
                        if (vertex->Y() < firstVertex->Y()) {
                                firstVertex = vertex;
                                firstVertexIndex = iVertex;
                        } else if (vertex->Y() == firstVertex->Y()) {
                                if (vertex->Z() < firstVertex->Z()) {
                                        firstVertex = vertex;
                                        firstVertexIndex = iVertex;
                                } else if (vertex->Z() == firstVertex->Z()) {
                                        throw GMDSException(
                                            "SubMapping::init two vertices of the model are at the same position.");
                                }
                        }
                }
        }
        std::vector<gmds::geom::GeomPoint*> orderedVertices_tmp = orderedVertices;
        std::vector<gmds::geom::GeomCurve*> orderedCurves_tmp = orderedCurves;
        int currentIndex = firstVertexIndex;
        for (int iVertex = 0; iVertex < orderedVertices.size(); iVertex++, currentIndex++) {
                if (currentIndex == orderedVertices.size()) {
                        currentIndex = 0;
                }
                orderedVertices[iVertex] = orderedVertices_tmp[currentIndex];
                orderedCurves[iVertex] = orderedCurves_tmp[currentIndex];
        }

        m_orderedDirectVertices[ASurf] = orderedVertices;
        m_orderedDirectCurves[ASurf] = orderedCurves;
        m_isOrderedDirectCurves[ASurf] = isOrderedDirectCurves;

        for (int iVertex = 0; iVertex < orderedVertices.size(); iVertex++) {
                // std::cout<<orderedVertices[iVertex]->getPoint()<<std::endl;
        }
        for (int iCurve = 0; iCurve < orderedCurves.size(); iCurve++) {
                // std::cout<<orderedCurves[iCurve]<<std::endl;
        }

        // compute the angle at each vertex seen from this surface
        std::cout << "---- SubMapping::init computing angles ----" << std::endl;

        std::map<gmds::geom::GeomPoint*, double> angleAtVertexFromSurface;

        // it works here because orderedCurves and orderedVertices have the same size
        for (int iVertex = 0; iVertex < orderedVertices.size(); iVertex++) {
                // vector of next curve
                gmds::math::Vector vext;
                orderedCurves[iVertex]->computeVector(*orderedVertices[iVertex], vext);
                // vector of previous curve
                gmds::math::Vector vint;
                orderedCurves[(iVertex + (orderedCurves.size() - 1)) % orderedCurves.size()]->computeVector(
                    *orderedVertices[iVertex], vint);
                // normal of the surface
                gmds::math::Vector vnormal;
                ASurf->computeNormal(orderedVertices[iVertex]->getPoint(), vnormal);

                double angle = vext.orientedAngle(vint, vnormal);
                if (angle < 0.) {
                        angle = gmds::submapping::PI + fabs(angle);
                }

                angleAtVertexFromSurface[orderedVertices[iVertex]] = angle;

                // std::cout<<iVertex<<" "<<orderedCurves.size()<<" "<<(iVertex-1)<<"
                // "<<(iVertex+(orderedCurves.size()-1))%orderedCurves.size()<<std::endl;
                // std::cout<<orderedCurves[iVertex]<<"
                // "<<orderedCurves[(iVertex+(orderedCurves.size()-1))%orderedCurves.size()]<<std::endl;
                // std::cout<<vext<<" "<<vint<<" "<<vnormal<<std::endl;
                // std::cout<<"angle "<<angle<<" "<<180.*(angle/gmds::submapping::PI)<<std::endl;
        }

        m_angleAtVertexFromSurface[ASurf] = angleAtVertexFromSurface;
}
/*----------------------------------------------------------------------------*/
void
SubMapping::exec(gmds::geom::GeomVolume* AVol)
{
        int glpErr = glp_init_env();
        switch (glpErr) {
        case 0:
                // initialization successful
                break;
        case 1:
                // environment is already initialized
                break;
        case 2:
                // initialization failed (insufficient memory)
                throw GMDSException("SubMapping::exec call to glp_init_env failed (insufficient memory)");
        case 3:
                // initialization failed (unsupported programming model)
                throw GMDSException("SubMapping::exec call to glp_init_env failed (unsupported programming model)");
        }

        init(AVol);
        std::vector<gmds::geom::GeomSurface*> surfaces;
        AVol->get(surfaces);
        for (int iSurf = 0; iSurf < surfaces.size(); iSurf++) {
                verticesClassification(surfaces[iSurf]);
                curvesClassification(surfaces[iSurf]);
        }
        intervalAssignment(AVol);

        exportVTKSubmappingInfo(AVol, "submappingInfo_beforeboundaryDiscretization");

        boundaryDiscretization(AVol);

        glpErr = glp_free_env();
}
/*----------------------------------------------------------------------------*/
void
SubMapping::exec(gmds::geom::GeomSurface* ASurf)
{
        int glpErr = glp_init_env();

        init(ASurf);
        verticesClassification(ASurf);
        curvesClassification(ASurf);
        intervalAssignment(ASurf);
        boundaryDiscretization(ASurf);

        glpErr = glp_free_env();
}
/*----------------------------------------------------------------------------*/
void
SubMapping::verticesClassification(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::verticesClassification ----" << std::endl;

        std::map<gmds::geom::GeomPoint*, int> verticesClassification;

        int ends = 0;
        int corners = 0;
        int reversals = 0;

        std::vector<gmds::geom::GeomPoint*> vertices = m_orderedDirectVertices[ASurf];
        // ASurf->get(vertices);

        for (int iVertex = 0; iVertex < vertices.size(); iVertex++) {
                gmds::geom::GeomPoint* vertex = vertices[iVertex];

                double angle = getAngle(ASurf, vertex);

                double alphaEnd = fabs(angle - END_VALUE);
                double alphaSide = fabs(angle - SIDE_VALUE);
                double alphaCorner = fabs(angle - CORNER_VALUE);
                double alphaReversal = fabs(angle - REVERSAL_VALUE);

                if (alphaEnd <= alphaSide && alphaEnd <= alphaCorner && alphaEnd <= alphaReversal) {
                        verticesClassification[vertex] = END;
                        ends++;
                } else if (alphaSide <= alphaEnd && alphaSide <= alphaCorner && alphaSide <= alphaReversal) {
                        verticesClassification[vertex] = SIDE;
                } else if (alphaCorner <= alphaEnd && alphaCorner <= alphaSide && alphaCorner <= alphaReversal) {
                        verticesClassification[vertex] = CORNER;
                        corners++;
                } else if (alphaReversal <= alphaEnd && alphaReversal <= alphaSide && alphaReversal <= alphaCorner) {
                        verticesClassification[vertex] = REVERSAL;
                        reversals++;
                }

                std::cout << vertex->getPoint() << " " << verticesClassification[vertex] << " " << angle << std::endl;
        }

        m_verticesClassification[ASurf] = verticesClassification;

        if (ends - corners - 2 * reversals != 4) {
                correctVerticesClassification(ASurf);
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::correctVerticesClassification(gmds::geom::GeomSurface* ASurf)
{
        std::vector<gmds::geom::GeomPoint*> vertices = m_orderedDirectVertices[ASurf];
        int nbVertices = vertices.size();

        glp_prob* lp;
        int *ia, *ja;
        double* ar;
        double *rho, *omega;
        int i, j, count;

        lp = glp_create_prob();

        /* GLP_MIN means "we want to minimize" */
        glp_set_obj_dir(lp, GLP_MIN);

        /* each row represents one constraint */
        glp_add_rows(lp, 1 + 3 * nbVertices);

        /* we want the first row to be equal to 4 */
        glp_set_row_bnds(lp, 1, GLP_FX, 4, 4);
        for (i = 0; i < nbVertices; i++) {
                /* the next rows represent the second constraints line of the linear problem i.e. Di */
                glp_set_row_bnds(lp, i + 2, GLP_LO, -1 * getClassification(ASurf, vertices[i]), 0);
        }
        for (; i < 2 * nbVertices; i++) {
                /* the third constraints line i.e. di */
                glp_set_row_bnds(lp, i + 2, GLP_LO, getClassification(ASurf, vertices[i - nbVertices]), 0);
        }
        for (; i < 3 * nbVertices; i++) {
                /* the fourth constraints line */
                glp_set_row_bnds(lp, i + 2, GLP_UP, 0, 1);
        }

        /* all classifications + Di for each vertices i + di for each vertices i */
        glp_add_cols(lp, 3 * nbVertices);
        for (i = 0; i < nbVertices; i++) {
                glp_set_col_kind(lp, i + 1, GLP_IV);
                glp_set_col_bnds(lp, i + 1, GLP_DB, -2, 1);
                glp_set_obj_coef(lp, i + 1, 0);
        }

        /* the next 14 lines are for coefficients generation */
        rho = (double*) malloc(nbVertices * sizeof(double));
        omega = (double*) malloc(nbVertices * sizeof(double));

        computeRhoAndOmega(ASurf, rho, omega);
        for (; i < 2 * nbVertices; i++) {
                glp_set_col_bnds(lp, i + 1, GLP_LO, 0, 0);
                glp_set_obj_coef(lp, i + 1, rho[i - nbVertices]);
        }
        for (; i < 3 * nbVertices; i++) {
                glp_set_col_bnds(lp, i + 1, GLP_LO, 0, 0);
                glp_set_obj_coef(lp, i + 1, omega[i - 2 * nbVertices]);
        }
        free(rho);
        free(omega);

        /* dynamic allocation of the matrix where ia represents the row index, ja the col index and ar the entry */
        ia = (int*) calloc(1 + ((1 + 3 * nbVertices) * (3 * nbVertices)), sizeof(int));
        ja = (int*) calloc(1 + ((1 + 3 * nbVertices) * (3 * nbVertices)), sizeof(int));
        ar = (double*) calloc(1 + ((1 + 3 * nbVertices) * (3 * nbVertices)), sizeof(double));

        count = 1;
        for (i = 0; i < 1 + 3 * nbVertices; i++) {
                for (j = 0; j < 3 * nbVertices; j++) {
                        ia[count] = i + 1;
                        ja[count] = j + 1;
                        ar[count] = 0;
                        count++;
                }
        }

        int rowLength = 3 * nbVertices;
        int indexRow = 1;

        // sum alphai == 4
        for (int icol = 0; icol < nbVertices; icol++) {
                ar[indexRow + icol] = 1;
        }
        indexRow += rowLength;

        // Di >= alphai - alphai_bar
        for (int irow = 0; irow < nbVertices; irow++) {
                // alphai
                ar[indexRow + irow] = -1;
                // Di
                ar[indexRow + nbVertices + irow] = 1;

                indexRow += rowLength;
        }

        // di >= alphai_bar - alphai
        for (int irow = 0; irow < nbVertices; irow++) {
                // alphai
                ar[indexRow + irow] = 1;
                // Di
                ar[indexRow + 2 * nbVertices + irow] = 1;

                indexRow += rowLength;
        }

        // Di + di <= 1
        for (int irow = 0; irow < nbVertices; irow++) {
                // Di
                ar[indexRow + nbVertices + irow] = 1;
                // di
                ar[indexRow + 2 * nbVertices + irow] = 1;

                indexRow += rowLength;
        }

        glp_load_matrix(lp, ((1 + 3 * nbVertices) * (3 * nbVertices)), ia, ja, ar);

        /* the next line disables the default terminal report */
        glp_term_out(GLP_ON);

        glp_iocp glpParams;
        glp_init_iocp(&glpParams);
        glpParams.presolve = GLP_ON;

        glp_write_lp(lp, NULL, "cplex.txt");

        int glpErr = 0;
        glpErr = glp_intopt(lp, &glpParams);
        switch (glpErr) {
        case 0:
                std::cout << "GLP OK" << std::endl;
                break;
        default:
                std::cout << "pb solving in GLP." << std::endl;
                throw GMDSException("SubMapping::correctVerticesClassification pb solving in GLPK.");
                break;
        }
        // glp_print_sol( lp, "poyop.txt");

        for (i = 0; i < nbVertices; i++) {
                /* we want the new classifications so we just need to stop at numberOfVertices (the classifications) */
                setClassification(ASurf, vertices[i], glp_mip_col_val(lp, i + 1));
        }

        /* free everything */
        free(ia);
        free(ja);
        free(ar);

        glp_delete_prob(lp);
}
/*----------------------------------------------------------------------------*/
void
SubMapping::computeRhoAndOmega(gmds::geom::GeomSurface* ASurf, double* rho, double* omega)
{
        std::vector<gmds::geom::GeomPoint*> vertices = m_orderedDirectVertices[ASurf];
        const int nbVertices = vertices.size();

        int i;
        int t1_x, t1_y, t2_x, t2_y;
        double n1_x, n1_y, n2_x, n2_y;
        double d1, d2, d;
        for (i = 0; i < nbVertices; i++) {
                if (i == 0) {
                        t1_x = vertices[nbVertices - 1]->X() - vertices[0]->X();
                        t1_y = vertices[nbVertices - 1]->Y() - vertices[0]->Y();
                        n1_x = vertices[nbVertices - 1]->Y() - vertices[0]->Y();
                        n1_y = vertices[0]->X() - vertices[nbVertices - 1]->X();
                        t2_x = vertices[1]->X() - vertices[0]->X();
                        t2_y = vertices[1]->Y() - vertices[0]->Y();
                        n2_x = vertices[1]->Y() - vertices[0]->Y();
                        n2_y = vertices[0]->X() - vertices[1]->X();
                } else if (i + 1 == nbVertices) {
                        t1_x = vertices[i - 1]->X() - vertices[i]->X();
                        t1_y = vertices[i - 1]->Y() - vertices[i]->Y();
                        n1_x = vertices[i - 1]->Y() - vertices[i]->Y();
                        n1_y = vertices[i]->X() - vertices[i - 1]->X();
                        t2_x = vertices[0]->X() - vertices[i]->X();
                        t2_y = vertices[0]->Y() - vertices[i]->Y();
                        n2_x = vertices[0]->Y() - vertices[i]->Y();
                        n2_y = vertices[i]->X() - vertices[0]->X();
                } else {
                        t1_x = vertices[i - 1]->X() - vertices[i]->X();
                        t1_y = vertices[i - 1]->Y() - vertices[i]->Y();
                        n1_x = vertices[i - 1]->Y() - vertices[i]->Y();
                        n1_y = vertices[i]->X() - vertices[i - 1]->X();
                        t2_x = vertices[i + 1]->X() - vertices[i]->X();
                        t2_y = vertices[i + 1]->Y() - vertices[i]->Y();
                        n2_x = vertices[i + 1]->Y() - vertices[i]->Y();
                        n2_y = vertices[i]->X() - vertices[i + 1]->X();
                }
                d1 = (t1_x * n1_y) - (n1_x * t1_y);
                d2 = (t2_x * n2_y) - (n2_x * t2_y);

                d = (d1 + d2) / 2;
                if (d > 0) {
                        rho[i] = 5 / 8 * (1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i])));
                        omega[i] =
                            3 / 8 * (1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i])));
                } else if (d < 0) {
                        rho[i] = 3 / 8 * (1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i])));
                        omega[i] =
                            5 / 8 * (1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i])));
                } else {
                        rho[i] = 1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i]));
                        omega[i] = 1 - 2 * (ceil(getAngle(ASurf, vertices[i])) - getAngle(ASurf, vertices[i]));
                }
                // temporary hack
                rho[i] = 1 - 2 * (fabs(getAngle(ASurf, vertices[i]) - nearbyint(getAngle(ASurf, vertices[i]))));
                omega[i] = 1 - 2 * (fabs(getAngle(ASurf, vertices[i]) - nearbyint(getAngle(ASurf, vertices[i]))));
                // rho[i] = 1;
                // omega[i] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::curvesClassification(gmds::geom::GeomSurface* ASurf)
{
        std::map<gmds::geom::GeomCurve*, gmds::submapping::CurveClassification> curvesClassification;

        std::vector<gmds::geom::GeomCurve*> orderedDirectCurves = m_orderedDirectCurves[ASurf];
        std::vector<gmds::geom::GeomPoint*> orderedDirectVertices = m_orderedDirectVertices[ASurf];

        curvesClassification[orderedDirectCurves[0]] = PLUS_I;

        for (int iCurve = 1; iCurve < orderedDirectCurves.size(); iCurve++) {
                switch (curvesClassification[orderedDirectCurves[iCurve - 1]]) {
                case PLUS_I:
                        switch (getClassification(ASurf, orderedDirectVertices[iCurve])) {
                        case END:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_J;
                                break;
                        case SIDE:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_I;
                                break;
                        case CORNER:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_J;
                                break;
                        case REVERSAL:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_I;
                                break;
                        default:
                                throw GMDSException("SubMapping::curvesClassification bad vertex classification");
                        }
                        break;
                case PLUS_J:
                        switch (getClassification(ASurf, orderedDirectVertices[iCurve])) {
                        case END:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_I;
                                break;
                        case SIDE:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_J;
                                break;
                        case CORNER:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_I;
                                break;
                        case REVERSAL:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_J;
                                break;
                        default:
                                throw GMDSException("SubMapping::curvesClassification bad vertex classification");
                        }
                        break;
                case MINUS_I:
                        switch (getClassification(ASurf, orderedDirectVertices[iCurve])) {
                        case END:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_J;
                                break;
                        case SIDE:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_I;
                                break;
                        case CORNER:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_J;
                                break;
                        case REVERSAL:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_I;
                                break;
                        default:
                                throw GMDSException("SubMapping::curvesClassification bad vertex classification");
                        }
                        break;
                case MINUS_J:
                        switch (getClassification(ASurf, orderedDirectVertices[iCurve])) {
                        case END:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_I;
                                break;
                        case SIDE:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_J;
                                break;
                        case CORNER:
                                curvesClassification[orderedDirectCurves[iCurve]] = MINUS_I;
                                break;
                        case REVERSAL:
                                curvesClassification[orderedDirectCurves[iCurve]] = PLUS_J;
                                break;
                        default:
                                throw GMDSException("SubMapping::curvesClassification bad vertex classification");
                        }
                        break;
                default:
                        throw GMDSException("SubMapping::curvesClassification bad curve classification");
                }
        }

        m_curvesClassification[ASurf] = curvesClassification;
}
/*----------------------------------------------------------------------------*/
void
SubMapping::intervalAssignment(gmds::geom::GeomVolume* AVol, double ASize)
{
        std::cout << "---- SubMapping::intervalAssignment ----" << std::endl;

        std::vector<gmds::geom::GeomCurve*> curves;
        AVol->get(curves);

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                double length = curves[iCurve]->length();
                int nbEdges = ceil(length / ASize);
                m_targetDiscretization[curves[iCurve]] = nbEdges;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::intervalAssignment(gmds::geom::GeomSurface* ASurf, double ASize)
{
        std::cout << "---- SubMapping::intervalAssignment ----" << std::endl;

        std::vector<gmds::geom::GeomCurve*> curves;
        ASurf->get(curves);

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                double length = curves[iCurve]->length();
                int nbEdges = ceil(length / ASize);
                m_targetDiscretization[curves[iCurve]] = nbEdges;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::boundaryDiscretization(gmds::geom::GeomVolume* AVol)
{
        std::cout << "---- SubMapping::boundaryDiscretization ----" << std::endl;

        std::vector<gmds::geom::GeomSurface*> surfaces;
        std::vector<gmds::geom::GeomCurve*> curves;
        AVol->get(surfaces);
        AVol->get(curves);

        boundaryDiscretization(surfaces, curves);
}
/*----------------------------------------------------------------------------*/
void
SubMapping::boundaryDiscretization(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::boundaryDiscretization ----" << std::endl;

        std::vector<gmds::geom::GeomSurface*> surfaces;
        std::vector<gmds::geom::GeomCurve*> curves;
        surfaces.push_back(ASurf);
        ASurf->get(curves);

        boundaryDiscretization(surfaces, curves);
}
/*----------------------------------------------------------------------------*/
void
SubMapping::boundaryDiscretization(std::vector<gmds::geom::GeomSurface*> ASurfaces,
                                   std::vector<gmds::geom::GeomCurve*> ACurves)
{
        std::cout << "---- SubMapping::boundaryDiscretization ----" << std::endl;

        glp_prob* lp;
        int *ia, *ja;
        double* ar;
        int count;

        int nbSurfaces = ASurfaces.size();
        int nbCurves = ACurves.size();

        std::map<gmds::geom::GeomCurve*, int> curves2ID;
        for (int i = 0; i < nbCurves; i++) {
                curves2ID[ACurves[i]] = i;
        }

        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);

        //
        int nbRows = 2 * nbSurfaces + 3 * nbCurves;
        // columns for ne, De, de, M and m
        int nbCols = 3 * nbCurves + 2;

        glp_add_rows(lp, nbRows);
        glp_add_cols(lp, nbCols);
        int irow = 1;
        int icol = 1;

        // this is for the i+i- / j+j- equalities
        for (int i = 0; i < 2 * nbSurfaces; i++, irow++) {
                glp_set_row_bnds(lp, irow, GLP_FX, 0, 0);
        }
        // this is for De - de - ne = - Ne
        for (int i = 0; i < nbCurves; i++, irow++) {
                glp_set_row_bnds(lp, irow, GLP_FX, -m_targetDiscretization[ACurves[i]], 0);
        }
        // this is for M - ne >= - Ne
        for (int i = 0; i < nbCurves; i++, irow++) {
                glp_set_row_bnds(lp, irow, GLP_LO, -m_targetDiscretization[ACurves[i]], 0);
        }
        // this is for m - ne =< - Ne
        for (int i = 0; i < nbCurves; i++, irow++) {
                glp_set_row_bnds(lp, irow, GLP_UP, 0, -m_targetDiscretization[ACurves[i]]);
        }

        // this is to enforce a minimum target size
        for (int i = 0; i < nbCurves; i++, icol++) {
                glp_set_col_kind(lp, icol, GLP_IV);
                glp_set_col_bnds(lp, icol, GLP_LO, m_targetDiscretization[ACurves[i]], 0);
                glp_set_obj_coef(lp, icol, 0);
        }
        // De >= 0
        for (int i = 0; i < nbCurves; i++, icol++) {
                glp_set_col_kind(lp, icol, GLP_IV);
                glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);
                double length = ACurves[i]->length();
                if (length == 0.) {
                        throw GMDSException("SubMapping::boundaryDiscretization a curve has length zero.");
                }
                double we = 1. / length;
                glp_set_obj_coef(lp, icol, we);
        }
        // de >= 0
        for (int i = 0; i < nbCurves; i++, icol++) {
                glp_set_col_kind(lp, icol, GLP_IV);
                glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);
                double length = ACurves[i]->length();
                if (length == 0.) {
                        throw GMDSException("SubMapping::boundaryDiscretization a curve has length zero.");
                }
                double we = 1. / length;
                glp_set_obj_coef(lp, icol, we);
        }
        // M
        glp_set_col_kind(lp, icol, GLP_IV);
        // glp_set_col_bnds( lp, icol, GLP_LO, 0, 0 );
        glp_set_col_bnds(lp, icol, GLP_FR, 0, 0);
        glp_set_obj_coef(lp, icol, 1);
        icol++;
        // m
        glp_set_col_kind(lp, icol, GLP_IV);
        // glp_set_col_bnds( lp, icol, GLP_LO, 0, 0 );
        glp_set_col_bnds(lp, icol, GLP_FR, 0, 0);
        glp_set_obj_coef(lp, icol, -1);
        icol++;

        ia = (int*) calloc(1 + nbRows * nbCols, sizeof(int));
        ja = (int*) calloc(1 + nbRows * nbCols, sizeof(int));
        ar = (double*) calloc(1 + nbRows * nbCols, sizeof(double));

        count = 1;
        for (int i = 0; i < nbRows; i++) {
                for (int j = 0; j < nbCols; j++) {
                        ia[count] = i + 1;
                        ja[count] = j + 1;
                        count++;
                }
        }

        // fill ar with zeros
        for (int i = 1; i < 1 + nbRows * nbCols; i++) {
                ar[i] = 0.;
        }

        // fill ar for the i+i- / j+j- equalities
        int indexRow = 1;
        for (int iSurf = 0; iSurf < nbSurfaces; iSurf++) {
                gmds::geom::GeomSurface* surf = ASurfaces[iSurf];
                std::vector<gmds::geom::GeomCurve*> curves = m_orderedDirectCurves[surf];

                for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                        gmds::geom::GeomCurve* curve = curves[iCurve];

                        switch (m_curvesClassification[surf][curve]) {
                        case PLUS_I:
                                ar[indexRow + curves2ID[curve]] = 1;
                                break;
                        case MINUS_I:
                                ar[indexRow + curves2ID[curve]] = -1;
                                break;
                        default:
                                break;
                        }
                }
                indexRow += nbCols;

                for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                        gmds::geom::GeomCurve* curve = curves[iCurve];

                        switch (m_curvesClassification[surf][curve]) {
                        case PLUS_J:
                                ar[indexRow + curves2ID[curve]] = 1;
                                break;
                        case MINUS_J:
                                ar[indexRow + curves2ID[curve]] = -1;
                                break;
                        default:
                                break;
                        }
                }
                indexRow += nbCols;
        }
        // De - de - ne = - Ne
        for (int iCurve = 0; iCurve < ACurves.size(); iCurve++) {
                // ne
                ar[indexRow + iCurve] = -1;
                // De
                ar[indexRow + nbCurves + iCurve] = 1;
                // de
                ar[indexRow + 2 * nbCurves + iCurve] = -1;

                indexRow += nbCols;
        }
        // M - ne >= - Ne
        for (int iCurve = 0; iCurve < ACurves.size(); iCurve++) {
                // ne
                ar[indexRow + iCurve] = -1;
                // M
                ar[indexRow + 3 * nbCurves] = 1;

                indexRow += nbCols;
        }
        // m - ne =< - Ne
        for (int iCurve = 0; iCurve < ACurves.size(); iCurve++) {
                // ne
                ar[indexRow + iCurve] = -1;
                // m
                ar[indexRow + 3 * nbCurves + 1] = 1;

                indexRow += nbCols;
        }

        glp_load_matrix(lp, nbRows * nbCols, ia, ja, ar);

        // glp_term_out( GLP_OFF );
        glp_term_out(GLP_ON);

        glp_iocp glpParams;
        glp_init_iocp(&glpParams);
        glpParams.presolve = GLP_ON;

        glp_write_lp(lp, NULL, "cplex.txt");

        int glpErr = 0;
        // glpErr = glp_simplex( lp, NULL );
        glpErr = glp_intopt(lp, &glpParams);
        switch (glpErr) {
        case 0:
                std::cout << "GLP OK" << std::endl;
                break;
        default:
                std::cout << "pb solving in GLP." << std::endl;
                break;
        }
        // glp_print_sol( lp, "poyop.txt");

        glpErr = glp_mip_status(lp);
        switch (glpErr) {
        case GLP_UNDEF:
                std::cout << " MIP solution is undefined" << std::endl;
                throw GMDSException("SubMapping::boundaryDiscretization MIP solution is undefined");
                break;
        case GLP_OPT:
                std::cout << " MIP solution is integer optimal" << std::endl;
                break;
        case GLP_FEAS:
                std::cout << " MIP solution is integer feasible" << std::endl;
                break;
        case GLP_NOFEAS:
                std::cout << "problem has no integer feasible solution" << std::endl;
                throw GMDSException("SubMapping::boundaryDiscretization problem has no integer feasible solution");
                break;
        default:
                throw GMDSException("SubMapping::boundaryDiscretization glp_intopt unknown return code.");
        }

        for (int iCurve = 0; iCurve < nbCurves; iCurve++) {
                m_computedDiscretization[ACurves[iCurve]] = glp_mip_col_val(lp, iCurve + 1);
        }

        free(ia);
        free(ja);
        free(ar);

        glp_delete_prob(lp);
}
/*----------------------------------------------------------------------------*/
int
SubMapping::getClassification(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint)
{
        return m_verticesClassification[ASurf][APoint];
}
/*----------------------------------------------------------------------------*/
void
SubMapping::setClassification(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint, const int AClassification)
{
        m_verticesClassification[ASurf][APoint] = AClassification;
}
/*----------------------------------------------------------------------------*/
double
SubMapping::getAngle(gmds::geom::GeomSurface* ASurf, gmds::geom::GeomPoint* APoint)
{
        return m_angleAtVertexFromSurface[ASurf][APoint];
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printVerticesClassification(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface*> surfaces;
        AVol->get(surfaces);

        for (int iSurf = 0; iSurf < surfaces.size(); iSurf++) {
                printVerticesClassification(surfaces[iSurf]);
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printCurvesClassification(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface*> surfaces;
        AVol->get(surfaces);

        for (int iSurf = 0; iSurf < surfaces.size(); iSurf++) {
                printCurvesClassification(surfaces[iSurf]);
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printVerticesClassification(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::printVerticesClassification ----" << std::endl;

        std::vector<gmds::geom::GeomPoint*> vertices = m_orderedDirectVertices[ASurf];

        for (int iVertex = 0; iVertex < vertices.size(); iVertex++) {
                std::cout << verticesClassificationStr[m_verticesClassification[ASurf][vertices[iVertex]]] << " "
                          << m_angleAtVertexFromSurface[ASurf][vertices[iVertex]] << " "
                          << m_verticesClassification[ASurf][vertices[iVertex]] << std::endl;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printCurvesClassification(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::printCurvesClassification ----" << std::endl;

        std::vector<gmds::geom::GeomCurve*> curves = m_orderedDirectCurves[ASurf];

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                std::cout << curvesClassificationStr[m_curvesClassification[ASurf][curves[iCurve]]] << " "
                          << m_curvesClassification[ASurf][curves[iCurve]] << std::endl;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printCurvesDiscretization(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::printCurvesDiscretization ----" << std::endl;

        std::vector<gmds::geom::GeomCurve*> curves = m_orderedDirectCurves[ASurf];

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                std::cout << m_targetDiscretization[curves[iCurve]] << " " << m_computedDiscretization[curves[iCurve]]
                          << std::endl;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::printSubmappingInfo(gmds::geom::GeomSurface* ASurf)
{
        std::cout << "---- SubMapping::printSubmappingInfo ----" << std::endl;

        std::vector<gmds::geom::GeomCurve*> curves = m_orderedDirectCurves[ASurf];

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                std::cout << curvesClassificationStr[m_curvesClassification[ASurf][curves[iCurve]]] << " "
                          << m_targetDiscretization[curves[iCurve]] << " " << m_computedDiscretization[curves[iCurve]]
                          << std::endl;
        }
}
/*----------------------------------------------------------------------------*/
void
SubMapping::exportVTKSubmappingInfo(gmds::geom::GeomVolume* AVol, const std::string& AFile)
{
        std::cout << "---- SubMapping::exportVTKSubmappingInfo ----" << std::endl;

        gmds::MeshModel model = DIM3 | N | F | R | F2N | R2N;
        gmds::IGMesh mesh(model);

        // we will represent the curves using degenerated faces.
        gmds::Variable<int>* curveDiscretizationVariable =
            mesh.newVariable<int>(GMDS_FACE, "curveDiscretizationVariable");

        std::vector<gmds::geom::GeomCurve*> curves;
        AVol->get(curves);

        for (int iCurve = 0; iCurve < curves.size(); iCurve++) {
                gmds::geom::FacetedCurve* current_curve = dynamic_cast<gmds::geom::FacetedCurve*>(curves[iCurve]);

                std::vector<gmds::Edge> edges;
                current_curve->getMeshEdges(edges);

                for (int iEdge = 0; iEdge < edges.size(); iEdge++) {
                        std::vector<gmds::Node> nodes = edges[iEdge].get<gmds::Node>();

                        gmds::Node n0 = mesh.newNode(nodes[0].getPoint());
                        gmds::Node n1 = mesh.newNode(nodes[1].getPoint());

                        gmds::Face f = mesh.newTriangle(n0, n1, n0);
                        if (m_computedDiscretization.find(current_curve) != m_computedDiscretization.end()) {
                                (*curveDiscretizationVariable)[f.getID()] = m_computedDiscretization[current_curve];
                        } else {
                                (*curveDiscretizationVariable)[f.getID()] = -1;
                        }
                }
        }

        // we will represent the vertices using small triangles, as seen by each adjacent surface
        gmds::Variable<int>* verticesAngleVariable = mesh.newVariable<int>(GMDS_FACE, "verticesAngleVariable");
        gmds::Variable<int>* verticesClassificationVariable =
            mesh.newVariable<int>(GMDS_FACE, "verticesClassificationVariable");

        std::vector<gmds::geom::GeomSurface*> surfaces;
        AVol->get(surfaces);

        for (int iSurface = 0; iSurface < surfaces.size(); iSurface++) {
                gmds::geom::GeomSurface* surf = surfaces[iSurface];

                std::vector<gmds::geom::GeomCurve*> orderedDirectCurves = m_orderedDirectCurves[surf];
                std::vector<gmds::geom::GeomPoint*> orderedDirectVertices = m_orderedDirectVertices[surf];

                for (int iVertex = 0; iVertex < orderedDirectVertices.size(); iVertex++) {
                        gmds::geom::GeomPoint* point = orderedDirectVertices[iVertex];
                        gmds::geom::GeomCurve* curve1 = orderedDirectCurves[iVertex];
                        gmds::geom::GeomCurve* curve2 =
                            orderedDirectCurves[(iVertex + orderedDirectVertices.size() - 1) %
                                                orderedDirectVertices.size()];

                        // get the adjacent curves vectors
                        gmds::math::Vector t1;
                        curve1->computeVector(*point, t1);
                        gmds::math::Vector t2;
                        curve2->computeVector(*point, t2);

                        gmds::math::Vector t_middle = t1 + t2;

                        double angle = getAngle(surf, point);

                        if (angle > SIDE_VALUE) {
                                t_middle = (-1) * t_middle;
                        } else {
                                if (angle == SIDE_VALUE) {
                                        // t_middle = gmds::math::Vector(t1.Y(),);
                                }
                        }

                        t_middle.normalize();

                        double length1 = curve1->length();
                        double length2 = curve2->length();
                        double length_middle = (length1 + length2) / 2.;

                        t1 = (1. / 10.) * length1 * t1;
                        t2 = (1. / 10.) * length2 * t2;
                        t_middle = (1. / 10.) * length_middle * t_middle;

                        gmds::Node n0 = mesh.newNode(point->getPoint());
                        gmds::Node n1 = mesh.newNode(n0.getPoint() + t1);
                        gmds::Node n2 = mesh.newNode(n0.getPoint() + t2);
                        gmds::Node n_middle = mesh.newNode(n0.getPoint() + t_middle);

                        gmds::Face f1 = mesh.newTriangle(n0, n1, n_middle);
                        gmds::Face f2 = mesh.newTriangle(n0, n_middle, n2);

                        (*verticesAngleVariable)[f1.getID()] = angle * (180. / gmds::submapping::PI);
                        (*verticesAngleVariable)[f2.getID()] = angle * (180. / gmds::submapping::PI);
                        (*verticesClassificationVariable)[f1.getID()] = m_verticesClassification[surf][point];
                        (*verticesClassificationVariable)[f2.getID()] = m_verticesClassification[surf][point];
                }
        }

        gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write(AFile, N | F);
}
/*----------------------------------------------------------------------------*/
}  // end namespace submapping
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
