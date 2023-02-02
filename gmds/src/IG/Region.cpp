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
 * Region.cpp
 *
 *  Created on: 20 may 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Region.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/Cell.h>
#include <GMDS/IG/RegionContainer.h>
#include <GMDS/Math/Hexahedron.h>
#include <GMDS/Math/Tetrahedron.h>
#include <GMDS/Math/Pyramid.h>
#include <GMDS/Math/Prism3.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    Region::
    Region()
    : Cell(0,GMDS_TETRA,NullID),m_regions_container(0),m_type_id(NullID)
    {
    }
    /*----------------------------------------------------------------------------*/
    Region::
    Region(IGMesh* AMesh, const ECellType AType, const TCellID& AID)
    : Cell(AMesh,AType,AID) //mesh, type and id are filled in
    {
        //============================================
        // we keep a reference on the face container
        if(AMesh!=0){
            m_regions_container = AMesh->m_regions_container;
            m_type_id = m_regions_container->getTypeID(AID);
        }
        else {
            m_regions_container = 0;
            m_type_id = NullID;
        }
        
    }
    /*----------------------------------------------------------------------------*/
    Region::
    Region(const Region& AReg)
    : Cell(AReg.m_owner,AReg.m_type,AReg.m_id)
    {
        if(m_owner!=0)
            m_regions_container = m_owner->m_regions_container;
        else
            m_regions_container = 0;
        
        m_type_id = AReg.m_type_id;
    }
    /*----------------------------------------------------------------------------*/
    Region::~Region()
    {}
    /*----------------------------------------------------------------------------*/
    bool Region::operator==(const Region& ARegion) const
    {
        return (m_owner == ARegion.m_owner && m_id==ARegion.m_id);
    }
    /*----------------------------------------------------------------------------*/
    bool Region::operator!=(const Region& ARegion) const
    {
        return (!(*this == ARegion));
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::getNbNodes() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val= 4;
        else if (m_type==GMDS_HEX)
            val= 8;
        else if (m_type==GMDS_PYRAMID)
            val= 5;
        else if (m_type==GMDS_PRISM3)
            val= 6;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2N)[m_type_id];
            val= cells.size();
        }
        return val;
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::getNbEdges() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=6;
        else if (m_type==GMDS_HEX)
            val=12;
        else if (m_type==GMDS_PYRAMID)
            val=8;
        else if (m_type==GMDS_PRISM3)
            val=9;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2E)[m_type_id];
            val=cells.size();
        }
        return val;
        
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::getNbFaces() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=4;
        else if (m_type==GMDS_HEX)
            val=6;
        else if (m_type==GMDS_PYRAMID || m_type==GMDS_PRISM3)
            val=5;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2F)[m_type_id];
            val=cells.size();
        }
        return val;
        
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::getNbRegions() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=4;
        else if (m_type==GMDS_PYRAMID || m_type==GMDS_PRISM3)
            val=5;
        else if (m_type==GMDS_HEX)
            val=6;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2R)[m_type_id];
            val=cells.size();
        }
        return val;
    }
    /*----------------------------------------------------------------------------*/
    math::Point
    Region::center() const
    {
        TCoord p_coords[3] = {0.0,0.0,0.0};
        
        std::vector<Node> nodes = this->get<Node>();
        int nb_nodes = nodes.size();
        
        for(int i=0; i<nb_nodes; i++)
        {
            Node n = nodes[i];
            p_coords[0] += n.X();
            p_coords[1] += n.Y();
            p_coords[2] += n.Z();
        }
        
        p_coords[0] = p_coords[0] / nodes.size();
        p_coords[1] = p_coords[1] / nodes.size();
        p_coords[2] = p_coords[2] / nodes.size();
        
        math::Point p(p_coords[0],p_coords[1],p_coords[2]);
        
        return p;
    }
    
    /*----------------------------------------------------------------------------*/
    TCoord Region::volume() const
    {
        TCoord vol=0.0;
        std::vector<Node> n = this->get<Node>();
        if(this->getType()==GMDS_TETRA){
            vol = math::Tetrahedron(n[0].getPoint(),
                                    n[1].getPoint(),
                                    n[2].getPoint(),
                                    n[3].getPoint()).getVolume();
        }
        else if(this->getType()==GMDS_HEX){
            vol = math::Hexahedron(n[0].getPoint(),
                                   n[1].getPoint(),
                                   n[2].getPoint(),
				   n[3].getPoint(),
				   n[4].getPoint(),
				   n[5].getPoint(),
				   n[6].getPoint(),
                                   n[7].getPoint()).getVolume();
	}
        else if(this->getType()==GMDS_PYRAMID){
            vol = math::Pyramid(n[0].getPoint(),
                                n[1].getPoint(),
                                n[2].getPoint(),
				n[3].getPoint(),
                                n[4].getPoint()).getVolume();
        }
        else if(this->getType()==GMDS_PRISM3){
            vol = math::Prism3(n[0].getPoint(),
                               n[1].getPoint(),
                               n[2].getPoint(),
			       n[3].getPoint(),
			       n[4].getPoint(),
			       n[5].getPoint()).getVolume();
        }
	else {
	  throw GMDSException("Region::volume can not be computed for this cell type.");
	}
        
        return vol;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<math::Point>
    Region::computeNGLLPoints(int ADegree) const
    {
        //	throw GMDSException("Region::computeNGLLPoints not implemented yet.");
        
        std::vector<math::Point> points;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
            {
                points.resize(27);
                
                // we get the 8 corners of the hex
                gmds::math::Point p[8];
                for(int iNode=0; iNode<nodes.size(); iNode++)
                {
                    p[iNode] = nodes[iNode].getPoint();
                }
                // Now we build an inner grid of  27 points
                
                // we build intermediate nodes to create our 27 inner points
                // middle of edge
                gmds::math::Point p01, p12, p23, p03,
                p45, p56, p67, p47, p04, p15, p26, p37;
                
                p01 = (p[0]+p[1])*0.5;
                p12 = (p[1]+p[2])*0.5;
                p23 = (p[2]+p[3])*0.5;
                p03 = (p[0]+p[3])*0.5;
                p45 = (p[4]+p[5])*0.5;
                p56 = (p[5]+p[6])*0.5;
                p67 = (p[6]+p[7])*0.5;
                p47 = (p[4]+p[7])*0.5;
                p04 = (p[0]+p[4])*0.5;
                p15 = (p[1]+p[5])*0.5;
                p26 = (p[2]+p[6])*0.5;
                p37 = (p[3]+p[7])*0.5;																	  //face center
                gmds::math::Point p0123, p4567, p0154, p1265, p2376, p0374;
                p0123 = (p[0]+p[1]+p[2]+p[3])*0.25;
                p4567 = (p[4]+p[5]+p[6]+p[7])*0.25;
                p0154 = (p[0]+p[1]+p[5]+p[4])*0.25;
                p1265 = (p[1]+p[2]+p[6]+p[5])*0.25;
                p2376 = (p[2]+p[3]+p[7]+p[6])*0.25;
                p0374 = (p[0]+p[3]+p[7]+p[4])*0.25;
                //hex center
                gmds::math::Point bary = (p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7])*0.125;
                points[0] = bary;
                points[1] = (bary+p[0])*0.5;
                points[2] = (bary+p[1])*0.5;
                points[3] = (bary+p[2])*0.5;
                points[4] = (bary+p[3])*0.5;
                points[5] = (bary+p[4])*0.5;
                points[6] = (bary+p[5])*0.5;
                points[7] = (bary+p[6])*0.5;
                points[8] = (bary+p[7])*0.5;
                points[9] = (bary+p01)*0.5;
                points[10]= (bary+p12)*0.5;
                points[11]= (bary+p23)*0.5;
                points[12]= (bary+p03)*0.5;
                points[13]= (bary+p45)*0.5;
                points[14]= (bary+p56)*0.5;
                points[15]= (bary+p67)*0.5;
                points[16]= (bary+p47)*0.5;
                points[17]= (bary+p04)*0.5;
                points[18]= (bary+p15)*0.5;
                points[19]= (bary+p26)*0.5;
                points[20]= (bary+p37)*0.5;
                points[21]= (bary+p0123)*0.5;
                points[22]= (bary+p4567)*0.5;
                points[23]= (bary+p0154)*0.5;
                points[24]= (bary+p1265)*0.5;
                points[25]= (bary+p2376)*0.5;
                points[26]= (bary+p0374)*0.5;
            }
                break;
            case GMDS_TETRA:
                for(int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].getPoint());
                }
                break;
            case GMDS_PYRAMID:
                for(int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].getPoint());
                }
                break;
            case GMDS_PRISM3:
                for(int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].getPoint());
                }
                break;
            default:
                throw GMDSException("Region::computeNGLLPoints not implemented yet for this cell type.");
                break;
        }
        
        return points;
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeScaledJacobian() const
    {
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),
                                     nodes[4].getPoint(),nodes[5].getPoint(),nodes[6].getPoint(),nodes[7].getPoint());
                return hex.computeScaledJacobian();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());
                return tet.computeScaledJacobian();
            }
                break;
            case GMDS_PYRAMID:
	    {
		math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());
		return pyr.computeScaledJacobian();
	    }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
                                   nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());
                return prism.computeScaledJacobian();
            }
                break;
            default:
                throw GMDSException("Region::computeScaledJacobian not implemented yet for this cell type.");
                break;
        }
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeNormalizedScaledJacobian() const
    {
        std::vector<Node> nodes = this->get<Node>();

        switch(this->getType()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),
                                     nodes[4].getPoint(),nodes[5].getPoint(),nodes[6].getPoint(),nodes[7].getPoint());
                return hex.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());
                return tet.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_PYRAMID:
            {
                math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());
                return pyr.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
                                   nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());
                return prism.computeNormalizedScaledJacobian();
            }
                break;
            default:
                throw GMDSException("Region::computeNormalizedScaledJacobian not implemented yet for this cell type.");
                break;
        }
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeMeanRatio() const
    {
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),
                                     nodes[4].getPoint(),nodes[5].getPoint(),nodes[6].getPoint(),nodes[7].getPoint());
                return hex.computeMeanRatio();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());
                return tet.computeMeanRatio();
            }
                break;
            case GMDS_PYRAMID:
	    {
		math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());
		return pyr.computeMeanRatio();
	    }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
                                   nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());
                return prism.computeMeanRatio();
            }
                break;
            default:
                throw GMDSException("Region::computeMeanRatio not implemented yet for this cell type.");
                break;
        }
    }
    /*----------------------------------------------------------------------------*/
    std::vector<std::vector<Node> >
    Region::getOrderedNodesFaces() const
    {
        std::vector<std::vector<Node> > orderedNodesFaces;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
                orderedNodesFaces.resize(6);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[7];
                orderedNodesFaces[1][1] = nodes[4];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[6];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[5];
                orderedNodesFaces[2][1] = nodes[4];
                orderedNodesFaces[2][2] = nodes[0];
                orderedNodesFaces[2][3] = nodes[1];
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[7];
                orderedNodesFaces[3][1] = nodes[6];
                orderedNodesFaces[3][2] = nodes[2];
                orderedNodesFaces[3][3] = nodes[3];
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[4];
                orderedNodesFaces[4][1] = nodes[7];
                orderedNodesFaces[4][2] = nodes[3];
                orderedNodesFaces[4][3] = nodes[0];
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[6];
                orderedNodesFaces[5][1] = nodes[5];
                orderedNodesFaces[5][2] = nodes[1];
                orderedNodesFaces[5][3] = nodes[2];
                break;
            case GMDS_TETRA:
                orderedNodesFaces.resize(4);
                orderedNodesFaces[0].resize(3);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[2];
                orderedNodesFaces[0][2] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[3];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[0];
                orderedNodesFaces[3][2] = nodes[3];
                break;
            case GMDS_PYRAMID:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[4];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[4];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[3];
                orderedNodesFaces[3][2] = nodes[4];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[0];
                orderedNodesFaces[4][2] = nodes[4];
                break;
            case GMDS_PRISM3:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[4];
                orderedNodesFaces[0][3] = nodes[3];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[1];
                orderedNodesFaces[1][1] = nodes[2];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[4];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[2];
                orderedNodesFaces[2][1] = nodes[0];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[2][3] = nodes[5];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[0];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[1];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[4];
                orderedNodesFaces[4][2] = nodes[5];
                break;
            default:
                throw GMDSException("Region::getOrderedNodesFaces not implemented for this region type");
                break;
        }
        
        return orderedNodesFaces;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<std::vector<TCellID> >
    Region::getOrderedNodesFacesIDs() const
    {
        std::vector<std::vector<TCellID> > orderedNodesFaces;
        
        std::vector<TCellID> nodes = this->getIDs<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
                orderedNodesFaces.resize(6);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[7];
                orderedNodesFaces[1][1] = nodes[4];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[6];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[5];
                orderedNodesFaces[2][1] = nodes[4];
                orderedNodesFaces[2][2] = nodes[0];
                orderedNodesFaces[2][3] = nodes[1];
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[7];
                orderedNodesFaces[3][1] = nodes[6];
                orderedNodesFaces[3][2] = nodes[2];
                orderedNodesFaces[3][3] = nodes[3];
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[4];
                orderedNodesFaces[4][1] = nodes[7];
                orderedNodesFaces[4][2] = nodes[3];
                orderedNodesFaces[4][3] = nodes[0];
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[6];
                orderedNodesFaces[5][1] = nodes[5];
                orderedNodesFaces[5][2] = nodes[1];
                orderedNodesFaces[5][3] = nodes[2];
                break;
            case GMDS_TETRA:
                orderedNodesFaces.resize(4);
                orderedNodesFaces[0].resize(3);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[2];
                orderedNodesFaces[0][2] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[3];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[0];
                orderedNodesFaces[3][2] = nodes[3];
                break;
            case GMDS_PYRAMID:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[4];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[4];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[3];
                orderedNodesFaces[3][2] = nodes[4];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[0];
                orderedNodesFaces[4][2] = nodes[4];
                break;
            case GMDS_PRISM3:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[4];
                orderedNodesFaces[0][3] = nodes[3];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[1];
                orderedNodesFaces[1][1] = nodes[2];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[4];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[2];
                orderedNodesFaces[2][1] = nodes[0];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[2][3] = nodes[5];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[0];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[1];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[4];
                orderedNodesFaces[4][2] = nodes[5];
                break;
            default:
                throw GMDSException("Region::getOrderedNodesFaces not implemented for this region type");
                break;
        }
        
        return orderedNodesFaces;
    }
    /*----------------------------------------------------------------------------*/
    bool
    Region::isFaceOrientedOutward(std::vector<Node> ANodes) const
    {
      std::vector<TCellID> ids(ANodes.size());
      for(size_t i=0; i<ANodes.size(); i++) {
	ids[i] = ANodes[i].getID();
      }

      try{
	return this->isFaceOrientedOutward(ids);
      } catch(GMDSException& e) {
	throw e;
      }
      
    }
    /*----------------------------------------------------------------------------*/
    bool
    Region::isFaceOrientedOutward(std::vector<TCellID> AIDs) const
    {
        std::vector<std::vector<TCellID> > orderedNodesFaces = this->getOrderedNodesFacesIDs();
        
        // first find the face
        for(int iFace=0; iFace<orderedNodesFaces.size(); iFace++) {
            
            if(AIDs.size() == orderedNodesFaces[iFace].size()) {
                
                int nbNodesMatched = 0;
                
                for(int iNode1=0; iNode1<orderedNodesFaces[iFace].size(); iNode1++) {
                    for(int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                        if(AIDs[iNode2] == orderedNodesFaces[iFace][iNode1]) {
                            nbNodesMatched++;
                        }
                    }
                }
                if(nbNodesMatched == AIDs.size()) {
                    // face is found, now check the orientation
                    if(orderedNodesFaces[iFace].size() < 3) {
                        throw GMDSException("Region::isFaceOrientedOutward face with less than 3 nodes.");
                    }
                    TCellID firstNode = orderedNodesFaces[iFace][0];
                    TCellID secondNode = orderedNodesFaces[iFace][1];
                    
                    for(int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                        if(AIDs[iNode2] == firstNode) {
                            if(AIDs[(iNode2+1)%AIDs.size()] == secondNode) {
                                return true;
                            } else  if (AIDs[(iNode2-1)%AIDs.size()] == secondNode) {
                                return false;
                            } else {
                                throw GMDSException("Region::isFaceOrientedOutward jumbled face.");
                            }
                        }
                    }
                }
            }
        }
        
        throw GMDSException("Region::isFaceOrientedOutward face not found.");
    }
   /*----------------------------------------------------------------------------*/
    std::vector<FakeFace>
    Region::getFakeFaces() const
    {
	std::vector<std::vector<Node> > orderedNodesFaces = this->getOrderedNodesFaces();
	std::vector<FakeFace> fakeFaces;

	for(size_t iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

		std::vector<gmds::TCellID> ids(orderedNodesFaces[iFace].size());
		for(size_t iNode=0; iNode<ids.size(); iNode++) {
			ids[iNode] = orderedNodesFaces[iFace][iNode].getID();
		}

		fakeFaces.push_back(FakeFace(ids));
	}

        return fakeFaces;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<FakeEdge>
    Region::getFakeEdges() const
    {
        std::vector<FakeEdge> fakeEdges;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->getType()) {
            case GMDS_HEX:
                fakeEdges.resize(12);
                fakeEdges[ 0] = FakeEdge(nodes[0].getID(),nodes[1].getID());
                fakeEdges[ 1] = FakeEdge(nodes[1].getID(),nodes[2].getID());
                fakeEdges[ 2] = FakeEdge(nodes[2].getID(),nodes[3].getID());
                fakeEdges[ 3] = FakeEdge(nodes[3].getID(),nodes[0].getID());
                fakeEdges[ 4] = FakeEdge(nodes[4].getID(),nodes[5].getID());
                fakeEdges[ 5] = FakeEdge(nodes[5].getID(),nodes[6].getID());
                fakeEdges[ 6] = FakeEdge(nodes[6].getID(),nodes[7].getID());
                fakeEdges[ 7] = FakeEdge(nodes[7].getID(),nodes[4].getID());
                fakeEdges[ 8] = FakeEdge(nodes[0].getID(),nodes[4].getID());
                fakeEdges[ 9] = FakeEdge(nodes[1].getID(),nodes[5].getID());
                fakeEdges[10] = FakeEdge(nodes[2].getID(),nodes[6].getID());
                fakeEdges[11] = FakeEdge(nodes[3].getID(),nodes[7].getID());
                break;
            case GMDS_TETRA:
                fakeEdges.resize(6);
		fakeEdges[ 0] = FakeEdge(nodes[0].getID(),nodes[1].getID());
		fakeEdges[ 1] = FakeEdge(nodes[1].getID(),nodes[2].getID());
		fakeEdges[ 2] = FakeEdge(nodes[2].getID(),nodes[0].getID());
		fakeEdges[ 3] = FakeEdge(nodes[0].getID(),nodes[3].getID());
		fakeEdges[ 4] = FakeEdge(nodes[1].getID(),nodes[3].getID());
		fakeEdges[ 5] = FakeEdge(nodes[2].getID(),nodes[3].getID());
                break;
            case GMDS_PYRAMID:
                fakeEdges.resize(8);
		fakeEdges[ 0] = FakeEdge(nodes[0].getID(),nodes[1].getID());
                fakeEdges[ 1] = FakeEdge(nodes[1].getID(),nodes[2].getID());
                fakeEdges[ 2] = FakeEdge(nodes[2].getID(),nodes[3].getID());
                fakeEdges[ 3] = FakeEdge(nodes[3].getID(),nodes[0].getID());
                fakeEdges[ 4] = FakeEdge(nodes[0].getID(),nodes[4].getID());
                fakeEdges[ 5] = FakeEdge(nodes[1].getID(),nodes[4].getID());
                fakeEdges[ 6] = FakeEdge(nodes[2].getID(),nodes[4].getID());
                fakeEdges[ 7] = FakeEdge(nodes[3].getID(),nodes[4].getID());
                break;
            case GMDS_PRISM3:
                fakeEdges.resize(9);
		fakeEdges[ 0] = FakeEdge(nodes[0].getID(),nodes[1].getID());
                fakeEdges[ 1] = FakeEdge(nodes[1].getID(),nodes[2].getID());
                fakeEdges[ 2] = FakeEdge(nodes[2].getID(),nodes[0].getID());
                fakeEdges[ 3] = FakeEdge(nodes[3].getID(),nodes[4].getID());
                fakeEdges[ 4] = FakeEdge(nodes[4].getID(),nodes[5].getID());
                fakeEdges[ 5] = FakeEdge(nodes[5].getID(),nodes[3].getID());
                fakeEdges[ 6] = FakeEdge(nodes[0].getID(),nodes[3].getID());
                fakeEdges[ 7] = FakeEdge(nodes[1].getID(),nodes[4].getID());
                fakeEdges[ 8] = FakeEdge(nodes[2].getID(),nodes[5].getID());
                break;
            default:
                throw GMDSException("Region::getFakeEdges not implemented for this region type");
                break;
        }
        
        return fakeEdges;
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Node>& ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Edge>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Face>& ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Region>& ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_regions_container->buildRegion(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Node>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Edge>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Face>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Region>& ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_regions_container->buildRegion(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");

        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        TInt nb_cells = getNbNodes();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_H2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        
        TInt nb_cells = getNbEdges();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_H2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_PY2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_PR2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        TInt nb_cells = getNbFaces();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size()){
                for(int c_id=0; c_id<ACells.size(); c_id++){
                    Face f = m_owner->get<Face>(ACells[c_id]);
                    std::cout<<"Face: "<<ACells[c_id]<<std::endl;
                    std::vector<Node> fn =f.get<Node>();
                    for(int ni=0; ni<fn.size(); ni++){
                        std::cout<<"\t Node "<<fn[ni].getID()<<": "<<fn[ni].getPoint()<<"\n";
                    }
                }
                throw GMDSException("Invalid number of adj. entities");
            }
            (*m_regions_container->m_H2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        
        TInt nb_cells = getNbRegions();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_H2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            
            (*m_regions_container->m_P2R)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    std::ostream & operator << (std::ostream & AStream, const Region & AR)
    {
        AStream<<"Region "<<AR.getID();;
        return AStream;
    }
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
