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
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/VectorND.h>
#include <GMDS/Math/Numerics.h>

#include <GMDS/Math/AxisAngleRotation.h>
/*---------------------------------------------------------------------------*/
// STL File Headers
#include <set>
/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include "CavitySurfacePaver.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;

/*---------------------------------------------------------------------------*/
BndPoint::BndPoint(const gmds::Node& AN, const math::Vector3d& AV)
:node(AN),inward(AV),type(UNDEFINED){;}
/*---------------------------------------------------------------------------*/
CavitySurfaceBnd::
CavitySurfaceBnd(const std::vector<gmds::Node>& ABnd,
                 CavitySurfacePatch* APatch)
{
    patch=APatch;
    
    points.reserve(ABnd.size());
    

    for(auto i=0; i<ABnd.size();i++){
        points.push_back(BndPoint(ABnd[i], math::Vector3d(0,0,0)));
    }
    
    // The orientation will be checked now, by considering the two
    // possible loop orientation and by keeping the one which maximize
    // the number of END comparing to REV+CORNER.
    // WARNING :If we have the same number we will have to find out
    // what to do
    
    // WARNING: it works since the patch triangle are oriented in a
    // consistent manner
    for(auto i=0; i<points.size();i++){
        math::Point cur = this->operator[](i).node.getPoint();
        math::Point nex = next(i).node.getPoint();
        math::Point pre = prev(i).node.getPoint();
        math::Vector3d n = patch->normal(cur);
        math::Vector3d d1(cur, nex);
        math::Vector3d d2(pre, cur);
        math::Vector3d in1 = n.cross(d1);
        math::Vector3d in2 = n.cross(d2);
//        math::Vector3d inv(-in_patch.X(), -in_patch.Y(), -in_patch.Z());
        points[i].inward = in1+in2;
        std::cout<<"Vec in "<< this->operator[](i).node.getID()
        <<": "<<points[i].inward
        <<std::endl;
        
    }
    
    initTypes();
    int end1 = 0, cor_rev1=0;
    for(auto p:points){
        if(p.type==BndPoint::END){
            end1++;
        }
        else if(p.type==BndPoint::CORNER ||
                p.type==BndPoint::REVERSAL){
            cor_rev1++;
        }
    }
    
    
//    math::Vector3d vc(0,0,0);
//    for(auto i=0; i<points.size();i++){
//        math::Point cur = points[i].node.getPoint();
//        math::Point pre = prev(i).node.getPoint();
//        math::Point nex = next(i).node.getPoint();
//        math::Vector3d v1(cur,pre);
//        math::Vector3d v2(cur,nex);
//        v1.normalize();
//        v2.normalize();
//        math::Vector3d n = v1.cross(v2);
//    
//        points[i].inward=v2.cross(n);
//    }
//    
//    initTypes();
//    int end1 = 0, cor_rev1=0;
//    for(auto p:points){
//        if(p.type==BndPoint::END){
//            end1++;
//        }
//        else if(p.type==BndPoint::CORNER ||
//                p.type==BndPoint::REVERSAL){
//            cor_rev1++;
//        }
//    }
    std::cout<<"(end1, cor1) =("<<end1<<", "<<cor_rev1<<")"<<std::endl;
    if(end1<cor_rev1){
        //We invers the inward vectors
        for(auto i =0; i<points.size();i++){
            math::Vector3d ref = points[i].inward;
            math::Vector3d inv(-ref.X(), -ref.Y(), -ref.Z());
            points[i].inward = inv;
        }
        initTypes();
        int end2 = 0, cor_rev2=0;
        for(auto p:points){
            if(p.type==BndPoint::END){
                end2++;
            }
            else if(p.type==BndPoint::CORNER ||
                    p.type==BndPoint::REVERSAL){
                cor_rev2++;
            }
        }
        std::cout<<"(end2, cor2) =("<<end2<<", "<<cor_rev2<<")"<<std::endl;

    }

}
/*---------------------------------------------------------------------------*/
int CavitySurfaceBnd::size() const{
    return points.size();
}
/*---------------------------------------------------------------------------*/
void CavitySurfaceBnd::clear() {
    points.clear();
}
/*---------------------------------------------------------------------------*/
double CavitySurfaceBnd::angle(const int AI) {
    BndPoint cur = this->operator[](AI);
    BndPoint pre = this->operator[](AI-1);
    BndPoint nex = this->operator[](AI+1);
    math::Vector3d u(cur.node.getPoint(),pre.node.getPoint());
    math::Vector3d u_inv(pre.node.getPoint(),cur.node.getPoint());
    math::Vector3d v(cur.node.getPoint(),nex.node.getPoint());
    
    math::Vector3d nu = u_inv.cross(pre.inward);
    math::Vector3d nv = v.cross(cur.inward);
    math::Vector3d n = nu+nv;
    n.normalize();
    return angle(u,v,n);
}
/*---------------------------------------------------------------------------*/
double CavitySurfaceBnd::angle(const math::Vector3d& AV1,
                               const math::Vector3d& AV2,
                               const math::Vector3d& ANormal)
{
    
    math::Vector3d v1 = AV1;
    math::Vector3d v2 = AV2;
    v1.normalize();
    v2.normalize();
    double d = v1.dot(v2);
    math::Vector3d v1xv2 = v1.cross(v2);
    double c = v1xv2.norm();
    if(v1xv2.dot(ANormal)>0)
        c = -c;
    
    return std::atan2(c, d)*math::Constants::INVPIDIV180;
}

/*---------------------------------------------------------------------------*/
BndPoint& CavitySurfaceBnd::operator[](const int AIndex)
{
    int index = AIndex;
    if(index<0 ){
        index =points.size()+index;
    }
    return points[index%points.size()];
}

/*---------------------------------------------------------------------------*/
BndPoint& CavitySurfaceBnd::prev(const int AIndex)
{
    return this->operator[](AIndex-1);
}
/*---------------------------------------------------------------------------*/
BndPoint& CavitySurfaceBnd::next(const int AIndex)
{
    return this->operator[](AIndex+1);
}
/*---------------------------------------------------------------------------*/
void CavitySurfaceBnd::initTypes()
{

    updateTypes(0, points.size()-1);
}
/*---------------------------------------------------------------------------*/
void CavitySurfaceBnd::updateTypes(const int& AI, const int& AJ)
{
    if(AI>=points.size() || AJ>=points.size())
        return;
    
    if(AI>=AJ)
        return;

    std::cout<<"---------- UPDATE PAVING CLASSIFICATION --------"<<std::endl;
    for(auto i = AI; i<=AJ; i++){
        Node cur_node = points[i].node;
        Node nex_node = next(i).node;
        
        Node pre_node = prev(i).node;
        
        math::Point cur = cur_node.getPoint();
        math::Point nex = nex_node.getPoint();
        
        math::Point pre=pre_node.getPoint();
        
        math::Vector3d from(cur,pre);
        math::Vector3d to  (cur,nex);
        from.normalize();
        to.normalize();
        math::Vector3d local_normal = to.cross(points[i].inward);
        
        
        local_normal.normalize();
        double a = angle(from,to,local_normal);
        
        std::cout<<"Node "<<points[i].node.getID();
        std::cout<<" with vec "<<local_normal<<". a= "<<a<<std::endl;
        std::cout<<"prev,cur,next ="<<pre_node.getID()<<", "
        <<cur_node.getID()<<", "<<nex_node.getID()<<std::endl;
        std::cout<<"And inward = "<<points[i].inward<<std::endl;
        std::cout<<"And to = "<<to<<std::endl;
        
        if(40<=a && a<= 130){
            std::cout<<": END"<<std::endl;
            points[i].type=BndPoint::END;
        }
        else if(std::abs(a)>130){
            std::cout<<": SIDE"<<std::endl;
            points[i].type=BndPoint::SIDE;
        }
        else if(-120<=a && a<=-40){
            std::cout<<": CORNER"<<std::endl;
            points[i].type=BndPoint::CORNER;
        }
        else if(std::abs(a)<40){
            std::cout<<": REVERSAL"<<std::endl;
            points[i].type=BndPoint::REVERSAL;
        }
    }
}


/*---------------------------------------------------------------------------*/
void CavitySurfaceBnd::
replace(const int AI, const int AJ,
        std::vector<Node>& AN,
        std::vector<math::Vector3d>& ADir)
{
    if(ADir.size()!=AN.size()+1){
        throw GMDSException("Wrong node/direction for paving bnd replacement");
    }
    
    std::vector<BndPoint> new_points;
    new_points.reserve(points.size());//almost the same size
    BndPoint pi = prev(AI);
    BndPoint new_pi(pi.node, ADir[0]);
    new_points.push_back(new_pi);

    for(auto i=0; i<AN.size();i++){
        BndPoint new_p(AN[i], ADir[i+1]);
        new_points.push_back(new_p);
    }
    BndPoint pj =next(AJ);
    int idx= AJ;
    while(pj.node.getID()!=pi.node.getID()){
        new_points.push_back(pj);
        pj = next(++idx);

    }
    //we switch the containers
    points = new_points;
    
    updateTypes(0, AN.size()+1);
}
/*---------------------------------------------------------------------------*/
CavitySurfacePatch::
CavitySurfacePatch(const std::vector<math::Triangle>& AT,
                   const std::vector<math::Vector3d>& AN)
:triangles(AT), normals(AN)
{}
/*---------------------------------------------------------------------------*/
math::Point CavitySurfacePatch::project(gmds::math::Point& AP) const
{
    math::Point proj;
    double dist = 100000;
    for(auto t:triangles){
        math::Point p = t.project(AP);
        if(p.distance2(AP)<dist){
            proj = p;
            dist = p.distance2(AP);
        }
    }
    return proj;
}

/*---------------------------------------------------------------------------*/
math::Vector3d CavitySurfacePatch::normal(gmds::math::Point& AP) const
{
    math::Point proj;
    math::Vector3d n;
    double dist = 100000;
    for(auto i=0; i<triangles.size();i++){
        math::Triangle t = triangles[i];
        math::Point p = t.project(AP);
        if(p.distance2(AP)<dist){
            proj = p;
            dist = p.distance2(AP);
            n = normals[i];
        }
    }
    return  n;
}
/*---------------------------------------------------------------------------*/
CavitySurfacePaver::
CavitySurfacePaver(gmds::IGMesh* AMesh,
                   const std::vector<gmds::Node>& ABnd,
                   const std::vector<gmds::math::Triangle>& ATri,
                   const std::vector<gmds::math::Vector3d>& AN)
:m_mesh(AMesh), m_patch(ATri, AN), m_bnd(ABnd, &m_patch)
{
    m_bnd.setPatch(&m_patch);
    m_initial_nodes = ABnd;
}
/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::execute()
{
    m_patch.write();
    //======================================================================
    //At the departure, some loop "edges" have an inward vectors and some
    // other not. Looking at those having, we search two consecutives with
    // quite close ones. Indeed, those vectors are not totally reliable. We
    // can not trust them if they were defined on a point located on a geom
    // curve. Either the point was a loop-curve intersection (and so the
    // vector is oriented along the curve), or a new point due to loop
    // splitting.
    
    
    //======================================================================
    // 1. points of the boundary are flagged as end, side, corner and
    //    reversal in an usual manner
    //======================================================================
    m_bnd.initTypes();
    write();

    while(m_bnd.size()>6){
        //(1) we detect TFI structure, i.e 4-sided loops with 4 END and same
        // number of S between them
        if(meshGridPatch())
            continue;
        
        // (2) smallest portion between corner-end or end-end are filled in
        // first
        int i1, i2;
        findSmallestRowToInsert(i1,i2);
        std::cout<<"**************** ROW BETWEEN "
        <<m_bnd[i1].node.getID()<<" AND "
        <<m_bnd[i2].node.getID()
        <<"***************"<<std::endl;
        
        std::vector<Node> row_nodes = insertRow(i1,i2);
        write();
        
        smooth(row_nodes);
        write();

    }
    
    finalize();
    write();

}
/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::smooth(std::vector<Node>& ANodes)
{
    //Only nodes in ANodes are going to be smoothed
    std::vector<math::Point> update_pnts;
    update_pnts.resize(m_free_nodes.size());
    for(auto i=0; i<m_free_nodes.size(); i++){
        Node n = m_free_nodes[i];
        std::vector<Face> nf = n.get<Face>();
        std::set<TCellID> rink;
        for(auto f:nf){
            std::vector<TCellID> fn_ids = f.getIDs<Node>();
            for(auto id:fn_ids){
                if(id!=n.getID())
                    rink.insert(id);
            }
        }
        int rink_size =rink.size();
        math::Vector3d v(0,0,0);
        for(auto id:rink){
            math::Point p = m_mesh->get<Node>(id).getPoint();
            v+=math::Vector3d(p.X(),p.Y(),p.Z());
        }
        v/=rink_size;
        update_pnts[i]=math::Point(v.X(),v.Y(),v.Z());
    }
    
    
    for(auto i=0; i<m_free_nodes.size(); i++){
        m_free_nodes[i].setPoint(update_pnts[i]);

    }
}
/*---------------------------------------------------------------------------*/
Face CavitySurfacePaver::addQuad(Node& AN1, Node& AN2,
                                 Node& AN3, Node& AN4)
{
    Face f = m_mesh->newQuad(AN1,AN2,AN3,AN4);
    AN1.add(f);
    AN2.add(f);
    AN3.add(f);
    AN4.add(f);
    m_faces.push_back(f);
    return f;
    
}
/*---------------------------------------------------------------------------*/
Face CavitySurfacePaver::addQuad(TCellID& AI1, TCellID& AI2,
                                 TCellID& AI3, TCellID& AI4)
{
    Face f = m_mesh->newQuad(AI1,AI2,AI3,AI4);
    m_mesh->get<Node>(AI1).add(f);
    m_mesh->get<Node>(AI2).add(f);
    m_mesh->get<Node>(AI3).add(f);
    m_mesh->get<Node>(AI4).add(f);
    m_faces.push_back(f);
    return f;
    
}
/*---------------------------------------------------------------------------*/
Face CavitySurfacePaver::addTriangle(Node& AN1, Node& AN2,
                                     Node& AN3)
{
    Face f = m_mesh->newTriangle(AN1,AN2,AN3);
    AN1.add(f);
    AN2.add(f);
    AN3.add(f);
    m_faces.push_back(f);
    return f;

}

/*---------------------------------------------------------------------------*/
Node CavitySurfacePaver::addNode(math::Point& AP)
{
   Node n = m_mesh->newNode(m_patch.project(AP));
    m_free_nodes.push_back(n);
    
    return n;
}
/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::finalize()
{
    std::vector<BndPoint>& pts = m_bnd.points;
    if(m_bnd.size()==3){
        //=================================================
        // 3-SIDED
        //=================================================
        //put a triangle
        Node n1 = pts.begin()->node;
        Node n2 = std::next(pts.begin(),1)->node;
        Node n3 = std::next(pts.begin(),2)->node;
        addTriangle(n1,n2,n3);
    }
    else if(m_bnd.size()==4){
        //=================================================
        // 4-SIDED
        //=================================================
        //put a quad
        Node n1 = pts.begin()->node;
        Node n2 = std::next(pts.begin(),1)->node;
        Node n3 = std::next(pts.begin(),2)->node;
        Node n4 = std::next(pts.begin(),3)->node;
        addQuad(n1,n2,n3,n4);
    }
    else if(m_bnd.size()==5){
        //=================================================
        // 5-SIDED
        //=================================================
        finalize5();
    }
    else if(m_bnd.size()==6){
        //=================================================
        // 6-SIDED
        //=================================================
        finalize6();
    }
}

/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::finalize5()
{
    std::vector<BndPoint>& pts = m_bnd.points;
    
    // put a quad and a triangle depending on the best quad we
    // can generate
    Node n1 = pts[0].node;
    Node n2 = pts[1].node;
    Node n3 = pts[2].node;
    Node n4 = pts[3].node;
    Node n5 = pts[4].node;
    
    
    double q[5] = {
        quality(n1,n2,n3,n4),
        quality(n1,n2,n3,n5),
        quality(n1,n2,n4,n5),
        quality(n1,n3,n4,n5),
        quality(n2,n3,n4,n5)
    };

    //A config is valid only if the generate triangle and quad are
    // oriented in a consistent manner. Otherwise, it means
    // an element will be inverted.
    std::vector<int> valid_config;
    //------ Config 0 ------
    //Quad(n1,n2,n3,n4) + Triangle(n1,n4,n5);
    math::Vector3d u1(n1.getPoint(),n2.getPoint());
    math::Vector3d v1(n1.getPoint(),n4.getPoint());
    math::Vector3d v2(n1.getPoint(),n5.getPoint());
    if(u1.cross(v1).dot((v1.cross(v2)))>0)
        valid_config.push_back(0);
    //------ Config 1 ------
    //Quad(n1,n2,n3,n5) + Triangle(n4,n5,n3);
    u1=math::Vector3d(n5.getPoint(),n1.getPoint());
    v1=math::Vector3d(n5.getPoint(),n3.getPoint());
    v2=math::Vector3d(n5.getPoint(),n4.getPoint());
    if(u1.cross(v1).dot((v1.cross(v2)))>0)
        valid_config.push_back(1);
    
    //------ Config 2 ------
    //Quad(n1,n2,n4,n5) + Triangle(n4,n2,n3);
    u1=math::Vector3d(n4.getPoint(),n5.getPoint());
    v1=math::Vector3d(n4.getPoint(),n2.getPoint());
    v2=math::Vector3d(n4.getPoint(),n3.getPoint());
    if(u1.cross(v1).dot((v1.cross(v2)))>0)
        valid_config.push_back(2);

    //------ Config 3 ------
    //Quad(n1,n3,n4,n5) + Triangle(n3,n1,n2)
    u1=math::Vector3d(n3.getPoint(),n4.getPoint());
    v1=math::Vector3d(n3.getPoint(),n1.getPoint());
    v2=math::Vector3d(n3.getPoint(),n2.getPoint());
    if(u1.cross(v1).dot((v1.cross(v2)))>0)
        valid_config.push_back(3);

    //------ Config 4 ------
    //Quad(n2,n3,n4,n5) + Triangle(n2,n5,n1)
    u1=math::Vector3d(n2.getPoint(),n3.getPoint());
    v1=math::Vector3d(n2.getPoint(),n5.getPoint());
    v2=math::Vector3d(n2.getPoint(),n1.getPoint());
    if(u1.cross(v1).dot((v1.cross(v2)))>0)
        valid_config.push_back(4);
    
    //Now among the valid configurations, we only keep
    // the best quad quality one
    bool found_best = false;
    int best_config = -1;
    for(auto i=0;i<valid_config.size() && !found_best;i++){
        double qi = q[valid_config[i]];
        bool best = true;
        for(auto j=0;j<valid_config.size()&& best;j++){
            if(i==j)
                continue;
            if(q[valid_config[j]]>qi){
                best =false;
            }
        }
        if(best){
            found_best=true;
            best_config=valid_config[i];
        }
    }
    //the best quad is the one with the best worst angle quality
    if(best_config==0){
        addQuad(n1,n2,n3,n4);
        addTriangle(n4,n5,n1);
        
    }
    else if(best_config==1){
        addQuad(n1,n2,n3,n5);
        addTriangle(n5,n4,n3);
    }
    else if(best_config==2){
        addQuad(n1,n2,n4,n5);
        addTriangle(n4,n3,n2);
        
    }
    else if(best_config==3){
        addQuad(n1,n3,n4,n5);
        addTriangle(n3,n2,n1);
        
    }
    else if(best_config==4){
        addQuad(n2,n3,n4,n5);
        addTriangle(n2,n1,n5);
    }
    
}

/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::finalize6()
{
    std::vector<BndPoint>& p = m_bnd.points;
    
    const int nb_patterns = 6;
    BndPoint::vertex_type pattern[nb_patterns][6]={
        {// CONFIG 0
            BndPoint::END, BndPoint::END,
            BndPoint::SIDE,BndPoint::END,
            BndPoint::END, BndPoint::SIDE
        },
        {// CONFIG 1
            BndPoint::END, BndPoint::SIDE,
            BndPoint::END, BndPoint::SIDE,
            BndPoint::END, BndPoint::SIDE
        },
        {// CONFIG 2
            BndPoint::SIDE, BndPoint::END,
            BndPoint::SIDE,BndPoint::END,
            BndPoint::END, BndPoint::END
        },
        {// CONFIG 3
            BndPoint::END, BndPoint::END,
            BndPoint::SIDE,BndPoint::SIDE,
            BndPoint::END, BndPoint::END
        },
        {// CONFIG 4
            BndPoint::END, BndPoint::END,
            BndPoint::END,BndPoint::END,
            BndPoint::END, BndPoint::CORNER
        },
        {// CONFIG 5
            BndPoint::END, BndPoint::END,
            BndPoint::END,BndPoint::END,
            BndPoint::END, BndPoint::SIDE
        }
    };
    
    auto found_pattern=false;
    auto pattern_id=0;
    auto pattern_start=0;
    for(auto i=0; i<nb_patterns && !found_pattern; i++){
        if(checkPattern6(pattern[i], pattern_start)){
            found_pattern=true;
            pattern_id=i;
        }
    }
    if(!found_pattern) {
        throw GMDSException("UNknow patter for a 6-sided patch");
    }
    
    Node n0 = p[pattern_start].node;
    Node n1 = p[(pattern_start+1)%6].node;
    Node n2 = p[(pattern_start+2)%6].node;
    Node n3 = p[(pattern_start+3)%6].node;
    Node n4 = p[(pattern_start+4)%6].node;
    Node n5 = p[(pattern_start+5)%6].node;
    if(pattern_id==0){
        addQuad(n0,n1,n2,n5);
        addQuad(n2,n3,n4,n5);
    }
    else if(pattern_id==1){
        math::Point p6 = (1.0/3.0)*(n1.getPoint()+
                                    n3.getPoint()+
                                    n5.getPoint());
        
        Node n6 = addNode(p6);
        addQuad(n0,n1,n6,n5);
        addQuad(n1,n2,n3,n6);
        addQuad(n3,n4,n5,n6);
    }
    else if(pattern_id==2){
        math::Point p6 = 0.25*(n1.getPoint()+
                               n3.getPoint()+
                               n4.getPoint()+
                               n5.getPoint());
        
        Node n6 = addNode(p6);
        addQuad(n0,n1,n2,n6);
        addQuad(n2,n3,n4,n6);
        addQuad(n4,n5,n0,n6);
   
    }
    else if(pattern_id==3){
        math::Point p05 = 0.5*(n0.getPoint()+
                               n5.getPoint());
        
        math::Point p6 = 0.25*(n0.getPoint()+
                               n1.getPoint()+
                               n2.getPoint()+ p05);
        math::Point p7 = 0.25*(n3.getPoint()+
                               n4.getPoint()+
                               n5.getPoint()+ p05);
        
        
        Node n6 = addNode(p6);
        Node n7 = addNode(p7);
        addQuad(n0,n1,n2,n6);
        addQuad(n2,n3,n7,n6);
        addQuad(n3,n4,n5,n7);
        addQuad(n5,n0,n6,n7);
    }
    else if(pattern_id==4 ||
            pattern_id==5){
        addQuad(n0,n1,n2,n5);
        addQuad(n2,n3,n4,n5);
    }
    
}
/*---------------------------------------------------------------------------*/
bool CavitySurfacePaver::checkPattern6(const BndPoint::vertex_type APattern[],
                                       int& AStartingIndex)
{
    std::vector<BndPoint>& p = m_bnd.points;
    for(auto i=0;i<6;i++){
        BndPoint::vertex_type config[6]={
            p[i].type,
            p[(i+1)%6].type,
            p[(i+2)%6].type,
            p[(i+3)%6].type,
            p[(i+4)%6].type,
            p[(i+5)%6].type
        };
        bool stop=false;
        for(auto j=0;j<6 && !stop;j++){
            if(config[j]!=APattern[j]){
                stop=true;
            }
        }
        if(!stop){
            //found the right pattern
            AStartingIndex=i;
            return true;
        }
    }
    
    return false;
}
/*---------------------------------------------------------------------------*/
bool CavitySurfacePaver::meshGridPatch()
{
    int index[4];
    if(!check4Sided(index))
        return false;
    
    std::vector<BndPoint>& pnts = m_bnd.points;
    // So we have our quad structure, is it regular? Same number of edges on
    // opposite sides?
    int side1 = index[1]-index[0];
    int side2 = index[2]-index[1];
    int side3 = index[3]-index[2];
    int side4 = index[0]+(pnts.size()-index[3]);
    
    if(side1!=side3 || side2!=side4)
        return false;
    
    //Se we have a quad structure and we distinguish 3 configurations:
    // (1) we have a line of quad along side 1
    // (2) we have a line of quad along side 2
    // (3) general case where innner points must be add
    // A 1x1 patch must not be encountered in theory
    if(side1==1){
        mesh1RowPatch(index[1], index[2]);
    }
    else if (side2==1){
        mesh1RowPatch(index[0], index[1]);
    }
    else{
        //it means that we have to insert innner points
        meshTFIPatch(index);
    }
    m_bnd.clear();
    return true;
}
/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::meshTFIPatch(const int AIndex[4]){
    std::vector<BndPoint>& pnts = m_bnd.points;
    int nb_pnts = pnts.size();
    
    int i0 = AIndex[0];
    int i1 = AIndex[1];
    int i2 = AIndex[2];
    int i3 = AIndex[3];
    
    const int s1 = i1-i0+1;
    const int s2 = i2-i1+1;
    
    gmds::TCellID t[s1][s2];
    for(int i=0; i<s1; i++){
        t[i][0] = pnts[(i0+i)%nb_pnts].node.getID();
    }
    for(int i=0; i<s1; i++){
        t[s1-1-i][s2-1] = pnts[(i2+i)%nb_pnts].node.getID();
    }
    for(int i=0; i<s2; i++){
        t[0][s2-1-i] = pnts[(i3+i)%nb_pnts].node.getID();
    }
    for(int i=0; i<s2; i++){
        t[s1-1][i] = pnts[(i1+i)%nb_pnts].node.getID();
    }

    // Now we create inside points using a simple algorithm
    // (to change to TFI one)
    math::Point p00 = m_mesh->get<Node>(t[0   ][0   ]).getPoint();
    math::Point p10 = m_mesh->get<Node>(t[s1-1][0   ]).getPoint();
    math::Point p01 = m_mesh->get<Node>(t[0   ][s2-1]).getPoint();
    math::Point p11 = m_mesh->get<Node>(t[s1-1][s2-1]).getPoint();
    for(auto i=1; i<s1-1; i++){
        double ei = (double)i/(double)(s1-1);
        for(auto j=1; j<s2-1; j++){
            double ej = (double)j/(double)(s2-1);
            math::Point pij=
            (1-ej)*m_mesh->get<Node>(t[i][0]).getPoint()+
            ej*m_mesh->get<Node>(t[i][s2-1]).getPoint()+
            (1-ei)*m_mesh->get<Node>(t[0][j]).getPoint()+
            ei*m_mesh->get<Node>(t[s1-1][j]).getPoint() -
            (1-ei)*(1-ej)*p00-
            (1-ei)*ej*p01-
            ei*(1-ej)*p10-
            ei*ej*p11;
            Node nij =  addNode(pij);
            t[i][j]=nij.getID();
            m_free_nodes.push_back(nij);
        }
    }
    
    //Finally quad elements are built on the grid structure
    for(auto i=1; i<s1; i++){
        for(auto j=1; j<s2; j++){
            addQuad(t[i  ][j  ],
                    t[i-1][j  ],
                    t[i-1][j-1],
                    t[i  ][j-1]);
        }
    }
}
/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::mesh1RowPatch(const int& AI, const int& AJ){
    std::vector<BndPoint>& p = m_bnd.points;
    int nb_pnts = p.size();
    int nb_quads = AJ-AI;
    
    for(int i=0; i<nb_quads;i++){
        int i0 = (AI+i  )%nb_pnts;
        int i1 = (AI+i+1)%nb_pnts;
        int i2 = (AJ+1+nb_quads-i-1)%nb_pnts;
        int i3 = (AJ+1+nb_quads-i  )%nb_pnts;
        addQuad(p[i0].node, p[i1].node, p[i2].node, p[i3].node);
    }

}
/*---------------------------------------------------------------------------*/
bool CavitySurfacePaver::check4Sided(int (&AIndex)[4])
{
    int index=0;
    std::vector<BndPoint>& pnts = m_bnd.points;
    for(auto i=0;i<pnts.size(); i++){
        if(pnts[i].type==BndPoint::END){
            //We cannot at a fifth END corner
            if(index>=4)
                return false;
            AIndex[index++]=i;
        }
    }
    //if we have 4 corners, it's ok
    return (index==4);
}

/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::findSmallestRowToInsert(int& AI, int& AJ)
{
    std::vector<int> end_index;
    std::vector<BndPoint>& pnts = m_bnd.points;
    int nb_pnts = pnts.size();
    for(auto i=0;i<nb_pnts; i++){
        if(pnts[i].type==BndPoint::END){
            end_index.push_back(i);
            std::cout<<"\t End point "<<i<<": "
            <<pnts[i].node.getID()<<std::endl;
        }
    }

    //we keep the id that begins the smallest row
    int smallest_id=-1;
    int smallest_size=10000;
    std::vector<int> smallest_start;
    if(end_index.empty() || end_index.size()==1){
        std::cout<<"Stop meshing";
        write();
        exit(0);
    }
    
    for(int i=0;i<end_index.size()-1;i++){
        int i1 = end_index[i];
        int i2 = end_index[i+1];
        
        if(i2-i1<smallest_size){
            smallest_start.clear();
            smallest_start.push_back(i);
            smallest_size = i2-i1;
            smallest_id   = i;
        }
        else if(i2-i1==smallest_size){
            smallest_start.push_back(i);
        }
    }
    //last interval to check
    int i1 = end_index[end_index.size()-1]; //last
    int i2 = end_index[0];//first
    if(nb_pnts-i1+i2<smallest_size){
        std::cout<<" last element"<<std::endl;
        smallest_size = nb_pnts-i1+i2;
        smallest_id   = end_index.size()-1;
        smallest_start.clear();
        smallest_start.push_back(smallest_id);
    }
    else  if(nb_pnts-i1+i2==smallest_size){
        smallest_start.push_back(end_index.size()-1);
    }
    if(smallest_start.size()==1){
        std::cout<<"() SIZE= "<<smallest_size<<std::endl;

        AI = end_index[smallest_id];
        AJ = (AI+smallest_size)%nb_pnts;
    }
    else{
        //We take the configuration with the angle the closest of 90
        //degree as possible
        int index_i = end_index[smallest_start[0]];
        int index_j = end_index[(index_i+smallest_size)%nb_pnts];
        double ai = std::abs(m_bnd.angle(index_i)-90);
        double aj = std::abs(m_bnd.angle(index_j)-90);
        double ref = (ai>aj)?ai:aj;
        int ref_id = index_i;
        //we keep the row that minimize the difference to 90 degree
        for(auto s=1; s<smallest_start.size();s++){
            index_i = end_index[smallest_start[s]];
            index_j = end_index[(index_i+smallest_size)%nb_pnts];
            ai = std::abs(m_bnd.angle(index_i)-90);
            aj = std::abs(m_bnd.angle(index_j)-90);
            double v = (ai>aj)?ai:aj;
            if(v<ref){
                ref=v;
                ref_id=index_i;
            }
        }
        std::cout<<"SIZE= "<<smallest_size<<std::endl;
        AI = ref_id;
        AJ = (AI+smallest_size)%nb_pnts;
    }
}

/*---------------------------------------------------------------------------*/
std::vector<Node> CavitySurfacePaver::insertRow(const int& AI, const int& AJ)
{
    double length_factor =1.6;
    std::vector<BndPoint>& pnts = m_bnd.points;
    int nb_pnts = pnts.size();
    
    int i = AI;
    int j = AJ;
    std::vector<Node> new_nodes;
    std::vector<math::Vector3d> new_dir;
    int row_size = (i>j)?(nb_pnts-i+j):(j-i);
    
    BndPoint prev2 = m_bnd.prev(i-1);
    BndPoint next2 = m_bnd.next(j+1);
    BndPoint prev = m_bnd.prev(i);
    BndPoint next = m_bnd.next(j);
    BndPoint ni = m_bnd[i];
    BndPoint nj = m_bnd[j];
    if(row_size==1){
        //We check if the length of the edge to insert is not to long
        //comparing to local edges
        bool insert_point = false;
        math::Point p0 =prev.node.getPoint();
        math::Point p1 =ni.node.getPoint();
        math::Point p2 =nj.node.getPoint();
        math::Point p3 =next.node.getPoint();
        double new_edge_length= math::Vector3d(p0,p3).norm();
        double e01= math::Vector3d(p0,p1).norm();
        double e23= math::Vector3d(p2,p3).norm();
        if(new_edge_length>e01*length_factor &&
           new_edge_length>e23*length_factor ){
            insert_point = true;
        }
        if(insert_point){
            //we split the biggest angle between in node ni or nj
            math::Vector3d v10(p1,p0);
            math::Vector3d v12(p1,p2);
            math::Vector3d v21(p2,p1);
            math::Vector3d v23(p2,p3);
            double angle_i = v10.angle(v12);
            double angle_j = v21.angle(v23);
            Node n_ref  = ni.node;
            Node n_next = nj.node;
            Node n_prev = prev.node;
            Node n_next2= next.node;
            Node n_prev2 = prev2.node;
            BndPoint::vertex_type prev_type = prev.type;
            math::Vector3d prev_inward = prev.inward;

            bool inverse_config=false;
            if(angle_j>angle_i){
                //We inverse the configuration
                inverse_config=true;
                n_ref = nj.node;
                n_next = ni.node;
                n_prev = next.node;
                prev_type = next.type;
                prev_inward = next.inward;

                n_next2 = prev.node;
                n_prev2 = next2.node;
            }
            //=======================
            //create point from prev
            math::Vector3d u1(n_prev.getPoint(), n_prev2.getPoint());
            math::Vector3d v1(n_prev.getPoint(), n_ref.getPoint());
            double length1 = 0.5*(u1.norm()+v1.norm());
            u1.normalize();
            v1.normalize();
            math::Vector3d dir1 = 0.5*(u1+v1);
            if(prev_type==BndPoint::CORNER ||
               prev_type==BndPoint::REVERSAL)
                dir1 = -dir1;
            else if(prev_type==BndPoint::SIDE && dir1.dot(prev_inward)<0){
                dir1 = -dir1;
            }
            dir1.normalize();
            dir1 *=length1;
            math::Point p1(n_prev.getPoint().X()+dir1.X(),
                           n_prev.getPoint().Y()+dir1.Y(),
                           n_prev.getPoint().Z()+dir1.Z());

            Node n1 = addNode(p1);
            
            //=======================
            //create point from prev
            math::Vector3d u2(n_ref.getPoint(), n_prev.getPoint());
            math::Vector3d v2(n_ref.getPoint(), n_next.getPoint());
            double length2= 0.5*(u2.norm()+v2.norm());
            u2.normalize();
            v2.normalize();
            math::Vector3d dir2 = 0.5*u2+0.5*v2;
            dir2.normalize();
            dir2 *=length2;
            math::Point p2(n_ref.getPoint().X()+dir2.X(),
                           n_ref.getPoint().Y()+dir2.Y(),
                           n_ref.getPoint().Z()+dir2.Z());
            Node n2 = addNode(p2);

            Face f1 = addQuad(n_prev,n1,n2,n_ref);
            Face f2 = addQuad(n_ref,n2,n_next2,n_next);
            
            math::Point x1 = 0.5*(n_prev.getPoint()+
                                  n1.getPoint());
            math::Point x2 = 0.5*(n1.getPoint()+
                                  n2.getPoint());
            math::Point x3 = 0.5*(n2.getPoint()+
                                  n_next2.getPoint());
            if(inverse_config){
                new_nodes.push_back(n2);
                new_nodes.push_back(n1);
                
              
                new_dir.push_back(math::Vector3d(f2.center(),x3));
                new_dir.push_back(math::Vector3d(f1.center(),x2));
                new_dir.push_back(math::Vector3d(f1.center(),x1));
            }
            else{
                new_nodes.push_back(n1);
                new_nodes.push_back(n2);
                new_dir.push_back(math::Vector3d(f1.center(),x1));
                new_dir.push_back(math::Vector3d(f1.center(),x2));
                new_dir.push_back(math::Vector3d(f2.center(),x3));
            }
            
        }
        else{
            //only one quad to create
            Face f = addQuad(prev.node,ni.node, nj.node,
                             next.node);
            
            math::Point x = 0.5*(prev.node.getPoint()+
                                 next.node.getPoint());
            
            new_dir.push_back(math::Vector3d(f.center(),x));
        }
    }
    else {
        math::Point p_prev = prev.node.getPoint();
        math::Point p_next = next.node.getPoint();
        // some inner points must be inserted too
        Node n0 = prev.node;

        for(auto x=i;x<i+row_size;x++){
            
            Node n1 =m_bnd[x].node;
            Node n2 =m_bnd[x+1].node;
            Node n3;
            if(x==i+row_size-1){
                n3 = next.node;
                Face f = addQuad(n0,n1,n2,n3);
                math::Point p03 = 0.5*(n0.getPoint() + n3.getPoint());
                new_dir.push_back(math::Vector3d(f.center(),p03));
                //we update n0 for the next quad building
            }
            else if(m_bnd[x+1].type==BndPoint::SIDE){
                math::Point p3 = createRowSideNode(x+1);
                n3 = addNode(p3);
                new_nodes.push_back(n3);
                Face f = addQuad(n0,n1,n2,n3);
                math::Point p03 = 0.5*(n0.getPoint() + n3.getPoint());
                new_dir.push_back(math::Vector3d(f.center(),p03));
                //we update n0 for the next quad building
                n0 = n3;
            }
            else if (m_bnd[x+1].type==BndPoint::CORNER){
                std::vector<math::Point> pjkl = createRowCornerNode(x+1);
                Node nj = addNode(pjkl[0]);
                Node nk = addNode(pjkl[1]);
                Node nl = addNode(pjkl[2]);
                new_nodes.push_back(nj);
                new_nodes.push_back(nk);
                new_nodes.push_back(nl);
                
                Face f = addQuad(n0,n1,n2,nl);
                math::Point p03 = 0.5*(n0.getPoint() + nl.getPoint());
                new_dir.push_back(math::Vector3d(f.center(),p03));
                f = addQuad(n2,nl,nk,nj);
                math::Point pjk = 0.5*(nj.getPoint() + nk.getPoint());
                new_dir.push_back(math::Vector3d(f.center(),pjk));
                math::Point pkl = 0.5*(nk.getPoint() + nl.getPoint());
                new_dir.push_back(math::Vector3d(f.center(),pkl));
                //we update n0 for the next quad building
                n0 = nl;
            }
            else{
                throw GMDSException("Row point type not yet implemented");
            }
            
       
        }

    }
    
    //the front must be updated
    m_bnd.replace(i,j, new_nodes, new_dir);
    //points are inserted and angles recomputed for altered points
    return new_nodes;
}
/*---------------------------------------------------------------------------*/
double CavitySurfacePaver::quality(const Node& AN1, const Node& AN2,
                                   const Node& AN3, const Node& AN4)
{
    math::Vector v12(AN1.getPoint(),AN2.getPoint());
    math::Vector v23(AN2.getPoint(),AN3.getPoint());
    math::Vector v34(AN3.getPoint(),AN4.getPoint());
    math::Vector v41(AN4.getPoint(),AN1.getPoint());
    
    double a[4] ={
        v12.angle(v41.opp()),
        v23.angle(v12.opp()),
        v34.angle(v23.opp()),
        v41.angle(v34.opp())
    };
    double q = a[0];
    for(auto i=1;i<4;i++){
        auto qi =a[i];
        if(qi<q){
            q=qi;
        }
    }
    
    return q;
}

/*---------------------------------------------------------------------------*/
void CavitySurfacePaver::write()
{
    static int i=0;
    std::string fileN = "paver_"+std::to_string(i);
    std::cout<<"***** WRITE "<<fileN<<std::endl;
    IGMesh m(MeshModel(DIM3 | F | N | F2N));
    
    Variable<math::Vector>* v;
    try{
        
        v= m_mesh->newVariable<math::Vector>(GMDS_NODE,"inward");
    }
    catch(GMDSException& e){
        v= m_mesh->getVariable<math::Vector>(GMDS_NODE,"inward");
    }
    
    for(int i=0; i<m_bnd.size();i++){
        BndPoint pi = m_bnd[i];
        math::Vector vi(pi.inward.X(),
                        pi.inward.Y(),
                        pi.inward.Z());
        std::cout<<"var for "<<pi.node.getID()<<": "<<vi<<std::endl;

        (*v)[pi.node.getID()]=vi;
    }
    VTKWriter<IGMesh> w(*m_mesh);
    w.write(fileN, DIM3 | F | N);
    i++;
    
}

/*---------------------------------------------------------------------------*/
math::Point CavitySurfacePaver::createRowSideNode(const int& AI)
{

    BndPoint cur = m_bnd[AI];
    BndPoint pre = m_bnd.prev(AI);
    BndPoint nex = m_bnd.next(AI);
    std::cout<<"========== SIDE NODE "<<cur.node.getID()
    <<" ==========="<<std::endl;
    
    math::Point c = cur.node.getPoint();
    math::Point p = pre.node.getPoint();
    math::Point n = nex.node.getPoint();

    double angle_deg =  m_bnd.angle(AI);
    if(angle_deg<0)
        angle_deg= 360-angle_deg;
    double angle =angle_deg*math::Constants::PIDIV180;
    std::cout<<"ANGLE = "<<angle_deg<<" ("<<angle<<")"<<std::endl;
    
    math::Vector3d in = pre.inward;
    math::Vector3d cn(c,n);
    math::Vector3d cp(c,p);
    std::cout<<"IN "<<in<<std::endl;
    std::cout<<"cn.norm = "<<cn.norm()<<std::endl;
    std::cout<<"cp.norm = "<<cp.norm()<<std::endl;
//    cn.normalize();
//    in.normalize();
    
    math::Vector3d normal = cn.cross(in);
    normal.normalize();
    std::cout<<"Normal "<<normal<<std::endl;
    
    math::AxisAngleRotation rot(normal,angle/2.0);
    std::cout<<"Rotation d'angle "<<(angle/2.0)*math::Constants::INVPIDIV180<<std::endl;;
    math::Vector3d v = rot*cn;
    
    std::cout<<"cn = "<<cn<<std::endl;
    std::cout<<"v  = "<<v<<std::endl;
    v.normalize();
    
    
    std::cout<<"ROW SIDE NODE VECTOR: "<<v<<std::endl;
    double v_length = 0.5*(cp.norm()+cn.norm());
    std::cout<<"length: "<<v_length<<std::endl;
    v *=v_length;
    std::cout<<"ROW SIDE NODE VECTOR: "<<v<<std::endl;;
    return math::Point(c.X()+v.X(),c.Y()+v.Y(),c.Z()+v.Z());

}
/*---------------------------------------------------------------------------*/
std::vector<math::Point>
CavitySurfacePaver::createRowCornerNode(const int& AI)
{
    BndPoint cur = m_bnd[AI];
    BndPoint pre = m_bnd.prev(AI);
    BndPoint nex = m_bnd.next(AI);
    
    
    std::cout<<"========== CORNER NODE "<<cur.node.getID()
    <<" ==========="<<std::endl;
    
    math::Point c = cur.node.getPoint();
    math::Point p = pre.node.getPoint();
    math::Point n = nex.node.getPoint();
    
    double angle_deg =  m_bnd.angle(AI);
    if(angle_deg<0)
        angle_deg= 360-angle_deg;
    double angle =angle_deg*math::Constants::PIDIV180;
    std::cout<<"ANGLE = "<<angle_deg<<" ("<<angle<<")"<<std::endl;
    
    math::Vector3d in = cur.inward;
    math::Vector3d cn(c,n);
    math::Vector3d cp(c,p);
    std::cout<<"IN "<<in<<std::endl;
    std::cout<<"cn.norm = "<<cn.norm()<<std::endl;
    std::cout<<"cp.norm = "<<cp.norm()<<std::endl;
    //    cn.normalize();
    //    in.normalize();
    
    math::Vector3d normal = cn.cross(in);
    normal.normalize();
    std::cout<<"Normal "<<normal<<std::endl;

    
    
    math::Vector3d vj = pre.inward;
    math::Vector3d vl = cur.inward;
    math::Vector3d vk = vj+vl;
    vj.normalize();
    vk.normalize();
    vl.normalize();
    
    
    double jl_length = 0.5*(cp.norm()+cn.norm());
    vj *=jl_length;
    vl *=jl_length;
    vk *=sqrt(2.0)*jl_length;
    
    math::Point pj(c.X()+vj.X(),c.Y()+vj.Y(),c.Z()+vj.Z());
    math::Point pk(c.X()+vk.X(),c.Y()+vk.Y(),c.Z()+vk.Z());
    math::Point pl(c.X()+vl.X(),c.Y()+vl.Y(),c.Z()+vl.Z());

    std::vector<math::Point> pjkl;
    pjkl.resize(3);
    pjkl[0] =pj;
    pjkl[1] =pk;
    pjkl[2] =pl;
    
    return pjkl;
}
/*---------------------------------------------------------------------------*/
std::vector<math::Point>
CavitySurfacePaver::createRowReversalNode(const int& AI)
{
    std::vector<math::Point> ns;
    return ns;
}
/*---------------------------------------------------------------------------*/
void CavitySurfacePatch::write()
{
    static int i=0;
    std::string fileN = "surface_patch_"+std::to_string(i);
    IGMesh m(MeshModel(DIM3 | F | N | F2N));
    for(auto t:triangles){
        Node n0 = m.newNode(t.getPoint(0));
        Node n1 = m.newNode(t.getPoint(1));
        Node n2 = m.newNode(t.getPoint(2));
        m.newTriangle(n0,n1,n2);
    }
    
    VTKWriter<IGMesh> w(m);
    w.write(fileN, DIM3 | F | N);
    i++;
    
}
/*---------------------------------------------------------------------------*/


