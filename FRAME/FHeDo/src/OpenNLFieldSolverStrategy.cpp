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
// STL File Headers
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// OpenNL File Headers
#include "OpenNL_psm.h"
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include "OpenNLFieldSolverStrategy.h"
#include "GMDS/Math/SHarmonicL4.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
/*---------------------------------------------------------------------------*/
OpenNLFieldSolverStrategy::OpenNLFieldSolverStrategy(){}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::solve(){
    /* Solve and get solution */
    nlSolve();
    NLint    nb_iter;
    NLdouble elapsed_time;
    
    nlGetIntegerv(NL_USED_ITERATIONS, &nb_iter);
    nlGetDoublev (NL_ELAPSED_TIME   , &elapsed_time);
    
    std::cout<<"... elapsed time "<<elapsed_time
    <<" in "<<nb_iter<<" iterations"<<std::endl;

}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::initializeAssembly(){
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, m_nb_unknowns);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::finalizeAssembly(){
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::clean()
{
    nlDeleteContext(nlGetCurrent());
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::setX()
{
    IGMesh::node_iterator it_node = m_mesh->nodes_begin();
    for (; !it_node.isDone(); it_node.next()) {
        Node n = it_node.value();
        int i = (*m_ordering)[n.getID()];
        math::SHarmonicL4 shi = (*m_harmonic_field)[n.getID()];
        for(int k=0;k<9;k++)
            nlSetVariable(9*i+k, shi[k]);
    }
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::getFeasibleSolution()
{
    IGMesh::node_iterator it_node = m_mesh->nodes_begin();
    for (; !it_node.isDone(); it_node.next())
    {
        Node n = it_node.value();
        
        if(!m_mesh->isMarked(n, m_markNodeLocked)){
            
            int local_id = 9*(*m_ordering)[n.getID()];
            double t[9] = {
                nlGetVariable(local_id),
                nlGetVariable(local_id+1),
                nlGetVariable(local_id+2),
                nlGetVariable(local_id+3),
                nlGetVariable(local_id+4),
                nlGetVariable(local_id+5),
                nlGetVariable(local_id+6),
                nlGetVariable(local_id+7),
                nlGetVariable(local_id+8)};

            math::Vector9d v(t);

            math::AxisAngleRotation r;
            
            math::SHarmonicL4 sh = math::SHarmonicL4::closest(v,r);
            
            math::Vector3d vx(1,0,0);
            math::Vector3d vy(0,1,0);
            math::Vector3d vz(0,0,1);
            vx= r*vx;
            vy= r*vy;
            vz= r*vz;
            
            (*m_harmonic_field)[n.getID()]=sh;
            (*m_chart_field)[n.getID()]=math::Chart(math::Vector(vx[0],vx[1],vx[2]),
                                                    math::Vector(vy[0],vy[1],vy[2]),
                                                    math::Vector(vz[0],vz[1],vz[2]));
            
        }
    }
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addLocalOptimConstraints(){
    int step = 9*m_nb_free_nodes+2*m_surface_nodes->size();
    
    for (IGMesh::node_iterator it_n = m_mesh->nodes_begin();
         !it_n.isDone(); it_n.next())
    {
        Node n = it_n.value();
        
        //      See what to do with fixed nodes
        //        if(m_mesh->isMarked(n,m_markNodeFixed))
        //            continue;
        
        TCellID n_id = n.getID();
        int i = (*m_ordering)[n_id];
        
        math::SHarmonicL4 prevSH = (*m_harmonic_field)[n_id];
        math::SHarmonicL4 ai = prevSH;
        math::SHarmonicL4 cx = math::EX()*ai;
        math::SHarmonicL4 cy = math::EY()*ai;
        math::SHarmonicL4 cz = math::EZ()*ai;
        //We create the smoothing constraint on the 9D rep vectors of i and j
        for(int k=0;k<9;k++){
            nlRowScaling(m_penalty);
            nlBegin(NL_ROW);
            nlCoefficient(9*i+k, 1.0);
            nlCoefficient(step+3*i  , -cx[k]);
            nlCoefficient(step+3*i+1, -cy[k]);
            nlCoefficient(step+3*i+2, -cz[k]); 
            nlRightHandSide(ai[k]);
            nlEnd(NL_ROW);
        }
    }
    
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addLockedTerms(){
    
    for (IGMesh::node_iterator it_n = m_mesh->nodes_begin();
         !it_n.isDone(); it_n.next()) {
        Node ni= it_n.value();
        if(m_mesh->isMarked(ni,m_markNodeLocked )){
            //Variable to be locked
            TCellID i   = ni.getID();
            int local_i = (*m_ordering)[i];
            
            math::SHarmonicL4 sh = (*m_harmonic_field)[ni.getID()];
            for(int k=0;k<9;k++){
                nlSetVariable(9*local_i+k, sh[k]);
                nlLockVariable(9*local_i+k);
            }
            
        }//if(m_mesh->isMarked(ni,m_markNodeFixed ))
        
    }//for (IGMesh::node_iterator it_n = m_mesh->nodes_begin()....

    nlBegin(NL_MATRIX);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addBoundaryConstraints(){
    for(unsigned int i_node=0; i_node<m_surface_nodes->size(); i_node++){
        Node ni = (*m_surface_nodes)[i_node];
        
        TCellID i   = ni.getID();
        int local_i = (*m_ordering)[i];
        
        math::SHarmonicL4 h0 = (*m_H0)[i];
        math::SHarmonicL4 h4 = (*m_H4)[i];
        math::SHarmonicL4 h8 = (*m_H8)[i];
        //We create the smoothing constraint on the 9D rep vectors of i and j
        
        for(int k=0;k<9;k++){
            
            nlBegin(NL_ROW);
            nlCoefficient(9*local_i+k                  ,  m_penalty);
            nlCoefficient(9*m_nb_free_nodes+2*local_i  , -m_penalty*h8[k]);
            nlCoefficient(9*m_nb_free_nodes+2*local_i+1, -m_penalty*h0[k]);
            nlRightHandSide(m_penalty*sqrt(7.0/12.0)*h4[k]);
            nlEnd(NL_ROW);
            
        }
    }
    
//    std::cout<<"Alignment constraint!!!"<<std::endl;
//    for(IGMesh::node_iterator it_n = m_mesh->nodes_begin();
//        !it_n.isDone(); it_n.next())
//    {
//        Node ni = it_n.value();
//        int local_i = (*m_ordering)[ni.getID()];
//        
//        math::Point pi = ni.getPoint();
//        if((!m_mesh->isMarked(ni, m_markNodeLocked)) &&
//           (!m_mesh->isMarked(ni, m_markNodeOnSurf)))
//        {
//            
//            //inner node so
//            Node closest_bnd_node =(*m_surface_nodes)[0];
//            double closest_dist = pi.distance(closest_bnd_node.getPoint());
//            
//            for(unsigned int i=1; i<m_surface_nodes->size(); i++){
//                Node bndi = (*m_surface_nodes)[i];
//                double di = pi.distance(bndi.getPoint());
//                if(di<closest_dist){
//                    closest_dist=di;
//                    closest_bnd_node=bndi;
//                }
//                
//            }//for(unsigned int i=1; i<m_surface_nodes->size(); i++)
//            
//            //We have the closest boundary node, we constraint the local
//            //frame onto the surface normal defined in that bnd node.
//            math::SHarmonicL4 h0 = (*m_H0)[closest_bnd_node.getID()];
//            math::SHarmonicL4 h4 = (*m_H4)[closest_bnd_node.getID()];
//            math::SHarmonicL4 h8 = (*m_H8)[closest_bnd_node.getID()];
//            //We create the smoothing constraint on the 9D rep vectors of i and j
//            
//            for(int k=0;k<9;k++){
//                
//                nlBegin(NL_ROW);
//                nlCoefficient(9*local_i+k                  ,  1);
//                nlCoefficient(9*m_nb_free_nodes+2*local_i  , -h8[k]);
//                nlCoefficient(9*m_nb_free_nodes+2*local_i+1, -h0[k]);
//                nlRightHandSide(sqrt(7.0/12.0)*h4[k]);
//                nlEnd(NL_ROW);
//                
//            }
//        }
//        
//    }//for(IGMesh::node_iterator it_n = m_mesh....
//    std::cout<<"\t > Initialized...OK"<<std::endl;
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addSmoothingTerms(){
    
    Variable<double>* cot_w = m_mesh->getVariable<double>(GMDS_EDGE,"cot_weight");
    

    for (IGMesh::edge_iterator it_e = m_mesh->edges_begin();
         !it_e.isDone(); it_e.next())
    {
        Edge e = it_e.value();
        std::vector<Node> e_nodes = e.get<Node>();
        Node ni = e_nodes[0];
        Node nj = e_nodes[1];
        int i = (*m_ordering)[ni.getID()];
        int j = (*m_ordering)[nj.getID()];

        double c=(*cot_w)[e.getID()];
        
        for(int k=0;k<9;k++){
            nlBegin(NL_ROW);
            nlCoefficient(9*i+k, c);
            nlCoefficient(9*j+k,-c);
            nlRightHandSide(0);
            nlEnd(NL_ROW);
        }
    }
}
/*----------------------------------------------------------------------------*/