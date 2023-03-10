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
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include "EigenFieldSolverStrategy.h"
/*---------------------------------------------------------------------------*/
using namespace fhedo;
/*---------------------------------------------------------------------------*/
EigenFieldSolverStrategy::EigenFieldSolverStrategy(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::solve(){
    //    static int compteur=0;
    //    int nb_lin = m_nb_equations;
    //    int nb_col = m_nb_unknowns;
    //    //    std::cout<<m_X<<std::endl;
    //    //    std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    //    //    std::cout<<m_A<<std::endl;
    //    //    std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    //    //    std::cout<<m_b<<std::endl;
    //    //    std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    //    Eigen::SparseMatrix<double> AT (nb_col,nb_lin);
    //    Eigen::SparseMatrix<double> ATA(nb_lin,nb_lin);
    //    AT= m_A.transpose();
    //    AT.makeCompressed();
    //    ATA=AT*m_A;
    //    //We finalize the Sparse Matrix DS
    //    ATA.makeCompressed();
    //
    //    std::cout<<"System solving for initialization"<<std::endl;
    //    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //    //    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver;
    //    solver.compute(ATA);
    //    Eigen::VectorXd prev_X = m_X;
    //
    //    solver.setMaxIterations(1);
    //    if(compteur==0)
    //        m_X = solver.solve(AT*m_b);
    //    else
    //        m_X = solver.solveWithGuess(AT*m_b, prev_X);
    //    std::cout<<"Movement: "<<(m_X-prev_X).norm()<<std::endl;
    //    std::cout<<" Solved in nb. iterations= "<<solver.iterations()<<", with error "<<solver.error()<<std::endl;
    //    //  std::cout<<"GET"<<m_X<<std::endl;
    //    compteur=1;

}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::initializeAssembly(){
    
    //We initialize the equation numbering variable
    //m_equation_numbering=0;
    
    
    //        m_A = Eigen::SparseMatrix<double>(m_nb_equations,m_nb_unknowns);
    //
    //        if(i==0)
    //            m_X = Eigen::VectorXd(m_nb_unknowns);
    //        else if (i==1){
    //            Eigen::VectorXd old_X = m_X;
    //            m_X = Eigen::VectorXd(m_nb_unknowns);
    //            for(int i=0; i <9*m_nb_free_nodes+ 2*m_nb_surface_nodes ; i++){
    //                m_X[i]=old_X[i];
    //            }
    //            for(int i=9*m_nb_free_nodes+ 2*m_nb_surface_nodes; i<m_nb_unknowns ; i++){
    //                m_X[i]=0;
    //            }
    //        }
    //
    //        m_b = Eigen::VectorXd(m_nb_equations);

}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::finalizeAssembly(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::clean(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::setX(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::getFeasibleSolution(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::addLocalOptimConstraints(){
    //
    //
    //            m_A.coeffRef(m_equation_numbering, 9*i+k) = m_lambda;
    //            m_A.coeffRef(m_equation_numbering, step+3*i  ) = -m_lambda*cx[k];
    //            m_A.coeffRef(m_equation_numbering, step+3*i+1) = -m_lambda*cy[k];
    //            m_A.coeffRef(m_equation_numbering, step+3*i+2) = -m_lambda*cz[k];
    //            m_b[m_equation_numbering]=m_lambda*ai[k];
    //            m_equation_numbering++;//next line

}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::addLockedTerms(){
    
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::addBoundaryConstraints(){
//    for(unsigned int i_node=0; i_node<m_surface_nodes.size(); i_node++){
//        Node ni = m_surface_nodes[i_node];
//        
//        TCellID i   = ni.getID();
//        int local_i = (*m_ordering)[i];
//        
//        math::SHarmonicL4 h0 = m_H0[i];
//        math::SHarmonicL4 h4 = m_H4[i];
//        math::SHarmonicL4 h8 = m_H8[i];
//        //We create the smoothing constraint on the 9D rep vectors of i and j
//        
//        for(int k=0;k<9;k++){
//            
//            // WITH EIGEN
//            //        m_A.coeffRef(m_equation_numbering, 9*local_i+k) = m_lambda;
//            //        m_A.coeffRef(m_equation_numbering,
//            //                     9*m_nb_free_nodes+2*local_i  ) = -m_lambda*h8[k];
//            //        m_A.coeffRef(m_equation_numbering,
//            //                     9*m_nb_free_nodes+2*local_i+1) = -m_lambda*h0[k];
//            //        m_b[m_equation_numbering]=m_lambda*sqrt(7.0/12.0)*h4[k];
//            //        m_equation_numbering++;
//            
//            nlBegin(NL_ROW);
//            nlCoefficient(9*local_i+k                  ,  m_lambda);
//            nlCoefficient(9*m_nb_free_nodes+2*local_i  , -m_lambda*h8[k]);
//            nlCoefficient(9*m_nb_free_nodes+2*local_i+1, -m_lambda*h0[k]);
//            nlRightHandSide(m_lambda*sqrt(7.0/12.0)*h4[k]);
//            nlEnd(NL_ROW);
//            
//        }
//    }
}
/*---------------------------------------------------------------------------*/
void EigenFieldSolverStrategy::addSmoothingTerms(){
    //
    //        std::vector<Node> e_nodes = e.get<Node>();
    //        if(!m_mesh->isMarked(e_nodes[0], m_markNodeFixed)&&
    //           !m_mesh->isMarked(e_nodes[1], m_markNodeFixed) )
    //        {
    //            int i = m_order_nodes[e_nodes[0].getID()];
    //            int j = m_order_nodes[e_nodes[1].getID()];
    //            //We create the smoothing constraint on the 9D rep vectors of i and j
    //            for(int k=0;k<9;k++){
    ////                m_A.coeffRef(m_equation_numbering+k, 9*i+k) = 1;
    ////                m_A.coeffRef(m_equation_numbering+k, 9*j+k) = -1;
    ////                m_b[m_equation_numbering+k]=0;
    ////                m_equation_numbering++;//next line
    //
    //                nlBegin(NL_ROW);
    //                nlCoefficient(9*i+k, 1.0);
    //                nlCoefficient(9*j+k,-1.0);
    //                nlRightHandSide(0);
    //                nlEnd(NL_ROW);
    //            }
    //        }
    //        else if( m_mesh->isMarked(e_nodes[0], m_markNodeFixed)&&
    //                (!m_mesh->isMarked(e_nodes[1], m_markNodeFixed)) ){
    //            int i = m_order_nodes[e_nodes[1].getID()];
    //            math::SHFrame fixed_f = (*m_frame_field)[e_nodes[0].getID()];
    //            math::SHVector a = fixed_f.a();
    //            for(int k=0;k<9;k++){
    ////                m_A.coeffRef(m_equation_numbering+k, 9*i+k) = 1;
    ////                m_b[m_equation_numbering+k]=a[k];
    ////                m_equation_numbering++;//next line
    //
    //
    //                nlBegin(NL_ROW);
    //                nlCoefficient(9*i+k, 1.0);
    //                nlRightHandSide(a[k]);
    //                nlEnd(NL_ROW);
    //
    //            }
    //        }
    //        else if((!m_mesh->isMarked(e_nodes[0], m_markNodeFixed))&&
    //                 m_mesh->isMarked(e_nodes[1], m_markNodeFixed) ){
    //            int i = m_order_nodes[e_nodes[0].getID()];
    //            math::SHFrame fixed_f = (*m_frame_field)[e_nodes[1].getID()];
    //            math::SHVector a = fixed_f.a();
    //            for(int k=0;k<9;k++){
    //                //                m_A.coeffRef(m_equation_numbering+k, 9*i+k) = 1;
    //                //                m_b[m_equation_numbering+k]=a[k];
    //                //                m_equation_numbering++;//next line
    //                nlBegin(NL_ROW);
    //                nlCoefficient(9*i+k, 1.0);
    //                nlRightHandSide(a[k]);
    //                nlEnd(NL_ROW);
    //
    //            }
    //        }

}
/*----------------------------------------------------------------------------*/