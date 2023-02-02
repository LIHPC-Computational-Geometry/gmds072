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
// FRAME File Headers
#include "CMLSurfEvalImpl.h"
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Ray.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
CMLSurfEvalImpl::CMLSurfEvalImpl(std::vector<Face>& AFaces,
                                 std::vector<math::Vector>& ANormals)
:m_tri(AFaces),m_normal(ANormals)
{};
/*---------------------------------------------------------------------------*/
CMLSurfEvalImpl::~CMLSurfEvalImpl()
{};
/*---------------------------------------------------------------------------*/
double CMLSurfEvalImpl::area()
{
    double a=0;
    for(auto t:m_tri){
        a+=t.area();
    }
    return a;
}
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::bounding_box(double min_xyz[3], double max_xyz[3])
{
    min_xyz[0]=LONG_MAX;
    min_xyz[1]=LONG_MAX;
    min_xyz[2]=LONG_MAX;
    
    max_xyz[0]=LONG_MIN;
    max_xyz[1]=LONG_MIN;
    max_xyz[2]=LONG_MIN;
    for(auto t:m_tri){
        std::vector<Node> n = t.get<Node>();
        for(auto ni:n){
            math::Point pi = ni.getPoint();
            
            // X COORDINATE
            if(pi.X()<min_xyz[0]){
                min_xyz[0]=pi.X();
            }
            else if(pi.X()>max_xyz[0]){
                max_xyz[0]=pi.X();
            }
            
            // Y COORDINATE
            if(pi.Y()<min_xyz[1]){
                min_xyz[1]=pi.Y();
            }
            else if(pi.Y()>max_xyz[1]){
                max_xyz[1]=pi.Y();
            }
            
            // Z COORDINATE
            if(pi.Z()<min_xyz[2]){
                min_xyz[2]=pi.Z();
            }
            else if(pi.Z()>max_xyz[2]){
                max_xyz[2]=pi.Z();
            }
        }
    }
}


//// FRANCK, since we gave you CAMAL, we renamed move_to_surface to move_to.  
////// Our CMLPave calls move_to, and your version of CMLPave calls move_to_surface,
////// So this change defines both so it will work for both of us.
////// -Matt Staten
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::move_to(double& x, double& y, double& z )
{
  move_to_surface(x,y,z);
}

/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::move_to_surface(double& x, double& y, double& z )
{
    math::Point p(x,y,z);
    std::vector<Node> cur_nodes = m_tri[0].get<Node>();
    math::Triangle t(cur_nodes[0].getPoint(),
                     cur_nodes[1].getPoint(),
                     cur_nodes[2].getPoint());
    math::Point proj_p = t.project(p);
    double min_dist = p.distance2(proj_p);
    
    for(auto i=1; i<m_tri.size();i++){
        cur_nodes = m_tri[i].get<Node>();
        math::Triangle ti(cur_nodes[0].getPoint(),
                          cur_nodes[1].getPoint(),
                          cur_nodes[2].getPoint());
        try{
            math::Point pi = ti.project(p);
            double di = p.distance2(pi);
            if(di<min_dist) {
                min_dist = di;
                proj_p   = pi;
            }
        }
        catch(GMDSException& e)
        {;}
        
    }
   // std::cout<<"Move to surface  ("<<x<<", "<<y<<", "<<z<<"): "<<proj_p<<std::endl;

    x=proj_p.X();
    y=proj_p.Y();
    z=proj_p.Z();
    
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::
move_to_surface_along_vector(double& x, double& y, double& z,
                             double &along_x, double &along_y, double &along_z)
{
    math::Point  from(x,y,z);
    math::Vector dir (along_x,along_y,along_z);
    math::Ray    ray (from,dir);
    
    for(auto i=0; i<m_tri.size();i++){
        std::vector<Node> cur_nodes = m_tri[i].get<Node>();
        math::Triangle ti(cur_nodes[0].getPoint(),
                          cur_nodes[1].getPoint(),
                          cur_nodes[2].getPoint());
        math::Point intersection;

        if(!ray.intersect3D(ti, intersection))
            continue;
        
        //We intersect
        x = intersection.X();
        y = intersection.Y();
        z = intersection.Z();
        return true;
    }
    return false;
}
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::move_to_surface(double& x, double& y, double& z,
                                      double& u_guess, double& v_guess)
{
    throw GMDSException("CMLSurfEvalImpl::move_to_surface() not implemented");
}
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::move_to_surface(double& x, double& y, double& z,
                                      int &vert_idx )
{
    move_to_surface(x,y,z);
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::normal_at(double x, double y, double z,
                                double& nx, double& ny, double& nz)
{
    math::Point p(x,y,z);
    int index = 0;
    std::vector<Node> cur_nodes = m_tri[0].get<Node>();
    math::Triangle t(cur_nodes[0].getPoint(),
                     cur_nodes[1].getPoint(),
                     cur_nodes[2].getPoint());
    math::Point proj_p = t.project(p);
    double min_dist = p.distance2(proj_p);
    
    for(auto i=1; i<m_tri.size();i++){
        cur_nodes = m_tri[i].get<Node>();
        math::Triangle ti(cur_nodes[0].getPoint(),
                          cur_nodes[1].getPoint(),
                          cur_nodes[2].getPoint());
        math::Point pi = t.project(p);
        double di = p.distance2(pi);
        if(di<min_dist) {
            min_dist = di;
            proj_p   = pi;
            index =i;
        }
        
    }
    
    math::Vector n = m_normal[index];
    nx=n.X();
    ny=n.Y();
    nz=n.Z();
  //  std::cout<<"normal at ("<<x<<", "<<y<<", "<<z<<"): "<<n<<std::endl;
    return true;
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::normal_at(double x, double y, double z,
                                double& u_guess, double& v_guess,
                                double& nx, double& ny, double& nz)
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::normal_at(double u, double v,
                                double& nx, double& ny, double& nz)
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::contained_in_surface( double u, double v )
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
int CMLSurfEvalImpl::num_poles()
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::is_planar()
{

    math::Vector n0 = m_normal[0];
    n0.normalize();
    double epsilon = 1e-5;
    
    for(auto i=1; i<m_tri.size();i++){
        math::Vector ni = m_normal[i];
        ni.normalize();
        double vi=std::abs(ni.dot(n0));
        if(vi>1+epsilon || vi<1-epsilon)
            return false;
    }
    
    return true;
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::is_parametric()
{
    return false;
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::is_periodic_in_u(double& u_period)
{
    return false;
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::is_periodic_in_v(double& v_period)
{
    return false;
}
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::get_param_range_u(double& u_low, double& u_high)
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::get_param_range_v(double& v_low, double& v_high)
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::get_geometry_sense()
{
    //by construction it must be true
    return true;
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::uv_from_position(double x, double y, double z,
                                       double& u, double& v)
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool CMLSurfEvalImpl::uv_from_position(double x, double y, double z,
                                       double& u, double& v,
                                       double& cx, double& cy, double& cz)
{
    throw GMDSException("Not yet implemented");
}

/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::position_from_uv(double u, double v,
                                       double& x, double& y, double& z)
{
    throw GMDSException("Not yet implemented");
}

/*---------------------------------------------------------------------------*/
void CMLSurfEvalImpl::distortion_at_uv(double u, double v,
                                       double du[3], double dv[3])
{
    throw GMDSException("Not yet implemented");
}
/*---------------------------------------------------------------------------*/


