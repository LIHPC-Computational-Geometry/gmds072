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
#ifndef SH_CAMAL_SURF_EVAL_IMPL_H_
#define SH_CAMAL_SURF_EVAL_IMPL_H_
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// CAMAL File Headers
#include <CMLSurfEval.hpp>
/*----------------------------------------------------------------------------*/
/** \class CMLSurfEvalImpl
 *  \brief Implementation of the surface representation for calling CAMAL
 *         paver
 */
class EXPORT_GMDS CMLSurfEvalImpl: public CMLSurfEval{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor from a vector of triangular faces that define the
     *         surface
     *
     */
    CMLSurfEvalImpl(std::vector<gmds::Face>& AFaces,
                    std::vector<gmds::math::Vector>& ANormals);

    ~CMLSurfEvalImpl();
    
    //! \brief Calculate the surface area.
    //!
    //! \return The surface area of the region.
    virtual double area();
    
    //! \brief Calculate the bounding box for a surface
    //!
    //! \param min_xyz The minimum coordinate of the bounding box.
    //! \param max_xyz The maximum coordinate of the bounding box.
    virtual void bounding_box(double min_xyz[3], double max_xyz[3]);
    
    //! \brief Move a point near the surface to the closest point
    //! on the surface.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    virtual void move_to_surface(double& x, double& y, double& z );
    virtual void move_to(double& x, double& y, double& z );
//// FRANCK, since we gave you CAMAL, we renamed move_to_surface to move_to.  
//// Our CMLPave calls move_to, and your version of CMLPave calls move_to_surface,
//// So this change defines both so it will work for both of us.
//// -Matt Staten

    //! \brief Move a point near the surface to the closest point
    //! on the surface, along a direction.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param along_x The x coordinate of the direction
    //! \param along_y The y coordinate of the direction
    //! \param along_z The z coordinate of the direction
    virtual bool  move_to_surface_along_vector(double& x, double& y, double& z,
                                               double &along_x, double &along_y, double &along_z);
    
    
    //! \brief Move a point near the surface to the closest point on the
    //! surface.
    //!
    //! Uses the guesses for the u,v parametric coordinates of the point.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param u_guess A guess for the u parametric value of the point
    //! \param v_guess A guess for the v parametric value of the point
    virtual void move_to_surface(double& x, double& y, double& z,
                                 double& u_guess, double& v_guess);
    
    //! \brief Move a point near the surface to the closest point on the
    //! surface.
    //!
    //! Uses the id of a point on a discrete surface as a starting location
    //! for its search.  This should just call basic move_to_surface if
    //! no vert_idx is available or for surfaces that are not discrete
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param vert_idx A guess for the id of a point on the discrete surface
    virtual void move_to_surface(double& x, double& y, double& z,
                                 int &vert_idx );
    
    //! \brief Get the surface normal at the closest point to x, y, z.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param nx The x component of the normal at the point
    //! \param ny The y component of the normal at the point
    //! \param nz The z component of the normal at the point
    //!
    //! \return \a true if normal has unit length (normalized),
    //! \a false otherwise
    virtual bool normal_at(double x, double y, double z,
                           double& nx, double& ny, double& nz);
    
    //! \brief Get the surface normal at the closest point to x, y, z.
    //!
    //! Uses the guesses for the u,v parametric coordinates of the point.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param u_guess A guess for the u parametric value of the point
    //! \param v_guess A guess for the v parametric value of the point
    //! \param nx The x component of the normal at the point
    //! \param ny The y component of the normal at the point
    //! \param nz The z component of the normal at the point
    //!
    //! \return \a true if normal has unit length (normalized),
    //! \a false otherwise
    virtual bool normal_at(double x, double y, double z,
                           double& u_guess, double& v_guess,
                           double& nx, double& ny, double& nz);
    
    
    //! \brief Get the surface normal at the closest point to u, v
    //!
    //! \param u The u parametric coordinate of the point
    //! \param v The v parametric coordinate of the point
    //! \param nx The x component of the normal at the point
    //! \param ny The y component of the normal at the point
    //! \param nz The z component of the normal at the point
    //!
    //! \return \a true if normal has unit length (normalized),
    //! \a false otherwise
    virtual bool normal_at(double u, double v,
                           double& nx, double& ny, double& nz);
    
    //! \brief for parametric surfaces, determine if a uv point
    //!        is inside of a surface.
    //!
    //! \param u The u parametric coordinate of the point
    //! \param v The v parametric coordinate of the point
    //!
    //! \return \a true if point is inside surface,
    //! \a false otherwise
    virtual bool contained_in_surface( double u, double v );
    
    //! \brief for parametric surfaces, count the number of singular
    //!                  points on this surface.
    //!
    //! \return \a int number of poles
    virtual int num_poles();
    
    //! \brief Get the planar status of the surface.
    //!
    //! \return \a true if surface is planar, \a false otherwise.
    virtual bool is_planar();
    
    //! \brief Get the planar status of the surface.
    //!
    //! \return \a true if surface is parametric, \a false otherwise.
    virtual bool is_parametric();
    
    //! \brief Get the periodic surface status in the u parameter direction.
    //!
    //! \param u_period The period if the surface is periodic in u.
    //!
    //! \return \a true if the surface is periodic in u.
    virtual bool is_periodic_in_u(double& u_period);
    
    //! \brief Get the periodic surface status in the v parameter direction.
    //!
    //! \param v_period The period if the surface is periodic in v.
    //!
    //! \return \a true if the surface is periodic in v.
    virtual bool is_periodic_in_v(double& v_period);
    
    //! \brief Get the periodic range for a periodic surface in u.
    //!
    //! \param u_low The minimum u parameter if periodic.
    //! \param u_high The maximum u parameter if periodic.
    virtual void get_param_range_u(double& u_low, double& u_high);
    
    //! \brief Get the periodic range for a periodic surface in v.
    //!
    //! \param v_low The minimum v parameter if periodic.
    //! \param v_high The maximum v parameter if periodic.
    virtual void get_param_range_v(double& v_low, double& v_high);
    
    //! \brief Get the geometry sense
    //!
    //! return true if the uv handed-ness agrees with the normal of the surface
    virtual bool get_geometry_sense();
    
    //! Get the u,v parametric coordinates on the surface closest to x,y,z.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param u The u parametric coordinate of the point closest to \a x,y,z.
    //! \param v The v parametric coordinate of the point closest to \a x,y,z.
    //!
    //! \return \a true if successful, \a false otherwise
    virtual bool uv_from_position(double x, double y, double z,
                                  double& u, double& v);
    
    //! \brief Get the u,v parametric coordinates on the surface closest
    //! to x,y,z.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param u The u parametric coordinate of the point closest to \a x,y,z.
    //! \param v The v parametric coordinate of the point closest to \a x,y,z.
    //! \param cx The x coordinate of the closest point to \a u,v
    //! on the surface
    //! \param cy The y coordinate of the closest point to \a u,v
    //! on the surface
    //! \param cz The z coordinate of the closest point to \a u,v
    //! on the surface
    //!
    //! \return \a true if successful, \a false otherwise
    virtual bool uv_from_position(double x, double y, double z,
                                  double& u, double& v,
                                  double& cx, double& cy, double& cz);
    
    //! \brief Get the xyz coordinates on the surface closest to u,v
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //! \param u The u parametric coordinate of the point
    //! \param v The v parametric coordinate of the point
    virtual void position_from_uv(double u, double v,
                                  double& x, double& y, double& z);
    
    //! \brief Get the distortion vectors in the parametric surface at u,v
    //!
    //! \param u The u parametric coordinate of the point
    //! \param v The v parametric coordinate of the point
    //! \param du The 3 components of distortion in u
    //! \param dv The 3 components of distortion in v
    virtual void distortion_at_uv(double u, double v, 
                                  double du[3], double dv[3]);
    
private:
    
    std::vector<gmds::Face> m_tri;
    std::vector<gmds::math::Vector> m_normal;
    
};

/*----------------------------------------------------------------------------*/
#endif /* SH_CAMAL_SURF_EVAL_IMPL_H_ */
/*----------------------------------------------------------------------------*/
