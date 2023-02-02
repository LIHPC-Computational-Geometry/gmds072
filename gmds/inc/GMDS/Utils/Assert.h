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
/** \file    Assert.h
 *  \author  F. LEDOUX
 *  \date    04/11/2016
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ASSERT_H_
#define GMDS_ASSERT_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonFlags.h>
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    /**
     * \enum  Assert termination mode
     * \brief Defines how assertions must be terminated, either with launching
     *        an exception or abording the program.
     */
    enum AssertionMode {
        ASSERT_MODE_THROW,
        ASSERT_MODE_ABORT
    };
    /*------------------------------------------------------------------------*/
    /** \brief Sets the assertion mode to \p AMode.
     *
     *  \param[in] AMode the assertion mode
     */
    void EXPORT_GMDS setAssertMode(AssertionMode AMode);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the current assertion mode
     */
    AssertionMode EXPORT_GMDS getAssertMode();

    /*------------------------------------------------------------------------*/
    /**
     * \brief Prints an assertion failure when \p ACondition is evaluated to 
     *        false. The termination is done according to the current assertion
     *        mode.
     *
     * \param[in] ACondition String representation of the condition
     * \param[in] AFileName  Name of the file where the assertion is raised
     * \param[in] ALineNum   Number Line where the assertion is raised
     */
    void EXPORT_GMDS assertFailed(const std::string& ACondition,
                                  const std::string& AFileName,
                                  int ALineNum);
    
    /*------------------------------------------------------------------------*/
    /** \brief Prints a range assertion failure when \p AVal is out of the range
     *         [\p AMin, \p AMax\]. The termination is done according to the
     *         current assertion mode.
     *
     * \param[in] AVal      the value to be checked
     * \param[in] AMin      Minimum value
     * \param[in] AMax      Maximum value
     * \param[in] AFileName Name of the file where the assertion is raised
     * \param[in] ALineNum  Number Line where the assertion is raised     
     */
    void EXPORT_GMDS assertRangeFailed(const double& AVal,
                                       const double& AMin,
                                       const double& AMax,
                                       const std::string& AFileName,
                                       int ALineNum);
}
/*----------------------------------------------------------------------------*/
// Two levels of assert,
// use gmds_XXX_assert()        for assertions called in release and debug mode
// use gmds_XXX__debug_assert() for assertions called only in debug mode
/*----------------------------------------------------------------------------*/
/** \brief Check that condition \p AX is verified (so evaluated to true)
 * 
 * \param[in] AX the boolean expression of the condition
 *
 * \see gmds::assert_failed()
 */
#define GMDS_ASSERT(AX) {                                \
        if(!(AX)) {                                      \
            gmds::assertFailed(#AX, __FILE__, __LINE__);\
        }                                                \
}
/*----------------------------------------------------------------------------*/
/** \brief Check that that \p AV is in [\p AMin, \p AMax]
 *
 * \param[in] AX   the value to be checked
 * \param[in] AMin minimum authorized value
 * \param[in] AMax maximum authorized value
 *
 * \see gmds::assert_range_failed()
 */
#define GMDS_RANGE_ASSERT(AX, AMin, AMax) {        \
        if(((AX) < (AMin)) || ((AX) > (AMax))) {   \
            gmds::assertRangeFailed(AX, AMin, AMax \
                __FILE__, __LINE__                 \
            );                                     \
        }                                          \
}
/*----------------------------------------------------------------------------*/
#ifdef GEO_DEBUG
#define geo_debug_assert(x) geo_assert(x)
#define geo_debug_range_assert(x, min_val, max_val) geo_range_assert(x, min_val, max_val)
#else
#define geo_debug_assert(x)
#define geo_debug_range_assert(x, min_val, max_val)
#endif
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_ASSERT_H_ */
/*----------------------------------------------------------------------------*/


