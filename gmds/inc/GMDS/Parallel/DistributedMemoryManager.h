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
/** \file    DistributedMemoryManager.h
 *  \author  F. LEDOUX
 *  \date    03/05/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_DISTRIBUTEDMEMORYMANAGER_H_
#define GMDS_DISTRIBUTEDMEMORYMANAGER_H_
/*----------------------------------------------------------------------------*/
// STL file headers
#include <string>
/*----------------------------------------------------------------------------*/
// GMDS file headers
#include "GMDS/Utils/CommonTypes.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** \class DistributedMemoryManager
 *
 *  \brief It is a singleton. It gives some general services for parallel
 *  	   implementation. The current implementation use mpi as underlying
 *  	   framework.
 */
/*----------------------------------------------------------------------------*/
class DistributedMemoryManager {

public:


	/*------------------------------------------------------------------------*/
    /** \brief  Destructor
     */
	~DistributedMemoryManager();

	/** Message severity level*/
	typedef enum {DEBUG1,DEBUG2,INFO,WARNING,ERROR} MessageLevel;

	/*------------------------------------------------------------------------*/
    /** \brief  Provides the unique instance of the manager
     */
	static DistributedMemoryManager* instance();

	/*------------------------------------------------------------------------*/
    /** \brief  Initialisation step necessary for mpi
     */
	void init(int &argc, char **&argv);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a barrier.
     */
	void barrier(int, const char*);
	/*------------------------------------------------------------------------*/
    /** \brief  Finalize
     */
	void finalize();
	/*------------------------------------------------------------------------*/
	/** \brief  Computes wall time
	 */
	double wTime () const;

	/*------------------------------------------------------------------------*/
	/** \brief gets the processor name
	 *
	 * @return the processor name
	 */
	std::string processorName() const;

	/*------------------------------------------------------------------------*/
	/** \brief Set the verbosity level of the messages (0,1,2,3)
	 */
	inline void setVertbosityLevel(int i){debug_verbose_level_ = i;}

	/// prints a message, same format as printf
	void message(MessageLevel lev, const char *fmt);

	/// abort a calculation
	void abort();

	static void setBufferSize (int b) {buffer_size_ = b;}

#ifdef GMDS_PARALLEL

	inline int rank() { return rank_; }

	inline int nbParts() { return nb_partitions_; }

	inline bool isMaster() { return rank_==0; }

	/// Envoie d'un message
	void send(int toProc, int tag,  void* tabTMP, int & tailleTMP);

	/// Réception d'un message
	void* receive(int fromProc, int tag, int & tailleTMP);


	/// collecte
    void allGather(const std::vector<int>& ASend, std::vector<int>& AGather);

#else
  /// gets the processor id
  inline int rank() { return 0; }

  /// gets the number of processors
  inline int nbParts() { return 1; }

  /// tells if it's processor 0
  inline bool  isMaster() { return 1; }
#endif

  inline void   resetBarrierSensor() { timeSpentOnBarriers_ = 0; }

  inline double getBarrierSensor()   { return timeSpentOnBarriers_; }
private:

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor
     */
	DistributedMemoryManager();

private:

	/** pointer on the unique instance of this manager */
	static DistributedMemoryManager *instance_;

	/** indicate the level of the verbosity of debug messages */
	int debug_verbose_level_;

	/** processor name */
	std::string procName_;

	double timeSpentOnBarriers_;
#ifdef GMDS_PARALLEL
	int rank_;
	int nb_partitions_;
#endif
	static TInt buffer_size_;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
//#include "GMDSPar/DistributedMemoryManager.t.h"
/*----------------------------------------------------------------------------*/
#endif /* GMDS_DISTRIBUTEDMEMORYMANAGER_H_*/
/*----------------------------------------------------------------------------*/

