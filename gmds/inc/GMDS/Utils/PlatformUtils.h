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
#ifndef GMDS_PLATFORM_UTILS_H_
#define GMDS_PLATFORM_UTILS_H_
/*----------------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/*!
 * \brief Variable d'environnement du nom \a name.
 *
 * Si aucune variable de nom \a name n'est définie,
 * la chaîne nulle est retournée.
 */
/*----------------------------------------------------------------------------*/
class platform{

public:
	static double getMemoryUsed(){
		double mem_used = 0.;

		FILE* f = ::fopen("/proc/self/stat","r");
		if (f){
			// See proc(1), category 'stat'
			int z_pid;
			char z_comm[4096];
			char z_state;
			int z_ppid;
			int z_pgrp;
			int z_session;
			int z_tty_nr;
			int z_tpgid;
			unsigned long z_flags;
			unsigned long z_minflt;
			unsigned long z_cminflt;
			unsigned long z_majflt;
			unsigned long z_cmajflt;
			unsigned long z_utime;
			unsigned long z_stime;
			long z_cutime;
			long z_cstime;
			long z_priority;
			long z_nice;
			long z_zero;
			long z_itrealvalue;
			long z_starttime;
			unsigned long z_vsize;
			long z_rss;
			fscanf(f,"%d %s %c %d %d %d %d %d",
					&z_pid,z_comm,&z_state,&z_ppid,&z_pgrp,&z_session,&z_tty_nr,&z_tpgid);

			fscanf(f,"%lu %lu %lu %lu %lu %lu %lu",
					&z_flags,&z_minflt,&z_cminflt,&z_majflt,&z_cmajflt,&z_utime,&z_stime);

			fscanf(f,"%ld %ld %ld %ld %ld %ld %ld %lu %ld",
					&z_cutime,&z_cstime,&z_priority,&z_nice,&z_zero,&z_itrealvalue,&z_starttime,&z_vsize,&z_rss);

//			mem_used = (double)(z_rss) * (double)getpagesize();
			::fclose(f);
		}

		return mem_used;
	}

	/*!
	 * \brief Mémoire utilisée en octets
	 *
	 * \return la mémoire utilisée ou un nombre négatif si inconnu
	 */

	static std::string getEnvironmentVariable(const std::string& name)
	{

		char* s = ::getenv(name.c_str());
		if (!s)
			return std::string();

		return std::string(s);
	}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_PLATFORM_UTILS_H_ */
/*----------------------------------------------------------------------------*/
