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
/** \file    Log.h
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LOG_H_
#define GMDS_LOG_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
/*---------------------------------------------------------------------------*/
#include <GMDS/Utils/LogStream.h>
/*---------------------------------------------------------------------------*/
namespace gmds {
    /*-----------------------------------------------------------------------*/
    /** \class Log
     *
     *  \brief Class defining a log manager
     */
    /*-----------------------------------------------------------------------*/
    class Log{
        
    public:
        
        static Log& mng();
        static LogLevel& reportingLevel();
        
        virtual ~Log();
        void addStream(LogStream& AStream);
        
        
        /*------------------------------------------------------------------*/
#ifdef DEBUG_GMDS
        /*------------------------------------------------------------------*/
        template<class T> friend inline const Log&
        operator<<(const Log& o,const T& AT)
        {
            for(unsigned int i=0; i<o.out_.size();i++)
            {
                if(o.out_[i].level()<=Log::reportingLevel())
                    o.out_[i] << AT;
            }
            return o;
        }
        
        void flush(){
            for(unsigned int i=0; i<out_.size();i++)
            {
                out_[i].flush();
            }
        }
        /*------------------------------------------------------------------*/
#else
        /*------------------------------------------------------------------*/
        template<class T> friend inline const Log&
        operator<<(const Log& o,const T& AT)
        {
            for(unsigned int i=0; i<o.out_.size();i++)
            {
                if(o.out_[i].level()<=Log::reportingLevel())
                    o.out_[i] << AT;
            }
            return o;

//            //do nothing
//            return o;
        }
        
        void flush(){
            for(unsigned int i=0; i<out_.size();i++)
            {
                out_[i].flush();
            }
        }
        /*------------------------------------------------------------------*/
#endif //DEBUG_GMDS
        /*------------------------------------------------------------------*/
        
    protected:
        
        Log();
        Log (const Log& ALog);
        Log& operator = (const Log&);
        
    private:
        
        static Log* mng_;
        std::vector<LogStream> out_;
        static LogLevel reporting_level_;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_LOG_H_ */
/*----------------------------------------------------------------------------*/
