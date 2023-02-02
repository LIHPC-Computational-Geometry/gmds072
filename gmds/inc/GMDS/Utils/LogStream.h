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
/** \file    LogStream.h
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LOG_STREAM_H_
#define GMDS_LOG_STREAM_H_
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
namespace gmds {
/*---------------------------------------------------------------------------*/
typedef enum {
    LOG_ERROR,  // 0
    LOG_WARNING,// 1
    LOG_INFO,   // 2
    LOG_DEBUG,  // 3
    LOG_DEBUG1, // 4
    LOG_DEBUG2  // 5
} LogLevel;
/*---------------------------------------------------------------------------*/
    class Log;
    /*---------------------------------------------------------------------------*/
    /** \class LogStream
     *
     *  \brief Class defining an output to be logged
     */
    /*----------------------------------------------------------------------------*/
    class LogStream{
        
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param AStream the stream the LogStream is connected to
         *  \param ALevel  the level of this output
         */
        LogStream(std::ostream* AStream, const LogLevel& ALevel =LOG_INFO)
        :m_stream(AStream),m_level(ALevel){;}
        
        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor.
         *
         *  \param ALog the LogStream to copy
         */
        LogStream(const LogStream& ALog) {
            m_stream = ALog.m_stream;
            m_level=ALog.m_level;
        }
        
        
        ~LogStream(){;}
        
        /*------------------------------------------------------------------------*/
        /** \brief Overloading of =
         *
         *  \param ALog the LogStream to copy
         */
        const LogStream& operator=(const LogStream& ALog){
            m_stream = ALog.m_stream;
            m_level  = ALog.m_level;
            return *this;
        }
        
        /*------------------------------------------------------------------------*/
        /** \brief give access to the stream
         *
         *  \return the stream
         */
        inline std::ostream& stream() const {return *m_stream;}
        
        /*------------------------------------------------------------------------*/
        /** \brief give consultation access to the level
         *
         *  \return the stream
         */
        inline LogLevel   level() const {return m_level;}
        /*------------------------------------------------------------------------*/
        /** \brief give modification access to the level
         *
         *  \return the stream
         */
        inline LogLevel& level() {return m_level;}
        
        /*----------------------------------------------------------------------------*/
        template<class T> friend inline const LogStream&
        operator<<(const LogStream& ALO,  const T& AT)
        { ALO.stream() << AT; return ALO;}
        
        void flush(){m_stream->flush();}
    protected:
        
        std::ostream* m_stream;
        LogLevel      m_level;
        
    };
    
    /*---------------------------------------------------------------------------*/
    /** \class FileLogStream
     *
     *  \brief Class defining an file stream to be logged
     */
    /*----------------------------------------------------------------------------*/
    class FileLogStream: public LogStream{
        
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param AFileName  the name of the file
         *  \param ALevel  the level of this output
         */
        FileLogStream(const std::string& AFileName,
                      const LogLevel& ALevel =LOG_INFO)
        :LogStream(NULL,ALevel) {
            m_file_name = AFileName;
            m_file_stream = new std::ofstream(m_file_name.c_str(),
                                              std::ios::out | std::ios::app);
            m_stream = m_file_stream;
        }
        
        
        
        ~FileLogStream() {
            m_file_stream->close();
            if(m_file_stream!=NULL)
                delete  m_file_stream;
        }
        
    private:
        std::ofstream* m_file_stream;
        std::string m_file_name;
        
        
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_LOG_STREAM_H_ */
/*----------------------------------------------------------------------------*/
