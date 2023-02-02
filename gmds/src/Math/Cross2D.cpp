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
/*
 * cross.cpp
 *
 *  Created on: Sep 05, 2014
 *      Author: franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/Math/Cross2D.h"
#include "GMDS/Math/Quaternion.h"
#include "GMDS/Math/Vector.h"
#include "GMDS/Math/Constants.h"
#include "GMDS/Math/Numerics.h"
/*-----------------------------------------------------------------*/
#include <iostream>
#include <math.h>
namespace gmds {
  /*-----------------------------------------------------------------*/
  namespace math{
    /*-----------------------------------------------------------------*/
    Cross2D::Cross2D() : m_angle(0),m_known_vectors(false){
        // m_vectors is then emppty
    }
      /*-----------------------------------------------------------------*/
      Cross2D::Cross2D(Vector& AV1, Vector& AV2)
      {
          if (AV1.dot(AV2) != 0.0)
              throw GMDSException("A CROSS2D object can only be built from 2 orthogonal vectors");
          
          Vector v(1,0,0);
          TCoord a = v.orientedAngle(AV1);
          if(a<0){
              a = Constants::PI2+a;
          }
          
          m_angle = modulo2PI(4*a);
          m_known_vectors = false;
      }
      
      /*-----------------------------------------------------------------*/
      Cross2D::
      Cross2D(const Vector& ARef)
      {
          Vector v(1,0,0);
          TCoord a = v.orientedAngle(ARef);
          if(a<0){
              a = Constants::PI2+a;
          }
          
          m_angle = a;
          m_known_vectors = false;
      }
      /*-----------------------------------------------------------------*/
      Cross2D::
      Cross2D(const TCoord& ARefAngle)
      {
          
      m_angle = modulo2PI(ARefAngle);
    
      m_known_vectors = false;
    }
    /*-----------------------------------------------------------------*/
    Cross2D::Cross2D(const Cross2D& AC)
    {
      m_angle = AC.m_angle;
      m_known_vectors = AC.m_known_vectors;
      if(m_known_vectors)
	m_vectors = AC.m_vectors;
    }
    /*-----------------------------------------------------------------*/
    Vector Cross2D::closestComponentVector(const Vector& AN) const
    {
      std::vector<Vector> vecs = (m_known_vectors)?m_vectors:componentVectors();

      Vector n = AN;
      n.normalize();
      
      Vector v = vecs[0];
      TCoord val = n.dot(v);

      for(unsigned int i=1;i<4;i++){
	Vector v_i = vecs[i];
	TCoord val_i=n.dot(v_i);
	if(val_i> val){
	  val = val_i;
	  v = v_i;
	}
      }

      return v;
    } 
    /*-----------------------------------------------------------------*/
    bool Cross2D::storeComponentVectors() const{
      return m_known_vectors;
    }
      /*-----------------------------------------------------------------*/
      std::vector<Vector> Cross2D::componentVectors() const
      {
          std::vector<Vector> v;
          if(m_known_vectors)
              v = m_vectors;
          else{
              v.resize(4);
              TCoord a = m_angle/4.0;
              for(unsigned int i=0;i<4;i++){
                  TCoord a_i = a +i*Constants::PIDIV2;
                  v[i] = Vector(cos(a_i),sin(a_i),0.0);
              }
          }
          return v;
      }
      /*-----------------------------------------------------------------*/
    void Cross2D::computeComponentVectors() 
    {
      if(!m_known_vectors){
	m_known_vectors = true;
	m_vectors.resize(4);
	TCoord a = m_angle/4.0;
	for(unsigned int i=0;i<4;i++){
	  TCoord a_i = a +i*Constants::PIDIV2;
	  m_vectors[i] = Vector(cos(a_i),sin(a_i),0.0);
	}
      }
    }
    /*-----------------------------------------------------------------*/
    Vector Cross2D::referenceVector() const
    {
      TCoord a = m_angle;
      return Vector(cos(a),sin(a),0.0);
    }
    /*-----------------------------------------------------------------*/
    Cross2D  Cross2D::mean(const vector<Cross2D> & ACrosses, 
			   const vector<TCoord> & AWeights,
			   const TInt ANbSteps)
    {   
      if (ACrosses.empty()){
	throw GMDSException("The mean of zero 2D crosses is undefined.");
      }

      if (AWeights.size() != ACrosses.size())	{
	throw GMDSException("There must exactly the same number of crosses and weights.");
      }

      TCoord ref_angle =  ACrosses[0].referenceAngle(); 
      TCoord prev_ref_angle =  0;
      Vector ref_vector =  ACrosses[0].referenceVector();
      TInt step_index=0;
      do{
	step_index++;
	prev_ref_angle =  ref_angle;
	TCoord pen_angle=0;
	for (unsigned int i = 0; i < ACrosses.size(); i++)
	  {
	    Vector vec_i = ACrosses[i].referenceVector();
	    const TCoord pen_i = ref_vector.orientedAngle(vec_i);
	    pen_angle += pen_i;//*AWeights[i];
	  }
	pen_angle /=ACrosses.size();
	
	ref_angle  = modulo2PI(ref_angle+pen_angle);
	ref_vector = Cross2D(ref_angle).referenceVector();
      }  
      while( (step_index<ANbSteps) && (ref_angle!=prev_ref_angle));

      TCoord new_angle = ref_angle;
      return Cross2D(new_angle);
      
    }
    /*-----------------------------------------------------------------*/
    Cross2D  Cross2D::mean(const Cross2D& AC1, const TCoord& AW1,
			   const Cross2D& AC2, const TCoord& AW2)
    {   
      Vector v1 = AC1.componentVectors()[0];      
      Vector v2 = AC2.closestComponentVector(v1);
      Vector v = AW1*v1 + AW2*v2;
      Vector v_ortho = v.cross(math::Vector(0,0,1));

      return Cross2D(v,v_ortho);
    }

    /*-----------------------------------------------------------------*/
    TCoord Cross2D::angle(const Cross2D& AC) const
    {
      Vector v1 = referenceVector();
      Vector v2 = AC.referenceVector();

      return v1.angle(v2);
    }
    /*-----------------------------------------------------------------*/
    int Cross2D::index(const Cross2D& AC1,
		       const Cross2D& AC2,
		       const Cross2D& AC3)
    {
      Vector vi = AC1.referenceVector();
      Vector vj = AC2.referenceVector();
      Vector vk = AC3.referenceVector();

      
      double wij = vi.orientedAngle(vj);
      double wjk = vj.orientedAngle(vk);
      double wki = vk.orientedAngle(vi);
      
      return (wij+wjk+wki)/Constants::PI2;
    
    }
    /*-----------------------------------------------------------------*/
    ostream & operator << (ostream & AStr, const Cross2D & AC)
    {
      AStr << "Cross2D (" << AC.referenceAngle() << ")";
      return AStr;
    }
    /*-----------------------------------------------------------------*/
    Cross2D operator+(const Cross2D& AC1, const Cross2D& AC2){
      TCoord a1 =  AC1.referenceAngle();
      Vector v1 =  AC1.referenceVector();
      Vector v2 =  AC2.referenceVector();
      TCoord penalty =v1.orientedAngle(v2);
      TCoord new_angle = modulo2PI(a1 + penalty/2.0);      
      return Cross2D(new_angle);
    }
  }
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
