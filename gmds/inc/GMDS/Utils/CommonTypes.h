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
 * CommonTypes.h
 *
 *  Created on: 6 janv. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_COMMONTYPES_H_
#define GMDS_COMMONTYPES_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonFlags.h>
#include <GMDS/Utils/Assert.h>
#include <GMDS/Utils/Exception.h>
/*----------------------------------------------------------------------------*/
// STL File headers
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <climits>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    typedef int TInt;
    /*----------------------------------------------------------------------------*/
    /** \type type used to define cell ids locally to a part */
    /*----------------------------------------------------------------------------*/
    typedef unsigned int TCellID;
    
    
    const TCellID NullID = UINT_MAX;
    const TCellID InfinityID = NullID-1;
    const TCellID MaxID = InfinityID-1;
    /*----------------------------------------------------------------------------*/
    /* \type TMark used to mark the different map entities. */
    /*----------------------------------------------------------------------------*/
    typedef std::int32_t TMark;
    /*----------------------------------------------------------------------------*/
    /* \type TCoord used for the coordinate of a point in the space. */
    /*----------------------------------------------------------------------------*/
    typedef double TCoord;
    /*----------------------------------------------------------------------------*/
    const TCoord  TCoord_Epsilon= 1e-4;//1e-11;
    
    /* \const GChunkSize of a chunk used for memory allocation*/
    /*----------------------------------------------------------------------------*/
    const TInt GChunkSize = 128;//4096;// 2048
    /*----------------------------------------------------------------------------*/
    /* \const GNbMaxThreads max number of threads*/
    /*----------------------------------------------------------------------------*/
    const TInt GNbMaxThreads = 16;
    /*----------------------------------------------------------------------------*/
    /** Types of the different usual cells. This type allows us to template element
     * building and editing */
    /*----------------------------------------------------------------------------*/
    typedef enum {
        GMDS_NODE, 				                        /* 0D*/
        GMDS_EDGE,							/* 1D*/
        GMDS_FACE, GMDS_QUAD, GMDS_TRIANGLE, GMDS_POLYGON, 		/* 2D*/
        GMDS_REGION, GMDS_HEX, GMDS_TETRA, GMDS_PRISM3, GMDS_PRISM5,	/* 3D*/
        GMDS_PRISM6, GMDS_PYRAMID, GMDS_POLYHEDRA			/* 3D*/
    } ECellType;
    /*----------------------------------------------------------------------------*/
    typedef enum {
        N2N = 1		  , N2E = 1 << 1  , N2F = 1 << 2  , N2R = 1 << 3  ,
        E2N = 1 << 4  , E2E = 1 << 5  , E2F = 1 << 6  , E2R = 1 << 7  ,
        F2N = 1 << 8  , F2E = 1 << 9  , F2F = 1 << 10 , F2R = 1 << 11 ,
        R2N = 1 << 12 , R2E = 1 << 13 , R2F = 1 << 14 , R2R = 1 << 15 ,
        DIM0 = 1 << 16,	DIM1 = 1 << 17,	DIM2 = 1 << 18,	DIM3 = 1 << 19,
        N = 1 << 20	  ,	E = 1 << 21   ,	F = 1 << 22   ,	R = 1 << 23
    } EMeshDefinition;
    /*----------------------------------------------------------------------------*/
    class EXPORT_GMDS MeshModel{
    public:
        MeshModel(const int meshDef):m_meshDef(meshDef){;}
        
        int getDef() const {return m_meshDef;}
        int getDim() const {
            if(has(DIM2))
                return 2;
            else if (has(DIM1))
                return 1;
            else if (has(DIM3))
                return 3;
            else
                return 0;
        }
        bool operator==(const MeshModel& AModel){return m_meshDef==AModel.m_meshDef;}
        
        bool has(const EMeshDefinition& e) const {return ((m_meshDef|e)==m_meshDef);}
        void add(const int& e)  {m_meshDef = m_meshDef|e;}
        
        /*** Get the model elements common to models A and B */
        static MeshModel intersection(const MeshModel& A, const MeshModel& B){
            int mesh_def=0;
            if(A.has(N) && B.has(N))
                mesh_def = mesh_def | N;
            if(A.has(E) && B.has(E))
                mesh_def = mesh_def | E;
            if(A.has(F) && B.has(F))
                mesh_def = mesh_def | F;
            if(A.has(R) && B.has(R))
                mesh_def = mesh_def | R;
            
            if(A.has(N2N) && B.has(N2N))
                mesh_def = mesh_def | N2N;
            if(A.has(N2E) && B.has(N2E))
                mesh_def = mesh_def | N2E;
            if(A.has(N2F) && B.has(N2F))
                mesh_def = mesh_def | N2F;
            if(A.has(N2R) && B.has(N2R))
                mesh_def = mesh_def | N2R;
            
            if(A.has(E2N) && B.has(E2N))
                mesh_def = mesh_def | E2N;
            if(A.has(E2E) && B.has(E2E))
                mesh_def = mesh_def | E2E;
            if(A.has(E2F) && B.has(E2F))
                mesh_def = mesh_def | E2F;
            if(A.has(E2R) && B.has(E2R))
                mesh_def = mesh_def | E2R;
            
            if(A.has(F2N) && B.has(F2N))
                mesh_def = mesh_def | F2N;
            if(A.has(F2E) && B.has(F2E))
                mesh_def = mesh_def | F2E;
            if(A.has(F2F) && B.has(F2F))
                mesh_def = mesh_def | F2F;
            if(A.has(F2R) && B.has(F2R))
                mesh_def = mesh_def | F2R;
            
            if(A.has(R2N) && B.has(R2N))
                mesh_def = mesh_def | R2N;
            if(A.has(R2E) && B.has(R2E))
                mesh_def = mesh_def | R2E;
            if(A.has(R2F) && B.has(R2F))
                mesh_def = mesh_def | R2F;
            if(A.has(R2R) && B.has(R2R))
                mesh_def = mesh_def | R2R;
            
            if(A.has(DIM0) && B.has(DIM0))
                mesh_def = mesh_def | DIM0;
            if(A.has(DIM1) && B.has(DIM1))
                mesh_def = mesh_def | DIM1;
            if(A.has(DIM2) && B.has(DIM2))
                mesh_def = mesh_def | DIM2;
            if(A.has(DIM3) && B.has(DIM3))
                mesh_def = mesh_def | DIM3;
            
            
            return MeshModel(mesh_def);
        }
        /*** Get the model elements that are in A and not in B */
        static MeshModel exclusion(const MeshModel& A, const MeshModel& B){
            int mesh_def=0;
            if(A.has(N) && !B.has(N))
                mesh_def = mesh_def | N;
            if(A.has(E) && !B.has(E))
                mesh_def = mesh_def | E;
            if(A.has(F) && !B.has(F))
                mesh_def = mesh_def | F;
            if(A.has(R) && !B.has(R))
                mesh_def = mesh_def | R;
            
            if(A.has(N2N) && !B.has(N2N))
                mesh_def = mesh_def | N2N;
            if(A.has(N2E) && !B.has(N2E))
                mesh_def = mesh_def | N2E;
            if(A.has(N2F) && !B.has(N2F))
                mesh_def = mesh_def | N2F;
            if(A.has(N2R) && !B.has(N2R))
                mesh_def = mesh_def | N2R;
            
            if(A.has(E2N) && !B.has(E2N))
                mesh_def = mesh_def | E2N;
            if(A.has(E2E) && !B.has(E2E))
                mesh_def = mesh_def | E2E;
            if(A.has(E2F) && !B.has(E2F))
                mesh_def = mesh_def | E2F;
            if(A.has(E2R) && !B.has(E2R))
                mesh_def = mesh_def | E2R;
            
            if(A.has(F2N) && !B.has(F2N))
                mesh_def = mesh_def | F2N;
            if(A.has(F2E) && !B.has(F2E))
                mesh_def = mesh_def | F2E;
            if(A.has(F2F) && !B.has(F2F))
                mesh_def = mesh_def | F2F;
            if(A.has(F2R) && !B.has(F2R))
                mesh_def = mesh_def | F2R;
            
            if(A.has(R2N) && !B.has(R2N))
                mesh_def = mesh_def | R2N;
            if(A.has(R2E) && !B.has(R2E))
                mesh_def = mesh_def | R2E;
            if(A.has(R2F) && !B.has(R2F))
                mesh_def = mesh_def | R2F;
            if(A.has(R2R) && !B.has(R2R))
                mesh_def = mesh_def | R2R;
            
            if(A.has(DIM0) && !B.has(DIM0))
                mesh_def = mesh_def | DIM0;
            if(A.has(DIM1) && !B.has(DIM1))
                mesh_def = mesh_def | DIM1;
            if(A.has(DIM2) && !B.has(DIM2))
                mesh_def = mesh_def | DIM2;
            if(A.has(DIM3) && !B.has(DIM3))
                mesh_def = mesh_def | DIM3;
            
            return MeshModel(mesh_def);
        }
        
        friend std::ostream& operator<<(std::ostream& stream, const gmds::MeshModel& model);
    private:
        int m_meshDef;
    };
    /*----------------------------------------------------------------------------*/
    typedef enum {
        size_undef = 0,
        size_1 = 1,
        size_2 = 2,
        size_3 = 3,
        size_4 = 4,
        size_5 = 5,
        size_6 = 6,
        size_7 = 7,
        size_8 = 8,
    }TabSize;
    /*----------------------------------------------------------------------------*/
    bool EXPORT_GMDS isZero(const TCoord a);
    bool EXPORT_GMDS areEquals(const TCoord a, const TCoord b);
    /*----------------------------------------------------------------------------*/
    template<int N> struct EXPORT_GMDS TabCellID {
    private:
        TCellID val[N];
        
    public:
        TabCellID(){
            for (int i=0;i<N;i++){
                val[i]=NullID;
            }
        }
        
        inline TInt size() const{
            return N;
        }
        void add(TCellID AID){
            for(unsigned int i=0;i<N;i++)
            {
                if(val[i]==NullID){
                    val[i]=AID;
                    return;
                }
            }
            throw GMDSException("Impossible to add an element in a full container");
        }
        void del(TCellID id){
            for(unsigned int i=0;i<N;i++)
                if(val[i]==id)
                    val[i]=NullID;
        }
        void replace(TCellID id1, TCellID id2){
            for(unsigned int i=0;i<N;i++)
                if(val[i]==id1)
                    val[i]=id2;
        }
        
        void values(std::vector<TCellID>& AVec){
            AVec.clear();
            for(unsigned int i=0;i<N;i++)
                if(val[i]!=NullID)
                    AVec.push_back(val[i]);
        }
        
        void allValues(std::vector<TCellID>& AVec){
            AVec.resize(N);
            for(unsigned int i=0;i<N;i++)
                AVec[i]=val[i];
        }
        void resize(TInt ASize){
            throw GMDSException("Irregular operation for this TabCellID specialization");
        }
        
        inline void operator=(const std::vector<TCellID>& AVec){
            if(AVec.size()!=N)
                throw GMDSException("Irregular operation for this TabCellID specialization");
            
            for(unsigned int i=0;i<N;i++)
                val[i]=AVec[i];
        }
        inline TCellID const& operator[](const TInt& i) const{
            return val[i];
        }
        inline TCellID& operator[](const TInt& i) {
            return val[i];
        }
    };
    
    /*----------------------------------------------------------------------------*/
    template<> struct EXPORT_GMDS TabCellID<size_undef> {
        
    private:
        //TCellID *val;
        std::vector<TCellID> val;
        //	size_t nb_vals;
    public:
        TabCellID(){
        }
        
        TabCellID(const TabCellID& ATab)
        {
            val.resize(ATab.size());
            for(unsigned int i=0;i<val.size();i++){
                val[i]=ATab.val[i];
            }
        }
        ~TabCellID(){
        }
        
        void serialize(std::ostream& stream){
            int s = val.size();
            stream.write((char*)&s,sizeof(int));
            stream.write((char*)&val[0],s*sizeof(TCellID));
        }
        
        void unserialize(std::istream& stream){
            val.clear();
            int s;
            stream.read((char*)&s,sizeof(int));
            resize(s);
            stream.read((char*)&val[0],s*sizeof(TCellID));
        }
        
        inline TInt size() const{
            return val.size();
        }
        
        void add(TCellID id){
            val.push_back(id);
        }
        
        void del(TCellID id){
            std::vector<TCellID> toKeep;
            TInt nbToRemove=0;
            for(unsigned int i=0;i<val.size();i++)
                if(val[i]!=id)
                    toKeep.push_back(val[i]);
                else
                    nbToRemove++;
            
            if(nbToRemove!=0){
                val.resize(toKeep.size());
                for(unsigned int i=0;i<toKeep.size();i++)
                    val[i]=toKeep[i];
            }
            
        }
        
        void replace(TCellID id1, TCellID id2){
            for(unsigned int i=0;i<val.size();i++)
                if(val[i]==id1)
                    val[i]=id2;
        }
        
        void resize(TInt ASize){
            if(ASize<0)
                throw GMDSException("Negative size is forbidden");
            val.resize(ASize);
        }
        inline void operator=(const std::vector<TCellID>& AVec){
            val.clear();
            val.resize(AVec.size());
            for(unsigned int i=0;i<AVec.size();i++)
                val[i]=AVec[i];
        }
        
        void values(std::vector<TCellID>& AVec){
            AVec.clear();
            AVec.reserve(val.size());
            for(unsigned int i=0;i<val.size();i++){
                if(val[i]!=NullID)
                    AVec.push_back(val[i]);
            }
        }
        
        void allValues(std::vector<TCellID>& AVec){
            AVec.resize(val.size());
            for(unsigned int i=0;i<val.size();i++)
                AVec[i]=val[i];
        }
        
        inline TCellID const& operator[](const TInt& i) const{
            return val[i];
        }
        inline TCellID& operator[](const TInt& i) {
            return val[i];
        }
        
    };
    /*----------------------------------------------------------------------------*/
    template<int N> EXPORT_GMDS
    std::ostream & operator << (std::ostream & AStream, const TabCellID<N> & ATab){
        for(int i=0;i<ATab.size();i++)
            AStream<<ATab[i]<<" ";
        return AStream;
    }
    /*----------------------------------------------------------------------------*/
    /*----------------------------------------------------------------------------*/
    /** \brief Function that returns the ids that are both in AS1 and AS2 but ABut
     *
     * @param AS1  A first set of ids
     * @param AS2  A secon set of ids
     * @param ABut An ID to avoid
     * @return the set of ids that are both in AS1 and AS2 but not equal to ABut
     */
    std::vector<TCellID> EXPORT_GMDS getCommonBut(const std::vector<TCellID>& AS1,
                                                  const std::vector<TCellID>& AS2, const TCellID ABut);
    /*----------------------------------------------------------------------------*/
    /** \brief Function that filters the contain of a multiset ASet to keep
     * 		   elements having at least ANb occurrences in ASet
     *
     * @param  ASet the set of elements to be filtered
     * @param  ANb the number of occurrences we target
     * @return A vector containing the elements of ASet having at least ANb
     * 		   occurrences in ASet
     */
    std::vector<TCellID>  EXPORT_GMDS keepFilter(const std::vector<TCellID>& ASet, const TInt ANb);
    /*----------------------------------------------------------------------------*/
    /** \class FakeEdge
     *
     *  \brief This class stores an edge relatively to its node ids.
     *  	   The node order is given starting from the lowest node id, N1.
     *
     */
    class EXPORT_GMDS FakeEdge{
        
    public:
        /*-----------------------------------------------------------------------*/
        /** \class ID
         *
         *  \brief defines an id for a fake edge
         */
        struct EdgeID{
            
            EdgeID(const TCellID& AID1=1, const TCellID& AID2=1)
            { id_tab_[0]=AID1; id_tab_[1]=AID2;}
            
            void set(const TCellID& AID1=1, const TCellID& AID2=1)
            { id_tab_[0]=AID1; id_tab_[1]=AID2;}
            
            bool operator==(const EdgeID& AID) const{
                return (id_tab_[0]==AID.id_tab_[0] && id_tab_[1]==AID.id_tab_[1] );
            }
            
            bool operator<(const EdgeID& AID) const{
                return ((id_tab_[0]<AID.id_tab_[0]) ||
                        (id_tab_[0]==AID.id_tab_[0] &&
                         id_tab_[1]<AID.id_tab_[1]) );
            }
            
            TCellID getID1() const {return id_tab_[0];}
            TCellID getID2() const {return id_tab_[1];}
            
        private:
            TCellID id_tab_[2];
        };
        
        FakeEdge(){
            init(NullID,NullID);
        }
        
        FakeEdge(const TCellID& AID1, const TCellID& AID2){
            init(AID1,AID2);
        }
        
        ~FakeEdge(){;}
        
        std::vector<TCellID> node_ids() const {
            std::vector<TCellID> node_ids;
            node_ids.push_back(id_.getID1());
            node_ids.push_back(id_.getID2());
            return node_ids;
        }
        
        TCellID first() const {return id_.getID1();}
        TCellID second() const {return id_.getID2();}
        void set(const TCellID& AID1, const TCellID& AID2){
            init(AID1,AID2);
        }
        
        EdgeID getID() const {return id_;}
        
        bool operator==(const FakeEdge& AF)
        {	return (id_==AF.id_);}
        
        bool operator!=(const FakeEdge& AF)
        {	return !(id_==AF.id_);}
        
        bool operator<(const FakeEdge& AF) const
        {	return (id_<AF.id_);}
        
    private:
        
        void init(const TCellID& AID1, const TCellID& AID2){
            if(AID1<AID2)
                id_.set(AID1,AID2);
            else
                id_.set(AID2,AID1);
            
        }
        
        
    private:
        
        EdgeID id_;
    };
    /*----------------------------------------------------------------------------*/
    /** \class FakeFace
     *
     *  \brief This class stores the node ids of a face. The node order is given
     *		   starting from the lowest node id, N1,  and traversing the nodes
     *		   in the direction of the N1 incident node having the lowest id.
     *
     */
    class EXPORT_GMDS FakeFace{
        
    public:
        
        
        /*-----------------------------------------------------------------------*/
        /** \class ID
         *
         *  \brief defines an id made of 3 id
         */
        struct FaceID{
            
            FaceID(){;}
            FaceID(const std::vector<TCellID>& AIDs){tab_ = AIDs;}
            
            void set(const std::vector<TCellID>& AIDs){tab_ = AIDs;}
            
            bool operator==(const FaceID& AID) const{
                if (tab_.size()!=AID.tab_.size())
                    return false;
                
                for(unsigned int i=0;i<tab_.size();i++)
                    if(tab_[i]!=AID.tab_[i])
                        return false;
                
                return true;
            }
            
            bool operator<(const FaceID& AID) const{
                if(tab_.size()<AID.tab_.size())	return true;
                if(tab_.size()>AID.tab_.size()) return false;
                /* now same size */
                for(unsigned int i=0;i<tab_.size();i++)
                    if(tab_[i]<AID.tab_[i])		  return true;
                    else if (tab_[i]>AID.tab_[i]) return false;
                /* this and AID are then equals */
                return false;
            }
            
            std::vector<TCellID> getNodeIDs() const {return tab_;}
        private:
            std::vector<TCellID> tab_;
        };
        
        FakeFace(const std::vector<TCellID>& ANodeIDs){
            init(ANodeIDs);
        }
        
        ~FakeFace(){;}
        
        FaceID getID() const {return id_;}
        
        std::vector<TCellID> node_ids() const {
            return id_.getNodeIDs();
        }
        
        
        void set(const std::vector<TCellID>& ANodeIDs){
            init(ANodeIDs);
        }
        
        bool operator==(const FakeFace& AF) const
        {	return (id_==AF.id_);}
        
        bool operator<(const FakeFace& AF) const
        {	return (id_<AF.id_);}
        
    private:
        
        void init_front(const std::vector<TCellID>& ANodeIDs, const TInt& AFirstItem,
                        std::vector<TCellID>& AOutNodeIDs){
            for(unsigned int i = AFirstItem;i<ANodeIDs.size();i++)
                AOutNodeIDs.push_back(ANodeIDs[i]);
            for(int i = 0;i<AFirstItem;i++)
                AOutNodeIDs.push_back(ANodeIDs[i]);
        }
        void init_back(const std::vector<TCellID>& ANodeIDs, const TInt& AFirstItem,
                       std::vector<TCellID>& AOutNodeIDs){
            for(int i = AFirstItem;i>=0;i--)
                AOutNodeIDs.push_back(ANodeIDs[i]);
            for(int i = ANodeIDs.size()-1;i>AFirstItem;i--)
                AOutNodeIDs.push_back(ANodeIDs[i]);
        }
        void init(const std::vector<TCellID>& ANodeIDs){
            /* ANodeIDs gives an ordered description of a face by nodes */
            std::vector<TCellID> node_ids;
            node_ids.reserve(ANodeIDs.size());
            
            /* look for the smallest node id*/
            TCellID lowest_id= NullID;
            TInt     lowest_item = 0;
            TInt nb_ids = ANodeIDs.size();
            for(int i=0;i<nb_ids;i++)
                if (ANodeIDs[i]<lowest_id)
                {
                    lowest_id = ANodeIDs[i];
                    lowest_item = i;
                }
            /* study of the surrounding node*/
            TCellID left_value, right_value;
            if(lowest_item==0){
                left_value = ANodeIDs[nb_ids-1];
                right_value = ANodeIDs[1];
            }
            else if(lowest_item==nb_ids-1){
                left_value = ANodeIDs[nb_ids-2];
                right_value = ANodeIDs[0];
            }
            else{
                left_value = ANodeIDs[lowest_item-1];
                right_value = ANodeIDs[lowest_item+1];
            }
            if(left_value<right_value)
                init_back(ANodeIDs,lowest_item,node_ids);
            else
                init_front(ANodeIDs,lowest_item,node_ids);
            
            id_.set(node_ids);
            
        }
        
        
    private:
        
        FaceID id_;
    };
    /*----------------------------------------------------------------------------*/
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
//std::ostream& operator<<(std::ostream& stream, const gmds::MeshModel& model);
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#endif /* GMDS_COMMONTYPES_H_ */
/*----------------------------------------------------------------------------*/
