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
/** \file    DistributedMemoryManager.cpp
 *  \author  F. LEDOUX
 *  \date    03/05/2009
 */
/*----------------------------------------------------------------------------*/
// STL file headers
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
// pour éviter les conflits entre les types de base du C et ceux utilisés par MPICH
#define MPICH_SKIP_MPICXX
#include "mpi.h"
#endif
/*----------------------------------------------------------------------------*/
// GMDS file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Parallel/DistributedMemoryManager.h>
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
// other file headers
#ifdef GMDS_PARALLEL
//#include "autopack.h"
#else
#ifndef MSC
//#include <sys/time.h>
#endif //MSC
#endif //PARALLEL
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
DistributedMemoryManager* DistributedMemoryManager::instance_ = 0;
/*----------------------------------------------------------------------------*/
int DistributedMemoryManager::buffer_size_ = 1024;
/*----------------------------------------------------------------------------*/
DistributedMemoryManager* DistributedMemoryManager::instance()
{
	if(!instance_)
		instance_= new DistributedMemoryManager();

	return instance_;
}
/*----------------------------------------------------------------------------*/
DistributedMemoryManager::~DistributedMemoryManager()
{
	message(INFO,"Finalizing...\n");
	finalize();
}
/*----------------------------------------------------------------------------*/
DistributedMemoryManager::DistributedMemoryManager()
{}
/*----------------------------------------------------------------------------*/
void DistributedMemoryManager::init(int &argc, char **&argv) {

#ifdef GMDS_PARALLEL
	int namelen;
	char name[1024];

	/* mpi initialization */
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_partitions_);

	MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
	MPI_Get_processor_name(name,&namelen);
	procName_ = name;
	debug_verbose_level_ = 1;

	/* autopack initialization */
	//  AP_init(&argc, &argv);
	//  AP_setparam(100000, 1, 10, 100);


	resetBarrierSensor();
#endif
}
/*----------------------------------------------------------------------------*/
double DistributedMemoryManager::wTime() const
{

#ifdef GMDS_PARALLEL
	return MPI_Wtime();
#else
	//#ifndef MSC
	//  struct timeval tp;
	//  struct timezone tzp;
	//  double timeval;
	//
	//  gettimeofday(&tp,&tzp);
	//
	//  timeval = (double) tp.tv_sec;
	//  timeval = timeval + (double) ((double) .000001 * (double) tp.tv_usec);
	//
	//  return(timeval);
	//#else
	return 0;
	//#endif
#endif
}
/*----------------------------------------------------------------------------*/
std::string DistributedMemoryManager::processorName() const
{
#ifdef GMDS_PARALLEL
	return procName_;
#else
	return std::string("localhost");
#endif
}
/*----------------------------------------------------------------------------*/
void DistributedMemoryManager:: message(DistributedMemoryManager::MessageLevel level, const char *fmt)
{
	char buff[1024];
//	va_list  args;
//	va_start (args, fmt);
//	vsprintf(buff, fmt, args);
//	va_end (args);

	switch(level)
	{
	case DEBUG1:
		if(debug_verbose_level_ > 1)
		{
			/* char logname[256];
	  sprintf(logname,"log-proc%d-%s.dat",myrank,procName);
	  log = fopen (logname,"a");
	  fprintf(log,"%s",buff);
	  fclose(log);*/
		}
		break;
	case DEBUG2:
		if(debug_verbose_level_ > 2)
		{
			/*	  char logname[256];
	  sprintf(logname,"log-proc%d-%s.dat",myrank,procName);
	  log = fopen (logname,"a");
	  fprintf(log,"%s",buff);
	  fclose(log);*/
		}
		break;
	case INFO:
		if(debug_verbose_level_ >= 0 && isMaster())
		{
			//	  fprintf(log,"%s",buff);
			fprintf(stdout,"%s",buff);
			//	  fflush(log);
		}
		if(debug_verbose_level_ > 2)
		{
			//	  fprintf(log,"%s",buff);
			//	  fflush(log);
		}
		break;
	case WARNING:
		fprintf(stdout,"Processor %d WARNING : %s",rank(),buff);
		fflush(stdout);
		break;
	case ERROR:
		fprintf(stdout,"FATAL ERROR : %s",buff);
		fflush(stdout);
		abort();
		break;
	}
}
/*----------------------------------------------------------------------------*/
void DistributedMemoryManager::abort()
{
#ifdef GMDS_PARALLEL
	MPI_Abort(MPI_COMM_WORLD, 1);
#else
	abort();
#endif
}
/*----------------------------------------------------------------------------*/
void DistributedMemoryManager::barrier(int line, const char *fn)
{
#ifdef GMDS_PARALLEL
	double t1 = wTime();
	//  message(DEBUG2,"BARRIER : Line %d in %s\n",line,fn);
	MPI_Barrier(MPI_COMM_WORLD);
	//  message(DEBUG2,"BARRIER PASSED : Line %d in %s\n",line,fn);
	timeSpentOnBarriers_ += (wTime()-t1);
#endif
}
/*----------------------------------------------------------------------------*/
void DistributedMemoryManager::finalize()
{
#ifdef GMDS_PARALLEL
	MPI_Finalize();
#endif
}
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::send (int toProc, int tag, void* tabTMP, int & tailleTMP)
{
#ifdef _DEBUG
	std::cout << rank()<<") DistributedMemoryManager::send("<<toProc<<", "<<tag<<", "<<tailleTMP<<")"<<std::endl;
#endif
	int ierror; /* retour d'erreur sur les fonctions MPI */

	ierror = MPI_Send (tabTMP,
			tailleTMP,
			MPI_BYTE,
			toProc,
			tag,
			MPI_COMM_WORLD);

	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication vers "<<toProc<<", lors de MPI_Send"<<std::endl;
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void* DistributedMemoryManager::receive (int fromProc, int tag, int & tailleTMP)
{
#ifdef _DEBUG
	std::cout << rank()<<") DistributedMemoryManager::receive("<<fromProc<<", "<<tag<<", "<<tailleTMP<<")"<<std::endl;
#endif
	void *TabTMP;
	int ierror;        /* retour d'erreur sur lib MPI */
	MPI_Status status; /* pour les communication MPI  */


	if (tailleTMP <= 0)
	{
		/* pour récupérer l'info sur la taille des tableaux qui vont suivre,
		 on attend qu'un message soit arrivé */
#ifdef _DEBUG
		std::cout << rank()<<") MPI_Probe(...)"<<std::endl;
#endif
		ierror = MPI_Probe (fromProc,
				tag,
				MPI_COMM_WORLD,
				&status);

		if (ierror != MPI_SUCCESS)
			std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<fromProc<<", lors de MPI_Probe"<<std::endl;

		/* récupération de la taille du message | tableau */
#ifdef _DEBUG
		std::cout << rank()<<") MPI_Get_count(...)"<<std::endl;
#endif
		ierror = MPI_Get_count (&status, MPI_BYTE, &tailleTMP);

		if (ierror != MPI_SUCCESS) {
			std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<fromProc<<", lors de MPI_Get_count"<<std::endl;
		}
	} /* end if tailleTMP <= 0 */

	/* allocation pour le tableau tampon */
	if (tailleTMP > 0) {
		TabTMP = new char [(size_t) tailleTMP];
	}
	else
		TabTMP = NULL;

	/* réception du message */
#ifdef _DEBUG
	std::cout << rank()<<") MPI_Recv(...)"<<std::endl;
#endif
	ierror = MPI_Recv (TabTMP,
			tailleTMP,
			MPI_BYTE,//MPI_PACKED,
			fromProc,
			tag,
			MPI_COMM_WORLD,
			&status);

	if (ierror != MPI_SUCCESS) {
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<fromProc<<", lors de MPI_Recv"<<std::endl;
	}

	return (TabTMP);
#ifdef _DEBUG
	std::cout << rank()<<") DistributedMemoryManager::receive("<<fromProc<<", "<<tag<<", "<<tailleTMP<<") terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::send(int toProc, std::vector<std::list<id> > & idPairs)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::send("<<toProc<<", idPairs)"<<std::endl;
#endif

	// Il faut calculer la taille du tableau à allouer
	int taille = sizeof(size_t); // pour le nombre de listes
	std::vector<std::list<id> >::iterator iter = idPairs.begin();
	for(; iter != idPairs.end(); ++iter){
		taille+=sizeof(size_t); // pour le nombre d'éléments dans la liste
		taille+=sizeof(id)*(*iter).size();
	}
#ifdef _DEBUG
	std::cout <<rank()<<") Emission d'un tableau de taille = "<<taille<<std::endl;
#endif

	char *TabTMP = new char [(size_t) taille];

	char* position = TabTMP;
	size_t nb = idPairs.size();
	memcpy((void*)position, (void*)&nb, sizeof(size_t)); position+=sizeof(size_t);

	// remplissage de TabTMP suivant idPairs
	iter = idPairs.begin();
	for(; iter != idPairs.end(); ++iter){
		std::list<id> ll = (*iter);
		nb = ll.size();
		memcpy((void*)position, (void*)&nb, sizeof(size_t)); position+=sizeof(size_t);

		std::list<id>::iterator iter2 = ll.begin();
		for(; iter2 != ll.end(); ++iter2){
			memcpy((void*)position, (void*)&(*iter2), sizeof(id)); position+=sizeof(id);
		}
	}

	send(toProc, 1, TabTMP, taille); // TODO faire une liste des tags utilisés

	delete [] TabTMP;
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::send("<<toProc<<", idPairs) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::sendPacked(const int& AToProc,
										  std::vector<id>&  AIDs,
										  std::vector<TInt>&  ANbs,
										  const int& ATag)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", ids->proc)"<<std::endl;
#endif
	int ierror;
	int vec1_size = AIDs.size();
	int vec2_size = ANbs.size();
	// Il faut calculer la taille du tableau à allouer
	int buf_size= sizeof(int); // pour le nombre d'éléments dans le premier vecteur
	buf_size+= sizeof(int); // pour le nombre d'éléments dans le second vecteur
	if(vec1_size!=0)
		buf_size+=sizeof(id)*vec1_size; // pour les éléments du premier vecteur
	if(vec2_size!=0)
		buf_size+=sizeof(TInt)*vec2_size; // pour les éléments du second vecteur

	// le tableau à allouer à la forme
	// longeur vec1; longeur vec2; vec1[0]; ... ; vec1[N]; vec2[0]; ...; vec2[N]
#ifdef _DEBUG
	std::cout <<rank()<<") Emission d'un tableau de taille = "<<buf_size<<" avec (v1,v2)="<<vec1_size<<", "<<vec2_size<<std::endl;
#endif

	/* allocated buffer to send data */
	char *buffer = new char [buf_size];
	/* position in this buffer during the packing process */
	int position = 0;
	MPI_Pack(&vec1_size, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&vec2_size, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);

	if(vec1_size!=0)
		MPI_Pack(&AIDs[0],vec1_size, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	if(vec2_size!=0)
		MPI_Pack(&ANbs[0],vec2_size, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);

	ierror = MPI_Send(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD);
//	MPI_Request* r;
//	ierror = MPI_Isend(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD,r);
	if (ierror != MPI_SUCCESS)
			std::cout<<"Problème sur "<<rank()<<" de communication vers "<<AToProc<<", lors de MPI_Send"<<std::endl;
	delete [] buffer;

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", ids->proc) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::receivePacked(const int& AFromProc,
											 std::vector<id>&  AIDs,
											 std::vector<TInt>&  ANbs,
											 const int& ATag)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receivePacked("<<AFromProc<<", ids->proc)"<<std::endl;
#endif

	int vec1_size, vec2_size;
	int buf_size;
	MPI_Status status;
	int ierror;        /* retour d'erreur sur lib MPI */

#ifdef _DEBUG
		std::cout << rank()<<") MPI_Probe(...)"<<std::endl;
#endif
	 ierror = MPI_Probe(AFromProc,ATag, MPI_COMM_WORLD, &status);
	 if (ierror != MPI_SUCCESS)
			std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Probe"<<std::endl;

		/* récupération de la taille du message | tableau */
#ifdef _DEBUG
	 std::cout << rank()<<") MPI_Get_count(...)"<<std::endl;
#endif
	 ierror = MPI_Get_count(&status, MPI_BYTE, &buf_size);

	 if (ierror != MPI_SUCCESS)
		 std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Get_count"<<std::endl;


		#ifdef _DEBUG
	std::cout <<" réception d'un tableau de taille = "<<buf_size<<" faite, conversion en 2 vecteurs"<<std::endl;
#endif
	 char* buffer = new char[buf_size];

	 ierror = MPI_Recv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&status);
//	 MPI_Request r;
//	 ierror = MPI_Irecv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&r);

	 if (ierror != MPI_SUCCESS) {
		 std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Recv"<<std::endl;
	 }

	 int position =0;
	 ierror = MPI_Unpack(buffer,buf_size,&position,&vec1_size, 1, MPI_INT,MPI_COMM_WORLD);
	 if (ierror != MPI_SUCCESS) {
		 std::cout<<"Problème sur "<<rank()<<" de depactage en provenance de "<<AFromProc<<", lors de MPI_Unpack"<<std::endl;
	 }

	 ierror = MPI_Unpack(buffer,buf_size,&position,&vec2_size, 1, MPI_INT,MPI_COMM_WORLD);
	 if (ierror != MPI_SUCCESS) {
		 std::cout<<"Problème sur "<<rank()<<" de depactage en provenance de "<<AFromProc<<", lors de MPI_Unpack"<<std::endl;
	 }

	  AIDs.clear();
	  ANbs.clear();
	  AIDs.reserve(vec1_size);
	  ANbs.reserve(vec2_size);
	  AIDs.resize(vec1_size);
	  ANbs.resize(vec2_size);

	  if(vec1_size!=0){
		  ierror = MPI_Unpack(buffer,buf_size,&position,&AIDs[0],vec1_size, MPI_LONG,MPI_COMM_WORLD);
		  if (ierror != MPI_SUCCESS) {
			  std::cout<<"Problème sur "<<rank()<<" de depactage en provenance de "<<AFromProc<<", lors de MPI_Unpack"<<std::endl;
		  }
	  }

	  if(vec2_size!=0){
		  ierror = MPI_Unpack(buffer,buf_size,&position,&ANbs[0],vec2_size, MPI_INT,MPI_COMM_WORLD);
		  if (ierror != MPI_SUCCESS) {
			  std::cout<<"Problème sur "<<rank()<<" de depactage en provenance de "<<AFromProc<<", lors de MPI_Unpack"<<std::endl;
		  }

	  }

	  delete [] buffer;

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receivePacked("<<AFromProc<<", ids->proc) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::sendPacked(const int& AToProc,
										  MigrationStruct& AData,
										  const int& ATag)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", migrationStruct)"<<std::endl;
#endif
	int ierror;
	int nb_regions = AData.region_ids.size();
	int nb_faces = AData.face_ids.size();
	int nb_edges = AData.edge_ids.size();
	int nb_nodes = AData.node_ids.size();
	int nb_partition_sharing_regions=0;
	int nb_partition_sharing_faces=0;
	int nb_partition_sharing_edges=0;
	int nb_partition_sharing_nodes=0;
	int nb_partition_group_regions=0;
	int nb_partition_group_faces=0;
	int nb_partition_group_edges=0;
	int nb_partition_group_nodes=0;
	// Il faut calculer la taille du buffer à allouer
	int buf_size= sizeof(int); 					// pour stocker le nombre de regions
	if(nb_regions!=0)	{
		buf_size+=sizeof(id)*nb_regions;  	// pour stocker les ids des regions
		buf_size+= sizeof(ECellType)*nb_regions; 	// pour stocker le type des regions
		for(unsigned int i=0;i<nb_regions;i++) 	// pour stocker les connectivités
		{
			buf_size+= sizeof(TInt)*AData.region_2N[i].size() + sizeof(int);
			// + 1 pour stocker le nombre de noeuds
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de partage
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.region_sharing.begin();
			it_sharing!=AData.region_sharing.end(); it_sharing++)
		{
			nb_partition_sharing_regions++;
			buf_size+= sizeof(id)*(it_sharing->second).size()+ 2*sizeof(int);
			// + 1 pour stocker la partition partageant les regions
			// + 1 pour stocker le nombre de regions partagées
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de groupes
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.region_groups.begin();
			it_group !=AData.region_groups.end(); it_group++)
		{
			nb_partition_group_regions++;
			buf_size+= sizeof(id)*(it_group->second).size()+ sizeof(char)*((it_group->first).size()+1) + 2*sizeof(int);
			// + 1 pour stocker la taille du nom de la partition partageant les regions
			// + 1 pour stocker le nombre de regions partagées
		}

	}
	buf_size+= sizeof(int); 					// pour stocker le nombre de faces
	if(nb_faces!=0)	{
		buf_size+=sizeof(id)*nb_faces;  	// pour stocker les ids des faces
		buf_size+= sizeof(ECellType)*nb_faces; 	// pour stocker le type des faces
		for(unsigned int i=0;i<nb_faces;i++) 	// pour stocker les connectivités
		{
			buf_size+= sizeof(TInt)*AData.face_2N[i].size() + sizeof(int);
			// + 1 pour stocker le nombre de noeuds
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de partage
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.face_sharing.begin();
			it_sharing!=AData.face_sharing.end(); it_sharing++)
		{
			nb_partition_sharing_faces++;
			buf_size+= sizeof(id)*(it_sharing->second).size()+ 2*sizeof(int);
			// + 1 pour stocker la partition partageant les faces
			// + 1 pour stocker le nombre de faces partagées
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de groupes
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.face_groups.begin();
			it_group !=AData.face_groups.end(); it_group++)
		{
			nb_partition_group_faces++;
			buf_size+= sizeof(id)*(it_group->second).size()+ sizeof(char)*((it_group->first).size()+1) + 2*sizeof(int);
			// + 1 pour stocker la taille du nom de la partition partageant les faces
			// + 1 pour stocker le nombre de faces partagées
		}

	}
	buf_size+= sizeof(int); 					// pour stocker le nombre d'aretes
	if(nb_edges!=0)	{
		buf_size+=sizeof(id)*nb_edges;  	// pour stocker les ids des aretes
		buf_size+=2*sizeof(TInt)*nb_edges;  	// pour stocker les connectivites
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de partage
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.edge_sharing.begin();
			it_sharing!=AData.edge_sharing.end(); it_sharing++)
		{
			nb_partition_sharing_edges++;
			buf_size+= sizeof(id)*(it_sharing->second).size()+ 2*sizeof(int);
			// + 1 pour stocker la partition partageant les aretes
			// + 1 pour stocker le nombre d'aretes partagees
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de groupes
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.edge_groups.begin();
			it_group !=AData.edge_groups.end(); it_group++)
		{
			nb_partition_group_edges++;
			buf_size+= sizeof(id)*(it_group->second).size()+ sizeof(char)*((it_group->first).size()+1) + 2*sizeof(int);
			// + 1 pour stocker la taille du nom de la partition partageant les aretes
			// + 1 pour stocker le nombre de aretes partagées
		}

	}
	buf_size+= sizeof(int); 		// pour stocker le nombre de noeuds
	if(nb_nodes!=0) {
		buf_size+=sizeof(id)*nb_nodes;  	// pour stocker les ids des noeuds
		buf_size+=3*sizeof(TCoord)*nb_nodes;   	// pour stocker les coords des noeuds
		buf_size+= sizeof(int); 				// pour stocker le nombre d'info de partage
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.node_sharing.begin();
			it_sharing!=AData.node_sharing.end(); it_sharing++)
		{
			nb_partition_sharing_nodes++;
			buf_size+= sizeof(id)*(it_sharing->second).size()+ 2*sizeof(int);
			// + 1 pour stocker la partition partageant les noeuds
			// + 1 pour stocker le nombre de noeuds partagés
		}
		buf_size+= sizeof(int); 	// pour stocker le nombre d'info de groupes
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.node_groups.begin();
			it_group !=AData.node_groups.end(); it_group++)
		{
			nb_partition_group_nodes++;
			buf_size+= sizeof(id)*(it_group->second).size()+ sizeof(char)*((it_group->first).size()+1) + 2*sizeof(int);
			// + 1 pour stocker la taille du nom de la partition partageant les noeuds
			// + 1 pour stocker le nombre de noeuds partagés
		}
	}

#ifdef _DEBUG
	std::cout <<rank()<<") Emission d'un tableau de taille = "<<buf_size<<std::endl;
#endif

	/* allocated buffer to send data */
	char *buffer = new char [buf_size];
	/* position in this buffer during the packing process */
	int position = 0;
	MPI_Pack(&nb_regions, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&nb_faces, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&nb_edges, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&nb_nodes, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);

	if(nb_regions!=0){
		MPI_Pack(&AData.region_ids [0],nb_regions, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
		MPI_Pack(&AData.region_type[0],nb_regions, MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);

		for(unsigned int i=0;i<nb_regions;i++)
		{
			TInt loc_nb_nodes 			   = AData.region_2N[i].size();
			std::vector<TInt> loc_node_ids = AData.region_2N[i];
			MPI_Pack(&loc_nb_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&loc_node_ids[0],loc_nb_nodes, MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_sharing_regions,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.region_sharing.begin();it_sharing!=AData.region_sharing.end(); it_sharing++)
		{
			int share_partition = it_sharing->first;
			int nb_shared_regions=(it_sharing->second).size();
			MPI_Pack(&share_partition,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_shared_regions,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> shared_regions= it_sharing->second;
			MPI_Pack(&shared_regions[0],nb_shared_regions, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_group_regions,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.region_groups.begin();it_group!=AData.region_groups.end(); it_group++)
		{
			std::string group_name = it_group->first;
			int group_name_size = group_name.size()+1;
			int nb_group_regions=(it_group->second).size();
			char* group_char = (char*)group_name.c_str();
			MPI_Pack(&group_name_size,1, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(group_char,group_name_size, MPI_CHAR, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_group_regions,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> group_regions= it_group->second;
			MPI_Pack(&group_regions[0],nb_group_regions, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

	}

	if(nb_faces!=0){
		MPI_Pack(&AData.face_ids [0],nb_faces, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
		MPI_Pack(&AData.face_type[0],nb_faces, MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);

		for(unsigned int i=0;i<nb_faces;i++)
		{
			TInt loc_nb_nodes 					= AData.face_2N[i].size();
			std::vector<TInt> loc_node_ids = AData.face_2N[i];
			MPI_Pack(&loc_nb_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&loc_node_ids[0],loc_nb_nodes, MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_sharing_faces,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.face_sharing.begin();it_sharing!=AData.face_sharing.end(); it_sharing++)
		{
			int share_partition = it_sharing->first;
			int nb_shared_faces =(it_sharing->second).size();
			MPI_Pack(&share_partition,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_shared_faces,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> shared_faces = it_sharing->second;
			MPI_Pack(&shared_faces[0],nb_shared_faces, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_group_faces,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.face_groups.begin();it_group!=AData.face_groups.end(); it_group++)
		{
			std::string group_name = it_group->first;
			int group_name_size = group_name.size()+1;
			int nb_group_faces=(it_group->second).size();
			char* group_char = (char*)group_name.c_str();
			MPI_Pack(&group_name_size,1, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(group_char,group_name_size, MPI_CHAR, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_group_faces,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> group_faces= it_group->second;
			MPI_Pack(&group_faces[0],nb_group_faces, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

	}
	if(nb_edges!=0){
		MPI_Pack(&AData.edge_ids[0],nb_edges, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
		MPI_Pack(&AData.edge_2N[0],2*nb_edges, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);

		MPI_Pack(&nb_partition_sharing_edges,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.edge_sharing.begin();it_sharing!=AData.edge_sharing.end(); it_sharing++)
		{
			int share_partition = it_sharing->first;
			int nb_shared_edges =(it_sharing->second).size();
			MPI_Pack(&share_partition,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_shared_edges,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> shared_edges= it_sharing->second;
			MPI_Pack(&shared_edges[0],nb_shared_edges, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_group_edges,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.edge_groups.begin();it_group!=AData.edge_groups.end(); it_group++)
		{
			std::string group_name = it_group->first;
			int group_name_size = group_name.size()+1;
			int nb_group_edges=(it_group->second).size();
			char* group_char = (char*)group_name.c_str();
			MPI_Pack(&group_name_size,1, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(group_char,group_name_size, MPI_CHAR, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_group_edges,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> group_edges= it_group->second;
			MPI_Pack(&group_edges[0],nb_group_edges, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

	}
	if(nb_nodes!=0){
		MPI_Pack(&AData.node_ids   [0],nb_nodes, MPI_LONG  , buffer, buf_size, &position,MPI_COMM_WORLD);
		MPI_Pack(&AData.node_coords[0],3*nb_nodes, MPI_DOUBLE, buffer, buf_size, &position,MPI_COMM_WORLD);

		MPI_Pack(&nb_partition_sharing_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<TInt, std::vector<id> >::iterator it_sharing;
		for(it_sharing = AData.node_sharing.begin();it_sharing!=AData.node_sharing.end(); it_sharing++)
		{
			int share_partition = it_sharing->first;
			int nb_shared_nodes =(it_sharing->second).size();
			MPI_Pack(&share_partition,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_shared_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> shared_nodes = it_sharing->second;
			MPI_Pack(&shared_nodes[0],nb_shared_nodes, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}

		MPI_Pack(&nb_partition_group_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
		std::map<std::string, std::vector<id> >::iterator it_group;
		for(it_group = AData.node_groups.begin();it_group!=AData.node_groups.end(); it_group++)
		{
			std::string group_name = it_group->first;
			int group_name_size = group_name.size()+1;
			int nb_group_nodes=(it_group->second).size();
			char* group_char = (char*)group_name.c_str();
			MPI_Pack(&group_name_size,1, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(group_char,group_name_size, MPI_CHAR, buffer, buf_size, &position,MPI_COMM_WORLD);
			MPI_Pack(&nb_group_nodes,1 , MPI_INT , buffer, buf_size, &position,MPI_COMM_WORLD);
			std::vector<id> group_nodes= it_group->second;
			MPI_Pack(&group_nodes[0],nb_group_nodes, MPI_LONG , buffer, buf_size, &position,MPI_COMM_WORLD);
		}
	}
//	MPI_Request* r;
//	ierror = MPI_Isend(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD,r);
	ierror = MPI_Send(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication vers "<<AToProc<<", lors de MPI_Send"<<std::endl;

	delete [] buffer;

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", ids->proc) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::sendPacked(const int& AToProc,
										  RefinementStruct& AData,
										  const int& ATag)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", refinementStruct)"<<std::endl;
#endif
	int ierror;
	int nb_infos = AData.nb_infos;
	// Il faut calculer la taille du buffer à allouer
	int buf_size= sizeof(int)+ sizeof(id)*4*nb_infos;

#ifdef _DEBUG
	std::cout <<rank()<<") Emission d'un tableau de taille = "<<buf_size<<std::endl;
#endif

	/* allocated buffer to send data */
	char *buffer = new char [buf_size];
	/* position in this buffer during the packing process */
	int position = 0;
	MPI_Pack(&nb_infos, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&AData.init_edge_id[0],nb_infos, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_node_id[0] ,nb_infos, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_edge1_id[0],nb_infos, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_edge2_id[0],nb_infos, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);

//	MPI_Request* r;
//	ierror = MPI_Isend(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD,r);
	ierror = MPI_Send(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication vers "<<AToProc<<", lors de MPI_Send"<<std::endl;

	delete [] buffer;

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", ids->proc) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::receivePacked(const int& AFromProc,
											 MigrationStruct& AData,
											 const int& ATag)
{

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receivePacked("<<AFromProc<<", MigrationStruct)"<<std::endl;
#endif

	int buf_size;
	MPI_Status status;
	int ierror;        /* retour d'erreur sur lib MPI */

#ifdef _DEBUG
		std::cout << rank()<<") MPI_Probe(...)"<<std::endl;
#endif
	ierror = MPI_Probe(AFromProc,ATag, MPI_COMM_WORLD, &status);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Probe"<<std::endl;

		/* récupération de la taille du message | tableau */
#ifdef _DEBUG
	std::cout << rank()<<") MPI_Get_count(...)"<<std::endl;
#endif
	ierror = MPI_Get_count(&status, MPI_BYTE, &buf_size);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Get_count"<<std::endl;


#ifdef _DEBUG
	std::cout <<" réception d'un tableau de taille = "<<buf_size<<std::endl;
#endif
	char* buffer = new char[buf_size];
//	MPI_Request r;
//	ierror = MPI_Irecv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&r);
	ierror = MPI_Recv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&status);
	if (ierror != MPI_SUCCESS) {
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Recv"<<std::endl;
	}
	int position =0;

	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_regions, 1, MPI_INT,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_faces, 1, MPI_INT,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_edges, 1, MPI_INT,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_nodes, 1, MPI_INT,MPI_COMM_WORLD);
	TInt nb_regions, nb_faces, nb_edges, nb_nodes;
	nb_regions = AData.nb_regions;
	nb_faces   = AData.nb_faces;
	nb_edges   = AData.nb_edges;
	nb_nodes   = AData.nb_nodes;
#ifdef _DEBUG
	std::cout<<rank()<<") "<<AData.nb_regions<<" regions received from "<<AFromProc<<std::endl;
	std::cout<<rank()<<") "<<AData.nb_faces<<" faces   received from "<<AFromProc<<std::endl;
	std::cout<<rank()<<") "<<AData.nb_edges<<" edges   received from "<<AFromProc<<std::endl;
	std::cout<<rank()<<") "<<AData.nb_nodes<<" nodes   received from "<<AFromProc<<std::endl;
#endif //_DEBUG
	if(nb_regions!=0){
		AData.region_ids.resize(nb_regions);
		AData.region_type.resize(nb_regions);
		AData.region_2N.resize(nb_regions);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.region_ids) [0], nb_regions, MPI_LONG,MPI_COMM_WORLD);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.region_type)[0], nb_regions, MPI_INT ,MPI_COMM_WORLD);
		for(unsigned int i=0;i<nb_regions;i++)
		{
			TInt loc_nb_nodes=0;
			MPI_Unpack(buffer,buf_size,&position,&loc_nb_nodes, 1, MPI_INT ,MPI_COMM_WORLD);
			(AData.region_2N[i]).resize(loc_nb_nodes);
			MPI_Unpack(buffer,buf_size,&position,&(AData.region_2N[i])[0], loc_nb_nodes, MPI_INT,MPI_COMM_WORLD);
		}

		TInt nb_partition_sharing_regions =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_sharing_regions, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_sharing_regions;i++)
		{
			TInt share_partition = 0;
			TInt nb_shared_regions =0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&share_partition, 1, MPI_INT,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_shared_regions, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> shared_regions;
			shared_regions.resize(nb_shared_regions);
			ierror = MPI_Unpack(buffer,buf_size,&position,&shared_regions[0], nb_shared_regions, MPI_LONG ,MPI_COMM_WORLD);
			AData.region_sharing[share_partition] = shared_regions;
		}

		TInt nb_partition_group_regions =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_group_regions, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_group_regions;i++)
		{
			TInt group_name_size=0;
			TInt nb_group_regions=0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_name_size, 1, MPI_INT,MPI_COMM_WORLD);
			char* group_name = new char[group_name_size];
			ierror = MPI_Unpack(buffer,buf_size,&position,group_name, group_name_size, MPI_CHAR,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_group_regions, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> group_regions;
			group_regions.resize(nb_group_regions);
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_regions[0], nb_group_regions, MPI_LONG ,MPI_COMM_WORLD);
			std::string group_string(group_name);
			AData.region_groups[group_name] = group_regions;

			delete group_name;
		}
	}
	if(nb_faces!=0){
		AData.face_ids.resize(nb_faces);
		AData.face_type.resize(nb_faces);
		AData.face_2N.resize(nb_faces);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.face_ids) [0], nb_faces, MPI_LONG,MPI_COMM_WORLD);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.face_type)[0], nb_faces, MPI_INT ,MPI_COMM_WORLD);
		for(unsigned int i=0;i<nb_faces;i++)
		{
			TInt loc_nb_nodes=0;
			MPI_Unpack(buffer,buf_size,&position,&loc_nb_nodes, 1, MPI_INT ,MPI_COMM_WORLD);
			(AData.face_2N[i]).resize(loc_nb_nodes);
			MPI_Unpack(buffer,buf_size,&position,&(AData.face_2N[i])[0], loc_nb_nodes, MPI_INT,MPI_COMM_WORLD);
		}

		TInt nb_partition_sharing_faces =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_sharing_faces, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_sharing_faces;i++)
		{
			TInt share_partition = 0;
			TInt nb_shared_faces =0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&share_partition, 1, MPI_INT,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_shared_faces, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> shared_faces;
			shared_faces.resize(nb_shared_faces);
			ierror = MPI_Unpack(buffer,buf_size,&position,&shared_faces[0], nb_shared_faces, MPI_LONG ,MPI_COMM_WORLD);
			AData.face_sharing[share_partition] = shared_faces;
		}

		TInt nb_partition_group_faces=0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_group_faces, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_group_faces;i++)
		{
			TInt group_name_size=0;
			TInt nb_group_faces=0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_name_size, 1, MPI_INT,MPI_COMM_WORLD);
			char* group_name = new char[group_name_size];
			ierror = MPI_Unpack(buffer,buf_size,&position,group_name, group_name_size, MPI_CHAR,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_group_faces, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> group_faces;
			group_faces.resize(nb_group_faces);
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_faces[0], nb_group_faces, MPI_LONG ,MPI_COMM_WORLD);
			std::string group_string(group_name);
			AData.face_groups[group_name] = group_faces;

			delete group_name;
		}

	}
	if(nb_edges!=0){
		AData.edge_ids.resize(nb_edges);
		AData.edge_2N.resize(2*nb_edges);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.edge_ids) [0], nb_edges, MPI_LONG,MPI_COMM_WORLD);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.edge_2N)[0], 2*nb_edges, MPI_INT ,MPI_COMM_WORLD);

		TInt nb_partition_sharing_edges =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_sharing_edges, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_sharing_edges;i++)
		{
			TInt share_partition = 0;
			TInt nb_shared_edges =0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&share_partition, 1, MPI_INT,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_shared_edges, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> shared_edges;
			shared_edges.resize(nb_shared_edges);
			ierror = MPI_Unpack(buffer,buf_size,&position,&shared_edges[0], nb_shared_edges, MPI_LONG ,MPI_COMM_WORLD);
			AData.edge_sharing[share_partition] = shared_edges;
		}

		TInt nb_partition_group_edges =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_group_edges, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_group_edges;i++)
		{
			TInt group_name_size=0;
			TInt nb_group_edges=0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_name_size, 1, MPI_INT,MPI_COMM_WORLD);
			char* group_name = new char[group_name_size];
			ierror = MPI_Unpack(buffer,buf_size,&position,group_name, group_name_size, MPI_CHAR,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_group_edges, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> group_edges;
			group_edges.resize(nb_group_edges);
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_edges[0], nb_group_edges, MPI_LONG ,MPI_COMM_WORLD);
			std::string group_string(group_name);
			AData.edge_groups[group_name] = group_edges;

			delete group_name;
		}

	}

	if(nb_nodes!=0){
		AData.node_ids.resize(nb_nodes);
		AData.node_coords.resize(nb_nodes*3);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.node_ids) [0], nb_nodes  , MPI_LONG  ,MPI_COMM_WORLD);
		ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.node_coords)[0], nb_nodes*3, MPI_DOUBLE,MPI_COMM_WORLD);

		TInt nb_partition_sharing_nodes =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_sharing_nodes, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_sharing_nodes;i++)
		{
			TInt share_partition = 0;
			TInt nb_shared_nodes =0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&share_partition, 1, MPI_INT,MPI_COMM_WORLD);
			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_shared_nodes, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> shared_nodes;
			shared_nodes.resize(nb_shared_nodes);
			ierror = MPI_Unpack(buffer,buf_size,&position,&shared_nodes[0], nb_shared_nodes, MPI_LONG ,MPI_COMM_WORLD);
			AData.node_sharing[share_partition] = shared_nodes;
		}

		TInt nb_partition_group_nodes =0;
		ierror = MPI_Unpack(buffer,buf_size,&position,&nb_partition_group_nodes, 1, MPI_INT,MPI_COMM_WORLD);
		for(int i=0;i<nb_partition_group_nodes;i++)
		{

			TInt group_name_size=0;
			TInt nb_group_nodes =0;
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_name_size, 1, MPI_INT,MPI_COMM_WORLD);
			char* group_name = new char[group_name_size];
			ierror = MPI_Unpack(buffer,buf_size,&position,group_name, group_name_size, MPI_CHAR,MPI_COMM_WORLD);

			ierror = MPI_Unpack(buffer,buf_size,&position,&nb_group_nodes, 1, MPI_INT,MPI_COMM_WORLD);
			std::vector<id> group_nodes;
			group_nodes.resize(nb_group_nodes);
			ierror = MPI_Unpack(buffer,buf_size,&position,&group_nodes[0], nb_group_nodes, MPI_LONG ,MPI_COMM_WORLD);
			std::string group_string(group_name);

			AData.node_groups[group_string] = group_nodes;
			delete group_name;
		}
	}

	delete [] buffer;
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::sendPacked(const int& AToProc,
										  RefinementStruct3D& AData,
										  const int& ATag)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", refinementStruct3D)"<<std::endl;
#endif
	int ierror;
	int nb_faces = AData.nb_faces;
	int nb_edges = AData.nb_edges;

	// Il faut calculer la taille du buffer à allouer

	TInt total_nb_sides = AData.new_subface_id.size();

	int buf_size= sizeof(int)	+				 /* nb faces */
				  sizeof(int)	+				 /* nb edges */
				  sizeof(int)	+				 /* total of subfaces*/
				  sizeof(id)*2*nb_faces + /* init face id + associated node id */
				  sizeof(TInt)*nb_faces +        /* number of sides of faces */
				  sizeof(id)*total_nb_sides + /* ids of new subfaces*/
				  sizeof(id)*3*nb_edges; /* 	fake edge to new node ids*/


#ifdef _DEBUG
	std::cout <<rank()<<") Emission d'un tableau de taille = "<<buf_size<<std::endl;
#endif

	/* allocated buffer to send data */
	char *buffer = new char [buf_size];
	/* position in this buffer during the packing process */
	int position = 0;
	MPI_Pack(&nb_faces, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&nb_edges, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&total_nb_sides, 1, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);
	MPI_Pack(&AData.init_face_id[0],nb_faces, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.init_face_nb_sides[0] ,nb_faces, MPI_INT, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_face_node_id[0],nb_faces, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_subface_id[0],total_nb_sides, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);
	MPI_Pack(&AData.new_edge_node_id[0],3*nb_edges, MPI_LONG, buffer, buf_size, &position,MPI_COMM_WORLD);

//	MPI_Request* r;
//	ierror = MPI_Isend(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD,r);
	ierror = MPI_Send(buffer,position,MPI_PACKED,AToProc,ATag,MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication vers "<<AToProc<<", lors de MPI_Send"<<std::endl;

	delete [] buffer;

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::sendPacked("<<AToProc<<", ids->proc) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::receivePacked(const int& AFromProc,
											 RefinementStruct3D& AData,
											 const int& ATag)
{

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receivePacked("<<AFromProc<<", RefinementStruct3D)"<<std::endl;
#endif

	int buf_size;
	MPI_Status status;
	int ierror;        /* retour d'erreur sur lib MPI */

#ifdef _DEBUG
		std::cout << rank()<<") MPI_Probe(...)"<<std::endl;
#endif
	ierror = MPI_Probe(AFromProc,ATag, MPI_COMM_WORLD, &status);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Probe"<<std::endl;

		/* récupération de la taille du message | tableau */
#ifdef _DEBUG
	std::cout << rank()<<") MPI_Get_count(...)"<<std::endl;
#endif
	ierror = MPI_Get_count(&status, MPI_BYTE, &buf_size);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Get_count"<<std::endl;


#ifdef _DEBUG
	std::cout <<" réception d'un tableau de taille = "<<buf_size<<std::endl;
#endif
	char* buffer = new char[buf_size];
//	MPI_Request r;
//	ierror = MPI_Irecv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&r);
	ierror = MPI_Recv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&status);
	if (ierror != MPI_SUCCESS) {
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Recv"<<std::endl;
	}
	int position =0;

	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_faces, 1, MPI_INT,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_edges, 1, MPI_INT,MPI_COMM_WORLD);
	TInt nb_faces = AData.nb_faces;
	TInt nb_edges = AData.nb_edges;
	TInt total_nb_subfaces=0;
	ierror = MPI_Unpack(buffer,buf_size,&position,&total_nb_subfaces, 1, MPI_INT,MPI_COMM_WORLD);

#ifdef _DEBUG
	std::cout<<rank()<<") "<<AData.nb_faces<<" faces received from "<<AFromProc<<std::endl;
#endif //_DEBUG
	AData.init_face_id.resize(nb_faces);
	AData.init_face_nb_sides.resize(nb_faces);
	AData.new_face_node_id.resize(nb_faces);
	AData.new_subface_id.resize(total_nb_subfaces);
	AData.new_edge_node_id.resize(3*nb_edges);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.init_face_id) [0], nb_faces, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.init_face_nb_sides) [0], nb_faces, MPI_INT,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_face_node_id) [0], nb_faces, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_subface_id) [0], total_nb_subfaces, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_edge_node_id) [0], 3*nb_edges, MPI_LONG,MPI_COMM_WORLD);

	delete [] buffer;
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::receivePacked(const int& AFromProc,
											 RefinementStruct& AData,
											 const int& ATag)
{

#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receivePacked("<<AFromProc<<", RefinementStruct)"<<std::endl;
#endif

	int buf_size;
	MPI_Status status;
	int ierror;        /* retour d'erreur sur lib MPI */

#ifdef _DEBUG
		std::cout << rank()<<") MPI_Probe(...)"<<std::endl;
#endif
	ierror = MPI_Probe(AFromProc,ATag, MPI_COMM_WORLD, &status);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Probe"<<std::endl;

		/* récupération de la taille du message | tableau */
#ifdef _DEBUG
	std::cout << rank()<<") MPI_Get_count(...)"<<std::endl;
#endif
	ierror = MPI_Get_count(&status, MPI_BYTE, &buf_size);
	if (ierror != MPI_SUCCESS)
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Get_count"<<std::endl;


#ifdef _DEBUG
	std::cout <<" réception d'un tableau de taille = "<<buf_size<<std::endl;
#endif
	char* buffer = new char[buf_size];
//	MPI_Request r;
//	ierror = MPI_Irecv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&r);
	ierror = MPI_Recv(buffer,buf_size,MPI_PACKED,AFromProc,ATag,MPI_COMM_WORLD,&status);
	if (ierror != MPI_SUCCESS) {
		std::cout<<"Problème sur "<<rank()<<" de communication en provenance de "<<AFromProc<<", lors de MPI_Recv"<<std::endl;
	}
	int position =0;

	ierror = MPI_Unpack(buffer,buf_size,&position,&AData.nb_infos, 1, MPI_INT,MPI_COMM_WORLD);
	TInt nb_infos = AData.nb_infos;
#ifdef _DEBUG
	std::cout<<rank()<<") "<<AData.nb_infos<<" infos received from "<<AFromProc<<std::endl;
#endif //_DEBUG
	AData.init_edge_id.resize(nb_infos);
	AData.new_node_id.resize(nb_infos);
	AData.new_edge1_id.resize(nb_infos);
	AData.new_edge2_id.resize(nb_infos);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.init_edge_id) [0], nb_infos, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_node_id) [0], nb_infos, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_edge1_id) [0], nb_infos, MPI_LONG,MPI_COMM_WORLD);
	ierror = MPI_Unpack(buffer,buf_size,&position,&(AData.new_edge2_id) [0], nb_infos, MPI_LONG,MPI_COMM_WORLD);

	delete [] buffer;
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::receive(int fromProc, std::vector<std::list<id> > & idPairs)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receive("<<fromProc<<", idPairs)"<<std::endl;
#endif

	int tailleTMP = 0;
	void *TabTMP = receive(fromProc, 1, tailleTMP);

#ifdef _DEBUG
	std::cout <<" réception d'un tableau de taille = "<<tailleTMP<<" faite, conversion en vecteur de listes"<<std::endl;
#endif

	char* position = (char*)TabTMP;
	size_t nb_listes;
	memcpy((void*)&(nb_listes), (void*)position, sizeof(size_t)); position+=sizeof(size_t);

	// remplir idPairs
	for (uint i=0; i<nb_listes; i++){
		std::list<id> ll;
		size_t taille_liste;
		memcpy((void*)&(taille_liste), (void*)position, sizeof(size_t)); position+=sizeof(size_t);

		for (uint j=0; j<taille_liste; j++){
			id elem;
			memcpy((void*)&(elem),(void*)position,  sizeof(id)); position+=sizeof(id);
			ll.push_back(elem);
		}
		idPairs.push_back(ll);
	}

	delete [] TabTMP;
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::receive("<<fromProc<<", idPairs) terminé"<<std::endl;
#endif
}
#endif
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
void DistributedMemoryManager::allGather(const std::vector<int>& ASend, std::vector<int>& AGather)
{
#ifdef _DEBUG
	std::cout <<rank()<<") DistributedMemoryManager::allGather(...)"<<std::endl;
#endif

	int localNbElt= ASend.size();
	int globalNbElt=0;



	MPI_Allreduce (&localNbElt, &globalNbElt, 1,
	                   MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

	std::cout<<rank()<<": "<<localNbElt<<" -> "<<globalNbElt<<std::endl;
	// Il faut créer un tableau C  à partir du vecteur C++;


//#ifdef _DEBUG
//	std::cout <<rank()<<") Emission d'un tableau de taille = "<<nbEltToSend<<std::endl;
//#endif
//
////	int *TabTMP = new int [nbEltToSend];
////
////	char* position = TabTMP;
////	size_t nb = idPairs.size();
////	memcpy((void*)position, (void*)&nb, sizeof(size_t)); position+=sizeof(size_t);
////
////	// remplissage de TabTMP suivant idPairs
////	iter = idPairs.begin();
////	for(; iter != idPairs.end(); ++iter){
////		std::list<id> ll = (*iter);
////		nb = ll.size();
////		memcpy((void*)position, (void*)&nb, sizeof(size_t)); position+=sizeof(size_t);
////
////		std::list<id>::iterator iter2 = ll.begin();
////		for(; iter2 != ll.end(); ++iter2){
////			memcpy((void*)position, (void*)&(*iter2), sizeof(id)); position+=sizeof(id);
////		}
////	}
////
//
//	ierror = MPI_Allgather(&nbEltToSend,
//				nbElt,
//				MPI_PACKED,
//				fromProc,
//				tag,
//				MPI_COMM_WORLD);
//
//	MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
//	                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
//	                  MPI_Comm comm)
//	int tailleTMP = 0;
//	void *TabTMP = receive(fromProc, 1, tailleTMP);
//
//#ifdef _DEBUG
//	std::cout <<" réception d'un tableau de taille = "<<tailleTMP<<" faite, conversion en vecteur de listes"<<std::endl;
//#endif
//
//	char* position = (char*)TabTMP;
//	size_t nb_listes;
//	memcpy((void*)&(nb_listes), (void*)position, sizeof(size_t)); position+=sizeof(size_t);
//
//	// remplir idPairs
//	for (uint i=0; i<nb_listes; i++){
//		std::list<id> ll;
//		size_t taille_liste;
//		memcpy((void*)&(taille_liste), (void*)position, sizeof(size_t)); position+=sizeof(size_t);
//
//		for (uint j=0; j<taille_liste; j++){
//			id elem;
//			memcpy((void*)&(elem),(void*)position,  sizeof(id)); position+=sizeof(id);
//			ll.push_back(elem);
//		}
//		idPairs.push_back(ll);
//	}
//
//	delete [] TabTMP;
//#ifdef _DEBUG
//	std::cout <<rank()<<") DistributedMemoryManager::receive("<<fromProc<<", idPairs) terminé"<<std::endl;
//#endif
}
#endif
/*----------------------------------------------------------------------------*/
