/*----------------------------------------------------------------------------*/
/** \file    VTKReader.t.h
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef VTKREADER_H_
#define VTKREADER_H_
/*----------------------------------------------------------------------------*/
// headers of STL files
#include <iostream>
#include <fstream>
#include <map>
/*----------------------------------------------------------------------------*/
// headers of VTK files
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IO/IReader.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/IOVTK/VTKReader_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
VTKReader<TMesh>::VTKReader(TMesh& AMesh)
:IReader<TMesh>(AMesh)
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
VTKReader<TMesh>::~VTKReader()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKReader<TMesh>::read(const std::string& AFileName)
{
    vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(AFileName.c_str());

	if(reader->CanReadFile(AFileName.c_str()))
	{
		read(reader);
	}
	else
	{
		vtkXMLPolyDataReader *reader2 = vtkXMLPolyDataReader::New();
		reader2->SetFileName(AFileName.c_str());
		if(reader2->CanReadFile(AFileName.c_str()))
		{
			read(reader2);
		}
		else
		{
			reader->Delete();
			reader2->Delete();
			throw GMDSException("Unable to read this file with VTK in GMDS");
		}
		reader2->Delete();
	}
	reader->Delete();

}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKReader<TMesh>::read(vtkXMLUnstructuredGridReader* AReader)
{

	vtkUnstructuredGrid* grid = 0;
    grid = AReader->GetOutput();

    grid->Update();

    int nb_regions = grid->GetNumberOfCells();
    int nb_nodes   = grid->GetNumberOfPoints();

    /* vtk points to store the mesh nodes */
    vtkPoints*  points  =  grid->GetPoints();
    vtkPointData* pointData = grid->GetPointData();
    vtkCellData* cellData = grid->GetCellData();



    /* these maps are created to store the link between node and cell position in
     * the vtk list and their gmds ids.  */
    std::map<int,int> VTK_to_GMDS_Node_Ids;
    std::map<int,int> VTK_to_GMDS_Cell_Ids;

    vtkIntArray* point_ids = dynamic_cast<vtkIntArray*>(pointData->GetArray("GMDS_IDS"));
    vtkIntArray* cell_ids  = dynamic_cast<vtkIntArray*>(cellData->GetArray("GMDS_IDS"));

    if(point_ids!=0)
    {
    	int max_VTK_node_id =point_ids->GetSize()-1;
    	/* in this case we need to get the max id */
    	for(int i=0;i<point_ids->GetSize();i++)
    		if(point_ids->GetValue(i)>max_VTK_node_id)
    			max_VTK_node_id = point_ids->GetValue(i);
    	/* now we can allocate the necessary space in the mesh node container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxNodeID(max_VTK_node_id);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_nodes;i++)
			VTK_to_GMDS_Node_Ids[i]=point_ids->GetValue(i);
    }
    else
    {
    	/* now we can allocate the necessary space in the mesh node container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxNodeID(grid->GetNumberOfPoints());

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_nodes;i++)
			VTK_to_GMDS_Node_Ids[i]=i;
    }

    if(cell_ids!=0)
    {
    	int max_VTK_cell_id =cell_ids->GetSize()-1;

    	/* in this case we need to get the max id */
    	for(int i=0;i<cell_ids->GetSize();i++)
    		if(cell_ids->GetValue(i)>max_VTK_cell_id)
    			max_VTK_cell_id = cell_ids->GetValue(i);
    	/* now we can allocate the necessary space in the mesh region container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxRegionID(max_VTK_cell_id);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_regions;i++)
			VTK_to_GMDS_Cell_Ids[i]=cell_ids->GetValue(i);
    }
    else
    {
    	/* now we can allocate the necessary space in the mesh region container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxRegionID(grid->GetNumberOfCells());

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_regions;i++)
			VTK_to_GMDS_Cell_Ids[i]=i;
    }

    /* we create an id list that is necessary for getting node ids from cells*/
	vtkIdList* cell_point_ids = vtkIdList::New();

	/* Creation of the GMDS nodes */
	for(int i=0;i<nb_nodes;i++)
	{
		double coord[3];
		points->GetPoint(i,coord);
		this->newNode(coord[0],coord[1],coord[2],VTK_to_GMDS_Node_Ids[i]);
	}

	/* Creation of the GMDS regions */
	for(int i=0; i<nb_regions;i++){

		grid->GetCellPoints(i,cell_point_ids);
		switch(grid->GetCellType(i)){
		case VTK_TETRA:
			this->newTet(VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(0)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(1)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(2)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(3)],
						 VTK_to_GMDS_Cell_Ids[i]);
			break;
		case VTK_HEXAHEDRON:
			this->newHex(VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(0)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(1)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(2)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(3)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(4)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(5)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(6)],
						 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(7)],
						 VTK_to_GMDS_Cell_Ids[i]);
			break;
		case VTK_PYRAMID:
			this->newPyramid(VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(0)],
							 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(1)],
							 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(2)],
							 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(3)],
							 VTK_to_GMDS_Node_Ids[cell_point_ids->GetId(4)],
							 VTK_to_GMDS_Cell_Ids[i]);
			break;
		default:
			throw GMDSException("Not supported type in the VTK Cell reader!");

		}
	}

	/* now we create the node and regions groups that are present in the VTK
	 * file*/
	for(int i=0; i<pointData->GetNumberOfArrays();i++)
	{
		vtkIntArray* array = dynamic_cast<vtkIntArray*>(pointData->GetArray(i));
		if(array!=NULL){
			std::string name = array->GetName();
			if(name!="GMDS_IDS")
			{
				typename TMesh::cloud& c = this->mesh_.newCloud(array->GetName());
				for(int j=0;j<array->GetSize();j++)
					if(array->GetValue(j)!=0)
						c.add(this->mesh_.template get<typename TMesh::node>(VTK_to_GMDS_Node_Ids[j]));
			}
		}
	}

	for(int i=0; i<cellData->GetNumberOfArrays();i++)
	{
		vtkIntArray* array = dynamic_cast<vtkIntArray*>(cellData->GetArray(i));
		if(array!=NULL){
			std::string name = array->GetName();
			if(name!="GMDS_IDS")
			{
				typename TMesh::volume& v = this->mesh_.newVolume(array->GetName());
				for(int j=0;j<array->GetSize();j++)
					if(array->GetValue(j)!=0)
						v.add(this->mesh_.template get<typename TMesh::region>(VTK_to_GMDS_Cell_Ids[j]));
			}
		}
	}

	cell_point_ids->Delete();
    this->updateMeshIDContainers();
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKReader<TMesh>::read(vtkXMLPolyDataReader* AReader)
{
	vtkPolyData* polydata= 0;

    //polydata->DeepCopy(AReader->GetOutput());

    polydata = AReader->GetOutput();

    polydata->Update();

    int nb_faces = polydata->GetNumberOfCells();
    int nb_nodes = polydata->GetNumberOfPoints();


    /* vtk points to store the mesh nodes */
    vtkPoints*  points  =  polydata->GetPoints();
    vtkPointData* pointData = polydata->GetPointData();
    vtkCellData* cellData = polydata->GetCellData();

    /* these maps are created to store the link between node and cell position in
     * the vtk list and their gmds ids.  */
    std::map<vtkIdType,int> VTK_to_GMDS_Node_Ids;
    std::map<vtkIdType,int> VTK_to_GMDS_Cell_Ids;

    vtkIntArray* point_ids = dynamic_cast<vtkIntArray*>(pointData->GetArray("GMDS_IDS"));
    vtkIntArray* cell_ids  = dynamic_cast<vtkIntArray*>(cellData->GetArray("GMDS_IDS"));

    if(point_ids!=0)
    {
    	int max_VTK_node_id =point_ids->GetSize()-1;
    	/* in this case we need to get the max id */
    	for(int i=0;i<point_ids->GetSize();i++)
    		if(point_ids->GetValue(i)>max_VTK_node_id)
    			max_VTK_node_id = point_ids->GetValue(i);
    	/* now we can allocate the necessary space in the mesh node container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxNodeID(max_VTK_node_id+1);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_nodes;i++)
			VTK_to_GMDS_Node_Ids[i]=point_ids->GetValue(i);
    }
    else
    {
    	/* now we can allocate the necessary space in the mesh node container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxNodeID(polydata->GetNumberOfPoints()+1);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_nodes;i++)
			VTK_to_GMDS_Node_Ids[i]=i;
    }

    if(cell_ids!=0)
    {
        int max_VTK_cell_id =cell_ids->GetSize()-1;
    	/* in this case we need to get the max id */
    	for(int i=0;i<cell_ids->GetSize();i++)
    		if(cell_ids->GetValue(i)>max_VTK_cell_id)
    			max_VTK_cell_id = cell_ids->GetValue(i);
    	/* now we can allocate the necessary space in the mesh region container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxFaceID(max_VTK_cell_id+1);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_faces;i++)
			VTK_to_GMDS_Cell_Ids[i]=cell_ids->GetValue(i);
    }
    else
    {
    	/* now we can allocate the necessary space in the mesh region container to allocate
    	 * node with id without any trouble.	 */
    	this->specifyMaxFaceID(polydata->GetNumberOfCells()+1);

    	/* VTK location -> gmds ids */
    	for(int i=0;i<nb_faces;i++)
			VTK_to_GMDS_Cell_Ids[i]=i;
    }

	/* Creation of the GMDS nodes */
	for(int i=0;i<nb_nodes;i++)
	{
		double coord[3];
		points->GetPoint(i,coord);
		this->newNode(coord[0],coord[1],coord[2],VTK_to_GMDS_Node_Ids[i]);
	}


	vtkIdList* cellPtIds = vtkIdList::New();
	for(int i=0; i<nb_faces;i++){
		/* we create an id list that is necessary for getting node ids from cells*/
		int id_gmds = VTK_to_GMDS_Cell_Ids[i];
		int cell_type =polydata->GetCellType(i);

		vtkCell* cell_i = polydata->GetCell(i);

		polydata->GetCellPoints(i,cellPtIds);

		switch(cell_type){
		case VTK_QUAD:
			this->newQuad(VTK_to_GMDS_Node_Ids[cellPtIds->GetId(0)],
						  VTK_to_GMDS_Node_Ids[cellPtIds->GetId(1)],
						  VTK_to_GMDS_Node_Ids[cellPtIds->GetId(2)],
						  VTK_to_GMDS_Node_Ids[cellPtIds->GetId(3)],
						  VTK_to_GMDS_Cell_Ids[i]);
			break;
		case VTK_TRIANGLE:{
			vtkIdType i0 = cellPtIds->GetId(0);
			vtkIdType i1 = cellPtIds->GetId(1);
			vtkIdType i2 = cellPtIds->GetId(2);

			this->newTriangle(VTK_to_GMDS_Node_Ids[cellPtIds->GetId(0)],
							  VTK_to_GMDS_Node_Ids[cellPtIds->GetId(1)],
							  VTK_to_GMDS_Node_Ids[cellPtIds->GetId(2)],
							  VTK_to_GMDS_Cell_Ids[i]);
		}
			break;
		case VTK_POLYGON:{
			std::vector<TCellID> ids;
			for(int k=0;k<cellPtIds->GetNumberOfIds();k++)
				ids.push_back(VTK_to_GMDS_Node_Ids[cellPtIds->GetId(k)]);

			this->newPolygon(ids,VTK_to_GMDS_Cell_Ids[i]);
			}
			break;
		default:
			throw GMDSException("Not supported type in the VTK Cell reader!");

		}

	}

	/* now we create the node and face groups that are present in the VTK
	 * file*/
	for(int i=0; i<pointData->GetNumberOfArrays();i++)
	{
		vtkIntArray* array = dynamic_cast<vtkIntArray*>(pointData->GetArray(i));
		if(array!=NULL){
			std::string name = array->GetName();
			if(name!="GMDS_IDS")
			{
				typename TMesh::cloud& c = this->mesh_.newCloud(array->GetName());
				for(int j=0;j<array->GetSize();j++)
					if(array->GetValue(j)!=0)
						c.add(this->mesh_.template get<typename TMesh::node>(VTK_to_GMDS_Node_Ids[j]));
			}
		}
	}

	for(int i=0; i<cellData->GetNumberOfArrays();i++)
	{
		vtkIntArray* array = dynamic_cast<vtkIntArray*>(cellData->GetArray(i));
		if(array!=NULL){
			std::string name = array->GetName();
			if(name!="GMDS_IDS")
			{
				typename TMesh::surface& s = this->mesh_.newSurface(array->GetName());
				for(int j=0;j<array->GetSize();j++)
					if(array->GetValue(j)!=0)
						s.add(this->mesh_.template get<typename TMesh::face>(VTK_to_GMDS_Cell_Ids[j]));
			}
		}
	}

	cellPtIds->Delete();

    this->updateMeshIDContainers();

}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //VTKREADER_H_
/*----------------------------------------------------------------------------*/
