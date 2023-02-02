/*----------------------------------------------------------------------------*/
/** \file    VTKWriter.t.h
 *  \author  F. LEDOUX
 *  \date    02/24/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef VTKWRITER_H_
#define VTKWRITER_H_
/*----------------------------------------------------------------------------*/
// headers of STL files
#include <iostream>
#include <fstream>
/*----------------------------------------------------------------------------*/
// headers of VTK files
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
/*----------------------------------------------------------------------------*/
// headers of GMDS
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/Node.h>
#include <GMDS/IG/Edge.h>
#include <GMDS/IG/Face.h>
#include <GMDS/IG/Region.h>
#include <GMDS/Utils/VariableManager.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/IOVTK/VTKWriter_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
VTKWriter<TMesh>::VTKWriter(TMesh& AMesh)
:mesh_(AMesh)
{

	if((AMesh.getDim())==2)
		mesh_dimension_ = 2;
	else 	if((AMesh.getDim())==3)
		mesh_dimension_ = 3;

}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
VTKWriter<TMesh>::~VTKWriter()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKWriter<TMesh>::write(const std::string& AFileName,const int& AMask)
{
	if (mesh_.getModel().has(R) && (AMask|R)==AMask )
		createUnstructuredGrid(AFileName);
	if (mesh_.getModel().has(F) && (AMask|F)==AMask )
		createPolyData(AFileName);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKWriter<TMesh>::createUnstructuredGrid(const std::string& AFileName)
{

	/* the vtk object that stores the mesh to export */
	vtkUnstructuredGrid* grid          = vtkUnstructuredGrid::New();

	/* vtk points to store the mesh nodes */
    vtkPoints*  points  =  vtkPoints::New();

    /* coordinates are double */
    points->SetDataTypeToDouble();
    /* and points are in 3D, thus 3 coordinates*/
    points->GetData()->SetNumberOfComponents(3);

    TInt nbNodes = mesh_.getNbNodes();
    /* memory allocation for the points */
    points->Allocate(nbNodes);

    /* definition of a pointer traversing the data structure that stores
     * points */
    vtkDoubleArray* pointsArray =
                            dynamic_cast<vtkDoubleArray*>(points->GetData());
    double* pointsPtr = pointsArray->WritePointer(0, 3*nbNodes);

    pointsArray->SetNumberOfTuples(nbNodes);


    /* points data structure is now ready to store mesh nodes */
    TInt maxID = mesh_.getMaxLocalID(0);

    /* this id array is useful to keep the connection between initial mesh nodes
     * and the corresponding new VTK points.
     */
    vtkIdTypeArray* IDLimaToVTK = vtkIdTypeArray::New();
    IDLimaToVTK->Allocate(maxID+1);
    IDLimaToVTK->SetNumberOfValues(maxID+1);

    vtkIdTypeArray* point_ids = vtkIdTypeArray::New();
    point_ids->Allocate(nbNodes);
    point_ids->SetName("GMDS_IDS");


    typename TMesh::node_iterator it_nodes = mesh_.nodes_begin();
    int index =0;
    for(;!it_nodes.isDone();it_nodes.next())
    {
    	Node n = it_nodes.value();
		IDLimaToVTK->SetValue(n.getID(),index);
		*pointsPtr++             = n.X();
		*pointsPtr++             = n.Y();
		*pointsPtr++             = n.Z();
		point_ids->InsertNextValue(n.getID());

		index++;
	}

	/* the points are now the points of the grid */
    grid->SetPoints(points);

    /* points was an intermediate data structure that can be deleted now*/
    points->Delete();

    grid->GetPointData()->AddArray(point_ids);
    /* we build regions nows (cells in VTK) */
    vtkCellArray* cellArray     = vtkCellArray::New();

    TInt nbRegions = mesh_.getNbRegions();
    int *vtkCellsTypes = new int [nbRegions];

    /* the number of nodes as the sum of region nodes */
    size_t nodeNumber = 0;

    typename TMesh::region_iterator it_regions     = mesh_.regions_begin();

	for(;!it_regions.isDone();it_regions.next())
        nodeNumber  += it_regions.value().getNbNodes();

	/* Memory Allocation for the grid cells */
    vtkIdTypeArray* idsArray    = vtkIdTypeArray::New();
    idsArray->Allocate(nbRegions + nodeNumber);

    vtkIdType*  cellsPtr = idsArray->GetPointer(0);

    vtkIdTypeArray* cell_ids = vtkIdTypeArray::New();
    cell_ids->Allocate(nbRegions);
    cell_ids->SetName("GMDS_IDS");

    index = 0;
	for(it_regions = mesh_.regions_begin(); !it_regions.isDone();it_regions.next())
    {
		Region r=it_regions.value();
		switch(r.getType()){
		case GMDS_HEX:
			vtkCellsTypes [index]   = VTK_HEXAHEDRON;
			{
				*cellsPtr++ = r.getNbNodes();
				std::vector<TCellID> node_ids = r.getIDs<Node>();
				for (unsigned int n = 0; n < node_ids.size(); n++)
					*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[n]);
			}
			break;
		case  GMDS_TETRA:
			vtkCellsTypes [index]   = VTK_TETRA;
			{
				*cellsPtr++ = r.getNbNodes();
				std::vector<TCellID> node_ids = r.getIDs<Node>();
				for (unsigned int n = 0; n < node_ids.size(); n++)
					*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[n]);
			}
			break;
		case  GMDS_PRISM3:
			vtkCellsTypes [index]   = VTK_WEDGE;
			{
				*cellsPtr++ = r.getNbNodes();
				std::vector<TCellID> node_ids = r.getIDs<Node>();
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[0]);
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[2]);
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[1]);
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[3]);
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[5]);
				*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[4]);
			}
			break;
		case  GMDS_PYRAMID:
			vtkCellsTypes [index]   = VTK_PYRAMID;
			{
				*cellsPtr++ = r.getNbNodes();
				std::vector<TCellID> node_ids = r.getIDs<Node>();
				for (unsigned int n = 0; n < node_ids.size(); n++)
					*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[n]);
			}
			break;
		default:
			throw GMDSException("This mesh has cell types not handled by VTK");
		}

		cell_ids->InsertNextValue(r.getID());
		index++;
    }

    idsArray->SetNumberOfValues(nbRegions+nodeNumber);
    cellArray->SetCells (nbRegions, idsArray);
    cellArray->Squeeze();


    grid->GetCellData()->AddArray(cell_ids);

    grid->SetCells (vtkCellsTypes, cellArray);
    /*==========================================================================*/
    /* Variables management*/
    /*==========================================================================*/
    //			VARIABLES ON NODES
    /*==========================================================================*/
    VariableManager* var_nodes = &(mesh_.m_node_variable_manager);
    TInt nb_nodes_variables = var_nodes->getNbVariables();
    for(TInt i_var=0;i_var<nb_nodes_variables;i_var++){
    	VariableItf* var_i = var_nodes->getVariable(i_var);
    	Variable<TInt>* var_int = dynamic_cast<Variable<TInt>* >(var_i);
    	if(var_int!=0)
    	{
    		vtkIntArray* cloud_array = vtkIntArray::New();
    		cloud_array->Allocate(nbNodes);
    		cloud_array->SetName(var_int->getName().c_str());

    		typename TMesh::node_iterator it;
    		int index =0;
    		for(it = mesh_.nodes_begin();!it.isDone();it.next())
    		{
    			Node n= it.value();
    			cloud_array->InsertNextValue((int)((*var_int)[n.getID()]));
    		}
    		grid->GetPointData()->AddArray(cloud_array);
    		cloud_array->Delete();
    	}
    	else{
    		Variable<double>* var_double = dynamic_cast<Variable<double>* >(var_i);
    		if(var_double!=0)
    		{

    			vtkDoubleArray* cloud_array = vtkDoubleArray::New();
    			cloud_array->Allocate(nbNodes);
    			cloud_array->SetName(var_double->getName().c_str());

    			typename TMesh::node_iterator it;
    			int index =0;
    			for(it = mesh_.nodes_begin();!it.isDone();it.next())
    			{
    				Node n= it.value();
    				cloud_array->InsertNextValue((double)((*var_double)[n.getID()]));
    			}
    			grid->GetPointData()->AddArray(cloud_array);
    			cloud_array->Delete();
    		}
        	else{
        		Variable<math::Vector>* var_vect = dynamic_cast<Variable<math::Vector>* >(var_i);
        		if(var_vect!=0)
        		{
        			vtkDoubleArray* cloud_array = vtkDoubleArray::New();
        			cloud_array->SetNumberOfComponents(3);
        			cloud_array->Allocate(nbNodes);
        			cloud_array->SetName(var_vect->getName().c_str());

        			typename TMesh::node_iterator it;
        			int index =0;
        			for(it = mesh_.nodes_begin();!it.isDone();it.next())
        			{
        				Node n= it.value();
        				math::Vector v = (*var_vect)[n.getID()];
        				double x = v.X();
        				double y = v.Y();
        				double z = v.Z();
        				cloud_array->InsertNextTuple3(x,y,z);
        			}
        			grid->GetPointData()->AddArray(cloud_array);
        			cloud_array->Delete();
        		}
        	}
    	}
    }
    /*==========================================================================*/
    //			VARIABLES ON REGIONS
    /*==========================================================================*/
    VariableManager* var_regions = &(mesh_.m_region_variable_manager);
    TInt nb_regions_variables = var_regions->getNbVariables();
    for(TInt i_var=0;i_var<nb_regions_variables;i_var++){
    	VariableItf* var_i = var_regions->getVariable(i_var);
    	Variable<TInt>* var_int = dynamic_cast<Variable<TInt>* >(var_i);
    	if(var_int!=0)
    	{
    		vtkIntArray* vol_array = vtkIntArray::New();
    		vol_array->Allocate(nbRegions);
    		vol_array->SetName(var_int->getName().c_str());

    		typename TMesh::region_iterator it;
    		int index =0;
    		for(it = mesh_.regions_begin();!it.isDone();it.next())
    		{
    			Region r= it.value();
    			vol_array->InsertNextValue((int)((*var_int)[r.getID()]));
    		}
    		grid->GetCellData()->AddArray(vol_array);
    		vol_array->Delete();
    	}
    	else{
    		Variable<double>* var_double = dynamic_cast<Variable<double>* >(var_i);
        	if(var_double!=0)
        	{
        		vtkDoubleArray* vol_array = vtkDoubleArray::New();
        		vol_array->Allocate(nbRegions);
        		vol_array->SetName(var_double->getName().c_str());

        		typename TMesh::region_iterator it;
        		int index =0;
        		for(it = mesh_.regions_begin();!it.isDone();it.next())
        		{
        			Region r = it.value();
        			vol_array->InsertNextValue((double)((*var_double)[r.getID()]));
        		}
        		grid->GetCellData()->AddArray(vol_array);
        		vol_array->Delete();
        	}
        	else
        	{
        		Variable<math::Vector>* var_vect = dynamic_cast<Variable<math::Vector>* >(var_i);
        		if(var_vect!=0)
        		{
        			vtkDoubleArray* vol_array = vtkDoubleArray::New();
        			vol_array->SetNumberOfComponents(3);
        			vol_array->Allocate(nbNodes);
        			vol_array->SetName(var_vect->getName().c_str());
        			typename TMesh::region_iterator it;
        			int index =0;
        			for(it = mesh_.regions_begin();!it.isDone();it.next())
        			{
        				Region r = it.value();
        				math::Vector v = (*var_vect)[r.getID()];
        				double x = v.X();
        				double y = v.Y();
        				double z = v.Z();
        				vol_array->InsertNextTuple3(x,y,z);
        			}
        			grid->GetCellData()->AddArray(vol_array);
        			vol_array->Delete();

        		}
        	}
    	}
    }

    /* now we create data array for recognizing clouds*/
//	std::vector<bool> nodes_in_cloud;
//	nodes_in_cloud.reserve(maxID+1);
//
//	typename Mesh<TMask>::clouds_iterator it_clouds = mesh_.clouds_begin();
//	for(;it_clouds!=mesh_.clouds_end();it_clouds++)
//	{
//
//		typename Mesh<TMask>::cloud& cl = *it_clouds;
//
//		for(unsigned int i=0; i<maxID+1;i++)
//			nodes_in_cloud[i]=0;
//
//		if(cl.size())
//		{
//
//
//		    std::vector<Node*>& nodes_in_cl = cl.cells();
//			for(unsigned int n_index=0; n_index<nodes_in_cl.size();n_index++)
//				nodes_in_cloud[nodes_in_cl[n_index]->getID()] = 1;
//
//
//		    vtkIntArray* cloud_array = vtkIntArray::New();
//		    cloud_array->Allocate(nbNodes);
//		    cloud_array->SetName(cl.name().c_str());
//
//			typename Mesh<TMask>::nodes_iterator it=mesh_.nodes_begin();
//
//			for(; !it->isDone();it->next())
//				cloud_array->InsertNextValue(nodes_in_cloud[it->currentItem()->getID()]);
//
//			grid->GetPointData()->AddArray(cloud_array);
//
//			cloud_array->Delete();
//		}
//	}
//
//    /* now we create data array for recognizing volumes*/
//	std::vector<bool> regions_in_volume;
//	int max_size = mesh_.getMaxLocalID(3);
//	regions_in_volume.reserve(max_size+1);
//
//	typename Mesh<TMask>::volumes_iterator it_vols = mesh_.volumes_begin();
//	for(;it_vols!=mesh_.volumes_end();it_vols++)
//	{
//
//		typename Mesh<TMask>::volume& vol = *it_vols;
//
//
//		for(unsigned int i=0; i<max_size+1;i++)
//			regions_in_volume[i]=0;
//
//		if(vol.size())
//		{
//
//
//		    std::vector<Region*>& regions_in_vol = vol.cells();
//			for(unsigned int r_index=0; r_index<regions_in_vol.size();r_index++)
//				regions_in_volume[regions_in_vol[r_index]->getID()] = 1;
//
//
//		    vtkIntArray* vol_array = vtkIntArray::New();
//		    vol_array->Allocate(nbRegions);
//		    vol_array->SetName(vol.name().c_str());
//
//			typename Mesh<TMask>::regions_iterator it=mesh_.regions_begin();
//			for(;!it->isDone();it->next())
//				vol_array->InsertNextValue(regions_in_volume[it->currentItem()->getID()]);
//
//			grid->GetCellData()->AddArray(vol_array);
//
//			vol_array->Delete();
//		}
//	}
//


    /* deletion of the intermediate data structures */
    delete [] vtkCellsTypes;
    idsArray->Delete();
    IDLimaToVTK->Delete();
    cellArray->Delete();

    point_ids->Delete();
    cell_ids->Delete();




    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetInput(grid);
    std::string file_name = AFileName+".vtu";
    writer->SetFileName(file_name.c_str());
    writer->SetDataModeToAscii();
    writer->Write();
    writer->Delete();


    /* deletion of the vtk object*/
    grid->Delete();
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void VTKWriter<TMesh>::createPolyData(const std::string& AFileName)
{

	/* the vtk object that stores the mesh to export */
	vtkPolyData* polydata          = vtkPolyData::New();

	/* vtk points to store the mesh nodes */
    vtkPoints*  points  =  vtkPoints::New();

    /* coordinates are double */
   points->SetDataTypeToDouble();
    /* and points are in 3D, thus 3 coordinates*/
    points->GetData()->SetNumberOfComponents(3);

    TInt nbNodes = mesh_.getNbNodes();
    /* memory allocation for the points */
    points->Allocate(nbNodes);

    /* definition of a pointer traversing the data structure that stores
     * points */
    vtkDoubleArray* pointsArray =
                            dynamic_cast<vtkDoubleArray*>(points->GetData());
    double* pointsPtr = pointsArray->WritePointer(0, 3*nbNodes);

    pointsArray->SetNumberOfTuples(nbNodes);


    /* points data structure is now ready to store mesh nodes */
    TInt maxID = mesh_.getMaxLocalID(0);

    /* this id array is useful to keep the connection between initial mesh nodes
     * and the corresponding new VTK points.
     */
    vtkIdTypeArray* IDLimaToVTK = vtkIdTypeArray::New();
    IDLimaToVTK->Allocate(maxID+1);
    IDLimaToVTK->SetNumberOfValues(maxID+1);

    vtkIdTypeArray* point_ids = vtkIdTypeArray::New();
    point_ids->Allocate(nbNodes);
    point_ids->SetName("GMDS_IDS");


    typename TMesh::node_iterator itn;
    int index =0;

    for(itn = mesh_.nodes_begin();!itn.isDone();itn.next())
    {
    	Node n= itn.value();

		IDLimaToVTK->SetValue(n.getID(),index);
		*pointsPtr++             = n.X();
		*pointsPtr++             = n.Y();
		*pointsPtr++             = n.Z();

		point_ids->InsertNextValue(n.getID());

		index++;
	}
    int nb_nodes = points->GetNumberOfPoints();
	/* the points are now the points of the polydata */
    polydata->SetPoints(points);

    /* points was an intermediate data structure that can be deleted now*/
    points->Delete();

    polydata->GetPointData()->AddArray(point_ids);
    /* we build face nows (cells in VTK) */
    vtkCellArray* cellArray     = vtkCellArray::New();

    TInt nbFaces = mesh_.getNbFaces();


    /* the number of nodes as the sum of face nodes */
    size_t nodeNumber = 0;
    typename TMesh::face_iterator itf;

    for(itf = mesh_.faces_begin();!itf.isDone();itf.next())
    {
        nodeNumber  += itf.value().getNbNodes();
	}

	/* Memory Allocation for the polydata cells */
    vtkIdTypeArray* idsArray    = vtkIdTypeArray::New();
    idsArray->Allocate(nbFaces + nodeNumber);

    vtkIdType*  cellsPtr = idsArray->GetPointer(0);

    vtkIdTypeArray* cell_ids = vtkIdTypeArray::New();
    cell_ids->Allocate(nbFaces);
    cell_ids->SetName("GMDS_IDS");

    index = 0;
    for(itf = mesh_.faces_begin();!itf.isDone();itf.next())
    {
		Face f=itf.value();

		*cellsPtr++ = f.getNbNodes();
		std::vector<TCellID> node_ids;
		f.getIDs<Node>(node_ids);
		for (unsigned int n = 0; n < node_ids.size(); n++)
			*cellsPtr++ = IDLimaToVTK->GetValue(node_ids[n]);

		cell_ids->InsertNextValue(f.getID());
		index++;
    }

    idsArray->SetNumberOfValues(nbFaces+nodeNumber);
    cellArray->SetCells (nbFaces, idsArray);
    cellArray->Squeeze();


    polydata->GetCellData()->AddArray(cell_ids);

    polydata->SetPolys(cellArray);


    /*==========================================================================*/
    /* Variables management*/
    /*==========================================================================*/
    //			VARIABLES ON NODES
    /*==========================================================================*/
    VariableManager* var_nodes = &(mesh_.m_node_variable_manager);
    TInt nb_nodes_variables = var_nodes->getNbVariables();
    for(TInt i_var=0;i_var<nb_nodes_variables;i_var++){
    	VariableItf* var_i = var_nodes->getVariable(i_var);
    	Variable<TInt>* var_int = dynamic_cast<Variable<TInt>* >(var_i);
    	if(var_int!=0)
    	{
    		vtkIntArray* cloud_array = vtkIntArray::New();
    		cloud_array->Allocate(nbNodes);
    		cloud_array->SetName(var_int->getName().c_str());

    		typename TMesh::node_iterator it;
    		int index =0;
    		for(it = mesh_.nodes_begin();!it.isDone();it.next())
    		{
    			Node n= it.value();
    			cloud_array->InsertNextValue((int)((*var_int)[n.getID()]));
    		}
    		polydata->GetPointData()->AddArray(cloud_array);
    		cloud_array->Delete();
    	}
    	else{
    		Variable<double>* var_double = dynamic_cast<Variable<double>* >(var_i);
    		if(var_double!=0)
    		{

    			vtkDoubleArray* cloud_array = vtkDoubleArray::New();
    			cloud_array->Allocate(nbNodes);
    			cloud_array->SetName(var_double->getName().c_str());

    			typename TMesh::node_iterator it;
    			int index =0;
    			for(it = mesh_.nodes_begin();!it.isDone();it.next())
    			{
    				Node n= it.value();
    				cloud_array->InsertNextValue((double)((*var_double)[n.getID()]));
    			}
    			polydata->GetPointData()->AddArray(cloud_array);
    			cloud_array->Delete();
    		}
    		else{
    			Variable<math::Vector>* var_vect = dynamic_cast<Variable<math::Vector>* >(var_i);
    			if(var_vect!=0)
    			{
    				vtkDoubleArray* cloud_array = vtkDoubleArray::New();
    				cloud_array->SetNumberOfComponents(3);
    				cloud_array->Allocate(nbNodes);
    				cloud_array->SetName(var_vect->getName().c_str());

    				typename TMesh::node_iterator it;
    				int index =0;
    				for(it = mesh_.nodes_begin();!it.isDone();it.next())
    				{
    					Node n= it.value();
    					math::Vector v = (*var_vect)[n.getID()];
    					double x = v.X();
    					double y = v.Y();
    					double z = v.Z();
    					cloud_array->InsertNextTuple3(x,y,z);
    				}
    				polydata->GetPointData()->AddArray(cloud_array);
    				cloud_array->Delete();
    			}
    		}
    	}
    }
    /*==========================================================================*/
    //			VARIABLES ON FACES
    /*==========================================================================*/
    VariableManager* var_faces = &(mesh_.m_face_variable_manager);
    TInt nb_faces_variables = var_faces->getNbVariables();
    for(TInt i_var=0;i_var<nb_faces_variables;i_var++){
    	VariableItf* var_i = var_faces->getVariable(i_var);
    	Variable<TInt>* var_int = dynamic_cast<Variable<TInt>* >(var_i);
    	if(var_int!=0)
    	{
    		vtkIntArray* surf_array = vtkIntArray::New();
    		surf_array->Allocate(nbFaces);
    		surf_array->SetName(var_int->getName().c_str());

    		typename TMesh::face_iterator it;
    		int index =0;
    		for(it = mesh_.faces_begin();!it.isDone();it.next())
    		{
    			Face f= it.value();
    			surf_array->InsertNextValue((int)((*var_int)[f.getID()]));
    		}
    		polydata->GetCellData()->AddArray(surf_array);
    		surf_array->Delete();
    	}
    	else{
    		Variable<double>* var_double = dynamic_cast<Variable<double>* >(var_i);
        	if(var_double!=0)
        	{
        		vtkDoubleArray* surf_array = vtkDoubleArray::New();
        		surf_array->Allocate(nbFaces);
        		surf_array->SetName(var_double->getName().c_str());

        		typename TMesh::face_iterator it;
        		int index =0;
        		for(it = mesh_.faces_begin();!it.isDone();it.next())
        		{
        			Face f= it.value();
        			surf_array->InsertNextValue((double)((*var_double)[f.getID()]));
        		}
        		polydata->GetCellData()->AddArray(surf_array);
        		surf_array->Delete();
        	}

        	else{
        		Variable<math::Vector>* var_vect = dynamic_cast<Variable<math::Vector>* >(var_i);
        		if(var_vect!=0)
        		{
        			vtkDoubleArray* surf_array = vtkDoubleArray::New();
        			surf_array->SetNumberOfComponents(3);
        			surf_array->Allocate(nbNodes);
        			surf_array->SetName(var_vect->getName().c_str());

        			typename TMesh::face_iterator it;
        			int index =0;
        			for(it = mesh_.faces_begin();!it.isDone();it.next())
        			{
        				Face f= it.value();
        				math::Vector v = (*var_vect)[f.getID()];
        				double x = v.X();
        				double y = v.Y();
        				double z = v.Z();
        				surf_array->InsertNextTuple3(x,y,z);
        			}
        			polydata->GetPointData()->AddArray(surf_array);
        			surf_array->Delete();
        		}
        	}
    	}
    }
    //    /* now we create data array for node variables*/
	std::vector<bool> nodes_in_cloud;
	nodes_in_cloud.resize(maxID+1);


	typename TMesh::clouds_iterator it_clouds = mesh_.clouds_begin();
	for(;it_clouds!=mesh_.clouds_end();it_clouds++)
	{

		typename TMesh::cloud& cl = *it_clouds;

		for(TInt i=0; i<maxID+1;i++)
			nodes_in_cloud[i]=0;

		if(cl.size())
		{
		    std::vector<TCellID>& nodes_in_cl = cl.cellIDs();
			for(unsigned int n_index=0; n_index<nodes_in_cl.size();n_index++)
				nodes_in_cloud[nodes_in_cl[n_index]] = 1;


		    vtkIntArray* cloud_array = vtkIntArray::New();
		    cloud_array->Allocate(nbNodes);
		    cloud_array->SetName(cl.name().c_str());

		    typename TMesh::node_iterator it = mesh_.nodes_begin();
			for(; !it.isDone();it.next())
				cloud_array->InsertNextValue(nodes_in_cloud[it.value().getID()]);

			polydata->GetPointData()->AddArray(cloud_array);

			cloud_array->Delete();
		}
	}



    /* now we create data array for recognizing surfaces*/
	std::vector<bool> faces_in_surface;
	TInt max_size = mesh_.getMaxLocalID(2);
	faces_in_surface.resize(max_size+1);


	typename TMesh::surfaces_iterator it_surfs = mesh_.surfaces_begin();
	for(;it_surfs!=mesh_.surfaces_end();it_surfs++)
	{

		typename TMesh::surface& surf = *it_surfs;

		for(TInt i=0; i<max_size+1;i++)
			faces_in_surface[i]=0;


		if(surf.size())
		{


		    std::vector<TCellID>& faces_in_surf = surf.cellIDs();
			for(unsigned int r_index=0; r_index<faces_in_surf.size();r_index++)
				faces_in_surface[faces_in_surf[r_index]] = 1;


		    vtkIntArray* surf_array = vtkIntArray::New();
		    surf_array->Allocate(nbFaces);
		    surf_array->SetName(surf.name().c_str());

		    typename  TMesh::face_iterator it = mesh_.faces_begin();
			for(; !it.isDone();it.next())
				surf_array->InsertNextValue(faces_in_surface[it.value().getID()]);

			polydata->GetCellData()->AddArray(surf_array);

			surf_array->Delete();
		}
	}



    /* deletion of the intermediate data structures */
    idsArray->Delete();
    IDLimaToVTK->Delete();
    cellArray->Delete();

    point_ids->Delete();
    cell_ids->Delete();

    vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
    writer->SetInput(polydata);
    writer->SetDataModeToBinary();
    std::string file_name = AFileName+".vtp";
    writer->SetFileName(file_name.c_str());
    //writer->SetDataModeToAscii();
    writer->Write();
    writer->Delete();


    /* deletion of the vtk object*/
    polydata->Delete();
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //VTKWRITER_H_
/*----------------------------------------------------------------------------*/
