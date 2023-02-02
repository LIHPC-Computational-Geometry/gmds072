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
/** \file    VTKReader.h
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VTKREADER_H_
#define GMDS_VTKREADER_H_
/*----------------------------------------------------------------------------*/
// headers of GMDS files
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <GMDS/IG/IG.h>
#include <GMDS/IO/IReader.h>
/*----------------------------------------------------------------------------*/
// headers of STL files
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <string>
/*----------------------------------------------------------------------------*/
namespace gmds{
	/*----------------------------------------------------------------------------*/
#include <GMDS/IO/VTKReader_def.h>
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	VTKReader<TMesh>::VTKReader(TMesh& AMesh)
		:IReader<TMesh>(AMesh), m_nb_imported_faces(0), m_nb_imported_regions(0)
	{}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	VTKReader<TMesh>::~VTKReader()
	{}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::read(const std::string& AFileName)
	{
		std::ifstream input(AFileName.c_str(), std::ios::in);
		if (!input)
			throw GMDSException("Impossible to read this VTK file");

		const int max_length = 60;
		char ch[max_length];
		int lg;
		bool fileTypeFound = false;
		bool isPolydata = false;
		do{
			input.getline(ch, max_length);
			lg = (int)(input.gcount());
			
			if(std::string(ch).find("DATASET UNSTRUCTURED_GRID") != std::string::npos) {
				isPolydata = false;
				fileTypeFound = true;
			}
			if(std::string(ch).find("DATASET POLYDATA") != std::string::npos) {
                                isPolydata = true;
				fileTypeFound = true;
                        }

		} while (
			!fileTypeFound && lg > 0);

		if (!fileTypeFound)
			throw GMDSException("Impossible to read this VTK file. Only unstructured grid are read.");

		std::string current_word;
		do{
			input >> current_word;
		} while (current_word != "POINTS");

		readNodes(input);

		do{
			input >> current_word;
		} while ( (!isPolydata && current_word != "CELLS") || (isPolydata && current_word != "POLYGONS"));

		readCells(input,isPolydata);

		readNodesData(AFileName);

		if (m_nb_imported_faces != 0 && m_nb_imported_regions == 0)
			readFacesData(AFileName);
		else if (m_nb_imported_faces == 0 && m_nb_imported_regions != 0)
			readRegionsData(AFileName);
		else
			throw GMDSException("Impossible to read a VTK file containing both faces and regions");

		input.close();


	}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::readNodes(std::ifstream& AIn)
	{
		TInt nb_nodes;
		std::string coord_type;
		AIn >> nb_nodes >> coord_type;
		std::cout << "Nb Vertices: " << nb_nodes << " of type " << coord_type << std::endl;

		for (int i = 0; i < nb_nodes; i++)
		{
			TCoord x, y, z;
			AIn >> x >> y >> z;
			this->mesh_.newNode(x, y, z);
		}
	}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::readNodesData(const std::string& AFileName)
	{
		std::ifstream AIn(AFileName.c_str(), std::ios::in);
		std::string current_word;


		do{
			AIn >> current_word;
		} while (!AIn.eof() && current_word != "POINT_DATA");

		TInt nb_values;
		AIn >> nb_values;
		std::cout << "POINT_DATA: " << nb_values << std::endl;


		//SCALARS VALUE
		std::ifstream scalar_stream(AFileName.c_str(), std::ios::in);
		while (!scalar_stream.eof() && current_word != "CELL_DATA")
		{
			do{
				scalar_stream >> current_word;
			} while (((current_word != "SCALARS") && (current_word != "FIELD")) && !scalar_stream.eof());
			if (current_word == "SCALARS")
			{
				std::string scalar_name, scalar_type, scalar_nb;
				scalar_stream >> scalar_name >> scalar_type >> scalar_nb;
				scalar_stream >> current_word;
				scalar_stream >> current_word;
				std::cout << "READ SCALAR " << scalar_name << " "
					<< scalar_type << " " << scalar_nb << std::endl;
				if (scalar_name != "GMDS_ID"){
					if (scalar_type == "int"){
						Variable<int>* var = this->mesh_.template newVariable<int>(GMDS_NODE, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++){
							int val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
					else if (scalar_type == "float"){
						Variable<double>* var = this->mesh_.template newVariable<double>(GMDS_NODE, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++){
							double val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
				}
			}
			if (current_word == "FIELD")
                        {
                                std::string field_name, field_nb;
                                scalar_stream >> field_name >> field_nb;
                                std::cout<<"READ FIELD "<<field_name<<" "<<field_nb<<std::endl;

				const int max_length = 1000;
	        	        char ch[max_length];
                		int lg;
		                
				std::streampos pos;
                		do{
				        bool curveFound = false;
				        bool vertexFound = false;

					pos = scalar_stream.tellg();
		                        scalar_stream.getline(ch, max_length);
                		        lg = (int)(scalar_stream.gcount());

		                        if(std::string(ch).find("Crb") != std::string::npos) {
                                		curveFound = true;
		                        }
					if(std::string(ch).find("Pt") != std::string::npos) {
                                		vertexFound = true;
		                        }
					

					if(curveFound || vertexFound) {
						scalar_stream.seekg(pos);
	
						std::string cloud_name, cloudComponentType;
						int cloudNumComponents, cloudNumTuples;
						scalar_stream >> cloud_name >> cloudNumComponents >> cloudNumTuples >> cloudComponentType;
						std::cout<<cloud_name << " " << cloudNumComponents << " " << cloudNumTuples << " " <<cloudComponentType <<std::endl;
						
						gmds::IGMesh::cloud& cl = this->mesh_.newCloud(cloud_name);

						for(unsigned int iTuple=0; iTuple<cloudNumTuples; iTuple++) {
						    int val;
						    scalar_stream >> val;
						    if(val == 1) {
						      cl.add(this->mesh_.template get<gmds::Node>(iTuple));
						    }
						}			
					}

		                } while (
		                        lg > 0);

				break;		
                        }
		};

		//VECTORS VALUE
		std::ifstream vector_stream(AFileName.c_str(), std::ios::in);
		do{
			vector_stream >> current_word;
		} while (!vector_stream.eof() && current_word != "POINT_DATA");

		while (!vector_stream.eof() && current_word != "CELL_DATA")
		{
			do{
				vector_stream >> current_word;
			} while (current_word != "VECTORS" && current_word != "CELL_DATA" &&  !vector_stream.eof());
			if (current_word == "VECTORS")
			{
				std::string vector_name, vector_type;
				vector_stream >> vector_name >> vector_type;

				std::cout << "READ VECTOR " << vector_name << " " << vector_type << std::endl;

				Variable<gmds::math::Vector >* var =
					this->mesh_.template newVariable<gmds::math::Vector >(GMDS_NODE, vector_name);

				var->setDomain(nb_values);
				for (int i = 0; i < nb_values; i++){
					double x, y, z;
					vector_stream >> x >> y >> z;
					math::Vector v(x, y, z);
					(*var)[i] = v;
				}

			}
		}
	}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::readCells(std::ifstream& AIn, bool AIsPolydata)
	{
		TInt nb_cells, nb_values;
		AIn >> nb_cells;
		AIn >> nb_values;

		std::cout << "Nb cells: " << nb_cells << std::endl;

		std::vector<std::vector<int> > nodes_of_cell;
		nodes_of_cell.resize(nb_cells);
		for (int i = 0; i < nb_cells; i++)
		{
			std::vector<int> nodes;
			int nb_nodes;
			AIn >> nb_nodes;
			nodes.resize(nb_nodes);
			for (int j = 0; j < nb_nodes; j++)
			{
				int val;
				AIn >> val;
				nodes[j] = val;
			}
			nodes_of_cell[i] = nodes;
		}

		if(AIsPolydata) {
			for (int i = 0; i < nb_cells; i++) {
				std::vector<int> nodes = nodes_of_cell[i];
				switch(nodes.size()) {
				case 3://TRIANGLE
					this->mesh_.newTriangle(nodes[0], nodes[1], nodes[2]);
                                        m_nb_imported_faces++;
					break;
				case 4://QUAD
					this->mesh_.newQuad(nodes[0], nodes[1], nodes[2], nodes[3]);
                                        m_nb_imported_faces++;
					break;
				default:
					throw GMDSException("VTKReader::readCells face type not handled");
					break;
				}
			}
		}
		else {

			std::string current_word;
			do{
				AIn >> current_word;
			} while (current_word != "CELL_TYPES");


			TInt nb_cell_types;
			AIn >> nb_cell_types;
			std::cout << "Nb cell types: " << nb_cell_types << std::endl;
			for (int i = 0; i < nb_cells; i++){
				int cell_type;
				AIn >> cell_type;
				std::vector<int> nodes = nodes_of_cell[i];
				switch (cell_type){
				case 5://TRIANGLE
				{
					   this->mesh_.newTriangle(nodes[0], nodes[1], nodes[2]);
					   m_nb_imported_faces++;
				}
					break;
				case 9://QUAD
				{
					   this->mesh_.newQuad(nodes[0], nodes[1], nodes[2], nodes[3]);
					   m_nb_imported_faces++;
				}
					break;
				case 10://TET
				{
						this->mesh_.newTet(nodes[0], nodes[1], nodes[2], nodes[3]);
						m_nb_imported_regions++;
				}
					break;
				case 12://HEX
				{
						this->mesh_.newHex(nodes[0], nodes[1], nodes[2], nodes[3],
							nodes[4], nodes[5], nodes[6], nodes[7]);
						m_nb_imported_regions++;
				}
					break;
				}
			}
		}
	

	}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::readFacesData(const std::string& AFileName)
	{
		std::ifstream AIn(AFileName.c_str(), std::ios::in);
		std::string current_word;


		do{
			AIn >> current_word;
		} while (!AIn.eof() && current_word != "CELL_DATA");

		TInt nb_values;
		AIn >> nb_values;
		std::cout << "CELL_DATA: " << nb_values << std::endl;


		//SCALARS VALUE
		std::ifstream scalar_stream(AFileName.c_str(), std::ios::in);
		while (!scalar_stream.eof())
		{
			do{
				scalar_stream >> current_word;
			} while (((current_word != "SCALARS") && (current_word != "FIELD")) && !scalar_stream.eof());
			if (current_word == "SCALARS")
			{
				std::string scalar_name, scalar_type, scalar_nb;
				scalar_stream >> scalar_name >> scalar_type >> scalar_nb;
				scalar_stream >> current_word;
				scalar_stream >> current_word;
				std::cout << "READ SCALAR " << scalar_name << " " << scalar_type << " " << scalar_nb << std::endl;
				if (scalar_name != "GMDS_ID")
				{
					if (scalar_type == "int")
					{
						Variable<int>* var = this->mesh_.template newVariable<int>(GMDS_FACE, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++)
						{
							int val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
					else if (scalar_type == "float")
					{
						Variable<double>* var = this->mesh_.template newVariable<double>(GMDS_FACE, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++)
						{
							double val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
				}
			}
                        if (current_word == "FIELD")
                        {
                                std::string field_name, field_nb;
                                scalar_stream >> field_name >> field_nb;
                                std::cout<<"READ FIELD "<<field_name<<" "<<field_nb<<std::endl;

				const int max_length = 1000;
	        	        char ch[max_length];
                		int lg;
		               
				std::streampos pos;
                		do{
				        bool surfaceFound = false;

					pos = scalar_stream.tellg();
		                        scalar_stream.getline(ch, max_length);
                		        lg = (int)(scalar_stream.gcount());

		                        if(std::string(ch).find("Surf") != std::string::npos) {
                                		surfaceFound = true;
		                        }

					if(surfaceFound) {
						scalar_stream.seekg(pos);
	
						std::string surf_name, surfComponentType;
						int surfNumComponents, surfNumTuples;
						scalar_stream >> surf_name >> surfNumComponents >> surfNumTuples >> surfComponentType;
						std::cout<<surf_name << " " << surfNumComponents << " " << surfNumTuples << " " <<surfComponentType <<std::endl;
						gmds::IGMesh::surface& surf = this->mesh_.newSurface(surf_name);

						for(unsigned int iTuple=0; iTuple<surfNumTuples; iTuple++) {
                                                        int val;
                                                        scalar_stream >> val;
							if(val == 1) {
								surf.add(this->mesh_.template get<gmds::Face>(iTuple));
							}
						}			
					}

		                } while (
		                        lg > 0);

				break;		
                        }
		};

		//VECTORS VALUE
		std::ifstream vector_stream(AFileName.c_str(), std::ios::in);

		do{
			vector_stream >> current_word;
		} while (!vector_stream.eof() && current_word != "CELL_DATA");

		while (!vector_stream.eof())
		{
			do{
				vector_stream >> current_word;
			} while (current_word != "VECTORS" && !vector_stream.eof());
			if (!vector_stream.eof())
			{
				std::string vector_name, vector_type;
				vector_stream >> vector_name >> vector_type;

				std::cout << "READ VECTOR " << vector_name << " " << vector_type << std::endl;

				Variable<gmds::math::Vector >* var =
					this->mesh_.template newVariable<gmds::math::Vector >(GMDS_FACE, vector_name);

				var->setDomain(nb_values);
				for (int i = 0; i < nb_values; i++)
				{
					double x, y, z;
					vector_stream >> x >> y >> z;
					math::Vector v(x, y, z);
					(*var)[i] = v;
				}

			}
		}
	}
	/*----------------------------------------------------------------------------*/
	template<typename TMesh>
	void VTKReader<TMesh>::readRegionsData(const std::string& AFileName)
	{
		std::ifstream AIn(AFileName.c_str(), std::ios::in);
		std::string current_word;


		do{
			AIn >> current_word;
		} while (!AIn.eof() && current_word != "CELL_DATA");

		TInt nb_values;
		AIn >> nb_values;
		std::cout << "CELL_DATA: " << nb_values << std::endl;


		//SCALARS VALUE
		std::ifstream scalar_stream(AFileName.c_str(), std::ios::in);
		do{
			scalar_stream >> current_word;
		} while (!scalar_stream.eof() && current_word != "CELL_DATA");

		while (!scalar_stream.eof())
		{
			do{
				scalar_stream >> current_word;
			} while (current_word != "SCALARS" && !scalar_stream.eof());
			if (current_word == "SCALARS")
			{
				std::string scalar_name, scalar_type, scalar_nb;
				scalar_stream >> scalar_name >> scalar_type >> scalar_nb;
				scalar_stream >> current_word;
				scalar_stream >> current_word;
				std::cout << "READ SCALAR " << scalar_name << " " << scalar_type << " " << scalar_nb << std::endl;
				if (scalar_name != "GMDS_ID")
				{
					if (scalar_type == "int")
					{
						Variable<int>* var = this->mesh_.template newVariable<int>(GMDS_REGION, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++)
						{
							int val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
					else if (scalar_type == "float")
					{
						Variable<double>* var = this->mesh_.template newVariable<double>(GMDS_REGION, scalar_name);

						var->setDomain(nb_values);
						for (int i = 0; i < nb_values; i++)
						{
							double val;
							scalar_stream >> val;
							(*var)[i] = val;
						}
					}
				}
			}
		};

		//VECTORS VALUE
		std::ifstream vector_stream(AFileName.c_str(), std::ios::in);
		do{
			vector_stream >> current_word;
		} while (!vector_stream.eof() && current_word != "CELL_DATA");

		while (!vector_stream.eof())
		{
			do{
				vector_stream >> current_word;
			} while (current_word != "VECTORS" && !vector_stream.eof());
			if (!vector_stream.eof())
			{
				std::string vector_name, vector_type;
				vector_stream >> vector_name >> vector_type;

				std::cout << "READ VECTOR " << vector_name << " " << vector_type << std::endl;

				Variable<gmds::math::Vector >* var =
					this->mesh_.template newVariable<gmds::math::Vector >(GMDS_REGION, vector_name);

				var->setDomain(nb_values);
				for (int i = 0; i < nb_values; i++)
				{
					double x, y, z;
					vector_stream >> x >> y >> z;
					math::Vector v(x, y, z);
					(*var)[i] = v;
				}

			}
		}
	}
	/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_VTKREADER_H_ */
/*----------------------------------------------------------------------------*/
