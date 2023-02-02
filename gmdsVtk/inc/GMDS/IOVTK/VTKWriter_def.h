/*----------------------------------------------------------------------------*/
/** \file    VTKWriter_def.h
 *  \author  F. LEDOUX
 *  \date    02/24/2009
 */
/*----------------------------------------------------------------------------*/
template<typename TMesh>
class VTKWriter{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param AMesh the mesh we want to write into a file.
     */
	VTKWriter(TMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~VTKWriter();

    /*------------------------------------------------------------------------*/
    /** \brief  Write the content of mesh_ into the file named AFileName.
     *
     *  \param AFileName name of the output file without its extension
     *  \param AMask	 we check it to find R and/or F to know which type of
     *  				 vtk file to create.
     */
	void write(const std::string& AFileName, const int& AMask);

protected:

    /*------------------------------------------------------------------------*/
    /** \brief  Method creating a VTK unstructured grid and saving in .vtu file.
     */
	void createUnstructuredGrid(const std::string& AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Method creating a VTK poly data and saving in .vtp file.
     */
	void createPolyData(const std::string& AFileName);

	/* a mesh */
	TMesh& mesh_;

	/** mesh dimension */
	int mesh_dimension_;

};
/*----------------------------------------------------------------------------*/
