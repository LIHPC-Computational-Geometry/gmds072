/*----------------------------------------------------------------------------*/
/** \file    VTKReader_def.h
 *  \author  F. LEDOUX
 *  \date    03/16/2009
 */
/*----------------------------------------------------------------------------*/
template<typename TMesh>
class VTKReader:public IReader<TMesh>
{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param AMesh the mesh we want to write into a file.
     */
	VTKReader(TMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~VTKReader();

    /*------------------------------------------------------------------------*/
    /** \brief  Read the content of the file named outputName_ and write it in
     *   		mesh_.
     */
	void read(const std::string& AFileName);

protected:

    /*------------------------------------------------------------------------*/
    /** \brief  Create a mesh  with the content of the reader AReader.
     *
     *  \param AReader a vtk reader giving access to the content of an XML file
     *  			   storing a vtk unstructured grid
     */
	void read(vtkXMLUnstructuredGridReader* AReader);

    /*------------------------------------------------------------------------*/
    /** \brief  Create a mesh  with the content of the reader AReader.
     *
     *  \param AReader a vtk reader giving access to the content of an XML file
     *  			   storing a vtk poly data
     */
	void read(vtkXMLPolyDataReader* AReader);

};
/*----------------------------------------------------------------------------*/
