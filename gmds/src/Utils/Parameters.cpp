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
 *  Parameters.cpp
 *
 *  Created on: April 12, 2016
 *  Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/Parameters.h>
/*----------------------------------------------------------------------------*/
#include <limits>
/*----------------------------------------------------------------------------*/
#include <dictionary.h>
#include <iniparser.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
Parameters::Parameters()
:m_initialized(false)
{}
/*----------------------------------------------------------------------------*/
Parameters::~Parameters()
{}
/*----------------------------------------------------------------------------*/
bool Parameters::add(const std::string& ASection,
                     const std::string& AName,
                     const ETypeParam AType)
{
    for(int e=0; e<m_entries.size(); e++){
        if(m_entries[e].section==ASection && m_entries[e].name==AName){
            return false;
        }
    }
    
    Entry e;
    e.section = ASection;
    e.name    = AName;
    e.type    = AType;
    
    m_entries.push_back(e);
    return true;
}
/*----------------------------------------------------------------------------*/
std::vector<Parameters::Entry> Parameters::getEntries() const
{
    return m_entries;
}
/*----------------------------------------------------------------------------*/
int Parameters::
find(const std::string& ASection, const std::string& AName)
{
    for(auto i=0; i<m_entries.size();i++){
        Entry e = m_entries[i];
        if(e.section==ASection && e.name==AName){
            //we found the right entry
            return i;
        }
    }
    return -1;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, double& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=DOUBLE_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=std::stod(m_values[i]);
    return true;
}/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, int& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=INT_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=std::stoi(m_values[i]);
    return true;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, std::string& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=STRING_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=m_values[i];
    return true;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, bool& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=BOOL_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    std::string str_val = m_values[i];
    if(str_val=="true"){
        AOut=true;
    }
    else if(str_val=="false"){
        AOut=false;
    }
    else{
        //Wrong value
        return false;
    }
    return true;
}
/*----------------------------------------------------------------------------*/
std::vector<std::string> Parameters::parseIni(const std::string& AFileName)
{
    std::vector<std::string> wrong;
    if(m_entries.empty()){
        //nothing to do
        return wrong;
    }
    m_values.resize(m_entries.size());
    
    //We create a dictionnary object to parse the read file
    dictionary* dico = iniparser_load(AFileName.c_str());
    
    //For each entry, we look a right type value in the file
    for(auto i=0; i<m_entries.size();i++){
        Entry ei = m_entries[i];
        std::string param_s = ei.section+":"+ei.name;
        const char* param_ch = param_s.c_str();
        // ----------- INT VALUE -----------------------------
        if(ei.type==INT_P){
	    int not_found = std::numeric_limits<int>::min();
            int v = iniparser_getint(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else{
                m_values[i]=std::to_string(v);
            }
        }
        // ----------- DOUBLE VALUE -----------------------------
        else if(ei.type==DOUBLE_P){
            double not_found = std::numeric_limits<double>::min();
            double v = iniparser_getdouble(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else{
                m_values[i]=std::to_string(v);
            }
        }  // ----------- BOOL VALUE -----------------------------
        else if(ei.type==BOOL_P){
            int not_found = 10;
            int v = iniparser_getboolean(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else if (v==1){
                m_values[i]="true";
            }
            else {
                m_values[i]="false";
            }
        }
        // ----------- STRING VALUE -----------------------------
        else {
            char* v = iniparser_getstr(dico,(char *)param_ch);
            if(v==NULL){
                wrong.push_back(param_s);
            }
            else{
                std::string v_s(v);
                m_values[i]=v_s;
            }
        }
    }
    
    //All the value have been found and initialized
    m_initialized=true;
    //We free the dictionnary
    iniparser_freedict(dico);
    
    return wrong;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
