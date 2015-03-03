#include "analyseStructure.h"
#include "Vec3D.h"
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>


void AnalysisInfo::getInputData()
{
   std::string line;
   std::getline (std::cin,line); //skip comments
   std::cin>>nmeasurements;
   std::cin.ignore(10000,'\n');
   std::getline (std::cin,line); //skip comments
   for (int i=0;i<nmeasurements;i++)
   {
      MeasurementEntry tempEntry;
      std::string line,temp;
      int itemp;
   
      //read by line 
      std::getline(std::cin, line);
      std::istringstream iss(line);
 
      //get filename
      iss>>tempEntry.filename;

      //get type (BON,ANG or DIH)
      iss>>temp;
      tempEntry.mtype = std::find(mtypestrings.begin(), mtypestrings.end(), temp) - mtypestrings.begin();

      //get atom indices
      while (iss >> itemp) 
      {    
         tempEntry.atomindices.push_back(itemp);
      }

      if (tempEntry.atomindices.size() != (tempEntry.mtype+2))
      {
         std::cout<<"Error! "<<temp<<" entries should have "<<tempEntry.mtype+2<<" atom indices"<<std::endl;
      }

      measurementEntries.push_back(tempEntry);
   }
}

void System::getInputData()
{
   std::string line,filetypeString;
   char coordFileName[50];
   std::getline (std::cin,line); //skip comments
   std::getline (std::cin,line); //skip comments
   std::cin>>filetypeString;
   if (filetypeString=="PDB") {std::cout<<"Warning, not using PBC for pdb filetype. (Assuming infinite cell)."<<std::endl;}
   coordFileType = std::find(filetypestrings.begin(), filetypestrings.end(), filetypeString) - filetypestrings.begin();
   std::cin.ignore(10000,'\n');
   std::cin>>coordFileName;
   if (filetypeString=="DCD")
   {
      inputCoordStream.open(coordFileName,std::ifstream::in | std::ios::binary);
   }
   else
   {
      inputCoordStream.open(coordFileName,std::ifstream::in);
   }
   std::cin.ignore(10000,'\n');
   std::cin>>natoms;
   std::cout<<"natoms"<<natoms<<std::endl;
   std::cin.ignore(10000,'\n');
   std::cin>>nframes;
   std::cout<<"nframes"<<nframes<<std::endl;
   std::cin.ignore(10000,'\n');
   std::cin>>nskip;
   std::cin.ignore(10000,'\n');
}

void System::initializeSystem()
{
   for (int i=0;i<natoms;i++)
   {
      atomEntries.push_back(emptyAtomEntry); 
   }
}

void System::readGroFrame(bool velocitiesPresent)
{
   std::string line;
   std::getline(inputCoordStream, line); //skip comments
   int natomstemp;
   inputCoordStream >> natomstemp;
   if (natomstemp!=natoms) {
      std::cout<<"Error! Natoms in gro file not equal to natoms in input file"<<std::endl;
   }
   inputCoordStream.ignore(10000,'\n');
   for (int i=0;i<natoms;i++)
   {
      std::getline(inputCoordStream, line);
      atomEntries.at(i).coords[0]=stof(line.substr(20,8));
      atomEntries.at(i).coords[1]=stof(line.substr(28,8));
      atomEntries.at(i).coords[2]=stof(line.substr(36,8));
      if (velocitiesPresent) 
      { 
        //TODO
      }
   }
   inputCoordStream >> cell[0] >> cell[1] >> cell[2];
   inputCoordStream.ignore(10000,'\n');
}

void System::readPdbFrame()
{
   std::string line;
   for (int i=0;i<natoms;i++)
   {
      std::getline(inputCoordStream, line);
      atomEntries.at(i).coords[0]=stof(line.substr(30,8));
      atomEntries.at(i).coords[1]=stof(line.substr(38,8));
      atomEntries.at(i).coords[2]=stof(line.substr(46,8));
   }
   cell[0]=10000.0;
   cell[1]=10000.0;
   cell[2]=10000.0;
}

void System::readDcdHeader()
{
   //format as given in http://code.zgib.net/mp3/doc/dcd.txt
   char onebyte;
   int idum,ntitle;
   int sizeDcdInt = 4;
   if (sizeof(int) != sizeDcdInt) {std::cout<<"Error size of int here != size of int in DCD file (4 bytes)"<<std::endl;}
   
   //first header section
   int headerByteSize = 92;
   int nIntegers = headerByteSize/sizeDcdInt;
   for(int i=0;i<nIntegers;i++)
   {
      inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum));
      //std::cout<<idum<<std::endl;
   }
   //std::cout<<"end of header, part 1"<<std::endl;
  
   //second header section
   inputCoordStream.read(reinterpret_cast<char*>(&ntitle),sizeof(ntitle));
   int lengthTitleLine=80;
   int nlines=(ntitle-sizeDcdInt)/lengthTitleLine;
   nIntegers = lengthTitleLine/sizeDcdInt;
   for(int i=0;i<nlines;i++)
   {
      for(int j=0;j<nIntegers;j++)
      {
         inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum)); 
      }
   }
   inputCoordStream.read(reinterpret_cast<char*>(&ntitle),sizeof(ntitle));
   //std::cout<<"end of header, part 2"<<std::endl;

   //third header section
   int nEntries = 4;
   for(int i=0;i<nEntries;i++)
   {
      inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum));
      //std::cout<<idum<<std::endl;
   }
   //std::cout<<"end of header, part 3"<<std::endl;
}

void System::readDcdFrame()
{
   double tempdble;
   float tempflt;
   int idum; //assumes size is 4 bytes, already checked in readDcdHeader
   double table[6];
   inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum)); 
std::cout<<"start cell"<<std::endl;
   for (int j=0;j<6;j++)
   { 
      inputCoordStream.read(reinterpret_cast<char*>(&tempdble), sizeof tempdble);
      std::cout<<tempdble<<std::endl;
      table[j]=tempdble;
   }
std::cout<<"end cell"<<std::endl;
   cell[0]=table[0];
   cell[1]=table[2];
   cell[2]=table[5];
   inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum)); 
   for (int j=0;j<3;j++)
   {
      inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum)); //should be natoms*4
      for (int i=0;i<natoms;i++)
      {
         //inputCoordStream.read(reinterpret_cast<char*>(&tempdble), sizeof tempdble);
         inputCoordStream.read(reinterpret_cast<char*>(&tempflt), sizeof tempflt);
         //std::cout<<tempflt<<std::endl;
         atomEntries.at(i).coords[j]=tempflt; 
      }
      inputCoordStream.read(reinterpret_cast<char*>(&idum),sizeof(idum)); //should be natoms*4
   }
}

double Calc::calcDist(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{
   double dist,dx,dy,dz;
   double x1,y1,z1,x2,y2,z2;

   if (indices.size()!=2) {std::cout<<"Error in calcDist!"<<std::endl;}
   int index1=indices.at(0)-1;
   int index2=indices.at(1)-1;
   x1=atomEntries.at(index1).coords[0];
   y1=atomEntries.at(index1).coords[1];
   z1=atomEntries.at(index1).coords[2];
   x2=atomEntries.at(index2).coords[0];
   y2=atomEntries.at(index2).coords[1];
   z2=atomEntries.at(index2).coords[2];

   dx = x1-x2;
   dy = y1-y2;
   dz = z1-z2;
   dx = dx - rint(dx/cell[0])*cell[0];
   dy = dy - rint(dy/cell[1])*cell[1];
   dz = dz - rint(dz/cell[2])*cell[2];
   dist=sqrtf(dx*dx + dy*dy + dz*dz);

   return dist;
}

double Calc::calcAngle(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{

   //calculate angle index1-index2-index3 (with index2 at corner of angle)

   double theta,dx,dy,dz;
   Vec3D vec21,vec23;
   
   if (indices.size()!=3) {std::cout<<"Error in calcAngle!"<<std::endl;}
   int index1=indices.at(0)-1;
   int index2=indices.at(1)-1;
   int index3=indices.at(2)-1;

   dx = atomEntries.at(index1).coords[0]-atomEntries.at(index2).coords[0];
   vec21[0] = dx - rint(dx/cell[0])*cell[0];
   dy = atomEntries.at(index1).coords[1]-atomEntries.at(index2).coords[1];
   vec21[1] = dy - rint(dy/cell[1])*cell[1];
   dz = atomEntries.at(index1).coords[2]-atomEntries.at(index2).coords[2];
   vec21[2] = dz - rint(dz/cell[2])*cell[2];

   dx = atomEntries.at(index3).coords[0]-atomEntries.at(index2).coords[0];
   vec23[0] = dx - rint(dx/cell[0])*cell[0];
   dy = atomEntries.at(index3).coords[1]-atomEntries.at(index2).coords[1];
   vec23[1] = dy - rint(dy/cell[1])*cell[1];
   dz = atomEntries.at(index3).coords[2]-atomEntries.at(index2).coords[2];
   vec23[2] = dz - rint(dz/cell[2])*cell[2];

   theta=acos(vec21*vec23/(vec21.len()*vec23.len()));
   //theta=theta/pi*180.0;
   theta=theta/3.14159*180.0;
   return theta;
}

double Calc::calcDihedral(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{

   //calculate dihedral index1-index2-index3-index4 

   double phi,dx,dy,dz,b;
   Vec3D vec12,vec23,vec34;
   Vec3D vectemp1,vectemp2,vectemp3;

   if (indices.size()!=4) {std::cout<<"Error in calcDihedral!"<<std::endl;}
   int index1=indices.at(0)-1;
   int index2=indices.at(1)-1;
   int index3=indices.at(2)-1;
   int index4=indices.at(3)-1;

   dx = atomEntries.at(index2).coords[0]-atomEntries.at(index1).coords[0];
   vec12[0] = dx - rint(dx/cell[0])*cell[0];
   dy = atomEntries.at(index2).coords[1]-atomEntries.at(index1).coords[1];
   vec12[1] = dy - rint(dy/cell[1])*cell[1];
   dz = atomEntries.at(index2).coords[2]-atomEntries.at(index1).coords[2];
   vec12[2] = dz - rint(dz/cell[2])*cell[2];

   dx = atomEntries.at(index3).coords[0]-atomEntries.at(index2).coords[0];
   vec23[0] = dx - rint(dx/cell[0])*cell[0];
   dy = atomEntries.at(index3).coords[1]-atomEntries.at(index2).coords[1];
   vec23[1] = dy - rint(dy/cell[1])*cell[1];
   dz = atomEntries.at(index3).coords[2]-atomEntries.at(index2).coords[2];
   vec23[2] = dz - rint(dz/cell[2])*cell[2];
 
   dx = atomEntries.at(index4).coords[0]-atomEntries.at(index3).coords[0];
   vec34[0] = dx - rint(dx/cell[0])*cell[0];
   dy = atomEntries.at(index4).coords[1]-atomEntries.at(index3).coords[1];
   vec34[1] = dy - rint(dy/cell[1])*cell[1];
   dz = atomEntries.at(index4).coords[2]-atomEntries.at(index3).coords[2];
   vec34[2] = dz - rint(dz/cell[2])*cell[2];
    
   vectemp1 = vec12.cross(vec23);
   vectemp2 = vec23.cross(vec34);

   b = vectemp1*vectemp2/(vectemp1.len()*vectemp2.len());

   phi=acos(b);

   //get sign
   vectemp3 = vectemp1.cross(vectemp2);
   double signcheck = vectemp3*vec23;
   if (signcheck > 0.0) signcheck = 1.0;
   else if (signcheck < 0.0) signcheck = -1.0;
   else signcheck = 0.0;

   phi=phi*signcheck;
   phi=phi/3.14159*180.0;
   return phi;

}

