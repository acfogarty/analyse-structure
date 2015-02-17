#include "analyseStructure.h"
#include "Vec3D.hpp"
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>


void AnalysisInfo::getInputData()
{
    std::string line;
    std::string temp;
    std::getline (std::cin,line); //skip comments
    std::cin>>nmeasurements;
    std::cin.ignore(10000,'\n');
    std::getline (std::cin,line); //skip comments
    for (int i=0;i<nmeasurements;i++)
    {
       MeasurementEntry tempEntry;
       std::string line,enumstring;
       int itemp;
       AnalysisInfo::measurementtypes itemp2;
     
       std::getline(std::cin, line);
       std::istringstream iss(line);
 
       //get filename
       iss>>tempEntry.filename;

       //get type (BON,ANG or DIH)
       iss>>temp;
       tempEntry.mtype = std::find(mtypestrings.begin(), mtypestrings.end(), temp) - mtypestrings.begin();

       //get atom indices
       while ( iss >> itemp) 
       {    
         tempEntry.atomindices.push_back(itemp);
       }

       if (tempEntry.atomindices.size() != (tempEntry.mtype+2))
       {
         std::cout<<"Error! "<<temp<<" entries should have "<<tempEntry.mtype+2<<" atom indices"<<std::endl;
       }

       measurementEntries.push_back(tempEntry);
    }
};

void System::getInputData()
{
    std::string line;
    std::getline (std::cin,line); //skip comments
    std::getline (std::cin,line); //skip comments
    char grofilename[50];
    std::cin>>grofilename;
    inputGroStream.open(grofilename,std::ifstream::in);
    std::cin.ignore(10000,'\n');
    std::cin>>natoms;
    std::cout<<"natoms"<<natoms<<std::endl;
    std::cin.ignore(10000,'\n');
    std::cin>>nframes;
    std::cout<<"nframes"<<nframes<<std::endl;
    std::cin.ignore(10000,'\n');
    std::cin>>nskip;
    std::cin.ignore(10000,'\n');
};

void System::initializeSystem()
{
    for (int i=0;i<natoms;i++)
    {
      atomEntries.push_back(emptyAtomEntry); 
    }
};

void System::readGroFrame(bool velocitiesPresent)
{
   std::string line;
   std::getline(inputGroStream, line); //skip comments
   int natomstemp;
   inputGroStream >> natomstemp;
   if (natomstemp!=natoms) {
     std::cout<<"Error! Natoms in gro file not equal to natoms in input file"<<std::endl;
   }
   inputGroStream.ignore(10000,'\n');
   for (int i=0;i<natoms;i++)
   {
      std::getline(inputGroStream, line);
      atomEntries.at(i).coords[0]=stof(line.substr(20,8));
      atomEntries.at(i).coords[1]=stof(line.substr(28,8));
      atomEntries.at(i).coords[2]=stof(line.substr(36,8));
      if (velocitiesPresent) 
      {
      }
   }
   inputGroStream >> cell[0] >> cell[1] >> cell[2];
   inputGroStream.ignore(10000,'\n');
};

double Calc::calcDist(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{
   double dist,dx,dy,dz;
   double x1,y1,z1,x2,y2,z2;

   if (indices.size()!=2) {std::cout<<"Error in calcDist!"<<std::endl;}
   int index1=indices.at(0);
   int index2=indices.at(1);
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
};

double Calc::calcAngle(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{

   //calculate angle index1-index2-index3 (with index2 at corner of angle)

   double theta,dx,dy,dz;
   Vec3D vec21,vec23;
   
   if (indices.size()!=3) {std::cout<<"Error in calcAngle!"<<std::endl;}
   int index1=indices.at(0);
   int index2=indices.at(1);
   int index3=indices.at(2);

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
};

double Calc::calcDihedral(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices)
{

   //calculate dihedral index1-index2-index3-index4 

   double phi,dx,dy,dz,b;
   Vec3D vec12,vec23,vec34;
   Vec3D vectemp1,vectemp2;

   if (indices.size()!=4) {std::cout<<"Error in calcDihedral!"<<std::endl;}
   int index1=indices.at(0);
   int index2=indices.at(1);
   int index3=indices.at(2);
   int index4=indices.at(3);

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
   phi=phi/3.14159*180.0;
   return phi;

};

