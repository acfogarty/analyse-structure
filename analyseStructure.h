#ifndef ANALYSESTRUCTURE_H
#define ANALYSESTRUCTURE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

struct AtomEntry
{
   double coords[3];
   double velocities[3];
};

struct MeasurementEntry
{
   std::string filename;
   int mtype; //ENUM
   std::vector<int> atomindices;
};

class Calc 
{
   public:
      static double calcDist(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices);
      //static double calcDist(std::vector<AtomEntry> atomEntries, std::vector<int> indices);
      static double calcAngle(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices);
      static double calcDihedral(std::vector<AtomEntry> atomEntries, double cell[3], std::vector<int> indices);

      const double pi = 3.14159;
};

class System 
{
   public:
      //System();
      //~System();
      void getInputData();
      void initializeSystem();
      void readGroFrame(bool velocitiesPresent);
      void readPdbFrame();
      void readDcdFrame();
      void readDcdHeader();

      std::vector<AtomEntry> atomEntries;
      double cell[3];
      int natoms,nframes,nskip;
      std::ifstream inputCoordStream;
      enum filetypeEnum {GRO, PDB, DCD};
      int coordFileType;
      std::vector<std::string> filetypestrings = {"GRO", "PDB", "DCD"};
      
   private:
      AtomEntry emptyAtomEntry; //should initialise with zeroes TODO
};

class AnalysisInfo
{
   public:
      int nmeasurements;
      void getInputData();
      std::vector<MeasurementEntry> measurementEntries;
      enum measurementtypeEnum {BON, ANG, DIH };
      std::vector<std::string> mtypestrings = {"BON", "ANG", "DIH"};
   private:
};


#endif
