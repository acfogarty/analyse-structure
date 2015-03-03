#include "analyseStructure.h"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>

main() 
{
   System system;
   AnalysisInfo analysisInfo;
   
   system.getInputData();
   system.initializeSystem();
   analysisInfo.getInputData();
 
   std::vector<double> avevalues,avevalues_sq,stdev;
   for (int i=0;i<analysisInfo.nmeasurements;i++)
   {
      avevalues.push_back(0.0);
      avevalues_sq.push_back(0.0);
      stdev.push_back(0.0);
   }
 
   int ncount=0;
 
   std::vector<std::ofstream*> outputFileStreams;
   outputFileStreams.resize(analysisInfo.nmeasurements);
   for (int j=0;j<analysisInfo.nmeasurements;j++)
   {
      outputFileStreams[j] = new std::ofstream(analysisInfo.measurementEntries[j].filename);
   }

   if (system.coordFileType==System::filetypeEnum::DCD) {system.readDcdHeader();}
 
   for (int i=0;i<system.nframes;i++)
   {
      switch(system.coordFileType)
      {
         case(System::filetypeEnum::GRO):
            system.readGroFrame(false);
            break;
         case(System::filetypeEnum::PDB):
            system.readPdbFrame();
            break;
         case(System::filetypeEnum::DCD):
            system.readDcdFrame();
            break;
         default:
            std::cout<<"Error in main!"<<std::endl;
            //TODO
      }
      if ((i%system.nskip) != 0) continue;
      ncount+=1;
      for (int j=0;j<analysisInfo.nmeasurements;j++)
      {
         *outputFileStreams[j] << i+1 << " ";
         //calculation the distances, angles and dihedrals specified in analysisInfo 
         double value;
         std::stringstream descriptionStream;
         switch(analysisInfo.measurementEntries[j].mtype)
         {
            case(AnalysisInfo::measurementtypeEnum::BON):
               value=Calc::calcDist(system.atomEntries, system.cell, analysisInfo.measurementEntries[j].atomindices);
               break;
            case(AnalysisInfo::measurementtypeEnum::ANG):
               value=Calc::calcAngle(system.atomEntries, system.cell, analysisInfo.measurementEntries[j].atomindices);
               break;
            case(AnalysisInfo::measurementtypeEnum::DIH):
               value=Calc::calcDihedral(system.atomEntries, system.cell, analysisInfo.measurementEntries[j].atomindices);
               break;
            default:
               std::cout<<"Error in main!"<<std::endl;
               //TODO
         }
         *outputFileStreams[j] << value << " ";
         avevalues[j]+=value;
         avevalues_sq[j]+=(value*value);
         *outputFileStreams[j] << "\n";
      } //loop over nmeasurements
   } //loop over nframes
   
   for (int j=0;j<analysisInfo.nmeasurements;j++)
   {
      outputFileStreams[j]->close();
      //get average and standard deviation
      avevalues[j]/=(double)ncount;
      avevalues_sq[j]/=(double)ncount;
      stdev[j]=sqrt(avevalues_sq[j]-avevalues[j]*avevalues[j]);
      //output
      std::stringstream descriptionStream;
      switch(analysisInfo.measurementEntries[j].mtype)
      {
         case(AnalysisInfo::measurementtypeEnum::BON):
            descriptionStream << "Bond     " << analysisInfo.measurementEntries[j].atomindices[0] << " " << analysisInfo.measurementEntries[j].atomindices[1] << "     ";
            break;
         case(AnalysisInfo::measurementtypeEnum::ANG):
            descriptionStream << "Angle    " << analysisInfo.measurementEntries[j].atomindices[0] << " " << analysisInfo.measurementEntries[j].atomindices[1] << " " << analysisInfo.measurementEntries[j].atomindices[2] << "   ";
            break;
         case(AnalysisInfo::measurementtypeEnum::DIH):
            descriptionStream << "Dihedral " << analysisInfo.measurementEntries[j].atomindices[0] << " " << analysisInfo.measurementEntries[j].atomindices[1] << " " << analysisInfo.measurementEntries[j].atomindices[2] << " " << analysisInfo.measurementEntries[j].atomindices[3] << " ";
            break;
      }
   std::cout<<"# "<<descriptionStream.str()<<" "<<avevalues[j]<<" stdev "<<stdev[j]<<std::endl;
   } //loop over nmeasurements
}
