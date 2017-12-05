//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: LXeUserEventInformation.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/src/LXeUserEventInformation.cc
/// \brief Implementation of the LXeUserEventInformation class
//
//
#include "G1UserEventInformation.hh"
#include <iomanip>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1UserEventInformation::G1UserEventInformation( G4int eventID) : fPhotonCount_Scint(0), fPhotonCount_Ceren(0) , fAbsorbedPhotonCount(0),  fPmt1DetCount(0), fPmt2DetCount(0), fPmt1AbsCount(0), fPmt2AbsCount(0), fEdep(0.)
  {
  
  fEventID = eventID;
  
  fProcessMap = new std::map < G4int, G4int> ();
  fPhotonEnergyVec = new std::vector <G4double> ();
  fDetectedPmtPhotonEnergyVec = new std::vector <G4double> ();
  fScinOrCherVec = new std::vector <G4int> ();
  
  (*fProcessMap)[0]=0; //Undefinied
  (*fProcessMap)[1]=0; //Transmission
  (*fProcessMap)[2]=0; //FresnelRefraction
  (*fProcessMap)[3]=0; //FresnelReflection
  (*fProcessMap)[4]=0; //TotalInternalReflection
  (*fProcessMap)[5]=0; //LambertianReflection
  (*fProcessMap)[6]=0; //LobeReflection
  (*fProcessMap)[7]=0; //SpikeReflection
  (*fProcessMap)[8]=0; //BackScattering
  (*fProcessMap)[9]=0; //Absorption
  (*fProcessMap)[10]=0; //Detection
  (*fProcessMap)[13]=0; //StepTooSmall
  
  
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1UserEventInformation::~G1UserEventInformation() {

delete fProcessMap;
delete fPhotonEnergyVec;
delete fDetectedPmtPhotonEnergyVec;
delete fScinOrCherVec;
}



void G1UserEventInformation::PrintEventResult()
{


G4cout<<"-------------------------------------------------"<<G4endl;

G4cout<<" \n Event "<<fEventID<<" Statistic: "<<G4endl;

G4cout<<"Deposited energy in this event: "<< std::setprecision(3) <<G4BestUnit(fEdep ,"Energy") <<G4endl;

G4cout<<"# of emitted scintillation photons: "<< fPhotonCount_Scint <<G4endl;

G4cout<<"# of emitted cherenkov photons : "<< fPhotonCount_Ceren <<G4endl;

G4cout<< " # of total emitted (scin + cheren ) photons :  "<< fPhotonCount_Ceren + fPhotonCount_Scint <<G4endl;

G4cout<<" # of detected photons from pmt1 surface :  "<< GetPmt1DetCount() <<G4endl;
G4cout<<" # of detected photons from pmt2 surface :  "<< GetPmt2DetCount() <<G4endl;
G4cout<<" # of total detected photons from pmt surface :  "<< GetTotalDetectedPhotonNumber() <<G4endl; 

G4cout<<" # of absorbed photons from pmt1(Negatif Side) surface :  "<< GetPmt1AbsCount() <<G4endl;
G4cout<<" # of absorbed photons from pmt2(Pozitif side) surface :  "<< GetPmt2AbsCount() <<G4endl;

G4cout<<" # of total absorbed photons from pmt surface :  "<< GetPmt1AbsCount() + GetPmt2AbsCount() <<G4endl;

G4cout<<" # of absorbed photon in all volume(bulk abs):  "<< GetAbsorbedPhotonCount() <<G4endl;
    
G4cout<<"---OpticalBoundaryProcess statistic---"<<G4endl;

for( std::map<const G4int, G4int>::iterator ii=fProcessMap->begin(); ii!=fProcessMap->end(); ++ii)

   	{
       //küçükten büyüğe doğru sıralıyor std::map(index numarasına göre)
       G4cout << (*ii).first << ": " << (*ii).second << G4endl;

   	}
 
 /*  
G4cout<<"------Photon PrimaryVertex Energy-------"<<G4endl;
G4cout<<"# of total photon in this event: "<<fPhotonEnergyVec->size()<<G4endl;

for (unsigned int i=0; i< fPhotonEnergyVec->size(); i++) //fPhotonEnergyVec->size() unsigned olduğu için mecbure unsigned olmalı i
	{
         G4cout<<i<<": "<<G4BestUnit((*fPhotonEnergyVec)[i],"Energy")<<G4endl;
	}   

G4cout<<"-------------------------------------------------"<<G4endl;

  */



/*
G4cout<<"------Pmt Detected Photon Energy-------"<<G4endl;
G4cout<<"# of detected photon from pmt surface "<<fDetectedPmtPhotonEnergyVec->size()<<G4endl;

for (unsigned int i=0; i< fDetectedPmtPhotonEnergyVec->size(); i++) //fPhotonEnergyVec->size() unsigned olduğu için mecbure unsigned olmalı i
	{
         G4cout<<i<<": "<<G4BestUnit((*fDetectedPmtPhotonEnergyVec)[i],"Energy")<<G4endl;
	}   

G4cout<<"-------------------------------------------------"<<G4endl;

  */

}





