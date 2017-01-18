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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1UserEventInformation::G1UserEventInformation() : fPhotonCount_Scint(0), fPhotonCount_Ceren(0) , fEdep(0)
  {
  
  fProcessMap = new std::map < G4int, G4int>();
  
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
}


void G1UserEventInformation::PrintProcessNumber() 
{


for( std::map<const G4int, G4int>::iterator ii=fProcessMap->begin(); ii!=fProcessMap->end(); ++ii)

   {
       //küçükten büyüğe doğru sıralıyor std::map(index numarasına göre)
       G4cout << (*ii).first << ": " << (*ii).second << G4endl;

   }



}


