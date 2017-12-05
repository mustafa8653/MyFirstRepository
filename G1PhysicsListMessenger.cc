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
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G1PhysicsListMessenger.hh"

#include "G1PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1PhysicsListMessenger::
  G1PhysicsListMessenger(G1PhysicsList* pPhys) 
  : G4UImessenger(), fPhysicsList(pPhys)
{

//G4UIdirectory (char *theCommandPath, G4bool commandsToBeBroadcasted=true)
  fG1Dir = new G4UIdirectory("/G1/");
  fG1Dir->SetGuidance("UI commands of this example");

  fPhysDir = new G4UIdirectory("/G1/phys/");
  fPhysDir->SetGuidance("PhysicsList control");
 
 
 //burda omittable true parametre girmese defaul deger var
 //G4UIcmdWithAnInteger (const char *theCommandPath, G4UImessenger *theMessenger)
  fVerboseCmd = new G4UIcmdWithAnInteger("/G1/phys/verbose",this);
  fVerboseCmd->SetGuidance("set verbose for physics processes");
  
 //void SetParameterName (const char *theName, G4bool omittable, G4bool currentAsDefault=false)
  fVerboseCmd->SetParameterName("verbose",true);
  
  //void SetDefaultValue (G4int defVal)
  fVerboseCmd->SetDefaultValue(1);
  fVerboseCmd->SetRange("verbose>=0");
  fVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
 
 
 //burda omittable false parametreye gerek var
  fCerenkovCmd =
           new G4UIcmdWithAnInteger("/G1/phys/cerenkovMaxPhotons",this);
  fCerenkovCmd->SetGuidance("set max nb of photons per step");
  fCerenkovCmd->SetParameterName("MaxNumber",false);
  fCerenkovCmd->SetRange("MaxNumber>=0");
  fCerenkovCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1PhysicsListMessenger::~G1PhysicsListMessenger()
{
  delete fVerboseCmd;
  delete fCerenkovCmd;
  delete fPhysDir;
  delete fG1Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  if( command == fVerboseCmd )
   {fPhysicsList->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));}

  if( command == fCerenkovCmd )
   {fPhysicsList->
              SetNbOfPhotonsCerenkov(fCerenkovCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
