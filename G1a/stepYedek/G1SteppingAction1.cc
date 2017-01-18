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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "G1SteppingAction.hh"


#include "G4Step.hh"
//#include "G4Event.hh"
//#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"

#include "G1EventAction.hh"
#include "G1DetectorConstruction.hh"

#include "G4RunManager.hh"

#include "G4TrackStatus.hh"  //fStopAndKill bunun içinde



#include "G4ParticleDefinition.hh" 
#include "G4ParticleTypes.hh"  //G4Gamma için

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::G1SteppingAction(G1EventAction* eventAction)
: G4UserSteppingAction(),fEventAction(eventAction),bc_400_logic(0)
  
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::~G1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1SteppingAction::UserSteppingAction(const G4Step* step)
{

/*
 static G4double totalEnergyDepositForAllStep=0;
  G4double energyDeposit=0;
  
  energyDeposit=step->GetTotalEnergyDeposit();
  
  if( step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()=="BC-400" ) 
  {
 
  G4cout<<step->GetTrack()->GetCurrentStepNumber()<<" inci stepde ,"<<step->GetPreStepPoint()->GetMaterial()->GetName()<<" materyeli içersinde depolanan enerji "<<energyDeposit<<G4endl;
  
  totalEnergyDepositForAllStep+=energyDeposit;
  
  G4cout<<"Bütüp stepler toplamında depolanan enerji: "<<totalEnergyDepositForAllStep<<G4endl;
  
  }

*/




/*
//kill photon when enter world volume

if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="World")
{
step->GetTrack()->SetTrackStatus(fStopAndKill);

G4cout<<"Killed world volume içinde"<<G4endl;

}


*/


/*
gamalar sintilator hacminden ayrıldığı anda track edilmesin (ikisinin farkını gör)
if(step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() && step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
if(step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() && step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="World")
step->GetTrack()->SetTrackStatus(fStopAndKill);
step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

*/


if(step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition() && step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
step->GetTrack()->SetTrackStatus(fStopAndKill);


if (!bc_400_logic) { 
    const G1DetectorConstruction* detectorConstruction
      = static_cast<const G1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    bc_400_logic = detectorConstruction->GetBc400();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
 if (volume != bc_400_logic) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep); 
  
  
  
  
 
  
  
  
 
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

