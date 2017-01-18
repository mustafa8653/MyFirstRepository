//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list o2f copyright holders.                             *
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
// $Id: B2EventAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file B2EventAction.cc
/// \brief Implementation of the B2EventAction class

#include "G1EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G1RunAction.hh"

#include "G1UserEventInformation.hh"

#include "HistoManager.hh" 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1EventAction::G1EventAction(G1RunAction *runAction, HistoManager* histo)
: G4UserEventAction(), fRunAction(runAction), fHistoManager(histo), fEdep(0.), fIndex(0.)
{




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1EventAction::~G1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1EventAction::BeginOfEventAction(const G4Event*)
{


G4cout<<"BeginOfEventAction started"<<G4endl;

//fpEventManager->SetUserInformation(new G1UserEventInformation); //fpEventManager protected attributes in G4UserEventAction class
G4EventManager::GetEventManager()->SetUserInformation(new G1UserEventInformation);
//anEvent->SetUserInformation(new G1UserEventInformation); bu olmaz cünkü const yetkisi yok

fEdep = 0.;
fIndex = 0.;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4cout<<"EndOfEventAction started"<<G4endl;



    G1UserEventInformation* eventInformation
    =(G1UserEventInformation*)anEvent->GetUserInformation();
    
    fEdep = eventInformation->GetEdep();
    
   // eventInformation->PrintProcessNumber();
    
    fRunAction->fillPerEvent(fEdep, fIndex);
    
    fHistoManager->FillHisto(0,fEdep);
    
    fHistoManager->FillNtuple1(fEdep);
    
   for( std::map< G4int, G4int>::iterator ii=eventInformation->GetProcessMap()->begin(); ii!=eventInformation->GetProcessMap()->end(); ++ii)

   {

   // G4cout << (*ii).first << ": " << (*ii).second << G4endl;
       
       G4int index= (*ii).first;
       G4int weight=(*ii).second;
       
       
        fHistoManager->FillHisto (1, index, weight);
        
  }
  
  
  
  
 
    
    
}  





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
