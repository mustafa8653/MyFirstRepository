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
: G4UserEventAction(), fRunAction(runAction), fHistoManager(histo), fEdep(0.), fIndex(0), fEmittedPhotonNumber(0), fTotalDetectedPhotonNumber(0), fPmt1DetCount(0), fPmt2DetCount(0), fPmt1AbsCount(0), fPmt2AbsCount(0)
{




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1EventAction::~G1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1EventAction::BeginOfEventAction(const G4Event* anEvent)
{


G4int eventID = anEvent->GetEventID();
//G4cout<<"BeginOfEventAction started"<<G4endl;

//fpEventManager->SetUserInformation(new G1UserEventInformation); //fpEventManager protected attributes in G4UserEventAction class
G4EventManager::GetEventManager()->SetUserInformation( new G1UserEventInformation ( eventID ));
//anEvent->SetUserInformation(new G1UserEventInformation); bu olmaz cünkü const yetkisi yok

fEdep = 0.;
fIndex = fEmittedPhotonNumber = fTotalDetectedPhotonNumber = fPmt1DetCount = fPmt2DetCount = fPmt1AbsCount = fPmt2AbsCount = 0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1EventAction::EndOfEventAction (const G4Event* anEvent)
{
     // G4cout<<"EndOfEventAction started"<<G4endl;

   /* if( anEvent->IsAborted() == true ) return; */
    
    G1UserEventInformation* eventInformation = (G1UserEventInformation*)anEvent->GetUserInformation();
    
 //eventInformation->PrintEventResult();
    
    
    fEdep = eventInformation->GetEdep();
   
    fEmittedPhotonNumber =  eventInformation->GetEmittedPhotonNumber();
    
    fPmt1DetCount = eventInformation->GetPmt1DetCount();
    fPmt2DetCount = eventInformation->GetPmt2DetCount();
    fPmt1AbsCount = eventInformation->GetPmt1AbsCount();
    fPmt2AbsCount = eventInformation->GetPmt2AbsCount();
    
    
   /*
    std::vector< double> pev = *(eventInformation->GetPhotonEnergyVec()); //PrimaryVertexPhotonEnergy,Emission spectrum
    std::vector< double> pdpv = *(eventInformation->GetDetectedPmtPhotonEnergyVec()); //PmtDetectedPhotonEnergy
    std::vector< int> csoc = *(eventInformation->GetScinOrCherVec()); //identify Scin or Cherenkov photon
    
   */ 
    
    
    //fHistoManager->FillNtuple2(  pev, pdpv, csoc);
    
    
     fHistoManager->FillNtuple3(fEdep, fEmittedPhotonNumber, fPmt1DetCount, fPmt2DetCount, fPmt1AbsCount, fPmt2AbsCount );
     fRunAction->fillPerEvent(fEdep, fIndex);
    
    
    
    
    
   for( std::map< G4int, G4int>::iterator ii = eventInformation->GetProcessMap()->begin(); ii != eventInformation->GetProcessMap()->end(); ++ii)

   {

   // G4cout << (*ii).first << ": " << (*ii).second << G4endl;
       
       G4int index = (*ii).first;
       G4int weight = (*ii).second;
       
       
          for(G4int i=0; i< weight; i++)
          {
           
              
              fHistoManager->FillNtuple1(index);     
                
              
           
           }
        
        
        
        
   }
  
  
  
              
              
              
        
  
 
    
    
}  





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
