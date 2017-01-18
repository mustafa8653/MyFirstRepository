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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"


#include "G1RunAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4ParameterManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo................oooOO0OOooo........oooOO0OOooo......

G1RunAction::G1RunAction(HistoManager* histo)
 : G4UserRunAction(), fTimer(0) , fHistoManager(histo), fEdep(0.), fIndex(0)
{

  
  fTimer = new G4Timer;
  
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1RunAction::~G1RunAction()
{
  delete fTimer;
  
 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1RunAction::BeginOfRunAction(const G4Run* )
{

G4cout << "BeginOfRunAction Started" << G4endl;

fTimer->Start();
  
  fEdep = 0.;
  fIndex = 0;
  
  fHistoManager->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1RunAction::fillPerEvent(G4double edep, G4int index)
                                  
{
  //accumulate statistic
  //
  fEdep += edep;
  fIndex += index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1RunAction::EndOfRunAction(const G4Run* run)
{

G4cout << "EndOfRunAction Started" << G4endl;



G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;


  fTimer->Stop();
  G4cout << "number of event = " << nofEvents
         << " " << *fTimer << G4endl;
         
         
 
         
  //save histograms
 // fHistoManager->PrintStatistic();
  fHistoManager->Save();        

  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Runnn-----------------------";
     
    
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------"
     << G4endl;
     
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << "  event"
     << G4endl;
     
  
   
  
         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



