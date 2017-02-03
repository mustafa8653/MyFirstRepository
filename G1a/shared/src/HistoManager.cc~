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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 92443 2015-09-01 13:56:16Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
 : fFactoryOn(false), v()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{


  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile("root/gate");
  if (! fileOpen) {
    G4cerr << "\n---> HistoManager::Book(): cannot open " 
           << analysisManager->GetFileName() << G4endl;
    return;
  }
  
  // Create histograms.
  // Histogram ids are generated automatically starting from 0.
  // The start value can be changed by:
  // analysisManager->SetFirstHistoId(1);  
  
  // id = 0
  analysisManager->CreateH1("EAbs", "Edep in absorber (MeV)", 100, 0., 1.5*MeV);
  // id = 1
  analysisManager->CreateH1("ProcessNumber", "NumberOfProcess", 13, 0.5, 13.5);
  // id = 2
  analysisManager->CreateH1("PmtPhotonDetection", "PmtPhotonCount", 200, 0., 600);
  

  // Create ntuples.
  // Ntuples ids are generated automatically starting from 0.
  // The start value can be changed by:
  // analysisManager->SetFirstMtupleId(1);  
  
  // Create 1st ntuple (id = 0)
  analysisManager->CreateNtuple("Ntuple1", "Ntuple");
  analysisManager->CreateNtupleDColumn("Eabs"); // column Id = 0
  analysisManager->CreateNtupleIColumn("PmtPhotonCount"); // column Id = 1
  analysisManager->FinishNtuple();

  // Create 2nd ntuple (id = 1)
  analysisManager->CreateNtuple("Ntuple2", "ProcessTypee");
  analysisManager->CreateNtupleIColumn("ProcessType"); // column Id = 0
  
  analysisManager->FinishNtuple();
  
  // Create 3 ntuple (id = 2)
  analysisManager->CreateNtuple("Ntuple3", "PhotonPrimaryVertexEnergy");
  analysisManager->CreateNtupleDColumn("PhotonEnergy", v); // column Id = 0
  analysisManager->FinishNtuple();
  
  fFactoryOn = true;       

  G4cout << "\n----> Output file is open in " 
         << analysisManager->GetFileName() << "." 
         << analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{
  if (! fFactoryOn) return;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->Write();
  analysisManager->CloseFile(); 
   
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
      
  delete G4AnalysisManager::Instance();
  fFactoryOn = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto (G4int ih, G4double xbin , G4double weight )
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void HistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  G4H1* h1 = analysisManager->GetH1(ih);
  if (h1) h1->scale(fac);
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple1(G4double edep, G4int pfc)
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 0)
  analysisManager->FillNtupleDColumn(0, 0, edep);
  analysisManager->FillNtupleIColumn(0, 1, pfc);
  analysisManager->AddNtupleRow(0);  
}  


void HistoManager::FillNtuple2(G4int index )
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 2st ntuple ( id = 1)
  analysisManager->FillNtupleIColumn(1, 0, index);
  analysisManager->AddNtupleRow(1);  
}  

void HistoManager::FillNtuple3(std::vector< double > &vector )
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 3 ntuple ( id = 2)
  v=vector;
  //G4cout<<"SEMELE"<<v[0]<<G4endl;
  //analysisManager->FillNtupleDColumn(2, 0,v); //bu yok cunku adresi vermisdik yaratırken
  analysisManager->AddNtupleRow(2); // AddNtupleRow(nTupleID) 
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void HistoManager::PrintStatistic()
{
  if (! fFactoryOn) return;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for ( G4int i=0; i<analysisManager->GetNofH1s(); ++i ) {
    G4String name = analysisManager->GetH1Name(i);
    G4H1* h1 = analysisManager->GetH1(i);
    
    G4String unitCategory;
    if (name[0U] == 'E' ) unitCategory = "Energy"; 
    if (name[0U] == 'L' ) unitCategory = "Length";
         // we use an explicit unsigned int type for operator [] argument
         // to avoid problems with windows compiler

    G4cout << name
           << ": mean = " << G4BestUnit(h1->mean(), unitCategory) 
           << " rms = " << G4BestUnit(h1->rms(), unitCategory ) 
           << G4endl;
  }
}

*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


