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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::G1SteppingAction()
: G4UserSteppingAction()
  
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::~G1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1SteppingAction::UserSteppingAction(const G4Step* step)
{
  
  
  /*
  
  //bir stepde ne kadar yeni track olduğunu belirtiyor
  //G4cout<<"GetNumberOfSecondariesInCurrentStep(): "<<step->GetNumberOfSecondariesInCurrentStep()<<G4endl;
  
  */
  
  
  
  
  /*
  //  step has to point pre and post...step refers to poststep..this two line is the same
  G4cout<<"Step kinetic energy: "<<step->GetTrack()->GetKineticEnergy ()<<G4endl;
  G4cout<<"post step:kinetic "<<step->GetPostStepPoint()->GetKineticEnergy()<<G4endl;
  */
  
  
  /*
  track nerden baslıo onu gösterio
  G4cout<<"Step#: "<<step->GetTrack()->GetVertexPosition ()<<G4endl;
  G4cout<<"Step#: "<<step->GetTrack()->GetVertexKineticEnergy ()<<G4endl;
 */
 
 //GetTrack
 G4Track *track=step->GetTrack();
 
 //Get pre and post step point
 G4StepPoint* preStepPoint=step->GetPreStepPoint();
 G4StepPoint* postStepPoint=step->GetPostStepPoint();
 
 
 //Step Number
 G4int currentStepNumber=track->GetCurrentStepNumber();
 G4int preStepNumber=track->GetCurrentStepNumber()-1;
 
 //step 0 için sıfır olacak tüm durumlar için bu değişken kullanılacak
 G4double initialStep=0;
 
 
 
 
 
 
 //StepLenth
 G4double stepLength=step->GetStepLength();
 
 
 //TrackLength
 G4double trackLength=track->GetTrackLength();
 
 //Deposited enegy in stepLength
 G4double dE=step->GetTotalEnergyDeposit();
 
 
 //PreStepPoint,,,  sadece step sıfır için kullanılıyor
 G4double preX=preStepPoint->GetPosition().x();
 G4double preY=preStepPoint->GetPosition().y();
 G4double preZ=preStepPoint->GetPosition().z();
 
 G4double postX=postStepPoint->GetPosition()[0];
 G4double postY=postStepPoint->GetPosition()[1];
 G4double postZ=postStepPoint->GetPosition()[2];
 
 
 //Kitetic energy
 G4double preEnergy=preStepPoint->GetKineticEnergy(); //just for initial step
 G4double postEnergy=postStepPoint->GetKineticEnergy();
 
 
 //PhysicalVolume
    G4String initialPV=preStepPoint->GetPhysicalVolume()->GetName();
    G4String postPV=postStepPoint->GetPhysicalVolume()->GetName(); //this equal to track->GetPhysicalVolume()->GetName()
 
 
 //ProcessName
G4String preProcName="initilStep";
G4String postProcName=postStepPoint->GetProcessDefinedStep()->GetProcessName();
 
 
 
 
 
 
 //
 G4int numberOfSecondariesInCurrentStep=step->GetNumberOfSecondariesInCurrentStep();
 
 
 
 
 
 //step refers to postStepPoint
 if(step->GetTrack()->GetCurrentStepNumber()==1)
 {
 
 //track başlangıcı
 
 //Header
 G4cout<<"************************************************************************************************************** "<<G4endl;
  G4cout<<"  G4Track Information \t Particle = "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<",\t TrackID= "<<step->GetTrack()->GetTrackID()<<",\t ParentID= "<<step->GetTrack()->GetParentID()<<G4endl;
  G4cout<<"************************************************************************************************************* "<<G4endl;
 
 
 //7 satır bosluk+diger karakter..12=3(spacenumber)+karakternumber
 G4cout<<"Step#"<<std::setw(8)<<"X(mm)"<<std::setw(8)<<"Y(mm)"<<std::setw(8)<<"Z(mm)"<<std::setw(12)<<"KinE(MeV)"<<std::setw(10)<<"dE(MeV)" <<std::setw(13)<< "stepLength" <<std::setw(14)<<"trackLength" <<std::setw(13)<<"nextVolume"<<std::setw(16)<<"ProcName"<<std::setw(13)<<"#Secondary"<<G4endl;
 
 //step 0
 G4cout<<std::setw(5)<<std::setprecision(3)<<preStepNumber<<std::setw(8)<<preX<<std::setw(8)<<preY<<std::setw(8)<<preZ<<std::setw(12)<<preEnergy<<std::setw(10)<<initialStep<<std::setw(13)<<initialStep<<std::setw(14)<<initialStep<<std::setw(13)<<initialPV<<std::setw(16)<<preProcName<<std::setw(13)<<initialStep<<G4endl;
 
 //step 1, step refers to poststep
 G4cout<<std::setw(5)<<std::setprecision(3)<<currentStepNumber<<std::setw(8)<<postX<<std::setw(8)<<postY<<std::setw(8)<<postZ<<std::setw(12)<<postEnergy<<std::setw(10)<<dE<<std::setw(13)<<stepLength<<std::setw(14)<<trackLength<<std::setw(13)<<postPV<<std::setw(16)<<postProcName<<std::setw(13)<<numberOfSecondariesInCurrentStep<<G4endl;
 
 
 
 
 
 }
 
 else{
 
 G4cout<<std::setw(5)<<std::setprecision(3)<<currentStepNumber<<std::setw(8)<<postX<<std::setw(8)<<postY<<std::setw(8)<<postZ<<std::setw(12)<<postEnergy<<std::setw(10)<<dE<<std::setw(13)<<stepLength<<std::setw(14)<<trackLength<<std::setw(13)<<postPV<<std::setw(16)<<postProcName<<std::setw(13)<<numberOfSecondariesInCurrentStep<<G4endl;
 

  
  }
  
  

//G4cout<<"Track id: "<<step->GetTrack()->GetTrackID()<<G4endl;


/*
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
     */ 
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

