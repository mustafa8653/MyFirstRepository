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

//#include "G4LogicalBorderSurface.hh"

#include "G4RunManager.hh"

#include "G4TrackStatus.hh"  //fStopAndKill bunun içinde



#include "G4ParticleDefinition.hh" 
#include "G4ParticleTypes.hh"  //G4Gamma için

#include "G4OpBoundaryProcess.hh"

#include "G4EventManager.hh"
#include "G1UserEventInformation.hh"


#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::G1SteppingAction(HistoManager* histo)
: G4UserSteppingAction(), fHistoManager(histo), fEventInformation(0), bc_400_logic(0), bc_skin_logic(0), alu_border_logic(0), pmt_border_logic(0)
  
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::~G1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1SteppingAction::UserSteppingAction(const G4Step* step)
{

fEventInformation = (G1UserEventInformation*)G4EventManager::GetEventManager()->GetUserInformation();

if(step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition()) 
{
      if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) //son adım track ediliyor
      step->GetTrack()->SetTrackStatus(fStopAndKill);
}




 const G4VUserDetectorConstruction* vdc = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 const G1DetectorConstruction* dc = static_cast<const G1DetectorConstruction*>(vdc);


 G4VPhysicalVolume* physBc = dc->GetBcPV();
 //G4VPhysicalVolume* physAlu = dc->GetAluPV();
 //G4VPhysicalVolume* physPmt = dc->GetPmtPV();
 G4VPhysicalVolume* physWorld = dc->GetWorldPV();

  
  if (!bc_400_logic)  
   bc_400_logic = dc->GetBc400();   
   
  
  
  










if( step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )
{
/*
  if(!bc_skin_logic)
  bc_skin_logic = dc->GetBcSkinLogic();
  
  if(!alu_border_logic)
 alu_border_logic = dc->GetAluBorderLogic();
   
  if(!pmt_border_logic)
  pmt_border_logic = dc->GetPmtBorderLogic();
  
  
G4OpticalSurface *opBcSurface = dynamic_cast < G4OpticalSurface* >
        (bc_skin_logic->GetSurface(bc_400_logic)->GetSurfaceProperty());
 
 
 
G4OpticalSurface *opAluSurface = dynamic_cast < G4OpticalSurface* >
        (alu_border_logic->GetSurface(physBc,physAlu)->GetSurfaceProperty());    
    
    
 G4OpticalSurface *opPmtSurface = dynamic_cast <G4OpticalSurface*>
        (pmt_border_logic->GetSurface (physBc,physPmt)->GetSurfaceProperty()); 
  
  
 */
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager(); //photonun yapabileceği processler. 5 tane
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector *pv = pm->GetProcessList(); //obje tutan vektor sınıfı
    
    //This loop just for boundary process
    for(G4int i=0; i<nprocesses; i++) //return all process for each photon step
    
    {
    
    G4VProcess* vp=(*pv)[i];
   
     	if( vp->GetProcessName() == "OpBoundary" ) //for all surface
     	{
   	 
   	 G4VPhysicalVolume* p1 =  step->GetPreStepPoint()->GetPhysicalVolume();
         G4VPhysicalVolume* p2 =  step->GetPostStepPoint()->GetPhysicalVolume(); 
         
         
                 if( (p1==physBc && p2==physWorld) || (p1==physWorld && p2==physBc) ) 
                 {
                 
                 G4OpBoundaryProcess* boundary = (G4OpBoundaryProcess*)vp; 
                 G4OpBoundaryProcessStatus status = boundary->GetStatus(); 
         	 CountProcess(status);
         	 
                 }
    	
    	}
    
    
    }
    
    
 
               
           /*   
         G4VPhysicalVolume* p1 =  step->GetPreStepPoint()->GetPhysicalVolume();
         G4VPhysicalVolume* p2 =  step->GetPostStepPoint()->GetPhysicalVolume(); 
  
  
        if(p1==physBc && p2==physWorld)
         G4cout<<opBcSurface->GetName()<<G4endl;
  
       if(p1==physBc && p2==physAlu)
           G4cout<<opAluSurface->GetName()<<G4endl; 
           
        if(p1==physBc && p2==physPmt)
         G4cout<<opPmtSurface->GetName()<<G4endl;   
        */
     
    
 }   
   
  
      
      
      // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
 if (volume != bc_400_logic) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventInformation->AddEdep(edepStep); 
  
 
  
  
  
 
  
  
}

void G1SteppingAction::CountProcess( G4int status ){

         
   
  
   
               switch(status){
		case 0 :
		 (*(fEventInformation->GetProcessMap()))[0]+=1;
		 fHistoManager->FillNtuple2(0.);
		    break;
		case 1 :
		(*(fEventInformation->GetProcessMap()))[1]+=1;
		fHistoManager->FillNtuple2(1.);
		   break;
		case 2 :
		(*(fEventInformation->GetProcessMap()))[2]+=1;
		fHistoManager->FillNtuple2(2.);
		//G4cout <<"FresnelRefraction "<<G4endl;
		
		   break;
		case 3 :
		//G4cout <<"FresnelReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[3]+=1;
		fHistoManager->FillNtuple2(3.);
		break; 
		
	        case 4 :
		//G4cout <<"TotalInternalReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[4]+=1;
		fHistoManager->FillNtuple2(4.);
		   break;  
		case 5 :
		//G4cout <<"LambertianReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[5]+=1;
		fHistoManager->FillNtuple2(5.);
		   break;  
		case 6 :
		
		//G4cout <<"LobeReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[6]+=1;
		fHistoManager->FillNtuple2(6.);
		   break;
		case 7 :
		//G4cout <<"SpikeReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[7]+=1;
		fHistoManager->FillNtuple2(7.);
		//G4cout<<"number SS: "<<fEventAction->CountProcessVector()->at(7);
		   break;
		case 8 :
		//G4cout <<"BackScattering"<<G4endl;
		(*(fEventInformation->GetProcessMap()))[8]+=1;
		fHistoManager->FillNtuple2(8.);
		   break;    
		case 9 :
		//G4cout <<"Absorption "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[9]+=1;
		fHistoManager->FillNtuple2(9.);
		   break;   
		case 10 :
		//G4cout <<"Detection"<<G4endl;
		(*(fEventInformation->GetProcessMap()))[10]+=1;
		fHistoManager->FillNtuple2(10.);
		   break; 
		case 13 :
		//G4cout <<"StepTooSmall "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[13]+=1;
		fHistoManager->FillNtuple2(13.);
		
		   break; 
		   
		     
		default :
		break;
               }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

