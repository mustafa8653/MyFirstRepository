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

#include "G4OpBoundaryProcess.hh"

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







//kill photon when enter world volume
/*
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


if(step->GetTrack()->GetDefinition()==G4Gamma::GammaDefinition()) 
{
      if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
      step->GetTrack()->SetTrackStatus(fStopAndKill);
}

/*
if(step->GetTrack()->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
{
const G4VProcess *vp=step->GetPostStepPoint()->GetProcessDefinedStep();
G4OpBoundaryProcess * bp=(G4OpBoundaryProcess*)vp;
G4cout <<"process name for photon: "<<vp->GetProcessName()<<G4endl;
G4cout<<"Status :"<<bp->GetStatus()<<G4endl;
}
*/

//For boundaryprocess



//burada reflectivity type ve refraction ön yüzey için...absorption 
//if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
if(step->GetTrack()->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
{
  G4VPhysicalVolume* p1 =  step->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* p2 =  step->GetPostStepPoint()->GetPhysicalVolume(); 
  
 if(p1!=p2){ //transport olanlar
    
    
    /*
       const G4String& pv1 =  p1->GetName();
       const G4String& pv2 =  p2->GetName(); 
  
         if(pv1=="BC-400" && pv2=="aluPV")
            G4cout<<"Bc400-Alu surface"<<G4endl;
         if(pv1=="BC-400" && pv2=="photocathodePV")
            //G4cout<<"Bc400-PhoCat surface"<<G4endl;
         if(pv1=="BC-400" && pv2=="World")
            G4cout<<"Bc400-World surface"<<G4endl;   
    */
    
    /*
    G4cout<<"Old MomentumDirection: "<<step->GetPreStepPoint()->GetMomentumDirection()<<G4endl;
    G4cout<<"Old Polarization: "<<step->GetPreStepPoint()->GetPolarization()<<G4endl;
    G4cout<<"New MomentumDirection: "<<step->GetTrack()->GetDynamicParticle()->GetMomentumDirection()<<G4endl;
    G4cout<<"New Polarization: "<<step->GetTrack()->GetDynamicParticle()->GetPolarization()<<G4endl;
    */
   
    
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager(); //photonun yapabileceği processler
    G4int nprocesses = pm->GetProcessListLength();
    //G4cout<<"number of OpticPhoton process: "<<nprocesses<<G4endl; //photon için tanmlı tüm processler
    
    G4ProcessVector* pv = pm->GetProcessList();
    
    
    
    for(G4int i=0; i<nprocesses; i++) //return all process for each step
    
    {
    //G4cout<<"Process : "<<i<<" "<<(*pv)[i]->GetProcessName()<<G4endl;
    
    G4VProcess* vp=(*pv)[i];
    G4OpBoundaryProcess* boundary = (G4OpBoundaryProcess*)vp; //G4OpBoundaryProcess can derive from all process
   
   
    
   
	   if(boundary->GetProcessName()=="OpBoundary" && boundary->GetProcessType()==fOptical) //photon için tanımlı birkaç process var..5 tane
	   {
   
            G4OpBoundaryProcessStatus status=boundary->GetStatus();
   
	     switch(status){
		case 0 :
		 fEventAction->CountProcessVector()->at(0)++;
		    break;
		case 1 :
		fEventAction->CountProcessVector()->at(1)++;
		   break;
		case 2 :
		fEventAction->CountProcessVector()->at(2)++;
		//G4cout <<"FresnelRefraction "<<G4endl;
		
		   break;
		case 3 :
		//G4cout <<"FresnelReflection "<<G4endl;
		fEventAction->CountProcessVector()->at(3)++;
		   break; 
		case 4 :
		//G4cout <<"TotalInternalReflection "<<G4endl;
		fEventAction->CountProcessVector()->at(4)++;
		   break;   
		case 6 :
		
		//G4cout <<"LobeReflection "<<G4endl;
		fEventAction->CountProcessVector()->at(6)++;
		   break;
		case 7 :
		//G4cout <<"SpikeReflection "<<G4endl;
		fEventAction->CountProcessVector()->at(7)++;
		//G4cout<<"number SS: "<<fEventAction->CountProcessVector()->at(7);
		   break;
		case 8 :
		//G4cout <<"BackScattering"<<G4endl;
		fEventAction->CountProcessVector()->at(8)++;
		   break;    
		case 9 :
		//G4cout <<"Absorption "<<G4endl;
		fEventAction->CountProcessVector()->at(9)++;
		   break;   
		case 10 :
		//G4cout <<"Detection"<<G4endl;
		fEventAction->CountProcessVector()->at(10)++;
		   break; 
		case 13 :
		//G4cout <<"StepTooSmall "<<G4endl;
		fEventAction->CountProcessVector()->at(13)++;
		   break;     
		default :
		break;
               }
    	}
   /*
    if( boundary->GetStatus()!=Undefined)
    { //undefined değilse
          
    //step->GetTrack()->SetTrackStatus(fStopAndKill);
    G4cout<<"boundary status: "<<boundary->GetStatus()<<G4endl;
    if(boundary->GetStatus()==LobeReflection)c++;
    //G4cout<<boundary->GetProcessName()<<G4endl;
    }
     
     */  
  }
      
      
      
 //    G4cout<<"dsfdsfdsf:   "<<countProcess.size()<<G4endl; 
      
    
    
}




}


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

