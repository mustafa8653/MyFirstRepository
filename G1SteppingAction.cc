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

#include "G4UImanager.hh"

#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::G1SteppingAction(HistoManager* histo)
: G4UserSteppingAction(), fHistoManager(histo), fEventInformation(0), fScin_logic(0), fScin_skin_logic(0), alu_border_logic(0), pmt_border_logic(0)
  
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1SteppingAction::~G1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1SteppingAction::UserSteppingAction(const G4Step* step)
{

 
      
 fEventInformation = (G1UserEventInformation*)G4EventManager::GetEventManager()->GetUserInformation();
 
 


 const G4VUserDetectorConstruction* vdc = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 const G1DetectorConstruction* dc = static_cast<const G1DetectorConstruction*>(vdc);


 //G4VPhysicalVolume* physBc = dc->GetBcPV();
 //G4VPhysicalVolume* physAlu = dc->GetAluPV();
// G4VPhysicalVolume* physPmt = dc->GetPmtPV();
 //G4VPhysicalVolume* physWorld = dc->GetWorldPV();

  
  /*
  if(step->GetTrack()->GetDefinition() == G4Gamma::GammaDefinition() && step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary )
     step->GetTrack()->SetTrackStatus(fStopAndKill);
  */
  
  if (!fScin_logic)  
   fScin_logic = dc->GetScinLV();   
   
  
  
  





if( step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) //particle is photon
{

 

  
    
  
    G4OpBoundaryProcessStatus boundaryStatus = Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL;
    
    
    //find the boundary process only once
    if(!boundary){  //bir kez çalışıyor sadece. oda başda
    //G4cout<<" boundarynooo"<<G4endl;
    
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager(); //photonun yapabileceği processler. 5 tane
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector *pv = pm->GetProcessList(); //obje tutan vektor sınıfı
    
                  for(G4int i=0; i<nprocesses; i++)
                  {
                     if((*pv)[i]->GetProcessName()=="OpBoundary")
                         {
                           boundary = (G4OpBoundaryProcess*)(*pv)[i];
                           
                              break;
                         }
                 }
    }


   if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary)
   {  //photon transport processi ise diyede yapılabilir.aksi taktirde Opabsorption da dahil oluyor..
   boundaryStatus = boundary->GetStatus();
    
      if(boundaryStatus == 10){
      
      fEventInformation->GetDetectedPmtPhotonEnergyVec()->push_back(step->GetTrack()->GetTotalEnergy() ) ;
                           
                           if(step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() ==0)
                           fEventInformation->IncPmt1DetPhotonCount();
                           
                           if(step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() ==1)
                           fEventInformation->IncPmt2DetPhotonCount();
      }
       
     if(boundaryStatus == 9){
      
      
       if(step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Pmt")  
       {         
                           if(step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() ==0)
                           fEventInformation->IncPmt1AbsPhotonCount();
                           
                           if(step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() ==1)
                           fEventInformation->IncPmt2AbsPhotonCount();
      
      }
      
      
      }
    
   
   
      CountProcess(boundaryStatus);
      //G4cout <<"statuss: "<<boundaryStatus<<G4endl;                     
   
   
   }   
 
  
  
  
//scintilator içinde absorblanan foton sayısı
  if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="OpAbsorption")
  fEventInformation->IncAbsorbedPhotonCount();
    
    //OpMie ve OpRaylei dedectorConstruction sınıfı içindeki parametreleri değiştirirsek aktif oluyor
    
    
 
               
           
     
    
 }   
   
  
      
      
      // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
 if (volume != fScin_logic) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
   
  
 if( step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() ) //optic fotonların üretiminde enerji korunumu yok
  fEventInformation->AddEdep(edepStep);
  
  
 
  
  
}

void G1SteppingAction::CountProcess( G4int status )
{

         
                switch(status){
		case 0 :
		 (*(fEventInformation->GetProcessMap()))[0]+=1;
		 
		    break;
		case 1 :
		(*(fEventInformation->GetProcessMap()))[1]+=1;
		
		   break;
		case 2 :
		(*(fEventInformation->GetProcessMap()))[2]+=1;
		
		//G4cout <<"FresnelRefraction "<<G4endl;
		
		   break;
		case 3 :
		//G4cout <<"FresnelReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[3]+=1;
		
		break; 
		
	        case 4 :
		//G4cout <<"TotalInternalReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[4]+=1;
		
		   break;  
		case 5 :
		//G4cout <<"LambertianReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[5]+=1;
		
		   break;  
		case 6 :
		
		//G4cout <<"LobeReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[6]+=1;
		
		   break;
		case 7 :
		//G4cout <<"SpikeReflection "<<G4endl;
		(*(fEventInformation->GetProcessMap()))[7]+=1;
		
		   break;
		case 8 :
		//G4cout <<"BackScattering"<<G4endl;
		(*(fEventInformation->GetProcessMap()))[8]+=1;
		
		   break;    
		case 9 :
		//G4cout <<"Absorption "<<G4endl; 
		(*(fEventInformation->GetProcessMap()))[9]+=1;
		
		   break;   
		case 10 :
		//G4cout <<"Semele"<<G4endl; //detection absorbtion içinde değil..yani detekte olan aynı zmanda absorbe gözükmüyor
		(*(fEventInformation->GetProcessMap()))[10]+=1;
		
		   break; 
		case 13 : //this step works also OpAbsorption process,not just for after reflection
		//G4cout <<"StepTooSmall semele"<<G4endl;
		(*(fEventInformation->GetProcessMap()))[13]+=1;
		
		
		   break; 
		   
		     
		default :
		
		break;
               }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

