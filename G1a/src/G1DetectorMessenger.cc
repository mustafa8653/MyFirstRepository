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
// $Id: B2aDetectorMessenger.cc 69706 2013-05-13 09:12:40Z gcosmo $
// 
/// \file B2aDetectorMessenger.cc
/// \brief Implementation of the B2aDetectorMessenger class

#include "G1DetectorMessenger.hh"
#include "G1DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1DetectorMessenger::G1DetectorMessenger(G1DetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fGDirectory = new G4UIdirectory("/G/");
  fGDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/G/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  
  
  
  fOpTypeCmd = new G4UIcmdWithAnInteger("/G/det/type",this);
  fOpTypeCmd->SetGuidance("Set the surface finish type.");
  fOpTypeCmd->SetParameterName("type",false);
  fOpTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fOpTypeCmd->SetToBeBroadcasted(false);

  fSigmaAlphaCmd = new G4UIcmdWithADouble("/G/det/sigmaAlpha",this);
  fSigmaAlphaCmd->SetGuidance("Set the sigmaAlpga of the surface .");
  fSigmaAlphaCmd->SetParameterName("sigmaAlpha",false);
  fSigmaAlphaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSigmaAlphaCmd->SetToBeBroadcasted(false);


  fReflectivityCmd = new G4UIcmdWithADouble("/G/det/reflectivity",this);
  fReflectivityCmd->SetGuidance("Set the reflectivity of the wrapping .");
  fReflectivityCmd->SetParameterName("reflectivity",false);
  fReflectivityCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fReflectivityCmd->SetToBeBroadcasted(false);
  
  fResolutionScaleCmd = new G4UIcmdWithADouble("/G/det/resolutionScale",this);
  fResolutionScaleCmd->SetGuidance("Set the resolutionScale .");
  fResolutionScaleCmd->SetParameterName("resolutionScale",false);
  fResolutionScaleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fResolutionScaleCmd->SetToBeBroadcasted(false);
  
  fDetecDimensionsCmd =
    new G4UIcmdWith3VectorAndUnit("/G/det/dimensions",this);
  fDetecDimensionsCmd->SetGuidance("Set the dimensions of the detector volume.");
  fDetecDimensionsCmd->SetParameterName("x","y","z",false);
  fDetecDimensionsCmd->SetDefaultUnit("cm");
  fDetecDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDetecDimensionsCmd->SetToBeBroadcasted(false);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1DetectorMessenger::~G1DetectorMessenger()
{
  
  
  delete fGDirectory;
  delete fDetDirectory;
  
  delete fOpTypeCmd;
  delete fSigmaAlphaCmd;
  delete fReflectivityCmd;
  delete fResolutionScaleCmd;
  delete fDetecDimensionsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
      
      
          if( command == fOpTypeCmd)
            {
              
              fDetectorConstruction->SetOpFinishType( fOpTypeCmd->GetNewIntValue( newValue ) );
           
            }
            
          if( command == fSigmaAlphaCmd )
          {
            fDetectorConstruction->SetSigmaAlpha( fSigmaAlphaCmd->GetNewDoubleValue(newValue) );
          }
      
	 if( command == fReflectivityCmd ) 
	  {
	    fDetectorConstruction->SetReflectivity( fReflectivityCmd->GetNewDoubleValue(newValue) );
	  }  
	  
	   if( command == fResolutionScaleCmd ) 
	  {
	    fDetectorConstruction->SetResolutionScale( fResolutionScaleCmd->GetNewDoubleValue(newValue) );
	  }  
	  
	   if( command == fDetecDimensionsCmd ) 
	  {
	    fDetectorConstruction->SetDetectorDimension( fDetecDimensionsCmd->GetNew3VectorValue(newValue) );
	  } 
  
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
