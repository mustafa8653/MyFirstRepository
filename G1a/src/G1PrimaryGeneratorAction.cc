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

#include "G1PrimaryGeneratorAction.hh"
#include "G4Box.hh"
//#include "G1PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4PhysicalVolumeStore.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1PrimaryGeneratorAction::G1PrimaryGeneratorAction()
//önce fBedBox initiliaze edilmeli fParticleGun dan önce yoksa hata veriyor
 : G4VUserPrimaryGeneratorAction(),fParticleGun(0),fBedBox(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //create a messenger for this class
 // fGunMessenger = new G1PrimaryGeneratorMessenger(this);


  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);
  
  
  fParticleGun->SetParticleEnergy(1115.546*keV);
  
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(G4UniformRand(),G4UniformRand(),G4UniformRand()));
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1PrimaryGeneratorAction::~G1PrimaryGeneratorAction()
{
  delete fParticleGun;
  //delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
/* 
//Source position
G4double sourcePositionx=0*cm;
G4double sourcePositiony=-(0.5*25+0.5*0.1+5)*cm;
G4double sourcePositionz=-0.5*91.5;
*/

G4VPhysicalVolume * physBc
      = G4PhysicalVolumeStore::GetInstance()->GetVolume("BC-400");


//world merkezine göre

//bu ikisinin bir farkı var mı? birisi refeerans döndürüyor..diğeri normal
G4ThreeVector vec=physBc->GetObjectTranslation();
//G4ThreeVector vec=physBc->GetTranslation();

//fParticleGun->SetParticlePosition(G4ThreeVector(sourcePositionx,sourcePositiony,sourcePositionz));

fParticleGun->SetParticlePosition(vec);

fParticleGun->GeneratePrimaryVertex(anEvent);
/*
if (!fBedBox)
  {
    G4LogicalVolume* bed_logic
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Bed");
    if ( bed_logic ) fBedBox = dynamic_cast<G4Box*>(bed_logic->GetSolid());
  }

if ( fBedBox ) {
   sourcePositiony= -10*fBedBox->GetYHalfLength();
    
  }  

 fParticleGun->SetParticlePosition(G4ThreeVector(sourcePositionx,sourcePositiony,sourcePositionz));

  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




//Bu iki fonksiyon sanırım geant içinde olusan photonları polarize ediyor

void G1PrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G1PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......