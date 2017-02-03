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
//o
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G1DetectorConstruction_h
#define G1DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G1DetectorMessenger.hh"
#include "G4OpticalSurface.hh"
#include <vector>

class G4VPhysicalVolume;
class G4VisAttributes;
class G4LogicalVolume;
class G4OpticalSurface;
class G4LogicalSkinSurface;
class G4LogicalBorderSurface;
class G4Material;

class G4Box;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    G1DetectorConstruction();
    virtual ~G1DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetBc400() const { return fBc_400_logic; }
    
    G4LogicalSkinSurface* GetBcSkinLogic() const { return fBc_skin_logic; }
    G4LogicalBorderSurface* GetAluBorderLogic() const { return fAlu_border_logic; }
    G4LogicalBorderSurface* GetPmtBorderLogic() const { return fPmt_border_logic; }
    
    G4VPhysicalVolume* GetBcPV() const { return fPhysBc; }
    G4VPhysicalVolume* GetAluPV() const { return fPhysAlu; }
    G4VPhysicalVolume* GetPmtPV() const { return fPhysPmt; }
    G4VPhysicalVolume* GetWorldPV() const { return fPhysWorld; }
    
    
    void SetDetectorDimension( G4ThreeVector tv );
    void SetOpFinishType( G4int type );
    void SetSigmaAlpha( G4double );
    void SetReflectivity( G4double reflectivity );
    void SetResolutionScale( G4double rs );
    
    void SetDefaults();
    void Print() const;
    
    void SetP();
    void SetG();
    void SetBP(); // PBP and GBP
    void SetFP(); // PFP and GFP
    
  private:
  
  
  
  void DefineMaterials();
   G4VPhysicalVolume* DefineVolumes();
  
  G1DetectorMessenger* fDetectorMessenger;
  
  
G4double fBc_400_size_x;
G4double fBc_400_size_y;
G4double fBc_400_size_z;


	//World Size
G4double fWorldSize_x;
G4double fWorldSize_y;
G4double fWorldSize_z;
  
  G4double fScintillationYield;
  G4double fResolutionScale;
  G4double fFastTimeConstant;
  G4double fYieldRatio;
  
  G4double fSigmaAlpha;
  G4double fReflectivity;
  
  G4Material *fWorld_mat;
  G4Material *fBc_400_mat;
  
  G4Box* fWorld_solid;
  G4Box* fBc_400_solid;
  
  G4LogicalVolume *fWorld_logic;
  G4LogicalVolume *fBc_400_logic;
  
  G4LogicalSkinSurface* fBc_skin_logic;
  G4LogicalBorderSurface* fAlu_border_logic;
  G4LogicalBorderSurface* fPmt_border_logic;
  
  
  G4VPhysicalVolume* fPhysWorld;
  G4VPhysicalVolume *fPhysBc;
  G4VPhysicalVolume *fPhysAlu;
  G4VPhysicalVolume *fPhysPmt;
  
  
  G4OpticalSurface* fOpBcSurface; 
  G4OpticalSurfaceFinish fSurfaceFinish;
  
  
  
  
   G4bool checkOverlaps;
    
   std::vector<G4VisAttributes*> fVisAttributes;




};


























//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*G1DetectorConstruction_h*/
