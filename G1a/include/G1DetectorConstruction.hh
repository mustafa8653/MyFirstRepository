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

#include <vector>

class G4VPhysicalVolume;
class G4VisAttributes;
class G4LogicalVolume;
class G4OpticalSurface;
class G4LogicalSkinSurface;
class G4LogicalBorderSurface;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    G1DetectorConstruction();
    virtual ~G1DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetBc400() const { return bc_400_logic; }
    
    G4LogicalSkinSurface* GetBcSkinLogic() const { return bc_skin_logic; }
    G4LogicalBorderSurface* GetAluBorderLogic() const { return alu_border_logic; }
    G4LogicalBorderSurface* GetPmtBorderLogic() const { return pmt_border_logic; }
    
    G4VPhysicalVolume* GetBcPV() const { return physBc; }
    G4VPhysicalVolume* GetAluPV() const { return physAlu; }
    G4VPhysicalVolume* GetPmtPV() const { return physPmt; }
    G4VPhysicalVolume* GetWorldPV() const { return physWorld; }
  private:
  void DefineMaterials();
   G4VPhysicalVolume* DefineVolumes();
  
  
  G4VPhysicalVolume *physBc;
  G4VPhysicalVolume *physAlu;
  G4VPhysicalVolume *physPmt;
  G4VPhysicalVolume* physWorld;
  
  G4LogicalVolume *bc_400_logic;
  
  
  G4LogicalSkinSurface* bc_skin_logic;
  G4LogicalBorderSurface* alu_border_logic;
  G4LogicalBorderSurface* pmt_border_logic;
  
   G4bool checkOverlaps;
    
   std::vector<G4VisAttributes*> fVisAttributes;



};


























//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*G1DetectorConstruction_h*/
