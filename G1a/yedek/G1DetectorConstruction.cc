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

#include "G1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"

#include "G4VisAttributes.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1DetectorConstruction::G1DetectorConstruction()
 :  G4VUserDetectorConstruction(), fWorld_mat(0), fBc_400_mat(0),
                                   fWorld_solid(0), fBc_400_solid(0),
                                   fWorld_logic(0), fBc_400_logic(0), fBc_skin_logic(0), fAlu_border_logic(0), fPmt_border_logic(0),
                                   fPhysWorld(0), fPhysBc(0), fPhysAlu(0), fPhysPmt(0),   
                                   checkOverlaps(true), fVisAttributes()

{

  SetDefaults();
  fDetectorMessenger = new G1DetectorMessenger(this);
  
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1DetectorConstruction::~G1DetectorConstruction()
{


      for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
         {
           delete fVisAttributes[i];
         }  
         
         delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G1DetectorConstruction::Construct()
{

   DefineMaterials();
   
  // Define volumes
   return DefineVolumes();

 
}


void G1DetectorConstruction::DefineMaterials(){

G4double a, z, density;
  G4int nelements;


/*
# of hydrogen atoms per cm3 5.23*10^{22}
# of carbon atoms per cm3 4.74*10^{22}

1 cm3 de kac gr hidrojen ve carbon var
(5.23*10^22/6.02*10^23)*1.00794*g/mole=
(4.74*10^22/6.02*10^23)*12.0107*g/mole=


*/
//BC-400 PlasticScintillator
   G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.00794*g/mole);
   G4Element* C = new G4Element("Carbon", "C", z=6 , a=12.0107*g/mole);

   fBc_400_mat = new G4Material("BC-400", density=1.032*g/cm3, nelements=2);
   fBc_400_mat->AddElement(H,8.5*perCent);
   fBc_400_mat->AddElement(C,91.5*perCent);

// Air

  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fWorld_mat = new G4Material("Air", density=1.29*g/cm3, nelements=2);
  fWorld_mat->AddElement(N, 70.*perCent);
  fWorld_mat->AddElement(O, 30.*perCent);




/*
// Water

 G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

*/
//Aluminium
//G4Material* Al = new G4Material(name="Aluminum", z=13., a = 26.98*g/mole, density = 2.700*g/cm3);

/*
//Stainless Steel

       G4NistManager* nist = G4NistManager::Instance();
        //G4Element* C  = nist->FindOrBuildElement("C");
        G4Element* Si =nist->FindOrBuildElement("Si");
        
        G4Element* Cr = nist->FindOrBuildElement("Cr");
	G4Element* Mn = nist->FindOrBuildElement("Mn");
	G4Element* Fe = nist->FindOrBuildElement("Fe");
	G4Element* Ni = nist->FindOrBuildElement("Ni");
	
	
	
	G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, nelements=6);
	StainlessSteel->AddElement(C, fractionmass=0.001);
	StainlessSteel->AddElement(Si, fractionmass=0.007);
	StainlessSteel->AddElement(Cr, fractionmass=0.18);
	StainlessSteel->AddElement(Mn, fractionmass=0.01);
	StainlessSteel->AddElement(Fe, fractionmass=0.712);
	StainlessSteel->AddElement(Ni, fractionmass=0.09);
*/
 //G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

G4VPhysicalVolume* G1DetectorConstruction::DefineVolumes(){

//World
// fWorld_mat = G4Material::GetMaterial("Air");

//Bc_400
 //fBc_400_mat = G4Material::GetMaterial("BC-400");



//Aluminum
//G4Material* alu_mat  = G4Material::GetMaterial("Aluminum");


//BC-400 property function of photon energy
 G4double photonEnergy[] =
            {3.1*eV,2.95*eV,2.93*eV,2.81*eV,2.7*eV,2.58*eV,2.48*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);


  G4double refractiveIndex1[] =
            {1.58,1.58,1.58,1.58,1.58,1.58,1.58 };

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
           {250*cm,250*cm,250*cm,250*cm,250*cm,250*cm,250*cm };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
            { 0.09,0.9,1.0,0.7,0.37,0.2,0.15 };
            
 // G4double scintilSlow[] =
   //         { 0.,0.,0.,0.,0.,0.,0. };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  

  

  G4MaterialPropertiesTable* scin_mat_prob = new G4MaterialPropertiesTable();

  scin_mat_prob->AddProperty("RINDEX",  photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
  scin_mat_prob->AddProperty("ABSLENGTH", photonEnergy, absorption,     nEntries)
        ->SetSpline(true);
  scin_mat_prob->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
        ->SetSpline(true);
  //scin_mat_prob->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
    //    ->SetSpline(true);



//Scintillator const property
  scin_mat_prob->AddConstProperty("SCINTILLATIONYIELD",fScintillationYield/MeV );
  scin_mat_prob->AddConstProperty("RESOLUTIONSCALE", fResolutionScale );
  scin_mat_prob->AddConstProperty("FASTTIMECONSTANT", fFastTimeConstant*ns );
 // scin_mat_prob->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  scin_mat_prob->AddConstProperty("YIELDRATIO", fYieldRatio);

  
  //Mie Scattering
  G4double bc400_Energy[] = {2.9*eV,3.0*eV,3.2*eV,3.7*eV,4.2*eV,5.0*eV,5.2*eV};
  G4double mie_BC400[] = {300*mm,300*mm,300*mm,300*mm,300*mm,300*mm,300*mm}; //küçültürsen bu process de devreye giriyor

const G4int nEntriesBC400 = sizeof(bc400_Energy)/sizeof(G4double);




  assert(sizeof(mie_BC400) == sizeof(bc400_Energy));
  
  G4double mie_BC400_const[3]={0.99,0.99,0.8};
  
  //spline olunca spline şekilde interpolasyon yapo.liner değil.forumdan bak
  scin_mat_prob->AddProperty("MIEHG",bc400_Energy,mie_BC400,nEntriesBC400)
        ->SetSpline(true);
  scin_mat_prob->AddConstProperty("MIEHG_FORWARD",mie_BC400_const[0]);
  scin_mat_prob->AddConstProperty("MIEHG_BACKWARD",mie_BC400_const[1]);
  scin_mat_prob->AddConstProperty("MIEHG_FORWARD_RATIO",mie_BC400_const[2]);

  //G4cout << "BC400 G4MaterialPropertiesTable" << G4endl;
  scin_mat_prob->DumpTable();


  fBc_400_mat->SetMaterialPropertiesTable(scin_mat_prob);
  
  //bunu öğren
  // Set the Birks Constant for the bc400 scintillator

  fBc_400_mat->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
  
  //Air
  
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
             

  G4MaterialPropertiesTable* air_mat_prob = new G4MaterialPropertiesTable();
  air_mat_prob->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  //G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  air_mat_prob->DumpTable();

  fWorld_mat->SetMaterialPropertiesTable(air_mat_prob);
  
  
  

// ------------- Volumes --------------
  
  
  
 //BC-400 PlasticScintillator


//World
 fWorld_solid = new G4Box("World", 0.5*fWorldSize_x, 0.5*fWorldSize_y, 0.5*fWorldSize_z);
 fWorld_logic= new G4LogicalVolume(fWorld_solid, fWorld_mat, "World", 0, 0, 0);

   fPhysWorld
    = new G4PVPlacement(0, G4ThreeVector(), fWorld_logic, "World", 0, false, checkOverlaps);




//Bc_400 solid,logic
 fBc_400_solid = new G4Box("BC-400", 0.5*fBc_400_size_x, 0.5*fBc_400_size_y, 0.5*fBc_400_size_z);
 fBc_400_logic = new G4LogicalVolume(fBc_400_solid, fBc_400_mat, "BC-400", 0, 0, 0);



  
  
  //place BC400 in world


fPhysBc = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fBc_400_logic,         //its logical volume
                      "BC-400",               //its name
                      fWorld_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking  





//photocathode side
G4double physCat_thick_z = 1.0e-5*cm;
G4Box* virtual_solid1 = new G4Box("photocathodePV", 0.5*fBc_400_size_x, 0.5*fBc_400_size_y, 0.5*physCat_thick_z);
G4LogicalVolume* virtual_logic1 = new G4LogicalVolume(virtual_solid1, fWorld_mat, "photocathodePV",0,0,0);

fPhysPmt =
 new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, -(0.5*fBc_400_size_z + 0.5*physCat_thick_z)),      //at (0,0,0)
                      virtual_logic1,         //its logical volume
                      "photocathodePV",               //its name
                      fWorld_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking  




//Aluminium side
G4double physAlu_thick_z = 1.0e-5*cm;
G4Box* virtual_solid2 = new G4Box("aluPV", 0.5*fBc_400_size_x, 0.5*fBc_400_size_y, 0.5*physAlu_thick_z);
G4LogicalVolume* virtual_logic2 = new G4LogicalVolume(virtual_solid2, fWorld_mat,"aluPV",0,0,0);

fPhysAlu =
 new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, 0.5*(fBc_400_size_z + physAlu_thick_z)),      //at (0,0,0)
                      virtual_logic2,         //its logical volume
                      "aluPV",               //its name
                      fWorld_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking  




//-------------------------SURFACES-----------------------------


//NOT:::the unified model applies only to dielectric_dielectric surface types


//Bc_400-air surface(cilt oluşturulacak)

  
  
  
  G4OpticalSurface* opBcSurface = new G4OpticalSurface("BC400-surface");
  opBcSurface->SetType(dielectric_dielectric);
  opBcSurface->SetModel(unified);
  opBcSurface->SetFinish( fSurfaceFinish );
  
  if( fSurfaceFinish == 0 ) G4cout<<"semele"<<G4endl;
  //burası sintilator-air yüzeyi için.ön taraf yani
  if(fSurfaceFinish == 3 || fSurfaceFinish == 2 || fSurfaceFinish == 5)
  opBcSurface->SetSigmaAlpha (0.1);


 
 
//hangi hacmi ciltliyosun.onun logic bilgisi burada fBc_400_logic
   fBc_skin_logic =
          new G4LogicalSkinSurface("BC400-surface", fBc_400_logic, opBcSurface);
          


  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
        (fBc_skin_logic->GetSurface(fBc_400_logic)->GetSurfaceProperty());

  if (opticalSurface) opticalSurface->DumpInfo();

 

 const G4int num = 2;
  G4double ephoton[num] = {2.9*eV, 5.2*eV};

//reflectivity refers to wrapping(arka yüzey), versus absorption rate,no refraction on painted surface,just reflection and absorption,wrapiing is smooth mirror(polish) only SS reflection.o yüzden burada(arka yüzey) reflection type belirtmeye gerek yok.zorunlu SS reflection
  
  G4double reflectivity[num] = {fReflectivity, fReflectivity};
  
  //yüzde 20 absorplandı fakat bunların kaçı detekte edildi onun olalığı...
  G4double efficiency[num]   = {0.5, 0.5};
  
 
  
  if( fSurfaceFinish == polishedbackpainted || fSurfaceFinish == groundbackpainted )  // PBP(2), GBP(5)
  G4double rindexforgap[num]   = {1.0, 1.0};  //refractivity refers to gap material, here air



if(fSurfaceFinish == 3 || fSurfaceFinish == 2 || fSurfaceFinish == 5) // G(3), PBP(2), GBP(5)
{

// bu kısım ön yüzey için (2,5 için).sintilator hava arası..reflectivity type refers to scintillator-air surface--not wrapping!!..
	G4double specularlobe[num] = {0.1, 0.1};
	G4double specularspike[num] = {0.7, 0.7};
	G4double backscatter[num] = {0.1, 0.1};
        //G4double ss[num] = {1.0,1.0};
}


  G4MaterialPropertiesTable *bc_mat_prob = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

bc_mat_prob->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
bc_mat_prob->AddProperty("EFFICIENCY", ephoton, efficiency, num);

if(fSurfaceFinish == 2 || fSurfaceFinish == 5) // PBP(2), GBP(5)
bc_mat_prob->AddProperty("RINDEX", ephoton, rindexforgap, num);

if(fSurfaceFinish == 3 || fSurfaceFinish == 2 || fSurfaceFinish == 5) // G(3), PBP(2), GBP(5)
{
bc_mat_prob->AddProperty("SPECULARLOBECONSTANT", ephoton, specularlobe, num);
bc_mat_prob->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularspike, num);
bc_mat_prob->AddProperty("BACKSCATTERCONSTANT", ephoton, backscatter, num);
//Lambertial 1-(diğer olasılıklar)..o yüzden belirtmeye gerek yok
}

 // G4cout << "bc_mat_prob G4MaterialPropertiesTable" << G4endl;
  bc_mat_prob->DumpTable();

  opBcSurface->SetMaterialPropertiesTable(bc_mat_prob);








//NOT:::the unified model applies only to dielectric_dielectric surface types

//Bc_400-Aluminium surface

G4OpticalSurface* opAluSurface = new G4OpticalSurface("Aluminum-surface");
  opAluSurface->SetType(dielectric_metal);
  opAluSurface->SetModel(glisur);
  opAluSurface->SetFinish(polished);


G4double alu_reflectivity[num] = {1.0, 1.0};
  
  //absorplama yok dolayısıyla efficicy verilen değer önemsiz
  G4double alu_efficiency[num]   = {0.0, 0.0};

   fAlu_border_logic = new G4LogicalBorderSurface("Aluminum-surface",fPhysBc,fPhysAlu,opAluSurface); 


//yüzeyin materyal özellikleri
 G4MaterialPropertiesTable *alu_mat_prop = new G4MaterialPropertiesTable();

alu_mat_prop->AddProperty("REFLECTIVITY", ephoton, alu_reflectivity, num);
alu_mat_prop->AddProperty("EFFICIENCY", ephoton, alu_efficiency, num);

//G4cout << "alu_mat_prop G4MaterialPropertiesTable" << G4endl;
  alu_mat_prop->DumpTable();

opAluSurface->SetMaterialPropertiesTable(alu_mat_prop);




//NOT:::the unified model applies only to dielectric_dielectric surface types
//Bc_400-PhotoCathode  surface



G4OpticalSurface* opPmtSurface = new G4OpticalSurface("Pmt-surface");
  opPmtSurface->SetType(dielectric_metal);
  opPmtSurface->SetModel(glisur);
  opPmtSurface->SetFinish(polished ); //so specular spike mandotary


//Hepsini absorplasın
G4double pmt_reflectivity[num] = {0.0, 0.0};
  
  
//  Absorplananların hepsi detekde edilsin verim yüzde yüz
  G4double pmt_efficiency[num]   = {1.0, 1.0};

  fPmt_border_logic = new G4LogicalBorderSurface("Pmt-surface",fPhysBc,fPhysPmt,opPmtSurface); 


//no refraction for metal--surface.surface finish is polish. so specular spike mandotary 
 G4MaterialPropertiesTable *pmt_mat_prop = new G4MaterialPropertiesTable();

pmt_mat_prop->AddProperty("REFLECTIVITY", ephoton, pmt_reflectivity, num);
pmt_mat_prop->AddProperty("EFFICIENCY", ephoton, pmt_efficiency, num);

G4cout << "pmt_mat_prob G4MaterialPropertiesTable" << G4endl;
  pmt_mat_prop->DumpTable();

opPmtSurface->SetMaterialPropertiesTable(pmt_mat_prop);

/*           
          G4Colour  white   ()              ;  // white
     G4Colour  white   (1.0, 1.0, 1.0) ;  // white
     G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
     G4Colour  black   (0.0, 0.0, 0.0) ;  // black
     G4Colour  red     (1.0, 0.0, 0.0) ;  // red
     G4Colour  green   (0.0, 1.0, 0.0) ;  // green
     G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
     G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
     G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
     G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow 
           
 */    





//VisAttributes of Bc_400

     G4VisAttributes * visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
          // Set the forced wireframe style 
     visAttributes->SetForceSolid(false);
          // Assignment of the visualization attributes to the logical volume
     fBc_400_logic->SetVisAttributes(visAttributes);
     fVisAttributes.push_back(visAttributes);



     

  
//VisAttributes of virtual_logic1 and virtual_logic2
visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
visAttributes->SetForceSolid(true);
virtual_logic1->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes);     
           
    
     visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));      
           visAttributes->SetForceSolid(true);
              
virtual_logic2->SetVisAttributes(visAttributes);
fVisAttributes.push_back(visAttributes); 


//VisAttributes of patient bed
visAttributes = new G4VisAttributes();
visAttributes->SetVisibility(false);
fWorld_logic->SetVisAttributes(visAttributes);
     fVisAttributes.push_back(visAttributes);


Print(); //buraya koyulmalı cünkü segmen fault verir. tüm değerlerin atanmış olması gerek cünkü..
return fPhysWorld;


}





void G1DetectorConstruction::SetReflectivity( G4double reflectivity )
{

this->fReflectivity = reflectivity;

G4RunManager::GetRunManager()->ReinitializeGeometry(); //gerekli yoksa aktif olmuyor

}

void G1DetectorConstruction::SetResolutionScale(G4double rs)
{

this->fResolutionScale = rs;

G4RunManager::GetRunManager()->ReinitializeGeometry();

}

void G1DetectorConstruction::SetDetectorDimension(G4ThreeVector tv)
{

this->fBc_400_size_x = tv[0];
this->fBc_400_size_y = tv[1];
this->fBc_400_size_z = tv[2];

G4RunManager::GetRunManager()->ReinitializeGeometry();  //gerekli yoksa aktif olmuyor

}

void G1DetectorConstruction::SetDefaults(){


fBc_400_size_x = 76*cm;
fBc_400_size_y = 25*cm;
fBc_400_size_z = 91.5*cm;


//World Size
fWorldSize_x = 300*cm;
fWorldSize_y = 300*cm;
fWorldSize_z = 300*cm;


//Scintillator const property
fScintillationYield = 10.0;
fResolutionScale = 1.0;
fFastTimeConstant = 2.4;
fYieldRatio = 1.0;

fSurfaceFinish = polishedbackpainted; // P(0)



fReflectivity = 0.8;

}

void G1DetectorConstruction::Print()
{
G4cout<<"----------------- Detector Property----------------"<<G4endl;
G4cout<< "Detector Dimension: "<<"x: "<<G4BestUnit(fBc_400_size_x,"Length")<<" y:"<<G4BestUnit(fBc_400_size_y, "Length")<<" z:"<<
G4BestUnit(fBc_400_size_z,"Length")<<G4endl;

G4cout<< "Detector Material: "<<fBc_400_mat->GetName()<<G4endl;
G4cout<< "Detector Mass: "<<G4BestUnit(fBc_400_logic->GetMass(),"Mass")<<G4endl;
G4cout<<"Bc-400 wrapping material reflectivity  coeff: "<<fReflectivity<<G4endl;
G4cout<<"Resolution Scale: "<<fResolutionScale<<G4endl;
G4cout<<"Surface Finish Type: "<<G4endl;

G4cout<<"-------------------------------------------"<<G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
