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
#include <fstream>

#include <sstream>
#include <iomanip>

#include "G1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4double N_A = (6.02e+23)*(1/mole);


G1DetectorConstruction::G1DetectorConstruction()
 :  G4VUserDetectorConstruction(), fWorld_mat(0), fScin_mat(0), fTrp_mat(0), fOptCement_mat(0),
                                   fWorld_solid(0), fScin_solid(0), fUnit_solid(0), fHall_solid(0), fLayer_solid(0),
                                   fWorld_logic(0), fScin_logic(0), fUnit_logic(0), fHall_logic(0), fLayer_logic(0), 
                                   fTrp_logic(0), fOptical_cement_logic(0), fOptical_cement_logic1(0),              
                                   fPhysWorld(0), fPhysScin(0), fPhysPmt(0), fPhysTrp(0), fPhysOpticalCement(0), 
                                   fPhysOpticalCement1(0) ,fPhysUnit(0), fPhysHall(0),
                                   fRepUnit(0), fRepLayer(0),
                                   checkOverlaps(true), fIsScinWrapped(false), fIsWrapHasTeoricReflectivity(false), fIsPmtHasTeoricQE(false), fVisAttributes()
                                   
                                   

{

  SetDefaults();
  fDetectorMessenger = new G1DetectorMessenger(this);
  fPhysTrp =new G4VPhysicalVolume*[2];
  fPhysPmt =new G4VPhysicalVolume*[2];
  
  fPhysOpticalCement =new G4VPhysicalVolume*[2]; //scin-lightGuide interface
   fPhysOpticalCement1 =new G4VPhysicalVolume*[2]; //lightGuide-pmt interface
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G1DetectorConstruction::~G1DetectorConstruction()
{

G4cout <<"Detector destruction..."<<G4endl;
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


void G1DetectorConstruction::DefineMaterials()
{

  G4double a, z, density;
  G4int nelements;
  
  G4double HydrogenAtomDensity, CarbonAtomDensity, OxygenAtomDensity;
  G4double HydrogenMassDensity, CarbonMassDensity, OxygenMassDensity;
  G4double HydrogenMassFraction,CarbonMassFraction, OxygenMassFraction;

   G4Element* H = new G4Element("Hydrogen", "H", z=1. , a=1.00794*g/mole);
   G4Element* C = new G4Element("Carbon", "C", z=6. , a=12.0107*g/mole);
   G4Element* N = new G4Element("Nitrogen", "N", z=7. , a=14.01*g/mole);
   G4Element* O = new G4Element("Oxygen"  , "O", z=8. , a=16.00*g/mole);
   
/*---------------------------------------------EJ-200 PlasticScintillator Material-----------------------------------------------*/
//Compounds are defined by mass fractions 


 HydrogenAtomDensity = (5.17e+22);
 CarbonAtomDensity = (4.69e+22);

HydrogenMassDensity = ConvertMassDensity(HydrogenAtomDensity,a=1.00794);
CarbonMassDensity = ConvertMassDensity(CarbonAtomDensity, a=12.0107);
//G4cout<<"Hyrogen Mass density: "<<HydrogenMassDensity/(g/cm*cm*cm)<<G4endl;

HydrogenMassFraction = HydrogenMassDensity/(HydrogenMassDensity+CarbonMassDensity);
CarbonMassFraction = CarbonMassDensity/(HydrogenMassDensity+CarbonMassDensity);
//G4cout <<"HydrogenMassFraction: "<<HydrogenMassFraction<<G4endl;
//G4cout<<"Scintilator density: "<<(HydrogenMassDensity+CarbonMassDensity)/(g/cm*cm*cm)<<G4endl; 1.023 olmalı
   
   fScin_mat = new G4Material("EJ-200", density=1.023*g/cm3, nelements=2);
   fScin_mat->AddElement(H,HydrogenMassFraction);
   fScin_mat->AddElement(C,CarbonMassFraction);

/*------------------------------------------------LightGuide Material---------------------------------------------------------*/
  
  HydrogenAtomDensity = 5.73e+22;
  CarbonAtomDensity = 3.58e+22;
  OxygenAtomDensity = 1.43e+22;
  
  HydrogenMassDensity = ConvertMassDensity(HydrogenAtomDensity,a=1.00794);
  CarbonMassDensity = ConvertMassDensity(CarbonAtomDensity,a=12.0107);
  OxygenMassDensity = ConvertMassDensity(OxygenAtomDensity, a=16.00);
  
  HydrogenMassFraction = HydrogenMassDensity/(HydrogenMassDensity+CarbonMassDensity+OxygenMassDensity);
  CarbonMassFraction = CarbonMassDensity/(HydrogenMassDensity+CarbonMassDensity+OxygenMassDensity);
  OxygenMassFraction = OxygenMassDensity/(HydrogenMassDensity+CarbonMassDensity+OxygenMassDensity);
  
  
  density = (HydrogenMassDensity+CarbonMassDensity+OxygenMassDensity);
  
  fTrp_mat = new G4Material("TrpM", density, nelements=3);
  fTrp_mat->AddElement(H, HydrogenMassFraction );
  fTrp_mat->AddElement(C, CarbonMassFraction);
  fTrp_mat->AddElement(O, OxygenMassFraction);
  
  //G4cout <<"LightGuide density:" <<density/(g/cm*cm*cm)<<G4endl;

/*------------------------OpticalCement Material---------------------------------------------------------------------------------*/

// (PolyVinylToluene, C_9H_10)
   G4NistManager* nistManager = G4NistManager::Instance();
  fOptCement_mat = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  





/*------------------------Air Material--------------------------------------------------------------------------------*/

  

  fWorld_mat = new G4Material("Air", density=1.29*g/cm3, nelements=2);
  fWorld_mat->AddElement(N, 70.*perCent);
  fWorld_mat->AddElement(O, 30.*perCent);

}



G4VPhysicalVolume* G1DetectorConstruction::DefineVolumes()
{

const G4int numberOfPhotons =51;
 
G4double photonWaveLength[numberOfPhotons];
G4double photonEnergy[numberOfPhotons];

G4double scintilFast[numberOfPhotons];
G4double scinRefractiveIndex[numberOfPhotons];
G4double absorption[numberOfPhotons];



ReadEmissionSpectrumFromFile("spectrum.txt", photonWaveLength, photonEnergy , scintilFast, scinRefractiveIndex, absorption);


// Add the material's optical properties

const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);


//  assert(sizeof(absorption) == sizeof(photonEnergy));
 // assert(sizeof(scintilFast) == sizeof(photonEnergy));

  //Scintillator property changes with photonEnergy
  G4MaterialPropertiesTable* scin_mat_prop = new G4MaterialPropertiesTable();

  scin_mat_prop->AddProperty("RINDEX",  photonEnergy, scinRefractiveIndex, nEntries)
        ->SetSpline(true);
  scin_mat_prop->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)
        ->SetSpline(true);
  scin_mat_prop->AddProperty("FASTCOMPONENT", photonEnergy, scintilFast, nEntries)
        ->SetSpline(true);
  //scin_mat_prop->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
    //    ->SetSpline(true);


//Scintillator const property
  scin_mat_prop->AddConstProperty("SCINTILLATIONYIELD", fScintillationYield/MeV );
  scin_mat_prop->AddConstProperty("RESOLUTIONSCALE", fResolutionScale );
  scin_mat_prop->AddConstProperty("FASTTIMECONSTANT", fFastTimeConstant*ns );
 // scin_mat_prop->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  scin_mat_prop->AddConstProperty("YIELDRATIO", fYieldRatio);

  /*------------------------------Mie Scattering---------------------------------------------------------------------------*/
  G4double EJ200_Energy[] = { 2.9*eV, 3.0*eV, 3.2*eV, 3.7*eV, 4.2*eV, 5.0*eV, 5.2*eV };
  G4double mie_EJ200[] = { 300*m, 300*m, 300*m, 300*m, 300*m, 300*m, 300*m }; //küçültürsen bu process de devreye giriyor

  const G4int nEntriesEJ200 = sizeof(EJ200_Energy)/sizeof(G4double);

 assert(sizeof(mie_EJ200) == sizeof(EJ200_Energy));
  
  G4double mie_EJ200_const[3]= { 0.99, 0.99, 0.8 };
  
  //spline olunca spline şekilde interpolasyon yapo.liner değil.forumdan bak
  scin_mat_prop->AddProperty("MIEHG",EJ200_Energy,mie_EJ200,nEntriesEJ200)
        ->SetSpline(true);
  scin_mat_prop->AddConstProperty("MIEHG_FORWARD",mie_EJ200_const[0]);
  scin_mat_prop->AddConstProperty("MIEHG_BACKWARD",mie_EJ200_const[1]);
  scin_mat_prop->AddConstProperty("MIEHG_FORWARD_RATIO",mie_EJ200_const[2]);

  G4cout<<"--------Scintillator G4MaterialPropertiesTable---------"<<G4endl;
  scin_mat_prop->DumpTable();
  

  fScin_mat->SetMaterialPropertiesTable(scin_mat_prop);
  
  //bunu öğren
  // Set the Birks Constant for the EJ200 scintillator
fScin_mat->GetIonisation()->SetBirksConstant( 0.126*mm/MeV );
  
  //fScin_mat->GetIonisation()->SetBirksConstant( 0.0126*g/MeV*cm*cm );
  
  /*------------------------------Air Refractive Index ---------------------------------------------------------------------------*/
  
   G4double airRefractiveIndex[numberOfPhotons];
  
   for(G4int i=0; i<numberOfPhotons; i++)
   airRefractiveIndex [i] =1.0002772;
   
  G4MaterialPropertiesTable* air_mat_prob = new G4MaterialPropertiesTable();
  air_mat_prob->AddProperty("RINDEX", photonEnergy, airRefractiveIndex, nEntries);

  G4cout << "-------Air G4MaterialPropertiesTable--------" << G4endl;
  air_mat_prob->DumpTable();
 
  fWorld_mat->SetMaterialPropertiesTable( air_mat_prob );
  
  
  /*---------------------------LightGuide Refractive Index------------------------------------------------------------------------------*/
  
  G4double trpRefractiveIndex[numberOfPhotons];
  
   for(G4int i=0; i<numberOfPhotons; i++)
   trpRefractiveIndex [i] =fTrpRefractiveIndex;
  
  
  G4MaterialPropertiesTable* trp_mat_prop = new G4MaterialPropertiesTable();
  trp_mat_prop->AddProperty("RINDEX", photonEnergy, trpRefractiveIndex, nEntries);
  trp_mat_prop->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true);
  G4cout << "-------LightGuide G4MaterialPropertiesTable--------" << G4endl;
  trp_mat_prop->DumpTable();
 
  fTrp_mat->SetMaterialPropertiesTable( trp_mat_prop );
 
 /*------------------------------OpticalCement Refractive Index---------------------------------------------------------------------------*/ 
  //OpticalCement refractive index
  G4double opticalCementRefractiveIndex[numberOfPhotons]; //between scintillator and lightguide
   
  for(G4int i=0; i<numberOfPhotons; i++)
   {
   opticalCementRefractiveIndex [i] =fCementRefractiveIndex;
   
   }
   
  G4MaterialPropertiesTable* optical_cement_mat_prop = new G4MaterialPropertiesTable();
  optical_cement_mat_prop->AddProperty("RINDEX", photonEnergy, opticalCementRefractiveIndex, nEntries);
 optical_cement_mat_prop->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true);
  G4cout << "-------OpticalCement G4MaterialPropertiesTable--------" << G4endl;
  optical_cement_mat_prop->DumpTable();
 
  fOptCement_mat->SetMaterialPropertiesTable( optical_cement_mat_prop );
  
  
  /*------------------------------SOLİD,LOGİC---------------------------------------------------------------------------*/
  
 //World
 fWorld_solid = new G4Box("World", 0.5*fWorldSize_x, 0.5*fWorldSize_y, 0.5*fWorldSize_z);
 fWorld_logic = new G4LogicalVolume(fWorld_solid, fWorld_mat, "World", 0, 0, 0);

   fPhysWorld = new G4PVPlacement(0, G4ThreeVector(), fWorld_logic, "World", 0, false, checkOverlaps);

//Unit

fUnit_solid = new G4Box("Unit", 0.5*fUnit_size_x, 0.5*fUnit_size_y, 0.5*fUnit_size_z);
fUnit_logic = new G4LogicalVolume(fUnit_solid, fWorld_mat, "Unit", 0, 0, 0);

//Layer
fLayer_solid= new G4Box("Layer", 0.5*fLayer_size_x, 0.5*fLayer_size_y, 0.5*fLayer_size_z);
fLayer_logic = new G4LogicalVolume(fLayer_solid, fWorld_mat, "Layer", 0, 0, 0);

//Hall
fHall_solid = new G4Box("Hall", 0.5*fHall_size_x, 0.5*fHall_size_y, 0.5*fHall_size_z);
fHall_logic = new G4LogicalVolume(fHall_solid, fWorld_mat, "Hall", 0, 0, 0);

//Ej_200 
 fScin_solid = new G4Box("EJ-200", 0.5*fScin_size_x, 0.5*fScin_size_y, 0.5*fScin_size_z);
 fScin_logic = new G4LogicalVolume(fScin_solid, fScin_mat , "EJ-200", 0, 0, 0);

 //Light guide ,pmt bağlantı kısmı kare olacak..dx1=dy1
G4double  dx1 = fTrp_Length_x/2; //bunlar pmt ağzına göre ayarlancak
G4double  dx2 = fScin_size_x/2.;
G4double  dy1 = fTrp_Length_y/2; //bunlar pmt ağzına göre ayarlancak
G4double  dy2 = fScin_size_y/2.;

G4Trd *trp_solid = new G4Trd ("LightGuide", dx1, dx2, dy1, dy2, 0.5*fTrp_Length_z); 
fTrp_logic = new G4LogicalVolume(trp_solid, fTrp_mat, "LightGuide",0,0,0);


//OpticalCement between scintillator and lightGuide
G4Box *optical_cement_solid = new G4Box ("OpticalCement", dx2, dy2, 0.5*fOptical_cement_size_z); 
fOptical_cement_logic = new G4LogicalVolume(optical_cement_solid, fOptCement_mat, "OpticalCement",0,0,0);
 
// OpticalCement between lightGuide and PMT

G4Tubs *optical_cement_solid1 = new G4Tubs ("OpticalCement1", fPmt_innerRadius, fPmt_outerRadius, 0.5*fOptical_cement_size_z, 0*deg, 360*deg); 
fOptical_cement_logic1 = new G4LogicalVolume (optical_cement_solid1, fOptCement_mat, "OpticalCement1",0,0,0); 
 
 
 //PMT 
 
G4Tubs* pmt_solid = new G4Tubs("Pmt", fPmt_innerRadius, fPmt_outerRadius, fPmt_Length_z/2., 0*deg, 360*deg);
  G4LogicalVolume* pmt_logic = new G4LogicalVolume(pmt_solid, fWorld_mat, "Pmt",0,0,0);
 
 
 
 /*------------------------------PLACEMENT---------------------------------------------------------------------------*/
 
 //place EJ-200 in unit
 fPhysScin = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,0),       //at (0,0,0)
                      fScin_logic,         //its logical volume
                      "EJ-200",               //its name
                      fUnit_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking  
                      
 //place pmt and ligt guides and optical cement in unit    
  
  G4int sign=-1;   
    
               for(int i=0; i<2; i++){   
      
        
            if(i==1) { //pozitif z tarafı
            sign=1; 
            G4RotationMatrix *rott = new G4RotationMatrix();
            rott->rotateX(180*deg);
            
            
              //place optical cement between scintillator and LightGuide
                      fPhysOpticalCement[i] = new G4PVPlacement(rott,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + 0.5*fOptical_cement_size_z)),      //at (0,0,0)
                      fOptical_cement_logic,         //its logical volume
                      "OpticalCement",               //its name
                      fUnit_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      i+1,                     //copy number
                      checkOverlaps);        //overlaps checking 
               
            
            //place trapezoid
            fPhysTrp[i] = new G4PVPlacement(rott,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + fOptical_cement_size_z + 0.5*fTrp_Length_z)),      //at (0,0,0)
                      fTrp_logic,         //its logical volume
                      "LightGuide",               //its name
                      fUnit_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      i,                     //copy number
                      checkOverlaps);        //overlaps checking 
             
            
            //place optical cement between LightGuide and Pmt
                      fPhysOpticalCement1[i] = new G4PVPlacement(rott,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + fOptical_cement_size_z + fTrp_Length_z +0.5*fOptical_cement_size_z)),    
                      fOptical_cement_logic1,         //its logical volume
                      "OpticalCement1",               //its name
                      fUnit_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      i+2,                     //copy number
                      checkOverlaps);        //overlaps checking 
            
            
             }   else
             {
             //negatif z tarafı
             
             //place optical cement between scintillator and LightGuide
             
                      fPhysOpticalCement[i] = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + 0.5*fOptical_cement_size_z)),      //at (0,0,0)
                      fOptical_cement_logic,         
                      "OpticalCement",               
                      fUnit_logic,                     
                      false,                 
                      i,                    
                      checkOverlaps);        
                  
                  
             //place LightGuide       
             fPhysTrp[i] = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + fOptical_cement_size_z + 0.5*fTrp_Length_z)),      //at (0,0,0)
                      fTrp_logic,         
                      "LightGuide",               
                      fUnit_logic,                     
                      false,                 
                      i,                     
                      checkOverlaps);        
                      
                      
                      //place optical cement between LightGuide and Pmt
            
            fPhysOpticalCement1[i] = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + fOptical_cement_size_z + fTrp_Length_z +0.5*fOptical_cement_size_z)),      //at (0,0,0)
                      fOptical_cement_logic1,         
                      "OpticalCement1",               
                      fUnit_logic,                     
                      false,                
                      i+1,                     
                      checkOverlaps);        
             }
                 
  
//place pmt
  fPhysPmt[i] = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, sign*(0.5*fScin_size_z + fOptical_cement_size_z + fTrp_Length_z + fOptical_cement_size_z + 0.5*fPmt_Length_z  )),      //at (0,0,0)
                      pmt_logic,         //its logical volume
                      "Pmt",               //its name
                      fUnit_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      i,                     //copy number
                      checkOverlaps);        //overlaps checking  


}

/*------------------------------Unit,Layer,Hall---------------------------------------------------------------------------*/
    //layer consist of units.index number for each unit   
 fRepUnit = new G4PVReplica("Unit", fUnit_logic, fLayer_logic, kXAxis, numberOfUnitInLayerX,fUnit_size_x );
 
 //Hall consist of layers.index number for layer
 fRepLayer = new G4PVReplica("Layer", fLayer_logic, fHall_logic, kYAxis, numberOfLayerAlongY,fUnit_size_y );


//place hall inside world
 new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,0),       //at (0,0,0)
                      fHall_logic,         //its logical volume
                      "Hall",               //its name
                      fWorld_logic,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking  



/*------------------------------SURFACES---------------------------------------------------------------------------*/

//if surfacefinish is polish and ground this means no wrapping.
//NOT:::the unified model applies only to dielectric_dielectric surface types

 //Scintillator yüzeyi,ciltleme için ->bütün yüzeyler kaplanıyor
 
  
  switch( fSurfaceFinish ){
   
    case 0 :
    SetP();
    break;
    
    case 1 :
    SetFP();
    break;
    
    case 2 :
    SetBP();
    //G4cout <<"case: 2"<<G4endl;
    break;
    
    case 3 :
    SetG();
    //G4cout <<"case: 3"<<G4endl;
    break;
    
    case 4 :
    SetFP();
    break;
    
    case 5 :
    SetBP();
    break;
    
    //polish
    default :
    break; 
  
   }  

/*------------------------------SURFACES INSIDE UNIT(6 surface)---------------------------------------------------------------------------*/
   //normalde scintillator, lightGuide ve cementin tüm yüzeyleri ciltlenmiş.iç yüzeyler için yeni bir yüzey tanımlayarak iç yüzeylerin ciltleme işini kaldırıyouz.toplam 6 yüzey..İç yüzeyler kaplanamayacağı için ya ground yada polish olmalı.
   
   
  G4OpticalSurface* OpInsideSurface = new G4OpticalSurface("EjTrp-surface"); 
  OpInsideSurface->SetType(dielectric_dielectric);
  OpInsideSurface->SetModel(unified);
  OpInsideSurface->SetFinish(fInsideSurfaceFinish ); 
   
   //dizi indisi sıfır ise negatif z tarafı, 1 ise pozitif z tarafı
   //negatif
  
          
          
          
   new G4LogicalBorderSurface("ScinNCement", fPhysScin , fPhysOpticalCement[0], OpInsideSurface ); 
   new G4LogicalBorderSurface("ScinNCement", fPhysOpticalCement[0], fPhysScin, OpInsideSurface ); //reverse side
   
   new G4LogicalBorderSurface("CementNTrp", fPhysOpticalCement[0] , fPhysTrp[0], OpInsideSurface );
   new G4LogicalBorderSurface("CementNTrp", fPhysTrp[0], fPhysOpticalCement[0] ,OpInsideSurface ); //reverse side
   
   //pozitif
   new G4LogicalBorderSurface("ScinPCement", fPhysScin , fPhysOpticalCement[1], OpInsideSurface ); 
   new G4LogicalBorderSurface("ScinPCement",  fPhysOpticalCement[1], fPhysScin , OpInsideSurface ); //reverse side
   
   new G4LogicalBorderSurface("CementPTrp", fPhysOpticalCement[1] , fPhysTrp[1], OpInsideSurface );
   new G4LogicalBorderSurface("CementPTrp", fPhysTrp[1],  fPhysOpticalCement[1] , OpInsideSurface ); //reverse side
   
   
   new G4LogicalBorderSurface("TrpNCement", fPhysTrp[0] , fPhysOpticalCement1[0], OpInsideSurface ); //negatif
   new G4LogicalBorderSurface("TrpNCement", fPhysOpticalCement1[0], fPhysTrp[0] ,  OpInsideSurface ); //reverse side
   
   new G4LogicalBorderSurface("TrpPCement", fPhysTrp[1] , fPhysOpticalCement1[1], OpInsideSurface ); //pozitif
   new G4LogicalBorderSurface("TrpPCement", fPhysOpticalCement1[1], fPhysTrp[1] ,  OpInsideSurface ); //pozitif
   
   
   
   if(fInsideSurfaceFinish == 3){
          
         SetInsideSurfaceGround(OpInsideSurface);
          G4cout <<"OLDUUUUUUU"<<G4endl;
          }
  
  
  
  /*------------------------------PMT SURFACE---------------------------------------------------------------------------*/
  G4OpticalSurface* opPmtSurface = new G4OpticalSurface("Pmt-surface");
  opPmtSurface->SetType(dielectric_metal);
  opPmtSurface->SetModel(glisur);
  opPmtSurface->SetFinish(polished ); 


   
  new G4LogicalBorderSurface("PmtNSurace", fPhysOpticalCement1[0], fPhysPmt[0] , opPmtSurface ); //negatif
  new G4LogicalBorderSurface("PmtPSurface", fPhysOpticalCement1[1], fPhysPmt[1], opPmtSurface ); //pozitif 
   

const G4int numberOfPhotonPmt = 36;

G4double pmt_pwl[numberOfPhotonPmt];
G4double pmt_pe[numberOfPhotonPmt];
G4double pmt_efficiency[numberOfPhotonPmt]; //QUANTUM EFFICIENCY
G4double pmt_reflectivity[numberOfPhotonPmt];

/*
--if a photon is absorbed but not detected it is counted for absorption process. if a photon absorbed and detected it is counted for detection process
--alu and dielectric_dielectric surface detection efficiency is equal to zero, so all detection process arise from pmt surface
--Absorplananların hepsi detekde edilsin verim yüzde yüz.. //all contribute for detection process
*/

//Read QE spectrum from file and get the efficiency value for each photon
  ReadPmtEfficiencyFromFile("efficiency.txt", pmt_pwl, pmt_pe, pmt_efficiency, pmt_reflectivity);
 
  

//no refraction for metal--surface.surface finish is polish. so specular spike mandotary 
 G4MaterialPropertiesTable *pmt_mat_prop = new G4MaterialPropertiesTable();

pmt_mat_prop->AddProperty("REFLECTIVITY", pmt_pe, pmt_reflectivity, numberOfPhotonPmt);
pmt_mat_prop->AddProperty("EFFICIENCY", pmt_pe, pmt_efficiency, numberOfPhotonPmt);

G4cout << "---------Pmt surface material property--------- " << G4endl;
  pmt_mat_prop->DumpTable();

opPmtSurface->SetMaterialPropertiesTable(pmt_mat_prop);

/*---------------------------------------VisAttributes------------------------------------------------------------------*/

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




//VisAttributes of Scintilator and
     G4VisAttributes * visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
          // Set the forced wireframe style 
    visAttributes->SetForceSolid(false);
          // Assignment of the visualization attributes to the logical volume
     fScin_logic->SetVisAttributes(visAttributes);
     fVisAttributes.push_back(visAttributes);



 //VisAttributes of trp_logic and
visAttributes = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
visAttributes->SetForceSolid(false);
fTrp_logic->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes);     
           
    
 //VisAttributes of opticalcement between scintillator and lightGuide 
visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); //magenta
visAttributes->SetForceSolid(false);

fOptical_cement_logic->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes);    
     
//VisAttributes of opticalcement between lightGuide and Pmt
visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));  // magenta ));
visAttributes->SetForceSolid(false);
fOptical_cement_logic1->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes);  


//visAttributes of pmt_logic

visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));      
           visAttributes->SetForceSolid(true);
              
pmt_logic->SetVisAttributes(visAttributes);
fVisAttributes.push_back(visAttributes); 



//VisAttributes of unit
visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
visAttributes->SetVisibility(false);

fUnit_logic->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes); 


//VisAttributes of layer
visAttributes = new G4VisAttributes();
visAttributes->SetVisibility(false);
fLayer_logic->SetVisAttributes(visAttributes);
//virtual_logic2->SetVisAttributes(visAttributes);
   fVisAttributes.push_back(visAttributes); 


//VisAttributes of Hall
visAttributes = new G4VisAttributes();
visAttributes->SetVisibility(false);

fHall_logic->SetVisAttributes(visAttributes);
     fVisAttributes.push_back(visAttributes);



//VisAttributes of World
visAttributes = new G4VisAttributes();
visAttributes->SetVisibility(false);
fWorld_logic->SetVisAttributes(visAttributes);
     fVisAttributes.push_back(visAttributes);



return fPhysWorld;


}


void G1DetectorConstruction::SetP()
{
/*

Snell's law is applied based on refractive index of 2 media. But before this reflectivity is used to determine whether photon is absorbed (Strictly speaking here reflectivity is not reflection coefficient. it is one minus absorbtion coefficient.  )

*/

//No wrapping
  G4OpticalSurface* ScinSurface = new G4OpticalSurface("OpPolishSurface");
  ScinSurface->SetType(dielectric_dielectric);
  ScinSurface->SetModel(unified);
  ScinSurface->SetFinish(fSurfaceFinish);

  
    new G4LogicalSkinSurface("ScintillatorSkinSurface", fScin_logic, ScinSurface);
    new G4LogicalSkinSurface("LightGuideSkinSurface", fTrp_logic, ScinSurface);
    new G4LogicalSkinSurface("OpticalCementSkinSurface", fOptical_cement_logic, ScinSurface); 
    new G4LogicalSkinSurface("OpticalCement1SkinSurface", fOptical_cement_logic1, ScinSurface);
         
/*
G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
        (Scin_skin_logic->GetSurface(fScin_logic)->GetSurfaceProperty());

 
  if (opticalSurface) 
  {
  G4cout<<"------Optical surface DumpInfo---------"<<G4endl;
  opticalSurface->DumpInfo();
  }

*/
  const G4int num = 2;
  G4double ephoton[num] = { 2.901*eV, 2.902*eV };

  //empirical reflectivity...default value 1.yüzeyde kirlenmeden dolayı(kirlenme madddesi absorbe yapabilir) absorbe olması ihtimaline karşılık.
  
  G4double reflectivity[num];
  G4double efficiency[num];
  
  for(G4int i=0; i<num; i++){
 
  reflectivity[i] = fReflectivity;
  efficiency[i] = fEfficiency;
  }
  


 G4MaterialPropertiesTable *ScinSurfaceMatProp = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

ScinSurfaceMatProp->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
ScinSurfaceMatProp->AddProperty("EFFICIENCY", ephoton, efficiency, num);

G4cout << "-------Ej_200 Polish surface material property--------" << G4endl;
ScinSurfaceMatProp->DumpTable();

ScinSurface->SetMaterialPropertiesTable(ScinSurfaceMatProp);


}



void G1DetectorConstruction::SetG()
{

//No  wrapping, rough surface
 G4OpticalSurface* ScinSurface = new G4OpticalSurface("OpGroundSurface");
  ScinSurface->SetType(dielectric_dielectric);
  ScinSurface->SetModel(unified);
  ScinSurface->SetFinish(fSurfaceFinish);


  ScinSurface->SetSigmaAlpha (fSigmaAlpha);

//hangi hacmi ciltliyosun.onun logic bilgisi burada fScin_logic
 
 
     new G4LogicalSkinSurface("ScintillatorSkinSurface", fScin_logic, ScinSurface);
     new G4LogicalSkinSurface("LightGuideSkinSurface", fTrp_logic, ScinSurface);
     new G4LogicalSkinSurface("OpticalCementSkinSurface", fOptical_cement_logic, ScinSurface); 
     new G4LogicalSkinSurface("OpticalCement1SkinSurface", fOptical_cement_logic1, ScinSurface);
 
   
       
   const G4int num = 2;
  G4double ephoton[num] = { 2.901*eV, 2.902*eV }; //en az iki tane lazım
  
  G4double reflectivity[num];
  G4double efficiency[num];
  
  G4double specularspike[num];
  G4double specularlobe[num];
  G4double backscatter[num];
  
 G4double SsRefProb = 1.0;
G4double SlRefProb = 0.0;
G4double BsRefProb = 0.0;
  
 for(G4int i=0; i<num; i++){
 
  reflectivity[i] = fReflectivity;
  efficiency[i] = fEfficiency;
  
  specularspike[i] = SsRefProb;
  specularlobe[i] = SlRefProb;
  backscatter[i] = BsRefProb;
  
  }

G4MaterialPropertiesTable *ScinSurfaceMatProp = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

ScinSurfaceMatProp->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
ScinSurfaceMatProp->AddProperty("EFFICIENCY", ephoton, efficiency, num);
ScinSurfaceMatProp->AddProperty("SPECULARLOBECONSTANT", ephoton, specularlobe, num);
ScinSurfaceMatProp->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularspike, num);
ScinSurfaceMatProp->AddProperty("BACKSCATTERCONSTANT", ephoton, backscatter, num);
//Lambertial 1-(diğer olasılıklar)..o yüzden belirtmeye gerek yok

G4cout << "-------Ej_200 Ground surface material property--------" << G4endl;
ScinSurfaceMatProp->DumpTable();

ScinSurface->SetMaterialPropertiesTable(ScinSurfaceMatProp);




}


void G1DetectorConstruction::SetBP() //for GBP and PBP
{

/*

BC-620/Ej-510 REflector paint for plastic scintillator.Reflectivity grafiğinden 400-600 nm aralığındaki fotonlar için reflectivity yaklaşık 0.96.yani sabit

Diffrence:
PBP----The Polish refers to wrapping. it implies wrapping is a perfectly smooth mirror with only specular spike reflection taking place.
GBP----The Ground refers to wrapping. it implies wrapping is a ground mirror with only lambertian reflection taking place.

Common for GBP and PBP:

The fSigmaAlpha value specified refers to the scintillator-air (ön yüzey) gap interface. Snell laws is apllied after sampling the facet normal, and if reflection takes place, one of the four (specularlobe,specularspike,backscatter,lambertian) takes place with respect ro facetnormal according to assigned probabilities.

*/
  
  //plastic scintillator and the light guides are wrapped in aluminized Mylar
  
  G4OpticalSurface* wrappingSurface = new G4OpticalSurface("WrappingSurface");
  wrappingSurface->SetType(dielectric_dielectric);
  wrappingSurface->SetModel(unified);
  wrappingSurface->SetFinish(fSurfaceFinish);
  
  //burası sintilator-air yüzeyi için.ön taraf yani
  wrappingSurface->SetSigmaAlpha (fSigmaAlpha);


 //hangi hacmi ciltliyosun.onun logic bilgisi burada fScin_logic
  
   
    
     new G4LogicalSkinSurface("ScintillatorSkinSurface", fScin_logic, wrappingSurface);
     new G4LogicalSkinSurface("LightGuideSkinSurface", fTrp_logic, wrappingSurface);
     new G4LogicalSkinSurface("OpticalCementSkinSurface", fOptical_cement_logic, wrappingSurface); 
     new G4LogicalSkinSurface("OpticalCement1SkinSurface", fOptical_cement_logic1, wrappingSurface);

  const G4int numberOfPhotonWrap = 23;

G4double wrap_pwl[numberOfPhotonWrap];
G4double wrap_pe[numberOfPhotonWrap];
G4double wrap_reflectivity[numberOfPhotonWrap]; //Reflection probability versus absorbtion.no refraction on painted surface ,wrapping is smooth mirror(polish) or only SS reflection.o yüzden burada(arka yüzey) reflection type belirtmeye gerek yok.zorunlu SS reflection
G4double wrap_efficiency[numberOfPhotonWrap];


// bu kısım ön yüzey için (2,5 için).sintilator hava arası..reflectivity type refers to scintillator-air surface--not wrapping!!..
G4double wrap_specularSpike[numberOfPhotonWrap];
G4double wrap_specularLobe[numberOfPhotonWrap];
G4double wrap_backScatter[numberOfPhotonWrap];

//G4double wrap_Lambertian[numberOfPhotonWrap];

G4double refractiveIndexForGap[numberOfPhotonWrap];

G4double refIndexForGap= 1.0002772;

G4double SsRefProb = 0.0;
G4double SlRefProb = 1.0;;
G4double BsRefProb = 0.0;
//G4double LamRefProb = 1-(SsRefProb+SlRefProb+BsRefProb);

ReadWrapReflectivityFromFile("wrappingReflectivity.txt", wrap_pwl, wrap_pe, wrap_reflectivity, wrap_efficiency);

for(G4int i=0; i< numberOfPhotonWrap; i++){

wrap_specularSpike[i] = SsRefProb;
wrap_specularLobe[i] =SlRefProb;
wrap_backScatter [i] =BsRefProb;
//wrap_Lambertian[i] = LamRefProb;

refractiveIndexForGap [i] =refIndexForGap;

}
  
  
  
  G4MaterialPropertiesTable *wrappingSurfaceMatProp = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

wrappingSurfaceMatProp->AddProperty("RINDEX", wrap_pe, refractiveIndexForGap, numberOfPhotonWrap);
wrappingSurfaceMatProp->AddProperty("REFLECTIVITY", wrap_pe, wrap_reflectivity, numberOfPhotonWrap);
wrappingSurfaceMatProp->AddProperty("EFFICIENCY", wrap_pe, wrap_efficiency, numberOfPhotonWrap);

wrappingSurfaceMatProp->AddProperty("SPECULARLOBECONSTANT", wrap_pe, wrap_specularLobe, numberOfPhotonWrap);
wrappingSurfaceMatProp->AddProperty("SPECULARSPIKECONSTANT", wrap_pe, wrap_specularSpike, numberOfPhotonWrap);
wrappingSurfaceMatProp->AddProperty("BACKSCATTERCONSTANT", wrap_pe, wrap_backScatter, numberOfPhotonWrap);
//Lambertial 1-(diğer olasılıklar)..o yüzden belirtmeye gerek yok


 G4cout << "-------Wrapping material(BackPainted) surface material property--------" << G4endl;
 wrappingSurfaceMatProp->DumpTable();

wrappingSurface->SetMaterialPropertiesTable(wrappingSurfaceMatProp);



}



void G1DetectorConstruction::SetFP() //for PFP and GFP
{

/*

BC-620/Ej-510 REflector paint for plastic scintillator.Reflectivity grafiğinden 400-600 nm aralığındaki fotonlar için reflectivity yaklaşık 0.96.yani sabit


The finish PFP (GFP) defines a volume with a painted surface, e.g. a scintillator coated with a reflecting material, where R is the probability for reflection by the paint. The GEANT4 code decides, with respect to R, if the optical photon is absorbed by the paint. If not absorbed the optical photon is SS (L) reflected. No refraction occurs and Snell’s law is not applied. Note: the paint is defined by the surface parameters; it is not required to define it as an additional geometrical volume.
it is perfectly coated scintillator.

GFP-------Only reflection or absorption; No refraction; Reflection probability set by Reflectivity. Only lambertian reflectin 

PFP-------Only reflection or absorption; No refraction; Reflection probability set by Reflectivity.if reflected only Specular Spike reflection

*/

G4OpticalSurface* wrappingSurface = new G4OpticalSurface("Trp-surface");
  wrappingSurface->SetType(dielectric_dielectric);
  wrappingSurface->SetModel(unified);
  wrappingSurface->SetFinish(fSurfaceFinish);

   
      new G4LogicalSkinSurface("ScintillatorSkinSurface", fScin_logic, wrappingSurface);
     new G4LogicalSkinSurface("LightGuideSkinSurface", fTrp_logic, wrappingSurface);
     new G4LogicalSkinSurface("OpticalCementSkinSurface", fOptical_cement_logic, wrappingSurface); 
     new G4LogicalSkinSurface("OpticalCement1SkinSurface", fOptical_cement_logic1, wrappingSurface);
          
          
const G4int numberOfPhotonWrap = 23;

G4double wrap_pwl[numberOfPhotonWrap];
G4double wrap_pe[numberOfPhotonWrap];
G4double wrap_reflectivity[numberOfPhotonWrap]; //Reflection probability versus absorbtion.no refraction
G4double wrap_efficiency[numberOfPhotonWrap];

ReadWrapReflectivityFromFile("wrappingReflectivity.txt",wrap_pwl, wrap_pe, wrap_reflectivity, wrap_efficiency);
  
  
  G4MaterialPropertiesTable *wrappingSurfaceMatProp = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

wrappingSurfaceMatProp->AddProperty("REFLECTIVITY", wrap_pe, wrap_reflectivity, numberOfPhotonWrap);
wrappingSurfaceMatProp->AddProperty("EFFICIENCY", wrap_pe, wrap_efficiency, numberOfPhotonWrap);

G4cout << "-------Wrapping Material(FrontPainted) surface material property--------" << G4endl;
wrappingSurfaceMatProp->DumpTable();

wrappingSurface->SetMaterialPropertiesTable(wrappingSurfaceMatProp);



}


void G1DetectorConstruction::SetInsideSurfaceGround(G4OpticalSurface *OpISG){

OpISG->SetSigmaAlpha(fInsideSigmaAlpha);

 const G4int num = 2;
  G4double ephoton[num] = { 2.901*eV, 2.902*eV };
  
  G4double reflectivity[num] = { fReflectivity, fReflectivity };
  G4double efficiency[num]   = { fEfficiency, fEfficiency };
  
G4double specularlobe[num] = { 0.3, 0.3 };
G4double specularspike[num] = { 0.0, 0.0 };
G4double backscatter[num] = { 0.7, 0.7 };
        //G4double ss[num] = {1.0,1.0};




G4MaterialPropertiesTable *InsideSurfaceMatProp = new G4MaterialPropertiesTable();
//AddProperty (const char *key, G4double *PhotonEnergies, G4double *PropertyValues, G4int NumEntries) ilk parametre keyfi değil..hata vermiyot ama etkiliyot

InsideSurfaceMatProp->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
InsideSurfaceMatProp->AddProperty("EFFICIENCY", ephoton, efficiency, num);
InsideSurfaceMatProp->AddProperty("SPECULARLOBECONSTANT", ephoton, specularlobe, num);
InsideSurfaceMatProp->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularspike, num);
InsideSurfaceMatProp->AddProperty("BACKSCATTERCONSTANT", ephoton, backscatter, num);
//Lambertial 1-(diğer olasılıklar)..o yüzden belirtmeye gerek yok

G4cout << "-------Ej_200 surface material property--------" << G4endl;
InsideSurfaceMatProp->DumpTable();

OpISG->SetMaterialPropertiesTable(InsideSurfaceMatProp);


}




void G1DetectorConstruction::SetSigmaAlpha( G4double sa)
{

this->fSigmaAlpha = sa;
G4RunManager::GetRunManager()->ReinitializeGeometry();

}


//dielectric_dielectric yüzeyi için
void G1DetectorConstruction::SetEfficiency( G4double efficiency )
{

this->fEfficiency = efficiency;

G4RunManager::GetRunManager()->ReinitializeGeometry(); //gerekli yoksa aktif olmuyor

}


//dielectric_dielectric yüzeyi için
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

this->fScin_size_x = tv[0];
this->fScin_size_y = tv[1];
this->fScin_size_z = tv[2];

G4RunManager::GetRunManager()->ReinitializeGeometry();  //gerekli yoksa aktif olmuyor

}

void G1DetectorConstruction::SetOpFinishType( G4int type )
{

  switch (type){
  
  case 0 :
  this->fSurfaceFinish = polished;
  break;
  
  case 1 :
  this->fSurfaceFinish = polishedfrontpainted;
  break;
  
  case 2 :
  this->fSurfaceFinish = polishedbackpainted;
  break;
  
  case 3 :
  this->fSurfaceFinish = ground;
  break;
  
  case 4 :
  this->fSurfaceFinish = groundfrontpainted;
  break;
  
  case 5 :
  this->fSurfaceFinish = groundbackpainted;
  break;
  
  default:
  break;
  
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
//G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void G1DetectorConstruction::SetDefaults()
{

/*
//Module 3
fCementRefractiveIndex = 1.43; //
fTrpRefractiveIndex = 1.43; //ScintillatorREfIndex =1.58..dosyadan okunuyor..lightGuide 1.502
*/


fCementRefractiveIndex = 1.57; //EJ-500
fTrpRefractiveIndex = 1.502; //ScintillatorREfIndex =1.58..dosyadan okunuyor..lightGuide 1.502


numberOfUnitInLayerX =1; //layer içinde x doğruldutunda
numberOfLayerAlongY =1;



fPmt_innerRadius = 0;
fPmt_outerRadius = 2.3*cm; //photocathode area(diameter 4.6cm) 

//fPmt_outerRadius = 3.5*cm;

fPmt_Length_z = 15*cm; //önemi yok



fOptical_cement_size_z = 0.0125*cm;

//fOptical_cement_size_z = 0.025*cm; //module 3 için

fTrp_Length_x = 10*cm;
fTrp_Length_y = 10*cm;
fTrp_Length_z = 10*cm;

//fTrp_Length_z = 0.05*cm; //module 3

//EJ-200
fScin_size_x = 10*cm;
fScin_size_y = 10*cm;
fScin_size_z = 100*cm;






//sıfır olursa her taraf sintilator materyali gibi oluyor.dolayısıyla fotonlar yansımıyor.wrapping thickness
fDist_scin_unit_x = 10e-05*cm;
fDist_scin_unit_y = 10e-05*cm;



fUnit_size_x = fScin_size_x+fDist_scin_unit_x;
fUnit_size_y = fScin_size_y+fDist_scin_unit_y;
fUnit_size_z = fScin_size_z + 2*fPmt_Length_z + 2*fTrp_Length_z + 4*fOptical_cement_size_z;


//x doğrultusunda yerleşiyor
fLayer_size_x = numberOfUnitInLayerX*fUnit_size_x;
fLayer_size_y = fUnit_size_y;
fLayer_size_z = fUnit_size_z;


//replicated geometry olduğu için buralarda boşluk olamaz
fHall_size_x = fLayer_size_x;
fHall_size_y = numberOfLayerAlongY*fLayer_size_y;
fHall_size_z = fLayer_size_z;

//World Size, 5cm away for each scin surface
fWorldSize_x = fHall_size_x+1*cm ;
fWorldSize_y = fHall_size_y+1*cm ;
fWorldSize_z = fHall_size_z+1*cm ;


//Scintillator const property
fScintillationYield = 6000.0;
//fScintillationYield = 10000.0;

fResolutionScale = 1.0;
fFastTimeConstant = 2.1;
fYieldRatio = 1.0;


//Surface Finish Type.Burası sintilatorun kaplanmış olup olmadığını gösterecek.
fSurfaceFinish = polishedbackpainted;


fInsideSurfaceFinish = polished; //ya polish ya da ground olur.
if( fInsideSurfaceFinish == 3 )
fInsideSigmaAlpha = 0.6;



   if ( fSurfaceFinish == 0 || fSurfaceFinish == 1 || fSurfaceFinish == 4 ) // P and FrontPainted(PFP,GFP), we dont need sigmaAlpha
	{
	fSigmaAlpha = 0.0;
	} 
   else // G and BackPainted(GBP,PBP) ,we nedd sigmaAlpha
	{

	fSigmaAlpha = 0.25; //bu değer 0 ile 1 arasında olmalı.yoksa backpaint durumlarında foton refraction yapım sintilatorden ayrılıyor.
	}



if( fSurfaceFinish == 0 || fSurfaceFinish == 3 ){ // for not wrapped scintillator , Polish and Ground
fIsScinWrapped = false;
fReflectivity = 1.0;
fEfficiency = 0.0;

} else // for wrapped scintillator. PFP,GFP,PBP,GBP
{
 
fIsScinWrapped = true;

fIsWrapHasTeoricReflectivity = true; //bunu değiştirip bakacağız
if(fIsWrapHasTeoricReflectivity == true)
fWrappTeoricReflectivity = 0.98;

}

//this is for wrapping surface:  PFP,PBP,GFP,GBP



fIsPmtHasTeoricQE = false;
if(fIsPmtHasTeoricQE == true)
fPmtTeoricQE = 1.0;




}

void G1DetectorConstruction::Print() const
{

G4cout<<"----------------- Detector Property----------------"<<G4endl;
G4cout<< "Detector Dimension: "<<"x: "<<G4BestUnit(fScin_size_x,"Length")<<" y:"<<G4BestUnit(fScin_size_y, "Length")<<" z:"<<
G4BestUnit(fScin_size_z,"Length")<<G4endl;

G4cout<<"Scintillation Yield: "<<fScintillationYield<<G4endl;
G4cout<< "Detector Material: "<<fScin_mat->GetName()<<G4endl;
G4cout<< "Detector Mass: "<<G4BestUnit(fScin_logic->GetMass(),"Mass")<<G4endl;
G4cout<<"Resolution Scale: "<<fResolutionScale<<G4endl;
G4cout<<"InnerSurface Finish Type: "<<fInsideSurfaceFinish<<G4endl;
G4cout<<"Scintilator Surface Finish Type: "<<fSurfaceFinish<<G4endl;




if( fIsScinWrapped == false ) //it not wrapped, Polish or ground
{

  G4cout <<"Scintilator is not wrapped"<<G4endl;
  
    if( fSurfaceFinish == 3 ) // it is ground
    G4cout<< "it is ground surface. Sigma alpha: "<<fSigmaAlpha<<G4endl;
    
G4cout<<"Sintillator surface reflectivity(maybe dirty): "<<fReflectivity<<G4endl;
G4cout<<"Sintillator surface efficiency: "<<fEfficiency<<G4endl;

}
else{  // it is wrapped. PFP,GFP,PBP,GBP

  G4cout <<"Scintilator is wrapped"<<G4endl;
  if( fSurfaceFinish == 1 || fSurfaceFinish == 4 ) // it is FP. no need sigma alpha
G4cout <<"it is FrontPainted. No need SigmaAlpha"<<G4endl;
  if( fSurfaceFinish == 2 || fSurfaceFinish == 5 ) // it is BP.
  G4cout<< "it is BackPainted. SigmaAlpha for surface between air gap and scintillator : "<<fSigmaAlpha<<G4endl;


}


if( fIsPmtHasTeoricQE == true ) 
 G4cout<<"Teoric Quantum Efficiency value of PMT : "<<fPmtTeoricQE<<G4endl; 


if( fIsScinWrapped == true && fIsWrapHasTeoricReflectivity == true)
  G4cout<<"Teroic Reflectivity value of wrapping: "<<fWrappTeoricReflectivity<<G4endl;





G4cout<<"-------------------------------------------"<<G4endl;

}






G4double G1DetectorConstruction::ConvertMassDensity (G4double ad, G4double mw){

G4double AtomDensity = ad*(1/cm*cm*cm);
G4double MolecularWeight = mw*(g/mole);
//hangi birimde ekranda görmek istiyorsan o birime böl.yoksa default km3 gösteriyor
//G4cout <<"AtomDensity"<<AtomDensity/(1/cm*cm*cm)<<G4endl;
G4double massDensity = (AtomDensity/N_A)*MolecularWeight;

return massDensity;

}



//read emission spectrum
void G1DetectorConstruction::ReadEmissionSpectrumFromFile(const G4String &filename, G4double* pwl , G4double* pe ,G4double* sf, G4double *ri, G4double* ab ){


  G4int i=0;
  G4String line;
  
  G4double waveLength;
  G4double relativeStrength;
  G4double refIndex;
  G4double abs;
  
  
std::ifstream myfile (filename.c_str());


  if (myfile.is_open())
  {
  
    while ( getline (myfile,line) )
    {
    
        std::istringstream iss(line);
        
        
        if ( (iss >> waveLength >> relativeStrength >> refIndex >> abs) ) 
         {
        
         pwl[i] = waveLength*nm;
         
         pe[i]=(1240*eV*nm)/pwl[i];
         
        // G4cout <<G4BestUnit(pe[i],"Energy")<<G4endl;
         //default olarak ekranda uzunluk mm energi MeV cinsinden cıkıyor
         sf [i] = relativeStrength; 
         ri[i] = refIndex;
         ab[i] = abs*cm;
       
       i++;
      
         
         } 
         
      
    }
    
     
     
    
    myfile.close();
  }

  else G4cout << "Unable to open file"; 


}

// Read PMT
void G1DetectorConstruction::ReadPmtEfficiencyFromFile(const G4String &filename, G4double* pwl, G4double* pe, G4double *eff, G4double *reff){

G4String line;
G4int i=0;

G4double waveLength;
G4double efficiency;
  
  
std::ifstream myfile (filename.c_str());


  if (myfile.is_open())
  {
  
    while ( getline (myfile,line) )
    {
    
        std::istringstream iss(line);
        
        
        if ( (iss >> waveLength >> efficiency) ) 
         {
        
         pwl[i] = waveLength*nm;
         
         pe[i]=(1240*eV*nm)/pwl[i];
         
         if( fIsPmtHasTeoricQE == true ){ 
         eff[i] = fPmtTeoricQE;
         
         }
         else{
         eff[i] = efficiency;
         }
         
         reff[i] = 0; //Hepsini absorplasın.No reflection on photoCathode surface.
         //G4cout <<"aaaaaaaaaaa:  "<<eff[i]<<G4endl;
        i++;
      
         
         } 
         
      
    }
    
     
     
    
    myfile.close();
  }

  else G4cout << "Unable to open file"; 


}


void G1DetectorConstruction::ReadWrapReflectivityFromFile(const G4String &filename, G4double* pwl, G4double* pe, G4double *reffcoeff, G4double * eff){

G4String line;
G4int i=0;

G4double waveLength;
G4double reflectivityCoefficient;
  
  
std::ifstream myfile (filename.c_str());


  if (myfile.is_open())
  {
  
    while ( getline (myfile,line) )
    {
    
        std::istringstream iss(line);
        
        
        if ( (iss >> waveLength >> reflectivityCoefficient) ) 
         {
        
         pwl[i] = waveLength*nm;
         pe[i]=(1240*eV*nm)/pwl[i];
         
         if(fIsWrapHasTeoricReflectivity == true){
         reffcoeff[i] = fWrappTeoricReflectivity;
          
         }
         else{
         
           reffcoeff[i] = reflectivityCoefficient;
         
         //G4cout <<"ReflectionCoefficient: "<<reffcoeff[i]<<G4endl;
         }
         
         eff[i] = 0; //no detected
        i++;
      
         
         } 
         
      
    }
    
     
     
    
    myfile.close();
  }

  else G4cout << "Unable to open file"; 


}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
