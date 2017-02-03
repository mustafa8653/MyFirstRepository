#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh" //MeV burdan geliyo

#include <vector>
#include <iostream>
#include <map>
#include <utility>
#include <string>
#include <iomanip> //setprecision burdan geliyo

#ifndef G1UserEventInformation_h
#define G1UserEventInformation_h 1

class G1UserEventInformation : public G4VUserEventInformation
{
public:
   G1UserEventInformation( G4int eventID );
  ~G1UserEventInformation();
  
 inline virtual void Print() const {};

  void IncPhotonCount_Scint() { fPhotonCount_Scint++; }
  void IncPhotonCount_Ceren() { fPhotonCount_Ceren++; }
  
  G4int GetPhotonCount_Scint() const { return fPhotonCount_Scint; }
  G4int GetPhotonCount_Ceren() const { return fPhotonCount_Ceren; }
  
  //Gets the total optical photon count produced
  G4int GetPhotonCount() { return fPhotonCount_Scint+fPhotonCount_Ceren; }
  
  //Get total nember of photon detected from pmt surface
  G4int GetPmtPhotonCount() { return (*fProcessMap)[10]; }
  
  void AddEdep(G4double edep) { fEdep += edep; }
  G4double GetEdep() { return fEdep; }
  
  
  
  void PrintEventResult();

std::map < G4int,G4int>* GetProcessMap()  { return fProcessMap; }

std::vector <G4double>* GetPhotonEnergyVec() { return fPhotonEnergyVec; }

private:

  G4int fEventID;
  
  G4int fPhotonCount_Scint;
  G4int fPhotonCount_Ceren;
  
  G4double fEdep;
   
  std::vector < G4double > *fPhotonEnergyVec;
  
  std::map < G4int,G4int> *fProcessMap;

};

#endif


