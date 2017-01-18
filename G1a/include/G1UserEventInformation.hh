#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <iostream>
#include <map>
#include <utility>
#include <string>

#ifndef G1UserEventInformation_h
#define G1UserEventInformation_h 1

class G1UserEventInformation : public G4VUserEventInformation
{
public:
   G1UserEventInformation();
  ~G1UserEventInformation();
  
 inline virtual void Print() const {};

  void IncPhotonCount_Scint() { fPhotonCount_Scint++; }
  void IncPhotonCount_Ceren() { fPhotonCount_Ceren++; }
  G4int GetPhotonCount_Scint() const { return fPhotonCount_Scint; }
  G4int GetPhotonCount_Ceren() const { return fPhotonCount_Ceren; }
  //Gets the total optical photon count produced
  G4int GetPhotonCount() { return fPhotonCount_Scint+fPhotonCount_Ceren; }
  
  void AddEdep(G4double edep) { fEdep += edep; }
  G4double GetEdep() { return fEdep; }
  
  void PrintProcessNumber(); 
  

std::map < G4int,G4int>* GetProcessMap()  { return fProcessMap; }

private:

  
  G4int fPhotonCount_Scint;
  G4int fPhotonCount_Ceren;
  
  G4double fEdep;
   

  std::map < G4int,G4int> *fProcessMap;

};

#endif


