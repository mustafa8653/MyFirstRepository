#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#ifndef PWUserEventInformation_h
#define PWUserEventInformation_h 1

class PWUserEventInformation : public G4VUserEventInformation
{
public:
  PWUserEventInformation();
  ~PWUserEventInformation();
  
  inline void Print()const{};

  void IncPhotonCount_Scint(){photonCount_Scint++;}
  void IncPhotonCount_Ceren(){photonCount_Ceren++;}
  void IncEDep(G4double dep){totE+=dep;}
  void IncAbsorption(){absorptionCount++;}
  void IncBoundaryAbsorption(){boundaryAbsorptionCount++;}
  void IncPMTHitCount(G4int i=1){pmtHitCount+=i;}

  void SetPosMax(const G4ThreeVector& p,G4double edep){posMax=p;edepMax=edep;}

  G4int GetPhotonCount_Scint()const {return photonCount_Scint;}
  G4int GetPhotonCount_Ceren()const {return photonCount_Ceren;}
  G4double GetEDep()const {return totE;}
  G4int GetAbsorptionCount()const {return absorptionCount;}
  G4int GetBoundaryAbsorptionCount() const {return boundaryAbsorptionCount;}
  G4int GetPMTHitCount()const {return pmtHitCount;}  

  G4ThreeVector GetPosMax(){return posMax;}
  G4double GetEDepMax(){return edepMax;}

  //Gets the total optical photon count produced
  G4int GetPhotonCount(){return photonCount_Scint+photonCount_Ceren;}

private:

  G4int pmtHitCount;
  G4int photonCount_Scint;
  G4int photonCount_Ceren;
  G4int absorptionCount;
  G4int boundaryAbsorptionCount;

  G4double totE;  

  //These only have meaning if totE > 0
  //If totE = 0 then these wont be set by EndOfEventAction

  G4ThreeVector posMax;
  G4double edepMax;

};

#endif




