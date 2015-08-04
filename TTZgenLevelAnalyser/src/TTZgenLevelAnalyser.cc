// -*- C++ -*-
//
// Package:    TTZgenLevelAnalyser
// Class:      TTZgenLevelAnalyser
// 
/**\class TTZgenLevelAnalyser TTZgenLevelAnalyser.cc ttZgenLevelAnalyser/TTZgenLevelAnalyser/src/TTZgenLevelAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexander D J Morton
//         Created:  Mon Aug  3 16:54:13 BST 2015
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <boost/shared_ptr.hpp>

// Gen level products
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TH1.h"

//
// class declaration
//

class TTZgenLevelAnalyser : public edm::EDAnalyzer {
public:
  explicit TTZgenLevelAnalyser(const edm::ParameterSet&);
  ~TTZgenLevelAnalyser();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  TH1D* h_topMass;
  TH1D* h_electronMass;
  TH1D* h_muonMass;
  TH1D* h_tauMass;
  TH1D* h_zMass;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTZgenLevelAnalyser::TTZgenLevelAnalyser(const edm::ParameterSet& iConfig)
  :  h_electronMass(0), h_muonMass(0), h_tauMass(0)
{
   //now do what ever initialization is needed

   //register your products
  //   produces<edm::HepMCProduct>();


}


TTZgenLevelAnalyser::~TTZgenLevelAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TTZgenLevelAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< LHEEventProduct > lheEventSrc;
  iEvent.getByLabel( "source", lheEventSrc) ;
  
  //  HepMC::GenEvent* Evt = new HepMC::GenEvent();

  for (int i = 0; i < lheEventSrc->hepeup().NUP; ++i) {
    
    int pgId = lheEventSrc->hepeup().IDUP[i];
    
    HepMC::FourVector p(lheEventSrc->hepeup().PUP[i][0], lheEventSrc->hepeup().PUP[i][1], lheEventSrc->hepeup().PUP[i][2], lheEventSrc->hepeup().PUP[i][3]);

    //lheEventSrc->hepeup().PUP[i][4]; generated mass

    //    std::cout << "i: " << i << std::endl;
    
    if ( fabs(pgId) == 6 ) {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Top." << std::endl;
      h_topMass->Fill(p.m());
    }

    else if ( fabs(pgId) == 11 ) {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Electron." << std::endl;
      h_electronMass->Fill(p.m());
    }

    else if ( fabs(pgId) == 13 ) {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Muon." << std::endl;
      h_muonMass->Fill(p.m());
    }

    else if ( fabs(pgId) == 15 ) {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Tauon." << std::endl;
      h_tauMass->Fill(p.m());
    }

    else if ( fabs(pgId) == 22 ) {
      //     std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Gamma." << std::endl;
    }

    else if ( fabs(pgId) == 23 ) {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: Z." << std::endl;
      h_zMass->Fill(p.m());
    }

    else {
      //      std::cout << "Status: " << lheEventSrc->hepeup().ISTUP[i] << "; pgId: " << pgId << std::endl;
    }

    //    std::cout << "Invariant Mass: " << p.m() << std::endl;
    
  }

}


// ------------ method called once each job just before starting event loop  ------------
void TTZgenLevelAnalyser::beginJob()
{
  edm::Service<TFileService> fs;

  h_topMass = fs->make<TH1D>( "h_topMass", "e+/e- inv. mass", 100, 0.0, 10.0 );
  h_electronMass = fs->make<TH1D>( "h_electronMass", "e+/e- inv. mass", 100, 0.0, 10.0 );
  h_muonMass = fs->make<TH1D>( "h_muonMass", "mu+/mu- inv. mass", 100, 60.0, 120.0);
  h_tauMass = fs->make<TH1D>( "h_tauMass", "tau+/tau- inv. mass", 100, 1400.0, 2400.0);
  h_zMass = fs->make<TH1D>( "h_zMass", "Z inv. mass", 100, 1400.0, 2400.0);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTZgenLevelAnalyser::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void TTZgenLevelAnalyser::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTZgenLevelAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TTZgenLevelAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TTZgenLevelAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTZgenLevelAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTZgenLevelAnalyser);
