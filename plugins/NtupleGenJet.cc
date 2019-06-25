// -*- C++ -*-
//
// Package:    validation/NtupleGenJet
// Class:      NtupleGenJet
// 
/**\class NtupleGenJet NtupleGenJet.cc validation/NtupleGenJet/plugins/NtupleGenJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lata Panwar
//         Created:  Wed, 14 Mar 2018 12:41:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtupleGenJet : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NtupleGenJet(const edm::ParameterSet&);
      ~NtupleGenJet();
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
       edm::EDGetTokenT <reco::GenParticleCollection> genparticlesToken;
	
        int genHiggs_n_=0;
        TH1D *nHiggs_histo;

        TH1F *pt_histo_H1;
        TH1F *eta_histo_H1;
	TH1F *phi_histo_H1;  
        TH1F *e_histo_H1;
        TH1F *mass_histo_H1;

        TH1F *pt_histo_H2;
        TH1F *eta_histo_H2;
        TH1F *phi_histo_H2;
        TH1F *e_histo_H2;
        TH1F *mass_histo_H2;

        TH1F *pt_histo_lead_pho;
        TH1F *eta_histo_lead_pho;
        TH1F *phi_histo_lead_pho;
        TH1F *e_histo_lead_pho;

        TH1F *pt_histo_sublead_pho;
        TH1F *eta_histo_sublead_pho;
        TH1F *phi_histo_sublead_pho;
        TH1F *e_histo_sublead_pho;

	TH1F *pt_histo_lead_b;
        TH1F *eta_histo_lead_b;
        TH1F *phi_histo_lead_b;
        TH1F *e_histo_lead_b;
        
	TH1F *pt_histo_sublead_b;
        TH1F *eta_histo_sublead_b;
        TH1F *phi_histo_sublead_b;
        TH1F *e_histo_sublead_b;

        TH1F *mass_HH;
       TLorentzVector p4_y1, p4_y2, p4_b, p4_bbar, p4_H1, p4_H2;
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
NtupleGenJet::NtupleGenJet(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   genparticlesToken	   = consumes <reco::GenParticleCollection> (std::string("genParticles"));
   edm::Service<TFileService> fs;
   nHiggs_histo = fs->make<TH1D>("N_higgs" , ";N_{H};Events;;" , 50 , 0 , 50 );
   pt_histo_H1  = fs->make<TH1F>("pT_H1"    , ";p_{T} of Higgs[GeV](#gamma#gamma);Events;;" , 100 , 0 , 500 );
   eta_histo_H1 = fs->make<TH1F>("eta_H1"   , ";#eta of Higgs(#gamma#gamma);Events;;" , 50 , -5 , 5 );
   phi_histo_H1 = fs->make<TH1F>("phi_H1"   , ";#phi of Higgs(#gamma#gamma);Events;;" , 50 , -5 , 5 );
   e_histo_H1  = fs->make<TH1F>("e_H1"      , ";energy of Higgs[GeV](#gamma#gamma);Events;;" , 100 , 0 , 500 );

   pt_histo_H2  = fs->make<TH1F>("pT_H2"  , ";p_{T} of Higgs[GeV](bb);Events;;" , 100 , 0 , 500 );
   eta_histo_H2 = fs->make<TH1F>("eta_H2" , ";#eta of Higgs(bb);Events;;" , 50 , -5 , 5 );
   phi_histo_H2 = fs->make<TH1F>("phi_H2" , ";#phi of Higgs(bb);Events;;" , 50 , -5 , 5 );
   e_histo_H2  = fs->make<TH1F>("e_H2"    , ";energy of Higgs[GeV](bb);Events;;" , 100 , 0 , 500 );

   mass_histo_H1=fs->make<TH1F>("mass_H1" , ";M_{#gamma#gamma} [GeV];Events;;" , 10000, 124.99, 125.05);
   mass_histo_H2=fs->make<TH1F>("mass_H2" , ";M_{bb~} [GeV];Events;;" , 10000, 124.95, 125.05);

   pt_histo_lead_pho = fs->make<TH1F>("pT_y1"  , ";p_{T} #gamma1[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_lead_pho= fs->make<TH1F>("eta_y1" , ";#eta #gamma1;Events;;" , 50 , -5 , 5 );
   phi_histo_lead_pho= fs->make<TH1F>("phi_y1" , ";#phi #gamma1;Events;;" , 50 , -5 , 5 );
   e_histo_lead_pho = fs->make<TH1F>("e_y1"   , ";energy #gamma1[GeV];Events;;" , 100 , 0 , 500 );

   pt_histo_sublead_pho = fs->make<TH1F>("pT_y2" , ";p_{T} #gamma2[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_sublead_pho= fs->make<TH1F>("eta_y2" , ";#eta #gamma2;Events;;" , 50 , -5 , 5 );
   phi_histo_sublead_pho= fs->make<TH1F>("phi_y2" , ";#phi #gamma2;Events;;" , 50 , -5 , 5 );
   e_histo_sublead_pho = fs->make<TH1F>("e_y2"   , ";energy #gamma2[GeV];Events;;" , 100 , 0 , 500 );

   pt_histo_lead_b = fs->make<TH1F>("pT_b1" , ";p_{T} of b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_lead_b=fs->make<TH1F>("eta_b1" , ";#eta b;Events;;" , 50 , -5 , 5 );
   phi_histo_lead_b=fs->make<TH1F>("phi_b1" , ";#phi b;Events;;" , 50 , -5 , 5 );
   e_histo_lead_b = fs->make<TH1F>("e_b1"   , ";energy b[GeV];Events;;" , 100 , 0 , 500 );

   pt_histo_sublead_b = fs->make<TH1F>("pT_b2" , ";p_{T} b~[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_sublead_b=fs->make<TH1F>("eta_b2" , ";#eta b~;Events;;" , 50 , -5 , 5 );
   phi_histo_sublead_b=fs->make<TH1F>("phi_b2" , ";#phi b~;Events;;" , 50 , -5 , 5 );
   e_histo_sublead_b = fs->make<TH1F>("e_b2" , ";energy b~[GeV];Events;;" , 100 , 0 , 500 );

   mass_HH=fs->make<TH1F>("mass_HH" , ";M_{HH} [GeV];Events;;" , 100, 0,1000);

   /*   pt_histo_add_b = fs->make<TH1F>("pT_b" , ";p_{T} of additional b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_add_b=fs->make<TH1F>("eta_b" , ";#eta of additional b;Events;;" , 50 , -5 , 5 );
   phi_histo_add_b=fs->make<TH1F>("phi_b" , ";#phi of additional b;Events;;" , 50 , -5 , 5 );*/
}



NtupleGenJet::~NtupleGenJet()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleGenJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genparticlesToken, genParticles);
   //   TLorentzVector p4_y1, p4_y2, p4_b, p4_bbar, p4_H1, p4_H2; 

   for(size_t i = 0; i < genParticles->size(); ++ i) {
     const reco::GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();
  
//     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     if (id == 25){
    int n = p.numberOfDaughters();
    if(n < 2 ) continue;
    const reco::Candidate * d1 = p.daughter( 0 );
    const reco::Candidate * d2 = p.daughter( 1 );

    if (std::abs(d1->pdgId())==22 && std::abs(d2->pdgId())==22){
     ++genHiggs_n_;
     nHiggs_histo->Fill(genHiggs_n_);

     pt_histo_H1->Fill(p.pt());
     eta_histo_H1->Fill(p.eta());
     phi_histo_H1 ->Fill(p.phi());
     e_histo_H1 ->Fill(p.energy());
     mass_histo_H1->Fill(p.mass());
 
    //// plotting for b's
    pt_histo_lead_pho->Fill(d1->pt());
    eta_histo_lead_pho->Fill(d1->eta());
    phi_histo_lead_pho->Fill(d1->phi());
    e_histo_lead_pho->Fill(d1->energy());
    p4_y1.SetPtEtaPhiE(d1->pt(), d1->eta(), d1->phi(), d1->energy());

    pt_histo_sublead_pho->Fill(d2->pt());
    eta_histo_sublead_pho->Fill(d2->eta());
    phi_histo_sublead_pho->Fill(d2->phi());
    e_histo_sublead_pho->Fill(d2->energy());
    p4_y2.SetPtEtaPhiE(d2->pt(), d2->eta(), d2->phi(), d2->energy());
    p4_H1 = p4_y1 + p4_y2;   
    }
    if (std::abs(d1->pdgId())==5 && std::abs(d2->pdgId())==5){

      ++genHiggs_n_;
      nHiggs_histo->Fill(genHiggs_n_);

      pt_histo_H2->Fill(p.pt());
      eta_histo_H2->Fill(p.eta());
      phi_histo_H2->Fill(p.phi());
      e_histo_H2 ->Fill(p.energy());
      mass_histo_H2->Fill(p.mass());
      //// plotting for b's                                                                                                                                                                                                                 
      if(d1->pdgId() == 5 && d2->pdgId() == -5){
      pt_histo_lead_b->Fill(d1->pt());
      eta_histo_lead_b->Fill(d1->eta());
      phi_histo_lead_b->Fill(d1->phi());
      e_histo_lead_b->Fill(d1->energy());
      p4_b.SetPtEtaPhiE(d1->pt(), d1->eta(), d1->phi(), d1->energy());

      pt_histo_sublead_b->Fill(d2->pt());
      eta_histo_sublead_b->Fill(d2->eta());
      phi_histo_sublead_b->Fill(d2->phi());
      e_histo_sublead_b->Fill(d2->energy());
      p4_bbar.SetPtEtaPhiE(d2->pt(), d2->eta(), d2->phi(), d2->energy());
      
      p4_H2 = p4_b + p4_bbar;
      }
      else if (d1->pdgId() == -5 && d2->pdgId() == 5){
	pt_histo_lead_b->Fill(d2->pt());
	eta_histo_lead_b->Fill(d2->eta());
	phi_histo_lead_b->Fill(d2->phi());
	e_histo_lead_b->Fill(d2->energy());
        p4_b.SetPtEtaPhiE(d2->pt(), d2->eta(), d2->phi(), d2->energy());

	pt_histo_sublead_b->Fill(d1->pt());
	eta_histo_sublead_b->Fill(d1->eta());
	phi_histo_sublead_b->Fill(d1->phi());
	e_histo_sublead_b->Fill(d1->energy());
        p4_bbar.SetPtEtaPhiE(d1->pt(), d1->eta(), d1->phi(), d1->energy());
	p4_H2 = p4_b + p4_bbar;
      }
    }
    mass_HH->Fill((p4_H1 + p4_H2).M());
     }
 
   }
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleGenJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleGenJet);
