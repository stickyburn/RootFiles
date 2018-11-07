#include <EventLoop/Worker.h>

#include "TSystem.h"
//#include "TMath.h" 

#include <CxAODReader_HH_bbtautau/AnalysisReader_hhbbtt.h>

#include <CxAODTools_HH_bbtautau/HHbbtautauLepHadSelection.h>
#include <CxAODTools_HH_bbtautau/HHbbtautauHadHadSelection.h>
#include <CxAODTools_HH_bbtautau/HHbbtautauLepHadJetSelection.h>
#include <CxAODTools_HH_bbtautau/HHbbtautauHadHadJetSelection.h>
#include <CxAODTools_HH_bbtautau/HHbbtautauBoostedSelection.h>
#include <CxAODTools_HH_bbtautau/HHbbtautauBoostedJetSelection.h>

#include "Mt2/mt2_bisect.h"

#include "BoostedJetTaggers/BoostedXbbTagger.h"

#define length(array) (sizeof(array)/sizeof(*(array)))

ClassImp(AnalysisReader_hhbbtt)

AnalysisReader_hhbbtt::AnalysisReader_hhbbtt() :
        AnalysisReader(),
  m_analysisType(""),
  m_doTruthTagging("false"),
  m_use2DbTagCut("false"),
  m_fillCR("false"),
  m_tree(nullptr),
  m_FF_File(nullptr),
  m_2DFF(nullptr),
  m_SelFJ(nullptr),
  m_SelDT(nullptr)
  //m_my_met(nullptr)
  //m_my_fatJets(nullptr),
  //m_my_ditauJets(nullptr)
  //m_tree_vbf(nullptr)
{
  m_trkBTagLimit = -0.3098; // added for Higgs Tagger

}

AnalysisReader_hhbbtt :: ~AnalysisReader_hhbbtt()
{
}

EL::StatusCode AnalysisReader_hhbbtt :: initializeSelection()
{

  m_config->getif<string>("analysisType", m_analysisType);
  m_config->getif<string>("analysisStrategy", m_analysisStrategy);

  // For 2D cut only use merged
  m_config->getif<bool>("use2DbTagCut", m_use2DbTagCut);
  if (m_use2DbTagCut){
    if(m_analysisStrategy != "Merged"){
      Error("initializeSelection()","Running with 2DbTagCut. Only possible for merged analysis, but %s is chosen. Exiting!",m_analysisStrategy.c_str());
      return EL::StatusCode::FAILURE;
    }
    Warning("initializeSelection()","Running with 2DbTagCut. b-tagging for calo jets disabled!!");
  }

  Info("initializeSelection()", "Initialize analysis '%s'.", m_analysisType.c_str());

  if (m_analysisType== "hadhad") {

    m_eventSelection = new HHbbtautauHadHadSelection(m_config);
    m_eventPostSelection = new HHbbtautauHadHadJetSelection(m_config);

    m_fillFunction = std::bind(&AnalysisReader_hhbbtt::fill_hadhad, this);

  } else if (m_analysisType=="lephad") { 

    m_eventSelection = new HHbbtautauLepHadSelection(m_config);
    m_eventPostSelection = new HHbbtautauLepHadJetSelection(m_config);

    m_fillFunction = std::bind(&AnalysisReader_hhbbtt::fill_lephad, this);

  }else if (m_analysisType == "boosted") {

    m_eventSelection = new HHbbtautauBoostedSelection(m_config);
    m_eventPostSelection = new HHbbtautauBoostedJetSelection(m_config);

    m_fillFunction = std::bind(&AnalysisReader_hhbbtt::fill_boosted, this);

  } else {
    Error("initializeSelection()", "Invalid analysis type %s", m_analysisType.c_str());
    return EL::StatusCode::FAILURE;
  }

  // TODO move to base class?
  bool writeMVATree = false;
  bool readMVA = false;
  m_config->getif< bool >("writeMVATree", writeMVATree);
  m_config->getif< bool >("readMVA", readMVA);
  // initialize MVATree
  m_tree = new MVATree_hhbbtt(writeMVATree, readMVA, m_analysisType, wk(), m_variations, false);

  //do truth tagging?
  m_config->getif<bool>("doTruthTagging", m_doTruthTagging);

  return EL::StatusCode::SUCCESS;
}

//Added this initialization since it was currectly getting called in generic reader and running on unwanted samples

EL::StatusCode AnalysisReader_hhbbtt::initializeCorrsAndSysts ()
{

  m_doTruthTagging&=m_isMC;
  Info("initializeSelection()","truth tagging1? %i",m_doTruthTagging);
  if (!m_isMC) return EL::StatusCode::SUCCESS;

  std::string comEnergy    = m_config->get<std::string>("COMEnergy");

  if ((comEnergy != "8TeV") || (comEnergy != "7TeV")) comEnergy = "13TeV";
  TString csname;
  std::string debugname;

  //if (m_analysisType == "hadhad") { csname = comEnergy + "_HadHad"; debugname = comEnergy + "_HadHad"; }
  if (m_analysisType == "hadhad" || m_analysisType == "boosted") { csname = comEnergy + "_ZeroLepton"; debugname = comEnergy + "_ZeroLepton"; }

  Info("initializeCorrsAndSysts()", "Initializing CorrsAndSysts for %s", debugname.c_str());

  m_corrsAndSysts = new CorrsAndSysts(csname);

  return EL::StatusCode::SUCCESS;
} // initializeCorrsAndSysts    


EL::StatusCode AnalysisReader_hhbbtt::initializeTools ()
{
  EL_CHECK("AnalysisReader_hhbbtt::initializeTools()", AnalysisReader::initializeTools());

  m_getFFInputs=false;
  m_config->getif<bool>("fillFakeHists", m_getFFInputs);
  m_applyFF=false;
  m_config->getif<bool>("ApplyQCDFF", m_applyFF);
  m_applyBoostedFF=false;
  m_config->getif<bool>("ApplyBoostedQCDFF", m_applyBoostedFF);
  bool m_FillSR=false;
  m_config->getif<bool>("FillSR", m_FillSR);
  bool m_FillSS=false;
  m_config->getif<bool>("FillSS", m_FillSS);

  if(m_applyFF) {
    std::string rootcore_path = gSystem->Getenv("ROOTCOREBIN");
    m_FF_File = TFile::Open((rootcore_path+"/data/FrameworkExe_HH_bbtautau/QCDFF.root").c_str(),"READ");
    if(m_FF_File==nullptr){
      std::cout << "File was not open" << std::endl;
      return EL::StatusCode::FAILURE;
    }
  }

  if(m_applyBoostedFF) {
    std::string rootcore_path = gSystem->Getenv("ROOTCOREBIN");
    m_FF_File = TFile::Open((rootcore_path+"/data/FrameworkExe_HH_bbtautau/FFBoosted.root").c_str(),"READ");
    if(m_FF_File==nullptr){
      std::cout << "File was not open" << std::endl;
      return EL::StatusCode::FAILURE;
    }
    TH1F *htemp1 =(TH1F*)m_FF_File->Get("h_N01");
    TH1F *htemp2 =(TH1F*)m_FF_File->Get("h_N02");
    m_weightFactor1=htemp1->Integral();
    m_weightFactor2=htemp2->Integral();
  }

  bool m_writeMVATree = false;
  bool m_readMVA = false;
  m_config->getif< bool >("writeMVA", m_writeMVATree);
  m_config->getif< bool >("readMVA", m_readMVA);

  m_tree = new MVATree_hhbbtt(m_writeMVATree, m_readMVA, m_analysisType, wk(), m_variations, false);
  if(m_writeMVATree || m_readMVA) m_tree->SetVariation(m_currentVar);

  m_FillSRSel260_300 = false;
  m_config->getif< bool >("FillSRSel260_300", m_FillSRSel260_300);
  m_FillSRSel400 = false;
  m_config->getif< bool >("FillSRSel400", m_FillSRSel400);
  m_FillSRSel500_600_700 = false;
  m_config->getif< bool >("FillSRSel500_600_700", m_FillSRSel500_600_700);
  m_FillSRSel800_900_1000 = false;
  m_config->getif< bool >("FillSRSel800_900_1000", m_FillSRSel800_900_1000);

  return EL::StatusCode::SUCCESS;
} //initializeTools

EL::StatusCode AnalysisReader_hhbbtt::initializeIsMC ()
{
  Info("initializeIsMC()", "Initialize isMC.");

  // get nominal event info
  // -------------------------------------------------------------
  const xAOD::EventInfo *eventInfo = m_eventInfoReader->getObjects("Nominal");

  if (!eventInfo) return EL::StatusCode::FAILURE;

  // get MC flag - different info on data/MC files
  // -----------------------------------------------
  m_isMC = Props::isMC.get(eventInfo);
  Info("initializeIsMC()", "isMC = %i", m_isMC);

  // COM energy
  std::string comEnergy = m_config->get<std::string>("COMEnergy");

  std::string xSectionFile = gSystem->Getenv("ROOTCOREBIN");
  xSectionFile      += "/data/FrameworkSub_HH_bbtautau/XSections_";
  xSectionFile      += comEnergy;
  xSectionFile      += ".txt";
  m_xSectionProvider = new XSectionProvider(xSectionFile);

  if (!m_xSectionProvider) {
    Error("initializeXSections()", "XSection provider not initialized!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
} // initializeIsMC



EL::StatusCode AnalysisReader_hhbbtt :: initializeSumOfWeights()
{
  Info("initializeSumOfWeights()", "Initialize sum of weights.");

  if (!m_isMC) {
    return EL::StatusCode::SUCCESS;
  }

  // which analysis
  std::string ana_read = "";
  if (m_analysisType == "hadhad") ana_read = "hadhad";
  else if (m_analysisType == "lephad") ana_read = "lephad";
  else if (m_analysisType == "boosted") ana_read = "boosted";
  else {
    Warning("initializeSumOfWeights()", "Invalid analysis type %s", m_analysisType.c_str());
    return EL::StatusCode::FAILURE;
  }

  // COM energy
  std::string comEnergy = m_config->get<std::string>("COMEnergy");
  std::string prodTag   = m_config->get<std::string>("prodTag");

  std::string sumOfWeightsFile = gSystem->Getenv("ROOTCOREBIN");
  sumOfWeightsFile += "/data/FrameworkSub_HH_bbtautau/yields.";
  sumOfWeightsFile += ana_read;
  sumOfWeightsFile += ".";
  sumOfWeightsFile += prodTag;
  sumOfWeightsFile += ".";
  sumOfWeightsFile += comEnergy;
  sumOfWeightsFile += ".txt";
  m_sumOfWeightsProvider = new sumOfWeightsProvider(sumOfWeightsFile);

  return EL::StatusCode::SUCCESS;
}

/*EL::StatusCode AnalysisReader_hhbbtt :: initializeOR2()
{
  m_overlapRemoval2 = new OverlapRemoval_HH_bbtautau( *m_config );
  EL_CHECK("AnalysisReader::initializeOR2()",m_overlapRemoval2->initialize());

  return EL::StatusCode::SUCCESS;
}*/

EL::StatusCode AnalysisReader_hhbbtt::setObjectsForOR(const xAOD::ElectronContainer * electrons,
				     const xAOD::PhotonContainer * photons,
				     const xAOD::MuonContainer * muons,
				     const xAOD::TauJetContainer * taus,
				     const xAOD::JetContainer * jets,
                 const xAOD::JetContainer* fatjets,
                 const xAOD::DiTauJetContainer* ditaus)
{
  
  if(electrons){
    for(const xAOD::Electron * elec : *electrons){
      Props::passPreSel.set(elec, Props::isHHLooseElectron.get(elec));
    }
  }
  
  if(muons){
    for(const xAOD::Muon * muon : *muons){
      Props::passPreSel.set(muon, Props::isHHLooseMuon.get(muon));
    }
  }
  
  if(jets){
    for(const xAOD::Jet * jet : *jets){
      Props::passPreSel.set(jet, Props::isVetoJet.get(jet));
    }
  }
  
  if(taus){
    for(const xAOD::TauJet * tau : *taus){
      Props::passPreSel.set(tau, Props::passTauSelector.get(tau));
    }
  }

  if (fatjets) {
    for (const xAOD::Jet* jet : *fatjets) {
      Props::passPreSel.set(jet,Props::isFatJet.get(jet));
    }
  }

  if (ditaus) {
    for (const xAOD::DiTauJet* ditau : *ditaus) {
      Props::passPreSel.set(ditau,Props::isDiTauJet.get(ditau));
    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader_hhbbtt :: fill_bTaggerHists(const xAOD::Jet* jet)
{
  // fills histograms per jet and jet flavour for a list of b-taggers
  std::string flav = "Data";

  if (m_isMC) {
    int label = Props::TruthLabelID.get(jet);

    flav = "L";

    if (label == 4) flav = "C";
    else if (label == 5) flav = "B";
  }

  float  SV1_IP3D = Props::SV1_IP3D.get(jet);
  float  MV2c10   = Props::MV2c10.get(jet);
  double BVar     = BTagProps::tagWeight.get(jet);

  if (m_isMC) {
    double BEff = BTagProps::eff.get(jet);
    m_histSvc->BookFillHist("eff_"   + flav, 110, -0.1, 1.1, BEff,  m_weight);
    m_histSvc->BookFillHist("pT_eff_"  + flav, 400,   0, 400, 110, 0.0, 1.1, jet->pt() / 1e3, BEff,  m_weight);
    m_histSvc->BookFillHist("eta_eff_" + flav, 100,  -5,   5, 110, 0.0, 1.1, jet->eta(), BEff,     m_weight);
  }
  m_histSvc->BookFillHist("btag_weight_"   + flav, 110, -1.1, 1.1, BVar,   m_weight);
  m_histSvc->BookFillHist("SV1_IP3D_" + flav, 110,  -30,  80, SV1_IP3D, m_weight);
  m_histSvc->BookFillHist("MV2c10_"   + flav, 110, -1.1, 1.1, MV2c10,   m_weight);

  

  return EL::StatusCode::SUCCESS;
}

void AnalysisReader_hhbbtt::compute_btagging (const std::vector<const xAOD::Jet*> &signalJets)
{
  if (m_trackJets) {
    m_bTagTool->setJetAuthor(m_trackJetReader->getContainerName());
    for (auto jet : *m_trackJets) {
      BTagProps::isTagged.set(jet, static_cast<decltype(BTagProps::isTagged.get(jet))>(m_bTagTool->isTagged(*jet)));
      BTagProps::tagWeight.set(jet, Props::MV2c10.get(jet));
    }
  }

  if (m_jets) {
    m_bTagTool->setJetAuthor(m_jetReader->getContainerName());
    for (auto jet : *m_jets) {
      if(!m_use2DbTagCut) BTagProps::isTagged.set(jet, static_cast<decltype(BTagProps::isTagged.get(jet))>(m_bTagTool->isTagged(*jet)));
      else BTagProps::isTagged.set(jet, false);

      BTagProps::tagWeight.set(jet, Props::MV2c10.get(jet));

      if (m_isMC && !m_use2DbTagCut) BTagProps::eff.set(jet, m_bTagTool->getEfficiency(*jet));
      else BTagProps::eff.set(jet, -999.);
    }

  }
  //Temporarily removed treatment of fat jets in lephad, but do not remove we might need it back
  /*
 if (m_fatJets) {
    for (auto fatjet : *m_fatJets) {
      using FatJetType = typename std::remove_pointer<decltype(fatjet)>::type;
      static FatJetType::ConstAccessor<vector<ElementLink<DataVector<xAOD::IParticle>>>> 
        GhostAccessor ("GhostAntiKt2TrackJet");
      if (Props::isFatJet.get(fatjet)) {
        // determine number of b-tags in a fat jet
        //  number of tracks will be the track jets with closest DR with the ghost matched track jets
        auto nTags            = 0;
        std::vector<const xAOD::Jet*> trackJetsInFatJet;

        if (m_trackJets) {
          // loop over the ghost associated track jets
          for (auto gTrackParticle : GhostAccessor(*fatjet)) {
            if (!gTrackParticle.isValid()) continue;
            const auto &temp = **gTrackParticle;
            auto gTrackJet = dynamic_cast<const xAOD::Jet*>(&temp);
            auto minDR = 100.0;
            const xAOD::Jet *minDRTrackJet = nullptr;

            for (auto trackJet : *m_trackJets) {
              auto Deta    = trackJet->eta() - gTrackJet->eta(),
                   Dphi = trackJet->phi() - gTrackJet->phi(),
                   DR   = sqrt(Deta * Deta + Dphi * Dphi);

              if (minDR > DR) {
                minDR = DR;
                minDRTrackJet = trackJet;
              }
            }
            if (minDRTrackJet) trackJetsInFatJet.push_back(minDRTrackJet);
          }
        }

        for (auto trackJet : trackJetsInFatJet)
          if (BTagProps::isTagged.get(trackJet)) ++nTags;

        Props::nTrackJets.set(fatjet, trackJetsInFatJet.size());
        Props::nBTags.set(fatjet, nTags);

        // finished b-tagging criteria

        // determine number of true b-jets in a fat jet with simple DR
        auto nBJets = 0;

        if (m_isMC) {
          for (auto trackJet : trackJetsInFatJet) {
            auto label = -999;
            if (Props::HadronConeExclTruthLabelID.exists(trackJet)) label = Props::HadronConeExclTruthLabelID.get(trackJet);
            else if (Props::TruthLabelID.exists(trackJet)) label = Props::TruthLabelID.get(trackJet);
            if (label == 5) ++nBJets;
          }
        }
        else nBJets = -1;  // Data
        Props::nTrueBJets.set(fatjet, nBJets);

        // finished true b-jet association criteria
      }
    }
 }//end if loop of m_fatjets
*/

} // compute_btagging

// temporarily removed fatjet_selection and implemented a different way to do it, but please do not permentantly remove

void AnalysisReader_hhbbtt::compute_TRF_tagging (const std::vector<const xAOD::Jet*>& signalJets)
{
  bool ttag_track_jets=false;
  if (m_trackJets) {
    m_bTagTool->setJetAuthor(m_trackJetReader->getContainerName());
    if(ttag_track_jets)    m_bTagTool->truth_tag_jets(m_eventInfo->eventNumber(),*m_trackJets,m_config);
    else{
      for (auto jet : *m_trackJets) {
	BTagProps::isTagged.set(jet, static_cast<decltype(BTagProps::isTagged.get(jet))>(m_bTagTool->isTagged(*jet)));
	BTagProps::tagWeight.set(jet, Props::MV2c10.get(jet));
      }
    }
  }

  if (m_jets) {//we assume signal jets are taken from the larger jet container
    m_bTagTool->setJetAuthor(m_jetReader->getContainerName());
    m_bTagTool->truth_tag_jets(m_eventInfo->eventNumber(),signalJets,m_config);
  }

  if (m_fatJets) {
    for (auto fatjet : *m_fatJets) {
      Props::nTrackJets.set(fatjet, 0);
      Props::nBTags.set(fatjet, 0);
      if (m_isMC) Props::nTrueBJets.set(fatjet, 0);
      else Props::nTrueBJets.set(fatjet, -1);
    }
  }

} // compute_TRF_tagging


/*
void AnalysisReader_hhbbtt::fatjet_selection ()
{
  if (m_fatJets) {
    for (auto fatjet : *m_fatJets) {
      if (Props::isFatJet.get(fatjet)) {
        // remove fatjet close to leptons with a simple DR cut
        auto fatjeteta = fatjet->eta();
        auto fatjetphi = fatjet->phi();

        for (auto electron : *m_electrons) {
          auto Deta    = electron->eta() - fatjeteta,
               Dphi = electron->phi() - fatjetphi,
               DR   = sqrt(Deta * Deta + Dphi * Dphi);

          if (DR < 1.0) {
            Props::isFatJet.set(fatjet, static_cast<int>(false));
            break;
          }
        }

        if (!Props::isFatJet.get(fatjet)) continue;

        for (auto muon : *m_muons) {
          auto Deta    = muon->eta() - fatjeteta,
               Dphi = muon->phi() - fatjetphi,
               DR   = sqrt(Deta * Deta + Dphi * Dphi);

          if (DR < 1.0) {
            Props::isFatJet.set(fatjet, static_cast<int>(false));
            break;
          }
        }

	// if (!Props::isFatJet.get(fatjet)) continue;

        // finished lepton isolation criteria
      }
    }
  }
} // fatjet_selection

*/



 //the function below defines 1 btag and 2 btag regions based on how many btagged jets there are 
void AnalysisReader_hhbbtt::tagjet_selection(
					   std::vector<const xAOD::Jet*> signalJets,
					   std::vector<const xAOD::Jet*> forwardJets,
					   std::vector<const xAOD::Jet*> &selectedJets,
					   int &tagcatExcl,
					   int &tagcatIncl)
{
  string tagStrategy, tagAlgorithm;

  m_config->getif<std::string>("tagStrategy", tagStrategy);   // AllSignalJets,Leading2SignalJets,LeadingSignalJets
  m_config->getif<std::string>("tagAlgorithm", tagAlgorithm); // FlavLabel,FlavTag

  selectedJets.clear();
  tagcatExcl = -1;
  tagcatIncl = -1;

  /////////////////////////////////////////////////////////////
  // **B-Tagging Selection**
  bool Lead1BTag   = false;
  bool Lead2BTag   = false;
  int  Ind1BTag    = -1;
  int  Ind2BTag    = -1;
  int  nbtag       = 0;
  int  jetidx      = 0;
  int  nSignalJet  = signalJets.size();
  int  nForwardJet = forwardJets.size();
  const xAOD::Jet *Jet1;
  const xAOD::Jet *Jet2;
  const xAOD::Jet *Jet3;

  if (m_isMC && tagAlgorithm == "FlavLabel") //use truth label to select b-jets (MC only!)
  {
    if (Props::TruthLabelID.get(signalJets.at(0)) == 5) Lead1BTag = true;
    if (Props::TruthLabelID.get(signalJets.at(1)) == 5) Lead2BTag = true;
  }
  else // if (tagAlgorithm == "FlavTag") use b-tagging to select b-jets (the only option for data)
  {
    //truth tagging (MC only)
    // if(m_isMC && m_doTruthTagging){
    //   if (BTagProps::isTruthTagged.get(signalJets.at(0)) == 1) Lead1BTag = true;
    //   if (BTagProps::isTruthTagged.get(signalJets.at(1)) == 1) Lead2BTag = true;
    // }
   
    //if direct tagging
    // else{
    if (BTagProps::isTagged.get(signalJets.at(0)) == 1) Lead1BTag = true; //isTagged property is set in compute_btagging function
    if (BTagProps::isTagged.get(signalJets.at(1)) == 1) Lead2BTag = true;
      //}
  }

  //set ind1btag to btagged jet with highest pt, ind1btag to btagged jet with 2nd highest pt
  for (const xAOD::Jet *jet : signalJets)
  {
    if (m_isMC && tagAlgorithm == "FlavLabel") //use truth label to select b-jets (MC only!)
    {
      if (Props::TruthLabelID.get(jet) == 5)
      {
        nbtag++;

        if (Ind1BTag < 0) Ind1BTag = jetidx;
        else if (Ind2BTag < 0) Ind2BTag = jetidx;
      }
    }
    else //if (tagAlgorithm == "FlavTag") use b-tagging to select b-jets (the only option for data)
    {
      //if truth tagging (MC only)
      // if(m_isMC && m_doTruthTagging){
      // 	if (BTagProps::isTruthTagged.get(jet) == 1)
      // 	      {
      // 	         nbtag++;
	    
      // 	         if (Ind1BTag < 0) Ind1BTag = jetidx;
      // 	         else if (Ind2BTag < 0) Ind2BTag = jetidx;
      // 	      }
      // }
      //if direct tagging 
      // else{
	      if (BTagProps::isTagged.get(jet) == 1)
	      {
	         nbtag++;
	    
	         if (Ind1BTag < 0) Ind1BTag = jetidx;
	         else if (Ind2BTag < 0) Ind2BTag = jetidx;
	      }
	      //}
    }
    jetidx++;
  }

  /////////////////////////////////////////////////////////////
  // **Tag-Category Definition**
  //define regions based on weather you want 0/1/2 btag categories in regions with 2 or more signal jets (first tagStrategy), or if you want 0/1/2 bag categories in regions with exactly 2 signal jets 
  //not sure the difference between tagcatExcl =0/1/2 when 0/1/2 btags, tagcatIncl is for what? 

  if ((tagStrategy == "AllSignalJets") || (tagStrategy == "LeadingSignalJets"))
  {
    if (signalJets.at(0)->pt() / 1000. > 45.)
    {
      tagcatIncl = 0; 

      if (nbtag == 0) tagcatExcl = 0;
      Jet1 = signalJets.at(0);
      Jet2 = signalJets.at(1);

      if (nSignalJet >= 3) Jet3 = signalJets.at(2);
      else if (nForwardJet >= 1) Jet3 = forwardJets.at(0);
    }

    if ((nbtag >= 1) && (signalJets.at(Ind1BTag)->pt() / 1000. > 45.))
    {
      tagcatIncl = 1;

      if (nbtag == 1) tagcatExcl = 1;

      // if(tagStrategy == "LeadingSignalJets") { if (!Lead1BTag) return EL::StatusCode::SUCCESS; }
      if (tagStrategy == "LeadingSignalJets") if (!Lead1BTag) tagcatExcl = -1; Jet1 = signalJets.at(Ind1BTag);

      if (Ind1BTag == 0) Jet2 = signalJets.at(1);
      else Jet2 = signalJets.at(0);

      if (nSignalJet >= 3)
      {
        if (Ind1BTag < 2) Jet3 = signalJets.at(2);
        else Jet3 = signalJets.at(1);
      }
      else if (nForwardJet >= 1) Jet3 = forwardJets.at(0);
    }

    if ((nbtag >= 2) && (signalJets.at(Ind1BTag)->pt() / 1000. > 45.))
    {
      tagcatIncl = 2;

      if (nbtag == 2) tagcatExcl = 2;

      if (tagStrategy == "LeadingSignalJets") if (!Lead1BTag) tagcatExcl = -1; 
      Jet1 = signalJets.at(Ind1BTag);
      Jet2 = signalJets.at(Ind2BTag);
      int jetidx3 = 0;

      if (nSignalJet >= 3) // may change with if/else to be consistent with 1tag
      {
        for (const xAOD::Jet *jet : signalJets)
        {
          if ((jetidx3 != Ind1BTag) && (jetidx3 != Ind2BTag)) { Jet3 = jet; break; }
          jetidx3++;
        }
      }
      else if (nForwardJet >= 1) Jet3 = forwardJets.at(0);
    }
  }

  /////////////////////////////////////////////////////////////
  if (tagStrategy == "Leading2SignalJets")
  {
    tagcatIncl = 0;

    if ((Lead1BTag == false) && (Lead2BTag == false) && (signalJets.at(0)->pt() / 1000. > 45.)) tagcatExcl = 0;

    if (((Lead1BTag == true) || (Lead2BTag == true)) && (signalJets.at(0)->pt() / 1000. > 45.)) tagcatIncl = 1;

    if ((((Lead1BTag == true) && (Lead2BTag == false)) || ((Lead1BTag == false) && (Lead2BTag == true))) && (signalJets.at(0)->pt() / 1000. > 45.)) tagcatExcl = 1;

    if ((Lead1BTag == true) && (Lead2BTag == true) && (signalJets.at(0)->pt() / 1000. > 45.)) { tagcatExcl = 2; tagcatIncl = 2; }

    Jet1 = signalJets.at(0);
    Jet2 = signalJets.at(1);

    if (nSignalJet >= 3) Jet3 = signalJets.at(2);
    else if (nForwardJet >= 1) Jet3 = forwardJets.at(0);
  }

  selectedJets.push_back(Jet1);
  selectedJets.push_back(Jet2);

  //set is tagged PROPERTY AGAIN (over wrtting is btag in compute btagging?) for leading /sublead signal jets chosen in taggins trategy loop

  BTagProps::isTagged.set(Jet1, static_cast<decltype(BTagProps::isTagged.get(Jet1))>(m_bTagTool->isTagged(*Jet1)));
  BTagProps::isTagged.set(Jet2, static_cast<decltype(BTagProps::isTagged.get(Jet2))>(m_bTagTool->isTagged(*Jet2)));

  if ((nSignalJet >= 3) || (nForwardJet >= 1)) selectedJets.push_back(Jet3);

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
} // tagjet_selection

EL::StatusCode AnalysisReader_hhbbtt::fill_nJetHistos(std::vector<const xAOD::Jet*> jets, string jetType) {
  int  nJet  = jets.size();
  bool isSig = (jetType == "Sig");

  m_histSvc->BookFillHist("N" + jetType + "Jets", 11, -0.5, 10.5, nJet, m_weight);

  for (const xAOD::Jet *jet : jets) {
    m_histSvc->BookFillHist("Pt" + jetType + "Jets", 100,  0, 100, jet->pt() / 1e3, m_weight);
    m_histSvc->BookFillHist("Eta" + jetType + "Jets", 100, -5,   5, jet->eta(),    m_weight);
  }

  if (isSig) {
    for (const xAOD::Jet *jet : jets) fill_bTaggerHists(jet);
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader_hhbbtt::setevent_flavour(std::vector<const xAOD::Jet*> selectedJets )
{
  int jet0flav = -999;
  int jet1flav = -999;
  //if(Props::ConeTruthLabelID.exists(selectedJets.at(0))){ //in rel20
  if (Props::HadronConeExclTruthLabelID.exists(selectedJets.at(0))) { // in rel20
     //jet0flav = Props::ConeTruthLabelID.get(selectedJets.at(0));
     //jet1flav = Props::ConeTruthLabelID.get(selectedJets.at(1));
     jet0flav = Props::HadronConeExclTruthLabelID.get(selectedJets.at(0));
     jet1flav = Props::HadronConeExclTruthLabelID.get(selectedJets.at(1));
  }
  if(jet0flav < 0 || jet1flav < 0){
     Error("AnalysisReader_hhbbtt::setevent_flavour","Failed to retrieve jet flavour! Exiting.");
     return EL::StatusCode::FAILURE;
  }

  m_histNameSvc -> set_eventFlavour(jet0flav, jet1flav);

  return EL::StatusCode::SUCCESS;
}

/*float AnalysisReader_hhbbtt::computeBTagSFWeight (std::vector<const xAOD::Jet*> &signalJets)
{
  using CP::CorrectionCode;
  using CP::SystematicCode;
  using CP::SystematicSet;

  float weight = 1;

  std::unordered_map<std::string, float> btageffweights;
  auto systs = m_bTagTool->affectingSystematics();

  // Set all weights to 1
  for (auto var : systs) {
    btageffweights[var.name()] = 1.0;
  }

  for (auto jet : signalJets) {

    float sf(1.0);
    auto  GeV = 1e3;

    // Code crashes without this!
    if (fabs(jet->eta()) > 2.5) continue;
    if (fabs(jet->pt()) < 20. * GeV) continue;

    //std::cout << "Seg fault" << std::endl;

    CorrectionCode result;
    if (static_cast<bool>(BTagProps::isTagged.get(jet))) result = m_bTagTool->getBTaggingEfficiencyTool().getScaleFactor(*jet, sf);
    else result = m_bTagTool->getBTaggingEfficiencyTool().getInefficiencyScaleFactor(*jet, sf);
    weight *= sf;

    int label(1000);
    if (Props::HadronConeExclTruthLabelID.exists(jet)) label = Props::HadronConeExclTruthLabelID.get(jet);
    else if (Props::TruthLabelID.exists(jet)) label = Props::TruthLabelID.get(jet);

    if (result != CorrectionCode::Ok) Warning("computeBTagSFWeight", "Get efficiency failed (jet details): eta %f, pt %f, sf %f, MV2c10 %f, flavor %d", jet->eta(), jet->pt(), sf, Props::MV2c10.get(jet), label);

    for (auto var : systs) {
      SystematicSet set;
      set.insert(var);
      auto sresult = m_bTagTool->applySystematicVariation(set);

      if (sresult != SystematicCode::Ok) {
        //	std::cout << var.name() << " apply systematic variation FAILED " << std::endl;
      }

      if (static_cast<bool>(BTagProps::isTagged.get(jet))) result = m_bTagTool->getBTaggingEfficiencyTool().getScaleFactor(*jet, sf);
      else result = m_bTagTool->getBTaggingEfficiencyTool().getInefficiencyScaleFactor(*jet, sf);

      if (result != CorrectionCode::Ok) {
        // std::cout << var.name() << " getScaleFactor FAILED" << std::endl;
      } else {
        //	std::cout << var.name() << " " <<m_hbbname[var.name()] << sf << std::endl;
        btageffweights[var.name()] *= sf;
      }
      // Info("computeBTagSFWeight", "pT: %f, eta: %f, label: %d, SF: %f, var: %s", jet->pt(), jet->eta(), label, sf, var.name().c_str());
    }

    // don't forget to switch back off the systematics...
    SystematicSet defaultSet;
    auto dummyResult = m_bTagTool->applySystematicVariation(defaultSet);

    if (dummyResult != SystematicCode::Ok) Warning("computeBTagSFWeight", "Problem disabling systematics setting!");
  }
  if (m_currentVar == "Nominal") {
    for (auto pair : btageffweights)
      m_weightSysts.push_back({pair.first, pair.second/weight});
  }
  // Info("computeBTagSFWeight", "event weight: %f", weight);
  return weight;
} */ 
  


std::pair<int,int> AnalysisReader_hhbbtt :: HiggsTag(const xAOD::Jet *fatjet, std::string wp) {
  // LASER: A simplified version of the Higgs tagger due to current limitations in CxAOD Framework
  // Input:
  //       1) const xAOD::Jet *fatjet [fatjet to tag]
  //       2) std::string wp [working point]
  // Output:
  //       1) std::pair result [result.first = pass/fail tagger, result.second = number of btagged track jets]
  // Result codes look like:
  //      -1: wrong configuration, kill everything
  //      -2: couldn't get associated track jets
  //      -3: less than 2 associated track jets
  //      -4: something went wrong in Laser's logic
  //       0: a fatjet that did not pass the mass/substructure tagger cuts
  //       1: a fatjet that did pass the mass/substructure tagger cuts
  // NB: No muon-in-jet correction applied (due to aforementioned limitations)
  
  std::pair<int,int> result;
  
  // Step 0: check working point and define cut values
  float lowmasscut; // low mass window cut
  float highmasscut; // high mass window cut
  float D2_a0, D2_a1, D2_a2, D2_a3, D2_a4; // D2 fit parameters (4th order polynomial)
  if (wp == "loose") {
    lowmasscut = 76.;
    highmasscut = 146.;
  }
  else if (wp == "medium") {
    lowmasscut = 93.;
    highmasscut = 134.;
  }
  else if (wp == "tight") {
    lowmasscut = 93.;
    highmasscut = 134.;
    D2_a0 = 11.3684217088;
    D2_a1 = -0.0834101931325;
    D2_a2 = 0.000244968552399;
    D2_a3 = -3.09799473883e-07;
    D2_a4 = 1.44703493877e-10;
  }
  else {
    if(m_debug) std::cout << "ERROR!!! Provided working point is not defined!  It must be either loose, medium, or tight!" << std::endl;
    result.first = -1;
    return result;
  }
  // Step 1: get the associated track jets
  std::vector<const xAOD::Jet*> trkJets, goodTrkJets;
  if (fatjet->getAssociatedObjects<xAOD::Jet>("GhostAntiKt2TrackJet", trkJets)) {
    if(m_debug) std::cout << "Debug HiggsTagger: Step 1 in trk jets loos" << std::endl;
    int ntrkjets = 0;
    for (const xAOD::Jet *trkJ : trkJets) {
      if (!trkJ) continue; // if the trackjet is not valid then skip it
      if (!(trkJ->pt() / 1000. > 10. && fabs(trkJ->eta()) < 2.5)) continue;
      goodTrkJets.push_back(trkJ);
      ntrkjets++;
      m_histSvc->BookFillHist("HiggsTagger_trkJets_pt", 400, 0, 2000, trkJ->pt()/1e3, m_weight);
      m_histSvc->BookFillHist("HiggsTagger_trkJets_MV2c10", 200, -1, 1, Props::MV2c10.get(trkJ), m_weight);
      m_histSvc->BookFillHist("HiggsTagger_trkJets_2d", 400, 0, 2000, 200, -1, 1, trkJ->pt()/1e3, Props::MV2c10.get(trkJ), 1.);
    }
    m_histSvc->BookFillHist("HiggsTagger_trkJets_n", 10, 0, 10, ntrkjets, m_weight);
  }
  else {
    result.first = -2;
    return result; // no associated track jets
  }
  
  // Step 2: count the number of btagged associated track jets
  if (goodTrkJets.size() < 2) {
    result.first = -3;
    return result;
  }
  std::sort(goodTrkJets.begin(), goodTrkJets.end(), EventSelection::sort_pt);
  const xAOD::Jet *trkJ_lead = goodTrkJets.at(0);
  const xAOD::Jet *trkJ_sublead = goodTrkJets.at(1);
  int nbtaggedtrkjets = 0;
  if (Props::MV2c10.get(trkJ_lead) > m_trkBTagLimit) nbtaggedtrkjets++;
  if (Props::MV2c10.get(trkJ_sublead) > m_trkBTagLimit) nbtaggedtrkjets++;
  m_histSvc->BookFillHist("HiggsTagger_trkJets_nbtagged", 10, 0, 10, nbtaggedtrkjets, m_weight);
  result.second = nbtaggedtrkjets; // how many btagged track jets are there
  
  // Step 3: apply mass window cut
  std::string btagregion;
  if (result.second == 0) btagregion = "0btag"; // number of btagged track jets defines which region
  if (result.second == 1) btagregion = "1btag";
  if (result.second >= 2) btagregion = "2btag";
  std::cout << "NI: Debug HiggsTagger: Step 3 result should not be negative" << result.second << "btagregion  "<< btagregion  << std::endl;
  float fatjetmass = fatjet->m() / 1000.;
  float fatjetD2 = Props::D2.get(fatjet);
  float fatjetC2 = Props::C2.get(fatjet);
  m_histSvc->BookFillHist("HiggsTagger_premasscut_m", 100, 0, 500, fatjetmass, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_premasscut_d2", 100, 0, 5, fatjetD2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_premasscut_c2", 20, 0, 1, fatjetC2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_premasscut_m", 100, 0, 500, fatjetmass, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_premasscut_d2", 100, 0, 5, fatjetD2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_premasscut_c2", 20, 0, 1, fatjetC2, m_weight);
  if (!(fatjetmass > lowmasscut && fatjetmass < highmasscut)) {
    result.first = 0;
    return result; // an unsuccessfully tagged fatjet, pass to the control region
  }
  m_histSvc->BookFillHist("HiggsTagger_postmasscut_d2", 100, 0, 5, fatjetD2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_postmasscut_c2", 20, 0, 1, fatjetC2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_postmasscut_d2", 100, 0, 5, fatjetD2, m_weight);
  m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_postmasscut_c2", 20, 0, 1, fatjetC2, m_weight);
  

  // Step 4: apply substructure cut if we are at the tight working point
  if (wp != "tight") {
    result.first = 1;
    return result; // a succesfully tagged fatjet at the loose or medium working points
  }
  else {
    float D2cut = D2_a0 + (D2_a1 * fatjet->pt() / 1000.) + (D2_a2 * pow(fatjet->pt() / 1000., 2)) + (D2_a3 * pow(fatjet->pt() / 1000., 3)) + (D2_a4 * pow(fatjet->pt()/ 1000., 4));
    if (fatjetD2 > D2cut) {
      result.first = 1;
      m_histSvc->BookFillHist("HiggsTagger_postd2cut_c2", 20, 0, 1, fatjetC2, m_weight);
      m_histSvc->BookFillHist("HiggsTagger_"+btagregion+"_postd2cut_c2", 20, 0, 1, fatjetC2, m_weight);
      return result; // a successfully tagged fatjet at the tight working point
    }
    else {
      result.first = 0;
      return result; // an unsuccessfully tagged fatjet, pass to the control region
    }
  }
  
  // If you got here, then something went wrong
  result.first = -3;
  return result; // something went wrong, start yelling
}//end of HiggsTag

EL::StatusCode AnalysisReader_hhbbtt :: fill_lephad()
{
  // m_tree->SetVariation(m_currentVar);
  ResultHHbbtautau selectionResult = ((HHbbtautauLepHadSelection*)m_eventSelection)->result();

  
  //do final event selection
  const xAOD::Electron* electron = selectionResult.el;
  const xAOD::Muon* muon = selectionResult.mu;
  const xAOD::MissingET* met = selectionResult.met;
  std::vector<const xAOD::TauJet*> taus = selectionResult.taus; 
  std::vector<const xAOD::Jet*> signalJets = selectionResult.signalJets;
  std::vector<const xAOD::Jet*> forwardJets = selectionResult.forwardJets;
  //  std::cout << "AT START Event  " << m_eventInfo->eventNumber() << " num btag signal jets " << signalJets.size()  <<std::endl;

  int nSignalJet = signalJets.size();
  int nForwardJet = forwardJets.size();
  int nJet = nSignalJet + nForwardJet;

  //  m_weight *= Props::leptonSF.get(m_eventInfo);

  /////////////////////////////////////////////////////////////
  // **SF** : lepton SF : comment out trigger scale factors for now
  ////m_weight *= Props::leptonSF.get(m_eventInfo);
 
  /////////////////////////////////////////////////////////////
  // **Selection & SF**: trigger
  // trigger selection : nominal case
  //double triggerSF_nominal = 1.;
  //if ( !m_triggerTool->getTriggerDecision(m_eventInfo, triggerSF_nominal, electron, 0, muon, 0, met, 0, m_puReweightingTool->GetRandomRunNumber(m_eventInfo), "Nominal") ) return EL::StatusCode::SUCCESS;
  //if (m_isMC) m_weight *= triggerSF_nominal;
  // handle systematics
  //if (m_isMC && m_currentVar=="Nominal") {
  //  for (size_t i = 0; i < m_triggerSystList.size(); i++) {
      // not computing useless systematics
  //    if (electron && m_triggerSystList.at(i).find("MUON_EFF_Trig") !=std::string::npos) continue;
  //    if (muon && m_triggerSystList.at(i).find("ELECTRON_EFF_Trig") !=std::string::npos) continue;
      // get decision + weight
  //    double triggerSF = 1.;
  //    if ( !m_triggerTool->getTriggerDecision(m_eventInfo, triggerSF, electron, 0, muon, 0, met, 0, m_puReweightingTool->GetRandomRunNumber(m_eventInfo), m_triggerSystList.at(i)) ) return EL::StatusCode::SUCCESS;
  //    if (triggerSF_nominal>0) m_weightSysts.push_back({m_triggerSystList.at(i), triggerSF/triggerSF_nominal});
  //    else Error("fill_1Lep()", "Nominal trigger SF=0!, The systematics will not be generated.");
  //  }
  //}

  //leptons
  TLorentzVector lepVec;
  if (muon) {
    lepVec = muon->p4();
  } else if (electron) {
    lepVec = electron->p4();
  } else {
    Error("fill_lephad", "Missing lepton!");
    return EL::StatusCode::FAILURE;
  }
  //met
  TLorentzVector metVec;
  metVec.SetPtEtaPhiM(met->met(), 0, met->phi(), 0);

  // Taus
  int ntaus=taus.size();  
  if (ntaus==0) return EL::StatusCode::SUCCESS;
  
  TLorentzVector tauVec = taus.at(0)->p4();
  if (taus.size()>1){
    if(taus.at(0)->pt() < taus.at(1)->pt())tauVec = taus.at(1)->p4();
  }


  //  std::cout << " lepvec pt "  <<lepVec.Pt()/1e3 <<  " met  "  << metVec.Pt()/1e3 <<std::endl;
 
  // exit if we don't have at least 2 jets
  if (!(nJet >= 2)) return EL::StatusCode::SUCCESS;

  //fatjet_selection();
  compute_btagging(signalJets); //sets the isTagged and TagWeight property for signal jets, and eff property for MC signal jets

  fill_nJetHistos(signalJets, "Sig");
  fill_nJetHistos(forwardJets, "Fwd");
  
  // cut on nSignalJet, exit if not at least 2 signal jets
  if (!(nSignalJet >= 2)) return EL::StatusCode::SUCCESS; 
  

  // leading jet pt cut. sorting is done by event selection, exit if leading jet pt isn't >45 GeV
  if (signalJets.at(0)->pt()/1000. < 45.) return EL::StatusCode::SUCCESS; 

  // **Jets Definition** : define jet pair/triplet for anti-QCD cut // what?
  TLorentzVector j1VecPresel, j2VecPresel, j3VecPresel;
  j1VecPresel = signalJets.at(0)->p4();
  j2VecPresel = signalJets.at(1)->p4();

  if (nSignalJet  >= 3) j3VecPresel = signalJets.at(2)->p4();
  else if (nForwardJet >= 1) j3VecPresel = forwardJets.at(0)->p4(); 

  // **Jets Definition** : jet pair for mBB (+ possibly third jet)
  TLorentzVector j1Vec, j2Vec, j3Vec;

  // **Jets Definition** : jet pair for mBB rescaling
  TLorentzVector j1VecRW,j2VecRW;


  //Nina commented out this chunk for debugging
  
   // **B-tagging-Jet Selection** : 
  int tagcatExcl=-1;
  int tagcatIncl=-1;
  std::vector<const xAOD::Jet*> selectedJets; // here selected jets are all jets
  selectedJets.clear(); 
  //here we are just identifiying 0/1/2 btagged categories in regions that have either exactly 2 or >2 jets as specified by tagStraegy


  //after this step our selected jets are chosen to be our Signal jets that are btagged based on tag stragy! 
  //std::cout << " STEP 4: about to define tagjet_selections  " << std::endl; 
  tagjet_selection(signalJets, forwardJets, selectedJets, tagcatExcl, tagcatIncl);
  
  if(tagcatExcl==-1) return EL::StatusCode::SUCCESS; // Here if tagcatExcl is not set the code exits
  
  

  //std::cout << " STEP 5: defined tagjet_selection  " << std::endl;

  if (m_isMC) m_weight *= computeBTagSFWeight(signalJets); // leave out Btag SF for now
  //if truth tagging

  
  
  if(m_isMC && m_doTruthTagging){
    //only 2tag for the moment
    if(tagcatExcl!=2) return EL::StatusCode::SUCCESS;
    m_weight *= BTagProps::truthTagEventWeight.get(m_eventInfo); //NI, we dont have truthTagEventWeight property why do we need this? its not in VHbb boosted on mono_VH, it is in VHbb
  }
 
  
  m_histNameSvc->set_nTag(tagcatExcl); //this means we make histograms for 0/1/2 tagged regions
  //  m_histNameSvc->set_nJet(nJet); // the tag strategy for now is set to tag 2 signal jets

  // tagcatIncl is not used for now
  /*
    tagcatIncl==0 --> 0ptag, !1ptag, !2ptag
    tagcatIncl==1 --> 0ptag, 1ptag, !2ptag
    tagcatIncl==2 --> 0ptag, 1ptag, 2ptag
   */
  
  //Nina added temporarily for code check 
  //j1Vec = signalJets.at(0)->p4(); 
  //j2Vec = signalJets.at(1)->p4();
  // if(signalJets.size() >= 3) j3Vec = signalJets.at(2)->p4(); //3rd jet

  j1Vec = selectedJets.at(0)->p4(); // j1Vec is leading jet in whole event, it was btagged in tagjet_selection
  j2Vec = selectedJets.at(1)->p4(); // j1vec is subleading jet in whole event

  if(selectedJets.size() >= 3) j3Vec = selectedJets.at(2)->p4(); //3rd jet

  if( nJet>=4 && m_fillCR ){
     m_histNameSvc->set_description("topCR"); //top CR
  }
   
  TLorentzVector HbbVec = j1Vec + j2Vec;
  TLorentzVector Hlephad = tauVec + lepVec;

  //  if (Hlephad.M()/1e3 < 5.){ // if Hlephad Mass is less than 5 GeV
  //m_histSvc -> BookFillHist("ptDiff",   100, 0,20,fabs(lepVec.Pt()-tauVec.Pt())/1e3, m_weight); 
  //m_histSvc -> BookFillHist("DRlephad",   100, 0,0.2,lepVec.DeltaR(tauVec), m_weight);
  // if(muon)m_histSvc -> BookFillHist("muORe",   2, 0,1 ,0, m_weight);
  //if(electron)m_histSvc -> BookFillHist("muORe",   2, 0,1 ,1, m_weight);
  //} 
 
  // **Selection** : mbb invariant mass window + mbb-sidebands
  // rescaled quantities are : j1VecRW, j2VecRW
  // double mbb_weight = 1;
  // int mbbRegion=0;
  // j1VecRW = j1Vec;
  // j2VecRW = j2Vec;
  // if(mbbwindow) 
  //   {
  //     if(HVec.M() < 95e3 || HVec.M() > 140e3) mbbRegion=1;
  //     else 
  // 	{  
  // 	  mbb_weight = 125.0 / (HVec.M()*1e-3); // Jets reweighting
  // 	  j1VecRW = (j1Vec)*mbb_weight;
  // 	  j2VecRW = (j2Vec)*mbb_weight;	 
  // 	}
  //   } 
  // if(mbbRegion==1 && m_histNameSvc->get_description()=="SR") m_histNameSvc->set_description("mBBcr");
  // else if (mbbRegion==1 && m_histNameSvc->get_description()=="topCR") return EL::StatusCode::SUCCESS;
  /*
  // reset MVA tree variables
  m_tree->Reset();
  
  m_tree->EventWeight = m_weight;
  m_tree->EventNumber = m_eventInfo->eventNumber();
  
  m_tree->dPhiLBmin = fabs( lepVec.DeltaPhi(j1Vec) );
  if (m_tree->dPhiLBmin > fabs( lepVec.DeltaPhi(j2Vec) )) {
    m_tree->dPhiLBmin = fabs( lepVec.DeltaPhi(j2Vec) );
  }
  m_tree->MET = met->met();
  */

  double bdt_pTL = lepVec.Pt();
  
  // fill jet histos and set nJet, flavour for histo names
  //fill_jetHistos(signalJets, forwardJets);
  ////fill_jetSelectedHistos(signalJets, forwardJets, selectedJets); //This line does not work, must debug!!
 
  //  m_tree->sample = m_histNameSvc->getFullSample();

  //HbbVec

  //  std::cout << " MyReader Event  " << m_eventInfo->eventNumber() << " weight before histograms are drawn "<< m_weight  << std::endl; 

  m_histSvc -> BookFillHist("HbbPt",   30, 0, 1000, HbbVec.Pt()/1e3, m_weight);
  m_histSvc -> BookFillHist("HbbM",   30, 0, 500, HbbVec.M()/1e3, m_weight);
  m_histSvc -> BookFillHist("HbbEta",   20, -5, 5, HbbVec.Eta(), m_weight);
  m_histSvc -> BookFillHist("HbbPhi",   20, -4, 4, HbbVec.Phi(), m_weight);

  //m_histSvc -> BookFillHist("dPhiVBB", 100, 0, 3.15, m_tree->dPhiVBB, m_weight);
  //m_histSvc -> BookFillHist("MET",   100, 0, 500, m_tree->MET/1e3, m_weight);
  //m_histSvc -> BookFillHist("pTV",   100, 0, 500, m_tree->pTV/1e3, m_weight);
  m_histSvc -> BookFillHist("pTL",   100, 0, 500, bdt_pTL/1e3, m_weight);

  //lephad vector
  m_histSvc -> BookFillHist("HlephadPt",   30, 0, 1000, Hlephad.Pt()/1e3, m_weight);
  m_histSvc -> BookFillHist("HlephadM",   30, 0, 500, Hlephad.M()/1e3, m_weight);
  m_histSvc -> BookFillHist("HlephadEta",   20, -5, 5, Hlephad.Eta(), m_weight);
  m_histSvc -> BookFillHist("HlephadPhi",   20, -4, 4, Hlephad.Phi(), m_weight);

  //lep vector
  m_histSvc -> BookFillHist("LepPt",   30, 0, 600, lepVec.Pt()/1e3, m_weight);
  m_histSvc -> BookFillHist("LepM",   10, 0, 1, lepVec.M()/1e3, m_weight);
  m_histSvc -> BookFillHist("LepEta",   20, -4, 4, lepVec.Eta(), m_weight);
  m_histSvc -> BookFillHist("LepPhi",   20, -4, 4, lepVec.Phi(), m_weight);

  //tau vector
  m_histSvc -> BookFillHist("TauPt",   30, 0, 700, tauVec.Pt()/1e3, m_weight);
  m_histSvc -> BookFillHist("TauM",   10, 0, 2, tauVec.M()/1e3, m_weight);
  m_histSvc -> BookFillHist("TauEta",   20, -3, 3, tauVec.Eta(), m_weight);
  m_histSvc -> BookFillHist("TauPhi",   20, -4, 4, tauVec.Phi(), m_weight);

  // lep and tau together
  m_histSvc -> BookFillHist("ptDiff",   20, 0,200,fabs(lepVec.Pt()-tauVec.Pt())/1e3, m_weight);
  m_histSvc -> BookFillHist("DRlephad",   5, 0, 5.0,lepVec.DeltaR(tauVec), m_weight);
  if(muon)m_histSvc -> BookFillHist("muORe",   2, 0,1 ,0, m_weight);
  if(electron)m_histSvc -> BookFillHist("muORe",   2, 0,1 ,1, m_weight);
 

  m_histSvc->BookFillHist("met", 30, 0, 100, metVec.Pt(), m_weight);


  // fill MVA tree
  //  if (m_histNameSvc->get_isNominal()) {    m_tree->Fill();
  //}

  return EL::StatusCode::SUCCESS;
} 



EL::StatusCode AnalysisReader_hhbbtt :: fill_hadhad()
{

  if(m_debug) std::cout << "FILLHADHAD: " << m_weight << " m_isMC " << m_isMC << std::endl;

  m_histNameSvc->set_description("PreSel");
  m_histNameSvc->set_nTag(0);

  if(m_isMC && fabs(Props::MCEventWeight.get(m_eventInfo)) > 10.){
    Info("execute ()", "Warning: mc weight too high %f",m_weight);
    return EL::StatusCode::SUCCESS;
  }

  //std::cout << "NowRunningOn: " << m_eventInfo->eventNumber() << " pileuprndRunNumber: " << Props::PileupRdmRun.get(m_eventInfo) << " randomrunumber: " << Props::RandomRunNumber.get(m_eventInfo) << st$

  bool passSel=true;
  bool passPreSel = true;

  SelectionContainers container;
  container.evtinfo   = m_eventInfo;
  container.met       = m_met;
  container.electrons = m_electrons;
  container.photons   = m_photons;
  container.muons     = m_muons;
  container.taus      = m_taus;
  container.jets      = m_jets;
  container.fatjets   = m_fatJets;
  container.trackjets = m_trackJets;
  //container.truthParticles = m_truthParts;

  //passPreSel &= m_eventSelection->passPreSelection(container,true);
  //passSel    &= m_eventSelection->passSelection(container,false);

  //if(!passSel || !passPreSel) return EL::StatusCode::SUCCESS;

  ResultHHbbtautau selectionResult = ((HHbbtautauHadHadSelection*)m_eventSelection)->result();
  ((HHbbtautauHadHadJetSelection*)m_eventPostSelection)->setResult(selectionResult);

  passPreSel &= m_eventPostSelection->passPreSelection(container,false);
  passSel &= m_eventPostSelection->passSelection(container,false);

  if(!passPreSel || !passSel) return EL::StatusCode::SUCCESS;

  ResultHHbbtautau postSelectionResult = ((HHbbtautauHadHadJetSelection*)m_eventPostSelection)->result();

  //fill after selection with final objects
  std::vector<const xAOD::TauJet*> tauvector = postSelectionResult.taus;
  std::vector<const xAOD::Jet*> jetvector = postSelectionResult.signalJets;
  std::vector<const xAOD::Jet*> fjetvector = postSelectionResult.forwardJets;
  const xAOD::MissingET* met=postSelectionResult.met;
  Trigger TriggerDecision=postSelectionResult.trigger;

  bool doTriggerStudy = false;
  if(doTriggerStudy) {
    if(Props::passHLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo.get(m_eventInfo)) m_histNameSvc->set_description("Trigger_DTT");
    fill_triggerPlots();
  }

  //-------------------------------------------------------------------------------------------------------

  if(m_debug) std::cout << "Total Number of Taus " << tauvector.size()<<std::endl;

 
  double TriggerSF=1;
  TriggerSF= Props::trigSF.get(m_eventInfo);
  
  if(m_debug) std::cout<<"TriggerSF "<<TriggerSF<<std::endl;

  m_weight*=TriggerSF;// apply trigger SF

  std::vector<const xAOD::Jet*> sjets;//vector containing all the Signal Jets
  std::vector<const xAOD::TauJet*> signalTaus;//vector containing all the Signal Taus
  std::vector<const xAOD::TauJet*> antiTaus;//vector containing all the AntiTaus

  int NumberSignalJet45=0;
  int NumberSignalJet50=0;
  int NumberSignalJet70=0;
  int NumberSignalJet80=0;

  for(unsigned int i=0; i<jetvector.size(); i++){
    if(!Props::isSignalJet.get(jetvector[i])) continue;
    if(jetvector[i]->pt()<=20e+3) continue;
    sjets.push_back(jetvector[i]);//fill vector of final signal jets with pt 20 cut
    if(jetvector[i]->pt()<=45e+3) continue;
    NumberSignalJet45++;
    if(jetvector[i]->pt()<=50e+3) continue;
    NumberSignalJet50++;
    if(jetvector[i]->pt()<=70e+3) continue;
    NumberSignalJet70++;
    if(jetvector[i]->pt()<=80e+3) continue;
    NumberSignalJet80++;
    if(m_debug) std::cout << "jet pt: " << jetvector[i]->pt() << std::endl;
  }
   
  if(m_debug) std::cout << "Number of signal Jets: " <<sjets.size()<< std::endl;

  int NumberSignalTau30=0;
  int NumberSignalTau40=0;
  int NumberSignalTau100=0;
  int NumberSignalTau140=0;
  int NumberSignalTau160=0;
  int NumberSignalTau180=0;
  int NumberSignalTau200=0; 

  for(unsigned int i=0; i<tauvector.size(); i++){
    if(!Props::isBDTMedium.get(tauvector[i])) continue;
    if(tauvector[i]->pt()<=20e+3) continue;
    signalTaus.push_back(tauvector[i]);//fill vector of final signal taus with pt 20 cut
    if(tauvector[i]->pt()<=30e+3) continue;
    NumberSignalTau30++;
    if(tauvector[i]->pt()<=40e+3) continue;
    NumberSignalTau40++;
    if(tauvector[i]->pt()<=100e+3) continue;
    NumberSignalTau100++;
    if(tauvector[i]->pt()<=140e+3) continue;
    NumberSignalTau140++;
    if(tauvector[i]->pt()<=160e+3) continue;
    NumberSignalTau160++;
    if(tauvector[i]->pt()<=180e+3) continue;
    NumberSignalTau180++;
    if(tauvector[i]->pt()<=200e+3) continue;
    NumberSignalTau200++;
  }

  if(m_debug) std::cout << "Number of signal BDTMedium Taus: " <<signalTaus.size()<< std::endl;

  for(unsigned int i=0; i<tauvector.size(); i++){
    if(!Props::isAntiTau.get(tauvector[i])) continue;
    if(!Props::isHHRandomTau.get(tauvector[i])) continue;
    if(tauvector[i]->pt()<=20e+3) continue;
    antiTaus.push_back(tauvector[i]);//fill vector of final antitaus with pt 20 cut
  }
  
  if(m_debug) std::cout << "Number of AntiTaus: " <<antiTaus.size()<< std::endl;

  if(m_debug) std::cout << "NumberSignalJet80: " << NumberSignalJet80 << " TriggerDecision: " << TriggerDecision << std::endl;

  m_PassTriggerJetSelection = false;
  if(TriggerDecision==eDTT || TriggerDecision==eSTT){
    m_PassTriggerSelection=true;
    if((TriggerDecision==eDTT && NumberSignalJet80>0) || (TriggerDecision==eSTT && NumberSignalJet45>0)){
      m_PassTriggerJetSelection=true;
    }
  }

  //btagging
  if(m_doTruthTagging && m_isMC){
    compute_TRF_tagging(sjets);
  }
  else{
    compute_btagging(sjets);// compute the b-tagging and assign the isTagged property for all the signal jets with pT>20 GeV
  }

  /*if(m_isMC){
    btagWeight=computeBTagSFWeight(sjets, m_jetReader->getContainerName());
    if(m_debug) std::cout << "BTagWeight: "<<btagWeight<< std::endl;

    if(m_doTruthTagging){
      BTagProps::truthTagEventWeight.set(m_eventInfo, btagWeight);
    }
  }

  m_weight *= btagWeight;// apply btag weight*/

  m_nBJets=0;
  std::vector<int> btag_index;
  for(unsigned int i=0; i<sjets.size();i++){
    if(BTagProps::isTagged.get(sjets.at(i))) {
      btag_index.push_back(i);
      m_nBJets++;
    }
  }

  m_PassLeadingJetPtCut=false;
  if(m_nBJets == 0 && sjets.at(0)->p4().Pt()>45e+3){
    m_PassLeadingJetPtCut = true;
  }else if(m_nBJets == 1 && sjets.at(0)->p4().Pt()>45e+3){
    if(btag_index.size() != 1 && m_debug) std::cout << "Something weird is going on with btagging (1tag) " << m_nBJets << std::endl;
    m_PassLeadingJetPtCut = true;
  }else if(m_nBJets >=2 && sjets.at(btag_index.at(0))->p4().Pt()>45e+3){
    if(m_debug) std::cout << "lead b-tag jet pt " << sjets.at(btag_index.at(0))->p4().Pt() << std::endl;
    m_PassLeadingJetPtCut = true;
  }

  if(m_debug) std::cout << "MMC fit status: " << Props::mmc_fit_status.get(m_mmc) << std::endl;
  //starting defining variables for the analysis
  if(Props::mmc_fit_status.get(m_mmc) == 1) m_MMCVec.SetPtEtaPhiM(Props::mmc_mlnu3p_4vect_pt.get(m_mmc),Props::mmc_mlnu3p_4vect_eta.get(m_mmc),Props::mmc_mlnu3p_4vect_phi.get(m_mmc),Props::mmc_mlnu3p_4vect_m.get(m_mmc));
  else if (Props::mmc_fit_status.get(m_mmc) == 0) m_MMCVec.SetPtEtaPhiM(-1,-1,-1,-1);
  if(m_debug) std::cout << "MMC " << m_MMCVec.M() << " nBJets " << m_nBJets << std::endl;

  m_METVec.SetPtEtaPhiM(met->met(),0,met->phi(),0);

  std::vector<const xAOD::TauJet*> finalTaus;
  m_region="";

  if(signalTaus.size()==2){
    finalTaus.push_back(signalTaus.at(0));
    finalTaus.push_back(signalTaus.at(1));
    m_region = "MM";
  }else if(signalTaus.size()==1 && antiTaus.size()==1){
    //Fill pT ordered finalTaus vector
    if(signalTaus.at(0)->pt()>antiTaus.at(0)->pt()){
      m_region = "MA";
      finalTaus.push_back(signalTaus.at(0));
      finalTaus.push_back(antiTaus.at(0));
    }else if(signalTaus.at(0)->pt()< antiTaus.at(0)->pt()){
      m_region = "AM";
      finalTaus.push_back(antiTaus.at(0));
      finalTaus.push_back(signalTaus.at(0));
    }
  }else if(antiTaus.size()==2){
    finalTaus.push_back(antiTaus.at(0));
    finalTaus.push_back(antiTaus.at(1));
    m_region = "AA";
  }

  m_histNameSvc->set_analysisType(HistNameSvc::AnalysisType::HHres);

  int pair_charge=Props::charge.get(finalTaus.at(0))*Props::charge.get(finalTaus.at(1));

  if(pair_charge<0){
    m_taupairCharge = "OS";
    m_region +="_"+m_taupairCharge;
  }
  else if(pair_charge>0) {
    m_taupairCharge = "SS";
    m_region += "_"+m_taupairCharge;
  }

  m_RegA = false;
  m_RegB = false;
  m_RegC = false;
  m_RegD = false;

  //Dictionary for regions conventions (for cutflow only)
  if(m_region=="MM_OS") m_RegA=true;
  if(m_region=="MM_SS") m_RegB=true;
  if(m_region=="AA_OS" || m_region=="MA_OS" || m_region=="AM_OS") m_RegC=true;
  if(m_region=="AA_SS" || m_region=="MA_SS" || m_region=="AM_SS") m_RegD=true;

  m_PassTauPtSelection = true;
  if(TriggerDecision==eDTT){//pt cut for 3p taus in case of DTT
    if(Props::nTracks.get(finalTaus.at(0))==3){
       if(finalTaus.at(0)->pt()/1e3<=50){
          m_PassTauPtSelection=false;
       }
    }
    if(Props::nTracks.get(finalTaus.at(1))==3){
       if(finalTaus.at(1)->pt()/1e3<=40){
          m_PassTauPtSelection=false;
       }
    }
  }
	
  // **Calculate** : b-tagging SF
  if (m_isMC) {

    float btagWeight=1;

    btagWeight = computeBTagSFWeight(sjets, m_jetReader->getContainerName());
    if(m_debug) std::cout << "BTagWeight: " << btagWeight << std::endl;

    if(m_doTruthTagging){
      BTagProps::truthTagEventWeight.set(m_eventInfo,btagWeight);
      if(m_debug) std::cout << "Truth tag event weight: " << BTagProps::truthTagEventWeight.get(m_eventInfo) << std::endl;
    }

    m_weight  *= btagWeight;
    //setevent_flavour(sjets); //Working

    float tauEff = Props::effSF.get(finalTaus.at(0))*Props::effSF.get(finalTaus.at(1));
    if(m_debug) std::cout << "Tau efficiency weight: " << tauEff << std::endl;
    m_weight *= tauEff;

    m_weight *=Props::trigSF.get(m_eventInfo);

  }

  bool PassCutFlowPres=fill_cutFlow(sjets,finalTaus);
  if(!PassCutFlowPres) return EL::StatusCode::SUCCESS;

  if(m_debug) std::cout << "Passed cut flow for region: " << m_region << " nTags: " << m_nBJets << std::endl;

  m_histNameSvc->set_nTag(m_nBJets);

  int reg_number = -1;

  if(m_region=="MM_SS"){
    reg_number=7.5;
  }else if(m_region=="MA_SS"){
    reg_number=6.5;
  }else if(m_region=="AM_SS"){
    reg_number=5.5;
  }else if(m_region=="AA_SS"){
    reg_number=4.5;
  }else if(m_region=="MM_OS"){
    reg_number=3.5;
  }else if(m_region=="MA_OS"){
    reg_number=2.5;
  }else if(m_region=="AM_OS"){
    reg_number=1.5;
  }else if(m_region=="AA_OS"){
    reg_number=0.5;
  }

  int sample_number=getSampleNumber();

  //std::cout << "sample_number: " << sample_number << std::endl;

  if(reg_number!=-1) m_histSvc->BookFillHist("2DEvntCat",20,0,20,8,0,8,sample_number,reg_number,m_weight);

  m_histNameSvc->set_description(m_region);

  std::string regionFF = get2DProngRegion(finalTaus);

  if(finalTaus.at(0)->pt()*1e-3<100 && finalTaus.at(1)->pt()*1e-3>100 && regionFF == "nProng11") std::cout << "Event found !" << std::endl;

  //applying QCD FF from AA_OS template for estimating QCD in signal region
  //and AA_SS is used to estimate the closure
  m_config->getif<bool>("applyFF", m_applyFF);
  if( (m_RegD || m_RegC ) && m_applyFF){
    EL_CHECK("AnalysisReader::fillhadhad()",applyFF(regionFF,finalTaus));
  }

  if(m_writeMVATree) m_tree->Reset();

  fill_xHistos(sjets,finalTaus);
  fill_jetHistos(jetvector);
  fill_tauHistos(tauvector);
  fill_sjetHistos(sjets);

  if(m_getFFInputs) fill_tauFFHistos(regionFF,finalTaus);

  if(m_readMVA){
     m_tree->ReadMVA();
     // fill histos with BDT output
     m_histSvc->BookFillHist("BDT", 100, -1, 1, m_tree->BDT, m_weight);
  }

  //fill MVA tree
  if(m_histNameSvc->get_isNominal() && m_writeMVATree) m_tree->Fill();

  //Apply SR CUTS
  bool passMbbCut = m_bb.M()>80 && m_bb.M()<135;
  bool passMMCCut = m_MMCVec.M()>85 && m_MMCVec.M()<135;
  if(m_region=="MM_OS" && passMbbCut && passMMCCut) fill_SRPlots();

  //Apply ZCR Cuts
  if(m_nBJets==1 && m_MMCVec.M()) fill_SRPlots();

  if(m_debug) std::cout << "Exiting successfully" << std::endl;

  return EL::StatusCode::SUCCESS;

}// end of fill had had

// Fill boosted analysis plots
EL::StatusCode AnalysisReader_hhbbtt :: fill_boosted()
{

  std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> dRMap;
  std::vector<const xAOD::DiTauJet*> my_ditauJets;

  // Get Objects being used in the analysis
  EL_CHECK("AnalysisReader::fill_boosted()",getObjects(my_ditauJets));

  if(m_debug) std::cout << "MyDiTauJets size: " << my_ditauJets.size() << std::endl;
  // Skip the event if there are zero fatjets and ditaujets
  if(my_ditauJets.size() == 0 || m_my_fatJets.size() < 2 ) return EL::StatusCode::SUCCESS;

  // Truth matching function for optimization studies
  if(m_isMC) fillTruthEffPlots(my_ditauJets);

  // Set description to "ALL" events in CxAOD
  m_histNameSvc->set_description("All");

  // Find number of BTagged track jets in fat jets
  findNBTagTJ();

  // Fill fatjet plots and variables
  fillFatJetPlots();

  // Fill ditaujet plots and variables
  fillDiTauJetPlots(my_ditauJets);

  // Fill fatJet/diTauJet dRMap
  dRMap = fillDRMap(my_ditauJets);

  // Get selected DiTauJet
  getSelDT(my_ditauJets);
  // Get selected FatJet
  getSelFJ(dRMap);

  //if(m_nBTagTJ==-1) std::cout << "m_nBTagTJ==-1" << std::endl;
  // If not selected FJ and DiTau are found return
  if(m_MatchFJets.empty()) return EL::StatusCode::SUCCESS;
  if(m_SelFJ == NULL || m_SelDT == NULL || m_nBTagTJ == -1) return EL::StatusCode::SUCCESS;

  if (m_isMC) m_weight *= Props::effSF.get(m_SelDT);

  // Fill full system variables
  fillFullSystemVars(dRMap);
  fillFullSystemPlots();

  m_histNameSvc->set_nTag(m_nBTagTJ);

  bool passSideBand= applySideBandSelection();
  bool passPreSel  = applyPreSelection();
  bool passSR      = applySRSelection();
  bool passQCDCR   = applyQCDSelection();
  bool passTTBarCR = applyTTBarSelection();

  if(passSR) fillSRPlots();
  else if(passQCDCR) fillQCDCRPlots();
  else if(passSideBand) fillSideBandPlots();
  else if(passTTBarCR) fillTTBarCRPlots();

  //if((passSR && passQCDCR) || (passSR && passTTBarCR) || (passQCDCR && passTTBarCR)) std::cout << "Warning SR and background CRs are not orthogonal" << std::endl;

  return EL::StatusCode::SUCCESS;


}//end of boosted

EL::StatusCode AnalysisReader_hhbbtt::getObjects(std::vector<const xAOD::DiTauJet*> &my_ditauJets)
{

  //Info("getObjects()", "On getObjects function '%s'.", m_analysisType.c_str());

  bool passSel=true;
  bool passPreSel = true;

  SelectionContainers container;
  container.evtinfo        = m_eventInfo;
  container.met            = m_met;
  container.electrons      = m_electrons;
  container.photons        = m_photons;
  container.muons          = m_muons;
  container.taus           = m_taus;
  container.jets           = m_jets;
  container.fatjets        = m_fatJets;
  container.trackjets      = m_trackJets;
  container.ditaus         = m_ditaus;
  container.truthParticles = m_truthParts;
  container.truthTaus      = m_truthParts;

  if(!m_truthTaus && m_isMC) std::cout<< "Warning: Truth taus not found" << std::endl;

  ResultHHbbtautau selectionResult = ((HHbbtautauBoostedSelection*)m_eventSelection)->result();
  ((HHbbtautauBoostedJetSelection*)m_eventPostSelection)->setResult(selectionResult);

  passPreSel &= m_eventPostSelection->passPreSelection(container,false);
  passSel &= m_eventPostSelection->passSelection(container,false);

  if(!passPreSel || !passSel) return EL::StatusCode::SUCCESS;

  if(m_debug){
     std::cout << Props::passHLT_j420_a10_lcw_L1J100.get(m_eventInfo) 
        << " " << Props::passHLT_j360_a10_lcw_L1J100.get(m_eventInfo) 
        << " " << Props::passHLT_j400_a10_lcw_L1J100.get(m_eventInfo) 
        << " " << Props::passHLT_j420_a10r_L1J100.get(m_eventInfo)
        << " " << Props::passHLT_j360_a10r_L1J100.get(m_eventInfo) 
        << " " << Props::passHLT_j400_a10r_L1J100.get(m_eventInfo) << std::endl;
  }

  ResultHHbbtautau postSelectionResult = ((HHbbtautauBoostedJetSelection*)m_eventPostSelection)->result();

  std::vector<const xAOD::Jet*>      fatJets  = postSelectionResult.fatJets; //includes 200 GeV cut, eta <2.4, mass >50 GeV
  //std::vector<const xAOD::Jet*>      jets     = postSelectionResult.signalJets;
  std::vector<const xAOD::DiTauJet*> diTaus   = postSelectionResult.ditaus;
  m_my_met      = postSelectionResult.met;

  m_METVec.SetPtEtaPhiM(m_my_met->met(),0,m_my_met->phi(),0);
  if(m_debug) std::cout << "MET: " << m_METVec.Pt()*1e-3 << " " << std::endl;

  my_ditauJets.clear();
  for(auto ditau:diTaus){
     if((ditau->pt()*1e-3>300 && ditau->pt()*1e-3<1500) && fabs(ditau->eta())<2.0) my_ditauJets.push_back(ditau);
  }

  m_my_fatJets.clear();
  for(auto fatjet:fatJets){
     if((fatjet->pt()*1e-3>300 && fatjet->pt()*1e-3<1500) && fatjet->m()*1e-3>50 && fabs(fatjet->eta())<2.0) m_my_fatJets.push_back(fatjet);
  }


  return EL::StatusCode::SUCCESS;
}

void AnalysisReader_hhbbtt::fillTruthEffPlots(std::vector<const xAOD::DiTauJet*> my_ditauJets)
{

  if(m_debug) std::cout << "Truth particles object size: " << m_truthTaus->size() << std::endl;

  if(m_truthTaus->size() ==0){
      //std::cout << "Error: truth tau vector cannot be larger than 2, or could be ?" << std::endl;
      return;
  }

  xAOD::TruthParticleContainer::const_iterator thItr = m_truthTaus->begin();
  xAOD::TruthParticleContainer::const_iterator thItrE = m_truthTaus->end();

  std::vector<const xAOD::TruthParticle*> taus;
  std::vector<const xAOD::TruthParticle*> bis;

  TLorentzVector ditau_vis4V(0,0,0,0);

  for( ; thItr != thItrE; ++thItr ){
    const xAOD::TruthParticle *truthP = *thItr;
    if(Props::TruthTauOrigin.exists(truthP)){
      //std::cout << "TruthTauOrigin: " << Props::TruthTauOrigin.get(truthP) << std::endl;
      if(Props::TruthTauOrigin.get(truthP) != 14) continue;
    }else{
      //std::cout << "Warning:variable TruthTauOrigin doesn't exist, check how are you filling it" << std::endl;
      continue;
    }
    TLorentzVector tau_vis4V;
    tau_vis4V.SetPtEtaPhiM(Props::pt_vis.get(truthP),Props::eta_vis.get(truthP),Props::phi_vis.get(truthP),Props::m_vis.get(truthP));
    ditau_vis4V += tau_vis4V;
  }

  if(m_debug) std::cout << "Truth DiTau: " << " pt: " << ditau_vis4V.Pt() << " eta: " << ditau_vis4V.Eta() << " phi: " << ditau_vis4V.Phi() << std::endl;

  for (auto ditaujet : my_ditauJets ) {
    if(ditau_vis4V.Pt()==0) continue;
    float dR=ditau_vis4V.DeltaR(ditaujet->p4());
    m_histNameSvc->set_description("All");
    m_histSvc->BookFillHist("VisdRTruthDiTaus", 50, 0, 5, dR, m_weight);

    if(Props::BDTScore.get(ditaujet)>0.75) m_histSvc->BookFillHist("VisdRTruthDiTaus_DR_BDTID_075", 50, 0, 5, dR, m_weight);
    if(Props::BDTScore.get(ditaujet)>0.7) m_histSvc->BookFillHist("VisdRTruthDiTaus_DR_BDTID_07", 50, 0, 5, dR, m_weight);
    if(Props::BDTScore.get(ditaujet)>0.65) m_histSvc->BookFillHist("VisdRTruthDiTaus_DR_BDTID_065", 50, 0, 5, dR, m_weight);
    if(Props::BDTScore.get(ditaujet)>0.6) m_histSvc->BookFillHist("VisdRTruthDiTaus_DR_BDTID_06", 50, 0, 5, dR, m_weight);

    if(dR < 1.0){
      m_histNameSvc->set_description("TruthMatch");
      m_histSvc -> BookFillHist("diTauJetbdtScore", 20, 0, 1, Props::BDTScore.get(ditaujet), m_weight);
      m_histSvc -> BookFillHist("diTauJetpt",   100, 0, 2000, ditaujet->p4().Pt()*1e-3, m_weight);
      m_histSvc -> BookFillHist("diTauJeteta",   50, -4, 4, ditaujet->p4().Eta(), m_weight);
      m_histSvc -> BookFillHist("diTauJetphi",   50, -4, 4, ditaujet->p4().Phi(), m_weight);
      //fillDiTauJetPlots(my_ditauJets);
    }
  }

  if(taus.size() >= 2 && bis.size() >=2){
    TLorentzVector diTauTruthVec = taus.at(0)->p4() +taus.at(1)->p4();
    TLorentzVector diJetTruthVec = bis.at(0)->p4() + bis.at(1)->p4();
    if(m_debug) std::cout << "Vis Tau Pt: " << Props::pt_vis.get(taus.at(0)) << std::endl;
    if(m_debug) std::cout << "Truth DiJet: " << " pt: " << diJetTruthVec.Pt() << " eta: " << diJetTruthVec.Eta() << " phi: " << diJetTruthVec.Eta() << std::endl;
  }

}

void AnalysisReader_hhbbtt::findNBTagTJ()
{

  int nJet=0;
  m_have1BTagTJ.clear();m_have2BTagTJ.clear();

  for (auto jet: m_my_fatJets) {

    if(Props::xbbResult_1tagNoMcut.get(jet)==BoostedXbbTagger::InvalidJet) m_have1BTagTJ.push_back(-1);
    else if(Props::xbbResult_1tagNoMcut.get(jet)==BoostedXbbTagger::MassPassBTagPassJSSPass) m_have1BTagTJ.push_back(1);
    else if(Props::xbbResult_1tagNoMcut.get(jet)==BoostedXbbTagger::MassPassBTagFailJSSPass) m_have1BTagTJ.push_back(0);
    else m_have1BTagTJ.push_back(0);

    if(Props::xbbResult_2tagNoMcut.get(jet)==BoostedXbbTagger::InvalidJet) m_have2BTagTJ.push_back(-1);
    else if(Props::xbbResult_2tagNoMcut.get(jet)==BoostedXbbTagger::MassPassBTagPassJSSPass) m_have2BTagTJ.push_back(1);
    else if(Props::xbbResult_2tagNoMcut.get(jet)==BoostedXbbTagger::MassPassBTagFailJSSPass) m_have2BTagTJ.push_back(0);
    else m_have2BTagTJ.push_back(0);

    /*std::cout << "Value: " << Props::xbbResult_2tagNoMcut.get(jet) << " " << Props::xbbResult_1tagNoMcut.get(jet) << std::endl;
    std::cout << "Value: " << m_have1BTagTJ.at(nJet) << " " << m_have2BTagTJ.at(nJet) << std::endl;

    if(Props::xbbResult_2tagNoMcut.get(jet)!=BoostedXbbTagger::MassPassBTagPassJSSPass) std::cout << "XbbResult 2: " << Props::xbbResult_2tagNoMcut.get(jet) << " Mass " << jet->p4().M()*1e-3 << std::endl;
    if(Props::xbbResult_1tagNoMcut.get(jet)!=BoostedXbbTagger::MassPassBTagPassJSSPass) std::cout << "XbbResult 1: " << Props::xbbResult_1tagNoMcut.get(jet) << " Mass " << jet->p4().M()*1e-3 << std::endl;

    if(Props::xbbResult_2tagNoMcut.get(jet)==BoostedXbbTagger::BTagPass || Props::xbbResult_1tagNoMcut.get(jet)==BoostedXbbTagger::BTagPass) std::cout << "BTag pass: " << std::endl;*/

    nJet++;
  }

}

void AnalysisReader_hhbbtt::fillFatJetPlots()
{

  auto nJet=0;
  m_TaggedFJets.clear();
  m_NonTaggedFJets.clear();

  std::vector<const xAOD::Jet*> trkJets;
  for (auto jet: m_my_fatJets) {
    trkJets.clear();
    m_histSvc -> BookFillHist("fatJet_pt",   100, 0, 2000, jet->p4().Pt()*1e-3, m_weight);
    m_histSvc -> BookFillHist("fatJet_eta",   50, -3, 3, jet->p4().Eta(), m_weight);
    m_histSvc -> BookFillHist("fatJet_phi",   50, -4, 4, jet->p4().Phi(), m_weight);
    m_histSvc -> BookFillHist("fatJet_m",   100, 0, 2000, jet->p4().M()*1e-3, m_weight);

    m_histSvc -> BookFillHist("XbbResult_1tagNoMcut", 11, -1, 10, Props::xbbResult_1tagNoMcut.get(jet), m_weight);
    m_histSvc -> BookFillHist("XbbResult_1tag68Mcut", 11, -1, 10, Props::xbbResult_1tag68Mcut.get(jet), m_weight);
    m_histSvc -> BookFillHist("XbbResult_1tag90Mcut", 11, -1, 10, Props::xbbResult_1tag90Mcut.get(jet), m_weight);
    m_histSvc -> BookFillHist("XbbResult_2tagNoMcut", 11, -1, 10, Props::xbbResult_2tagNoMcut.get(jet), m_weight);
    m_histSvc -> BookFillHist("XbbResult_2tag68Mcut", 11, -1, 10, Props::xbbResult_2tag68Mcut.get(jet), m_weight);
    m_histSvc -> BookFillHist("XbbResult_2tag90Mcut", 11, -1, 10, Props::xbbResult_2tag90Mcut.get(jet), m_weight);

    if(m_have1BTagTJ.at(nJet) == 1 || m_have2BTagTJ.at(nJet) == 1) m_TaggedFJets.push_back(jet);
    else if(m_have1BTagTJ.at(nJet) < 0 && m_have2BTagTJ.at(nJet) < 0) m_NonTaggedFJets.push_back(jet);

    if (jet->getAssociatedObjects<xAOD::Jet>("GhostAntiKt2TrackJet", trkJets)) {
      int ntrkjets = 0;
      for (const xAOD::Jet *trkJ : trkJets) {
        if (!trkJ) continue; // if the trackjet is not valid then skip it
        if (!(trkJ->pt() / 1000. > 10. && fabs(trkJ->eta()) < 2.5)) continue;
        ntrkjets++;
        m_histSvc->BookFillHist("trkJets_pt", 400, 0, 2000, trkJ->pt()/1e3, m_weight);
        m_histSvc->BookFillHist("trkJets_MV2c10", 200, -1, 1, Props::MV2c10.get(trkJ), m_weight);
        m_histSvc->BookFillHist("trkJets_2d", 400, 0, 2000, 200, -1, 1, trkJ->pt()/1e3, Props::MV2c10.get(trkJ), 1.);
      }
      m_histSvc->BookFillHist("trkJets_n", 10, 0, 10, ntrkjets, m_weight);
    }

    nJet++;

  }

}

void AnalysisReader_hhbbtt::fillDiTauJetPlots(std::vector<const xAOD::DiTauJet*> my_ditauJets)
{

  TLorentzVector lead_subjet, sublead_subjet, ditau_subjets, ditau_subjets_met;

  for (auto ditau: my_ditauJets){
    if(!ditau || ditau==NULL) return;
    if(!Props::BDTScore.exists(ditau)) continue;
    m_histSvc -> BookFillHist("diTauJet_bdtScore", 20, 0, 1, Props::BDTScore.get(ditau), m_weight);
    m_histSvc -> BookFillHist("subjet_lead_pt", 100, 0, 2000, Props::subjet_lead_pt.get(ditau)*1e-3,m_weight);
    m_histSvc -> BookFillHist("subjet_lead_eta", 50, -4, 4, Props::subjet_lead_eta.get(ditau), m_weight);
    m_histSvc -> BookFillHist("subjet_lead_phi", 50, -4, 4, Props::subjet_lead_phi.get(ditau),m_weight);
    m_histSvc -> BookFillHist("subjet_lead_e", 100, 0, 2000, Props::subjet_lead_e.get(ditau)*1e-3,m_weight);
    m_histSvc -> BookFillHist("subjet_subl_pt", 100, 0, 2000, Props::subjet_subl_pt.get(ditau)*1e-3,m_weight);
    m_histSvc -> BookFillHist("subjet_subl_eta", 50, -4, 4, Props::subjet_subl_eta.get(ditau), m_weight);
    m_histSvc -> BookFillHist("subjet_subl_phi", 50, -4, 4, Props::subjet_subl_phi.get(ditau),m_weight);
    m_histSvc -> BookFillHist("subjet_subl_e", 100, 0, 2000, Props::subjet_subl_e.get(ditau)*1e-3,m_weight);
    m_histSvc -> BookFillHist("n_subjets", 10, 0, 10, Props::n_subjets.get(ditau),m_weight);
    m_histSvc -> BookFillHist("diTauJet_pt",   100, 0, 2000, ditau->p4().Pt()*1e-3, m_weight);
    m_histSvc -> BookFillHist("diTauJet_eta",   50, -4, 4, ditau->p4().Eta(), m_weight);
    m_histSvc -> BookFillHist("diTauJet_phi",   50, -4, 4, ditau->p4().Phi(), m_weight);
    m_histSvc -> BookFillHist("diTauJet_m",   100, 0, 2000, ditau->p4().M()*1e-3, m_weight);
    if(Props::n_subjets.get(ditau)>1){
      lead_subjet.SetPtEtaPhiE(Props::subjet_lead_pt.get(ditau),Props::subjet_lead_eta.get(ditau),Props::subjet_lead_phi.get(ditau),Props::subjet_lead_e.get(ditau));
      sublead_subjet.SetPtEtaPhiE(Props::subjet_subl_pt.get(ditau),Props::subjet_subl_eta.get(ditau),Props::subjet_subl_phi.get(ditau),Props::subjet_subl_e.get(ditau));
      ditau_subjets = lead_subjet+sublead_subjet;
      m_histSvc -> BookFillHist("diTauJet_subJets_M",100,0,2000,ditau_subjets.M()*1e-3,m_weight);
      ditau_subjets_met = lead_subjet+sublead_subjet+m_METVec;
      m_histSvc -> BookFillHist("diTauJet_subJetsMET_M",100,0,2000,ditau_subjets_met.M()*1e-3,m_weight);
    }
  }

}

std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> AnalysisReader_hhbbtt::fillDRMap(std::vector<const xAOD::DiTauJet*> my_ditauJets)
{

  std::pair <const xAOD::Jet*, const xAOD::DiTauJet*> par;
  std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> dRMap;
  dRMap.clear();
  //std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> m_dRMap;
  for (auto ditau: my_ditauJets) {
    if(!ditau || ditau==NULL) continue;
    for (auto fatjet: m_my_fatJets){
      if(!fatjet || fatjet==NULL) continue;
       par.first = fatjet; par.second = ditau;
       float dR = ditau->p4().DeltaR(fatjet->p4());
       m_histSvc->BookFillHist("deltaR", 100, 0, 5, dR, m_weight);
       dRMap[par] = dR;
    }
  }

  return dRMap;
}

void AnalysisReader_hhbbtt::getSelDT(std::vector<const xAOD::DiTauJet*> my_ditauJets)
{

  std::vector<float> BDTScore;
  int Ns = 0;
  for (auto ditau: my_ditauJets){
    if(Props::BDTScore.exists(ditau)) BDTScore.push_back(Props::BDTScore.get(ditau));
    Ns++;
  }

  m_histSvc->BookFillHist("Nditaus", 5, 0, 5, Ns, m_weight);

  if(BDTScore.size()==0){
    m_SelDT=NULL;
    return;
  }

  float ditaujet_bdtMax = *max_element(BDTScore.begin(), BDTScore.end());
  if(m_debug) std::cout << "bdt_max: " << ditaujet_bdtMax << std::endl;

  int index=-1, j=0;
  for (auto ditau: my_ditauJets){
    if(!Props::BDTScore.exists(ditau)) continue;
    if(ditaujet_bdtMax==Props::BDTScore.get(ditau)) index=j; // Select the diTAu object if it has the Max BDT score!!
    j++;
  }

  // Select the ditaujet with the higher BDT score
  if(index>-1) m_SelDT = my_ditauJets.at(index);
  else m_SelDT=NULL;

  if(m_debug) std::cout << "index: " << index << std::endl;
}

void AnalysisReader_hhbbtt::getSelFJ(std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> dRMap)
{
  std::pair <const xAOD::Jet*, const xAOD::DiTauJet*> par;
  par.second = m_SelDT;
  m_SelFJets.clear();
  m_MatchFJets.clear();
  int Ns=0;
  for (auto fatjet: m_my_fatJets){
     par.first = fatjet;
     float dR=dRMap[par];
     if(m_debug) std::cout << "dR(FJ,DTJ): " <<dR<<std::endl; 
     if(dR>0.5) m_SelFJets.push_back(fatjet);
     else m_MatchFJets.push_back(fatjet);
     Ns++;
  }

  m_histSvc->BookFillHist("Nfjets", 5, 0, 5, Ns, m_weight);

  if(m_SelFJets.size()==1) m_SelFJ=m_SelFJets.at(0);
  else if(m_SelFJets.size()>1){
    if(m_debug) std::cout << "more than one FJ found taking the leading one for now..." << std::endl;
    m_SelFJ=m_SelFJets.at(0);
  }else if(m_SelFJets.size()==0){
    m_SelFJ=NULL;
    return;
  }

  if(m_TaggedFJets.size()==0 && m_NonTaggedFJets.size()==0){
    m_nBTagTJ = 0;
    return;
  }else if(m_TaggedFJets.size()==0 && m_NonTaggedFJets.size()!=0){
    m_nBTagTJ = -1;
    return;
  }

  if(m_debug){
  std::cout << "@ tagged fjet vector size " <<m_TaggedFJets.size()<< std::endl; 
  std::cout << "@ tagged fjet vector size " <<m_have1BTagTJ.size()<< " " << m_have1BTagTJ.at(0) << std::endl; 
  std::cout << "@ tagged fjet vector size " <<m_have2BTagTJ.size()<< " " << m_have2BTagTJ.at(0) << std::endl;
  }

  int ijet=0, index=-1;
  for(auto fjet: m_my_fatJets){
    //std::cout << m_have1BTagTJ.at(ijet) << " " << m_have2BTagTJ.at(ijet) << std::endl;
    if((m_have1BTagTJ.at(ijet) == 1 || m_have2BTagTJ.at(ijet) == 1) && fjet == m_SelFJ) index = ijet;
    ijet++;
  }

  if(index == -1){
     m_nBTagTJ = 0;
     return;
  }else if(index>=0){
     // This is the number of B-tagged track jets inside a fat jet
     if(m_have2BTagTJ.at(index)==1) m_nBTagTJ = 2;
     else if(m_have1BTagTJ.at(index)==1) m_nBTagTJ = 1;
  }


}

void AnalysisReader_hhbbtt::fillFullSystemVars(std::map< std::pair <const xAOD::Jet*, const xAOD::DiTauJet*>, float> dRMap)
{

  // DPhi between selected DiTauJet and MET
  m_dPhi = m_SelDT->p4().DeltaPhi(m_METVec);
  m_dPhiFJ = m_SelFJ->p4().DeltaPhi(m_METVec);
  
  m_dEta = fabs(m_SelDT->p4().Eta()-m_SelFJ->p4().Eta());

  m_dR   = m_SelDT->p4().DeltaR(m_SelFJ->p4());

  TLorentzVector corrDTau, metVec;
  metVec.SetPtEtaPhiM(m_my_met->met(),m_SelDT->p4().Eta(),m_my_met->phi(),0);
  m_corrDTau = m_SelDT->p4()+metVec;
  m_corrDTau2 = m_SelDT->p4()+m_METVec;

  m_fullMass = (m_SelDT->p4()+m_SelFJ->p4()).M();

  /*std::pair <const xAOD::Jet*, const xAOD::DiTauJet*> par;
  par.second = m_SelDT;
  std::cout << "Point 1 " << m_my_fatJets.size() << std::endl;
  for(auto fatjet: m_my_fatJets){
     par.first = fatjet;
     float dR=dRMap[par];
     std::cout << "dR: " << dR << std::endl;
     if(dR<0.5) m_SelFJ2 = fatjet;
  }
  std::cout << "Point 2" << std::endl;*/

  if(m_debug) std::cout << "Match FJets vector size: " << m_MatchFJets.size() << std::endl;
  if(!m_MatchFJets.empty()){
     m_SelFJ2 = m_MatchFJets.at(0);

    const xAOD::Jet* LeadFJ;
    const xAOD::Jet* SubLFJ;

    if(m_SelFJ->pt()>m_SelFJ2->pt()){
      LeadFJ = m_SelFJ;
      SubLFJ = m_SelFJ2;
    }else{
      LeadFJ = m_SelFJ2;
      SubLFJ = m_SelFJ;
    }

    float mLeadFJ = LeadFJ->m()*1e-3;
    float mSubLFJ = SubLFJ->m()*1e-3;

    float C = (mLeadFJ-120)/20;
    float D = (mSubLFJ-100)/25;

    float E = (mLeadFJ-120)/30;
    float F = (mSubLFJ-100)/35;

    float G = (mLeadFJ-120)/40;
    float H = (mSubLFJ-100)/45;

    m_Xhh2 = sqrt(C*C+D*D);
    m_Xhh3 = sqrt(E*E+F*F);
    m_Xhh4 = sqrt(G*G+H*H);

    m_MCut2 = sqrt(pow(mLeadFJ-124,2)+pow(mSubLFJ-115,2));
  }


  float mFJ = m_SelFJ->m()*1e-3;
  float mDT = m_SelDT->m()*1e-3;

  float A = (mFJ-124)/(0.1*mFJ);
  float B = (mDT-80)/(0.5*mDT);

  m_Xhh  = sqrt(A*A+B*B);

  m_MCut = sqrt(pow(mFJ-124,2)+pow(mDT-80,2));

  m_MTW_Max = sqrt( mTsqr(m_SelDT->p4(),m_METVec) );
  m_MTW_Clos = sqrt( mTsqr(m_SelFJ->p4(),m_METVec) );

  m_ElCutX = sqrt(pow(m_MTW_Max-255,2)+pow(m_MTW_Clos,2));
  m_ElCutY = sqrt(pow(m_MTW_Max,2)+pow(m_MTW_Clos-255,2));

  if(m_debug) std::cout << "Sel DiTau M: " << m_SelDT->p4().M() << " Corr DiTauM: " << m_corrDTau.M() << std::endl;

}

void AnalysisReader_hhbbtt::fillFullSystemPlots()
{

  m_histSvc->BookFillHist("diTauMassvsfjJetMass", 100, 0, 500, 100, 0, 500, m_SelDT->p4().M()*1e-3, m_SelFJ->p4().M()*1e-3, m_weight);
  m_histSvc->BookFillHist("CorrDiTauMassvsfjJetMass", 100, 0, 500, 100, 0, 500, m_corrDTau.M()*1e-3, m_SelFJ->p4().M()*1e-3, m_weight);
  m_histSvc->BookFillHist("CorrDiTauMassvsfjJetMass", 100, 0, 500, 100, 0, 500, m_corrDTau.M()*1e-3, m_SelFJ->p4().M()*1e-3, m_weight);
  m_histSvc->BookFillHist("dPhidTauMET", 100, -5, 5, m_dPhi,m_weight);
  m_histSvc->BookFillHist("dPhifJetMET", 100, -5, 5, m_dPhiFJ,m_weight);
  m_histSvc->BookFillHist("dPhidTauMETdPhifJetMET", 100, -5, 5, 100, -5, 5, m_dPhi, m_dPhiFJ, m_weight);
  m_histSvc->BookFillHist("dEta",100,0,5,m_dEta,m_weight);
  m_histSvc->BookFillHist("MET", 100, 0, 1500, m_METVec.Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("nBTagTJ", 4, 0, 4, m_nBTagTJ, m_weight);
  m_histSvc->BookFillHist("FullXMass",100,200,8000,m_fullMass*1e-3,m_weight);
  m_histSvc->BookFillHist("deltaR", 100, 0, 5, m_dR, m_weight);
  m_histSvc->BookFillHist("Xhh",40,0,20,m_Xhh,m_weight);
  m_histSvc->BookFillHist("MCut",100,0,200,m_MCut,m_weight);
  m_histSvc->BookFillHist("MTW_Max",100,0,2000,m_MTW_Max,m_weight);
  m_histSvc->BookFillHist("MTW_Clos",100,0,2000,m_MTW_Clos,m_weight);
  m_histSvc->BookFillHist("MTW2D", 400, 0, 2000, 400, 0, 2000, m_MTW_Clos, m_MTW_Max, m_weight);
  m_histSvc->BookFillHist("SelfatJet_pt",   100, 500, 3000, m_SelFJ->p4().Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("SelfatJet_eta",   50, -3, 3, m_SelFJ->p4().Eta(), m_weight);
  m_histSvc->BookFillHist("SelfatJet_phi",   50, -4, 4, m_SelFJ->p4().Phi(), m_weight);
  m_histSvc->BookFillHist("SelfatJet_m",   100, 0, 2000, m_SelFJ->p4().M()*1e-3, m_weight);
  m_histSvc->BookFillHist("SelditauJet_pt",   100, 400, 3000, m_SelDT->p4().Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("SelditauJet_eta",   50, -3, 3, m_SelDT->p4().Eta(), m_weight);
  m_histSvc->BookFillHist("SelditauJet_phi",   50, -4, 4, m_SelDT->p4().Phi(), m_weight);
  m_histSvc->BookFillHist("SelditauJet_m",   50, 0, 2000, m_SelDT->p4().M()*1e-3, m_weight);
  m_histSvc->BookFillHist("SeldiTauJet_bdtScore", 20, 0, 1, Props::BDTScore.get(m_SelDT), m_weight);
  m_histSvc->BookFillHist("Selsubjet_lead_pt", 100, 0, 2000, Props::subjet_lead_pt.get(m_SelDT)*1e-3,m_weight);
  m_histSvc->BookFillHist("Selsubjet_lead_eta", 50, -4, 4, Props::subjet_lead_eta.get(m_SelDT), m_weight);
  m_histSvc->BookFillHist("Selsubjet_lead_phi", 50, -4, 4, Props::subjet_lead_phi.get(m_SelDT),m_weight);
  m_histSvc->BookFillHist("Selsubjet_lead_e", 100, 0, 2000, Props::subjet_lead_e.get(m_SelDT)*1e-3,m_weight);
  m_histSvc->BookFillHist("Selsubjet_subl_pt", 100, 0, 2000, Props::subjet_subl_pt.get(m_SelDT)*1e-3,m_weight);
  m_histSvc->BookFillHist("Selsubjet_subl_eta", 50, -4, 4, Props::subjet_subl_eta.get(m_SelDT), m_weight);
  m_histSvc->BookFillHist("Selsubjet_subl_phi", 50, -4, 4, Props::subjet_subl_phi.get(m_SelDT),m_weight);
  m_histSvc->BookFillHist("Selsubjet_subl_e", 100, 0, 2000, Props::subjet_subl_e.get(m_SelDT)*1e-3,m_weight);
  m_histSvc->BookFillHist("Seln_subjets", 10, 0, 10, Props::n_subjets.get(m_SelDT),m_weight);

  TLorentzVector lead_subjet, sublead_subjet, ditau_subjets, ditau_subjets_met;
  if(Props::n_subjets.get(m_SelDT) >1){
     lead_subjet.SetPtEtaPhiE(Props::subjet_lead_pt.get(m_SelDT),Props::subjet_lead_eta.get(m_SelDT),Props::subjet_lead_phi.get(m_SelDT),Props::subjet_lead_e.get(m_SelDT));
     sublead_subjet.SetPtEtaPhiE(Props::subjet_subl_pt.get(m_SelDT),Props::subjet_subl_eta.get(m_SelDT),Props::subjet_subl_phi.get(m_SelDT),Props::subjet_subl_e.get(m_SelDT));
     ditau_subjets = lead_subjet+sublead_subjet;
     m_histSvc -> BookFillHist("SeldiTauJet_subJets_M",100,0,2000,ditau_subjets.M()*1e-3,m_weight);
     ditau_subjets_met = lead_subjet+sublead_subjet+m_METVec;
     m_histSvc -> BookFillHist("SeldiTauJet_subJetsMET_M",100,0,2000,ditau_subjets_met.M()*1e-3,m_weight);
     m_histSvc -> BookFillHist("DiTauSubJMassvsfjJetMass", 100, 0, 500, 100, 0, 500, ditau_subjets.M()*1e-3, m_SelFJ->p4().M()*1e-3, m_weight);
  }

  std::vector<const xAOD::Jet*> trkJets, seltrkJets;
  if (m_SelFJ->getAssociatedObjects<xAOD::Jet>("GhostAntiKt2TrackJet", trkJets)) {
    for (const xAOD::Jet *trkJ : trkJets) {
      if (!trkJ) continue; // if the trackjet is not valid then skip it
      if (!(trkJ->pt() / 1000. > 10. && fabs(trkJ->eta()) < 2.5)) continue;
      seltrkJets.push_back(trkJ);
    }
  }
  m_nTJ=seltrkJets.size();
  m_histSvc->BookFillHist("nTJ", 10, 0, 10, seltrkJets.size(), m_weight);
  for(auto jet:seltrkJets){
     m_histSvc->BookFillHist("Seltj_pt", 100, 0, 2000, jet->pt()*1e-3,m_weight);
     m_histSvc->BookFillHist("Seltj_eta", 50, -4, 4, jet->eta(), m_weight);
     m_histSvc->BookFillHist("Seltj_phi", 50, -4, 4, jet->phi(),m_weight);
     m_histSvc->BookFillHist("Seltj_e", 100, 0, 2000, jet->e()*1e-3,m_weight);
     m_histSvc->BookFillHist("Seltj_MV2c10", 100, -1, 1, Props::MV2c10.get(jet),m_weight);
  }

  if(!m_MatchFJets.empty()){
    const xAOD::Jet* LeadFJ;
    const xAOD::Jet* SubLFJ;

    if(m_SelFJ->pt()>m_SelFJ2->pt()){
      LeadFJ = m_SelFJ;
      SubLFJ = m_SelFJ2;
    }else{
      LeadFJ = m_SelFJ2;
      SubLFJ = m_SelFJ;
    }

    float mLeadFJ = LeadFJ->m()*1e-3;
    float mSubLFJ = SubLFJ->m()*1e-3;

    m_histSvc->BookFillHist("FatJetMasses2DXhh",100,0,250,100,0,250,mLeadFJ,mSubLFJ,m_Xhh2);
    m_histSvc->BookFillHist("FatJetMasses2DMCut",100,0,250,100,0,250,mLeadFJ,mSubLFJ,m_MCut2);
    m_histSvc->BookFillHist("FatJetMasses2D",100,0,250,100,0,250,mLeadFJ,mSubLFJ,m_weight);

    m_histSvc->BookFillHist("Xhh2",40,0,20,m_Xhh2,m_weight);
    m_histSvc->BookFillHist("MCut2",100,0,200,m_MCut2,m_weight);
  }

  m_histSvc->BookFillHist("FJDT2DMassXhh",100,0,250,100,0,250,m_SelDT->m()*1e-3,m_SelFJ->m()*1e-3,m_Xhh);
  m_histSvc->BookFillHist("FJDT2DMassCut",100,0,250,100,0,250,m_SelDT->m()*1e-3,m_SelFJ->m()*1e-3,m_MCut);
  m_histSvc->BookFillHist("FJDT2D",100,0,250,100,0,250,m_SelDT->m()*1e-3,m_SelFJ->m()*1e-3,m_weight);

  //float 
  //m_histSvc->BookFillHist("FatJetMasses2D2",100,0,250,100,0,250,m_)

}

bool AnalysisReader_hhbbtt::applyPreSelection()
{

  /*if(Props::n_subjets.get(m_SelDT) == 0) return false;  
  //m_histNameSvc->set_description("PreSel_PostDiTauJetPtCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 1 " << m_SelDT->p4().Pt() << std::endl;

  if(m_SelDT->p4().Pt()*1e-3<450) return false;
  //m_histNameSvc->set_description("PreSel_PostDiTauJetPtCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 2 " << m_SelDT->p4().Pt() << std::endl;

  if(m_SelFJ->p4().Pt()*1e-3<600) return false;

  //m_histNameSvc->set_description("PreSel_PostFatJetPtCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 3 " << m_SelFJ->p4().Pt() << std::endl;*/

  if(m_dR<1.0) return false;

  m_histNameSvc->set_description("PreSel_PostDRCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 4 " << m_dR << std::endl;

  if(m_dEta>1.7) return false;

  m_histNameSvc->set_description("PreSel_PostDEtaCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 4 " << m_dR << std::endl;

  return true;
}

bool AnalysisReader_hhbbtt::applySRSelection()
{

  if(fabs(m_dPhi)>TMath::Pi()/2) return false;
  //m_histNameSvc->set_description("SR_PostDPhiCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applySRSelection 1 " << m_dPhi << std::endl;

  if(Props::BDTScore.get(m_SelDT)<0.6) return false;
  //m_histNameSvc->set_description("SR_PostDiTauBDTScoreCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applySRSelection 2 " << Props::BDTScore.get(m_SelDT) << std::endl;

  if(m_nBTagTJ == 0) return false;
  m_histNameSvc->set_description("SR_PostnBTagTJCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applySRSelection 3 " << m_nBTagTJ << std::endl;

  if(Props::n_subjets.get(m_SelDT) > 4 ) return false;

  m_histNameSvc->set_description("SR_PostNSubJetsCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 5 " << Props::n_subjets.get(m_SelDT) << std::endl;

  if(m_Xhh2 > 1 ) return false;

  m_histNameSvc->set_description("SR_PostEllipseCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 6 " << m_Xhh2 << std::endl;


  return true;

}

bool AnalysisReader_hhbbtt::applyQCDSelection()
{

  //if(fabs(m_dPhi)>TMath::Pi()/2) return false;
  //m_histNameSvc->set_description("QCDCR_PostDPhiCut");
  //fillFullSystemPlots();

  //if(m_debug) std::cout << "@ applyQCDCRSelection 1 " << m_dPhi << std::endl;

  if(Props::BDTScore.get(m_SelDT)>=0.6) return false;
  //m_histNameSvc->set_description("QCDCR_PostDiTauBDTScoreCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyQCDCRSelection 2 " << Props::BDTScore.get(m_SelDT) << std::endl;

  //if(m_nBTagTJ > 0) return false;
  //m_histNameSvc->set_description("QCDCR_PostnBTagTJCut");
  //fillFullSystemPlots();

  //if(m_debug) std::cout << "@ applyQCDCRSelection 3 " << m_nBTagTJ << std::endl;

  /*if(m_ElCutX <= 85){
    m_histNameSvc->set_description("QCDCR_PostEllipseXCut");
    fillFullSystemPlots();
    if(m_debug) std::cout << "@ applyQCDCRSelection 4 " << m_ElCutX << std::endl;
  }


  if(m_ElCutY <= 85){
    m_histNameSvc->set_description("QCDCR_PostEllipseYCut");
    fillFullSystemPlots();
    if(m_debug) std::cout << "@ applyQCDCRSelection 5 " << m_ElCutY << std::endl;
  }*/

  if(m_Xhh2 <= 1 || m_Xhh3 > 1 ) return false;

  m_histNameSvc->set_description("QCDCR_PostEllipseCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyPreSelection 6 " << m_Xhh2 << std::endl;

  return true;
}

bool AnalysisReader_hhbbtt::applySideBandSelection()
{

  //if(fabs(m_dPhi)>TMath::Pi()/2) return false;
  //m_histNameSvc->set_description("QCDCR_PostDPhiCut");
  //fillFullSystemPlots();

  //if(m_debug) std::cout << "@ applyQCDCRSelection 1 " << m_dPhi << std::endl;

  //m_histNameSvc->set_description("QCDCR_PostDiTauBDTScoreCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applySideBandSelection 2 " << Props::BDTScore.get(m_SelDT) << std::endl;

  //if(m_nBTagTJ > 0) return false;
  //m_histNameSvc->set_description("QCDCR_PostnBTagTJCut");
  //fillFullSystemPlots();

  //if(m_debug) std::cout << "@ applyQCDCRSelection 3 " << m_nBTagTJ << std::endl;

  /*if(m_ElCutX <= 85){
    m_histNameSvc->set_description("QCDCR_PostEllipseXCut");
    fillFullSystemPlots();
    if(m_debug) std::cout << "@ applyQCDCRSelection 4 " << m_ElCutX << std::endl;
  }


  if(m_ElCutY <= 85){
    m_histNameSvc->set_description("QCDCR_PostEllipseYCut");
    fillFullSystemPlots();
    if(m_debug) std::cout << "@ applyQCDCRSelection 5 " << m_ElCutY << std::endl;
  }*/

  if(m_Xhh3 <= 1 || m_Xhh4 > 1 ) return false;

  if(m_applyBoostedFF){

    m_histNameSvc->set_sample("MultiJet");
    m_histNameSvc->set_description("SR_PostEllipseCut");

    float temp_weight=m_weight;

    //m_histNameSvc->set_nTag(1);
    if(m_isMC) m_weight *= -1;
    else { 
     if(m_nBTagTJ==1) m_weight *= m_weightFactor1;
     else if(m_nBTagTJ==2) m_weight *= m_weightFactor2;
     else if(m_nBTagTJ==0) m_weight = temp_weight;
    }
    fillFullSystemPlots();

    m_weight = temp_weight;

    /*m_histNameSvc->set_nTag(2);
    if(m_isMC) m_weight *= -1;
    else m_weight *= m_weightFactor2;
    fillFullSystemPlots();

    m_weight = temp_weight;*/

  }

  if(Props::BDTScore.get(m_SelDT)<0.6) {
    m_histNameSvc->set_description("SBAnti2Tau_PostEllipseCut");
    if(!m_applyBoostedFF) fillFullSystemPlots();
  }else if(Props::BDTScore.get(m_SelDT)>=0.6){
    m_histNameSvc->set_description("SBMedium2Tau_PostEllipseCut");
    if(!m_applyBoostedFF) fillFullSystemPlots();
  }

  if(m_debug) std::cout << "@ applySideBandSelection 6 " << m_Xhh2 << std::endl;

  return true;
}

bool AnalysisReader_hhbbtt::applyTTBarSelection()
{

  //if(fabs(m_dPhi)<TMath::Pi()/2) return false;
  //m_histNameSvc->set_description("TTBarCR_PostDPhiCut");
  //fillFullSystemPlots();

  //if(m_debug) std::cout << "@ applyTTBarCRSelection 1 " << m_dPhi << std::endl;

  if(Props::BDTScore.get(m_SelDT)>=0.6) return false;
  //m_histNameSvc->set_description("QCDCR_PostDiTauBDTScoreCut");
  //fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyTTBarCRSelection 1 " << Props::BDTScore.get(m_SelDT) << std::endl;

  if(m_nBTagTJ == 0) return false;

  m_histNameSvc->set_description("TTBarCR_PostnBTagTJCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyTTBarCRSelection 2 " << m_nBTagTJ << std::endl;

  if(m_SelFJ->p4().M()*1e-3 < 140.0 || m_SelFJ->p4().M()*1e-3 > 200.0) return false;
  if(m_SelDT->p4().M()*1e-3 < 120.0 || m_SelDT->p4().M()*1e-3 > 220.0) return false;

  m_histNameSvc->set_description("TTBarCR_Post2DSquareMCut");
  fillFullSystemPlots();

  if(m_debug) std::cout << "@ applyTTBarCRSelection 3 " << m_SelFJ->p4().M()*1e-3 << " " << m_SelDT->p4().M()*1e-3 << std::endl;

  return true;
}

void AnalysisReader_hhbbtt::fillSRPlots()
{
}

void AnalysisReader_hhbbtt::fillQCDCRPlots()
{
}

void AnalysisReader_hhbbtt::fillSideBandPlots()
{
}

void AnalysisReader_hhbbtt::fillTTBarCRPlots()
{
}

void AnalysisReader_hhbbtt::fill_xHistos(std::vector<const xAOD::Jet*> jets, std::vector<const xAOD::TauJet*> taus)
{

  if(m_debug) std::cout << "On fill_xHistos function" << std::endl;

  if(m_region=="MM_OS" && m_nBJets == 2 && m_debug) std::cout << "Weight: " << m_weight << std::endl;

  std::vector<TLorentzVector> leadjets;
  leadjets.push_back(jets.at(0)->p4());leadjets.push_back(jets.at(1)->p4());

  m_bb=jets.at(0)->p4()+jets.at(1)->p4();
  TLorentzVector bb2 = m_bb;
  TLorentzVector tautau=taus.at(0)->p4()+taus.at(1)->p4();

  TLorentzVector X=m_bb+tautau;
  m_Xmmc=m_bb+m_MMCVec;

  //if(fabs(m_MMCVec.M()*1e-3-125)<25 && !m_isMC) std::cout << "m_MMCVec.M() : " << m_MMCVec.M() << " Xmmc.M() : " << Xmmc.M() << std::endl;

  m_deltaRJJ=jets.at(0)->p4().DeltaR(jets.at(1)->p4());
  m_deltaPhiJJ=jets.at(0)->p4().DeltaPhi(jets.at(1)->p4());
  m_deltaEtaJJ=jets.at(0)->eta()-jets.at(1)->eta();

  m_deltaRTT=taus.at(0)->p4().DeltaR(taus.at(1)->p4());
  m_deltaPhiTT=taus.at(0)->p4().DeltaPhi(taus.at(1)->p4());
  m_deltaEtaTT=taus.at(0)->eta()-taus.at(1)->eta();

  m_MTW_Max = sqrt( mTsqr(taus.at(0)->p4(),m_METVec) );
  m_mTtot = sqrt( mTsqr(taus.at(0)->p4(),taus.at(1)->p4()) + mTsqr(taus.at(0)->p4(),m_METVec) + mTsqr(taus.at(1)->p4(),m_METVec) );

  TLorentzVector met_tau = m_METVec + taus.at(0)->p4() + taus.at(1)->p4();

  m_mt2 = getMT2(leadjets,met_tau);

  TLorentzVector tauCloseMet=taus.at(0)->p4();
  if(fabs(taus.at(1)->phi()-m_METVec.Phi())< fabs(taus.at(0)->phi()-m_METVec.Phi())){
    tauCloseMet=taus.at(1)->p4();
  }

  if(m_debug) std::cout << "mT2: " << m_mt2 << " mmc: " << m_MMCVec.M() << " mmc pT: " << Props::mmc_mlnu3p_4vect_pt.get(m_mmc) << std::endl;

  m_MTW_Clos = sqrt( mTsqr(tauCloseMet,m_METVec) );

  float A = sin(m_METVec.Phi()-taus.at(0)->phi())/sin(taus.at(1)->phi()-taus.at(0)->phi());
  float B = sin(taus.at(1)->phi()-m_METVec.Phi())/sin(taus.at(1)->phi()-taus.at(0)->phi());
  m_METCentrality = (A+B)/sqrt(A*A+B*B);

  float pTBins[13]={0.,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0,180.0,220.0,300.0,500.0};
  float METBins[10]={0.,12.,24.,36.,48.,60.,80.,120.,200.,400.};

  m_dRDiJetDiTau=m_bb.DeltaR(m_MMCVec);
  m_dPhiDiJetDiTau=m_bb.DeltaPhi(m_MMCVec);

  //scaling the dijet and ditau mass to the higgs mass
  float scalebb = 125.0/m_bb.M();
  m_bbScaled=DiJetVec2*scalebb;
  float scalemmc = 125.0/m_MMCVec.M();
  m_MMCVecScaled= m_MMCVec*scalemmc;

  m_XmmcScaled=m_bbScaled+m_MMCVecScaled;

  m_histSvc->BookFillHist("pTBB", 12, pTBins, m_bb.Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("pTTT", 12, pTBins, tautau.Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("mBB", 23, 20, 250, m_bb.M()*1e-3, m_weight);
  m_histSvc->BookFillHist("mTauTau", 23, 20, 250, tautau.M()*1e-3, m_weight);
  m_histSvc->BookFillHist("mmc",23,20,250,m_MMCVec.M(),m_weight);
  m_histSvc->BookFillHist("mmcpT", 12, pTBins, m_MMCVec.Pt(), m_weight);
  m_histSvc->BookFillHist("mX", 100, 0, 1200, X.M()*1e-3, m_weight);
  m_histSvc->BookFillHist("mXmmc", 100, 0, 1200, m_Xmmc.M()*1e-3, m_weight);
  m_histSvc->BookFillHist("MET", 9,METBins, m_METVec.Pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("MTW_Max", 30, 0, 600, m_MTW_Max, m_weight);
  m_histSvc->BookFillHist("MTW_Clos", 20, 0, 300, m_MTW_Clos, m_weight);
  m_histSvc->BookFillHist("m_mTtot",30,0,600,m_mTtot, m_weight);
  m_histSvc->BookFillHist("m_mt2",21,80,500,m_mt2*1e-3, m_weight);
  m_histSvc->BookFillHist("METPhiCentrality", 25, -1.5, 1.5, m_METCentrality, m_weight);

  m_histSvc->BookFillHist("DeltaRJJ", 20, 0, 5, m_deltaRJJ, m_weight);
  m_histSvc->BookFillHist("DeltaRTauTau", 20, 0, 5, m_deltaRTT, m_weight);
  m_histSvc->BookFillHist("DeltaEtaTauTau", 20, -5, 5, m_deltaEtaTT, m_weight);
  m_histSvc->BookFillHist("DeltaEtaJJ", 20, -5, 5, m_deltaEtaJJ, m_weight);
  m_histSvc->BookFillHist("DeltaPhiJJ", 20, -5, 5, m_deltaPhiJJ, m_weight);
  m_histSvc->BookFillHist("DeltaPhiTauTau", 20, -5, 5, m_deltaPhiTT, m_weight);

  m_histSvc->BookFillHist("2DFinalTTBB", 100, 0, 500, 100, 0, 500, m_bb.Pt()/1e3, tautau.Pt()/1e3, m_weight);

  std::string hist_name = m_histNameSvc->getFullHistName("mmc");
  if(m_debug) std::cout << "Full Hist Name: " << hist_name << std::endl;

  if(m_writeMVATree){
    m_tree->diTauVisM = tautau.M()/1e3;
    m_tree->diTauVisPt = tautau.Pt()/1e3;
    m_tree->diTauVisEta = tautau.Eta();
    m_tree->diTauVisPhi = tautau.Phi();

    m_tree->diTauDR   = m_deltaRTT;
    m_tree->diTauDEta = m_deltaEtaTT;
    m_tree->diTauDPhi = m_deltaPhiTT;

   m_tree->diJetM = m_bb.M()/1e3;
   m_tree->diJetPt = m_bb.Pt()/1e3;
   m_tree->diJetEta = m_bb.Eta();
   m_tree->diJetPhi = m_bb.Phi();

   m_tree->diJetDR = m_deltaRJJ;
   m_tree->diJetDEta = m_deltaEtaJJ;
   m_tree->diJetDPhi = m_deltaPhiJJ;
   
   m_tree->diTauMMCM = m_MMCVec.M();
   m_tree->diTauMMCPt = m_MMCVec.Pt();
   m_tree->diTauMMCEta = m_MMCVec.Eta();
   m_tree->diTauMMCPhi = m_MMCVec.Phi();
     
   m_tree->diHiggsM  = m_Xmmc.M();
   m_tree->diHiggsPt = m_Xmmc.Pt();

   m_tree->diHiggsMScaled =m_XmmcScaled.M();

   m_tree->diJetdiTauDR = m_deltaRDiJetDiTau;
   m_tree->diJetdiTauDPhi = m_deltaPhiDiJetDiTau;

   m_tree->MTW_Max =m_MTW_Max;
   m_tree->MTW_Clos = m_MTW_Clos;
   m_tree->METCentrality = m_METCentrality;
   
   m_tree->MET = m_METVec.Pt()/1e3;

  }

  


}

void AnalysisReader_hhbbtt::fill_SRPlots()
{
       //Filling histograms after signal selection - cut-based
       if(m_FillSRSel260_300){
         m_histSvc->BookFillHist("Sel_M260300_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
         m_histSvc->BookFillHist("Sel_M260300_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
       }
       if(m_FillSRSel400){
         if(m_bb.Pt()>100){//DiJet pt cut for M400
            m_histSvc->BookFillHist("Sel_M400_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
            m_histSvc->BookFillHist("Sel_M400_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
         }
       }
       if(m_FillSRSel500_600_700){
         if(m_bb.Pt()>150){//DiJet pt cut for M500600700
            m_histSvc->BookFillHist("Sel_M500600700_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
            m_histSvc->BookFillHist("Sel_M500600700_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
         }
       }
       if(m_FillSRSel800_900_1000){
         if(m_bb.Pt()>200){//DiJet pt cut for M8009001000
            m_histSvc->BookFillHist("Sel_M8009001000_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
            m_histSvc->BookFillHist("Sel_M8009001000_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
         }
       }


   //Filling histograms after signal selection - BDT cut for MVA studies
   if(m_readMVA){
     if(m_tree->BDT>0){
       m_histSvc->BookFillHist("BDT0_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT0_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.1){
       m_histSvc->BookFillHist("BDT1_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT1_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.2){
       m_histSvc->BookFillHist("BDT2_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT2_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.3){
       m_histSvc->BookFillHist("BDT3_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT3_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.4){
       m_histSvc->BookFillHist("BDT4_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT4_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.5){
       m_histSvc->BookFillHist("BDT5_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT5_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.6){
       m_histSvc->BookFillHist("BDT6_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT6_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.7){
       m_histSvc->BookFillHist("BDT7_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT7_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.8){
       m_histSvc->BookFillHist("BDT8_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT8_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
     if(m_tree->BDT>0.9){
       m_histSvc->BookFillHist("BDT9_DiHiggsM", 300/2, 0, 3000, m_Xmmc.M(), m_weight);
       m_histSvc->BookFillHist("BDT9_DiHiggsMScaled", 300/2, 0, 3000, m_XmmcScaled.M(), m_weight);
     }
   }

}

void AnalysisReader_hhbbtt::fill_eventLVars()
{
  if(m_writeMVATree){
   m_tree->EventWeight = m_weight;
   m_tree->EventNumber = m_eventInfo->eventNumber();
   m_tree->sample = m_histNameSvc->getFullSample();
  }
}

bool AnalysisReader_hhbbtt::fill_cutFlow(std::vector<const xAOD::Jet*> jets, std::vector<const xAOD::TauJet*> taus)
{
  //NEW cut-flow-----------------------------------------------------------
  bool passCutFlow=false;

  m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 0, 1);//start
  m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,0, m_weight);

  if(m_RegA || m_RegB || m_RegC || m_RegD){
    m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 1, 1);
    m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,1, m_weight);
    
    if(m_PassTriggerSelection){//trigger DTT or STT
      m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 2, 1);
      m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,2, m_weight);
      if(m_debug) std::cout << "Tau 0 pT: " << taus.at(0)->pt()/1e3 << " Tau 1 pT: " << taus.at(1)->pt()/1e3 << std::endl;
      if(taus.at(0)->pt()/1e3>40 && taus.at(1)->pt()/1e3>20){//tau pt 40,20
        m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 3, 1);
        m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,3, m_weight);
     
    	if(jets.size()>1){
          m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 4, 1);//jet pt 20,20
          m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,4, m_weight);
          if(m_debug) std::cout << "MMC: " << m_MMCVec.M() << std::endl;
          if(m_MMCVec.M()>60){
            m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 5, 1);//mmc>60
            m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,5, m_weight);
            if(m_debug) std::cout << "PassTriggerJetSelection: " << m_PassTriggerJetSelection << std::endl;
            if(m_PassTriggerJetSelection){
              m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 6, 1);//jet pt cut 80 for DTT, 45 for STT
              m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,6, m_weight);
              if(m_debug) std::cout << "m_PassLeadingJetPtCut: " << m_PassLeadingJetPtCut << std::endl;
              if(m_PassLeadingJetPtCut){//jet pt cut 45
                m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 7, 1);
                m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,7, m_weight);
                if(m_debug) std::cout << "m_PassTauPtSelection: " << m_PassTauPtSelection << std::endl;
                //if(m_PassTauPtSelection){//tau pt selection for 3p - Note 240717 - Removed from preselection 
                //  m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 8, 1);
                //  m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,8, m_weight);
                  if(m_debug) std::cout << "m_nBJets: " << m_nBJets << std::endl;
                  if(m_nBJets==0){
                    if(m_debug) std::cout << "At 0 tag" << std::endl;
                    m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 9, 1);//0Tag
                    m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0.,20.,9, m_weight);
                  }else if(m_nBJets==1){
                    if(m_debug) std::cout << "At 1 tag" << std::endl;
                    m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 10, 1);//1Tag
                    m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0.,20.,10, m_weight);
                  }else if(m_nBJets==2){
                    if(m_debug) std::cout << "At 2 tag" << std::endl;
                    m_histSvc->BookFillHist("CutFlow_New", 20, 0.,20., 11, 1);//2Tag
                    m_histSvc->BookFillHist("CutFlow_New_Weighted", 20, 0., 20.,11, m_weight);
                  }

                  if(m_nBJets<3){
                    passCutFlow=true;
                  }
                //}
              }
            }
          }
	     }
      }
    }
  }

  if(m_debug) std::cout << "PassCutFlow: " << passCutFlow << std::endl;

  return passCutFlow;
}

void AnalysisReader_hhbbtt::fill_jetHistos(std::vector<const xAOD::Jet*> jets)
{
  //Leading jet histos
  m_histSvc->BookFillHist("Jet0Eta",50,-5,5, jets.at(0)->eta(), m_weight);
  m_histSvc->BookFillHist("Jet0Phi",50,-5,5, jets.at(0)->phi(), m_weight);
  m_histSvc->BookFillHist("Jet0Pt",100,0,1000, jets.at(0)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("Jet0isTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(0)), m_weight);
  m_histSvc->BookFillHist("Jet0tagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(0)), m_weight);

  //Sub-Leading jet histos
  m_histSvc->BookFillHist("Jet1Eta",50,-5,5, jets.at(1)->eta(), m_weight);
  m_histSvc->BookFillHist("Jet1Phi",50,-5,5, jets.at(1)->phi(), m_weight);
  m_histSvc->BookFillHist("Jet1Pt",100,0,1000, jets.at(1)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("Jet1isTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(1)), m_weight);
  m_histSvc->BookFillHist("Jet1tagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(1)), m_weight);

  //Inclusive jet histos
  for(unsigned int i=0;i<jets.size();i++){
     m_histSvc->BookFillHist("JetEta",50,-5,5, jets.at(i)->eta(), m_weight);
     m_histSvc->BookFillHist("JetPhi",50,-5,5, jets.at(i)->phi(), m_weight);
     m_histSvc->BookFillHist("JetPt",100,0,1000, jets.at(i)->pt()*1e-3, m_weight);
     m_histSvc->BookFillHist("JetisTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(i)), m_weight);
     m_histSvc->BookFillHist("JettagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(i)), m_weight);
  }

  m_histSvc->BookFillHist("NJets", 10, 0, 10, jets.size(), m_weight);
  m_histSvc->BookFillHist("NBJets", 10, 0, 10, m_nBJets, m_weight);

  if(m_writeMVATree){
    m_tree->Jet1Pt = jets.at(0)->pt()/1e3;
    m_tree->Jet1Eta = jets.at(0)->eta();
    m_tree->Jet1Phi = jets.at(0)->phi();
    m_tree->Jet1M = jets.at(0)->m()/1e3;

    m_tree->Jet2Pt = jets.at(1)->pt()/1e3;
    m_tree->Jet2Eta = jets.at(1)->eta();
    m_tree->Jet2Phi = jets.at(1)->phi();
    m_tree->Jet2M = jets.at(1)->m()/1e3;

    m_tree->NJets = jets.size();
    m_tree->NJetsbtagged = m_nBJets;
  }

}

void AnalysisReader_hhbbtt::fill_tauHistos(std::vector<const xAOD::TauJet*> taus)
{
  m_histNameSvc->set_nProng(-1);
  m_histSvc->BookFillHist("Tau0nTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(0)),m_weight);
  m_histSvc->BookFillHist("Tau1nTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(1)),m_weight);

  //Leading jet histos
  m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(0)));
  m_histSvc->BookFillHist("Tau0Eta",30,-3,3, taus.at(0)->eta(), m_weight);
  m_histSvc->BookFillHist("Tau0Pt",100,0,1000, taus.at(0)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("Tau0Phi",30,-5,5, taus.at(0)->phi(), m_weight);
  m_histSvc->BookFillHist("Tau0BDTScore",20,0.35,1,Props::BDTScore.get(taus.at(0)), m_weight);

  //Sub-Leading jet histos
  m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(1)));
  m_histSvc->BookFillHist("Tau1Eta",30,-3,3, taus.at(1)->eta(), m_weight);
  m_histSvc->BookFillHist("Tau1Pt",100,0,1000, taus.at(1)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("Tau1Phi",30,-5,5, taus.at(1)->phi(), m_weight);
  m_histSvc->BookFillHist("Tau1BDTScore",20,0.35,1,Props::BDTScore.get(taus.at(1)), m_weight);

  int nTaus=0,nSignalTaus=0,nAntiTaus=0;
  //Inclusive jet histos
  for(unsigned int i=0;i<taus.size();i++){
     m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(i)));
     m_histSvc->BookFillHist("TauEta",50,-5,5, taus.at(i)->eta(), m_weight);
     m_histSvc->BookFillHist("TauPhi",50,-5,5, taus.at(i)->phi(), m_weight);
     m_histSvc->BookFillHist("TauPt",100,0,1000, taus.at(i)->pt()*1e-3, m_weight);
     m_histSvc->BookFillHist("TauBDTScore",20,0.35,1, Props::BDTScore.get(taus.at(i)), m_weight);
     m_histNameSvc->set_nProng(-1);
     m_histSvc->BookFillHist("TaunTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(i)),m_weight);
     if(Props::isBDTMedium.get(taus.at(i))) nSignalTaus++;
     if(!Props::isBDTMedium.get(taus.at(i))) nAntiTaus++;
  }

  m_histSvc->BookFillHist("NTaus", 10, 0, 10, taus.size(), m_weight);
  m_histSvc->BookFillHist("NSignalTaus", 10, 0, 10, nSignalTaus, m_weight);
  m_histSvc->BookFillHist("NAntiTaus", 10, 0, 10, nAntiTaus, m_weight);

  if(m_writeMVATree){
    m_tree->Tau1Pt = taus.at(0)->pt()/1e3;
    m_tree->Tau1Eta = taus.at(0)->eta();
    m_tree->Tau1Phi = taus.at(0)->phi();
    m_tree->Tau1M = taus.at(0)->m()/1e3;

    m_tree->Tau2Pt = taus.at(1)->pt()/1e3;
    m_tree->Tau2Eta = taus.at(1)->eta();
    m_tree->Tau2Phi = taus.at(1)->phi();
    m_tree->Tau2M = taus.at(1)->m()/1e3;
  }

}

void AnalysisReader_hhbbtt::fill_sjetHistos(std::vector<const xAOD::Jet*> jets)
{
  m_histSvc->BookFillHist("NSignalJets", 10, 0, 10, jets.size(), m_weight);

  //Leading jet histos
  m_histSvc->BookFillHist("SignalJet0Eta",50,-5,5, jets.at(0)->eta(), m_weight);
  m_histSvc->BookFillHist("SignalJet0Phi",50,-5,5, jets.at(0)->phi(), m_weight);
  m_histSvc->BookFillHist("SignalJet0Pt",100,0,1000, jets.at(0)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("SignalJet0isTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(0)), m_weight);
  m_histSvc->BookFillHist("SignalJet0tagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(0)), m_weight);

  //Sub-Leading jet histos
  m_histSvc->BookFillHist("SignalJet1Eta",50,-5,5, jets.at(1)->eta(), m_weight);
  m_histSvc->BookFillHist("SignalJet1Phi",50,-5,5, jets.at(1)->phi(), m_weight);
  m_histSvc->BookFillHist("SignalJet1Pt",100,0,1000, jets.at(1)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("SignalJet1isTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(1)), m_weight);
  m_histSvc->BookFillHist("SignalJet1tagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(1)), m_weight);

  int nBJets=0;
  //Inclusive jet histos
  for(unsigned int i=0;i<jets.size();i++){
     m_histSvc->BookFillHist("SignalJetEta",50,-5,5, jets.at(i)->eta(), m_weight);
     m_histSvc->BookFillHist("SignalJetPhi",50,-5,5, jets.at(i)->phi(), m_weight);
     m_histSvc->BookFillHist("SignalJetPt",100,0,1000, jets.at(i)->pt()*1e-3, m_weight);
     m_histSvc->BookFillHist("SignalJetisTagged",2,-0.5,1.5, BTagProps::isTagged.get(jets.at(i)), m_weight);
     m_histSvc->BookFillHist("SignalJettagWeight",20,-2.0,2.0, BTagProps::tagWeight.get(jets.at(i)), m_weight);
     if(BTagProps::isTagged.get(jets.at(i))) nBJets++;
  }

  m_histSvc->BookFillHist("NJets", 10, 0, 10, jets.size(), m_weight);
  m_histSvc->BookFillHist("NBJets", 10, 0, 10, nBJets, m_weight);

}

std::string AnalysisReader_hhbbtt::get2DProngRegion(std::vector<const xAOD::TauJet*> taus)
{
   std::string string_out="";

   if(Props::nTracks.get(taus.at(0)) == 1 && Props::nTracks.get(taus.at(1)) == 3){
       string_out="nProngs13";
   }else if(Props::nTracks.get(taus.at(0)) == 3 && Props::nTracks.get(taus.at(1)) == 1){
       string_out="nProngs31";
   }else if(Props::nTracks.get(taus.at(0)) == 1 && Props::nTracks.get(taus.at(1)) == 1){
       string_out="nProngs11";
   }else if(Props::nTracks.get(taus.at(0)) == 3 && Props::nTracks.get(taus.at(1)) == 3){
       string_out="nProngs33";
   }

   return string_out;
}

void AnalysisReader_hhbbtt::fill_tauFFHistos(std::string regionFF,std::vector<const xAOD::TauJet*> taus)
{

  int NBinspT1=8;
  int NBinspT2=6;
  int NBinsEta2=10;
  
  float FF_pT1_Edges[NBinspT1+1]={40.0,50.0,60.0,70.0,80.0,100.0,140.0,240.0,1000.0};
  float FF_pT2_Edges[NBinspT2+1]={20.0,30.0,40.0,50.0,60.0,100.0,1000.0};
  float FF_eta_Edges[NBinsEta2+1]={-3.0,-2.4,-1.8,-1.2,-0.6,0.0,0.6,1.2,1.8,2.4,3.0};

  m_histSvc->BookFillHist("2DFFPt1Pt2", NBinspT1, FF_pT1_Edges, NBinspT2, FF_pT2_Edges, taus.at(0)->pt()/1e3, taus.at(1)->pt()/1e3, m_weight);
  m_histSvc->BookFillHist("2DFFPt1Eta1", NBinspT1, FF_pT1_Edges, NBinsEta2, FF_eta_Edges, taus.at(0)->pt()/1e3, taus.at(1)->eta(), m_weight);

  int tauTTag0=fabs(Props::TATTruthMatch.get(taus.at(0)));
  int tauTTag1=fabs(Props::TATTruthMatch.get(taus.at(1)));


  if( (m_isMC && tauTTag0 == 15 && tauTTag1 == 15) || !m_isMC){
     //if(m_region=="AM_OS" && regionFF=="nProngs11") std::cout  << "Tau 0 Pt: " << taus.at(0)->pt()/1e3 << "Tau 1 Pt: " << taus.at(1)->pt()/1e3 << std::endl;
     m_histSvc->BookFillHist(("2DFFPt1Pt2_"+regionFF).c_str(), NBinspT1, FF_pT1_Edges, NBinspT2, FF_pT2_Edges, taus.at(0)->pt()/1e3, taus.at(1)->pt()/1e3, m_weight);
  }

  m_histNameSvc->set_nProng(-1);
  m_histSvc->BookFillHist("FinalTau0nTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(0)),m_weight);
  m_histSvc->BookFillHist("FinalTau1nTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(1)),m_weight);

  m_histNameSvc->set_nProng(-1);
  m_histSvc->BookFillHist("2DFinalTauPt1Pt2", 100, 0, 500, 100, 0, 500, taus.at(1)->pt()/1e3, taus.at(0)->pt()/1e3, m_weight);


  if(m_isMC){

    m_histNameSvc->set_description(m_taupairCharge);
    m_histSvc->BookFillHist("BDTScore2D",50,0.35,1.0,50,0.35,1.0,Props::BDTScore.get(taus.at(0)),Props::BDTScore.get(taus.at(1)),1);

    //m_histNameSvc->set_description("AllReg");

    int tauTTag0=fabs(Props::TATTruthMatch.get(taus.at(0)));
    int tauTTag1=fabs(Props::TATTruthMatch.get(taus.at(0)));
    if(tauTTag0 == 15){
      m_histSvc->BookFillHist("BDTScore_FCompTau",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag0 == 1 || tauTTag0 == 2 || tauTTag0 == 3 ){
      m_histSvc->BookFillHist("BDTScore_FCompL",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag0 == 4){
      m_histSvc->BookFillHist("BDTScore_FCompC",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag0 == 5){
      m_histSvc->BookFillHist("BDTScore_FCompB",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag0 == 13 || tauTTag0 == 11 || tauTTag0 == 0){
      m_histSvc->BookFillHist("BDTScore_FCompOthers",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag0 == 21){
      m_histSvc->BookFillHist("BDTScore_FCompGluon",50,0.0,1,Props::BDTScore.get(taus.at(0)), 1);
    }
    if(tauTTag1 == 15){
      m_histSvc->BookFillHist("BDTScore_FCompTau",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    } 
    if(tauTTag1 == 1 || tauTTag1 == 2 || tauTTag1 == 3 ){ 
      m_histSvc->BookFillHist("BDTScore_FCompL",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    } 
    if(tauTTag1 == 4){
      m_histSvc->BookFillHist("BDTScore_FCompC",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    } 
    if(tauTTag1 == 5){
      m_histSvc->BookFillHist("BDTScore_FCompB",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    } 
    if(tauTTag1 == 13 || tauTTag1 == 11 || tauTTag1 == 0){ 
      m_histSvc->BookFillHist("BDTScore_FCompOthers",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    } 
    if(tauTTag1 == 21){ 
      m_histSvc->BookFillHist("BDTScore_FCompGluon",50,0.0,1,Props::BDTScore.get(taus.at(1)), 1);
    }

    m_histNameSvc->set_description(m_region);
    m_histSvc->BookFillHist("BDTScore2D",50,0.35,1.0,50,0.35,1.0,Props::BDTScore.get(taus.at(0)),Props::BDTScore.get(taus.at(1)),1);
  }

  //Leading jet histos
  m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(0)));
  m_histSvc->BookFillHist("FinalTau0Eta",30,-3,3, taus.at(0)->eta(), m_weight);
  m_histSvc->BookFillHist("FinalTau0Pt",100,0,1000, taus.at(0)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("FinalTau0Phi",30,-5,5, taus.at(0)->phi(), m_weight);
  m_histSvc->BookFillHist("FinalTau0BDTScore",20,0.35,1,Props::BDTScore.get(taus.at(0)), m_weight);
  if(m_isMC){
    int tauTTag=fabs(Props::TATTruthMatch.get(taus.at(0)));
    if(tauTTag == 15) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompTau",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 1) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompD",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 2) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompU",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 3) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompS",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 4) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompC",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 5) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompB",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 13) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompMuon",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 11) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompElectron",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 21) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompGluon",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
    if(tauTTag == 0) m_histSvc->BookFillHist("FinalTau0BDTScore_FCompNoMatched",20,0.35,1,Props::BDTScore.get(taus.at(0)), 1);
  }

  //Sub-Leading jet histos
  m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(1)));
  m_histSvc->BookFillHist("FinalTau1Eta",30,-3,3, taus.at(1)->eta(), m_weight);
  m_histSvc->BookFillHist("FinalTau1Pt",100,0,1000, taus.at(1)->pt()*1e-3, m_weight);
  m_histSvc->BookFillHist("FinalTau1Phi",30,-5,5, taus.at(1)->phi(), m_weight);
  m_histSvc->BookFillHist("FinalTau1BDTScore",20,0.35,1,Props::BDTScore.get(taus.at(1)), m_weight);
  if(m_isMC){
    int tauTTag=fabs(Props::TATTruthMatch.get(taus.at(1)));
    if(tauTTag == 15) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompTau",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 1) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompD",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 2) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompU",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 3) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompS",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 4) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompC",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 5) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompB",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 13) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompMuon",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 11) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompElectron",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 21) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompGluon",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
    if(tauTTag == 0) m_histSvc->BookFillHist("FinalTau1BDTScore_FCompNoMatched",20,0.35,1,Props::BDTScore.get(taus.at(1)), 1);
  }

  //Inclusive jet histos
  for(unsigned int i=0;i<taus.size();i++){
     m_histNameSvc->set_nProng(Props::nTracks.get(taus.at(i)));
     m_histSvc->BookFillHist("FinalTauEta",50,-5,5, taus.at(i)->eta(), m_weight);
     m_histSvc->BookFillHist("FinalTauPhi",50,-5,5, taus.at(i)->phi(), m_weight);
     m_histSvc->BookFillHist("FinalTauPt",100,0,1000, taus.at(i)->pt()*1e-3, m_weight);
     m_histSvc->BookFillHist("FinalTauBDTScore",20,0.35,1, Props::BDTScore.get(taus.at(i)), m_weight);
     m_histNameSvc->set_nProng(-1);
     m_histSvc->BookFillHist("FinalTaunTracks",5,-0.5,4.5,Props::nTracks.get(taus.at(i)),m_weight);
  }

}

EL::StatusCode AnalysisReader_hhbbtt::applyFF(std::string region,std::vector<const xAOD::TauJet*> finalTaus )
{

    int NBinspT1=8,NBinspT2=6;
    float FF_pT1_Edges[NBinspT1+1]={40.0,50.0,60.0,70.0,80.0,100.0,140.0,240.0,1000.0};
    float FF_pT2_Edges[NBinspT2+1]={20.0,30.0,40.0,50.0,60.0,100.0,1000.0};

    bool useNTracksInfo=true;
    if(useNTracksInfo){
      m_2DFF = (TH2F*)m_FF_File->Get(("fakes_"+region).c_str());
      if(m_debug) std::cout << "Getting 2D hist fakes for region " << region << std::endl;
    }else{ m_2DFF = (TH2F*)m_FF_File->Get(("fakes_"+region).c_str());}

    bool norder=true;
    float x,y;

    float tauPt1=finalTaus.at(0)->p4().Pt()*1e-3;
    float tauPt2=finalTaus.at(1)->p4().Pt()*1e-3;

    if(tauPt1< tauPt2) norder=false;
    
    if(norder) {
       x=tauPt1;
       y=tauPt2;
    }else{
       x=tauPt2;
       y=tauPt1;
    }

    Int_t binx = 0;
    if(m_debug) std::cout << "x (Pt_1): " << x << std::endl;
    binx=m_2DFF->GetXaxis()->FindBin(x);
    if(m_debug) std::cout << "x Bin: " << binx << std::endl;

    float fakef=1.0;
    for(unsigned int j=0;j<NBinspT2+1;j++){
       if(y > FF_pT2_Edges[j] && y< FF_pT2_Edges[j+1]){
          fakef=m_2DFF->GetBinContent(binx,j+1);
          //m_subtract = m_mcBkgA->GetBinContent(binx,j+1);
          if(m_debug) std::cout << "Fake factor: " << fakef << std::endl;
       }
    }

    m_weight *= fakef;

    if(m_isMC) m_weight *= -1;

    m_histNameSvc->set_sample("fakes");
    //if(m_region == "AA_OS") m_histNameSvc->set_description("MM_OS");
    //else if(m_region == "AA_SS") m_histNameSvc->set_description("MM_SS");
    if(m_RegC) m_histNameSvc->set_description("MM_OS");
    else if(m_RegD) m_histNameSvc->set_description("MM_SS");

    return EL::StatusCode::SUCCESS;

}

void AnalysisReader_hhbbtt::fill_triggerPlots()
{
    std::vector<float> L1_EMPhi  = Props::L1_EMPhi.get(m_eventInfo);
    std::vector<float> L1_EMEt   = Props::L1_EMEt.get(m_eventInfo);
    std::vector<float> L1_EMEta  = Props::L1_EMEta.get(m_eventInfo);
    std::vector<float> L1_JetEta = Props::L1_JetEta.get(m_eventInfo);
    std::vector<float> L1_JetPt  = Props::L1_JetPt.get(m_eventInfo);
    std::vector<float> L1_JetPhi = Props::L1_JetPhi.get(m_eventInfo);

    std::vector<float> tauPt;
    std::vector<float> tauPhi;
    std::vector<float> tauEta;

    std::vector<float> jetPt = L1_JetPt;

    std::vector<float> tauPtUnsorted;

    for(int j=0;j<jetPt.size();j++) {
      std::cout << "L1_JetPt: " << jetPt.at(j) << " L1_JetEta: " << L1_JetEta.at(j) << " L1_JetPhi: "<< L1_JetPhi.at(j)<< std::endl;
    }
    //std::cout << "here again" << std::endl;
    for(int j=0;j<L1_EMEt.size();j++) {
      if((L1_EMPhi.at(j)==0) && (L1_EMEt.at(j)==0) && (L1_EMEta.at(j)==0)) continue;
      tauPt.push_back(L1_EMEt.at(j));
      tauPhi.push_back(L1_EMPhi.at(j));
      tauEta.push_back(L1_EMEta.at(j));
      tauPtUnsorted.push_back(L1_EMEt.at(j));
      std::cout << "L1_TauPt: " << L1_EMEt.at(j) << " L1_TauEta: " << L1_EMEta.at(j) << " L1_TauPhi: "<< L1_EMPhi.at(j)<< std::endl;
    }

    if(tauPhi.size()<2 && L1_JetPhi.size()<2) return;

    float deltaR=9999.9;
    int j=0;

    std::cout << "L1PtJet size before " << L1_JetPt.size() << std::endl;

    for(int k=0;k<tauPt.size();k++) {
      //std::cout << "tauPt: " << tauPt.at(k) << " tauPhi: " << tauPhi.at(k) << " tauEta: " << tauEta.at(k) << std::endl;
      //for(vector<float>::iterator it = L1_JetPt.begin(); it != L1_JetPt.end(); ) {
      for(int j=0;j<L1_JetPt.size();j++){
        //std::cout << "L1_JetPt: " << L1_JetPt.at(j) << " L1_JetPhi: " << L1_JetPhi.at(j) << " L1_JetEta: " << L1_JetEta.at(j) << std::endl;
        deltaR=sqrt(pow(L1_JetPhi.at(j)-tauPhi.at(k),2)+pow(L1_JetEta.at(j)-tauEta.at(k),2));
        //std::cout << "Delta R: " << deltaR << " (k,j) " << k << " , " << j << std::endl;
        if(deltaR<0.4) {
          std::cout << "Ovelap found with Delta R: " << deltaR << std::endl;
          std::cout << "Matched jet pT: " << L1_JetPt.at(j) << " jetPt size: " << jetPt.size() << std::endl;
          if(jetPt.size() ==0) break;
          for(vector<float>::iterator it = jetPt.begin(); it != jetPt.end(); ) {
             //std::cout << "just erasing " << std::endl;
             if(*it == L1_JetPt.at(j)) it = jetPt.erase(it);
             else ++it;
          }
	}
      }
      if(jetPt.size() ==0) break;
      //++it;
    }

    std::cout << "L1PtJet size after " << jetPt.size() << " TauPt size: " << tauPt.size() << std::endl;

    std::sort(tauPt.begin(),tauPt.end());

    std::vector<int> tau_index, jet_index;
    vector<float>::iterator it1;

    for(int j=tauPt.size()-1;j>=tauPt.size()-2;j--) {
      if(j<0) {/*std::cout << "j: " << j << std::endl;*/ break;}
      std::cout << "tauPt: " << tauPt.at(j) << " j: " << j <<std::endl;
      it1=find(tauPtUnsorted.begin(),tauPtUnsorted.end(),tauPt.at(j));
      int pos = distance(tauPtUnsorted.begin(), it1);
      tau_index.push_back(pos);
      std::cout << "Position found at: " << pos << std::endl;
    }

    //std::cout << "After tau loop" << std::endl;
    int tauv_size=tauPt.size()-1;
    m_histSvc->BookFillHist("L1_lead_tau_pt",50,0,250,tauPt.at(tauv_size)*1e-3, 1.0);
    m_histSvc->BookFillHist("L1_sublead_tau_pt",50,0,250,tauPt.at(tauv_size-1)*1e-3, 1.0);

    TVector2 tau0(tauEta.at(tau_index.at(0)),tauPhi.at(tau_index.at(0)));
    TVector2 tau1(tauEta.at(tau_index.at(1)),tauPhi.at(tau_index.at(1)));
    float dPhi=tau0.DeltaPhi(tau1);
    float dR=(tau0-tau1).Mod();
    m_histSvc->BookFillHist("L1_dphitautau",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    m_histSvc->BookFillHist("L1_dRtautau",25,0,5,dR, 1.0);

    if(jetPt.size()>1){
      std::sort(jetPt.begin(),jetPt.end());
      for(int j=jetPt.size()-1;j>=jetPt.size()-2;j--) {
         if(j<0) {/*std::cout << "j: " << j << std::endl;*/ break;}
         std::cout << "jetPt: " << jetPt.at(j) << std::endl;
         it1=find(L1_JetPt.begin(),L1_JetPt.end(),jetPt.at(j));
         int pos = distance(L1_JetPt.begin(), it1);
         jet_index.push_back(pos);
         std::cout << "jet Position found at: " << pos << std::endl;
      }
    }else if(jetPt.size()<2) return;

    int jetv_size=jetPt.size()-1;
    m_histSvc->BookFillHist("L1_lead_jet_pt",50,0,250,jetPt.at(jetv_size)*1e-3, 1.0);
    m_histSvc->BookFillHist("L1_sublead_jet_pt",50,0,250,jetPt.at(jetv_size-1)*1e-3, 1.0);

    TVector2 jet0(L1_JetEta.at(jet_index.at(0)),L1_JetPhi.at(jet_index.at(0)));
    TVector2 jet1(L1_JetEta.at(jet_index.at(1)),L1_JetPhi.at(jet_index.at(1)));

    dPhi=jet0.DeltaPhi(jet1);
    dR=(jet0-jet1).Mod();
    m_histSvc->BookFillHist("L1_dphijetjet",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    m_histSvc->BookFillHist("L1_dRjetjet",25,0,5,dR, 1.0);
    
    dPhi=jet0.DeltaPhi(tau0);
    m_histSvc->BookFillHist("L1_dphijet0tau0",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    dPhi=jet0.DeltaPhi(tau1);
    m_histSvc->BookFillHist("L1_dphijet0tau1",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    dPhi=jet1.DeltaPhi(tau0);
    m_histSvc->BookFillHist("L1_dphijet1tau0",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    dPhi=jet1.DeltaPhi(tau1);
    m_histSvc->BookFillHist("L1_dphijet1tau1",25,-TMath::Pi(),TMath::Pi(),dPhi, 1.0);
    dR=(jet0-tau0).Mod();
    m_histSvc->BookFillHist("L1_dRjet0tau0",25,0,5,dR, 1.0);
    dR=(jet0-tau1).Mod();
    m_histSvc->BookFillHist("L1_dRjet0tau1",25,0,5,dR, 1.0);
    dR=(jet1-tau0).Mod();
    m_histSvc->BookFillHist("L1_dRjet1tau0",25,0,5,dR, 1.0);
    dR=(jet1-tau1).Mod();
    m_histSvc->BookFillHist("L1_dRjet1tau1",25,0,5,dR, 1.0);

    //std::cout << "Here at last" << std::endl;

    return;

    //int jetv_size=jetPt.size()-1;
    //m_histSvc->BookFillHist("L1_lead_tau_pt",50,0,250,tauPt.at(jetv_size)*1e-3, 1.0);
    //m_histSvc->BookFillHist("L1_sublead_tau_pt",50,0,250,tauPt.at(jetv_size-1)*1e-3, 1.0);

    /*std::map<int,int> indexMap;

    for(int j=0;j<tauPhi.size();j++) {
       for(int k=0;k<tauPhi.size();k++){
          if(j==k) continue;
          
       }
    }

    for(int j=0;j<L1_JetPhi.size();j++) {
       for(int k=0;k<L1_JetPhi.size();k++){
          if(j==k) continue;
       }
    }

    
    for(int j=0;j<L1_JetPhi.size();j++) {
       for(int k=0;k<tauPhi.size();k++){
          
       }
    }*/

}


float AnalysisReader_hhbbtt::getMT2(std::vector<TLorentzVector> bjets, TLorentzVector METv_new)
{
  mt2_bisect::mt2 mt2_event;
  double pa[3] = {bjets.at(0).M(), bjets.at(0).Px(), bjets.at(0).Py()};  // visible particle #1: bjets.at(0)
  double pb[3] = {bjets.at(1).M(), bjets.at(1).Px(), bjets.at(1).Py()};  // visible particle #2: bjets.at(1)

  double pmiss[3] = {0.0, METv_new.Px(), METv_new.Py()};  // MET_new = MET + tau + tau
  double mn = 80398.0;  // mass of W in MeV
  mt2_event.set_momenta(pa, pb, pmiss);
  mt2_event.set_mn(mn);

  return mt2_event.get_mt2();

}

float AnalysisReader_hhbbtt::mTsqr(TLorentzVector A, TLorentzVector B)
{
  float mTsqr = 2*A.Pt()*B.Pt()*(1-cos(A.Phi()-B.Phi()))*1e-6;
  return mTsqr;
}

int AnalysisReader_hhbbtt::getSampleNumber()
{
  int sample_number=-1;
  std::map<std::string,int> sample_nmap;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M260"]=0;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M300"]=1;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M400"]=2;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M500"]=3;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M600"]=4;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M700"]=5;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M800"]=6;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M900"]=7;
  sample_nmap["RS_G_hh_bbtt_hh_c10_M1000"]=8;
  sample_nmap["data"]=9;

  std::map<std::string, int>::iterator it = sample_nmap.begin();
  for( ; it != sample_nmap.end(); ++it ){
    std::string sampleName=m_histNameSvc->getFullSample();
    std::string testSName =it->first;
    //std::cout << "testSName: " << testSName << " sampleName: " << sampleName << std::endl;
    if(sampleName == testSName) sample_number=it->second;
  }

  return sample_number;
}
