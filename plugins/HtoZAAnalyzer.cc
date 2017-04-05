#include <cp3_llbb/HtoZAAnalysis/interface/HtoZAAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/HtoZAAnalysis/interface/Categories.h>
#include <cp3_llbb/HtoZAAnalysis/interface/GenStatusFlags.h>

#include <cp3_llbb/Framework/interface/EventProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#include <TH1F.h>
#include <TCanvas.h>

#define HtoZA_GEN_DEBUG (true)

void HtoZAAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    edm::ParameterSet newconfig = edm::ParameterSet(config);
    newconfig.addUntrackedParameter("m_analyzer_name", this->m_name);
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", newconfig);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", newconfig);
}



void HtoZAAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const AnalyzersManager&, const CategoryManager&) {

    //Reset event
    leptons.clear();
    ll.clear();
    jj.clear();
    lljj.clear();
    

    const JetsProducer& alljets = producers.get<JetsProducer>(m_jets_producer);
    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>(m_electrons_producer);
    const MuonsProducer& allmuons = producers.get<MuonsProducer>(m_muons_producer);
    const EventProducer& fwevent = producers.get<EventProducer>("event");
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    //const METProducer& pf_met = producers.get<METProducer>(m_met_producer);


    if(!event.isRealData()){
    
    //BR thing? Maybe in HH it applies to W;
    //Does it have to be taken into account for Z as well?
    size_t n_taus = 0;
    const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");

#if HtoZA_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gp, &print_mother_chain](size_t p) {
    
        if (gp.pruned_mothers_index[p].empty()) {
            std::cout << std::endl;
            return;
        }
 
        // Get the index in gp.pruned_p4 of the particle that is a mother
        size_t index = gp.pruned_mothers_index[p][0];
            std::cout << " <- #" << index << "(" << gp.pruned_pdg_id[index] << ")" << std::endl;
            print_mother_chain(index);
    };
#endif


    std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gp](size_t particle_index, size_t mother_index) -> bool {
      // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
      if (gp.pruned_mothers_index[particle_index].empty())
          return false;

      size_t index = gp.pruned_mothers_index[particle_index][0];

      if (index == mother_index) {
          return true;
      }

      if(pruned_decays_from(index, mother_index))
          return true;

      return false;
    };

    std::function<bool(size_t, size_t, bool)> pruned_decays_from_pdg_id = [&pruned_decays_from_pdg_id, &gp](size_t particle_index, uint64_t pdg_id, bool direct) -> bool {
      // Iterator over all pruned particles to find if the particle `particle_index` decays from a particle with pdg id == pdg_id
      if (gp.pruned_mothers_index[particle_index].empty())
          return false;

      size_t index = gp.pruned_mothers_index[particle_index][0];

      if (std::abs(gp.pruned_pdg_id[index]) == pdg_id) {
          return true;
      }

      if (!direct && pruned_decays_from_pdg_id(index, pdg_id, direct))
          return true;

      return false;
    };


    // Construct signal gen info
    std::cout << "New event" << std::endl;
    
    gen_iH = -1;
    gen_iA = -1;
    gen_iA_afterFSR = -1;
    gen_iZ = -1;
    gen_iZ_afterFSR = -1;
    gen_ih = -1;
    gen_ih_afterFSR = -1;
    gen_iB = gen_iBbar = -1;
    gen_iB_afterFSR = gen_iBbar_afterFSR = -1;
    gen_iW = -1;
    gen_iW_afterFSR = -1;
    gen_iLminus = gen_iLplus = -1;
    gen_iLminus_afterFSR = gen_iLplus_afterFSR = -1;
    gen_iNu1 = gen_iNu2 = -1;

    for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
        std::bitset<15> flags (gp.pruned_status_flags[ip]);

        if (!flags.test(8))
            continue;

        int64_t pdg_id = gp.pruned_pdg_id[ip];
        int8_t status = gp.pruned_status[ip];

#if HtoZA_GEN_DEBUG
        //std::cout << "[" << ip << "] pdg id: " << pdg_id << " flags: " << flags << " p = " << gp.pruned_p4[ip] << std::endl;
        std::cout << "[" << ip << "] pdg id: " << pdg_id << " flags: " << flags << " p = " << gp.pruned_p4[ip] << " status: " << (int) status << std::endl;
        //std::cout << "gp.pruned_pdg_id.size(): " << gp.pruned_pdg_id.size() << std::endl;
        //std::cout << "gp.pruned_p4.size(): " << gp.pruned_p4.size() << std::endl;
        print_mother_chain(ip);
#endif

        auto p4 = gp.pruned_p4[ip];
        //std::cout << "gen_iZ: " << (int) gen_iZ << ", gen_iZ_afterFSR: " << (int) gen_iZ_afterFSR << std::endl;
        if (std::abs(pdg_id) == 35 || std::abs(pdg_id) == 39) {
            ASSIGN_HtoZA_GEN_INFO_NO_FSR(H, "H");
            //std::cout << "ONE H FOUND WITH gen_iH: " << (int) gen_iH << std::endl;
        } else if (pdg_id == 36) {
            ASSIGN_HtoZA_GEN_INFO(A, "A");
        } else if (pdg_id == 23) {
            ASSIGN_HtoZA_GEN_INFO(Z, "Z"); 
            std::cout << " AFTER ASSIGNEMENT: gen_iZ: " << (int) gen_iZ << std::endl;
        } else if (pdg_id == 25) {
            ASSIGN_HtoZA_GEN_INFO(h, "SM Higgs");
        }

        // Only look for A and Z decays (not for the H decay because it's often off-shell) -- check it
        if ((gen_iA == -1) || (gen_iZ == -1))
            continue;

        // And if the particles come directly from Z and A decays (and H???)
        bool from_Z_decay = pruned_decays_from(ip, gen_iZ);
        bool from_A_decay = pruned_decays_from(ip, gen_iA);

        // Only keep particles coming from Z and A decays
        if (!from_Z_decay && !from_A_decay)
            continue;

        if (pdg_id == 5) {
            ASSIGN_HtoZA_GEN_INFO(B, "B");
        } else if (pdg_id == -5) {
            ASSIGN_HtoZA_GEN_INFO(Bbar, "Bbar");        
        }

        // Ignore B decays: why false?????
        if (pruned_decays_from_pdg_id(ip, 5, false))
            continue;

        // Count the number of tau coming directly from a W or a Z
        if ((std::abs(pdg_id) == 15) && (pruned_decays_from_pdg_id(ip, 24, true) || pruned_decays_from_pdg_id(ip, 23, true))) {
            n_taus++;
        }

        // Electron, muon, tauon
        if ((pdg_id == 11) || (pdg_id == 13) || (pdg_id == 15)) {
            ASSIGN_HtoZA_GEN_INFO(Lminus, "L-");
            //std::cout << "AFTER ASSIGNEMENT: gen_iLminus = " << (int) gen_iLminus << std::endl;
        } else if ((pdg_id == -11) || (pdg_id == -13) || (pdg_id == -15)) {
            ASSIGN_HtoZA_GEN_INFO(Lplus, "L+");
            //std::cout << "AFTER ASSIGNEMENT: gen_iLplus = " << (int) gen_iLplus << std::endl;
        } else if (std::abs(pdg_id) == 24) {
            ASSIGN_HtoZA_GEN_INFO(W, "W");
        } else if ((std::abs(pdg_id) == 12) || (std::abs(pdg_id) == 14) || (std::abs(pdg_id) == 16)) {
            ASSIGN_HtoZA_GEN_INFO_2_NO_FSR(Nu1, Nu2, "neutrinos");
          }
    } // End of loop

    if (gen_iH == -1 ) {
        offshellH_counter++;
        std::cout << "NO H IN THE EVENT!!!!!!!!!!!" << std::endl;
    }
        std::cout << "offshellH_counter: " << offshellH_counter << std::endl;
        events_counter++;
        std::cout << "events_counter: " << events_counter << std::endl;

    std::cout << "gen_iLminus: " << (int) gen_iLminus << std::endl;
    std::cout << "gen_iLplus: " << (int) gen_iLplus << std::endl;


      // Get the invariant mass of Z+A
      if ((gen_iZ != -1) && (gen_iA != -1)) {
          gen_mZA = (gen_Z + gen_A).M();
          gen_costhetastar = getCosThetaStar_CS(gen_Z, gen_A);
      }

 #if HtoZA_GEN_DEBUG
      std::cout << "*************************" << std::endl;
      PRINT_PARTICULE(H);
      PRINT_RESONANCE(Z, A);
      PRINT_RESONANCE(B, Bbar);
      PRINT_RESONANCE(Lminus, Lplus);
      PRINT_RESONANCE_NO_FSR(Nu1, Nu2);
      std::cout << "*************************" << std::endl;

      //Rebuild resonances for consistency checks
      auto LplusLminus = gen_Lplus + gen_Lminus;
      std::cout << "    gen_(L+ L-).M() = " << LplusLminus.M() << std::endl;

      auto LplusLminus_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR;
      std::cout << "    gen_(L+ L-)_afterFSR.M() = " << LplusLminus_afterFSR.M() << std::endl;

      auto BBbar = gen_B + gen_Bbar;
      std::cout << "    gen_(B Bbar).M() = " << BBbar.M() << std::endl;

      auto BBbar_afterFSR = gen_B_afterFSR + gen_Bbar_afterFSR;
      std::cout << "    gen_(B Bbar)_afterFSR.M() = " << BBbar_afterFSR.M() << std::endl;

      auto LLBB = gen_Lplus + gen_Lminus + gen_B + gen_Bbar;
      std::cout << "    gen_(L+ L- B Bbar).M() = " << LLBB.M() << std::endl;

      auto LLBB_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR + gen_B_afterFSR + gen_Bbar_afterFSR;
      std::cout << "    gen_(L+ L- B Bbar)_afterFSR.M() = " << LLBB_afterFSR.M() << std::endl;
      std::cout << "+++++++++++++++++++++++++++" << std::endl;

#endif


        // ***** ***** *****
        // Matching
        // ***** ***** *****
        for (auto p4: alljets.gen_p4) {
            gen_deltaR_jet_B.push_back(deltaR(p4, gen_B));
            gen_deltaR_jet_Bbar.push_back(deltaR(p4, gen_Bbar));
            gen_deltaR_jet_B_afterFSR.push_back(deltaR(p4, gen_B_afterFSR));
            gen_deltaR_jet_Bbar_afterFSR.push_back(deltaR(p4, gen_Bbar_afterFSR));
        }
        for (auto p4: allelectrons.gen_p4) {
            gen_deltaR_electron_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_electron_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_electron_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_electron_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }
        for (auto p4: allmuons.gen_p4) {
            gen_deltaR_muon_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_muon_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_muon_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_muon_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }
    }


    LorentzVector null_p4(0., 0., 0., 0.);
    float event_weight = fwevent.weight;
    float tmp_count_has2leptons = 0.;
    float tmp_count_has2leptons_elel = 0.;
    float tmp_count_has2leptons_elmu = 0.;
    float tmp_count_has2leptons_muel = 0.;
    float tmp_count_has2leptons_mumu = 0.;
    float tmp_count_has2leptons_1lljj = 0.;
    float tmp_count_has2leptons_elel_1lljj = 0.;
    float tmp_count_has2leptons_elmu_1lljj = 0.;
    float tmp_count_has2leptons_muel_1lljj = 0.;
    float tmp_count_has2leptons_mumu_1lljj = 0.;
    float tmp_count_has2leptons_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_elel_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_elmu_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_muel_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_mumu_1lljj_2btagM = 0.;

    // ***** ***** *****
    // Trigger Matching
    // ***** ***** *****

    // the actual trigger matching to dilepton HLT paths happens only once we have a dilepton candidate to consider in the event
    
    // ********** 
    // Leptons and dileptons
    // ********** 



    // Define a function that is true whenever the electron passes the trigger
    static auto electron_pass_HLT_ID = [&allelectrons, this](size_t index) {
        auto electron = allelectrons.products[index];

        // Use POG HLT-safe id
        
        bool result = allelectrons.ids[index][m_electron_hlt_safe_wp_name];

        // Add dxy and dz cuts described at https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
        if (electron->isEB()) {
            result &= std::abs(allelectrons.dz[index]) < 0.1;
            result &= std::abs(allelectrons.dxy[index]) < 0.05;
        } else {
            result &= std::abs(allelectrons.dz[index]) < 0.2;
            result &= std::abs(allelectrons.dxy[index]) < 0.1;
        }

        return result;
    };


    // Fill lepton structures
    // allelectrons is a vector with the reco electrons but it also contains gen info (in allelectrons.gen_p4)
    for (unsigned int ielectron = 0; ielectron < allelectrons.p4.size(); ielectron++)
    {
        if (allelectrons.p4[ielectron].Pt() > m_subleadingElectronPtCut && fabs(allelectrons.p4[ielectron].Eta()) < m_electronEtaCut)
        {
            // some selection
            // Ask for medium ID
            if (!allelectrons.ids[ielectron][m_electron_medium_wp_name])
                continue;

            HtoZA::Lepton ele;
            ele.p4 = allelectrons.p4[ielectron];
            ele.charge = allelectrons.charge[ielectron];
            ele.idx = ielectron;
            ele.isMu = false;
            ele.isEl = true;
            ele.ele_hlt_id = electron_pass_HLT_ID(ielectron);

            ele.gen_matched = allelectrons.matched[ielectron];
            ele.gen_p4 = ele.gen_matched ? allelectrons.gen_p4[ielectron] : null_p4;
            ele.gen_DR = ele.gen_matched ? ROOT::Math::VectorUtil::DeltaR(ele.p4, ele.gen_p4): -1.;
            ele.gen_DPtOverPt = ele.gen_matched ? (ele.p4.Pt() - ele.gen_p4.Pt()) / ele.p4.Pt() : -10.;
            ele.hlt_leg1 = false;
            ele.hlt_leg2 = false;

            ele.sc_eta = allelectrons.products[ielectron]->superCluster()->eta();

            leptons.push_back(ele);
        }
    } // end of loop over electrons

    for (unsigned int imuon = 0; imuon < allmuons.p4.size(); imuon++)
    {
        if (allmuons.p4[imuon].Pt() > m_subleadingMuonPtCut && fabs(allmuons.p4[imuon].Eta()) < m_muonEtaCut)
        {
            // Ask for tight ID & tight ISO
            if (!allmuons.isTight[imuon] || allmuons.relativeIsoR04_deltaBeta[imuon] >= m_muonTightIsoCut)
                continue;

            HtoZA::Lepton mu;
            mu.p4 = allmuons.p4[imuon];
            mu.charge = allmuons.charge[imuon];
            mu.idx = imuon;
            mu.isMu = true;
            mu.isEl = false;
            mu.gen_matched = allmuons.matched[imuon];
            mu.gen_p4 = mu.gen_matched ? allmuons.gen_p4[imuon] : null_p4;
            mu.gen_DR = mu.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mu.p4, mu.gen_p4) : -1.;
            mu.gen_DPtOverPt = mu.gen_matched ? (mu.p4.Pt() - mu.gen_p4.Pt()) / mu.p4.Pt() : -10.;
            mu.hlt_leg1 = false;
            mu.hlt_leg2 = false;
        
            leptons.push_back(mu);
        }
    } // end of loop over muons



    // sort leptons by pt (ignoring flavour, id and iso)
    std::sort(leptons.begin(), leptons.end(), [](const HtoZA::Lepton& lep1, const HtoZA::Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt();});

    for (unsigned int ilep1 = 0; ilep1 < leptons.size(); ilep1++)
    {
        if ((leptons[ilep1].isMu && leptons[ilep1].p4.Pt() < m_leadingMuonPtCut) || (leptons[ilep1].isEl && leptons[ilep1].p4.Pt() < m_leadingElectronPtCut))
            continue;

        for (unsigned int ilep2 = ilep1+1; ilep2 < leptons.size(); ilep2++)
        {
            HtoZA::Dilepton dilep;
            dilep.p4 = leptons[ilep1].p4 + leptons[ilep2].p4;
            // indices for dilepton: simply the pair (ilep1, ilep2)
            dilep.idxs = std::make_pair(leptons[ilep1].idx, leptons[ilep2].idx);
            dilep.ilep1 = ilep1;
            dilep.ilep2 = ilep2;
            // if the two leptons have opposite sign
            dilep.isOS = leptons[ilep1].charge * leptons[ilep2].charge < 0;
            dilep.isPlusMinus = leptons[ilep1].charge > 0 && leptons[ilep2].charge < 0;
            dilep.isMinusPlus = leptons[ilep1].charge < 0 && leptons[ilep2].charge > 0;
            dilep.isMuMu = leptons[ilep1].isMu && leptons[ilep2].isMu;
            dilep.isElEl = leptons[ilep1].isEl && leptons[ilep2].isEl;
            dilep.isElMu = leptons[ilep1].isEl && leptons[ilep2].isMu;
            dilep.isMuEl = leptons[ilep1].isMu && leptons[ilep2].isEl;
            // if the two leptons have same flavour
            dilep.isSF = dilep.isMuMu || dilep.isElEl;
            dilep.DR_l_l = ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4);
            dilep.DPhi_l_l = fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            dilep.ht_l_l = leptons[ilep1].p4.Pt() + leptons[ilep2].p4.Pt();
            // We are selecting all the possible combinations of dileptons in an event. But which one is the
            // one that triggered the event? You set this by calling the function "matchOfflineLepton()".
            if (!hlt.paths.empty()) {
                matchOfflineLepton(hlt, dilep);
                dilep.hlt_idxs = std::make_pair(leptons[dilep.ilep1].hlt_idx, leptons[dilep.ilep2].hlt_idx);
            }
            dilep.gen_matched = leptons[ilep1].gen_matched && leptons[ilep2].gen_matched;
            dilep.gen_p4 = dilep.gen_matched ? leptons[ilep1].gen_p4 + leptons[ilep2].gen_p4 : null_p4;
            dilep.gen_DR = dilep.gen_matched ? ROOT::Math::VectorUtil::DeltaR(dilep.p4, dilep.gen_p4) : -1.;
            dilep.gen_DPtOverPt = dilep.gen_matched ? (dilep.p4.Pt() - dilep.gen_p4.Pt()) / dilep.p4.Pt() : -10.;


            if (event.isRealData()) {
                dilep.trigger_efficiency = 1.;
                dilep.trigger_efficiency_downVariated = 1.;
                dilep.trigger_efficiency_upVariated = 1.;
            } else {
                fillTriggerEfficiencies(leptons[ilep1], leptons[ilep2], dilep); 
            }
            
            // Some selection
            // Note that ID and isolation criteria are in both electron and muon loops
            if (!dilep.isOS)
                continue;
            // Should I have also a "if (!dilep.isSF) continue;" ??

            // FIXME L1 EMTF bug mitigation -- cut the overlap on data if it's a run affected by the bug
            // On MC, apply the fraction of lumi the bug was not present
            if (dilep.isMuMu && isCSCWithOverlap(leptons[ilep1], leptons[ilep2])) {
                if (event.isRealData() && fwevent.run < 278167) {
                    continue;
                } else if (!event.isRealData()) {
                    dilep.trigger_efficiency *= 0.5265;
                    dilep.trigger_efficiency_downVariated *= 0.5265;
                    dilep.trigger_efficiency_upVariated *= 0.5265;
                }
            }


            // Throw event if there is no matched dilepton trigger path (only on data)
            if (event.isRealData() && !((leptons[dilep.ilep1].hlt_leg1 && leptons[dilep.ilep2].hlt_leg2)
                                  || (leptons[dilep.ilep1].hlt_leg2 && leptons[dilep.ilep2].hlt_leg1))) {
                continue;
            }

            // Counters
            tmp_count_has2leptons = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_elel = event_weight;
            if (dilep.isElMu)
                tmp_count_has2leptons_elmu = event_weight;
            if (dilep.isMuEl)
                tmp_count_has2leptons_muel = event_weight;
            if (dilep.isMuMu)
                tmp_count_has2leptons_mumu = event_weight;

            // Fill
            ll.push_back(dilep);
        } //end loop with lep2
    } //end loop with lep1

    // have the ll collection sorted by ht:
    // "sort" outputs the sorted vector
    std::sort(ll.begin(), ll.end(), [&](HtoZA::Dilepton& a, HtoZA::Dilepton& b){return a.ht_l_l > b.ht_l_l;});

    // Keep only the first ll candidate
    if (ll.size() > 1) {
        ll.resize(1); // discard all the elements beyond the first one
    }
    

    // ***** 
    // Jets and dijets 
    // *****

    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        float correctionFactor = m_applyBJetRegression ? alljets.regPt[ijet] / alljets.p4[ijet].Pt() : 1.;

        if ((alljets.p4[ijet].Pt() * correctionFactor > m_jetPtCut) && (fabs(alljets.p4[ijet].Eta()) < m_jetEtaCut))
        {
        
            if (!alljets.passLooseID[ijet])
                continue;

            HtoZA::Jet myjet;
            myjet.p4 = alljets.p4[ijet] * correctionFactor;
            myjet.idx = ijet;

            myjet.CSV = alljets.getBTagDiscriminant(ijet, "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            myjet.CMVAv2 = alljets.getBTagDiscriminant(ijet, "pfCombinedMVAV2BJetTags");
            float mybtag = alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName); // e.g. m_jet_bDiscrName = 'CSV'
            //myjet.btagL = mybtag > m_jet_bDiscrCut_loose;
            myjet.btag_M = mybtag > m_jet_bDiscrCut_medium;
            //myjet.btagT = mybtag > m_jet_bDiscrCut_tight;
            myjet.gen_matched_bParton = (std::abs(alljets.partonFlavor[ijet]) == 5);
            myjet.gen_matched_bHadron = (alljets.hadronFlavor[ijet]) == 5;
            myjet.gen_matched = alljets.matched[ijet];
            myjet.gen_p4 = myjet.gen_matched ? alljets.gen_p4[ijet] : null_p4;
            myjet.gen_DR = myjet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjet.p4, myjet.gen_p4) : -1.;
            myjet.gen_DPtOverPt = myjet.gen_matched ? (myjet.p4.Pt() - myjet.gen_p4.Pt()) / myjet.p4.Pt() : -10.;
            myjet.gen_b = (alljets.hadronFlavor[ijet]) == 5; // redundant with gen_matched_bHadron defined above
            myjet.gen_c = (alljets.hadronFlavor[ijet]) == 4;
            //light jets:
            myjet.gen_l = (alljets.hadronFlavor[ijet]) < 4;

            // We don't want a lepton close to the jet, right??????
            bool isThereACloseSelectedLepton = false;
            for (auto& mylepton: leptons) {
                if (ROOT::Math::VectorUtil::DeltaR(myjet.p4, mylepton.p4) < m_minDR_l_j_Cut) {
                    isThereACloseSelectedLepton = true;
                    break;
                }
            }

            if (isThereACloseSelectedLepton)
                continue;


            jets.push_back(myjet);
        } //end if 
    }// end loop over jets

    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ijet1 = 0; ijet1 < jets.size(); ijet1++)
    {
        for (unsigned int ijet2 = ijet1 + 1; ijet2 < jets.size(); ijet2++)
        {
            HtoZA::Dijet myjj;
            myjj.p4 = jets[ijet1].p4 + jets[ijet2].p4;
            myjj.idxs = std::make_pair(jets[ijet1].idx, jets[ijet2].idx);
            myjj.ijet1 = ijet1;
            myjj.ijet2 = ijet2;
            // Commented things
            // What's btag_M?????????????????????????????????
            myjj.btag_MM = jets[ijet1].btag_M && jets[ijet2].btag_M;
            // Commented things
            myjj.sumCSV = jets[ijet1].CSV + jets[ijet2].CSV;
            myjj.sumCMVAv2 = jets[ijet1].CMVAv2 + jets[ijet2].CMVAv2;
            myjj.DR_j_j = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet2].p4);
            myjj.DPhi_j_j = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[ijet1].p4, jets[ijet2].p4));
            myjj.ht_j_j = jets[ijet1].p4.Pt() + jets[ijet2].p4.Pt();
            myjj.gen_matched_bbPartons = jets[ijet1].gen_matched_bParton && jets[ijet2].gen_matched_bParton;
            myjj.gen_matched_bbHadrons = jets[ijet1].gen_matched_bHadron && jets[ijet2].gen_matched_bHadron;
            myjj.gen_matched = jets[ijet1].gen_matched && jets[ijet2].gen_matched;
            myjj.gen_p4 = myjj.gen_matched ? jets[ijet1].gen_p4 + jets[ijet2].gen_p4 : null_p4;
            myjj.gen_DR = myjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjj.p4, myjj.gen_p4) : -1.;
            myjj.gen_DPtOverPt = myjj.gen_matched ? (myjj.p4.Pt() - myjj.gen_p4.Pt()) / myjj.p4.Pt() : -10.;
            myjj.gen_bb = (jets[ijet1].gen_b && jets[ijet2].gen_b);
            myjj.gen_bc = (jets[ijet1].gen_b && jets[ijet2].gen_c) || (jets[ijet1].gen_c && jets[ijet2].gen_b);
            myjj.gen_bl = (jets[ijet1].gen_b && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_b);
            myjj.gen_cc = (jets[ijet1].gen_c && jets[ijet2].gen_c);
            myjj.gen_cl = (jets[ijet1].gen_c && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_c);
            myjj.gen_ll = (jets[ijet1].gen_l && jets[ijet2].gen_l);
            jj.push_back(myjj);
        } //end loop wtih jet2
    } //end loop with jet1

    // have the jj collection sorted by ht
    std::sort(jj.begin(), jj.end(), [&](HtoZA::Dijet& a, HtoZA::Dijet& b){return a.p4.Pt() > b.p4.Pt();});


    // ********** 
    // lljj, llbb
    // **********
    for (unsigned int ill = 0; ill < ll.size(); ill++)
    {
        for (unsigned int ijj = 0; ijj < jj.size(); ijj++)
        {
            unsigned int ijet1 = jj[ijj].ijet1;
            unsigned int ijet2 = jj[ijj].ijet2;
            unsigned int ilep1 = ll[ill].ilep1;
            unsigned int ilep2 = ll[ill].ilep2;
            HtoZA::DileptonDijet mylljj;
            mylljj.p4 = ll[ill].p4 + jj[ijj].p4;
            mylljj.lep1_p4 = leptons[ilep1].p4;
            mylljj.lep2_p4 = leptons[ilep2].p4;
            mylljj.jet1_p4 = jets[ijet1].p4;
            mylljj.jet2_p4 = jets[ijet2].p4;
            mylljj.ll_p4 = ll[ill].p4;
            mylljj.jj_p4 = jj[ijj].p4;
            // gen info
            mylljj.gen_matched = ll[ill].gen_matched && jj[ijj].gen_matched;
            mylljj.gen_p4 = mylljj.gen_matched ? ll[ill].gen_p4 + jj[ijj].gen_p4 : null_p4;
            mylljj.gen_DR = mylljj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mylljj.p4, mylljj.gen_p4) : -1.;
            mylljj.gen_DPhi = mylljj.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(mylljj.p4, mylljj.gen_p4)) : -1.;
            mylljj.gen_DPtOverPt = mylljj.gen_matched ? (mylljj.p4.Pt() - mylljj.gen_p4.Pt()) / mylljj.p4.Pt() : -10.;
            mylljj.gen_lep1_p4 = leptons[ilep1].gen_p4;
            mylljj.gen_lep2_p4 = leptons[ilep2].gen_p4;
            mylljj.gen_jet1_p4 = jets[ijet1].gen_p4;
            mylljj.gen_jet2_p4 = jets[ijet2].gen_p4;
            mylljj.gen_ll_p4 = ll[ill].gen_p4;
            mylljj.gen_jj_p4 = jj[ijj].gen_p4;
            // blind copy of the jj content
            mylljj.ijet1 = jj[ijj].ijet1;
            mylljj.ijet2 = jj[ijj].ijet2;
            mylljj.btag_MM = jj[ijj].btag_MM;
            mylljj.sumCSV = jj[ijj].sumCSV;
            mylljj.sumCMVAv2 = jj[ijj].sumCMVAv2;
            mylljj.DR_j_j = jj[ijj].DR_j_j;
            mylljj.DPhi_j_j = jj[ijj].DPhi_j_j;
            mylljj.ht_j_j = jj[ijj].ht_j_j;
            mylljj.gen_matched_bbPartons = jj[ijj].gen_matched_bbPartons;
            mylljj.gen_matched_bbHadrons = jj[ijj].gen_matched_bbHadrons;
            mylljj.gen_bb = jj[ijj].gen_bb;
            mylljj.gen_bc = jj[ijj].gen_bc;
            mylljj.gen_bl = jj[ijj].gen_bl;
            mylljj.gen_cc = jj[ijj].gen_cc;
            mylljj.gen_cl = jj[ijj].gen_cl;
            mylljj.gen_ll = jj[ijj].gen_ll;
            // blind copy of the ll content
            mylljj.ilep1 = ll[ill].ilep1;
            mylljj.ilep2 = ll[ill].ilep2;
            mylljj.isOS = ll[ill].isOS;
            mylljj.isPlusMinus = ll[ill].isPlusMinus;
            mylljj.isMinusPlus = ll[ill].isMinusPlus;
            mylljj.isMuMu = ll[ill].isMuMu;
            mylljj.isElEl = ll[ill].isElEl;
            mylljj.isElMu = ll[ill].isElMu;
            mylljj.isMuEl = ll[ill].isMuEl;
            mylljj.isSF = ll[ill].isSF;
            mylljj.DR_l_l = ll[ill].DR_l_l;
            mylljj.DPhi_l_l = ll[ill].DPhi_l_l;
            mylljj.ht_l_l = ll[ill].ht_l_l;
            mylljj.trigger_efficiency = ll[ill].trigger_efficiency;
            mylljj.trigger_efficiency_downVariated = ll[ill].trigger_efficiency_downVariated;
            mylljj.trigger_efficiency_upVariated = ll[ill].trigger_efficiency_upVariated;
        
            float DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2;
            DR_j1l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep1].p4);
            DR_j1l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep2].p4);
            DR_j2l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep1].p4);
            DR_j2l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep2].p4);
            mylljj.maxDR_l_j = std::max({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            mylljj.minDR_l_j = std::min({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            mylljj.DR_ll_jj = ROOT::Math::VectorUtil::DeltaR(ll[ill].p4, jj[ijj].p4);
            mylljj.DPhi_ll_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, jj[ijj].p4));
            mylljj.MT_fullsystem = mylljj.p4.Mt();
            mylljj.melaAngles = getMELAAngles(ll[ill].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4);
            mylljj.visMelaAngles = getMELAAngles(ll[ill].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4); // only take the visible part of the H(ww) candidate - WHAT??????????????
            
            // Compute MT2. See https://arxiv.org/pdf/1309.6318v1.pdf and https://arxiv.org/pdf/1411.4312v5.pdf
            // No need to compute MT in this case, right? It uses invisible momenta that in this case we don't have


            //Counters
            tmp_count_has2leptons_1lljj = event_weight;
            if (mylljj.isElEl)
                tmp_count_has2leptons_elel_1lljj = event_weight;
            if (mylljj.isElMu)
                tmp_count_has2leptons_elmu_1lljj = event_weight;
            if (mylljj.isMuEl)
                tmp_count_has2leptons_muel_1lljj = event_weight;
            if (mylljj.isMuMu)
                tmp_count_has2leptons_mumu_1lljj = event_weight;
            if (mylljj.btag_MM)
            {
                tmp_count_has2leptons_1lljj_2btagM = event_weight;
                if (mylljj.isElEl)
                    tmp_count_has2leptons_elel_1lljj_2btagM = event_weight;
                if (mylljj.isElMu)
                    tmp_count_has2leptons_elmu_1lljj_2btagM = event_weight;
                if (mylljj.isMuEl)
                    tmp_count_has2leptons_muel_1lljj_2btagM = event_weight;
                if (mylljj.isMuMu)
                    tmp_count_has2leptons_mumu_1lljj_2btagM = event_weight;
            }
            // Fill
            lljj.push_back(mylljj);
        } //end loop over jj
    } //end loop over ll

    std::sort(lljj.begin(), lljj.end(), [&](HtoZA::DileptonDijet& a, const HtoZA::DileptonDijet& b){ return a.sumCMVAv2 > b.sumCMVAv2; });

    // Keep only the first candidate
    if (lljj.size() > 1){
        lljj.resize(1);
    }


    /*
    TH1F *mass_lljj = new TH1F("h1", "m_lljj", 1500, 0, 2000);

    if(lljj.size() != 0) {
        mass_lljj->Fill(lljj.at(0).p4.M());
        std::cout << "lljj.at(0).p4.M(): " << lljj.at(0).p4.M() << std::endl;
    } else { std::cout << "lljj is empty!" << std::endl;
      }

    
    TCanvas c1("canvas", "canvas", 1024, 768);
    c1.cd();
    mass_lljj->Draw();
    c1.SaveAs("m_lljj.pdf");
    */

    // ***** ***** *****
    // Event variables
    // ***** ***** *****

    // HT: the two selected leptons - if present - plus all selected jets
    HT = 0;
    if (lljj.size() > 0)
        // take the first lljj since it's been resized: it has only one entry! :-)
        HT += lljj[0].lep1_p4.Pt() + lljj[0].lep2_p4.Pt();
    for (unsigned int ijet=0; ijet < jets.size(); ijet++) {
        HT += jets[ijet].p4.Pt();
    }

    nJetsL = jets.size();  //We selected the jets to be all loose
    if (! doingSystematics()) {
        nBJetsM = 0;
        for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
            if (jets[ijet].btag_M)
                nBJetsM++;
        }


        nMuonsT = 0; // We selected the muons as tight
        nElectronsM = 0;  // We selected the electrons as medium
        for (unsigned int ilepton = 0; ilepton < leptons.size(); ilepton++)
        {
            if (leptons[ilepton].isMu) {
                nMuonsT++;
            }

            if (leptons[ilepton].isEl) {
                nElectronsM++;
            }
        }

        count_has2leptons += tmp_count_has2leptons;
        count_has2leptons_elel += tmp_count_has2leptons_elel;
        count_has2leptons_elmu += tmp_count_has2leptons_elmu;
        count_has2leptons_muel += tmp_count_has2leptons_muel;
        count_has2leptons_mumu += tmp_count_has2leptons_mumu;
        count_has2leptons_1lljj += tmp_count_has2leptons_1lljj;
        count_has2leptons_elel_1lljj += tmp_count_has2leptons_elel_1lljj;
        count_has2leptons_elmu_1lljj += tmp_count_has2leptons_elmu_1lljj;
        count_has2leptons_muel_1lljj += tmp_count_has2leptons_muel_1lljj;
        count_has2leptons_mumu_1lljj += tmp_count_has2leptons_mumu_1lljj;
        count_has2leptons_1lljj_2btagM += tmp_count_has2leptons_1lljj_2btagM;
        count_has2leptons_elel_1lljj_2btagM += tmp_count_has2leptons_elel_1lljj_2btagM;
        count_has2leptons_elmu_1lljj_2btagM += tmp_count_has2leptons_elmu_1lljj_2btagM;
        count_has2leptons_muel_1lljj_2btagM += tmp_count_has2leptons_muel_1lljj_2btagM;
        count_has2leptons_mumu_1lljj_2btagM += tmp_count_has2leptons_mumu_1lljj_2btagM;
    } //end if (! doingSystematics())

}

void HtoZAAnalyzer::endJob(MetadataManager& metadata) {
    if (! doingSystematics()) {
        metadata.add(this->m_name + "_count_has2leptons", count_has2leptons);
        metadata.add(this->m_name + "_count_has2leptons_elel", count_has2leptons_elel);
        metadata.add(this->m_name + "_count_has2leptons_elmu", count_has2leptons_elmu);
        metadata.add(this->m_name + "_count_has2leptons_muel", count_has2leptons_muel);
        metadata.add(this->m_name + "_count_has2leptons_mumu", count_has2leptons_mumu);
        metadata.add(this->m_name + "_count_has2leptons_1lljj", count_has2leptons_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_elel_1lljj", count_has2leptons_elel_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_elmu_1lljj", count_has2leptons_elmu_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_muel_1lljj", count_has2leptons_muel_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_mumu_1lljj", count_has2leptons_mumu_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_1lljj_2btagM", count_has2leptons_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_elel_1lljj_2btagM", count_has2leptons_elel_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_elmu_1lljj_2btagM", count_has2leptons_elmu_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_muel_1lljj_2btagM", count_has2leptons_muel_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_mumu_1lljj_2btagM", count_has2leptons_mumu_1lljj_2btagM);
    }
}





