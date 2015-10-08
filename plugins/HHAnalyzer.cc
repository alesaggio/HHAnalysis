#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Categories.h>

#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#define HHANADEBUG 0

using namespace HH;

void HHAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", config);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", config);
    manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, subleading as muon", config);
    manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, subleading as electron", config);
}


void HHAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const CategoryManager&) {

    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& allmuons = producers.get<MuonsProducer>("muons");
    const METProducer& met = producers.get<METProducer>("met");
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    //const METProducer& met = producers.get<METProducer>("puppimet"); 
    float mh = event.isRealData() ? 125.02 : 125.0;

    // ********** 
    // Leptons and dileptons
    // ********** 

    leptons.clear();
    ll.clear();

    for (unsigned int ielectron = 0; ielectron < allelectrons.p4.size(); ielectron++)
    {
        if (allelectrons.ids[ielectron][m_electron_tight_wp_name]
            && allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut
            && allelectrons.p4[ielectron].Pt() > m_electronPtCut
            && abs(allelectrons.p4[ielectron].Eta()) < m_electronEtaCut) 
            {
            electrons.push_back(ielectron);
            Lepton ele = { allelectrons.p4[ielectron], allelectrons.charge[ielectron], ielectron, false, true };
            leptons.push_back(ele);
        }
    }//end of loop on electrons

    for (unsigned int imuon = 0; imuon < allmuons.p4.size(); imuon++)
    {
        if (allmuons.isTight[imuon] 
            && allmuons.relativeIsoR04_withEA[imuon] < m_muonIsoCut 
            && allmuons.p4[imuon].Pt() > m_muonPtCut 
            && abs(allmuons.p4[imuon].Eta()) < m_muonEtaCut)
            {
            muons.push_back(imuon);
            Lepton mu = { allmuons.p4[imuon], allmuons.charge[imuon], imuon, true, false };
            leptons.push_back(mu);
        }
    }//end of loop on muons
           
    std::sort(leptons.begin(), leptons.end(), [](const Lepton& lep1, const Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt(); });     

    for (unsigned int ilep1 = 0; ilep1 < leptons.size(); ilep1++)
    {
        leptons_p4.push_back(leptons[ilep1].p4);
        leptons_idx.push_back(leptons[ilep1].idx);
        leptons_isMu.push_back(leptons[ilep1].isMu);
        leptons_isEl.push_back(leptons[ilep1].isEl);
        for (unsigned int ilep2 = ilep1+1; ilep2 < leptons.size(); ilep2++)
        {
            Dilepton dilep = { leptons[ilep1].p4 + leptons[ilep2].p4, std::make_pair(ilep1, ilep2), leptons[ilep1].charge * leptons[ilep2].charge, (leptons[ilep1].isMu && leptons[ilep2].isMu), (leptons[ilep1].isEl && leptons[ilep2].isEl), (leptons[ilep1].isEl && leptons[ilep2].isMu), (leptons[ilep1].isMu && leptons[ilep2].isEl) };
            ll.push_back(dilep); 
            ll_p4.push_back(dilep.p4);
            llmet_p4.push_back(dilep.p4 + met.p4);
            ll_idx.push_back(dilep.idxs);
            ll_isMuMu.push_back(dilep.isMuMu);
            ll_isElEl.push_back(dilep.isElEl);
            ll_isElMu.push_back(dilep.isElMu);
            ll_isMuEl.push_back(dilep.isMuEl);
            ll_DR.push_back(ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4));
            ll_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            float dphi = ROOT::Math::VectorUtil::DeltaPhi(dilep.p4, met.p4);
            ll_DPhi_met.push_back(dphi);
            ll_MT.push_back((dilep.p4 + met.p4).M());
            ll_MT_formula.push_back(std::sqrt(2 * dilep.p4.Pt() * met.p4.Pt() * (1-std::cos(dphi))));
            float mindphi = std::min(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, met.p4), ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep2].p4, met.p4));
            ll_minDPhi_lmet.push_back(mindphi);
            float maxdphi = std::max(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, met.p4), ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep2].p4, met.p4));
            ll_maxDPhi_lmet.push_back(maxdphi);
            ll_projectedMet.push_back(mindphi >= M_PI ? met.p4.Pt() : met.p4.Pt() * std::sin(mindphi));
        }
    }

    // ***** 
    // Jets and dijets 
    // ***** 

    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        if ((alljets.p4[ijet].Pt() > m_jetPtCut) 
            && (abs(alljets.p4[ijet].Eta()) < m_jetEtaCut))
            {
            jets_p4.push_back(alljets.p4[ijet]);
            jets_idx.push_back(ijet);
            if (alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName) > m_jet_bDiscrCut){
                bjets_p4.push_back(alljets.p4[ijet]);
                bjets_idx.push_back(ijet);
            }
        }
    }

    float diffWithMh = 14000;
    unsigned int dijetCounter = 0;
    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ijet1 = 0; ijet1 < jets_p4.size(); ijet1++)
    {  
        for (unsigned int ijet2 = ijet1+1; ijet2 < jets_p4.size(); ijet2++)
        {
            LorentzVector jj = jets_p4[ijet1] + jets_p4[ijet2];
            jj_p4.push_back(jj);
            jj_idx.push_back(std::make_pair(ijet1, ijet2));
            jj_DR.push_back(ROOT::Math::VectorUtil::DeltaR(jets_p4[ijet1], jets_p4[ijet2]));
            jj_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(jets_p4[ijet1], jets_p4[ijet2]));
            jj_DPhi_met.push_back(ROOT::Math::VectorUtil::DeltaPhi(jj, met.p4));
            jj_minDPhi_jmet.push_back(std::min(ROOT::Math::VectorUtil::DeltaPhi(jets_p4[ijet1], met.p4), ROOT::Math::VectorUtil::DeltaPhi(jets_p4[ijet2], met.p4)));
            jj_maxDPhi_jmet.push_back(std::max(ROOT::Math::VectorUtil::DeltaPhi(jets_p4[ijet1], met.p4), ROOT::Math::VectorUtil::DeltaPhi(jets_p4[ijet2], met.p4)));
            if (abs(jj.M() - mh) < diffWithMh) {
                h_dijet_idx = dijetCounter;
                diffWithMh = abs(jj.M() - mh);
            }
            dijetCounter++;
        }
    }

    diffWithMh = 14000;
    dijetCounter = 0;
    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ibjet1 = 0; ibjet1 < bjets_p4.size(); ibjet1++)
    {  
        for (unsigned int ibjet2 = ibjet1+1; ibjet2 < bjets_p4.size(); ibjet2++)
        {
            LorentzVector bb = bjets_p4[ibjet1] + bjets_p4[ibjet2];
            bb_p4.push_back(bb);
            bb_idx.push_back(std::make_pair(ibjet1, ibjet2));
            bb_DR.push_back(ROOT::Math::VectorUtil::DeltaR(bjets_p4[ibjet1], bjets_p4[ibjet2]));
            bb_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(bjets_p4[ibjet1], bjets_p4[ibjet2]));
            bb_DPhi_met.push_back(ROOT::Math::VectorUtil::DeltaPhi(bb, met.p4));
            bb_minDPhi_jmet.push_back(std::min(ROOT::Math::VectorUtil::DeltaPhi(bjets_p4[ibjet1], met.p4), ROOT::Math::VectorUtil::DeltaPhi(bjets_p4[ibjet2], met.p4)));
            bb_maxDPhi_jmet.push_back(std::max(ROOT::Math::VectorUtil::DeltaPhi(bjets_p4[ibjet1], met.p4), ROOT::Math::VectorUtil::DeltaPhi(bjets_p4[ibjet2], met.p4)));
            if (abs(bb.M() - mh) < diffWithMh) {
                h_dibjet_idx = dijetCounter;
                diffWithMh = abs(bb.M() - mh);
            }
            dijetCounter++;
        }
    }
            
    // ********** 
    // lljj, llbb, +met
    // ********** 
    for (unsigned int ill = 0; ill < ll_p4.size(); ill++)
    {
        for (unsigned int ijj = 0; ijj < jj_p4.size(); ijj++)
        {
            LorentzVector lljj = ll_p4[ill] + jj_p4[ijj];
            lljj_p4.push_back(lljj);
            lljj_idx.push_back(std::make_pair(ill, ijj));
            lljj_DR.push_back(ROOT::Math::VectorUtil::DeltaR(ll_p4[ill], jj_p4[ijj]));
            lljj_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(ll_p4[ill], jj_p4[ijj]));
            float DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2;
            DR_j1l1 = ROOT::Math::VectorUtil::DeltaR(jets_p4[jj_idx[ijj].first], leptons_p4[ll_idx[ill].first]);
            DR_j1l2 = ROOT::Math::VectorUtil::DeltaR(jets_p4[jj_idx[ijj].first], leptons_p4[ll_idx[ill].second]);
            DR_j2l1 = ROOT::Math::VectorUtil::DeltaR(jets_p4[jj_idx[ijj].second], leptons_p4[ll_idx[ill].first]);
            DR_j2l2 = ROOT::Math::VectorUtil::DeltaR(jets_p4[jj_idx[ijj].second], leptons_p4[ll_idx[ill].second]);
            float maxDR = std::max({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            lljj_maxDR_lj.push_back(maxDR);
            float minDR = std::min({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            lljj_minDR_lj.push_back(minDR);
            lljjmet_p4.push_back(lljj + met.p4);
            lljjmet_DR.push_back(ROOT::Math::VectorUtil::DeltaR(lljj, met.p4));
            lljjmet_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(ll_p4[ill] + met.p4, jj_p4[ijj]));
            lljjmet_cosThetaStar_CS.push_back(getCosThetaStar_CS(ll_p4[ill] + met.p4, jj_p4[ijj]));
        }
    }

    for (unsigned int ill = 0; ill < ll_p4.size(); ill++)
    {
        for (unsigned int ibb = 0; ibb < bb_p4.size(); ibb++)
        {
            LorentzVector llbb = ll_p4[ill] + bb_p4[ibb];
            llbb_p4.push_back(llbb);
            llbb_idx.push_back(std::make_pair(ill, ibb));
            llbb_DR.push_back(ROOT::Math::VectorUtil::DeltaR(ll_p4[ill], jj_p4[ibb]));
            llbb_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(ll_p4[ill], jj_p4[ibb]));
            float DR_b1l1, DR_b1l2, DR_b2l1, DR_b2l2;
            DR_b1l1 = ROOT::Math::VectorUtil::DeltaR(bjets_p4[bb_idx[ibb].first], leptons_p4[ll_idx[ill].first]);
            DR_b1l2 = ROOT::Math::VectorUtil::DeltaR(bjets_p4[bb_idx[ibb].first], leptons_p4[ll_idx[ill].second]);
            DR_b2l1 = ROOT::Math::VectorUtil::DeltaR(bjets_p4[bb_idx[ibb].second], leptons_p4[ll_idx[ill].first]);
            DR_b2l2 = ROOT::Math::VectorUtil::DeltaR(bjets_p4[bb_idx[ibb].second], leptons_p4[ll_idx[ill].second]);
            float maxDR = std::max({DR_b1l1, DR_b1l2, DR_b2l1, DR_b2l2});
            llbb_maxDR_lb.push_back(maxDR);
            float minDR = std::min({DR_b1l1, DR_b1l2, DR_b2l1, DR_b2l2});
            llbb_minDR_lb.push_back(minDR);
            llbbmet_p4.push_back(llbb + met.p4);
            llbbmet_DR.push_back(ROOT::Math::VectorUtil::DeltaR(llbb, met.p4));
            llbbmet_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(llbb, met.p4));
            llbbmet_DPhi.push_back(ROOT::Math::VectorUtil::DeltaPhi(ll_p4[ill] + met.p4, bb_p4[ibb]));
            llbbmet_cosThetaStar_CS.push_back(getCosThetaStar_CS(ll_p4[ill] + met.p4, bb_p4[ibb]));
        }
    }

    // ***** ***** *****
    // Event variables
    // ***** ***** *****
    nJets = jets_p4.size();
    nBJets = bjets_p4.size();
    nMuons = muons.size();
    nElectrons = electrons.size();
    nLeptons = leptons.size();



    if (!event.isRealData())
    {
// ***** ***** *****
// Get the MC truth information on the hard process
// ***** ***** *****
// from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
//    enum StatusBits {
//0      kIsPrompt = 0,
//1      kIsDecayedLeptonHadron,
//2      kIsTauDecayProduct,
//3      kIsPromptTauDecayProduct,
//4      kIsDirectTauDecayProduct,
//5      kIsDirectPromptTauDecayProduct,
//6      kIsDirectHadronDecayProduct,
//7      kIsHardProcess,
//8      kFromHardProcess,
//9      kIsHardProcessTauDecayProduct,
//10      kIsDirectHardProcessTauDecayProduct,
//11      kFromHardProcessBeforeFSR,
//12      kIsFirstCopy,
//13      kIsLastCopy,
//14      kIsLastCopyBeforeFSR
//    };

        const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");
    
        BRANCH(gen_B1, LorentzVector);
        BRANCH(gen_B2, LorentzVector);
        BRANCH(gen_B1FSR, LorentzVector);
        BRANCH(gen_B2FSR, LorentzVector);
        BRANCH(gen_Nu1, LorentzVector);
        BRANCH(gen_Nu2, LorentzVector);
        BRANCH(gen_L1, LorentzVector);
        BRANCH(gen_L2, LorentzVector);
        BRANCH(gen_LL, LorentzVector);
        BRANCH(gen_BB, LorentzVector);
        BRANCH(gen_BBFSR, LorentzVector);
        BRANCH(gen_L1FSRNu, LorentzVector);
        BRANCH(gen_L2FSRNu, LorentzVector);
        BRANCH(gen_L1FSR, LorentzVector);
        BRANCH(gen_L2FSR, LorentzVector);
        BRANCH(gen_LLFSR, LorentzVector);
        BRANCH(gen_NuNu, LorentzVector);
        BRANCH(gen_LLNuNu, LorentzVector);
        BRANCH(gen_LLFSRNuNu, LorentzVector);
        BRANCH(gen_LLFSRNuNuBB, LorentzVector);
        gen_B1.SetPxPyPzE(0.,0.,0.,0.);
        gen_B2.SetPxPyPzE(0.,0.,0.,0.);
        gen_B1FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_B2FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_Nu1.SetPxPyPzE(0.,0.,0.,0.);
        gen_Nu2.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2.SetPxPyPzE(0.,0.,0.,0.);
        gen_LL.SetPxPyPzE(0.,0.,0.,0.);
        gen_BB.SetPxPyPzE(0.,0.,0.,0.);
        gen_BBFSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1FSRNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2FSRNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_NuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLNuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSRNuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSRNuNuBB.SetPxPyPzE(0.,0.,0.,0.);
    
        BRANCH(gen_iX, int);
        BRANCH(gen_iH1, int);
        BRANCH(gen_iH2, int);
        BRANCH(gen_iV1, int);
        BRANCH(gen_iV2, int);
        BRANCH(gen_iB1, int);
        BRANCH(gen_iB2, int);
        BRANCH(gen_iL1, int);
        BRANCH(gen_iL2, int);
        BRANCH(gen_iNu1, int);
        BRANCH(gen_iNu2, int);
        BRANCH(gen_iG1, std::vector<int>);
        BRANCH(gen_iG2, std::vector<int>);
        BRANCH(gen_iGlu1, std::vector<int>);
        BRANCH(gen_iGlu2, std::vector<int>);
        gen_iX = -1;
        gen_iH1 = gen_iH2 = -1;
        gen_iV1 = gen_iV2 = -1;
        gen_iB1 = gen_iB2 = -1;
        gen_iL1 = gen_iL2 = -1;
        gen_iNu1 = gen_iNu2 = -1;
        gen_iG1.clear();
        gen_iG2.clear();
        gen_iGlu1.clear();
        gen_iGlu2.clear();
        bool isThereFSRforL1 = false;
        bool isThereFSRforL2 = false;
        bool isThereFSRforB1 = false;
        bool isThereFSRforB2 = false;
    
        // Construct the signal gen-level info
        for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp.pruned_status_flags[ip]);
            if (!flags.test(13)) continue; // take the last copies
            if (flags.test(8)) // first look at the hard process
            {
                if (HHANADEBUG)
                    std::cout << "\tip= " << ip << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[ip] << "\tflags= " << flags << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[ip].Pt() << ", " << gp.pruned_p4[ip].Eta() << ", " << gp.pruned_p4[ip].Phi() << ", " << gp.pruned_p4[ip].E() << ", " << gp.pruned_p4[ip].M() << ")" << std::endl;
                if (abs(gp.pruned_pdg_id[ip]) == 25)
                {
                    if (gen_iH1 == -1)
                        gen_iH1 = ip;
                    else if (gen_iH2 == -1)
                        gen_iH2 = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 5)
                {
                    if (gen_iB1 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforB1 = true;
                        gen_iB1 = ip;
                    }
                    else if (gen_iB2 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforB2 = true;
                        gen_iB2 = ip;
                    }
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 11 || abs(gp.pruned_pdg_id[ip]) == 13)
                {
                // if the lepton is not in the hard process itself, then look for FSR photons later on
                    if (gen_iL1 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforL1 = true;
                        gen_iL1 = ip;
                    }
                    else if (gen_iL2 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforL2 = true;
                        gen_iL2 = ip;
                    }
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 12 || abs(gp.pruned_pdg_id[ip]) == 14 || abs(gp.pruned_pdg_id[ip]) == 16)
                {
                    if (gen_iNu1 == -1)
                        gen_iNu1 = ip;
                    else if (gen_iNu2 == -1)
                        gen_iNu2 = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 35 || abs(gp.pruned_pdg_id[ip]) == 39)
                {
                    if (gen_iX == -1)
                        gen_iX = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 23 || abs(gp.pruned_pdg_id[ip]) == 24)
                {
                    if (gen_iV1 == -1)
                        gen_iV1 = ip;
                    else if (gen_iV2 == -1)
                        gen_iV2 = ip;
                }
            } // end if coming from hard process
        } // end of loop over pruned gen particles
        // new loop to find FSR photons
        if (isThereFSRforL1 || isThereFSRforL2)
        {
            for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
                if (gp.pruned_pdg_id[ip] != 22) continue;
                std::bitset<15> flags (gp.pruned_status_flags[ip]);
                if (!flags.test(13)) continue; // take the last copies
                for (unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size(); imother++)
                {
                    if (isThereFSRforL1)
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iL1].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iL1][jmother])
                                gen_iG1.push_back(ip);
                        }
                    }
                    if (isThereFSRforL2) 
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iL2].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iL2][jmother])
                                gen_iG2.push_back(ip);
                        }
                    }
                }
            } // end of loop over pruned gen particles
        } // end of FSR loop
        // new loop to find FSR gluons
        if (isThereFSRforB1 || isThereFSRforB2)
        {
            for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
                if (gp.pruned_pdg_id[ip] != 21) continue;
                std::bitset<15> flags (gp.pruned_status_flags[ip]);
                if (!flags.test(13)) continue; // take the last copies
                for (unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size(); imother++)
                {
                    if (isThereFSRforB1)
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iB1].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iB1][jmother])
                                gen_iGlu1.push_back(ip);
                        }
                    }
                    if (isThereFSRforB2) 
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iB2].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iB2][jmother])
                                gen_iGlu2.push_back(ip);
                        }
                    }
                }
            } // end of loop over pruned gen particles
        } // end of FSR loop
    
        if (HHANADEBUG)
        {
            if (gen_iX != -1)
                std::cout << "\tiX= " << gen_iX << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iX] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iX]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iX].Pt() << ", " << gp.pruned_p4[gen_iX].Eta() << ", " << gp.pruned_p4[gen_iX].Phi() << ", " << gp.pruned_p4[gen_iX].E() << ", " << gp.pruned_p4[gen_iX].M() << ")" << std::endl;
            std::cout << "\tiH1= " << gen_iH1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iH1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iH1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iH1].Pt() << ", " << gp.pruned_p4[gen_iH1].Eta() << ", " << gp.pruned_p4[gen_iH1].Phi() << ", " << gp.pruned_p4[gen_iH1].E() << ", " << gp.pruned_p4[gen_iH1].M() << ")" << std::endl;
            std::cout << "\tiH2= " << gen_iH2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iH2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iH2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iH2].Pt() << ", " << gp.pruned_p4[gen_iH2].Eta() << ", " << gp.pruned_p4[gen_iH2].Phi() << ", " << gp.pruned_p4[gen_iH2].E() << ", " << gp.pruned_p4[gen_iH2].M() << ")" << std::endl;
            std::cout << "\tiB1= " << gen_iB1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iB1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iB1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iB1].Pt() << ", " << gp.pruned_p4[gen_iB1].Eta() << ", " << gp.pruned_p4[gen_iB1].Phi() << ", " << gp.pruned_p4[gen_iB1].E() << ", " << gp.pruned_p4[gen_iB1].M() << ")" << std::endl;
            std::cout << "\tiB2= " << gen_iB2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iB2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iB2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iB2].Pt() << ", " << gp.pruned_p4[gen_iB2].Eta() << ", " << gp.pruned_p4[gen_iB2].Phi() << ", " << gp.pruned_p4[gen_iB2].E() << ", " << gp.pruned_p4[gen_iB2].M() << ")" << std::endl;
            if (gen_iV1 != -1)
                std::cout << "\tiV1= " << gen_iV1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iV1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iV1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iV1].Pt() << ", " << gp.pruned_p4[gen_iV1].Eta() << ", " << gp.pruned_p4[gen_iV1].Phi() << ", " << gp.pruned_p4[gen_iV1].E() << ", " << gp.pruned_p4[gen_iV1].M() << ")" << std::endl;
            if (gen_iV2 != -1)
                std::cout << "\tiV2= " << gen_iV2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iV2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iV2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iV2].Pt() << ", " << gp.pruned_p4[gen_iV2].Eta() << ", " << gp.pruned_p4[gen_iV2].Phi() << ", " << gp.pruned_p4[gen_iV2].E() << ", " << gp.pruned_p4[gen_iV2].M() << ")" << std::endl;
            std::cout << "\tiL1= " << gen_iL1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iL1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iL1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iL1].Pt() << ", " << gp.pruned_p4[gen_iL1].Eta() << ", " << gp.pruned_p4[gen_iL1].Phi() << ", " << gp.pruned_p4[gen_iL1].E() << ", " << gp.pruned_p4[gen_iL1].M() << ")" << std::endl;
            std::cout << "\tiL2= " << gen_iL2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iL2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iL2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iL2].Pt() << ", " << gp.pruned_p4[gen_iL2].Eta() << ", " << gp.pruned_p4[gen_iL2].Phi() << ", " << gp.pruned_p4[gen_iL2].E() << ", " << gp.pruned_p4[gen_iL2].M() << ")" << std::endl;
            std::cout << "\tiNu1= " << gen_iNu1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iNu1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iNu1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iNu1].Pt() << ", " << gp.pruned_p4[gen_iNu1].Eta() << ", " << gp.pruned_p4[gen_iNu1].Phi() << ", " << gp.pruned_p4[gen_iNu1].E() << ", " << gp.pruned_p4[gen_iNu1].M() << ")" << std::endl;
            std::cout << "\tiNu2= " << gen_iNu2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iNu2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iNu2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iNu2].Pt() << ", " << gp.pruned_p4[gen_iNu2].Eta() << ", " << gp.pruned_p4[gen_iNu2].Phi() << ", " << gp.pruned_p4[gen_iNu2].E() << ", " << gp.pruned_p4[gen_iNu2].M() << ")" << std::endl;
        
            for (unsigned int iG = 0; iG < gen_iG1.size(); iG++)
            {
                std::cout << "\tiG1[" << iG << "]= " << gen_iG1[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iG1[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iG1[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iG1[iG]].Pt() << ", " << gp.pruned_p4[gen_iG1[iG]].Eta() << ", " << gp.pruned_p4[gen_iG1[iG]].Phi() << ", " << gp.pruned_p4[gen_iG1[iG]].E() << ", " << gp.pruned_p4[gen_iG1[iG]].M() << ")" << std::endl;
            }
            for (unsigned int iG = 0; iG < gen_iG2.size(); iG++)
            {
                std::cout << "\tiG2[" << iG << "]= " << gen_iG2[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iG2[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iG2[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iG2[iG]].Pt() << ", " << gp.pruned_p4[gen_iG2[iG]].Eta() << ", " << gp.pruned_p4[gen_iG2[iG]].Phi() << ", " << gp.pruned_p4[gen_iG2[iG]].E() << ", " << gp.pruned_p4[gen_iG2[iG]].M() << ")" << std::endl;
            }
            for (unsigned int iGlu = 0; iGlu < gen_iGlu1.size(); iGlu++)
            {
                std::cout << "\tiGlu1[" << iGlu << "]= " << gen_iGlu1[iGlu] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iGlu1[iGlu]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iGlu1[iGlu]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iGlu1[iGlu]].Pt() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].Eta() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].Phi() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].E() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].M() << ")" << std::endl;
            }
            for (unsigned int iGlu = 0; iGlu < gen_iGlu2.size(); iGlu++)
            {
                std::cout << "\tiGlu2[" << iGlu << "]= " << gen_iGlu2[iGlu] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iGlu2[iGlu]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iGlu2[iGlu]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iGlu2[iGlu]].Pt() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].Eta() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].Phi() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].E() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].M() << ")" << std::endl;
            }
        } // end if HHANADEBUG
    
    
        gen_B1 += gp.pruned_p4[gen_iB1];
        gen_B2 += gp.pruned_p4[gen_iB2];
        gen_Nu1 += gp.pruned_p4[gen_iNu1];
        gen_Nu2 += gp.pruned_p4[gen_iNu2];
        gen_L1 += gp.pruned_p4[gen_iL1];
        gen_L2 += gp.pruned_p4[gen_iL2];
        gen_LL += gen_L1 + gen_L2;
        gen_NuNu += gen_Nu1 + gen_Nu2;
        gen_B1FSR += gen_B1;
        gen_B2FSR += gen_B2;
        gen_L1FSR += gen_L1;
        gen_L2FSR += gen_L2;
        for (unsigned int iG = 0; iG < gen_iG1.size(); iG++)
            gen_L1FSR += gp.pruned_p4[gen_iG1[iG]];
        for (unsigned int iG = 0; iG < gen_iG2.size(); iG++)
            gen_L2FSR += gp.pruned_p4[gen_iG2[iG]];
        for (unsigned int iGlu = 0; iGlu < gen_iGlu1.size(); iGlu++)
            gen_B1FSR += gp.pruned_p4[gen_iGlu1[iGlu]];
        for (unsigned int iGlu = 0; iGlu < gen_iGlu2.size(); iGlu++)
            gen_B2FSR += gp.pruned_p4[gen_iGlu2[iGlu]];
        gen_BB += gen_B1 + gen_B2;
        gen_BBFSR += gen_B1FSR + gen_B2FSR;
        gen_L1FSRNu += gen_L1FSR + gen_Nu1;
        gen_L2FSRNu += gen_L2FSR + gen_Nu2;
        gen_LLFSR += gen_L1FSR + gen_L2FSR;
        gen_LLNuNu += gen_LL + gen_NuNu;
        gen_LLFSRNuNu += gen_LLFSR + gen_NuNu;
        gen_LLFSRNuNuBB += gen_LLFSRNuNu + gen_BBFSR;
    
        if (HHANADEBUG)
        {
            std::cout << "\tgen_B1.M()= " << gen_B1.M() << std::endl;
            std::cout << "\tgen_B2.M()= " << gen_B2.M() << std::endl;
            std::cout << "\tgen_B1FSR.M()= " << gen_B1FSR.M() << std::endl;
            std::cout << "\tgen_B2FSR.M()= " << gen_B2FSR.M() << std::endl;
            std::cout << "\tgen_Nu1.M()= " << gen_Nu1.M() << std::endl;
            std::cout << "\tgen_Nu2.M()= " << gen_Nu2.M() << std::endl;
            std::cout << "\tgen_L1.M()= " << gen_L1.M() << std::endl;
            std::cout << "\tgen_L2.M()= " << gen_L2.M() << std::endl;
            std::cout << "\tgen_LL.M()= " << gen_LL.M() << std::endl;
            std::cout << "\tgen_L1FSRNu.M()= " << gen_L1FSRNu.M() << std::endl;
            std::cout << "\tgen_L2FSRNu.M()= " << gen_L2FSRNu.M() << std::endl;
            std::cout << "\tgen_L1FSR.M()= " << gen_L1FSR.M() << std::endl;
            std::cout << "\tgen_L2FSR.M()= " << gen_L2FSR.M() << std::endl;
            std::cout << "\tgen_LLFSR.M()= " << gen_LLFSR.M() << std::endl;
            std::cout << "\tgen_NuNu.M()= " << gen_NuNu.M() << std::endl;
            std::cout << "\tgen_LLNuNu.M()= " << gen_LLNuNu.M() << std::endl;
            std::cout << "\tgen_LLFSRNuNu.M()= " << gen_LLFSRNuNu.M() << std::endl;
            std::cout << "\tgen_BB.M()= " << gen_BB.M() << std::endl;
            std::cout << "\tgen_BBFSR.M()= " << gen_BBFSR.M() << std::endl;
            std::cout << "\tgen_LLFSRNuNuBB.M()= " << gen_LLFSRNuNuBB.M() << std::endl;
        }

        // ***** ***** *****
        // Matching
        // ***** ***** *****
        BRANCH(gen_deltaR_jet_B1, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B2, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B2FSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2FSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2FSR, std::vector<float>);    
    
        for (auto p4: alljets.gen_p4) {
            gen_deltaR_jet_B1.push_back(deltaR(p4, gen_B1));
            gen_deltaR_jet_B2.push_back(deltaR(p4, gen_B2));
            gen_deltaR_jet_B1FSR.push_back(deltaR(p4, gen_B1FSR));
            gen_deltaR_jet_B2FSR.push_back(deltaR(p4, gen_B2FSR));
        }
        for (auto p4: allelectrons.gen_p4) {
            gen_deltaR_electron_L1.push_back(deltaR(p4, gen_L1));
            gen_deltaR_electron_L2.push_back(deltaR(p4, gen_L2));
            gen_deltaR_electron_L1FSR.push_back(deltaR(p4, gen_L1FSR));
            gen_deltaR_electron_L2FSR.push_back(deltaR(p4, gen_L2FSR));
        }
        for (auto p4: allmuons.gen_p4) {
            gen_deltaR_muon_L1.push_back(deltaR(p4, gen_L1));
            gen_deltaR_muon_L2.push_back(deltaR(p4, gen_L2));
            gen_deltaR_muon_L1FSR.push_back(deltaR(p4, gen_L1FSR));
            gen_deltaR_muon_L2FSR.push_back(deltaR(p4, gen_L2FSR));
        }
    } // end of if !event.isRealData()

    // ***** ***** *****
    // Trigger
    // ***** ***** *****
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

}

