#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/HtoZAAnalysis/interface/Categories.h>

#include <regex>

// ***** ***** *****
// Dilepton categories
// ***** ***** *****

static std::regex s_mumu_hlt_regex("^HLT_Mu.*_(Tk)?Mu");
static std::regex s_elel_hlt_regex("^HLT_Ele.*_Ele");
static std::regex s_muel_elmu_hlt_regex("^HLT_Mu.*_Ele");

const std::vector<HtoZA::Lepton>& DileptonCategory::getLeptons(const AnalyzersManager& analyzers) const {
    const HtoZAAnalyzer& hZA_analyzer = analyzers.get<HtoZAAnalyzer>(m_analyzer_name);
    return hZA_analyzer.leptons;
}

const std::vector<HtoZA::Dilepton>& DileptonCategory::getDileptons(const AnalyzersManager& analyzers) const {
    const HtoZAAnalyzer& hZA_analyzer = analyzers.get<HtoZAAnalyzer>(m_analyzer_name);
    return hZA_analyzer.ll;
}

const std::vector<HtoZA::DileptonDijet>& DileptonCategory::getDileptonDijets(const AnalyzersManager& analyzers) const {
    const HtoZAAnalyzer& hZA_analyzer = analyzers.get<HtoZAAnalyzer>(m_analyzer_name);
    return hZA_analyzer.lljj;
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****

void MuMuCategory::configure(const edm::ParameterSet& conf) {
    DileptonCategory::configure(conf);

    m_leadingLeptonPtCut = conf.getUntrackedParameter<double>("mumu_leadingLeptonPtCut");
    m_subleadingLeptonPtCut = conf.getUntrackedParameter<double>("mumu_subleadingLeptonPtCut");
}

bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return (muons.p4.size() >= 2);
};

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HtoZA::Lepton>& leptons = getLeptons(analyzers);
    const std::vector<HtoZA::Dilepton>& ll = getDileptons(analyzers);
    const std::vector<HtoZA::DileptonDijet>& lljj = getDileptonDijets(analyzers);

    if (ll.empty())
        return false;

    if (lljj.empty())
        return false;

    // Only look at the first dilepton pair
    return ll[0].isMuMu &&
        (leptons[ll[0].ilep1].p4.Pt() > m_leadingLeptonPtCut) &&
        (leptons[ll[0].ilep2].p4.Pt() > m_subleadingLeptonPtCut);
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger", "HLT_Mu*");
};

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) {
        if (std::regex_search(path, s_mumu_hlt_regex)) {
            manager.pass_cut("fire_trigger");
            break;
        }
    }
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****

void ElElCategory::configure(const edm::ParameterSet& conf) {
    DileptonCategory::configure(conf);

    m_leadingLeptonPtCut = conf.getUntrackedParameter<double>("elel_leadingLeptonPtCut");
    m_subleadingLeptonPtCut = conf.getUntrackedParameter<double>("elel_subleadingLeptonPtCut");
}

bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return (electrons.p4.size() >= 2);
};

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HtoZA::Lepton>& leptons = getLeptons(analyzers);
    const std::vector<HtoZA::Dilepton>& ll = getDileptons(analyzers);
    const std::vector<HtoZA::DileptonDijet>& lljj = getDileptonDijets(analyzers);

    if (ll.empty())
        return false;

    if (lljj.empty())
        return false;

    // Only look at the first dilepton pair
    return ll[0].isElEl &&
        (leptons[ll[0].ilep1].p4.Pt() > m_leadingLeptonPtCut) &&
        (leptons[ll[0].ilep2].p4.Pt() > m_subleadingLeptonPtCut);
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger", "HLT_Ele*");
};

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) {
        if (std::regex_search(path, s_elel_hlt_regex)) {
            manager.pass_cut("fire_trigger");
            break;
        }
    }
}


/*
// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****

void ElMuCategory::configure(const edm::ParameterSet& conf) {
    DileptonCategory::configure(conf);

    m_leadingLeptonPtCut = conf.getUntrackedParameter<double>("elmu_leadingLeptonPtCut");
    m_subleadingLeptonPtCut = conf.getUntrackedParameter<double>("elmu_subleadingLeptonPtCut");
}

bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return ((electrons.p4.size() + muons.p4.size()) >= 2);
};

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HtoZA::Lepton>& leptons = getLeptons(analyzers);
    const std::vector<HtoZA::Dilepton>& ll = getDileptons(analyzers);
    const std::vector<HtoZA::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);

    if (ll.empty())
        return false;

    if (llmetjj.empty())
        return false;

    // Only look at the first dilepton pair
    return ll[0].isElMu &&
        (leptons[ll[0].ilep1].p4.Pt() > m_leadingLeptonPtCut) &&
        (leptons[ll[0].ilep2].p4.Pt() > m_subleadingLeptonPtCut);
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger", "HLT_Mu*Ele*");
};

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) {
        if (std::regex_search(path, s_muel_elmu_hlt_regex)) {
            manager.pass_cut("fire_trigger");
            break;
        }
    }
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****

void MuElCategory::configure(const edm::ParameterSet& conf) {
    DileptonCategory::configure(conf);

    m_leadingLeptonPtCut = conf.getUntrackedParameter<double>("muel_leadingLeptonPtCut");
    m_subleadingLeptonPtCut = conf.getUntrackedParameter<double>("muel_subleadingLeptonPtCut");
}

bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return ((electrons.p4.size() + muons.p4.size()) >= 2);
};

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HtoZA::Lepton>& leptons = getLeptons(analyzers);
    const std::vector<HtoZA::Dilepton>& ll = getDileptons(analyzers);
    const std::vector<HtoZA::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);

    if (ll.empty())
        return false;

    if (llmetjj.empty())
        return false;

    // Only look at the first dilepton pair
    return ll[0].isMuEl &&
        (leptons[ll[0].ilep1].p4.Pt() > m_leadingLeptonPtCut) &&
        (leptons[ll[0].ilep2].p4.Pt() > m_subleadingLeptonPtCut);
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger", "HLT_Mu*Ele*");
};

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) {
        if (std::regex_search(path, s_muel_elmu_hlt_regex)) {
            manager.pass_cut("fire_trigger");
            break;
        }
    }
}
*/
