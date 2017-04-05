#include <cp3_llbb/HtoZAAnalysis/interface/Types.h>
#include <vector>

namespace {
    struct dictionnary {
        HtoZA::Lepton dummy;
        std::vector<HtoZA::Lepton> dummy2;
        HtoZA::Dilepton dummy3;
        std::vector<HtoZA::Dilepton> dummy4;
        //HtoZA::Met dummy5;
        //std::vector<HtoZA::Met> dummy6;
        //HtoZA::DileptonMet dummy7;
        //std::vector<HtoZA::DileptonMet> dummy8;
        std::vector< std::vector<int> > dummy9;
        HtoZA::Jet dummy10;
        std::vector<HtoZA::Jet> dummy11;
        HtoZA::Dijet dummy12;
        std::vector<HtoZA::Dijet> dummy13;
        //HtoZA::DileptonMetDijet dummy14;
        //std::vector<HtoZA::DileptonMetDijet> dummy15;
        HtoZA::DileptonDijet dummy19;
        std::vector<HtoZA::DileptonDijet> dummy18;
        std::pair<int8_t, int8_t>  dummy16;
        HtoZA::MELAAngles dummy17;
    };
}
