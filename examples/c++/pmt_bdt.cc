#include <stdio.h>
#include <iostream>
#include <vector>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
//#include <RAT/DS/RunStore.hh>
//#include <RAT/DS/Run.hh>
//#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>

#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>

#include <glob.h>
#include <algorithm>
#include <numeric>

using namespace std;

int ENDpmt2dom(int pmt_id){
    return (int)pmt_id/31;
}

int GetDomX(int dom_id) {
    int dom_plane = dom_id/18;
    if (dom_plane == 0){ return 1; }
    if (dom_plane == 1){ return 0; }
    if (dom_plane == 2){ return 1; }
    if (dom_plane == 3){ return 2; }
    if (dom_plane == 4){ return 2; }
    return 0;
}

int GetDomY(int dom_id) {
    int dom_plane = dom_id/18;
    if (dom_plane == 0){ return 2; }
    if (dom_plane == 1){ return 1; }
    if (dom_plane == 2){ return 0; }
    if (dom_plane == 3){ return 0; }
    if (dom_plane == 4){ return 1; }
    return 0;
}

int GetDomZ(int dom_id) {
    return (int)dom_id % 18;
}

void process(std::string pattern, std::string outfilename){

    TFile* xsec = new TFile("xsec_graphs.root", "read");
    TGraph* gr_Si28 = (TGraph*)xsec->Get("nu_mu_Si28/tot_cc");
    TGraph* gr_O16 = (TGraph*)xsec->Get("nu_mu_O16/tot_cc");
    TGraph* gr_H1 = (TGraph*)xsec->Get("nu_mu_H1/tot_cc");
    xsec->Close();
    double xsec_weight = (gr_Si28->Eval(3.5) + 2.*gr_O16->Eval(3.5))/(gr_O16->Eval(3.5) + 2.*gr_H1->Eval(3.5));

    glob_t glob_result;
    int result = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(result != 0){
        std::cerr << "Couldn't find files! Exiting!" << std::endl;
        exit(1);
    }


    RAT::DSReader *dsreader = new RAT::DSReader(pattern.c_str());
    for(int i = 0; i < (int)glob_result.gl_pathc; i++)
        dsreader->Add(glob_result.gl_pathv[i]);
    const unsigned int nevents = dsreader->GetT()->GetEntries();
    std::cout << "NEvents : " << nevents << std::endl;

    // ttree at the dom level
    TH1D* hevts = new TH1D("nevts", "Total Evts", 1, 0, 1);
    TTree* tout = new TTree("doms", "doms");

    int dom_x, dom_y, dom_z;
    int npmts, npe, pe_min, pe_spread, pe_max;
    float pe_mean, pe_rms;
    float t_min, t_spread, t_mean, t_rms;
    int ndoms;
    int event_id;
    float vtxX, vtxY, vtxZ;
    float rock_wgt;

    tout->Branch("event_id", &event_id, "event_id/I");
    tout->Branch("vtxX", &vtxX, "vtxX/F");
    tout->Branch("vtxY", &vtxY, "vtxY/F");
    tout->Branch("vtxZ", &vtxZ, "vtxZ/F");
    tout->Branch("rock_wgt", &rock_wgt, "rock_wgt/F");
    tout->Branch("dom_x", &dom_x, "dom_x/I");
    tout->Branch("dom_y", &dom_y, "dom_y/I");
    tout->Branch("dom_z", &dom_z, "dom_z/I");
    tout->Branch("npmts", &npmts, "npmts/I");
    tout->Branch("npe", &npe, "npe/I");
    tout->Branch("pe_min", &pe_min, "pe_min/I");
    tout->Branch("pe_spread", &pe_spread, "pe_spread/I");
    tout->Branch("pe_rms", &pe_rms, "pe_rms/F");
    tout->Branch("t_min", &t_min, "t_min/F");
    tout->Branch("t_spread", &t_spread, "t_spread/F");
    tout->Branch("t_mean", &t_mean, "t_mean/F");
    tout->Branch("t_rms", &t_rms, "t_rms/F");

    // Loop over all triggered events
    double nwgt_events = 0;
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        if(rds->GetEVCount() == 0) continue;

        RAT::DS::MC *mc = rds->GetMC();
        RAT::DS::MCSummary *mcs = mc->GetMCSummary();
        RAT::DS::EV *ev = rds->GetEV(0);

        int mcpcount = mc->GetMCParticleCount();
        vtxX = -100000.;
        vtxY = -100000.;
        vtxZ = -100000.;
        for (int pid = 0; pid < mcpcount; pid++) {
            RAT::DS::MCParticle *particle = mc->GetMCParticle(pid);

            if (abs(particle->GetPDGCode()) == 13) {
                TVector3 mcpos = particle->GetPosition();
                vtxX = mcpos.X();
                vtxY = mcpos.Y();
                vtxZ = mcpos.Z();
            }
        }
        rock_wgt = 1.;
        if(((pattern.find("rockbed") != string::npos) || (pattern.find("equalx") != string::npos) || (pattern.find("aframe") != string::npos)) && (pattern.find("cosmic") == string::npos)){
            // up-weight rock events manually
            if(vtxY > -22000. && vtxY < -12000. && (pattern.find("short") == string::npos)){
                rock_wgt = xsec_weight;
            }
            if(vtxY > -22500. && vtxY < -12500. && (pattern.find("short") != string::npos)){
                rock_wgt = (float)xsec_weight;
            }
        }
        nwgt_events += rock_wgt;

        event_id = (int)iev;
        std::map<int, int> dom_pmts;
        std::map<int, std::vector<int>> dom_pes;
        std::map<int, std::vector<float>> dom_times;

        for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            // Total npe detected by a PMT
            int npe_i = mcpmt->GetMCPhotonCount();
            int pmt_i = mcpmt->GetID();
            int dom_i = ENDpmt2dom(pmt_i);
            // For End, type=0 is 8'', type=1 is 12''
            int type = mcpmt->GetType();
            if(type == 0) continue; // Skip the HQE PMTs in this example
            RAT::DS::PMT *pmt = ev->GetOrCreatePMT(pmt_i);
            if(npe_i >= 1) {
                dom_pes[dom_i].push_back(npe_i);
                dom_pmts[dom_i] += 1;
                dom_times[dom_i].push_back((float)pmt->GetTime());
            }
        }


        for(auto &dom: dom_pmts){
            dom_x = GetDomX(dom.first);
            dom_y = GetDomY(dom.first);
            dom_z = GetDomZ(dom.first);
            npmts = dom.second;
            std::vector<int> dom_pe = dom_pes[dom.first];
            std::vector<float> dom_time = dom_times[dom.first];

            std::sort(dom_pe.begin(), dom_pe.end());
            std::sort(dom_time.begin(), dom_time.end());

            npe = std::accumulate(dom_pe.begin(), dom_pe.end(), 0);
            t_mean = std::accumulate(dom_time.begin(), dom_time.end(), 0)/dom_time.size();
            pe_mean = npe/dom_pe.size();
            pe_min = dom_pe.at(0);
            pe_max = dom_pe.at(dom_pe.size()-1);
            pe_spread = pe_max - pe_min;
            t_spread = dom_time.at(dom_time.size()-1) - dom_time.at(0);
            t_min = dom_time.at(0);
            pe_rms = 0.;
            t_rms = 0.;
            for(int i = 0; i < (int)dom_pe.size(); i++){
                pe_rms += (TMath::Power(dom_pe[i]-pe_mean, 2));
                t_rms += (TMath::Power(dom_time[i]-t_mean, 2));
            }
            pe_rms = TMath::Sqrt(pe_rms)/dom_pe.size();
            t_rms = TMath::Sqrt(t_rms)/dom_time.size();

            tout->Fill();

        }
    }
    hevts->SetBinContent(1, nwgt_events);
    std::cout << "Writing to File!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    tout->Write();
    hevts->Write();
    outFile->Close();
}


int main(int argc, char *argv[]){

    if(argc == 3){
        std::string input_pattern(argv[1]);
        std::string outfilename(argv[2]);
        process(input_pattern, outfilename);
    }
    else{
        std::cout << "Wrong number of arguments." << std::endl;
    }

    return 0;
}
