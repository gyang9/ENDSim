#include <stdio.h>
#include <iostream>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/Root.hh>

#include <TCanvas.h>
#include <TH1D.h>

using namespace std;

void process(std::string filename){

    TH1D* ht = new TH1D("","",600,-20,100);
    TH1D* hq = new TH1D("","",250,-5,20);

    RAT::DSReader *dsreader = new RAT::DSReader(filename.c_str());
    const unsigned int nevents = dsreader->GetT()->GetEntries();

    int ch_total = 0;
    int ch_det = 0;

    // Loop over all triggered events
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        RAT::DS::MC *mc = rds->GetMC();

        // Use the MCSummary to get total # ch photons produced
        RAT::DS::MCSummary *mcs = mc->GetMCSummary();
        ch_total += mcs->GetNumCerenkovPhoton();

        // Loop over all hit PMTs
        for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            // Total npe detected by a PMT
            int npes = mcpmt->GetMCPhotonCount();

            // For End, type=0 is 8'', type=1 is 12''
            int type = mcpmt->GetType();
            if(type == 0) continue; // Skip the HQE PMTs in this example

            // Loop over all PEs detected by the PMT
            for(int ipe = 0; ipe < npes; ipe++){
                double time = mcpmt->GetMCPhoton(ipe)->GetFrontEndTime();
                double charge = mcpmt->GetMCPhoton(ipe)->GetCharge();
                bool isDark = mcpmt->GetMCPhoton(ipe)->IsDarkHit();
                if(isDark) continue; // Skip the Dark hits in this example
                ht->Fill(time);
                hq->Fill(charge);
                std::string process = mcpmt->GetMCPhoton(ipe)->GetCreatorProcess(); 
                if(strcmp(process.c_str(), "Cerenkov")==0){
                    ch_det += 1;
                }
            }
        }
    }

    std::cout << "Total Created Cherenkov photons/event: " << ch_total/nevents << std::endl;
    std::cout << "Total Detected Chernekov photons/event: " << ch_det/nevents << std::endl;

    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    ht->GetXaxis()->SetTitleFont(132);
    ht->GetXaxis()->SetLabelFont(132);
    ht->GetYaxis()->SetTitleFont(132);
    ht->GetYaxis()->SetLabelFont(132);
    ht->GetXaxis()->SetTitle("Hit times (ns)");
    ht->GetYaxis()->SetTitle("Arb. Units");
    ht->GetYaxis()->SetTitleOffset(1.25);
    ht->SetStats(0);

    ht->GetXaxis()->SetRangeUser(-5.0, 45.0);
    ht->SetLineColor(kBlack);
    ht->Draw("E HIST");

    c1->Update();

    c1->SaveAs("timing.png");

    TCanvas *c2 = new TCanvas("c2","c2",800,600);

    hq->GetXaxis()->SetTitleFont(132);
    hq->GetXaxis()->SetLabelFont(132);
    hq->GetYaxis()->SetTitleFont(132);
    hq->GetYaxis()->SetLabelFont(132);
    hq->GetXaxis()->SetTitle("Charge (pC)");
    hq->GetYaxis()->SetTitle("Arb. Units");
    hq->GetYaxis()->SetTitleOffset(1.25);
    hq->SetStats(0);

    hq->SetLineColor(kBlack);
    hq->Draw("E HIST");

    c2->Update();

    c2->SaveAs("charge.png");

}

int main(int argc, char *argv[]){

    if(argc == 2){
        std::string input_filename(argv[1]);
        process(input_filename);
    }
    else{
        std::cout << "Wrong number of arguments." << std::endl;
    }

    return 0;
}

