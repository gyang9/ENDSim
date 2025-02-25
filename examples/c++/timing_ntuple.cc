#include <iostream>
#include <string>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

// Proces the input ntuple file
void process(const char* filename){

    TH1D* ht = new TH1D("","",600,-20,100);
    TH1D* ht_digit = new TH1D("","",600,-20,100);

    TH1D* hpe_time = new TH1D("","",600,-20,100);
    TH1D* hpe_ch_time = new TH1D("","",600,-20,100);
    TH1D* hpe_sc_time = new TH1D("","",600,-20,100);

    TFile *f = TFile::Open(filename);
    TTree *t = (TTree*)f->Get("output"); 
    TTree *meta = (TTree*)f->Get("meta"); 

    const double OFFSET = 67.7;

    int mcNPE;
    int mcNhits;
    std::vector<double> *hitPMTTime = new std::vector<double>();
    std::vector<double> *hitPMTDigitizedTime = new std::vector<double>();
    std::vector<double> *mcPETime = new std::vector<double>();
    std::vector<int>  *mcPEProcess = new std::vector<int>();

    std::vector<int> *hitPMTID = new std::vector<int>();
    std::vector<double> *pmtX = new std::vector<double>();
    std::vector<double> *pmtY = new std::vector<double>();
    std::vector<double> *pmtZ = new std::vector<double>();

    t->SetBranchAddress("hitPMTID", &hitPMTID);
    t->SetBranchAddress("mcNPE", &mcNPE);
    t->SetBranchAddress("mcNhits", &mcNhits);  
    t->SetBranchAddress("hitPMTTime", &hitPMTTime);
    t->SetBranchAddress("mcPETime", &mcPETime);
    t->SetBranchAddress("mcPEProcess", &mcPEProcess);
    t->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);

    meta->SetBranchAddress("pmtX", &pmtX);
    meta->SetBranchAddress("pmtY", &pmtY);
    meta->SetBranchAddress("pmtZ", &pmtZ);

    for(int iEV = 0; iEV < t->GetEntries(); iEV++){

        t->GetEntry(iEV);
        meta->GetEntry(iEV);

        // Loop over the hit PMTs (triggered events)
        for(size_t i = 0; i < hitPMTTime->size(); i++){
            double time = (*hitPMTTime)[i];
            double digitized_time = (*hitPMTDigitizedTime)[i];
            // Example of how to grab the PMT position
            int id = (*hitPMTID)[i];
            double x = (*pmtX)[id];
            double y = (*pmtY)[id];
            double z = (*pmtZ)[id];
            ht->Fill(time);
            ht_digit->Fill(digitized_time - OFFSET);
        }

        // Loop over the total number of PEs
        // mcPETime->size() == mcNPE
        for(size_t i = 0; i < mcPETime->size(); i++){
            double time = (*mcPETime)[i];
            int process = (*mcPEProcess)[i];
            hpe_time->Fill(time);
            // These PE were created by Cherenkov photons
            if(process){
                hpe_ch_time->Fill(time);
            }
            // These PE were create by scint. or reemitted photons
            else{
                hpe_sc_time->Fill(time);
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    double m = std::max(ht->GetMaximum(), ht_digit->GetMaximum());
    ht->GetYaxis()->SetRangeUser(0, 1.1*m);
    ht->GetXaxis()->SetRangeUser(-5, 30);
    ht->SetStats(0);
    ht->SetLineColor(kBlack);
    ht->GetXaxis()->SetTitleFont(132);
    ht->GetXaxis()->SetLabelFont(132);
    ht->GetYaxis()->SetTitleFont(132);
    ht->GetYaxis()->SetLabelFont(132);
    ht->GetXaxis()->SetTitle("Time (ns)");
    ht->GetYaxis()->SetTitle("Counts / 0.2 ns");
    ht->Draw("");

    ht_digit->SetLineColor(kRed);
    ht_digit->Draw("same");

    TLegend *t1 = new TLegend(0.5, 0.5, 0.85, 0.85);
    t1->SetTextFont(132);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->AddEntry(ht, "PMT hit times");
    t1->AddEntry(ht_digit, "Digitizer hit times");
    t1->Draw("same");

    c1->Update();

    TCanvas *c2 = new TCanvas("c2","c2",800,600);

    hpe_time->GetXaxis()->SetRangeUser(0, 30);
    hpe_time->SetStats(0);
    hpe_time->SetLineColor(kBlack);
    hpe_time->GetXaxis()->SetTitleFont(132);
    hpe_time->GetXaxis()->SetLabelFont(132);
    hpe_time->GetYaxis()->SetTitleFont(132);
    hpe_time->GetYaxis()->SetLabelFont(132);
    hpe_time->GetXaxis()->SetTitle("Time (ns)");
    hpe_time->GetYaxis()->SetTitle("Counts / 0.2 ns");
    hpe_time->SetLineColor(kBlack);
    hpe_time->Draw("same");

    hpe_ch_time->SetLineColor(kBlue);
    hpe_ch_time->Draw("same");

    hpe_sc_time->SetLineColor(kRed);
    hpe_sc_time->Draw("same");

    TLegend *t2 = new TLegend(0.5, 0.5, 0.85, 0.85);
    t2->SetTextFont(132);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->AddEntry(hpe_time, "All PE hit times");
    t2->AddEntry(hpe_ch_time, "Cher. PE hit times");
    t2->AddEntry(hpe_sc_time, "Scint. PE hit times");
    t2->Draw("same");

    c2->Update();

    c1->SaveAs("pmt_hit_times.png");
    c2->SaveAs("pe_hit_times.png");
}

int main(int argc, char** argv){

    if(argc == 2){
        process(argv[1]);
    }
    else{
        std::cout << "Wrong number of arguments." << std::endl; 
    }

    return 0;
}
