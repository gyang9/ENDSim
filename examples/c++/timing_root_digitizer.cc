#include <stdio.h>
#include <iostream>
#include <string>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Digit.hh>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

using namespace std;

void process(std::string filename){

    TH1D* hwaveform = new TH1D("","",512,0,511);

    RAT::DSReader *dsreader = new RAT::DSReader(filename.c_str());
    const unsigned int nevents = dsreader->GetT()->GetEntries();

    // Loop over all triggered events
    for(size_t iev = 0; iev < nevents; iev++){

        // Read in the event and get the digitizer information
        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        if(!rds->ExistEV()) continue;
        RAT::DS::EV *ev = rds->GetEV(0);
        RAT::DS::MC *mc = rds->GetMC();

        if(!ev->DigitizerExists()) continue;

        RAT::DS::Digit digitizer = ev->GetDigitizer();

        // Get the digitizer dynamic range and the number of bits
        double dynamic_range = digitizer.GetDynamicRange(); // in mV
        int nbits = digitizer.GetNBits();
        double voltage_res = dynamic_range/(1 << nbits);

        // Digitizer sampling rate and number of samples
        float sampling_rate = digitizer.GetSamplingRate();
        double time_step = 1.0/sampling_rate;
        int nsamples = digitizer.GetNSamples();

        // Loop over the true set of hit PMTs
        for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            int npe = mcpmt->GetMCPhotonCount();
            int pmtID = mcpmt->GetID();
        }

        // Loop over the hit PMTs that are built into a triggered event
        for(int iPMT = 0; iPMT < ev->GetPMTCount(); iPMT++){

            // Get the waveform for each PMT
            RAT::DS::PMT* pmt = ev->GetPMT(iPMT);
            int pmtID = pmt->GetID();
            std::vector<UShort_t> waveform = digitizer.GetWaveform(pmtID);

            double charge = pmt->GetDigitizedCharge();
            int crossings = pmt->GetNCrossings();
            double t_over_thresh = pmt->GetTimeOverThreshold();
            double pedestal = pmt->GetPedestal();
            double peaks = -1.0*pmt->GetPeakVoltage();

            for(size_t sample = 0; sample < waveform.size(); sample++){
                double time = sample*time_step;
                // The waveform contains the ADC counts for each sample
                int adc = int(waveform[sample]);
                // Convert to voltage
                double voltage = (adc - pedestal)*voltage_res;
                hwaveform->SetBinContent(sample, voltage);
            }
        }
    }


    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    hwaveform->GetXaxis()->SetRangeUser(20, 140);
    hwaveform->GetXaxis()->SetTitleFont(132);
    hwaveform->GetXaxis()->SetLabelFont(132);
    hwaveform->GetYaxis()->SetTitleFont(132);
    hwaveform->GetYaxis()->SetLabelFont(132);
    hwaveform->GetXaxis()->SetTitle("Sample");
    hwaveform->GetYaxis()->SetTitle("Voltage (mV)");
    hwaveform->GetYaxis()->SetTitleOffset(1.25);
    hwaveform->SetStats(0);

    hwaveform->SetLineColor(kBlack);
    hwaveform->Draw("");

    c1->Update();

    c1->SaveAs("waveform.png");
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

