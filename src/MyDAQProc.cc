#include "MyDAQProc.hh"
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPhoton.hh>
#include <RAT/Log.hh>

namespace END {

MyDAQProc::MyDAQProc() : RAT::Processor("mydaq") {
    fEventCounter = 0;
}

MyDAQProc::~MyDAQProc() {}

RAT::Processor::Result MyDAQProc::DSEvent(RAT::DS::Root *ds) {
    if (ds->ExistEV()) {
        ds->PruneEV();
    }
    RAT::DS::EV *ev = ds->AddNewEV();
    ev->SetID(fEventCounter++);

    RAT::DS::MC *mc = ds->GetMC();
    if (!mc) {
        return RAT::Processor::OK;
    }

    for (int i = 0; i < mc->GetMCPMTCount(); i++) {
        RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(i);
        if (mcpmt->GetMCPhotonCount() > 0) {
            RAT::DS::PMT *pmt = ev->GetOrCreatePMT(mcpmt->GetID());
            
            double time = mcpmt->GetMCPhoton(0)->GetFrontEndTime();
            double charge = 0;

            for (int j = 0; j < mcpmt->GetMCPhotonCount(); j++) {
                if (time > mcpmt->GetMCPhoton(j)->GetFrontEndTime()) {
                    time = mcpmt->GetMCPhoton(j)->GetFrontEndTime();
                }
                charge += mcpmt->GetMCPhoton(j)->GetCharge();
            }
            
            pmt->SetTime(time);
            pmt->SetCharge(charge);
        }
    }

    return RAT::Processor::OK;
}

} // namespace END