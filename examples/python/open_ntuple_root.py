import ROOT
import sys

def pretty_plot(h):

    h.SetStats(0)
    h.GetXaxis().SetTitleFont(132)
    h.GetXaxis().SetLabelFont(132)
    h.GetYaxis().SetTitleFont(132)
    h.GetYaxis().SetLabelFont(132)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitle("Counts")


def open_ntuple(filename):

    # Arb. time offset to align peaks
    # in true time and digitized time
    OFFSET = 67.7

    ht1 = ROOT.TH1D("","",500,-20,80)
    ht2 = ROOT.TH1D("","",500,-20,80)

    hq1 = ROOT.TH1D("","",120,-2,10)
    hq2 = ROOT.TH1D("","",120,-2,10)

    f = ROOT.TFile.Open(filename) # Input the ntuple file
    tree = f.Get("output")
    meta = f.Get("meta")

    for i in range(tree.GetEntries()):

        try:
            tree.GetEntry(i)
            meta.GetEntry(i)
        except AttributeError:
            print ("Error with opening file!")
            sys.exit(1)

        pmtid = list(tree.hitPMTID)
        # Hit time of first PE for each PMT
        pmttime = list(tree.hitPMTTime)
        # Hit time from processed waveform
        pmtdigittime = list(tree.hitPMTDigitizedTime)
        pmtcharge = list(tree.hitPMTCharge)
        pmtdigitcharge = list(tree.hitPMTDigitizedCharge)

        pmtx = list(meta.pmtX)
        pmty = list(meta.pmtY)
        pmtz = list(meta.pmtZ)

        # Loop over the PMT hits
        for i in range(len(pmtid)):
            ht1.Fill(pmttime[i])
            hq1.Fill(pmtcharge[i])
            # Information from the processed waveforms
            ht2.Fill(pmtdigittime[i] - OFFSET)
            hq2.Fill(pmtdigitcharge[i])

            # How to get the PMT positions
            pid = pmtid[i]
            pmt_x = pmtx[pid]
            pmt_y = pmty[pid]
            pmt_z = pmtz[pid]


    return ht1, ht2, hq1, hq2


if __name__=='__main__':


    ht1, ht2, hq1, hq2 = open_ntuple(sys.argv[1])

    for h in [ht1, hq1, ht2, hq2]:
        pretty_plot(h)

    c1 = ROOT.TCanvas("c1","c1",800,600)

    ht1.GetXaxis().SetTitle("Time (ns)")
    ht1.SetLineColor(ROOT.kBlack)
    ht1.Draw("")
    ht2.SetLineColor(ROOT.kRed)
    ht2.Draw("same")

    c1.Update()

    c2 = ROOT.TCanvas("c2","c2",800,600)

    hq1.GetXaxis().SetTitle("Charge (pC)")
    hq1.SetLineColor(ROOT.kBlack)
    hq1.Draw("")
    hq2.SetLineColor(ROOT.kRed)
    hq2.Draw("same")

    c2.Update()

    input()

