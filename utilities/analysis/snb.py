import ROOT
import numpy as np

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gROOT.SetBatch(True)

gLeaked=[]
def Clone(obj, *args):
    ret = obj.Clone(*args)
    gLeaked.append(ret)
    return ret

def New(cls, *args):
    ret = cls(*args)
    gLeaked.append(ret)
    return ret


outfolder = "/home/nitish/public_html/shared/end/efficiency/domPE/"

fvs =  {
       #  "pen_small_equalx_short" : np.array([[-22.5, 22.5], [-22.5, 22.5], [-20, 64.5]])
       "aframe_spacing5m" : np.array([[-23, 23], [-23, 23], [-20, 70]])
       }

def mass(V, conf="long"):
    if "rockbed" not in conf and "equalx" not in conf and "aframe" not in conf:
        return (V[0][1] - V[0][0])*(V[1][1] - V[1][0])*(V[2][1] - V[2][0])*997./1.E6
    else:
        return (V[0][1] - V[0][0])*(V[2][1] - V[2][0])*((V[1][1] - (V[1][0]+10))*997. + (10)*2700.)/1.E6

frates = ROOT.TFile("/home/nitish/public_html/shared/uboone/flugg_study/for_milind/plots_pretty/data/event_rates_lp3_depth_numu.root", "read")
rate_perPOT_perKT = frates.Get("h_numu_lp3_me").Integral()/10.

rate_perPOT = {conf: rate_perPOT_perKT*mass(fvs[conf], conf) for conf in fvs}
print([mass(fvs[conf], conf) for conf in fvs])
pots = {}
cols = {"pen_small": ROOT.kBlack, "trans" : ROOT.kRed, "aframe_spacing5m" : ROOT.kGreen+2, "pen_small_equalx_short": ROOT.kViolet}
legends = {"pen_small" : "2m Pentagon (3D)", "trans" : "Vertical (2D)", "long": "Horizontal (2D)", "pen_small_equalx_short": "2.5m Pentagon (3D), Wider Bottom, 3m Z-spacing", "aframe_spacing5m": "A Frame, 5m Spacing"}

cosmic_surface = {
                    #  "pen_small_equalx_short" : np.array([[-30., 30.], [-30., 100.]]),
                    "aframe_spacing5m" : np.array([[-30., 30.], [-30., 80.]])
                 }
def area(surface):
    return (surface[0][1]-surface[0][0])*(surface[1][1]-surface[1][0])
cosmic_rate = {conf: 1.5*area(cosmic_surface[conf]) for conf in cosmic_surface}

hdoms_perPOT = {}
hdoms_perDay = {}
hdoms_c_perDay = {}
for conf in fvs:
    fin = ROOT.TFile("hists_numu_%s_domPE.root"% conf, "read")
    htot = fin.Get("fid_tot_xy")
    htot.SetDirectory(0)

    pot = htot.Integral()/rate_perPOT[conf]

    pots[conf] = pot

    hdom = fin.Get("hnDoms")
    hdom.SetDirectory(0)

    hcumdom = hdom.GetCumulative()
    for i in range(1, hcumdom.GetNbinsX()+1):
        hcumdom.SetBinContent(i, htot.Integral()-hcumdom.GetBinContent(i))

    hdom_perPOT = Clone(hcumdom, "dom_perPOT_%s"%conf)
    #  hdom_perPOT.SetBinContent(1, 0)
    hdom_perPOT.Scale(1./pot)
    hdom_perPOT.GetYaxis().SetTitle("Cumulative Triggered Events per 10^{20} POT")
    #  hdom_perPOT.GetXaxis().SetTitle("N DOMs with 5L0 Trigger")
    hdom_perPOT.GetXaxis().SetTitle("N DOMs with >=3 PE")
    #  hdom_perPOT.GetYaxis().SetRangeUser(0, hdom_perPOT.GetMaximum()*1.5)
    hdom_perPOT.SetLineWidth(2)
    hdom_perPOT.SetLineColor(cols[conf])
    hdom_perPOT.SetDirectory(0)
    hdoms_perPOT[conf] = hdom_perPOT

    hdom_perDay = Clone(hcumdom, "dom_perDay_%s"%conf)
    hdom_perDay.Scale(0.025/pot)
    hdom_perDay.GetYaxis().SetTitle("Cumulative Triggered Events per Day")
    hdom_perDay.GetXaxis().SetTitle("N DOMs with >=3 PE")
    hdom_perDay.SetLineWidth(2)
    hdom_perDay.SetLineColor(ROOT.kRed)
    hdom_perDay.SetDirectory(0)
    hdoms_perDay[conf] = hdom_perDay

    fin2 = ROOT.TFile("hists_cosmics_%s_domPE.root"%conf, "read")
    htot2 = fin2.Get("fid_tot_xy")
    htot2.SetDirectory(0)
    lt = htot2.Integral()/cosmic_rate[conf] # in seconds
    print(lt)
    #  livetime[conf] = lt

    hdom_c = fin2.Get("hnDoms")
    hdom_c.SetDirectory(0)

    hcumdom_c = hdom_c.GetCumulative()
    for i in range(1, hcumdom_c.GetNbinsX()+1):
        hcumdom_c.SetBinContent(i, htot2.Integral()-hcumdom_c.GetBinContent(i))
        print(hcumdom_c.GetBinContent(i)*86400./lt)

    print(86400*0.937*1.E-5/lt)

    hdom_c_perDay = Clone(hcumdom_c, "dom_cosmics_perDay_%s"%conf)
    hdom_c_perDay.Scale(86400*0.937*1.E-5/lt)
    hdom_c_perDay.GetYaxis().SetTitle("Cumulative Cosmic Triggered Events per Day")
    hdom_c_perDay.GetXaxis().SetTitle("N DOMs with >=3 PE")
    hdom_c_perDay.SetLineWidth(2)
    hdom_c_perDay.SetLineColor(ROOT.kBlue)
    hdom_c_perDay.SetDirectory(0)
    hdoms_c_perDay[conf] = hdom_c_perDay


    print(conf)
    print("-------")

print(pots)
for conf in fvs:
    c = ROOT.TCanvas("c", "", 2000, 800)
    c.Divide(2, 1)
    leg = ROOT.TLegend(0.45, 0.6, 0.85, 0.85)
    leg.SetBorderSize(0)
    c.cd(1)
    hdoms_c_perDay[conf].Draw("hist same")
    hdoms_perDay[conf].Draw("hist same")
    leg.AddEntry(hdoms_c_perDay[conf], "Cosmics", "le")
    leg.AddEntry(hdoms_perDay[conf], "Beam", "le")
    ROOT.gPad.SetLogy()
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.125)
    hratio = Clone(hdoms_perDay[conf], "%s_ratio"%conf)
    hratio.Divide(hdoms_c_perDay[conf])
    hratio.GetYaxis().SetTitle("Signal to Noise Ratio")
    hratio.Draw("hist same")
    #  ROOT.gPad.SetLogy()

    c.cd(1)
    leg.Draw()
    c.Print(outfolder+"ndoms_metric_%s_cosmics.pdf" % conf)
