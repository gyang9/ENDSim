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

fvs = {"pen_small" : np.array([[-22, 22], [-22, 22], [-20, 90]]),
       #  "pen_small_equalx" : np.array([[-22, 22], [-22, 22], [-20, 90]]),
       "pen_small_equalx_short" : np.array([[-22.5, 22.5], [-22.5, 22.5], [-20, 64.5]]),
       "aframe_spacing5m" : np.array([[-23, 23], [-23, 23], [-20, 70]])
       #  "trans" : np.array([[-20, 26], [-20, 90], [-20, 20]]),
       #  "long" : np.array([[-20, 26], [-20, 20], [-20, 90]])
      }

def mass(V, conf="long"):
    if "rockbed" not in conf and "equalx" not in conf and "aframe" not in conf:
        return (V[0][1] - V[0][0])*(V[1][1] - V[1][0])*(V[2][1] - V[2][0])*997./1.E6
    #  if "rockbed" not in conf:
    #      return (V[0][1] - V[0][0])*(V[1][1] - V[1][0])*(V[2][1] - V[2][0])*997./1.E6
    else:
        return (V[0][1] - V[0][0])*(V[2][1] - V[2][0])*((V[1][1] - (V[1][0]+10))*997. + (10)*2700.)/1.E6

#  mass = lambda V : (V[0][1] - V[0][0])*(V[1][1] - V[1][0])*(V[2][1] - V[2][0])/997.

frates = ROOT.TFile("/home/nitish/public_html/shared/uboone/flugg_study/for_milind/plots_pretty/data/event_rates_lp3_depth_numu.root", "read")
rate_perPOT_perKT = frates.Get("h_numu_lp3_me").Integral()/10.

rate_perPOT = {conf: rate_perPOT_perKT*mass(fvs[conf], conf) for conf in fvs}
print([mass(fvs[conf], conf) for conf in fvs])
print(rate_perPOT)
pots = {}
livetime = {}
cols = {"pen_small": ROOT.kBlack, "trans" : ROOT.kRed, "aframe_spacing5m" : ROOT.kGreen+2, "pen_small_equalx_short": ROOT.kViolet}
legends = {"pen_small" : "2m Pentagon (3D)", "trans" : "Vertical (2D)", "long": "Horizontal (2D)", "pen_small_equalx_short": "2.5m Pentagon (3D), Wider Bottom, 3m Z-spacing", "aframe_spacing5m": "A Frame, 5m Spacing"}

hdoms_perPOT = {}
hdoms_perDay = {}
for conf in fvs:
    fin = ROOT.TFile("hists_numu_%s_domPE.root"% conf, "read")
    htot = fin.Get("fid_tot_xy")
    htot.SetDirectory(0)

    pot = htot.Integral()/rate_perPOT[conf]
    if ("rockbed" in conf or "equalx" in conf) and "short" not in conf:
        pot = pots["pen_small"]
    #  if "rockbed" in conf:
    #      pot = pots["pen_small"]

    pots[conf] = pot
    livetime[conf] = pot/0.025 # pot -> day conversion

    hdom = fin.Get("hnDoms")
    hdom.SetDirectory(0)

    hcumdom = hdom.GetCumulative()
    for i in range(1, hcumdom.GetNbinsX()+1):
        hcumdom.SetBinContent(i, htot.Integral()-hcumdom.GetBinContent(i))
        print(hcumdom.GetBinContent(i)/pot)

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
    #  hdom_perDay.SetBinContent(1, 0)
    hdom_perDay.Scale(0.025/pot)
    hdom_perDay.GetYaxis().SetTitle("Cumulative Triggered Events per Day")
    #  hdom_perDay.GetXaxis().SetTitle("N DOMs with 5L0 Trigger")
    hdom_perDay.GetXaxis().SetTitle("N DOMs with >=3 PE")
    #  hdom_perDay.GetYaxis().SetRangeUser(0, hdom_perDay.GetMaximum()*1.5)
    hdom_perDay.SetLineWidth(2)
    hdom_perDay.SetLineColor(cols[conf])
    hdom_perDay.SetDirectory(0)
    hdoms_perDay[conf] = hdom_perDay

    print(conf)
    print("-------")
    #  for i in range(1, hdom_perPOT.GetNbinsX()+1):
    #      print(i, hdom_perPOT.GetBinContent(i))
    #  print(hdom_perPOT.GetBinContent(1))
    #  print(hdom_perPOT.GetBinContent(2))
    #  print("-------")

print(pots)
c = ROOT.TCanvas("c", "", 2000, 800)
c.Divide(2, 1)
leg = ROOT.TLegend(0.45, 0.6, 0.85, 0.85)
leg.SetBorderSize(0)
for conf in fvs:
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.125)
    hdoms_perPOT[conf].Draw("hist same")
    leg.AddEntry(hdoms_perPOT[conf], legends[conf], "le")
    #  ROOT.gPad.SetLogy()

    c.cd(2)
    hdoms_perDay[conf].Draw("hist same")
    #  ROOT.gPad.SetLogy()
leg.Draw()
c.Print(outfolder+"ndoms_metric_aframe.pdf")
