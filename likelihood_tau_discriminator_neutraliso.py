import ROOT
import math
import numpy as np
import copy
from CMGTools.H2TauTau.proto.plotter.ROCPlotter import *

ROOT.gROOT.SetBatch(True) # don't disply canvas

create_maps = False
create_discriminator = True
create_plots = True

#settings for track maps
nbins_pt = 22
nbins_dr = 60
pt_binning = np.array([0.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2] + [1.0+39.0*0.7**x for x in reversed(range(10))])

#settings for discriminator hists
upperlim = 20.0
nbins_discrim = 1000

file0 = ROOT.TFile("/eos/user/s/swozniew/TauPOG_files/ntuples/DYJetsToLL_M50_LO.root")
tree = file0.Get("tree")
c1 = ROOT.TCanvas()

if create_maps:
    file1 = ROOT.TFile("photon_maps.root", "RECREATE")
    for charge in [["P", 1.0], ["M", -1.0]]:
        for PUwin in [[0, 1], [30, 35], [60, 70]]:
            print "Creating maps for PU {down} to {up} and charge {ch} ...".format(down=PUwin[0], up=PUwin[1], ch=charge[0])
            
            #count event yields
            n_sig = ROOT.TH1I("n_sig", "n_sig", 1, 0, 2)
            tree.Draw("1>>n_sig", "tau_gen_match==5&&n_true_interactions>={down}&&n_true_interactions<={up}&&tau_charge=={ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[1]))
            n_bkg = ROOT.TH1I("n_bkg", "n_bkg", 1, 0, 2)
            tree.Draw("1>>n_bkg", "tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>={down}&&n_true_interactions<={up}&&tau_charge=={ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[1]))
            print "Number of signal events:     %i"%n_sig.GetBinContent(1)
            print "Number of background events: %i"%n_bkg.GetBinContent(1)
            
            #fill maps
            histsig = ROOT.TH2D("sig_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), "sig_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), nbins_pt, pt_binning, nbins_dr, 0.0, 0.5)
            histbkg = ROOT.TH2D("bkg_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), "bkg_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), nbins_pt, pt_binning, nbins_dr, 0.0, 0.5)
            tree.Draw("tau_iso_ph_dr:tau_iso_ph_pt>>sig_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), "(tau_gen_match==5&&n_true_interactions>={down}&&n_true_interactions<={up}&&tau_charge=={ch})".format(down=PUwin[0], up=PUwin[1], ch=charge[1]))
            tree.Draw("tau_iso_ph_dr:tau_iso_ph_pt>>bkg_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0]), "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>={down}&&n_true_interactions<={up}&&tau_charge=={ch})".format(down=PUwin[0], up=PUwin[1], ch=charge[1]))
            
            #normalize maps
            histsig.Scale(1.0/n_sig.GetBinContent(1))
            histbkg.Scale(1.0/n_bkg.GetBinContent(1))
            print "Expected integral number of charged tracks in the regarded phase space for real taus: %f"%histsig.Integral()
            print "Expected integral number of charged tracks in the regarded phase space for fake taus: %f"%histbkg.Integral()
            
            #save maps
            histsig.Write()
            histbkg.Write()
    file1.Close()

if create_discriminator:
    nevents_sig = [61016.0, 476480.0, 143537.0, 60690.0, 471014.0, 142300.0] #number of events per PU class (determined with first part of script)
    nevents_bkg = [46440.0, 368436.0, 113992.0, 43208.0, 345605.0, 107060.0]

    #load likelihood density
    file1 = ROOT.TFile("photon_maps.root")
    LDsig = []
    LDbkg = []
    for charge in ["P", "M"]:
        for PUwin in [[0, 1], [30, 35], [60, 70]]:
            LDsig.append(file1.Get("sig_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
            LDbkg.append(file1.Get("bkg_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
    
    #calculate the 0 track likelihoods
    print "Determine log likelihood difference between signal and fake for an event without tracks ..."
    Ldelta_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(1, nbins_pt):
        for j in range(1, nbins_dr):
            bin_i = LDsig[0].GetBin(i, j)
            for l in range(6):
                Ldelta_0[l] += math.log(1.0-LDbkg[l].GetBinContent(bin_i)) - math.log(1.0-LDsig[l].GetBinContent(bin_i))
    print "... Positive charge: PU 0 to 1: %f; PU 30 to 35: %f; PU 60 to 70 %f --- Negative charge: PU 0 to 1: %f; PU 30 to 35: %f; PU 60 to 70 %f"%tuple(Ldelta_0)
    #reset Ldelta_0
    for i in range(3):
        if Ldelta_0[i] > Ldelta_0[i+3]:
            Ldelta_0[i] = Ldelta_0[i] - Ldelta_0[i+3]
            Ldelta_0[i+3] = 0.0
        else:
            Ldelta_0[i+3] = Ldelta_0[i+3] - Ldelta_0[i]
            Ldelta_0[i] = 0.0

    #calculate iso and fill histograms
    file2 = ROOT.TFile("LLR_discrim_hists_neutraliso.root", "RECREATE")
    discrim_sig1 = ROOT.TH1I("LR_sig_ph_PU0to1", "LR_sig_ph_PU0to1", nbins_discrim, 0.0, upperlim)
    discrim_sig2 = ROOT.TH1I("LR_sig_ph_PU30to35", "LR_sig_ph_PU30to35", nbins_discrim, 0.0, upperlim)
    discrim_sig3 = ROOT.TH1I("LR_sig_ph_PU60to70", "LR_sig_ph_PU60to70", nbins_discrim, 0.0, upperlim)
    discrim_sig = [discrim_sig1, discrim_sig2, discrim_sig3]
    discrim_bkg1 = ROOT.TH1I("LR_bkg_ph_PU0to1", "LR_bkg_ph_PU0to1", nbins_discrim, 0.0, upperlim)
    discrim_bkg2 = ROOT.TH1I("LR_bkg_ph_PU30to35", "LR_bkg_ph_PU30to35", nbins_discrim, 0.0, upperlim)
    discrim_bkg3 = ROOT.TH1I("LR_bkg_ph_PU60to70", "LR_bkg_ph_PU60to70", nbins_discrim, 0.0, upperlim)
    discrim_bkg = [discrim_bkg1, discrim_bkg2, discrim_bkg3]

    nevents = tree.GetEntries()
    print "Process %i events..."%nevents
    i=0
    for event in tree:
        #if i==100000:
        #    break
        i+=1
        if (100.0*i/nevents)%10-(100.0*(i-1)/nevents)%10 < 0.0:
            print "Processed %i0%%"%int(10*i/nevents)
        
        is_true_tau = False
        pu_class = -1
        is_positive = True
        #classify true tau or fake
        if event.tau_gen_match==5:
            is_true_tau = True
        elif event.tau_gen_match==6 and event.gen_jet_pt>18 and abs(event.gen_jet_eta)<2.3:
            is_true_tau = False
        else:
            continue
        #classify PU
        if event.n_true_interactions>=0 and event.n_true_interactions<=1:
            pu_class = 0
        elif event.n_true_interactions>=30 and event.n_true_interactions<=35:
            pu_class = 1
        elif event.n_true_interactions>=60 and event.n_true_interactions<=70:
            pu_class = 2
        else:
            continue
        #classify charge
        if event.tau_charge==1.0:
            is_positive = True
            pu_class += 3 #evaluate with map created with opLLR_discrim_hists_neutralisopositely charged taus
        elif event.tau_charge==-1.0:
            is_positive = False
        else:
            continue
        #read charged tracks and determine bin numbers
        ph_pt = event.tau_iso_ph_pt
        ph_dr = event.tau_iso_ph_dr
        track_bins = []
        for track in zip(ph_pt, ph_dr):
            track_bins.append([LDsig[pu_class].GetXaxis().FindBin(track[0]), LDsig[pu_class].GetYaxis().FindBin(track[1])]) #[int(track[0]/(25.0/nbins))+1, int(track[1]/(0.5/nbins))+1])
        #calculate resulting log likelihood ratio
        llr = Ldelta_0[pu_class]
        for track in track_bins:
            exp_sig = LDsig[pu_class].GetBinContent(LDsig[pu_class].GetBin(*track))
            exp_bkg = LDbkg[pu_class].GetBinContent(LDbkg[pu_class].GetBin(*track))
            if exp_sig == 0.0 and exp_bkg == 0.0:
                continue
            elif exp_sig == 0.0:
                #llr = 999.0
                #break
                exp_sig = 1.0/nevents_sig[pu_class]
                #llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
            elif exp_bkg == 0.0:
                #llr = 0.0
                #break
                exp_bkg = 1.0/nevents_bkg[pu_class]    
            #else:
            llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
        if llr>upperlim*0.99:
            llr=upperlim*0.99
        if is_true_tau:
            discrim_sig[pu_class%3].Fill(llr)
        else:
            discrim_bkg[pu_class%3].Fill(llr)
    
    #create hists with standard isolation
    histsig1 = ROOT.TH1F("sig_ph_PU0to1", "sig_ph_PU0to1", nbins_discrim, 0.0, 20.0)
    histbkg1 = ROOT.TH1F("bkg_ph_PU0to1", "bkg_ph_PU0to1", nbins_discrim, 0.0, 20.0)
    histsig2 = ROOT.TH1F("sig_ph_PU30to35", "sig_ph_PU30to35", nbins_discrim, 0.0, 20.0)
    histbkg2 = ROOT.TH1F("bkg_ph_PU30to35", "bkg_ph_PU30to35", nbins_discrim, 0.0, 20.0)
    histsig3 = ROOT.TH1F("sig_ph_PU60to70", "sig_ph_PU60to70", nbins_discrim, 0.0, 20.0)
    histbkg3 = ROOT.TH1F("bkg_ph_PU60to70", "bkg_ph_PU60to70", nbins_discrim, 0.0, 20.0)
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>sig_ph_PU0to1", "(tau_gen_match==5&&n_true_interactions>=0&&n_true_interactions<=1)")
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>bkg_ph_PU0to1", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=0&&n_true_interactions<=1)")
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>sig_ph_PU30to35", "(tau_gen_match==5&&n_true_interactions>=30&&n_true_interactions<=35)")
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>bkg_ph_PU30to35", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=30&&n_true_interactions<=35)")
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>sig_ph_PU60to70", "(tau_gen_match==5&&n_true_interactions>=60&&n_true_interactions<=70)")
    tree.Draw("(tau_neutralIsoPtSum*(tau_neutralIsoPtSum<20.0)+19.99*(tau_neutralIsoPtSum>=20.0))>>bkg_ph_PU60to70", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=60&&n_true_interactions<=70)")
    histsig_std = [histsig1, histsig2, histsig3]
    histbkg_std = [histbkg1, histbkg2, histbkg3]
    
    for j in range(3):
        file2.cd()
        discrim_sig[j].Write()
        discrim_bkg[j].Write()
        histsig_std[j].Write()
        histbkg_std[j].Write()
        
    file2.Close()
    file1.Close()
file0.Close()

if create_plots:
    #load hists and create rebinned versions for plotting
    rebinning_factor=20
    PUs = ["PU0to1", "PU30to35", "PU60to70"]
    file2 = ROOT.TFile("LLR_discrim_hists_neutraliso.root")
    discrim_sig_load = []
    discrim_bkg_load = []
    stdiso_sig_load = []
    stdiso_bkg_load = []
    discrim_sig_rebinned = []
    discrim_bkg_rebinned = []
    for i in range(3):
        discrim_sig_load.append(file2.Get("LR_sig_ph_%s"%PUs[i]))
        discrim_sig_rebinned.append(ROOT.TH1F("LR_sig_%s_rebinned"%PUs[i], "LR_sig_ph_%s"%PUs[i], nbins_discrim/rebinning_factor, 0.0, upperlim))
        discrim_sig_rebinned[i].SetBinContent(0, discrim_sig_load[i].GetBinContent(0))
        discrim_sig_rebinned[i].SetBinContent(nbins_discrim/rebinning_factor+1, discrim_sig_load[i].GetBinContent(nbins_discrim+1))
        discrim_bkg_load.append(file2.Get("LR_bkg_ph_%s"%PUs[i]))
        discrim_bkg_rebinned.append(ROOT.TH1F("LR_bkg_%s_rebinned"%PUs[i], "LR_sig_ph_%s"%PUs[i], nbins_discrim/rebinning_factor, 0.0, upperlim))
        discrim_bkg_rebinned[i].SetBinContent(0, discrim_sig_load[i].GetBinContent(0))
        discrim_bkg_rebinned[i].SetBinContent(nbins_discrim/rebinning_factor+1, discrim_sig_load[i].GetBinContent(nbins_discrim+1))
        for j in range(1, nbins_discrim/rebinning_factor+1):
            content_sig = 0
            content_bkg = 0
            for k in range(rebinning_factor):
                content_sig += discrim_sig_load[i].GetBinContent(j*rebinning_factor-k)
                content_bkg += discrim_bkg_load[i].GetBinContent(j*rebinning_factor-k)
            discrim_sig_rebinned[i].SetBinContent(j, content_sig)
            discrim_bkg_rebinned[i].SetBinContent(j, content_bkg)
        stdiso_sig_load.append(file2.Get("sig_ph_%s"%PUs[i]))
        stdiso_bkg_load.append(file2.Get("bkg_ph_%s"%PUs[i]))
    
    #create plots
    colors = [1, 2, 4]
    discrim_sig_rebinned[0].GetXaxis().SetTitle("Relative log likelihood ratio")
    discrim_sig_rebinned[0].GetYaxis().SetTitle("Fraction of events")
    discrim_sig_rebinned[0].Draw("hist") #for the range and labels
    for i in range(3):
        discrim_sig_rebinned[i].Scale(1.0/discrim_sig_rebinned[i].Integral())
        discrim_bkg_rebinned[i].Scale(1.0/discrim_bkg_rebinned[i].Integral())
        discrim_sig_rebinned[i].SetLineColor(colors[i])
        discrim_bkg_rebinned[i].SetLineColor(colors[i])
        discrim_bkg_rebinned[i].SetLineStyle(2)
        discrim_sig_rebinned[i].SetLineWidth(3)
        discrim_bkg_rebinned[i].SetLineWidth(3)
        discrim_sig_rebinned[i].Draw("histsame")
        discrim_bkg_rebinned[i].Draw("histsame")
    legend1 = ROOT.TLegend(0.5, 0.76, 0.9, 0.9)
    legend1.AddEntry(discrim_sig_rebinned[0], "real taus", "l")
    legend1.AddEntry(discrim_bkg_rebinned[0], "fake taus", "l")
    legend2 = ROOT.TLegend(0.5, 0.5, 0.9, 0.71)
    legend2.AddEntry(discrim_sig_rebinned[0], "PU 0 to 1", "l")
    legend2.AddEntry(discrim_sig_rebinned[1], "PU 30 to 35", "l")
    legend2.AddEntry(discrim_sig_rebinned[2], "PU 60 to 70", "l")
    legend1.Draw()
    legend2.Draw()
    c1.SaveAs("ph_lldiscrim.png")
    
    rocs = []
    for i in range(3):
        rocs.append(histsToRoc(stdiso_sig_load[i], stdiso_bkg_load[i]))
        rocs[-1].name = "standard iso %s"%PUs[i]
        rocs[-1].title = "standard iso %s"%PUs[i]
    for i in range(3):
        rocs.append(histsToRoc(discrim_sig_load[i], discrim_bkg_load[i]))
        rocs[-1].name = "LLR iso %s"%PUs[i]
        rocs[-1].title = "LLR iso %s"%PUs[i]
    #makeROCPlot(rocs, "ph_LLR_ROCs", xmin=0.9, logy=True)
    makeROCPlot(rocs, "ph_LLR_ROCs_fullRange", xmin=0.9, ymin=0.5, logy=False)
    #makeROCPlot(rocs, "ph_LLR_ROCs_WPs", ymax=0.2, logy=False)

    file2.Close()

