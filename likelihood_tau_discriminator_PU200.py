import ROOT
import math
import copy
import numpy as np
from CMGTools.H2TauTau.proto.plotter.ROCPlotter import *

ROOT.gROOT.SetBatch(True) # don't disply canvas

create_maps = True
create_discriminator = True
create_plots = True

#settings for track maps
nbins_pt = 17
nbins_dz = 4
nbins_dr = 4
#pt_binning = np.array([0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2] + [1.0+39.0*0.7**x for x in reversed(range(10))])
#pt_binning = np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.4, 3.8, 4.4, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0])
pt_binning = np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.6, 2.0, 2.6, 3.4, 4.4, 7.0, 10.0, 15.0, 20.0, 30.0])
dz_binning = np.array([0.0, 0.05, 0.1, 0.2, 0.4])#([0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4])#([0.4/nbins_dz*x for x in range(nbins_dz+1)])
dr_binning = np.array([0.5/nbins_dr*x for x in range(nbins_dr+1)])

#settings for discriminator hists
upperlim = 100.0
nbins_discrim = 2000

file0 = ROOT.TFile("/eos/user/s/swozniew/TauPOG_files/ntuples/PU200/Histo_signal_pu200.root")
file00 = ROOT.TFile("/eos/user/s/swozniew/TauPOG_files/ntuples/PU200/Histo_QCD_pu200.root")
tree_sig = file0.Get("tree")
tree_bkg = file00.Get("tree")
c1 = ROOT.TCanvas()
sig_selection = "genTauPt>20&&genTauEta<3.0&&genTauEta>-3.0&&tauPt>20&&tauEta<3.0&&tauEta>-3.0&&taupfTausDiscriminationByDecayModeFinding==1&&genTauMatch==1"
bkg_selection = "jetPt>20&&jetEta<3.0&&jetEta>-3.0&&tauPt>20&&tauEta<3.0&&tauEta>-3.0&&taupfTausDiscriminationByDecayModeFinding==1&&genJetMatch==1"        
        

if create_maps:
    file1 = ROOT.TFile("track_maps_PU200.root", "RECREATE")
    for fold in [["EVEN", 0], ["ODD", 1]]:
        print "Creating maps for PU 200 and fold {fold} ...".format(fold=fold[0])
        
        #count event yields
        n_sig = ROOT.TH1I("n_sig", "n_sig", 1, 0, 2)
        tree_sig.Draw("1>>n_sig", "{sel}&&event%2=={fold}".format(sel=sig_selection, fold=fold[1]))
        n_bkg = ROOT.TH1I("n_bkg", "n_bkg", 1, 0, 2)
        tree_bkg.Draw("1>>n_bkg", "{sel}&&event%2=={fold}".format(sel=bkg_selection, fold=fold[1]))
        print "Number of signal events:     %i"%n_sig.GetBinContent(1)
        print "Number of background events: %i"%n_bkg.GetBinContent(1)
        
        #fill maps
        histsig = ROOT.TH3D("sig_{fold}".format(fold=fold[0]), "sig_{fold}".format(fold=fold[0]), nbins_pt, pt_binning, nbins_dz, dz_binning, nbins_dr, dr_binning)
        histbkg = ROOT.TH3D("bkg_{fold}".format(fold=fold[0]), "bkg_{fold}".format(fold=fold[0]), nbins_pt, pt_binning, nbins_dz, dz_binning, nbins_dr, dr_binning)
        tree_sig.Draw("DR:abs(dz):pt>>sig_{fold}".format(fold=fold[0]), "{sel}&&event%2=={fold}".format(sel=sig_selection, fold=fold[1]))
        tree_bkg.Draw("DR:abs(dz):NEW_pt>>bkg_{fold}".format(fold=fold[0]), "{sel}&&event%2=={fold}".format(sel=bkg_selection, fold=fold[1]))
        
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
    nevents_sig = [2665.0, 2629.0] #number of events per fold (determined with first part of script)
    nevents_bkg = [9039.0, 9076.0]

    #load likelihood density
    file1 = ROOT.TFile("track_maps_PU200.root")
    LDsig = []
    LDbkg = []
    for fold in ["EVEN", "ODD"]:
        LDsig.append(file1.Get("sig_{fold}".format(fold=fold)))
        LDbkg.append(file1.Get("bkg_{fold}".format(fold=fold)))
    
    #calculate the 0 track likelihoods
    print "Determine log likelihood difference between signal and fake for an event without tracks ..."
    Ldelta_0 = [0.0, 0.0]
    for i in range(1, LDsig[0].GetNbinsX()+1):
        for j in range(1, LDsig[0].GetNbinsY()+1):
            for k in range(1, LDsig[0].GetNbinsZ()+1):
                bin_i = LDsig[0].GetBin(i, j, k)
                for l in range(2):
                    Ldelta_0[l] += math.log(1.0-LDbkg[l].GetBinContent(bin_i)) - math.log(1.0-LDsig[l].GetBinContent(bin_i))
    print "... Even fold: %f --- Odd fold: %f"%tuple(Ldelta_0)
    #reset Ldelta_0
    if Ldelta_0[0] > Ldelta_0[1]:
        Ldelta_0[0] = Ldelta_0[0] - Ldelta_0[1]
        Ldelta_0[1] = 0.0
    else:
        Ldelta_0[1] = Ldelta_0[1] - Ldelta_0[0]
        Ldelta_0[0] = 0.0

    #calculate iso and fill histograms
    file2 = ROOT.TFile("LLR_discrim_hists_PU200.root", "RECREATE")
    discrim_sig = ROOT.TH1I("LR_sig_PU200", "LR_sig_PU200", nbins_discrim, 0.0, upperlim)
    discrim_bkg = ROOT.TH1I("LR_bkg_PU200", "LR_bkg_PU200", nbins_discrim, 0.0, upperlim)

    nevents = tree_sig.GetEntries()
    print "Process %i signal events..."%nevents
    i=0
    for event in tree_sig:
        #if i==100000:
        #    break
        i+=1
        if (100.0*i/nevents)%10-(100.0*(i-1)/nevents)%10 < 0.0:
            print "Processed %i0%%"%int(10*i/nevents)
        
        #filter signal
        if not (event.genTauPt > 20 \
                and event.genTauEta < 3.0 \
                and event.genTauEta > -3.0 \
                and event.tauPt > 20 \
                and event.tauEta < 3.0 \
                and event.tauEta > -3.0 \
                and event.taupfTausDiscriminationByDecayModeFinding == "\x01" \
                and event.genTauMatch == 1):
            continue
        #classify fold
        fold = 1 - (event.event%2) #evaluate on opposite fold
        #read charged tracks and determine bin numbers
        ch_pt = event.pt
        ch_dz = event.dz
        ch_dr = event.DR
        track_bins = []
        for track in zip(ch_pt, ch_dz, ch_dr):
            track_bins.append([LDsig[fold].GetXaxis().FindBin(track[0]), LDsig[fold].GetYaxis().FindBin(abs(track[1])), LDsig[fold].GetZaxis().FindBin(track[2])]) #[int(track[0]/(20.0/nbins))+1, int(abs(track[1])/(0.4/nbins))+1, int(track[2]/(0.5/nbins))+1])
        #calculate resulting log likelihood ratio
        llr = Ldelta_0[fold]
        for track in track_bins:
            exp_sig = LDsig[fold].GetBinContent(LDsig[fold].GetBin(*track))
            exp_bkg = LDbkg[fold].GetBinContent(LDbkg[fold].GetBin(*track))
            if exp_sig == 0.0 and exp_bkg == 0.0:
                continue
            elif exp_sig == 0.0:
                #llr = 999.0
                #break
                exp_sig = 1.0/nevents_sig[fold]
                #llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
            elif exp_bkg == 0.0:
                #llr = 0.0
                #break
                exp_bkg = 1.0/nevents_bkg[fold]    
            #else:
            llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
        if llr>upperlim*0.99:
            llr=upperlim*0.99
        if llr<0.0:
            #for track in zip(ch_pt, ch_dz, ch_dr):
            #    print track
            llr=0.0
        discrim_sig.Fill(llr)
    nevents = tree_bkg.GetEntries()
    print "Process %i background events..."%nevents
    i=0
    for event in tree_bkg:
        #if i==100000:
        #    break
        i+=1
        if (100.0*i/nevents)%10-(100.0*(i-1)/nevents)%10 < 0.0:
            print "Processed %i0%%"%int(10*i/nevents)
        
        #filter signal
        if not (event.jetPt > 20 \
                and event.jetEta < 3.0 \
                and event.jetEta > -3.0 \
                and event.tauPt > 20 \
                and event.tauEta < 3.0 \
                and event.tauEta > -3.0 \
                and event.taupfTausDiscriminationByDecayModeFinding == "\x01" \
                and event.genJetMatch == 1):
            continue
        #classify fold
        fold = 1 - (event.event%2) #evaluate on opposite fold
        #read charged tracks and determine bin numbers
        ch_pt = event.NEW_pt
        ch_dz = event.dz
        ch_dr = event.DR
        track_bins = []
        for track in zip(ch_pt, ch_dz, ch_dr):
            track_bins.append([LDsig[fold].GetXaxis().FindBin(track[0]), LDsig[fold].GetYaxis().FindBin(abs(track[1])), LDsig[fold].GetZaxis().FindBin(track[2])]) #[int(track[0]/(20.0/nbins))+1, int(abs(track[1])/(0.4/nbins))+1, int(track[2]/(0.5/nbins))+1])
        #calculate resulting log likelihood ratio
        llr = Ldelta_0[fold]
        for track in track_bins:
            exp_sig = LDsig[fold].GetBinContent(LDsig[fold].GetBin(*track))
            exp_bkg = LDbkg[fold].GetBinContent(LDbkg[fold].GetBin(*track))
            if exp_sig == 0.0 and exp_bkg == 0.0:
                continue
            elif exp_sig == 0.0:
                #llr = 999.0
                #break
                exp_sig = 1.0/nevents_sig[fold]
                #llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
            elif exp_bkg == 0.0:
                #llr = 0.0
                #break
                exp_bkg = 1.0/nevents_bkg[fold]    
            #else:
            llr += math.log(exp_bkg) - math.log(1.0-exp_bkg) - math.log(exp_sig) + math.log(1.0-exp_sig)
        #print llr
        if llr>upperlim*0.99:
            llr=upperlim*0.99
        if llr<0.0:
            llr=0.0
        discrim_bkg.Fill(llr)
    
    #create hists with standard isolation
    histsig = ROOT.TH1F("sig_PU200", "sig_PU200", nbins_discrim, 0.0, 100.0)
    histbkg = ROOT.TH1F("bkg_PU200", "bkg_PU200", nbins_discrim, 0.0, 100.0)
    histsig_w = ROOT.TH1F("sig_PU200_w", "sig_PU200_w", nbins_discrim, 0.0, 100.0)
    histbkg_w = ROOT.TH1F("bkg_PU200_w", "bkg_PU200_w", nbins_discrim, 0.0, 100.0)
    tree_sig.Draw("(tauChargedIsoPtSum*(tauChargedIsoPtSum<100.0)+99.99*(tauChargedIsoPtSum>=100.0))>>sig_PU200", sig_selection)
    tree_bkg.Draw("(tauChargedIsoPtSum*(tauChargedIsoPtSum<100.0)+99.99*(tauChargedIsoPtSum>=100.0))>>bkg_PU200", bkg_selection)
    tree_sig.Draw("(cand_Iso_Weight1*(cand_Iso_Weight1<100.0)+99.99*(cand_Iso_Weight1>=100.0))>>sig_PU200_w", sig_selection)
    tree_bkg.Draw("(cand_Iso_Weight1*(cand_Iso_Weight1<100.0)+99.99*(cand_Iso_Weight1>=100.0))>>bkg_PU200_w", bkg_selection)
    
    file2.cd()
    discrim_sig.Write()
    discrim_bkg.Write()
    histsig.Write()
    histbkg.Write()
    histsig_w.Write()
    histbkg_w.Write()
    
    file2.Close()
    file1.Close()
file0.Close()

if create_plots:
    #load hists and create rebinned versions for plotting
    rebinning_factor=5
    PUs = ["PU200"]
    file2 = ROOT.TFile("LLR_discrim_hists_PU200.root")
    discrim_sig_load = []
    discrim_bkg_load = []
    stdiso_sig_load = []
    stdiso_bkg_load = []
    stdiso_sig_w_load = []
    stdiso_bkg_w_load = []
    discrim_sig_rebinned = []
    discrim_bkg_rebinned = []
    for i in range(1):
        discrim_sig_load.append(file2.Get("LR_sig_%s"%PUs[i]))
        discrim_sig_rebinned.append(ROOT.TH1F("LR_sig_%s_rebinned"%PUs[i], "LR_sig_%s"%PUs[i], nbins_discrim/rebinning_factor, 0.0, upperlim))
        discrim_sig_rebinned[i].SetBinContent(0, discrim_sig_load[i].GetBinContent(0))
        discrim_sig_rebinned[i].SetBinContent(nbins_discrim/rebinning_factor+1, discrim_sig_load[i].GetBinContent(nbins_discrim+1))
        discrim_bkg_load.append(file2.Get("LR_bkg_%s"%PUs[i]))
        discrim_bkg_rebinned.append(ROOT.TH1F("LR_bkg_%s_rebinned"%PUs[i], "LR_sig_%s"%PUs[i], nbins_discrim/rebinning_factor, 0.0, upperlim))
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
        stdiso_sig_load.append(file2.Get("sig_%s"%PUs[i]))
        stdiso_bkg_load.append(file2.Get("bkg_%s"%PUs[i]))
        stdiso_sig_w_load.append(file2.Get("sig_%s_w"%PUs[i]))
        stdiso_bkg_w_load.append(file2.Get("bkg_%s_w"%PUs[i]))
    
    #create plots
    colors = [1, 2, 4]
    discrim_sig_rebinned[0].GetXaxis().SetTitle("Relative log likelihood ratio")
    discrim_sig_rebinned[0].GetYaxis().SetTitle("Fraction of events")
    discrim_sig_rebinned[0].Draw("hist") #for the range and labels
    for i in range(1):
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
    #legend2 = ROOT.TLegend(0.5, 0.5, 0.9, 0.71)
    #legend2.AddEntry(discrim_sig_rebinned[0], "PU 0 to 1", "l")
    #legend2.AddEntry(discrim_sig_rebinned[1], "PU 30 to 35", "l")
    #legend2.AddEntry(discrim_sig_rebinned[2], "PU 60 to 70", "l")
    legend1.Draw()
    #legend2.Draw()
    c1.SaveAs("PU200_lldiscrim.png")
    
    rocs = []
    for i in range(1):
        rocs.append(histsToRoc(stdiso_sig_load[i], stdiso_bkg_load[i]))
        rocs[-1].name = "standard iso %s"%PUs[i]
        rocs[-1].title = "standard iso %s"%PUs[i]
    for i in range(1):
        rocs.append(histsToRoc(stdiso_sig_w_load[i], stdiso_bkg_w_load[i]))
        rocs[-1].name = "reweighted iso %s"%PUs[i]
        rocs[-1].title = "reweighted iso %s"%PUs[i]
    for i in range(1):
        rocs.append(histsToRoc(discrim_sig_load[i], discrim_bkg_load[i]))
        rocs[-1].name = "LLR iso %s"%PUs[i]
        rocs[-1].title = "LLR iso %s"%PUs[i]
    #makeROCPlot(rocs, "LLR_ROCs", xmin=0.9, logy=True)
    #makeROCPlot(rocs, "LLR_ROCs_fullRange", xmin=0.0, logy=True)
    makeROCPlot(rocs, "PU200_LLR_ROCs_WPs", ymax=0.04, logy=False)

    file2.Close()

