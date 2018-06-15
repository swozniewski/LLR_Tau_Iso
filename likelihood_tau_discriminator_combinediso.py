import ROOT
import math
import copy
from CMGTools.H2TauTau.proto.plotter.ROCPlotter import *

ROOT.gROOT.SetBatch(True) # don't disply canvas

create_discriminator = True
create_plots = True

#settings for discriminator hists
upperlim = 40.0
nbins_discrim = 200

file0 = ROOT.TFile("/eos/user/s/swozniew/TauPOG_files/ntuples/DYJetsToLL_M50_LO.root")
tree = file0.Get("tree")
c1 = ROOT.TCanvas()

if create_discriminator:
    nevents_sig = [61016.0, 476480.0, 143537.0, 60690.0, 471014.0, 142300.0] #number of events per PU class (determined with first part of script)
    nevents_bkg = [46440.0, 368436.0, 113992.0, 43208.0, 345605.0, 107060.0]

    #load likelihood density
    file1 = ROOT.TFile("combined_maps.root")
    LDsig_trk = []
    LDbkg_trk = []
    LDsig_ph = []
    LDbkg_ph = []
    for charge in ["P", "M"]:
        for PUwin in [[0, 1], [30, 35], [60, 70]]:
            LDsig_trk.append(file1.Get("sig_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
            LDbkg_trk.append(file1.Get("bkg_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
            LDsig_ph.append(file1.Get("sig_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
            LDbkg_ph.append(file1.Get("bkg_ph_PU{down}to{up}_{ch}".format(down=PUwin[0], up=PUwin[1], ch=charge[0])))
    
    #calculate the 0 track likelihoods
    print "Determine log likelihood difference between signal and fake for an event without tracks ..."
    Ldelta_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(1, LDsig_trk[0].GetNbinsX()+1):
        for j in range(1, LDsig_trk[0].GetNbinsY()+1):
            for k in range(1, LDsig_trk[0].GetNbinsZ()+1):
                bin_i = LDsig_trk[0].GetBin(i, j, k)
                for l in range(6):
                    Ldelta_0[l] += math.log(1.0-LDbkg_trk[l].GetBinContent(bin_i)) - math.log(1.0-LDsig_trk[l].GetBinContent(bin_i))
    for i in range(1, LDsig_ph[0].GetNbinsX()+1):
        for j in range(1, LDsig_ph[0].GetNbinsY()+1):
            bin_i = LDsig_ph[0].GetBin(i, j)
            for l in range(6):
                Ldelta_0[l] += math.log(1.0-LDbkg_ph[l].GetBinContent(bin_i)) - math.log(1.0-LDsig_ph[l].GetBinContent(bin_i))
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
    file2 = ROOT.TFile("LLR_discrim_hists_combined.root", "RECREATE")
    discrim_sig1 = ROOT.TH1I("LR_sig_PU0to1", "LR_sig_PU0to1", nbins_discrim, 0.0, upperlim)
    discrim_sig2 = ROOT.TH1I("LR_sig_PU30to35", "LR_sig_PU30to35", nbins_discrim, 0.0, upperlim)
    discrim_sig3 = ROOT.TH1I("LR_sig_PU60to70", "LR_sig_PU60to70", nbins_discrim, 0.0, upperlim)
    discrim_sig = [discrim_sig1, discrim_sig2, discrim_sig3]
    discrim_bkg1 = ROOT.TH1I("LR_bkg_PU0to1", "LR_bkg_PU0to1", nbins_discrim, 0.0, upperlim)
    discrim_bkg2 = ROOT.TH1I("LR_bkg_PU30to35", "LR_bkg_PU30to35", nbins_discrim, 0.0, upperlim)
    discrim_bkg3 = ROOT.TH1I("LR_bkg_PU60to70", "LR_bkg_PU60to70", nbins_discrim, 0.0, upperlim)
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
            pu_class += 3 #evaluate with map created with oppositely charged taus
        elif event.tau_charge==-1.0:
            is_positive = False
        else:
            continue
        #read charged tracks and determine bin numbers
        ch_pt = event.tau_iso_ch_pt
        ch_dz = event.tau_iso_ch_dz
        ch_dr = event.tau_iso_ch_dr
        ph_pt = event.tau_iso_ph_pt
        ph_dr = event.tau_iso_ph_dr
        track_bins = []
        photon_bins = []
        for track in zip(ch_pt, ch_dz, ch_dr):
            track_bins.append([LDsig_trk[pu_class].GetXaxis().FindBin(track[0]), LDsig_trk[pu_class].GetYaxis().FindBin(abs(track[1])), LDsig_trk[pu_class].GetZaxis().FindBin(track[2])]) #[int(track[0]/(20.0/nbins))+1, int(abs(track[1])/(0.4/nbins))+1, int(track[2]/(0.5/nbins))+1])
        for photon in zip(ph_pt, ph_dr):
            photon_bins.append([LDsig_ph[pu_class].GetXaxis().FindBin(photon[0]), LDsig_ph[pu_class].GetYaxis().FindBin(photon[1])])
        #calculate resulting log likelihood ratio
        llr = Ldelta_0[pu_class]
        for track in track_bins:
            exp_sig = LDsig_trk[pu_class].GetBinContent(LDsig_trk[pu_class].GetBin(*track))
            exp_bkg = LDbkg_trk[pu_class].GetBinContent(LDbkg_trk[pu_class].GetBin(*track))
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
        for photon in photon_bins:
            exp_sig = LDsig_ph[pu_class].GetBinContent(LDsig_ph[pu_class].GetBin(*photon))
            exp_bkg = LDbkg_ph[pu_class].GetBinContent(LDbkg_ph[pu_class].GetBin(*photon))
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
    histsig1 = ROOT.TH1F("sig_PU0to1", "sig_PU0to1", nbins_discrim, 0.0, 40.0)
    histbkg1 = ROOT.TH1F("bkg_PU0to1", "bkg_PU0to1", nbins_discrim, 0.0, 40.0)
    histsig2 = ROOT.TH1F("sig_PU30to35", "sig_PU30to35", nbins_discrim, 0.0, 40.0)
    histbkg2 = ROOT.TH1F("bkg_PU30to35", "bkg_PU30to35", nbins_discrim, 0.0, 40.0)
    histsig3 = ROOT.TH1F("sig_PU60to70", "sig_PU60to70", nbins_discrim, 0.0, 40.0)
    histbkg3 = ROOT.TH1F("bkg_PU60to70", "bkg_PU60to70", nbins_discrim, 0.0, 40.0)
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>sig_PU0to1", "(tau_gen_match==5&&n_true_interactions>=0&&n_true_interactions<=1)")
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>bkg_PU0to1", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=0&&n_true_interactions<=1)")
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>sig_PU30to35", "(tau_gen_match==5&&n_true_interactions>=30&&n_true_interactions<=35)")
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>bkg_PU30to35", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=30&&n_true_interactions<=35)")
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>sig_PU60to70", "(tau_gen_match==5&&n_true_interactions>=60&&n_true_interactions<=70)")
    tree.Draw("(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits<20.0)+19.99*(tau_byCombinedIsolationDeltaBetaCorrRaw3Hits>=20.0))>>bkg_PU60to70", "(tau_gen_match==6&&gen_jet_pt>18&&abs(gen_jet_eta)<2.3&&n_true_interactions>=60&&n_true_interactions<=70)")
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
    rebinning_factor=5
    PUs = ["PU0to1", "PU30to35", "PU60to70"]
    file2 = ROOT.TFile("LLR_discrim_hists_combined.root")
    discrim_sig_load = []
    discrim_bkg_load = []
    stdiso_sig_load = []
    stdiso_bkg_load = []
    discrim_sig_rebinned = []
    discrim_bkg_rebinned = []
    for i in range(3):
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
    c1.SaveAs("combined_lldiscrim.png")
    
    rocs = []
    for i in range(3):
        rocs.append(histsToRoc(stdiso_sig_load[i], stdiso_bkg_load[i]))
        rocs[-1].name = "standard iso %s"%PUs[i]
        rocs[-1].title = "standard iso %s"%PUs[i]
    for i in range(3):
        rocs.append(histsToRoc(discrim_sig_load[i], discrim_bkg_load[i]))
        rocs[-1].name = "LLR iso %s"%PUs[i]
        rocs[-1].title = "LLR iso %s"%PUs[i]
    #makeROCPlot(rocs, "combined_LLR_ROCs", xmin=0.9, logy=True)
    #makeROCPlot(rocs, "combined_LLR_ROCs_fullRange", xmin=0.0, logy=True)
    makeROCPlot(rocs, "combined_LLR_ROCs_WPs", ymax=0.2, logy=False)

    file2.Close()

