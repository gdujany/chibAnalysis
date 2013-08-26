#!/usr/bin/python
from ROOT import gSystem, gROOT,  gStyle, TFile, TCanvas, TH1, TH1D, TTree, THStack, TLine, TPad, TPaveText
import ROOT, sys
sys.path.append('..')
sys.path.append('../Polarization')
from pyUtils import *
from polarization import *

TH1.SetDefaultSumw2(True) #Cosi` in automatico per tutti gli isto dovrebbe abilitare il calcolo degli errori come si deve (come somma dei pesi)



ptFileLabel = 'Upsilon2SPt' # 'Upsilon2SPt', 'Upsilon2SPt', 'Upsilon3SPt'
ptFileLabel1 = ptFileLabel#'UpsilonPt'
ptFileLabel2 = ptFileLabel#'Upsilon3SPt'

mc_ptSpectrum = '3S' # '1S', '2S', '3S', 'flat' 


output_name = 'agreementDataMC_'+mc_ptSpectrum#+'_rw'
outputFile_name = output_name+'.ps'

 
def makeHistos():

    cuts = Cuts()
    cuts = cuts.loadFromFile("../cuts.txt")
    # cuts.photon_pt_cut = 0.7

    
    rootupla_Chib1_Sel = '../../store/ChiB_1P_1_'+ptFileLabel1+'.root'
    rootupla_Chib2_Sel = '../../store/ChiB_1P_2_'+ptFileLabel2+'.root'
    rootupla_2012Data  = '../../store/2012_AllData.root'

    # Histograms
    histos_SR1 = dict()
    histos_SR2 = dict()
    histos_LSB = dict()
    histos_RSB = dict()
    histos_MC = dict()

    end = 40

    histos_SR1['photon_pt'] = TH1D("photon_pt_SR1","#chi_{b1};E_{T} [GeV];",99,0.2,3.5)
    histos_SR1['photon_eta'] = TH1D("photon_eta_SR1", "#chi_{b1};#eta;",99,-1.,1.)
    histos_SR1['dimuon_pt'] = TH1D("dimuon_pt_SR1","#chi_{b1};p_{T} [GeV];",99,5.,end)
    histos_SR1['dimuon_rapidity'] = TH1D("dimuon_rapidity_SR1","#chi_{b1};y;",99,-1.3,1.3) 
    histos_SR1['chib_pt'] = TH1D("chib_pt_SR1","#chi_{b1};p_{T} [GeV];",99,5.,end)

    histos_SR2['photon_pt'] = TH1D("photon_pt_SR2","#chi_{b2};E_{T} [GeV];",99,0.2,3.5)
    histos_SR2['photon_eta'] = TH1D("photon_eta_SR2", "#chi_{b2};#eta;",99,-1.,1.)
    histos_SR2['dimuon_pt'] = TH1D("dimuon_pt_SR2","#chi_{b2};p_{T} [GeV];",99,5.,end)
    histos_SR2['dimuon_rapidity'] = TH1D("dimuon_rapidity_SR2","#chi_{b2};y;",99,-1.3,1.3) 
    histos_SR2['chib_pt'] = TH1D("chib_pt_SR2","#chi_{b2};p_{T} [GeV];",99,5.,end)
        
    histos_LSB['photon_pt'] = TH1D("photon_pt_LSB","LSB;E_{T} [GeV];",99,0.2,3.5)
    histos_LSB['photon_eta'] = TH1D("photon_eta_LSB", "LSB;#eta;",99,-1.,1.)
    histos_LSB['dimuon_pt'] = TH1D("dimuon_pt_LSB","LSB;p_{T} [GeV];",99,5.,end)
    histos_LSB['dimuon_rapidity'] = TH1D("dimuon_rapidity_LSB","LSB;y;",99,-1.3,1.3)
    histos_LSB['chib_pt'] = TH1D("chib_pt_LSB","LSB;p_{T} [GeV];",99,5.,end)
    
    histos_RSB['photon_pt'] = TH1D("photon_pt_RSB","RSB;E_{T} [GeV];",99,0.2,3.5)
    histos_RSB['photon_eta'] = TH1D("photon_eta_RSB", "RSB;#eta;",99,-1.,1.)
    histos_RSB['dimuon_pt'] = TH1D("dimuon_pt_RSB","RSB;p_{T} [GeV];",99,5.,end)
    histos_RSB['dimuon_rapidity'] = TH1D("dimuon_rapidity_RSB","RSB;y;",99,-1.3,1.3)
    histos_RSB['chib_pt'] = TH1D("chib_pt_RSB","RSB;p_{T} [GeV];",99,5.,end)

    histos_MC['photon_pt'] = TH1D("photon_pt_MC","MC #chi_{b1,2};E_{T} [GeV];",99,0.2,3.5)
    histos_MC['photon_eta'] = TH1D("photon_eta_MC", "MC #chi_{b1,2};#eta;",99,-1.,1.)
    histos_MC['dimuon_pt'] = TH1D("dimuon_pt_MC","MC #chi_{b1,2};p_{T} [GeV];",99,5.,end)
    histos_MC['dimuon_rapidity'] = TH1D("dimuon_rapidity_MC","MC #chi_{b1,2};y;",99,-1.3,1.3)
    histos_MC['chib_pt'] = TH1D("chib_pt_MC","MC #chi_{b1,2};p_{T} [GeV];",99,5.,end)


    tree = TTree()
    chibTreeName = 'rootuple/chibTree'

    def fillMCHistos(inputFileName, num_chib, histos, weight = 1):
        inputFile = TFile.Open(inputFileName,"READ")
        inputFile.GetObject(chibTreeName, tree)
        for event in tree:
            if(event.Y1S_nsigma < 2.5 and event.probFit1S > 0.02 and event.dz < cuts.dz_cut and event.dimuon_pt > cuts.upsilon_pt_lcut and event.dimuon_pt < cuts.upsilon_pt_hcut and event.photon_pt > cuts.photon_pt_cut and abs(event.photon_eta) < cuts.photon_eta_cut and  abs(event.dimuon_rapidity) < cuts.upsilon_rapidity_cut) and abs(event.muonP_eta) < cuts.muon_eta_cut and abs(event.muonM_eta) < cuts.muon_eta_cut and event.muonP_pt > cuts.muon_pt_cut and event.muonM_pt > cuts.muon_pt_cut: 
                PtWeight = givePtweight(pt=event.chib_p4.Pt(), outPt=mc_ptSpectrum, inPt='2S')
                ups_dir, mu_dir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4,'hx')
                weight_polarization = angDist(ups_dir, mu_dir, num_chib, None)
                for key, histo in histos.items():
                    histo.Fill(getattr(event, key), weight * PtWeight * weight_polarization)
        inputFile.Close()

    print 'Filling MC histograms'
    fillMCHistos(rootupla_Chib1_Sel, 1, histos_MC, 1)
    fillMCHistos(rootupla_Chib2_Sel, 2, histos_MC, 0.61)
   
    print 'Filling Data histograms'
    inputFile = TFile.Open(rootupla_2012Data,"READ")
    inputFile.GetObject(chibTreeName, tree)
    for event in tree:
        for key in histos_SR1.keys():
            if(event.Y1S_nsigma < 2.5 and event.probFit1S > 0.02 and event.dz < cuts.dz_cut and event.dimuon_pt > cuts.upsilon_pt_lcut and event.dimuon_pt < cuts.upsilon_pt_hcut and event.photon_pt > cuts.photon_pt_cut and abs(event.photon_eta) < cuts.photon_eta_cut and  abs(event.dimuon_rapidity) < cuts.upsilon_rapidity_cut) and abs(event.muonP_eta) < cuts.muon_eta_cut and abs(event.muonM_eta) < cuts.muon_eta_cut and event.muonP_pt > cuts.muon_pt_cut and event.muonM_pt > cuts.muon_pt_cut:
                if event.rf1S_chib_mass > 9.83 and event.rf1S_chib_mass < 9.9 :
                    histos_SR1[key].Fill(getattr(event, key))
                if event.rf1S_chib_mass > 9.9 and event.rf1S_chib_mass < 9.93 :
                    histos_SR2[key].Fill(getattr(event, key))
                elif event.rf1S_chib_mass > 9.7 and event.rf1S_chib_mass < 9.83:
                    histos_LSB[key].Fill(getattr(event, key))
                elif event.rf1S_chib_mass > 9.93 and event.rf1S_chib_mass < 10.1:
                    histos_RSB[key].Fill(getattr(event, key))
    inputFile.Close()

        
    outputFile = TFile(output_name+'.root','RECREATE')
    for histo in histos_SR1.values() + histos_SR2.values() + histos_MC.values() + histos_LSB.values() + histos_RSB.values():
        histo.Write()
    outputFile.Close()

def makePlots():

    # gROOT.SetStyle('Plain') # white background
    # gStyle.SetOptStat(222211)

    ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style

    inputFile = TFile(output_name+'.root','READ')

    histos_SR1 = dict()
    histos_SR2 = dict()
    histos_LSB = dict()
    histos_RSB = dict()
    histos_MC = dict()

    keys = ['photon_pt', 'photon_eta', 'dimuon_pt', 'dimuon_rapidity', 'chib_pt'] 

    for key in keys:
        histos_SR1[key] =  inputFile.Get(key+'_SR1')
        histos_SR2[key] =  inputFile.Get(key+'_SR2')
        histos_LSB[key] =  inputFile.Get(key+'_LSB')
        histos_RSB[key] = inputFile.Get(key+'_RSB')
        histos_MC[key] =  inputFile.Get(key+'_MC')

    titles = dict(
        photon_pt = 'photon E_{T};E_{T}(#gamma) [GeV];Entries / 100 MeV',
        photon_eta = 'photon_eta;#eta(#gamma);Entries / 0.05',
        dimuon_pt = '#Upsilon p_{T};p_{T}(#Upsilon) [GeV];Entries / GeV',
        dimuon_rapidity = '#Upsilon rapidity(#Upsilon);y(#Upsilon);Entries / 0.1',
        chib_pt = '#chi_{b} p_{T};p_{T}(#chi_{b}) [GeV];Entries / GeV',
        )

    titles_chib = dict(
        photon_pt = 'photon E_{T};E_{T}(#gamma) [GeV];A.U.',
        photon_eta = 'photon_eta;#eta(#gamma);A.U.',
        dimuon_pt = '#Upsilon p_{T};p_{T}(#Upsilon) [GeV];A.U.',
        dimuon_rapidity = '#Upsilon rapidity(#Upsilon);y(#Upsilon);A.U.',
        chib_pt = '#chi_{b} p_{T};p_{T}(#chi_{b}) [GeV];A.U.',
        )
    
    logY = dict(
        photon_pt = 0,
        photon_eta = 0,
        dimuon_pt = 0,
        dimuon_rapidity = 0,
        chib_pt = 0,
        )

    maxY = dict(
        photon_pt = 750,
        photon_eta = 350,
        dimuon_pt = 800,
        dimuon_rapidity = 400,
        chib_pt = 800,
        )

    integralData1 = histos_SR1[keys[0]].Integral()
    integralData2 = histos_SR2[keys[0]].Integral()
    integralData = integralData1 + integralData2
  
    for histo in histos_LSB.values() + histos_RSB.values()  + histos_MC.values():
        if(histo.Integral() != 0):
            histo.Scale(1./histo.Integral()) 
            histo.Rebin(3)
    for histo in histos_SR1.values() + histos_SR2.values():
        histo.Rebin(3)
    
    fBG = 0.48 #0.28
    fLSB = 0.5

    histos_BG = dict()
    for key in keys:
        histos_BG[key] = histos_RSB[key].Clone()
        histos_BG[key].SetNameTitle(key+'_BG', 'bgd from side-bands')
        histos_BG[key].Scale(1-fLSB)
        histos_BG[key].Add(histos_LSB[key], fLSB)

    histos_Signal1 = dict()
    histos_Signal2 = dict()
    histos_Signal = dict()
    histos_SR = dict()
    for key in keys:
        histos_Signal1[key] = histos_SR1[key].Clone()
        histos_Signal1[key].SetNameTitle(key+'_signal1', 'Signal1')
        histos_Signal1[key].Add(histos_BG[key], -1*integralData1*fBG)
        
        histos_Signal2[key] = histos_SR2[key].Clone()
        histos_Signal2[key].SetNameTitle(key+'_signal2', 'Signal2')
        histos_Signal2[key].Add(histos_BG[key], -1*integralData2*fBG)
        
        histos_Signal[key] = histos_Signal1[key].Clone()
        histos_Signal[key].Add(histos_Signal2[key])
        histos_Signal[key].SetNameTitle(key+'_signal','data bkg-subtracted')

        histos_SR[key] = histos_SR1[key].Clone()
        histos_SR[key].Add(histos_SR2[key])
        histos_SR[key].SetNameTitle(key+'_SR','data')
        
    for histo in histos_Signal.values() + histos_BG.values():
        if(histo.Integral() != 0):
            histo.Scale(1./histo.Integral()) 

    # Normalize SR to their original value, MC to be comparable with signal background-subtracted and background such that Signal+background = SR
    for histo in histos_MC.values():
        histo.Scale(integralData * (1-fBG))
    for histo in histos_BG.values():
        histo.Scale(integralData * fBG)
    for histo in histos_LSB.values():
        histo.Scale(integralData * fBG * fLSB)
    for histo in histos_RSB.values():
        histo.Scale(integralData * fBG * (1-fLSB))
    
    chi2Prob = dict()
    kolmogorovProb = dict()
    h_ratio = dict()
    for key in histos_Signal.keys():
        h_ratio[key] = histos_SR[key].Clone()
        histo_sum_MC_BG = histos_MC[key].Clone()
        histo_sum_MC_BG.Add(histos_BG[key])
        h_ratio[key].Divide(histo_sum_MC_BG)
        h_ratio[key].SetNameTitle(key+'_ratio',titles[key])
        h_ratio[key].SetMarkerStyle(20)
        h_ratio[key].SetMarkerSize(1)
        h_ratio[key].GetYaxis().SetLabelSize(0.18)
        h_ratio[key].GetYaxis().SetTitle('DATA/MC')
        h_ratio[key].GetYaxis().SetTitleSize(0.22)
        h_ratio[key].GetYaxis().SetTitleOffset(0.2)
        h_ratio[key].GetYaxis().SetNdivisions(005)
        h_ratio[key].GetXaxis().SetLabelSize(0.16)
        h_ratio[key].GetXaxis().SetTitleSize(0.16)
        h_ratio[key].SetMaximum(1.5)
        h_ratio[key].SetMinimum(0.667)
        chi2Prob[key] = histos_SR[key].Chi2Test(histo_sum_MC_BG)#, 'UW')#, 'NORM')
        kolmogorovProb[key] = histos_SR[key].KolmogorovTest(histo_sum_MC_BG)

    for histo in histos_SR.values():
        histo.SetLineColor(ROOT.kBlack)
        histo.SetLineWidth(2)
        histo.SetMarkerColor(ROOT.kBlack)
        histo.SetMarkerStyle(20)

    for histo in histos_SR1.values()+histos_Signal1.values():
        histo.SetLineColor(ROOT.kBlue)
        histo.SetLineWidth(2)
        histo.SetMarkerColor(ROOT.kBlue)
        histo.SetMarkerStyle(22)

    for histo in histos_SR2.values()+histos_Signal2.values():
        histo.SetLineColor(ROOT.kRed)
        histo.SetLineWidth(2)
        histo.SetMarkerColor(ROOT.kRed)
        histo.SetMarkerStyle(23)

    for histo in histos_Signal.values():
        histo.SetLineColor(ROOT.kBlack)
        histo.SetLineWidth(2)
        # histo.SetMarkerColor(ROOT.kBlack)
        # histo.SetMarkerStyle(21)

    for histo in histos_BG.values():
        histo.SetLineColor(ROOT.kRed)
        histo.SetLineWidth(2)
        histo.SetFillColor(ROOT.kRed)
        # histo.SetMarkerColor(ROOT.kRed)
        # histo.SetMarkerStyle(22)

    # for histo in histos_RSB.values():
    #     histo.SetLineColor(ROOT.kViolet)
    #     histo.SetLineWidth(2)  
        
                
    for histo in histos_MC.values():
        histo.SetLineColor(ROOT.kBlue)
        histo.SetLineWidth(2)
        histo.SetFillColor(ROOT.kBlue)
        # histo.SetMarkerColor(ROOT.kBlue)
        # histo.SetMarkerStyle(20)

    
    outputFile = TFile(output_name+'_plots.root','RECREATE')
    for histo in histos_Signal.values() + histos_MC.values() + histos_BG.values() + histos_SR.values():
        histo.Write()
    outputFile.Close()
        
    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+"[") #apre file .ps

    for key in keys:
        pad1 = TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0)
        pad1.cd()
        pad1.SetLogy(logY[key])
        hs = THStack('hs', titles[key])
        #hs.Add(histos_SR[key])
        #hs.Add(histos_Signal[key])
        hs.Add(histos_BG[key])
        hs.Add(histos_MC[key])
        #hs.Add(histos_LSB[key])
        #hs.Add(histos_RSB[key])
        #maxHisto = max(hs.GetMaximum(), histos_SR[key].GetMaximum())
        hs.SetMaximum(maxY[key])#maxHisto)
        hs.Draw('HIST')
        histos_SR[key].Draw('same')
        leg=pad1.BuildLegend(.65,.93,.95,.7)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("SAME")
        pad2 = TPad("pad2", "The pad 20% of the height",0.0,0.01,1.0,0.2)
        pad2.cd()
        h_ratio[key].Draw()
        xmin = h_ratio[key].GetBinLowEdge(1)
        xmax = h_ratio[key].GetBinLowEdge(h_ratio[key].GetNbinsX())+h_ratio[key].GetBinWidth(h_ratio[key].GetNbinsX())
        line0 = TLine(xmin,1,xmax,1)
        line0.SetLineColor(ROOT.kBlue)
        line0.SetLineWidth(2)
        line0.Draw('same')
        canvas.cd()
        pad1.Draw()
        pad2.Draw()
        canvas.Update()
        canvas.Print(outputFile_name)
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'.png')
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'.pdf')
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'.root')

        # Make comparisison between chib1 and chib2
        histos_SR1[key].Scale(1./histos_SR1[key].Integral()) 
        histos_SR2[key].Scale(1./histos_SR2[key].Integral())
        histos_Signal1[key].Scale(1./histos_Signal1[key].Integral()) 
        histos_Signal2[key].Scale(1./histos_Signal2[key].Integral())
        # histos_SR1[key].Rebin(3)
        # histos_SR2[key].Rebin(3)
        # histos_Signal1[key].Rebin(3)
        # histos_Signal2[key].Rebin(3)
        hsc = THStack('hsc', titles_chib[key])
        hsc.Add(histos_Signal1[key])
        hsc.Add(histos_Signal2[key])
        #hsc.SetMaximum(maxY[key])
        hsc.Draw('nostack')
        legc = canvas.BuildLegend(.65,.93,.95,.7)
        legc.SetFillStyle(0)
        legc.SetBorderSize(0)
        legc.Draw("SAME")
        print 'kolmogorov prob SR     '+key+' '+str( histos_SR1[key].KolmogorovTest(histos_SR2[key]))
        print 'kolmogorov prob Signal '+key+' '+str( histos_Signal1[key].KolmogorovTest(histos_Signal2[key]))
        cfr = 2 if key not in ['dimuon_pt', 'chib_pt'] else 5
        kolmogorovProb_chi = round(histos_Signal1[key].KolmogorovTest(histos_Signal2[key]),cfr)
        pvtxt = TPaveText(.2, .0, 0.5, 0.1,"NDC")
        pvtxt.SetFillStyle(0)
        pvtxt.SetBorderSize(0)
        pvtxt.SetTextSize(0.04)
        pvtxt.AddText('Kolmogorov prob = '+str(kolmogorovProb_chi))
        pvtxt.Draw()
        canvas.Update()
        canvas.Print(outputFile_name)
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'_comparison'+'.png')
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'_comparison'+'.pdf')
        canvas.Print('plots/'+key+'_'+mc_ptSpectrum+'_comparison'+'.root')

   
    canvas.Clear()
    canvas.Print(outputFile_name+"]") #chiude file .ps

    for key in keys:
        print '\n'
        print key+' Kolmogorov probability: '+str(chi2Prob[key])
        print key+' chiSqure probability: '+str(kolmogorovProb[key])
    #print 'Kolmogorov probability with sidebands: '+str(histos_MC['dimuon_pt'].KolmogorovTest(histos_LSB['dimuon_pt'])) # 
    
    inputFile.Close()

if __name__ == '__main__':
    makeHistos()
    makePlots()
