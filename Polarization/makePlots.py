#!/usr/bin/python

from ROOT import gSystem, gROOT,  gStyle, TFile, TCanvas, TH1, TH1D, TTree, THStack
import ROOT, sys 
from polarization.py import *

gROOT.SetStyle('Plain') # white background
gStyle.SetOptStat(222211)
TH1.SetDefaultSumw2(True) #Cosi` in automatico per tutti gli isto dovrebbe abilitare il calcolo degli errori come si deve (come somma dei pesi)

outputFile_name = 'CutsPlots.ps'

def makeHistos():

    rootupla_Chib1_Sel = '../../store/ChiB_1P_1_UpsilonPt.root'
    rootupla_Chib2_Sel = '../../store/ChiB_1P_2_UpsilonPt.root'
    rootupla_2012Data  = '../../store/2012_AllData.root'

    # Histograms
    histos_chib1 = dict()
    for frame in ('hx', 'cs'):
        for helicity in ('1', '2', 'u'):
            histos_chib[frame+'_'+helicity] = dict()

    for frame in ('hx', 'cs'):
        for helicity in ('1', '2', 'u'): 
            histos_chib1['invm1S'] = TH1D('invm1S_chib1'+'_'+frame+'_'+helicity,'chib1;m_{#mu#mu#gamma} - m_{#mu#mu} + m_{#Upsilon}^{PDG};',80,9.7,10.1)
            histos_chib1['Y1S_nsigma'] = TH1D("Y1S_nsigma_chib1"+'_'+frame+'_'+helicity,"chib1; Number of theoric standad deviation from #Upsilon(1S);",100,0.,20.)
            histos_chib1['probFit1S'] = TH1D("probFit1S_chib1"+'_'+frame+'_'+helicity, "chib1;fit probability based on #chi^{2} over number dof;",100,0.,1.)
            histos_chib1['photon_pt'] = TH1D("photon_pt_chib1"+'_'+frame+'_'+helicity,"chib1;p_{T} [GeV];",100,0.,5.)
            histos_chib1['photon_eta'] = TH1D("photon_eta_chib1"+'_'+frame+'_'+helicity, "chib1;#eta;",100,-3.,3.)
            histos_chib1['dimuon_pt'] = TH1D("dimuon_pt_chib1"+'_'+frame+'_'+helicity,"chib1;p_{T} [GeV];",100,0.,40.)
            histos_chib1['dimuon_rapidity'] = TH1D("dimuon_rapidity_chib1"+'_'+frame+'_'+helicity,"chib1;y;",100,-1.8,1.8)
            histos_chib1['dz'] = TH1D("dz_chib1"+'_'+frame+'_'+helicity,"chib1;dz [cm];",100, 0., 0.5)


    tree = TTree()
    chibTreeName = 'rootuple/chibTree'

    def fillHistos(inputFileName, histos):
        inputFile = TFile.Open(inputFileName,"READ")
        inputFile.GetObject(chibTreeName, tree)
        for event in tree:
            for key, histo in histos.items():
                if( key == 'Y1S_nsigma' or key =='probFit1S' or key =='dz' or 
                   (event.Y1S_nsigma < 2.5 and event.probFit1S > 0.02 and event.dz < 0.1 and 
                    event.photon_pt > 0.5 and event.dimuon_pt > 7 and abs(event.photon_eta) < 1. and 
                    abs(event.dimuon_rapidity) < 1.25 and event.invm1S > 9.87 and event.invm1S < 9.92)):# and event.pi0_abs_mass > 0.025)):
                    upsdir, leptondir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4,frame='hx'):
                    weight = angDist(upsdir, leptondir, chib_state=1, helicity=0)
                    histo.Fill(getattr(event, key), weight)
        inputFile.Close()

    print 'Filling MC histograms'
    fillHistos(rootupla_Chib1_Sel, histos_MC)
    #fillHistos(rootupla_Chib2_Sel, histos_MC)
    

    outputFile = TFile('CutsPlots.root','RECREATE')
    for histo in  histos_MC.values():
        histo.Write()
    outputFile.Close()

def makePlots():

    inputFile = TFile('CutsPlots.root','READ')

    histos_data = dict()
    histos_MC = dict()

    keys = ['Y1S_nsigma', 'probFit1S', 'photon_pt', 'photon_eta', 'dimuon_pt', 'dimuon_rapidity', 'dz', 'invm1S'] 

    for key in keys:
        histos_data[key] =  inputFile.Get(key+'_data')
        histos_MC[key] =  inputFile.Get(key+'_MC')

    titles = dict(
        Y1S_nsigma = 'Y1S_nsigma;Number of theoric standad deviation from #Upsilon(1S);',
        probFit1S = 'probFit1S;fit probability based on #chi^{2} over number dof;',
        photon_pt = 'photon p_{T};p_{T} [GeV];',
        photon_eta = 'photon_eta;#eta;',
        dimuon_pt = '#Upsilon p_{T};p_{T} [GeV];',
        dimuon_rapidity = '#Upsilon rapidity;y;',
        dz = '|dz|;|dz| [cm];',
        invm1S = 'Q-Value;m_{#mu#mu#gamma} - m_{#mu#mu} + m_{#Upsilon}^{PDG};'
        )
    
    logY = dict(
        Y1S_nsigma = 0,
        probFit1S = 1,
        photon_pt = 0,
        photon_eta = 0,
        dimuon_pt = 0,
        dimuon_rapidity = 0,
        dz = 0,
        invm1S = 0
        )

    for histo in histos_data.values() + histos_MC.values():
        if(histo.Integral() != 0):
            histo.Scale(1./histo.Integral()) 
        

    for histo in histos_data.values():
        histo.SetLineColor(ROOT.kGreen+2)
        histo.SetLineWidth(2)
        # histo.SetMarkerColor(ROOT.kGreen+2)
        # histo.SetMarkerStyle(21)
                
    for histo in histos_MC.values():
        histo.SetLineColor(ROOT.kBlue)
        histo.SetLineWidth(2)
        # histo.SetMarkerColor(ROOT.kBlue)
        # histo.SetMarkerStyle(20)
        
    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+"[") #apre file .ps

    for key in keys:
        canvas.SetLogy(logY[key])
        hs = THStack('hs', titles[key])
        hs.Add(histos_data[key])
        hs.Add(histos_MC[key])
        hs.Draw('nostack')
        leg=canvas.BuildLegend(.8,.95,.92,.8);
        leg.SetFillColor(0);
        leg.Draw("SAME");
        canvas.Update()
        canvas.Print(outputFile_name)
        canvas.Print('CutsPlots/'+key+'.png')

    canvas.Print(outputFile_name+"]") #chiude file .ps
    inputFile.Close()

if __name__ == '__main__':
    makeHistos()
    makePlots()
