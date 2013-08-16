#!/usr/bin/python

from ROOT import gSystem
gSystem.Load('My_double_CB/My_double_CB_cxx')
from ROOT import My_double_CB

from ROOT import RooDataSet, RooRealVar, RooArgSet, RooFormulaVar, RooGenericPdf, RooCmdArg, RooStats
from ROOT import RooCBShape, RooAddPdf, RooArgList, RooPlot, RooDataHist, RooFitResult, RooAbsPdf, RooGaussian
from ROOT import RooFit, gROOT, TStyle, gStyle, gPad
from ROOT import TFile, TCanvas, TPad, TH1F, TGraphErrors, TPad, TLegend, TPaveText, TMultiGraph, TGraphErrors, TMath
from ROOT import TH1D, TH1F, TTree, RooHistPdf, TLine
import ROOT, sys, getopt
#import ChicUtils
from array import array
#from numpy import savetxt
from pyUtils import *


def doMCFit(inputfile_name, mass_chib, cuts, output_name='ChiB',  plotTitle = "#Chi_{b}", fittedVariable='qValue', returnOnlyFitResult = False, useOtherSignalParametrization = False, drawPulls = False, legendOnPlot=True): 

    mass_error = 0.15 
 
    print "Creating DataSet from file "+str(inputfile_name)
    dataSet = makeRooDataset(inputfile_name)

    if(fittedVariable == 'refittedMass'):
        x_var = 'rf1S_chib_mass'
        output_suffix = '_refit'
        x_axis_label= 'm_{#mu^{+} #mu^{-} #gamma} [GeV]'
    else:
        x_var = 'invm1S'
        output_suffix = '_qValue'
        x_axis_label = 'm_{#gamma #mu^{+} #mu^{-}} - m_{#mu^{+} #mu^{-}} + m^{PDG}_{#Upsilon}  [GeV]'

    cuts_str = str(cuts)
    #cuts_str = quality_cut +"photon_pt > 0.5 && abs(photon_eta) < 1.0 &&  abs(dimuon_rapidity) < 1.0 && dimuon_pt>5.0 && pi0_abs_mass > 0.025 &&  fabs(dz) < 0.1 "#&& numPrimaryVertices < 16"
    data = dataSet.reduce( RooFit.Cut(cuts_str) )
    
    print 'Creating pdf'
    x=RooRealVar(x_var,'m(#mu #mu #gamma) - m(#mu #mu) + m_{#Upsilon}',9.8,9.96)#9.7,10.1,'GeV')
    numBins = 32 # define here so that if I change it also the ndof change accordingly
    x.setBins(numBins)

    
    # Double sided Crystal Ball
    mean=RooRealVar("#mu","mean ChiB",mass_chib, mass_chib-mass_error,mass_chib+mass_error,"GeV")
    sigma=RooRealVar("#sigma","sigma ChiB",0.006, 0,0.2,'GeV')
    a1 = RooRealVar('#alpha1', '#alpha1', 0.75, 0, 3)
    a2 = RooRealVar('#alpha2', '#alpha2', 1.6, 0, 3)
    n1 = RooRealVar('n1', 'n1', 2.8)#, 1.0,4.0) # 2 per 2S
    n2 = RooRealVar('n2', 'n2', 3)#, 1.,4.0) # 1 per 2S
    parameters = RooArgSet(mean, sigma, a1, a2, n1, n2)
    
    cb_pdf = My_double_CB('chib', 'chib', x, mean, sigma, a1, n1, a2, n2)
    #cb_pdf = RooCBShape('chib', 'chib', x, mean, sigma, a1, n1)

    # ndof
    floatPars = parameters.selectByAttrib("Constant",ROOT.kFALSE)
    ndof = numBins - floatPars.getSize() - 1

    if useOtherSignalParametrization: # In this case I redefine cb_pdf
        n1 = RooRealVar('n1', 'n1', 2)
        cb = RooCBShape('cb1', 'cb1', x, mean, sigma, a1, n1)
        # I use a2 as the sigma of my second CB  
        a2 = RooRealVar('#alpha2', '#alpha2', 0.5, 0, 3)
        gauss = RooCBShape('cb2', 'cb2',x, mean, a2, a1, n1)
        # I use n2 as the ratio of cb1 with respect to cb2 
        n2 = RooRealVar('n2', 'n2',0.,1.)
        cb_pdf = RooAddPdf('chib','chib',RooArgList(cb, gauss),RooArgList(n2))
        #cb_pdf = cb
    
    print 'Fitting to data'
    fit_region = x.setRange("fit_region",9.8,9.96)
    result = cb_pdf.fitTo(data, RooFit.Save(), RooFit.Range("fit_region"))
    
    # a1_val = a1.getVal()
    # a2_val = a2.getVal()
    # a1 = RooRealVar('#alpha1', '#alpha1', a1_val)
    # a2 = RooRealVar('#alpha2', '#alpha2', a2_val)

    # n1 = RooRealVar('n1', 'n1', 1.7,1.0,5.)
    # n2 = RooRealVar('n2', 'n2', 1.7,0.,4.)
    # cb_pdf = My_double_CB('chib', 'chib', x, mean, sigma, a1, n1, a2, n2)
    # result = cb_pdf.fitTo(data, RooFit.Save())

    
    
    if returnOnlyFitResult:
        return result
    
    # define frame
    frame = x.frame()
    frame.SetNameTitle("fit_resonance","Fit Resonanace")
    frame.GetXaxis().SetTitle(x_axis_label )
    frame.GetYaxis().SetTitle( "Events/5 MeV " )
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.SetLineWidth(1)
    frame.SetTitle(plotTitle)
    
    # plot things on frame
    data.plotOn(frame, RooFit.MarkerSize(0.7))
    cb_pdf.plotOn(frame, RooFit.LineWidth(2))
    

    # chiSquare legend
    chi2 = frame.chiSquare()
    probChi2 = TMath.Prob(chi2*ndof, ndof)
    chi2 = round(chi2,2)
    probChi2 = round(probChi2,2)
    leg = TLegend(0.3,0,.10,.10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    cb_pdf.paramOn(frame, RooFit.Layout(0.17,0.56,0.93))
    leg.AddEntry(0,'#chi^{2} ='+str(chi2),'')
    leg.AddEntry(0,'Prob #chi^{2} = '+str(probChi2),'')
    leg.SetTextSize(0.04)
    frame.addObject(leg)

    if legendOnPlot:
        legend = TLegend(.08,.5,.55,.7)
        #legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        #legend.AddEntry(0,'CMS','')
        legend.AddEntry(0,str(cuts.upsilon_pt_lcut)+' GeV < p_{T}(#Upsilon) < '+str(cuts.upsilon_pt_hcut)+' GeV','')
        #legend.AddEntry(0,'p_{T}(#Upsilon)<'+str(cuts.upsilon_pt_hcut),'')
        frame.addObject(legend)
        
        title = TLegend(.8,.75,.9,.93)
        title.SetTextSize(0.15)
        title.SetFillStyle(0)
        title.SetBorderSize(0)
        title.AddEntry(0,plotTitle,'')
        frame.addObject(title)
    
     # Canvas
    c1=TCanvas(output_name+output_suffix,output_name+output_suffix)
    frame.Draw()

    if drawPulls:
        #c1=TCanvas(output_name+output_suffix,output_name+output_suffix,700, 625)
        hpull = frame.pullHist()
        framePulls = x.frame()
        framePulls.SetTitle(';;Pulls')
        framePulls.GetYaxis().SetLabelSize(0.18)
        framePulls.GetYaxis().SetTitle('Pulls')
        framePulls.GetYaxis().SetTitleSize(0.18)
        framePulls.GetYaxis().SetTitleOffset(0.15)
        framePulls.GetYaxis().SetNdivisions(005)
        framePulls.GetXaxis().SetLabelSize(0.16)
        framePulls.GetXaxis().SetTitle('')
        line0 = TLine(9.8, 0, 9.96, 0)
        line0.SetLineColor(ROOT.kBlue)
        line0.SetLineWidth(2)
        framePulls.addObject(line0)
        framePulls.addPlotable(hpull,"P") 
        framePulls.SetMaximum(5)
        framePulls.SetMinimum(-5)
        pad1 = TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0)
        pad1.cd()
        frame.Draw()
        pad2 = TPad("pad2", "The pad 20% of the height",0.0,0.01,1.0,0.2)
        pad2.cd()
        framePulls.Draw()
        c1.cd()
        pad1.Draw()
        pad2.Draw()
    #c1.SaveAs(output_name+output_suffix+'.png')
    print 'Chi2 = '+str(frame.chiSquare())

    pdf_parameters = CB_parameters(mean=mean.getVal(), sigma=sigma.getVal(), n1=n1.getVal(), n2=n2.getVal(), a1=a1.getVal(), a2=a2.getVal(), s_mean=mean.getError(), s_sigma=sigma.getError(), s_n1=n1.getError(), s_n2=n2.getError(), s_a1=a1.getError(), s_a2=a2.getError(), chiSquare = frame.chiSquare())
    #pdf_parameters.saveToFile("CB_"+output_name+output_suffix+".txt")
    return pdf_parameters, c1
    


# MAIN
if __name__ == '__main__':

    #gROOT.SetBatch()        # don't pop up canvases
    #gROOT.SetStyle('Plain') # white background
    #gStyle.SetOptStat(222211)
    
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style


    rootupla_Chib1_Sel = '../store/ChiB_1P_1_Upsilon2SPt.root'
    rootupla_Chib2_Sel = '../store/ChiB_1P_2_Upsilon2SPt.root'
    
    mass_chib1_1P =  9.89278 
    mass_chib2_1P =  9.91221 
    leg = 0

    ptBins = [7,11,16,20,40]
    
    # GRAB OPTIONS
    # Defaults
    fittedVariable='qValue'
    num_chib = 1
    ptBin = None
    makeAll = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fac:b:',['refit', 'all', 'chib=', 'ptBin='])
    except getopt.GetoptError:
        print './dataFit.py [-f] [-c <1, 2>] [-b <1,2,3,4>]'
        print '-f = --refit  -a = --all   -c = --chib   -b = --ptBin' 
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            fittedVariable='refittedMass'
        if opt in ('-a', '--all'):
            makeAll = True
        if opt in ('-c', '--chib'):
            num_chib = arg
        if opt in ('-b', '--ptBin'):
            ptBin = int(arg)
    
    cuts = Cuts()
    cuts = cuts.loadFromFile("cuts.txt")

    def make(num_chib=1, ptBin=None, fittedVariable='qValue'):

        if num_chib in ('2', 2):
            print 'processing chiB_1P_2'
            inputfile_name = rootupla_Chib2_Sel
            mass_chib = mass_chib2_1P
            output_name = 'ChiB2'
            plotTitle = "#chi_{b2}"
        else: # num_chib == 1
            print 'processing ChiB_1P_1'
            inputfile_name = rootupla_Chib1_Sel
            mass_chib = mass_chib1_1P
            output_name = 'ChiB1'
            plotTitle = "#chi_{b1}"   

        if ptBin == None:
            ptBin_label = ''
            cuts.upsilon_pt_lcut = ptBins[0]
            cuts.upsilon_pt_hcut = ptBins[-1]
            ptBin_addTitle = ''
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            cuts.upsilon_pt_lcut = ptBins[ptBin-1]
            cuts.upsilon_pt_hcut = ptBins[ptBin]
            ptBin_addTitle = ''#str(ptBins[ptBin-1])+' GeV <p_{T#Upsilon} < '+str(ptBins[ptBin])+' Gev'

        print str(cuts)
        if fittedVariable == 'refittedMass': cuts.varLabel = 'rf1S_'

        pdf_parameters, c1 = doMCFit(inputfile_name=inputfile_name, mass_chib=mass_chib, cuts=cuts, output_name=output_name+ptBin_label, plotTitle=plotTitle+ptBin_addTitle, fittedVariable=fittedVariable,useOtherSignalParametrization = False, drawPulls = True)
        pdf_parameters.saveToFile("txt/CB_"+c1.GetName()+".txt")   
        c1.SaveAs('plots/'+c1.GetName()+'.png')
        c1.SaveAs('plots/'+c1.GetName()+'.pdf')
        c1.SaveAs('plots/'+c1.GetName()+'.root')
    

    if makeAll:
        for num_chib in (1,2):
            for ptBin in [None] + range(1,len(ptBins)):
                for fittedVariable in ('refittedMass', 'qvalue'):
                    make(num_chib=num_chib, ptBin=ptBin, fittedVariable=fittedVariable) 
    else:
        make(num_chib=num_chib, ptBin=ptBin, fittedVariable=fittedVariable)

