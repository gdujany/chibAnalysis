#!/usr/bin/python
from ROOT import gSystem
gSystem.Load('My_double_CB/My_double_CB_cxx')
from ROOT import My_double_CB

from ROOT import RooDataSet, RooRealVar, RooArgSet, RooFormulaVar, RooGenericPdf, RooCmdArg, RooStats
from ROOT import RooCBShape, RooAddPdf, RooArgList, RooPlot, RooDataHist, RooFitResult, RooAbsPdf, RooGaussian, RooChebychev
from ROOT import RooFit, gROOT, TStyle, gStyle, gPad
from ROOT import TFile, TCanvas, TH1F, TGraphErrors, TPad, TLegend, TPaveText, TMultiGraph, TGraphErrors, TMath, TLine
from ROOT import TH1D, TH1F, TTree, RooHistPdf
import ROOT
import sys, getopt
#import ChicUtils
from array import array
#from numpy import savetxt
from pyUtils import *


def doDataFit(Chib1_parameters,Chib2_parameters, cuts, inputfile_name = None, RooDataSet = None, ptBin_label='', plotTitle = "#chi_{b}",fittedVariable='qValue', printSigReso = False, noPlots = False, useOtherSignalParametrization = False, useOtherBackgroundParametrization = False, massFreeToChange = False, sigmaFreeToChange = False, legendOnPlot=True, drawPulls=False, titleOnPlot=False, cmsOnPlot=True, printLegend=True):

    if RooDataSet != None:
        dataSet = RooDataSet 
    elif inputfile_name != None:
        print "Creating DataSet from file "+str(inputfile_name)
        dataSet = makeRooDataset(inputfile_name)
    else:
        raise ValueError('No dataset and no inputfile passed to function doDataFit')
    
    if(fittedVariable == 'refittedMass'):
        x_var = 'rf1S_chib_mass'
        output_suffix = '_refit'
        x_axis_label= 'm_{#mu^{+} #mu^{-} #gamma} [GeV]'
    else:
        x_var = 'invm1S'
        output_suffix = '_qValue'
        x_axis_label = 'm_{#gamma #mu^{+} #mu^{-}} - m_{#mu^{+} #mu^{-}} + m^{PDG}_{#Upsilon}  [GeV]'
    
    cuts_str = str(cuts)
    #cuts_str = quality_cut + "photon_pt > 0.5 && abs(photon_eta) < 1.0 && ctpv < 0.01  && abs(dimuon_rapidity) < 1.3 && pi0_abs_mass > 0.025 &&  abs(dz) < 0.5"
    data = dataSet.reduce( RooFit.Cut(cuts_str) )
    
    print 'Creating pdf'
    x=RooRealVar(x_var, 'm(#mu #mu #gamma) - m(#mu #mu) + m_{#Upsilon}',9.7,10.1,'GeV')
    numBins = 80 # define here so that if I change it also the ndof change accordingly
    x.setBins(numBins)
    
    # cristal balls
    mean_1 = RooRealVar("mean_1","mean ChiB1",Chib1_parameters.mean,"GeV")
    sigma_1 = RooRealVar("sigma_1","sigma ChiB1",Chib1_parameters.sigma,'GeV')
    a1_1 = RooRealVar('#alpha1_1', '#alpha1_1', Chib1_parameters.a1)
    n1_1 = RooRealVar('n1_1', 'n1_1', Chib1_parameters.n1)
    a2_1 = RooRealVar('#alpha2_1', '#alpha2_1',Chib1_parameters.a2)
    n2_1 = RooRealVar('n2_1', 'n2_1', Chib1_parameters.n2)
    parameters = RooArgSet(a1_1, a2_1, n1_1, n2_1)
    
    mean_2 = RooRealVar("mean_2","mean ChiB2",Chib2_parameters.mean,"GeV")
    sigma_2 = RooRealVar("sigma_2","sigma ChiB2",Chib2_parameters.sigma,'GeV')
    a1_2 = RooRealVar('#alpha1_2', '#alpha1_2', Chib2_parameters.a1)
    n1_2 = RooRealVar('n1_2', 'n1_2', Chib2_parameters.n1)
    a2_2 = RooRealVar('#alpha2_2', '#alpha2_2', Chib2_parameters.a2)
    n2_2 = RooRealVar('n2_2', 'n2_2', Chib2_parameters.n2)
    parameters.add(RooArgSet( a1_2, a2_2, n1_2, n2_2))

    if massFreeToChange:
        # scale_mean = RooRealVar('scale_mean', 'Scale that multiplies masses found with MC', 0.8,1.2)
        # mean_1_fixed = RooRealVar("mean_1_fixed","mean ChiB1",Chib1_parameters.mean,"GeV")
        # mean_2_fixed = RooRealVar("mean_2_fixed","mean ChiB2",Chib2_parameters.mean,"GeV")
        # mean_1 = RooFormulaVar("mean_1",'@0*@1', RooArgList(scale_mean, mean_1_fixed))
        # mean_2 = RooFormulaVar("mean_2",'@0*@1', RooArgList(scale_mean, mean_2_fixed))
        variazione_m = 0.05 # 50 MeV
        diff_m_12 = RooRealVar('diff_m_12', 'Difference between masses chib1 and chib2',0.021,'GeV') # 21 MeV from PDG
        mean_1=RooRealVar("mean_1","mean ChiB1",Chib1_parameters.mean,Chib1_parameters.mean-variazione_m,Chib1_parameters.mean+variazione_m ,"GeV")
        mean_2=RooFormulaVar('mean_2', '@0+@1',RooArgList(mean_1, diff_m_12))
        # mean_2=RooRealVar("mean_2","mean ChiB2",Chib2_parameters.mean,Chib2_parameters.mean-variazione_m,Chib2_parameters.mean+variazione_m ,"GeV")
        parameters.add(mean_1)
    else:
        parameters.add(RooArgSet(mean_1, mean_2))
        
    
    chib1_pdf = My_double_CB('chib1', 'chib1', x, mean_1, sigma_1, a1_1, n1_1, a2_1, n2_1)
    chib2_pdf = My_double_CB('chib2', 'chib2', x, mean_2, sigma_2, a1_2, n1_2, a2_2, n2_2)
    
    if sigmaFreeToChange:
        scale_sigma = RooRealVar('scale_sigma', 'Scale that multiplies sigmases found with MC', 1, 1.1)#1.01
        sigma_1_fixed = RooRealVar("sigma_1","sigma ChiB1",Chib1_parameters.sigma,'GeV')
        sigma_2_fixed = RooRealVar("sigma_2","sigma ChiB2",Chib2_parameters.sigma,'GeV')
        sigma_1 = RooFormulaVar("sigma_1",'@0*@1', RooArgList(scale_sigma, sigma_1_fixed))
        sigma_2 = RooFormulaVar("sigma_2",'@0*@1', RooArgList(scale_sigma, sigma_2_fixed))
        parameters.add(scale_sigma)
    else:
        parameters.add(RooArgSet(sigma_1, sigma_2))

    chib1_pdf = My_double_CB('chib1', 'chib1', x, mean_1, sigma_1, a1_1, n1_1, a2_1, n2_1)
    chib2_pdf = My_double_CB('chib2', 'chib2', x, mean_2, sigma_2, a1_2, n1_2, a2_2, n2_2)

    if useOtherSignalParametrization: # In this case I redefine cb_pdf
        cb1 = RooCBShape('cb1', 'cb1', x, mean_1, sigma_1, a1_1, n1_1)
        cb2 = RooCBShape('cb2', 'cb2', x, mean_2, sigma_2, a1_2, n1_2)
        # I use a2 as the sigma of my gaussian 
        gauss1 = RooCBShape('gauss1', 'gauss1',x, mean_1, a2_1, a1_1, n1_1)
        gauss2 = RooCBShape('gauss2', 'gauss2',x, mean_2, a2_2, a1_2, n1_2)
        # I use n2 as the ratio of cb with respect to gauss 
        chib1_pdf = RooAddPdf('chib1','chib1',RooArgList(cb1, gauss1),RooArgList(n2_1))
        chib2_pdf = RooAddPdf('chib2','chib2',RooArgList(cb2, gauss2),RooArgList(n2_2))
        
    
    #background
    q01S_Start = 9.5
    alpha   =   RooRealVar("#alpha","#alpha",1.5,0.2,3.5)
    beta    =   RooRealVar("#beta","#beta",-2.5,-7.,0.)
    q0      =   RooRealVar("q0","q0",q01S_Start)#,9.5,9.7)
    delta   =   RooFormulaVar("delta","TMath::Abs(@0-@1)",RooArgList(x,q0))
    b1      =   RooFormulaVar("b1","@0*(@1-@2)",RooArgList(beta,x,q0))
    signum1 =   RooFormulaVar( "signum1","( TMath::Sign( -1.,@0-@1 )+1 )/2.", RooArgList(x,q0) )
    
    
    background = RooGenericPdf("background","Background", "signum1*pow(delta,#alpha)*exp(b1)", RooArgList(signum1,delta,alpha,b1) )

    if useOtherBackgroundParametrization: # in thies case I redefine background
        a0 = RooRealVar('a0','a0',1.,-1.,1.) #,0.5,0.,1.)
        a1 = RooRealVar('a1','a1',0.1,-1.,1.) #-0.2,0.,1.)
        #a2 = RooRealVar('a2','a2',-0.1,1.,-1.)
        background = RooChebychev('background','Background',x,RooArgList(a0,a1))
        parameters.add(RooArgSet(a0, a1))
    else:
        parameters.add(RooArgSet(alpha, beta, q0))

    #together
    chibs = RooArgList(chib1_pdf,chib2_pdf,background)    
    
    # ndof
    floatPars = parameters.selectByAttrib("Constant",ROOT.kFALSE)
    ndof = numBins - floatPars.getSize() - 1

    # # Here I have as parameters N1, N2, and N_background
    # n_chib1 = RooRealVar("n_chib1","n_chib1",1250, 0, 50000)
    # n_chib2 =  RooRealVar("n_chib2","n_chib2",825, 0, 50000)
    # n_background = RooRealVar('n_background','n_background',4550, 0, 50000)
    # ratio_list = RooArgList(n_chib1, n_chib2, n_background)
    # modelPdf = RooAddPdf('ModelPdf', 'ModelPdf', chibs, ratio_list)

    # Here I have as parameters N_12, ratio_12, N_background
    n_chib = RooRealVar("n_chib","n_chib",2075, 0, 100000)
    ratio_21 = RooRealVar("ratio_21","ratio_21",0.6, 0, 1)
    n_chib1 = RooFormulaVar("n_chib1","@0/(1+@1)",RooArgList(n_chib, ratio_21))
    n_chib2 = RooFormulaVar("n_chib2","@0/(1+1/@1)",RooArgList(n_chib, ratio_21))
    n_background = RooRealVar('n_background','n_background',4550, 0, 50000)
    ratio_list = RooArgList(n_chib1, n_chib2, n_background)
    parameters.add(RooArgSet(n_chib1, n_chib2, n_background))
    modelPdf = RooAddPdf('ModelPdf', 'ModelPdf', chibs, ratio_list)
    
    print 'Fitting to data'
    fit_region = x.setRange("fit_region",9.7,10.1)
    result=modelPdf.fitTo(data,RooFit.Save(), RooFit.Range("fit_region"))
    
        
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
    chib1P_set = RooArgSet(chib1_pdf)
    modelPdf.plotOn(frame,RooFit.Components(chib1P_set), RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(2), RooFit.LineWidth(1))
    chib2P_set = RooArgSet(chib2_pdf)
    modelPdf.plotOn(frame, RooFit.Components(chib2P_set),RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(2), RooFit.LineWidth(1))
    background_set =  RooArgSet(background)
    modelPdf.plotOn(frame,RooFit.Components(background_set), RooFit.LineColor(ROOT.kBlack), RooFit.LineStyle(2), RooFit.LineWidth(1))
    modelPdf.plotOn(frame, RooFit.LineWidth(2))
    frame.SetName("fit_resonance")  

    # Make numChib object
    numChib = NumChib(numChib=n_chib.getVal(), s_numChib=n_chib.getError(), ratio_21=ratio_21.getVal(), s_ratio_21=ratio_21.getError(), numBkg=n_background.getVal(), s_numBkg=n_background.getError(), corr_NB=result.correlation(n_chib, n_background),corr_NR=result.correlation(n_chib, ratio_21) , name='numChib'+output_suffix+ptBin_label,q0=q0.getVal(),s_q0=q0.getError(),alpha=alpha.getVal(),s_alpha=alpha.getError(), beta=beta.getVal(), s_beta=beta.getError(), chiSquare=frame.chiSquare())
    #numChib.saveToFile('numChib'+output_suffix+'.txt')

    if noPlots:
        chi2 = frame.chiSquare()
        del frame
        return numChib, chi2
    
    # Legend
    parameters_on_legend = RooArgSet()#n_chib, ratio_21, n_background)
    if massFreeToChange:
        #parameters_on_legend.add(scale_mean)
        parameters_on_legend.add(mean_1)
        #parameters_on_legend.add(mean_2)
    if sigmaFreeToChange:
        parameters_on_legend.add(scale_sigma)
    if massFreeToChange or sigmaFreeToChange:
        modelPdf.paramOn(frame, RooFit.Layout(0.1,0.6,0.2),RooFit.Parameters(parameters_on_legend))
    
    if printLegend: #chiquadro, prob, numchib...
        leg = TLegend(0.48,0.75,0.97,0.95)
        leg.SetBorderSize(0)
        #leg.SetTextSize(0.04)
        leg.SetFillStyle(0)
        chi2 = frame.chiSquare()
        probChi2 = TMath.Prob(chi2*ndof, ndof)
        chi2 = round(chi2,2)
        probChi2 = round(probChi2,2)
        leg.AddEntry(0,'#chi^{2} = '+str(chi2),'')
        leg.AddEntry(0,'Prob #chi^{2} = '+str(probChi2),'')
        N_bkg, s_N_bkg = roundPair(numChib.numBkg, numChib.s_numBkg)
        leg.AddEntry(0,'N_{bkg} = '+str(N_bkg)+' #pm '+str(s_N_bkg),'')
        N_chib12, s_N_chib12 = roundPair(numChib.numChib, numChib.s_numChib)
        leg.AddEntry(0,'N_{#chi_{b}} = '+str(N_chib12)+' #pm '+str(s_N_chib12),'')
        Ratio = numChib.calcRatio()
        s_Ratio = numChib.calcRatioError()
        Ratio, s_Ratio = roundPair(Ratio, s_Ratio)
        leg.AddEntry(0,'N_{2}/N_{1} = '+str(Ratio)+' #pm '+str(s_Ratio),'')

        if printSigReso: # Add Significance
            Sig =  numChib.calcSignificance()
            s_Sig = numChib.calcSignificanceError()
            Sig, s_Sig = roundPair(Sig, s_Sig)
            leg.AddEntry(0,'Sig = '+str(Sig)+' #pm '+str(s_Sig),'')
            if(Chib1_parameters.sigma>Chib2_parameters.sigma):
                Reso = Chib1_parameters.sigma * 1000 # So it's in MeV and not in GeV
                s_Reso = Chib1_parameters.s_sigma * 1000 # So it's in MeV and not in GeV
            else:
                Reso = Chib2_parameters.sigma * 1000 # So it's in MeV and not in GeV
                s_Reso = Chib2_parameters.s_sigma * 1000 # So it's in MeV and not in GeV
            Reso, s_Reso =roundPair(Reso, s_Reso)
            leg.AddEntry(0,'Reso = '+str(Reso)+' #pm '+str(s_Reso)+' MeV','')
            #N1 = numChib.numChib1
            #s_N1 = numChib.s_numChib1
            #N1, s_N1 = roundPair(N1, s_N1)
            #leg.AddEntry(0,'N_{1} = '+str(N1)+' #pm '+str(s_N1),'')
            #N2 = numChib.numChib2
            #s_N2 = numChib.s_numChib2
            #N2, s_N2 = roundPair(N2, s_N2)
            #leg.AddEntry(0,'N_{2} = '+str(N2)+' #pm '+str(s_N2),'')

        frame.addObject(leg)

    if legendOnPlot:  #  < pT <
        legend = TLegend(.06,.75,.53,.93)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        #legend.AddEntry(0,'CMS','')
        legend.AddEntry(0,str(cuts.upsilon_pt_lcut)+' GeV < p_{T}(#Upsilon) < '+str(cuts.upsilon_pt_hcut)+' GeV','')
        #legend.AddEntry(0,'p_{T}(#Upsilon)<'+str(cuts.upsilon_pt_hcut),'')
        frame.addObject(legend)

    if titleOnPlot:
        titleLegend = TLegend(.06,.4,.55,.73)
       
        #titleLegend.SetTextSize(0.03)
        titleLegend.SetFillStyle(0)
        titleLegend.SetBorderSize(0)
        titleLegend.AddEntry(0,plotTitle,'')
        frame.addObject(titleLegend)

    if cmsOnPlot:
        if printLegend:
            pvtxt = TPaveText(.1,.55,.55,.73,"NDC")
        else:
            pvtxt = TPaveText(0.5,0.75,0.97,0.9,"NDC") #(.06,.4,.55,.73)
        pvtxt.AddText('CMS Preliminary')
        pvtxt.AddText('pp, #sqrt{s} = 8 TeV')
        pvtxt.AddText('L = 20.7 fb^{-1}')
        pvtxt.SetFillStyle(0)
        pvtxt.SetBorderSize(0)
        pvtxt.SetTextSize(0.04)
        frame.addObject(pvtxt)
    
    # Canvas
    c1=TCanvas('Chib12_1P'+output_suffix+ptBin_label,'Chib12_1P'+output_suffix+ptBin_label)
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
        line0 = TLine(9.7, 0, 10.1, 0)
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
    #c1.SaveAs('Chib12_1P'+output_suffix+'.png')
    print 'Chi2 = '+str(frame.chiSquare())
    

    # print ratio background/all in the signal refion
    signal_region = x.setRange("signal_region",9.87,9.92)
    pdf_integral = modelPdf.createIntegral(RooArgSet(x), RooFit.Range('signal_region')).getVal() * (n_chib.getVal() + n_background.getVal())
    bkg_integral = background.createIntegral(RooArgSet(x), RooFit.Range('signal_region')).getVal() * n_background.getVal()

    print 'Ratio bkg/all in signal region = '+str(bkg_integral/pdf_integral)

    return numChib, c1
    
# MAIN
if __name__ == '__main__':

    #gROOT.SetBatch()        # don't pop up canvases
    #gROOT.SetStyle('Plain') # white background
    #gStyle.SetOptStat(222211)

    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style
    
    inputfile_name = "../store/2012_AllData.root"

    ptBins = [7,11,16,20,40]
    

    # GRAB OPTIONS
    # Defaults
    fittedVariable='qValue'
    ptBin = None
    makeAll = False
    pas_label=''
    printLegend = True
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fab:',['refit', 'all', 'ptBin=','pas'])
    except getopt.GetoptError:
        print './dataFit.py [-f] [-b <1, 2, 3, 4>] [--pas]'
        print '-f = --refit   -a = -all  -b = --ptBin'
        print '--pas makes figure for pas without right legend'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            fittedVariable='refittedMass'
        if opt in ('-a', '--all'):
            makeAll = True
        if opt in ('-b', '--ptBin'):
            ptBin = int(arg)
        if opt ==  '--pas':
            pas_label='_forPAS'
            printLegend = False

          
  
    cuts = Cuts()
    cuts = cuts.loadFromFile("cuts.txt")
  
    def make(ptBin=None, fittedVariable='qValue', inputfile_name=None, RooDataSet=None):

        if ptBin == None:
            ptBin_label = ''
            cuts.upsilon_pt_lcut = ptBins[0]
            cuts.upsilon_pt_hcut = ptBins[-1]
            ptBin_addTitle = ''
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            cuts.upsilon_pt_lcut = ptBins[ptBin-1]
            cuts.upsilon_pt_hcut = ptBins[ptBin]
            ptBin_addTitle = str(ptBins[ptBin-1])+' GeV <p_{T#Upsilon} < '+str(ptBins[ptBin])+' Gev'

        Chib1_parameters = CB_parameters()
        Chib2_parameters = CB_parameters()
        if(fittedVariable=='refittedMass'):
            print 'Fitting kinematic refitted mass'
            fittedVariable='refittedMass'
            Chib1_parameters = Chib1_parameters.loadFromFile("txt/CB_ChiB1{0}_refit.txt".format(ptBin_label))
            Chib2_parameters = Chib2_parameters.loadFromFile("txt/CB_ChiB2{0}_refit.txt".format(ptBin_label))
        else:
            print 'Fitting Q-value + m_PDG'
            fittedVariable='qValue'
            Chib1_parameters = Chib1_parameters.loadFromFile("txt/CB_ChiB1{0}_qValue.txt".format(ptBin_label))
            Chib2_parameters = Chib2_parameters.loadFromFile("txt/CB_ChiB2{0}_qValue.txt".format(ptBin_label))


        if(fittedVariable == 'refittedMass'): cuts.varLabel = 'rf1S_'
        print str(cuts)

        numChib, c1 = doDataFit(inputfile_name=inputfile_name,RooDataSet=RooDataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, plotTitle = "#chi_{b} "+ptBin_addTitle,ptBin_label=ptBin_label, useOtherSignalParametrization = False, useOtherBackgroundParametrization = False, massFreeToChange = False, sigmaFreeToChange = False, printLegend=printLegend)
        numChib.saveToFile('txt/'+numChib.name+'.txt')
        c1.SaveAs('plots/'+c1.GetName()+pas_label+'.png')
        c1.SaveAs('plots/'+c1.GetName()+pas_label+'.pdf')
        c1.SaveAs('plots/'+c1.GetName()+pas_label+'.root')


    if makeAll:
        print "Creating DataSet from file "+str(inputfile_name)
        dataSet = makeRooDataset(inputfile_name)
        for ptBin in [None] + range(1,len(ptBins)):
            for fittedVariable in ('refittedMass', ):#'qValue'
                make(ptBin=ptBin, fittedVariable=fittedVariable, RooDataSet=dataSet) 
    else:
        make(ptBin=ptBin, fittedVariable=fittedVariable, inputfile_name=inputfile_name)
