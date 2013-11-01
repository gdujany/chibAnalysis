#!/usr/bin/python
from ROOT import gSystem, gROOT, gStyle, TFile, TTree, TH1, TH1D, TCanvas, THStack, TF1
import ROOT, sys, getopt, math
sys.path.append('Polarization')
from pyUtils import *
from polarization import *
from array import array

ptBins = [7,11,16,20,40]

def makeEffHistos(sel_file_name, gen_file_name, cuts=Cuts(), chib_state=1, frame='hx', outputFile_name='effHistos/chib.root', ptSpectrum_rw='2S', ptSpectrum_file='2S'):
    """
    Make file with histograms containing variables used to compute efficiency and ratio efficiency(heliciy)/efficiency(unpolarized) with error
    """
    # I increase the cut to avoid effects on the acceptance due to the difference between generated and reconstructed pT
    cuts.upsilon_pt_lcut = 5
    cuts.upsilon_pt_hcut = 50

    # inFile_ptSpectra = TFile('ptSpectra.root','READ')
    # h_ptSpectrum = dict()
    # for key in ['1S', '2S', '3S', 'flat']:
    #     h_ptSpectrum[key] = inFile_ptSpectra.Get(key)

    chib_state = int(chib_state)
    
    gen_file = TFile.Open(gen_file_name,'READ')
    sel_file = TFile.Open(sel_file_name,'READ')
    
    gen_tree = TTree()
    gen_file.GetObject('rootuple/GenParticlesTree', gen_tree)
    sel_tree = TTree()
    sel_file.GetObject('rootuple/chibTree', sel_tree)

    variables = ['gen_w', 'gen_w2', 'gen_wu', 'sel_w', 'sel_w2', 'sel_wu']
    keys = ['N_noWeight']
    for variable in variables:
        for helicity in range(chib_state+1)+['u']:
            keys.append(variable+'_'+str(helicity))

    histos = dict()
    for key in keys:
        histos[key] = TH1D(key, key, 70, 5., 40.)
    
    print 'Computing number generated particles'
    for event in gen_tree:
        if event.Upsilon_pt > cuts.upsilon_pt_lcut and event.Upsilon_pt < cuts.upsilon_pt_hcut and event.photon_pt > cuts.photon_pt_cut and abs(event.photon_eta) < cuts.photon_eta_cut and abs(event.Upsilon_rapidity) < cuts.upsilon_rapidity_cut and abs(event.muP_p4.Eta()) < cuts.muon_eta_cut and abs(event.muM_p4.Eta()) < cuts.muon_eta_cut and event.muP_p4.Pt() > cuts.muon_pt_cut and event.muM_p4.Pt() > cuts.muon_pt_cut:
            ups_pt = event.Upsilon_pt
            chib_pt = event.chib_p4.Pt()
            histos['N_noWeight'].Fill(ups_pt)
            ups_dir, mu_dir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4,frame)
            weigh_ptSpectrum = givePtweight(pt=chib_pt, outPt=ptSpectrum_rw, inPt=ptSpectrum_file)#h_ptSpectrum[ptSpectrum_rw].GetBinContent(h_ptSpectrum[ptSpectrum_rw].FindBin(event.chib_p4.Pt())) / h_ptSpectrum[ptSpectrum_file].GetBinContent(h_ptSpectrum[ptSpectrum_file].FindBin(event.chib_p4.Pt()))
            weight = dict()
            for helicity in range(chib_state+1):
                weight[helicity] = angDist(ups_dir, mu_dir, chib_state, helicity) * weigh_ptSpectrum#*(4*math.pi)**2
            weight['u'] = angDist(ups_dir, mu_dir, chib_state, None) * weigh_ptSpectrum
            for hel in weight.keys():
                histos['gen_w_'+str(hel)].Fill(ups_pt, weight[hel])
                histos['gen_w2_'+str(hel)].Fill(ups_pt, weight[hel]**2)
                histos['gen_wu_'+str(hel)].Fill(ups_pt, weight[hel]*weight['u'])
    

    print 'Computing number selected particles'
    for event in sel_tree:
        if event.dimuon_pt > cuts.upsilon_pt_lcut and event.dimuon_pt < cuts.upsilon_pt_hcut and event.photon_pt > cuts.photon_pt_cut and abs(event.photon_eta) < cuts.photon_eta_cut and abs(event.dimuon_rapidity) < cuts.upsilon_rapidity_cut and abs(event.dz)<cuts.dz_cut and abs(event.muonP_eta) < cuts.muon_eta_cut and abs(event.muonM_eta) < cuts.muon_eta_cut and event.muonP_pt > cuts.muon_pt_cut and event.muonM_pt > cuts.muon_pt_cut and event.Y1S_nsigma<2.5 and event.probFit1S>0.02:
            ups_pt = event.Upsilon_p4.Pt()
            chib_pt = event.chib_p4.Pt()
            ups_dir, mu_dir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4, frame)
            weigh_ptSpectrum = givePtweight(pt=chib_pt, outPt=ptSpectrum_rw, inPt=ptSpectrum_file) #h_ptSpectrum[ptSpectrum_rw].GetBinContent(h_ptSpectrum[ptSpectrum_rw].FindBin(event.chib_p4.Pt())) / h_ptSpectrum[ptSpectrum_file].GetBinContent(h_ptSpectrum[ptSpectrum_file].FindBin(event.chib_p4.Pt()))
            weight = dict()
            for helicity in range(chib_state+1):
                weight[helicity] = angDist(ups_dir, mu_dir, chib_state, helicity) * weigh_ptSpectrum #*(4*math.pi)**2
            weight['u'] = angDist(ups_dir, mu_dir, chib_state, None) * weigh_ptSpectrum
            for hel in weight.keys():
                histos['sel_w_'+str(hel)].Fill(ups_pt, weight[hel])
                histos['sel_w2_'+str(hel)].Fill(ups_pt, weight[hel]**2)
                histos['sel_wu_'+str(hel)].Fill(ups_pt, weight[hel]*weight['u'])
    
    outputFile = TFile(outputFile_name, 'recreate')
    for histo in histos.values():
        histo.Write()
    outputFile.Close()
    #inFile_ptSpectra.Close()


def makePlotsEff(inFile1_name = 'effHistos/chib1_hx_refit.root', inFile2_name = 'effHistos/chib2_hx_refit.root', outputFile_name = 'effHistos/plots/plotsEff'):
    '''
    Make plots of efficiency vs pt upsilon for both chib1 and chib2 and save the histo in a root file and in a ps
    '''
    inFile1 = TFile(inFile1_name,'read')
    inFile2 = TFile(inFile2_name,'read')
       
    gROOT.SetStyle('Plain') # white background
    gStyle.SetOptStat(222211)
    gStyle.SetOptFit(1111)
    gStyle.SetStatY(0.89) # Set y-position (fraction of pad size)
    gStyle.SetStatX(0.5) # Set x-position (fraction of pad size)
    gStyle.SetStatW(0.2) # Set width of stat-box (fraction of pad size)
    gStyle.SetStatH(0.2) # Set height of stat-box (fraction of pad size)

    
    def makeHistoEff(numChib=1):
        gStyle.SetOptStat(0000)
        h_eff = dict()
        keys = ['gen_w_u', 'gen_w2_u',  'sel_w_u', 'sel_w2_u', 'N_noWeight']
 
        inFile = inFile2
        if numChib==1: inFile = inFile1
        for key in keys:
            h_eff[key] = inFile.Get(key)

        lowHedges = array('d',[i/10. for i in range(60, 110, 5)]+range(11,16)+range(16,20,2)+range(20,40,4))
        lowHedges = array('d',range(6,10)+range(10,14,2)+range(14,20,3)+range(20,40,5))
        lowHedges = array('d',range(6,18)+range(17,23,2)+range(23,32,3)+range(32,36,3)) # ok for full stat
        lowHedges = array('d',range(7,16)+range(16,20,2)+range(20,32,3)+range(32,41,4))
        #lowHedges = array('d',range(7, 41,2))
        
        for key in keys:
            h_eff[key] = h_eff[key].Rebin(len(lowHedges)-1,key+'_rebinned', lowHedges)
        
        h_eff['eff'] = h_eff['gen_w_u'].Clone()
        h_eff['eff'].SetNameTitle('chib'+str(numChib),'Efficiency chib'+str(numChib)+';p_{T#Upsilon} [GeV];')
        for i in range(h_eff['eff'].GetNbinsX()+2):
            try:
                efficiency = Efficiency(**dictEffArgs(gen_w = h_eff['gen_w_u'].GetBinContent(i), 
                                                      gen_w2 = h_eff['gen_w2_u'].GetBinContent(i),
                                                      sel_w = h_eff['sel_w_u'].GetBinContent(i),
                                                      sel_w2 = h_eff['sel_w2_u'].GetBinContent(i),
                                                      N_noWeight =h_eff['N_noWeight'].GetBinContent(i)))
                h_eff['eff'].SetBinContent(i, efficiency.eff)
                h_eff['eff'].SetBinError(i, efficiency.s_eff)
            except ZeroDivisionError:                
                h_eff['eff'].SetBinContent(i, 0) ### set here value if denominator is 0
                h_eff['eff'].SetBinError(i, 0)

        return h_eff['eff']
    
    h_chib1 = makeHistoEff(1)
    h_chib2 = makeHistoEff(2)

    h_ratio = h_chib1.Clone()
    h_ratio.SetNameTitle('ratio','#epsilon_{1}/#epsilon_{2};p_{T}(#Upsilon) [GeV];')
    for i in range(h_ratio.GetNbinsX()+2):
            try:
                ratio = h_chib1.GetBinContent(i)/h_chib2.GetBinContent(i)
                s_ratio = ratio*math.sqrt((h_chib1.GetBinError(i)/h_chib1.GetBinContent(i))**2 + (h_chib2.GetBinError(i)/h_chib2.GetBinContent(i))**2)
                h_ratio.SetBinContent(i, ratio)
                h_ratio.SetBinError(i, s_ratio)
            except ZeroDivisionError:
                h_ratio.SetBinContent(i, 0) ### set here value if denominator is 0
                h_ratio.SetBinError(i, 0)

    fa = TF1("fa","[0]",6,40)
    fa.SetParameter(0, 1)
    fa.SetParName(0,"Constant");
    
    outputFile = TFile(outputFile_name+'.root', 'recreate')
    h_chib1.Write()
    h_chib2.Write()
    outputFile.Close()

    h_chib1.SetLineColor(ROOT.kBlue)
    h_chib1.SetLineWidth(2)
    h_chib2.SetLineColor(ROOT.kRed)
    h_chib2.SetLineWidth(2)
    h_ratio.SetLineColor(ROOT.kGreen+2)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMinimum(0.6)
    h_ratio.SetMaximum(1.2)

    f_eff1 = TF1('f_eff1', '1++x++x^2',6, 40)
    f_eff1.SetLineColor(h_chib1.GetLineColor())
    f_eff1.SetLineStyle(2)
    f_eff1.SetLineWidth(1)
    f_eff2 = f_eff1.Clone()
    f_eff2.SetName('f_eff2')
    f_eff2.SetLineColor(h_chib2.GetLineColor())

    h_chib1.Fit('f_eff1', 'I')
    h_chib2.Fit('f_eff2', 'I')

    f_ratio = TF1('f_ratio','f_eff1/f_eff2', 6, 40)
    f_ratio.SetLineColor(ROOT.kGreen+2)
    f_ratio.SetLineStyle(2)
    f_ratio.SetLineWidth(1)

    #h_ratio.GetListOfFunctions().Add(f_ratio)
    #h_ratio.Fit('f_ratio')
    


    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+".ps[") #apre file .ps
    for histo in (h_chib1, h_chib2):
        canvas.SetLogy(0)
        # histo.Fit("1++x++x^2")
        # histo.GetFunction("1++x++x^2").SetLineColor(histo.GetLineColor())
        # histo.GetFunction("1++x++x^2").SetLineStyle(2)
        # histo.GetFunction("1++x++x^2").SetLineWidth(1)
        histo.Draw()
        canvas.Update()
        histo.Draw()
        canvas.Update()
        canvas.Print(outputFile_name+'.ps')
        canvas.Print(outputFile_name+'_'+histo.GetName()+'.pdf')
        canvas.Print(outputFile_name+'_'+histo.GetName()+'.root')
    
    # canvas.SetLogy(0)
    # h_ratio.Draw()
    # f_ratio.Draw('same')
    # canvas.Update()
    # canvas.Print(outputFile_name+'.ps')  
    # canvas.Print(outputFile_name+'_'+'ratio.pdf')
    # canvas.Print(outputFile_name+'_'+'ratio.root')

    h_ratio.Fit("fa")
    h_ratio.SetMaximum(1.4)
    h_ratio.Draw()
    f_ratio.Draw('same')
    canvas.Update()
    canvas.Print(outputFile_name+'.ps')
    canvas.Print(outputFile_name+'_'+'ratio.pdf')
    canvas.Print(outputFile_name+'_'+'ratio.root')
    
    # h_ratio.Fit("1++x")
    # h_ratio.Draw()
    # canvas.Update()
    # canvas.Print(outputFile_name+'.ps') 
    # #canvas.Print(outputFile_name+'_'+'ratio_lin.pdf')
    # #canvas.Print(outputFile_name+'_'+'ratio_lin.root')


    # hs = THStack('hs', 'Efficiencies;p_{T#Upsilon} [GeV];')
    # hs.Add(h_chib1)
    # hs.Add(h_chib2)
    # hs.Draw('nostack')
    # leg=canvas.BuildLegend(.70,1.,.94,.92);
    # leg.SetFillColor(0);
    # leg.Draw("SAME");
    # canvas.Update()    
    # canvas.Print(outputFile_name+'.ps')
    canvas.Print(outputFile_name+".ps]") #chiude file .ps





def makeDictPolEff(inFile_name,  chib_state):
    '''
    make dictionnary containing the ratioPolUnpol efficiencies (pairs value and error), the key is a string:
    equal to the helicity ('u', '0', ...) for the "all data" case, and helicity_ptMin_ptMax for the ptBins
    '''

    inFile = TFile(inFile_name,'read')
    h_chib = dict()
       
    variables = ['gen_w', 'gen_w2', 'gen_wu', 'sel_w', 'sel_w2', 'sel_wu']
    keys = ['N_noWeight']
    for variable in variables:
        for helicity in ['u']+range(chib_state+1):
            keys.append(variable+'_'+str(helicity))
    for key in keys:
        h_chib[key] = inFile.Get(key)

    def getRatioPolUnpolEff(chib_state, helicity, ptMin, ptMax):
        binMin = h_chib['gen_w_u'].FindBin(ptMin)
        binMax = h_chib['gen_w_u'].FindBin(ptMax)
        
        ratioPolUnpolEff = RatioEff_Polar_Unpolar(
            gen_w = h_chib['gen_w_'+str(helicity)].Integral(binMin, binMax),
            gen_w2 = h_chib['gen_w2_'+str(helicity)].Integral(binMin, binMax),
            gen_wu = h_chib['gen_wu_'+str(helicity)].Integral(binMin, binMax),
            sel_w = h_chib['sel_w_'+str(helicity)].Integral(binMin, binMax),
            sel_w2 = h_chib['sel_w2_'+str(helicity)].Integral(binMin, binMax),
            sel_wu = h_chib['sel_wu_'+str(helicity)].Integral(binMin, binMax),
            u_gen_w = h_chib['gen_w_u'].Integral(binMin, binMax),
            u_gen_w2 = h_chib['gen_w2_u'].Integral(binMin, binMax),
            u_sel_w = h_chib['sel_w_u'].Integral(binMin, binMax),
            u_sel_w2 = h_chib['sel_w2_u'].Integral(binMin, binMax),
            N_noWeight = h_chib['N_noWeight'].Integral(binMin, binMax),
            label = 'chib'+str(chib_state)+' hel = '+str(helicity)+' bin = '+str(ptMin)+' - '+str(ptMax),
            )

        return (ratioPolUnpolEff.CalcRatioEff(), ratioPolUnpolEff.CalcRatioEffError())
    
    # dictionnary containing the ratioPolUnpol efficiencies, the key is a string:
    # equal to the helicity ('u', '0', ...) for the all data case, and helicity_ptMin_ptMax for the ptBins
    ratioPolUnpolEff_dict = dict()
    for hel in ['u']+range(chib_state+1):
       ratioPolUnpolEff_dict[str(hel)] = getRatioPolUnpolEff(chib_state, hel, ptMin=ptBins[0], ptMax=ptBins[-1])
       for ptBin in range(1,len(ptBins)):
           ratioPolUnpolEff_dict[str(hel)+'_'+str(ptBin-1)+'_'+str(ptBin)] = getRatioPolUnpolEff(chib_state, hel, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        
    return ratioPolUnpolEff_dict


def findMinMaxPol(ratioPolUnpolEff1, ratioPolUnpolEff2, ptBin):
    addKey = '' 
    if ptBin != None: addKey = '_'+str(ptBin-1)+'_'+str(ptBin)
    polSists = []
    for hel1 in ['u', 0, 1]:
            for hel2 in ['u', 0, 1, 2]:
                if hel1 == 'u' and hel2 =='u':
                    continue
                r1 = ratioPolUnpolEff1[str(hel1)+addKey][0]
                r2 = ratioPolUnpolEff2[str(hel2)+addKey][0]
                polSist = r1/r2 - 1
                polSists.append(polSist)
    polSists.sort()
    return polSists[0], polSists[-1]


def makePolarizationTable(ratioPolUnpolEff1, ratioPolUnpolEff2,  outFile_name, printErrors=True, printAllData=True):
    '''
    make the nice table with ratio efficiencies
    '''

    # Now I print the table
    pol_label = {'u' : 'Unpol.', 0 : '0', 1 : '$\\pm 1$', 2 : '$\\pm 2$'}
    outFile = open(outFile_name,'w')
    numColumns = len(ptBins)
    if printAllData:
        numColumns += 1
    outFile.write('\\begin{center} \n \\begin{tabular}{p{0.15\linewidth}*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')
    outFile.write(' &  \\multicolumn{'+str(numColumns-1)+'}{c}{$p_T(\\Upsilon)$ [GeV]} \\\\ \n')
    outFile.write('\\centering Polarization scenario \\\\($m_{\\chi_{b1}}$, $m_{\\chi_{b2}}$) ')
    for ptBin in range(1, len(ptBins)):
        outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
    if printAllData:
        outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
    outFile.write('\\\\ \n \\midrule \n')
    for hel1 in ['u', 0, 1]:
        for hel2 in ['u', 0, 1, 2]:
            if hel1 == 'u' and hel2 =='u':
                continue
            outFile.write('\\centering ('+pol_label[hel1]+', '+pol_label[hel2]+')')
            for ptBin in range(1, len(ptBins)):
                r1 = ratioPolUnpolEff1[str(hel1)+'_'+str(ptBin-1)+'_'+str(ptBin)][0]
                r2 = ratioPolUnpolEff2[str(hel2)+'_'+str(ptBin-1)+'_'+str(ptBin)][0]
                s_r1 = ratioPolUnpolEff1[str(hel1)+'_'+str(ptBin-1)+'_'+str(ptBin)][1]
                s_r2 = ratioPolUnpolEff2[str(hel2)+'_'+str(ptBin-1)+'_'+str(ptBin)][1]
                ratioRatioEff = r1/r2
                s_ratioRatioEff = ratioRatioEff * sqrt((s_r1/r1)**2 + (s_r2/r2)**2)
                #val, err = roundPair(ratioRatioEff, s_ratioRatioEff, 1)
                val = '{:.3f}'.format(ratioRatioEff)
                err = '{:.3f}'.format(s_ratioRatioEff)
                outFile.write(' & '+val+' ')
                if printErrors:
                    outFile.write('$\\pm$ '+err+' ')
            if printAllData:
                r1 = ratioPolUnpolEff1[str(hel1)][0]
                r2 = ratioPolUnpolEff2[str(hel2)][0]
                s_r1 = ratioPolUnpolEff1[str(hel1)][1]
                s_r2 = ratioPolUnpolEff2[str(hel2)][1]
                ratioRatioEff = r1/r2
                s_ratioRatioEff = ratioRatioEff * sqrt((s_r1/r1)**2 + (s_r2/r2)**2)
                val, err = roundPair(ratioRatioEff, s_ratioRatioEff, 1)
                outFile.write(' & '+val+' ')
                if printErrors:
                    outFile.write('$\\pm$ '+err+' ')
            outFile.write('\\\\ \n')
    outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
    outFile.close()


def makeEfficiency(inFile_name, ptMin, ptMax):
    '''
    return efficiency from histograms file
    '''
    inFile = TFile(inFile_name,'read')
    h_eff = dict()
    keys = ['gen_w_u', 'gen_w2_u',  'sel_w_u', 'sel_w2_u', 'N_noWeight']
    for key in keys:
        h_eff[key] = inFile.Get(key)

    ## TMP rebinno gli istogrammi come i plots delle efficienze cosi` sono piu` coerente per un check
    # lowHedges = array('d',range(6,16)+range(16,20,2)+range(20,32,3)+range(32,41,4))
    #for key in keys:
        #     h_eff[key] = h_eff[key].Rebin(len(lowHedges)-1,key+'_rebinned', lowHedges)
        #h_eff[key] = h_eff[key].Rebin(2)

    binMin = h_eff['gen_w_u'].FindBin(ptMin)
    binMax = h_eff['gen_w_u'].FindBin(ptMax)-1
       
    return Efficiency(**dictEffArgs(gen_w = h_eff['gen_w_u'].Integral(binMin, binMax),
                                    gen_w2 = h_eff['gen_w2_u'].Integral(binMin, binMax),
                                    sel_w = h_eff['sel_w_u'].Integral(binMin, binMax),
                                    sel_w2 = h_eff['sel_w2_u'].Integral(binMin, binMax),
                                    N_noWeight = h_eff['N_noWeight'].Integral(binMin, binMax)
                                    )
                        )

def makeEfficiencyXMethod(h_eff, h_ptY, ptMin, ptMax):

    binEdges = []
    for i in range(1, h_eff.GetXaxis().GetNbins()+2):
        binEdges.append(h_eff.GetBinLowEdge(i))
    h_ptY = h_ptY.Rebin(len(binEdges)-1,'ptY',array('d',binEdges))

    binMin = h_eff.FindBin(ptMin)
    binMax = h_eff.FindBin(ptMax)-1
    sumEff = 0
    for binNum in range(binMin, binMax+1):
        if h_eff.GetBinContent(binNum) != 0:
            sumEff +=  h_ptY.GetBinContent(binNum)/h_eff.GetBinContent(binNum)
    return Efficiency(eff = h_ptY.Integral(binMin, binMax)/sumEff , s_eff = 0.000001)
    

def makeEfficiencyTable(inFile1_name, inFile2_name, outFile_name, printAllData=False):
    '''
    make table of ratio efficiencies in the various bins
    '''

    outFile = open(outFile_name,'w')
    numColumns = len(ptBins)
    if printAllData:
        numColumns += 1
    outFile.write('\\begin{center} \n \\begin{tabular}{c*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')
    outFile.write('$p_T(\\Upsilon)$ [GeV]')
    for ptBin in range(1, len(ptBins)):
        outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
    if printAllData:
        outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
    outFile.write('\\\\ \n \\midrule \n')
    outFile.write('$\\frac{\\epsilon_1}{\\epsilon_2}$')
    for ptBin in range(1, len(ptBins)):
        eff1 = makeEfficiency(inFile_name=inFile1_name, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        eff2 = makeEfficiency(inFile_name=inFile2_name, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        ratio = eff1.eff/eff2.eff
        s_ratio = ratio * sqrt((eff1.s_eff/eff1.eff)**2 + (eff2.s_eff/eff2.eff)**2)
        val, err = roundPair(ratio, s_ratio)
        outFile.write(' & '+val+' $\\pm$ '+err+' ')
    if printAllData:
        eff1 = makeEfficiency(inFile_name=inFile1_name, ptMin=ptBins[0], ptMax=ptBins[-1])
        eff2 = makeEfficiency(inFile_name=inFile2_name, ptMin=ptBins[0], ptMax=ptBins[-1])
        ratio = eff1.eff/eff2.eff
        s_ratio = ratio * sqrt((eff1.s_eff/eff1.eff)**2 + (eff2.s_eff/eff2.eff)**2)
        val, err = roundPair(ratio, s_ratio)
        outFile.write(' & '+val+' $\\pm$ '+err+' ')
    outFile.write('\\\\ \n')
    outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
    outFile.close()


def makeEfficiencyTableXMethod(varLabel, ptY_label, outFile_name, ptFileLabel='Upsilon2SPt', ptSpectrum_rw='2S', printAllData=False, useData=False, useMC=False):
    '''
    make table of ratio efficiencies in the various bins
    '''
    inFileEff_name = 'effHistos/plots/plotsEff_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
    inFileEff = TFile(inFileEff_name,'read')
    h_eff1 = inFileEff.Get('chib1')
    h_eff2 = inFileEff.Get('chib2')

    if useData or useMC:
        inFilePtY_name = 'AgreementDataMC/agreementDataMC_'+ptY_label+'_plots.root'
        ptYFile = TFile(inFilePtY_name, 'read')
        if useData:
            h_ptY = ptYFile.Get('dimuon_pt_signal') 
        else:
            h_ptY = ptYFile.Get('dimuon_pt_MC') 
            
        h_ptY1 = h_ptY
        h_ptY2 = h_ptY
    else:
        inFilePtY1_name = 'effHistos/chib1_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
        inFilePtY2_name = 'effHistos/chib2_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        ptYFile1 = TFile(inFilePtY1_name, 'read')
        ptYFile2 = TFile(inFilePtY2_name, 'read')  
        h_ptY1 = ptYFile1.Get('sel_w_u')
        h_ptY2 = ptYFile2.Get('sel_w_u')
        

    outFile = open(outFile_name,'w')
    numColumns = len(ptBins)
    if printAllData:
        numColumns += 1
    outFile.write('\\begin{center} \n \\begin{tabular}{c*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')
    outFile.write('$p_T(\\Upsilon)$ [GeV]')
    for ptBin in range(1, len(ptBins)):
        outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
    if printAllData:
        outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
    outFile.write('\\\\ \n \\midrule \n')
    outFile.write('$\\frac{\\epsilon_1}{\\epsilon_2}$')
    for ptBin in range(1, len(ptBins)):
        eff1 = makeEfficiencyXMethod(h_eff1, h_ptY1, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        eff2 = makeEfficiencyXMethod(h_eff2, h_ptY2, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        ratio = eff1.eff/eff2.eff
        s_ratio = ratio * sqrt((eff1.s_eff/eff1.eff)**2 + (eff2.s_eff/eff2.eff)**2)
        val, err = roundPair(ratio, s_ratio)
        outFile.write(' & '+val+' $\\pm$ '+err+' ')
    if printAllData:
        eff1 = makeEfficiencyXMethod(h_eff1, h_ptY1, ptMin=ptBins[0], ptMax=ptBins[-1])
        eff2 = makeEfficiencyXMethod(h_eff2, h_ptY1, ptMin=ptBins[0], ptMax=ptBins[-1])
        ratio = eff1.eff/eff2.eff
        s_ratio = ratio * sqrt((eff1.s_eff/eff1.eff)**2 + (eff2.s_eff/eff2.eff)**2)
        val, err = roundPair(ratio, s_ratio)
        outFile.write(' & '+val+' $\\pm$ '+err+' ')
    outFile.write('\\\\ \n')
    outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
    outFile.close()
    

# MAIN
if __name__ == '__main__':
 

    # GRAB OPTIONS
    # Defaults
    fittedVariable='qValue'
    num_chib = 1
    ptBin = None
    makeAll = False
    frame = 'hx'
    ptSpectrum_rw = None
    ptSpectrum_file = '2S'
    makeHistos = False
    makePlots = False
    makeTable = False
    isMakeSistematici = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fac:sb:',['refit', 'all', 'chib=', 'ptBin=', 'cs', 'hx', 'filePt=', 'reweighPt=', 'makeHistos', 'makePlots', 'makeTable', 'makeTables', 'sistematici'])
    except getopt.GetoptError:
        print './efficiency.py [-f] [-c <1, 2>] [-b <1,2,3,4>] [--cs --hx] [--reweighPt <1S, 2S, 3S, flat>] [--filePt <1S, 2S, 3S, flat, full>] [--makeHistos] [--makeTable] [--sistematici]'
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
        if opt == '--cs':
            frame = 'cs'
        if opt == '--hx':
            frame = 'hx'
        if opt == '--reweighPt':
            ptSpectrum_rw = arg
        if opt == '--filePt':
            ptSpectrum_file = arg
        if opt == '--makeHistos':
            makeHistos = True
        if opt == '--makePlots':
            makePlots = True
        if opt in ('--makeTable', '--makeTables'):
            makeTable = True
        if opt in ('-s', '--sistematici'):
            isMakeSistematici = True
        if opt == '--histosForPlot':
            isHistosForPlot = True
        

    if ptSpectrum_file in ('1S', '1s'):
        ptFileLabel = 'UpsilonPt'
        ptSpectrum_file = '1S'
    elif ptSpectrum_file in ('2S', '2s'):
        ptFileLabel = 'Upsilon2SPt'
        ptSpectrum_file = '2S'
    elif ptSpectrum_file in ('3S', '3s'):
        ptFileLabel = 'Upsilon3SPt'
        ptSpectrum_file = '3S'
    elif ptSpectrum_file in ('flat', 'Flat'):
        ptFileLabel = 'FlatPt'
        ptSpectrum_file = 'flat'
    elif ptSpectrum_file in ('full', 'Full'):
        ptFileLabel = 'FullStat'

    ptY_label = ptSpectrum_file

    if ptSpectrum_rw == None:
        ptSpectrum_rw = ptSpectrum_file
        
    if ptSpectrum_rw in ('1S', '1s'):
        ptSpectrum_rw = '1S'
    elif ptSpectrum_rw in ('2S', '2s'):
        ptSpectrum_rw = '2S'
    elif ptSpectrum_rw in ('3S', '3s'):
        ptSpectrum_rw = '3S'
    elif ptSpectrum_rw in ('flat', 'Flat'):
        ptSpectrum_rw = 'flat'


    rootupla_Chib1_Gen = '../store/GenParticles_ChiB_1P_1_'+ptFileLabel+'.root'
    rootupla_Chib1_Sel = '../store/ChiB_1P_1_'+ptFileLabel+'.root'
    rootupla_Chib2_Gen = '../store/GenParticles_ChiB_1P_2_'+ptFileLabel+'.root'
    rootupla_Chib2_Sel = '../store/ChiB_1P_2_'+ptFileLabel+'.root'

    if num_chib in ('2', 2):
        sel_file_name = rootupla_Chib2_Sel
        gen_file_name = rootupla_Chib2_Gen
    else: # num_chib == 1
        sel_file_name = rootupla_Chib1_Sel
        gen_file_name = rootupla_Chib1_Gen
        
  
    varLabel = ('qValue', 'refit')[fittedVariable == 'refittedMass']


    def makeEfficiencyFile(num_chib, varLabel, ptBin, ptFileLabel='Upsilon2SPt', ptSpectrum_rw = '2S'):
        
        if ptBin == None:
            ptBin_label = ''
            ptMin = ptBins[0]
            ptMax = ptBins[-1]
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            ptMin = ptBins[ptBin-1]
            ptMax = ptBins[ptBin]
        

        inFile_name = 'effHistos/chib'+str(num_chib)+'_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' # I take hx as frame because in these plots I look only the unpolarized case so it's the same
        outFile_name = 'txt/eff_chib'+str(num_chib)+'_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'.txt'
        efficiency = makeEfficiency(inFile_name, ptMin, ptMax)
        efficiency.saveToFile(outFile_name)

    def makeEfficiencyFileXMethod(num_chib, varLabel, ptBin, ptFileLabel='Upsilon2SPt', ptY_label='2S', ptSpectrum_rw = '2S', useData=False, useMC=False):
        
        if ptBin == None:
            ptBin_label = ''
            ptMin = ptBins[0]
            ptMax = ptBins[-1]
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            ptMin = ptBins[ptBin-1]
            ptMax = ptBins[ptBin]
       

        inFileEff_name = 'effHistos/plots/plotsEff_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        inFileEff = TFile(inFileEff_name,'read')
        h_eff1 = inFileEff.Get('chib1')
        h_eff2 = inFileEff.Get('chib2')
    
        if useData or useMC:
            inFilePtY_name = 'AgreementDataMC/agreementDataMC_'+ptY_label+'_plots.root'
            ptYFile = TFile(inFilePtY_name, 'read')
            if useData:
                h_ptY = ptYFile.Get('dimuon_pt_signal') 
                outFile1_name = 'txt/eff_chib1_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_xSignal.txt'
                outFile2_name = 'txt/eff_chib2_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_xSignal.txt'
            else:
                h_ptY = ptYFile.Get('dimuon_pt_MC') 
                outFile1_name = 'txt/eff_chib1_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_xMC.txt'
                outFile2_name = 'txt/eff_chib2_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_xMC.txt'
            h_ptY1 = h_ptY
            h_ptY2 = h_ptY
        else:
            inFilePtY1_name = 'effHistos/chib1_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
            inFilePtY2_name = 'effHistos/chib2_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
            ptYFile1 = TFile(inFilePtY1_name, 'read')
            ptYFile2 = TFile(inFilePtY2_name, 'read')  
            h_ptY1 = ptYFile1.Get('sel_w_u')
            h_ptY2 = ptYFile2.Get('sel_w_u')
            outFile1_name = 'txt/eff_chib1_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_x.txt'
            outFile2_name = 'txt/eff_chib2_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_rw+ptBin_label+'_x.txt'        
  
              
        efficiency = makeEfficiencyXMethod(h_eff1, h_ptY1, ptMin, ptMax)
        efficiency.saveToFile(outFile1_name)
        efficiency = makeEfficiencyXMethod(h_eff2, h_ptY2, ptMin, ptMax)
        efficiency.saveToFile(outFile2_name)
      


 
#################################
    if makeHistos:
        print 'Making histos'
        cuts = Cuts()
        cuts = cuts.loadFromFile("cuts.txt")
        outputFile_name = 'effHistos/chib'+str(num_chib)+'_'+frame+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        makeEffHistos(sel_file_name = sel_file_name, gen_file_name = gen_file_name, cuts=cuts, chib_state=int(num_chib), outputFile_name=outputFile_name, ptSpectrum_rw=ptSpectrum_rw, ptSpectrum_file=ptSpectrum_file, frame=frame)


    elif makePlots:
        print 'Making plots'
        inFile1_name = 'effHistos/chib1_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' # I take hx as frame because in these plots I look only the unpolarized case so it's the same
        inFile2_name = 'effHistos/chib2_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        outputFile_name = 'effHistos/plots/plotsEff_'+frame+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel
        
        makePlotsEff(inFile1_name, inFile2_name, outputFile_name)
    
    elif makeTable:
        print 'Making table'
        varLabel = 'refit' # Faccio le tabelle solo per il refit (che e` la baseline) tanto non cambia molto
        for frame in ('hx', 'cs'):
            inFile1_name = 'effHistos/chib1_'+frame+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
            inFile2_name = 'effHistos/chib2_'+frame+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
            outFile_name = 'tables/ratio_efficiencies_'+frame+'.tex'

            ################################################################################
            ### !!!!!!!! SCOMMENTARE PER AVERE ANCHE POLARIZZAZIONI !!!!!!!!!!!!!!!!!!!!!!!!
            ################################################################################

            ratioPolUnpolEff1 = makeDictPolEff(inFile_name=inFile1_name,  chib_state=1)
            ratioPolUnpolEff2 = makeDictPolEff(inFile_name=inFile2_name,  chib_state=2)
            makePolarizationTable(ratioPolUnpolEff1, ratioPolUnpolEff2, outFile_name)
            ################################################################################

        ptFileLabel = 'Upsilon2SPt'
        for ptSpectrum_rw in ['1S', '2S', '3S']: 
            inFile1_name = 'effHistos/chib1_hx'+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
            inFile2_name = 'effHistos/chib2_hx'+'_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
            outFile_name = 'tables/efficiencies_'+ptFileLabel+'_'+ptSpectrum_rw+'.tex'

            makeEfficiencyTable(inFile1_name, inFile2_name, outFile_name)
            makeEfficiencyTableXMethod(varLabel='refit', ptY_label=ptSpectrum_rw, ptSpectrum_rw=ptSpectrum_rw, outFile_name='tables/efficiencies_'+ptFileLabel+'_'+ptSpectrum_rw+'_x.tex')
            makeEfficiencyTableXMethod(varLabel='refit', ptY_label=ptSpectrum_rw, ptSpectrum_rw=ptSpectrum_rw, outFile_name='tables/efficiencies_'+ptFileLabel+'_'+ptSpectrum_rw+'_xMC.tex', useMC=True)
            makeEfficiencyTableXMethod(varLabel='refit', ptY_label=ptSpectrum_rw, ptSpectrum_rw=ptSpectrum_rw, outFile_name='tables/efficiencies_'+ptFileLabel+'_'+ptSpectrum_rw+'_xSignal.tex', useData=True)
        
        # chib1, chib2 = 1S, 2S
        inFile1_name = 'effHistos/chib1_hx'+'_'+'Upsilon2SPt_1S'+'_'+varLabel+'.root' 
        inFile2_name = 'effHistos/chib2_hx'+'_'+'Upsilon2SPt_2S'+'_'+varLabel+'.root'
        outFile_name = 'tables/efficiencies_'+'1S_2S'+'.tex'
        makeEfficiencyTable(inFile1_name, inFile2_name, outFile_name)

        # chib1, chib2 = 2S, 3S
        inFile1_name = 'effHistos/chib1_hx'+'_'+'Upsilon2SPt_2S'+'_'+varLabel+'.root' 
        inFile2_name = 'effHistos/chib2_hx'+'_'+'Upsilon2SPt_3S'+'_'+varLabel+'.root'
        outFile_name = 'tables/efficiencies_'+'2S_3S'+'.tex'
        makeEfficiencyTable(inFile1_name, inFile2_name, outFile_name)

        
  
      
    elif isMakeSistematici:     
        # Max distance due to polarization
        inFile1_name_hx = 'effHistos/chib1_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
        inFile2_name_hx = 'effHistos/chib2_hx_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        inFile1_name_cs = 'effHistos/chib1_cs_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root' 
        inFile2_name_cs = 'effHistos/chib2_cs_'+ptFileLabel+'_'+ptSpectrum_rw+'_'+varLabel+'.root'
        
        ratioPolUnpolEff1_hx = makeDictPolEff(inFile_name=inFile1_name_hx, chib_state=1)
        ratioPolUnpolEff2_hx = makeDictPolEff(inFile_name=inFile2_name_hx, chib_state=2)
        ratioPolUnpolEff1_cs = makeDictPolEff(inFile_name=inFile1_name_cs, chib_state=1)
        ratioPolUnpolEff2_cs = makeDictPolEff(inFile_name=inFile2_name_cs, chib_state=2)
        
        inFiles = dict()
        inFiles['1_default'] = inFile1_name_hx
        inFiles['2_default'] = inFile2_name_hx
        inFiles['1_1S'] = 'effHistos/chib1_hx_Upsilon2SPt_1S'+'_'+varLabel+'.root'
        inFiles['2_1S'] = 'effHistos/chib2_hx_Upsilon2SPt_1S'+'_'+varLabel+'.root'
        inFiles['1_2S'] = 'effHistos/chib1_hx_Upsilon2SPt_2S'+'_'+varLabel+'.root'
        inFiles['2_2S'] = 'effHistos/chib2_hx_Upsilon2SPt_2S'+'_'+varLabel+'.root'
        inFiles['1_3S'] = 'effHistos/chib1_hx_Upsilon2SPt_3S'+'_'+varLabel+'.root'
        inFiles['2_3S'] = 'effHistos/chib2_hx_Upsilon2SPt_3S'+'_'+varLabel+'.root'

        def sistematicoPtSpectrum(ptBin, txtFile):
            if ptBin == None:
                ptMin = ptBins[0]
                ptMax = ptBins[-1]
            else:
                ptMin = ptBins[ptBin-1]
                ptMax = ptBins[ptBin]
            
            eff_sist = dict()
            for key in ['1_default', '2_default', '1_1S', '2_1S', '1_2S', '2_2S', '1_3S', '2_3S']:
                eff_sist[key] = makeEfficiency(inFiles[key], ptMin, ptMax)

            eRatios = dict()
            eRatios['default'] = eff_sist['1_default'].eff/eff_sist['2_default'].eff
            eRatios['1S_2S'] = eff_sist['1_1S'].eff/eff_sist['2_2S'].eff
            eRatios['2S_3S'] = eff_sist['1_2S'].eff/eff_sist['2_3S'].eff
            eRatios['1S_1S'] = eff_sist['1_1S'].eff/eff_sist['2_1S'].eff
            eRatios['3S_3S'] = eff_sist['1_3S'].eff/eff_sist['2_3S'].eff
            
            for key in ['default', '1S_2S', '2S_3S', '1S_1S', '3S_3S']:
                txtFile.write(key+' = '+str(eRatios[key])+'\n')
            

            return max([abs(eRatios['default'] - eRatios[key]) for key in ['3S_3S', '1S_1S', '1S_2S', '2S_3S']])/eRatios['default'] #errore relativo


        txtFile = open('varieEfficienze.txt', 'w')
        for varLabel in ('refit',):# 'qValue'):
            for ptBin in [None] + range(1,len(ptBins)):
                polar_hx =  findMinMaxPol(ratioPolUnpolEff1_hx, ratioPolUnpolEff2_hx, ptBin)
                polar_cs = findMinMaxPol(ratioPolUnpolEff1_cs, ratioPolUnpolEff2_cs, ptBin)
                txtFile.write('\nptBin = '+str(ptBin)+'\n')
                systPtSpectrum = sistematicoPtSpectrum(ptBin, txtFile)
                ptBin_label = ''
                if ptBin != None: 
                    ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
                sistFile_name = 'txt/sistematici_'+varLabel+ptBin_label+'.txt'
                sistematici = Sistematici()
                sistematici = sistematici.loadFromFile(sistFile_name)
                sistematici.polar_hx = polar_hx
                sistematici.polar_cs = polar_cs
                sistematici.ptDistro = systPtSpectrum
                sistematici.saveToFile(sistFile_name)
        

    else:
        if makeAll:
            for num_chib in (1, 2):
                for ptBin in [None] + range(1,len(ptBins)):
                    for varLabel in ('refit', ):#'qValue'):
                        for pt_label in ['1S', '2S', '3S']: 
                            makeEfficiencyFile(num_chib, varLabel, ptBin, ptSpectrum_rw=pt_label)
                            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_rw=pt_label, ptY_label=pt_label)
                            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_rw=pt_label, ptY_label=pt_label, useData=True)
                            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_rw=pt_label, ptY_label=pt_label, useMC=True)

        else:
            makeEfficiencyFile(num_chib, varLabel, ptBin, ptSpectrum_rw=ptSpectrum_label)
            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_label=ptSpectrum_label, ptY_label=ptSpectrum_label)
            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_label=ptSpectrum_label, ptY_label=ptSpectrum_label, useData=True)
            makeEfficiencyFileXMethod(num_chib, varLabel='refit', ptBin=ptBin, ptSpectrum_label=ptSpectrum_label, ptY_label=ptSpectrum_label, useMC=True)
            
