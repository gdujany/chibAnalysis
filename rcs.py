#!/usr/bin/python
from __future__ import division
from pyUtils import *
import sys, getopt, ROOT
from array import array
from ROOT import TGraphAsymmErrors, TCanvas, TMultiGraph, TLegend, TPaveText

# Branchin ratios
BF_chib1 = 33.9 #35.
s_BF_chib1 = 2.2 #8.
BF_chib2 = 19.1 #22.
s_BF_chib2 = 1.2 #4.

ptBins = [7,11,16,20,40]


def scrPrint(rcs):
    
    val, err = roundPair(rcs.ratioEfficiency(),rcs.s_ratioEfficiency())
    s_rel_ratioEfficiency = rcs.s_ratioEfficiency()/rcs.ratioEfficiency()
    print 'eff_chib1/eff_chib2 = '+str(val)
    print 'Errore =              '+str(err)
    print 'Errore relativo =     '+str(round(s_rel_ratioEfficiency*100,1))+'%'+'\n'

    val, err = roundPair(rcs.ratioN(),rcs.s_ratioN())
    s_rel_ratioN = rcs.s_ratioN()/rcs.ratioN()
    print 'N2/N1 =           '+str(val)
    print 'Errore =          '+str(err)
    print 'Errore relativo = '+str(round(s_rel_ratioN*100,1))+'%'+'\n'

    val, err = roundPair(rcs.R(),rcs.s_R())
    s_rel_R = rcs.s_R()/rcs.R()
    print 'R =               '+str(val)
    print 's_R =             '+str(err)
    print 'Errore relativo = '+str(round(s_rel_R*100,1))+'%'+'\n'


    cs_ratios = dict(
        rcs = rcs.rcs(),
        s_stat = rcs.s_stat(),
        s_sist = rcs.rcs() * rcs.s_rel_sist(),
        s_BR = rcs.s_BR(),
        s_polar_hx_l = rcs.rcs() * (rcs.s_rel_polar_hx())[0],
        s_polar_hx_h = rcs.rcs() * rcs.s_rel_polar_hx()[1],
        s_polar_cs_l = rcs.rcs() * rcs.s_rel_polar_cs()[0],
        s_polar_cs_h = rcs.rcs() * rcs.s_rel_polar_cs()[1],
        )
    
    cs_ratios = roundDict(cs_ratios)

    print 'rcs         =  '+cs_ratios['rcs']
    print 'Statistico  =  '+cs_ratios['s_stat']
    print 'Sistematico =  '+cs_ratios['s_sist']
    print 'BR          =  '+cs_ratios['s_BR']
    print 'HX          =  +'+cs_ratios['s_polar_hx_h']+' '+cs_ratios['s_polar_hx_l']
    print 'CS          =  +'+cs_ratios['s_polar_cs_h']+' '+cs_ratios['s_polar_cs_l']


def makeDictRCS(varLabel, ptLabel):
    rcs = dict()
    numChibs = dict()
    chib1_par = dict()
    chib2_par = dict()

    for ptBin in [None]+range(1, len(ptBins)):
        ptBin_label = '' if ptBin == None else '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
        numChib = NumChib()
        numChib = numChib.loadFromFile('txt/numChib_'+varLabel+ptBin_label+'.txt')
        chib1_eff = Efficiency()
        chib1_eff = chib1_eff.loadFromFile('txt/eff_chib1'+'_'+varLabel+'_'+'Upsilon2SPt'+'_'+ptLabel+ptBin_label+'.txt')
        chib2_eff = Efficiency()
        chib2_eff = chib2_eff.loadFromFile('txt/eff_chib2'+'_'+varLabel+'_'+'Upsilon2SPt'+'_'+ptLabel+ptBin_label+'.txt')
        sistematici = Sistematici()
        sistematici = sistematici.loadFromFile('txt/sistematici_'+varLabel+ptBin_label+'.txt')
        rcs[ptBin] = RCS(numChib=numChib, chib1_eff=chib1_eff, chib2_eff=chib2_eff, sistematici=sistematici)
        numChibs[ptBin] = numChib

        chib1_params = CB_parameters()
        chib2_params = CB_parameters()
        chib1_params = chib1_params.loadFromFile('txt/CB_ChiB1'+ptBin_label+'_'+varLabel+'.txt')
        chib2_params = chib2_params.loadFromFile('txt/CB_ChiB2'+ptBin_label+'_'+varLabel+'.txt')    
        chib1_par[ptBin] = chib1_params
        chib2_par[ptBin] = chib2_params
        chib_par_dict = {1 : chib1_par, 2 : chib2_par}

    return rcs, numChibs, chib_par_dict
        
  
#### Table RCS       

def makeTableRCS(rcs, printAllData=False):
    '''
    Make table rcs, take as imput the dictionary of rcs as made by makeDictRCS(varLabel, printAllData=True)
    '''
    with open('tables/rcs_'+varLabel+'.tex', 'w') as outFile:
        outFile.write('\\begin{center} \n \\begin{tabular}{p{0.1\columnwidth}crr}\n\\toprule\n')
        outFile.write('$p_T(\\Upsilon)$ [GeV]  &  $\\sigma(\\chi_{b2})/\\sigma(\\chi_{b1})$  &  HX  &  CS \\\\ \n \\midrule \n')
        for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
            if ptBin != None:
                outFile.write('\\multirow{2}*{'+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+'} & ')
            else:
                outFile.write('\\multirow{2}*{'+str(ptBins[0])+' - '+str(ptBins[-1])+'} & ')

            cs_ratio = dict(
                rcs = rcs[ptBin].rcs(),
                s_stat = rcs[ptBin].s_stat(),
                s_sist = rcs[ptBin].rcs() * rcs[ptBin].s_rel_sist(),
                s_BR = rcs[ptBin].s_BR(),
                s_polar_hx_l = abs(rcs[ptBin].rcs() * (rcs[ptBin].s_rel_polar_hx())[0]),
                s_polar_hx_h = rcs[ptBin].rcs() * rcs[ptBin].s_rel_polar_hx()[1],
                s_polar_cs_l = abs(rcs[ptBin].rcs() * rcs[ptBin].s_rel_polar_cs()[0]),
                s_polar_cs_h = rcs[ptBin].rcs() * rcs[ptBin].s_rel_polar_cs()[1],
                )
            
            cs_ratio = roundDict(cs_ratio)
            outFile.write('\\multirow{2}*{'+cs_ratio['rcs']+' $\\pm$ '+cs_ratio['s_stat']+' (stat.) $\\pm$ '+cs_ratio['s_sist']+' (syst.) $\\pm$ '+cs_ratio['s_BR']+' (BR)} &')
            outFile.write('+'+cs_ratio['s_polar_hx_h']+'  &  +'+cs_ratio['s_polar_cs_h']+' \\\\ \n')
            outFile.write(' &  &  $-$'+cs_ratio['s_polar_hx_l']+'  &  $-$'+cs_ratio['s_polar_cs_l']+' \\\\ \n')
        outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
  
##### Table R   

def makeTableR(rcs, printAllData=False):
    '''
    Make table R, take as imput the dictionary of rcs as made by makeDictRCS(varLabel, printAllData=True)
    '''  

    with open('tables/R_'+varLabel+'.tex', 'w') as outFile:
        outFile.write('\\begin{center} \n \\begin{tabular}{ccrr}\n\\toprule\n')
        outFile.write('$p_T(\\Upsilon)$ [GeV]  &  $\\frac{\\sigma(\\chi_{b2}) \\mathcal{B}(\\chi_{b2})}{\\sigma(\\chi_{b1}) \\mathcal{B}(\\chi_{b1})}$  &  HX  &  CS \\\\ \n \\midrule \n')
        for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
            if ptBin != None:
                outFile.write('\\multirow{2}*{'+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+'} & ')
            else:
                outFile.write('\\multirow{2}*{'+str(ptBins[0])+' - '+str(ptBins[-1])+'} & ')

            R = dict(
                val = rcs[ptBin].R(),
                s_stat = rcs[ptBin].s_R(),
                s_sist = rcs[ptBin].R() * rcs[ptBin].s_rel_sist(),
                s_polar_hx_l = abs(rcs[ptBin].R() * (rcs[ptBin].s_rel_polar_hx())[0]),
                s_polar_hx_h = rcs[ptBin].R() * rcs[ptBin].s_rel_polar_hx()[1],
                s_polar_cs_l = abs(rcs[ptBin].R() * rcs[ptBin].s_rel_polar_cs()[0]),
                s_polar_cs_h = rcs[ptBin].R() * rcs[ptBin].s_rel_polar_cs()[1],
                )
            
            R = roundDict(R)
            outFile.write('\\multirow{2}*{'+R['val']+' $\\pm$ '+R['s_stat']+' (stat.) $\\pm$ '+R['s_sist']+' (syst.)} &')
            outFile.write('+'+R['s_polar_hx_h']+'  &  +'+R['s_polar_cs_h']+' \\\\ \n')
            outFile.write(' &  &  $-$'+R['s_polar_hx_l']+'  &  $-$'+R['s_polar_cs_l']+' \\\\ \n')
        outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')


    #### Table N2/N1 
def makeTableN_21(numChibs, printAllData=False):
    '''
    Make table N2/N1, take as imput the dictionary of numChibs as made by makeDictRCS(varLabel, printAllData=True)
    '''

    with open('tables/N21_'+varLabel+'.tex', 'w') as outFile:
        righe = ['numChib1', 'numChib2', 'ratio_21', 'numBkg', 'alpha', 'beta']
        nomi = dict(numChib1 = '$N_{\\chi_{b1}}$',
                    numChib2 = '$N_{\\chi_{b2}}$',
                    ratio_21 = '$N_{\\chi_{b2}}/N_{\\chi_{b1}}$',
                    numBkg = '$N_{bkg}$',
                    alpha = '$\\lambda$',
                    beta = '$\\nu$ $[GeV^{-1}]$')
        numColumns = len(ptBins)
        if printAllData:
            numColumns += 1
        outFile.write('\\begin{center} \n \\begin{tabular}{c*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')
        outFile.write('$p_T(\\Upsilon)$ [GeV] ')
        for ptBin in range(1, len(ptBins)):
            outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
        if printAllData:
            outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
        outFile.write('\\\\ \n \\midrule \n')
        for key in righe:
            outFile.write(nomi[key])
            for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
                precision = 2 if key == 'ratio_21' else 1
                val, err = roundPair(getattr(numChibs[ptBin],key), getattr(numChibs[ptBin],'s_'+key), precision)
                outFile.write(' & '+val+' $\pm$ '+err)    
            outFile.write('\\\\ \n')
        outFile.write('$\\chi^2/ndof$')
        for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
            outFile.write(' & '+ str(round(numChibs[ptBin].chiSquare,2)))
        outFile.write('\\\\ \n')
        outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
   

  
#### Table chib parms from MC   
def makeTableChibParams(chib_par_dict, printAllData=False):
    '''
    Make table chib parms from MC, take as imput the dictionary chib_par_dict as made by makeDictRCS(varLabel, printAllData=True)
    '''
    numRows = len(ptBins)
    with open('tables/chibParameters_'+varLabel+'.tex', 'w') as outFile:
        outFile.write('\\begin{center} \n \\begin{tabular}{*{7}{c}}\n\\toprule\n')
        outFile.write(' & $p_T(\\Upsilon)$ [GeV] & m [GeV] & $\\sigma$ [GeV] & $\\alpha_1$ & $\\alpha_2$ & $\\chi^2/ndof$  \\\\ \n')
        for numChib in [1, 2]:
            chib_params = chib_par_dict[numChib]
            outFile.write('\\midrule \n')
            outFile.write('\\multirow{'+str(numRows-1)+'}*{$\\chi_{b'+str(numChib)+'}$} \n')
            for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
                if ptBin != None:
                    outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin]))
                else:
                    outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1]))
                for key in ['mean', 'sigma', 'a1', 'a2']:
                    val, err =  roundPair(getattr(chib_params[ptBin],key), getattr(chib_params[ptBin],'s_'+key), 1)
                    outFile.write(' & '+val+' $\pm$ '+err)
                outFile.write(' & '+str(round(chib_params[ptBin].chiSquare,2)))
                outFile.write('\\\\ \n')
        outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
              
    
##### Final Plot
def makePlotRCS(rcs):
    
    #x = array('d', [(ptBins[ptBin-1]+ptBins[ptBin])/2. for ptBin in range(1, len(ptBins))])
    x = array('d', [8.71, 12.94, 17.51, 26.02])
    exl = array('d',[(x[ptBin-1]-ptBins[ptBin-1]) for ptBin in range(1, len(ptBins))])
    exh = array('d',[(ptBins[ptBin]-x[ptBin-1]) for ptBin in range(1, len(ptBins))])
    
    y_rcs = array('d',[rcs[ptBin].rcs() for ptBin in range(1, len(ptBins))])
    ey_rcs = array('d',[rcs[ptBin].s_stat() for ptBin in range(1, len(ptBins))])
    
    y_R = array('d',[rcs[ptBin].R() for ptBin in range(1, len(ptBins))])
    ey_R = array('d',[rcs[ptBin].s_R() for ptBin in range(1, len(ptBins))])

    syst_rcs = array('d',[(rcs[ptBin].rcs()*rcs[ptBin].s_rel_sist()) for ptBin in range(1, len(ptBins))])
    syst_R = array('d',[(rcs[ptBin].R()*rcs[ptBin].s_rel_sist()) for ptBin in range(1, len(ptBins))])

    g_rcs = TGraphAsymmErrors(len(x), x, y_rcs, exl, exh, ey_rcs, ey_rcs)
    g_rcs_syst = TGraphAsymmErrors(len(x), x, y_rcs, exl, exh, syst_rcs, syst_rcs)
    
    g_R = TGraphAsymmErrors(len(x), x, y_R, exl, exh, ey_R, ey_R)
    g_R_syst = TGraphAsymmErrors(len(x), x, y_R, exl, exh, syst_R, syst_R)

    g_rcs.SetMarkerStyle(20)
    g_rcs_syst.SetMarkerStyle(20)
    g_R.SetMarkerStyle(21)
    g_R_syst.SetMarkerStyle(21)

    g_rcs.SetFillColor(0)
    g_rcs_syst.SetFillColor(ROOT.kGreen+2)
    g_R.SetFillColor(0)
    g_R_syst.SetFillColor(ROOT.kGreen+2)


    canvas = TCanvas('canvas','canvas')
    mg = TMultiGraph('mg', 'mg')
    mg.Add(g_rcs_syst, 'p2')
    mg.Add(g_rcs, 'psame')
    mg.Add(g_R_syst, 'p2')
    mg.Add(g_R, 'psame')

    mg.SetTitle(';p_{T}(#Upsilon) [GeV];ratio')
    mg.Draw('ap')

    mg.GetYaxis().SetRangeUser(0,1.5)
    mg.GetYaxis().SetTitleSize(0.03)
    #mg.GetYaxis().SetTitleOffset(1.5)
    mg.GetYaxis().SetTitleFont(42)
    mg.GetYaxis().SetTitleSize(0.04)

    ll = TLegend(.2,.17,.4,.35)
    ll.AddEntry(g_rcs,'#sigma_{#chi_{b2}}/#sigma_{#chi_{b1}}','LP')
    ll.AddEntry(g_R,'#frac{#sigma(#chi_{b2})}{#sigma(#chi_{b1})} #frac{B(#chi_{b2} #rightarrow #Upsilon #gamma)}{B(#chi_{b1} #rightarrow #Upsilon #gamma)}','LP')
    ll.SetFillStyle(0)
    ll.SetBorderSize(0)
    ll.SetTextSize(0.03)
    ll.SetTextFont(42)
    ll.Draw()

    pvtxt = TPaveText(.7,.8,0.92,0.93,"NDC")
#pvtxt.AddText('CMS Preliminary 2011')
    pvtxt.AddText('CMS')
    pvtxt.AddText('pp, #sqrt{s} = 8 TeV')
    pvtxt.AddText('L = 20.7 fb^{-1}')
    pvtxt.SetFillStyle(0)
    pvtxt.SetBorderSize(0)
    pvtxt.SetTextSize(0.03)

    
    pvtxt2 = TPaveText(.65,.18,0.75,0.33,"NDC")
    #pvtxt2.AddText('p_{T}(#Upsilon) > 7 GeV')
    pvtxt2.AddText('|y(#Upsilon)| <1.25')
    pvtxt2.AddText('|#eta(#gamma)| <1')
    pvtxt2.AddText('Unpolarized') 
    pvtxt2.SetTextSize(0.03)
    pvtxt2.SetFillStyle(0)
    pvtxt2.SetBorderSize(0)
    pvtxt2.SetTextFont(42)
    pvtxt2.SetTextAlign(11)

    pvtxt.Draw()
    pvtxt2.Draw()

    canvas.SaveAs('plots/rcs.pdf')
    canvas.SaveAs('plots/rcs.root')


if __name__ == '__main__':

    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style

    varLabel = 'qValue'
    ptBin = None
    isMakeTable = False
    isMakePlot = False
    ptSpectrum = '2S'
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fb:',['refit', 'ptBin=', 'makeTable', 'makePlot', 'pt='])
    except getopt.GetoptError:
        print './rcs.py [-f] [-b <1,2, 3, 4>] [--makeTable] [--pt]'
        print '-f = --refit   -b = --ptBin'  
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            varLabel = 'refit'
        if opt in ('-b', '--ptBin'):
            ptBin = int(arg)
        if opt == '--makeTable':
            isMakeTable = True
        if opt == '--makePlot':
            isMakePlot = True
        if opt == '--pt':
            ptSpectrum = arg

    if ptSpectrum in ('1S', '1s'):
        ptLabel = '1S'
    elif ptSpectrum in ('2S', '2s'):
        ptLabel = '2S'
    elif ptSpectrum in ('3S', '3S'):
        ptLabel = '3S'

    if isMakeTable:
        for varLabel in ('refit', ):#'qValue'):
            rcs, numChibs, chib_par_dict = makeDictRCS(varLabel, ptLabel)
            makeTableRCS(rcs)
            makeTableR(rcs)
            makeTableN_21(numChibs)
            makeTableChibParams(chib_par_dict)
           
  

    if isMakePlot:
        for varLabel in ('refit', ):#'qValue'):
            rcs, numChibs, chib_par_dict = makeDictRCS(varLabel, ptLabel)
            makePlotRCS(rcs)
        
    else:
        
        if ptBin == None:
            ptBin_label = ''
            ptMin = ptBins[0]
            ptMax = ptBins[-1]
            ptBin_print = str(ptBins[0])+' GeV <p_{T#Upsilon} < '+str(ptBins[-1])+' Gev/c'
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            ptMin = ptBins[ptBin-1]
            ptMax = ptBins[ptBin]
            ptBin_print = str(ptBins[ptBin-1])+' GeV <p_{T#Upsilon} < '+str(ptBins[ptBin])+' Gev/c'
           
        Chib1_eff_file = 'txt/eff_chib1'+'_'+varLabel+'_'+'Upsilon2SPt'+'_'+ptLabel+ptBin_label+'.txt'
        Chib2_eff_file = 'txt/eff_chib2'+'_'+varLabel+'_'+'Upsilon2SPt'+'_'+ptLabel+ptBin_label+'.txt'
        numChib_file = 'txt/numChib_'+varLabel+ptBin_label+'.txt'
        sistematici_file = 'txt/sistematici_'+varLabel+ptBin_label+'.txt'
        
        messaggi = {'qValue' : 'Using Q-value + m_PDG, ',
                    'refit' : 'Using kinematic refitted mass, ',
                    }
        
        print messaggi[varLabel]+ptBin_print
        
        
        chib1_eff = Efficiency()
        chib1_eff = chib1_eff.loadFromFile(Chib1_eff_file)
        
        chib2_eff = Efficiency()
        chib2_eff = chib2_eff.loadFromFile(Chib2_eff_file)
        
        numChib = NumChib()
        numChib = numChib.loadFromFile(numChib_file)
        
        sistematici = Sistematici()
        sistematici = sistematici.loadFromFile(sistematici_file)
        
        rcs = RCS(numChib=numChib, chib1_eff=chib1_eff, chib2_eff=chib2_eff, sistematici=sistematici)
        scrPrint(rcs)
        
        
        
    
