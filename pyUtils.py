from __future__ import division
import pickle
from math import sqrt, log10, floor
from ROOT import TFile, TTree, RooRealVar, RooArgSet, RooDataSet

def roundPair(val, err, sig=2):
    try:
        cfr = sig-int(floor(log10(err)))-1
    except ValueError:
        cfr = 2
    #return round(val, cfr), round(err, cfr)
    try:
        return ('{:.'+str(cfr)+'f}').format(val), ('{:.'+str(cfr)+'f}').format(err)
    except ValueError:
        if cfr > 0:
            return str(round(val, cfr)), str(round(err, cfr))
        else:
            return str(int(round(val, cfr))), str(int(round(err, cfr)))
    #     return ('{:.2f}').format(val), ('{:.2f}').format(err)

def roundList(ll, sig=1):
    cfr_list = sorted([sig-int(floor(log10(abs(err))))-1 for err in ll if err !=0.])
    cfr = cfr_list[-1]
    return [('{:.'+str(cfr)+'f}').format(val) for val in ll]

def roundDict(dd, sig=1, cfr_fixed=None):
    cfr_list = sorted([sig-int(floor(log10(abs(err))))-1 for err in dd.values() if err !=0.])
    cfr = cfr_list[-1]
    if cfr_fixed: cfr = cfr_fixed
    return dict([(key,('{:.'+str(cfr)+'f}').format(val)) for key,val in dd.items()])


def makeRooDataset(inputfile_name):
    chibTree_name = 'rootuple/chibTree'
    inputfile = TFile.Open(inputfile_name,"READ")
    tree = TTree()
    inputfile.GetObject(chibTree_name, tree)
    
    invm1S = RooRealVar("invm1S", "invm1S", 9.7, 10.1)#9.5, 11.5
    invm2S = RooRealVar("invm2S", "invm2S", 9.5, 13.0)
    invm3S = RooRealVar("invm3S", "invm3S", 9.5, 13.0)
    chib_mass = RooRealVar("chib_mass", "chib_mass", 5., 15.0)
    chib_pt = RooRealVar("chib_pt", "chib_pt", 0., 300.)
    chib_eta = RooRealVar("chib_eta", "chib_eta", -5., 5.)
    chib_phi = RooRealVar("chib_phi", "chib_phi", -3.2, 3.2)
    dimuon_mass = RooRealVar("dimuon_mass","dimuon_mass", 8.0, 12.0)
    dimuon_rapidity = RooRealVar("dimuon_rapidity", "dimuon_rapidity", -5.0, 5.0)
    dimuon_pt = RooRealVar("dimuon_pt","dimuon_pt", 0.0, 100.0)
    photon_eta = RooRealVar("photon_eta","photon_eta", -5.0, 5.0)
    photon_pt = RooRealVar("photon_pt","photon_pt", 0.0, 100.0)
    ctpv = RooRealVar("ctpv","ctpv", -5.0, 5.0)
    ctpv_error = RooRealVar("ctpv_err","ctpv_err", -5.0, 5.0)
    pi0_abs_mass = RooRealVar("pi0_abs_mass","pi0_abs_mass", 0.0, 2.0)
    Y1S_nsigma = RooRealVar("Y1S_nsigma","Y1S_nsigma",0.0,30.0)
    Y2S_nsigma = RooRealVar("Y2S_nsigma","Y2S_nsigma",0.0,30.0)
    Y3S_nsigma = RooRealVar("Y3S_nsigma","Y3S_nsigma",0.0,35.0)
    conv_vertex = RooRealVar("conv_vertex", "conv_vertex", 0.0, 70.0)
    dz = RooRealVar("dz","dz", 0., 0.6)
    numPrimaryVertices = RooRealVar("numPrimaryVertices","numPrimaryVertices",0,60)
    rf1S_chib_mass = RooRealVar("rf1S_chib_mass","rf1S_chib_mass",9.5,60.)
    probFit1S = RooRealVar("probFit1S","probFit1S",0.,1.)
    rf1S_chib_pt = RooRealVar("rf1S_chib_pt", "rf1S_chib_pt", 0., 300.)
    rf1S_chib_eta = RooRealVar("rf1S_chib_eta", "rf1S_chib_eta", -5., 5.)
    rf1S_dimuon_mass = RooRealVar("rf1S_dimuon_mass","rf1S_dimuon_mass", 8.0, 12.0)
    rf1S_dimuon_rapidity = RooRealVar("rf1S_dimuon_rapidity", "rf1S_dimuon_rapidity", -5.0, 5.0)
    rf1S_dimuon_pt = RooRealVar("rf1S_dimuon_pt","rf1S_dimuon_pt", 0.0, 230.0)
    rf1S_photon_eta = RooRealVar("rf1S_photon_eta","rf1S_photon_eta", -5.0, 5.0)
    rf1S_photon_pt = RooRealVar("rf1S_photon_pt","rf1S_photon_pt", 0.0, 100.0)
    rf2S_chib_mass = RooRealVar("rf2S_chib_mass","rf2S_chib_mass",9.5,60.)
    probFit2S = RooRealVar("probFit2S","probFit2S",0.,1.)
    rf3S_chib_mass = RooRealVar("rf3S_chib_mass","rf3S_chib_mass",9.5,60.)
    probFit3S = RooRealVar("probFit3S","probFit3S",0.,1.)

    muonP_pt = RooRealVar("muonP_pt","muonP_pt",0.,100.)
    muonM_pt = RooRealVar("muonM_pt","muonM_pt",0.,100.)
    muonP_eta = RooRealVar("muonP_eta","muonP_eta",-5.,5.)
    muonM_eta = RooRealVar("muonM_eta","muonM_eta",-5.,5.)

    dataArgSet = RooArgSet(invm1S, invm2S, invm3S, dimuon_mass, dimuon_rapidity, dimuon_pt, photon_eta, photon_pt, ctpv)
    dataArgSet.add(RooArgSet(chib_mass, chib_pt, chib_eta, chib_phi))
    dataArgSet.add(RooArgSet(ctpv_error, pi0_abs_mass, Y1S_nsigma, Y2S_nsigma, Y3S_nsigma, conv_vertex, dz, numPrimaryVertices))
    dataArgSet.add(RooArgSet(rf1S_chib_mass, probFit1S, rf2S_chib_mass, probFit2S, rf3S_chib_mass, probFit3S))
    dataArgSet.add(RooArgSet(rf1S_chib_pt, rf1S_chib_eta, rf1S_dimuon_mass, rf1S_dimuon_rapidity, rf1S_dimuon_pt, rf1S_photon_eta, rf1S_photon_pt))
    dataArgSet.add(RooArgSet(muonP_pt, muonM_pt, muonP_eta, muonM_eta))
    dataSet = RooDataSet("chibds","Chib RooDataSet", tree, dataArgSet)
    return dataSet



def dNdpT(pt, pt2m, beta) :
        return pt * (1 + pt**2/(pt2m*(beta-2)))**(-beta)

def givePtweight(pt, outPt = '2S', inPt='2S'):
    
    if inPt == outPt:
        return 1
    
    pt2m = {'1S' : 51.,
            '2S' : 75.,
            '3S' : 75.4}

    beta = {'1S' : 3.24,
            '2S' : 2.84,
            '3S' : 3.31}

    return dNdpT(pt, pt2m[outPt], beta[outPt])/dNdpT(pt, pt2m[inPt], beta[inPt])
    


class Cuts:
    """Class where I implement Cuts that can be used in the various fits and in computing efficiences"""
    def __init__(self, photon_pt_cut = 0., photon_eta_cut = 1, upsilon_pt_lcut = 7, upsilon_pt_hcut = 40, upsilon_rapidity_cut = 1.25,  pi0_abs_mass_cut = 0.,  dz_cut = 0.1, muon_pt_cut = 2.5, muon_eta_cut = 2.4, varFit = 'qValue',str_cut = ''):
        self.photon_pt_cut = photon_pt_cut
        self.photon_eta_cut = photon_eta_cut
        self.upsilon_pt_lcut = upsilon_pt_lcut
        self.upsilon_pt_hcut = upsilon_pt_hcut
        self.upsilon_rapidity_cut = upsilon_rapidity_cut
        self.pi0_abs_mass_cut = pi0_abs_mass_cut
        self.dz_cut = dz_cut
        self.muon_pt_cut = muon_pt_cut
        self.muon_eta_cut = muon_eta_cut
        self.str_cut = str_cut
        if(varFit == 'refit'): self.varLabel = 'rf1S_'
        else: self.varLabel = ''
       

    def __str__(self):
        #return "probFit1S > 0.02  && Y1S_nsigma < 2.5 && {varLabel}photon_pt > {photon_pt_cut} && abs({varLabel}photon_eta) < {photon_eta_cut} && {varLabel}dimuon_pt > {dimuon_pt_lcut} && {varLabel}dimuon_pt < {dimuon_pt_hcut} && abs({varLabel}dimuon_rapidity) < {dimuon_rapidity_cut} && pi0_abs_mass > {pi0_abs_mass_cut} &&  abs(dz) < {dz_cut}".format(photon_pt_cut=self.photon_pt_cut, photon_eta_cut=self.photon_eta_cut, dimuon_pt_lcut=self.upsilon_pt_lcut, dimuon_pt_hcut=self.upsilon_pt_hcut, dimuon_rapidity_cut=self.upsilon_rapidity_cut, pi0_abs_mass_cut=self.pi0_abs_mass_cut, dz_cut=self.dz_cut, varLabel = self.varLabel)+self.str_cut
        return "probFit1S > 0.02  && Y1S_nsigma < 2.5 && {varLabel}photon_pt > {photon_pt_cut} && abs({varLabel}photon_eta) < {photon_eta_cut} && {varLabel}dimuon_pt > {dimuon_pt_lcut} && {varLabel}dimuon_pt < {dimuon_pt_hcut} && abs({varLabel}dimuon_rapidity) < {dimuon_rapidity_cut} && pi0_abs_mass > {pi0_abs_mass_cut} &&  abs(dz) < {dz_cut} && muonP_pt > {muon_pt_cut} && muonM_pt > {muon_pt_cut} && abs(muonP_eta) < {muon_eta_cut} && abs(muonM_eta) < {muon_eta_cut}".format(photon_pt_cut=self.photon_pt_cut, photon_eta_cut=self.photon_eta_cut, dimuon_pt_lcut=self.upsilon_pt_lcut, dimuon_pt_hcut=self.upsilon_pt_hcut, dimuon_rapidity_cut=self.upsilon_rapidity_cut, pi0_abs_mass_cut=self.pi0_abs_mass_cut, dz_cut=self.dz_cut, muon_pt_cut=self.muon_pt_cut, muon_eta_cut=self.muon_eta_cut, varLabel = self.varLabel)+self.str_cut
  
    def saveToFile(self, fileName = "cuts.txt"):
        pickle.dump(self, open( fileName, "wb" ))
        
    def loadFromFile(self,fileName = "cuts.txt"):
        return pickle.load( open( fileName, "rb" ) )
        

class CB_parameters:
    """Class where I implemen methods to save the Crystal ball parameters found in the MC-Fit and to be used in the data-Fit"""
    def __init__(self, mean=0, sigma=1, n1=3.8, n2=2.0, a1=1, a2=1, s_mean=0, s_sigma=0, s_n1=0, s_n2=0, s_a1=0, s_a2=0, chiSquare=0):
        self.mean = mean
        self.s_mean = s_mean
        self.sigma = sigma
        self.s_sigma = s_sigma
        self.n1 = n1
        self.s_n1 = s_n1
        self.a1 = a1
        self.s_a1 = s_a1
        self.n2 = n2
        self.s_n2 = s_n2
        self.a2 = a2
        self.s_a2 = s_a2
        self.chiSquare = chiSquare

    def saveToFile(self, fileName = "CB_parameters.txt"):
           pickle.dump(self, open( fileName, "wb" ))
            
    def loadFromFile(self,fileName = "CB_parameters.txt"):
        return pickle.load( open( fileName, "rb" ) )


class Efficiency:
    '''
    Class to store and evaluate information about a single efficiency (e1 or e2)
    '''
    def __init__(self, eff=None, s_eff=None, N_gen=None, N_sel=None, s2_N_gen=None, s2_N_sel=None, cov_N_gen_sel=None):
        if((eff==None or s_eff==None) and (N_gen==None or N_sel==None) and (s2_N_gen==None or s2_N_sel==None or cov_N_gen_sel==None)):
            self.eff = 0
            self.s_eff = 1
        elif(eff==None and s_eff==None and (s2_N_gen==None or s2_N_sel==None or cov_N_gen_sel==None)):
            self.eff = N_sel/N_gen
            self.s_eff = sqrt( self.eff * (1-self.eff) / N_gen )
        elif(eff==None or s_eff==None):
            self.eff = N_sel/N_gen
            if N_sel != 0:
                self.s_eff = N_sel/N_gen * sqrt( s2_N_sel/N_sel**2 + s2_N_gen/N_gen**2 - 2*cov_N_gen_sel/(N_sel*N_gen) )
            else:
                self.s_eff = 0
        else:
            self.eff = eff
            self.s_eff = s_eff
  
    def s_rel_eff(self):
        return self.s_eff/self.eff 

    def __str__(self):
        return 'efficiency = '+str(self.eff)+'\n'+\
            'error = '+str(self.s_eff)+'\n'+\
            'relative error = '+str(self.s_rel_eff())

    def saveToFile(self, fileName = "efficiency.txt"):
        pickle.dump(self, open( fileName, "wb" ))
        
    def loadFromFile(self,fileName = "efficiency.txt"):
        return pickle.load( open( fileName, "rb" ) )
 


def dictEffArgs(gen_w, gen_w2, sel_w, sel_w2, N_noWeight):
    return dict(N_gen = gen_w,
                N_sel = sel_w,
                s2_N_gen = gen_w2 - gen_w**2/N_noWeight,
                s2_N_sel = sel_w2 - sel_w**2/N_noWeight,
                cov_N_gen_sel = sel_w2 - sel_w*gen_w/N_noWeight)



class RatioEff_Polar_Unpolar:
    '''
    Class to store and evaluate information about ratio efficiency in the polarized case over the efficiency in the unpolarized case (es. e1/e1_unpolarized or e2/e2_unpolarized)
    '''
    def __init__(self, gen_w, gen_w2, gen_wu, sel_w, sel_w2, sel_wu, u_gen_w, u_gen_w2, u_sel_w, u_sel_w2, N_noWeight, label='Default'):
        self.gen_w = gen_w
        self.gen_w2 = gen_w2
        self.gen_wu = gen_wu
        self.sel_w = sel_w
        self.sel_w2 = sel_w2
        self.sel_wu = sel_wu
        self.u_gen_w = u_gen_w
        self.u_gen_w2 = u_gen_w2
        self.u_sel_w = u_sel_w
        self.u_sel_w2 = u_sel_w2
        self.N_noWeight = N_noWeight
        self.label = label

    def eff(self):
        ''' Return efficiency polarized'''
        return Efficiency(**dictEffArgs(gen_w=self.gen_w, gen_w2=self.gen_w2, sel_w=self.sel_w, sel_w2=self.sel_w2, N_noWeight=self.N_noWeight))

    def u_eff(self):
        ''' Return efficiency unpolarized'''
        return Efficiency(**dictEffArgs(gen_w=self.u_gen_w, gen_w2=self.u_gen_w2, sel_w=self.u_sel_w, sel_w2=self.u_sel_w2, N_noWeight=self.N_noWeight))

    def rel_cov_pol_unpol(self):
        ''' Return relative covariance between polarized and unpolarized, 
        s_ab/(ab), exacly the thing needed in the error propagation formula'''
        cov_gen_ugen = self.gen_wu - self.gen_w*self.u_gen_w/self.N_noWeight
        cov_sel_usel = self.sel_wu - self.sel_w*self.u_sel_w/self.N_noWeight
        cov_gen_usel = self.sel_wu - self.gen_w*self.u_sel_w/self.N_noWeight
        cov_ugen_sel = self.sel_wu - self.sel_w*self.u_gen_w/self.N_noWeight

        return cov_sel_usel/(self.sel_w * self.u_sel_w) + cov_gen_ugen/(self.gen_w * self.u_gen_w) - cov_gen_usel/(self.gen_w * self.u_sel_w) - cov_ugen_sel/(self.sel_w * self.u_gen_w)

    def CalcRatioEff(self):
        eff_pol = self.eff()
        eff_unpol = self.u_eff()
        return eff_pol.eff / eff_unpol.eff

    def CalcRatioEffError(self):
        eff_pol = self.eff()
        eff_unpol = self.u_eff()
        if self.gen_w == self.u_gen_w: #The unpolarized case, the error is 0,
            return 0 # if I don't do this I can have ~0 but negative resulting in an error because of sqrt
        else:
            return self.CalcRatioEff()*sqrt((eff_pol.s_eff/eff_pol.eff)**2 + (eff_unpol.s_eff/eff_unpol.eff)**2 - 2*self.rel_cov_pol_unpol())

        



class NumChib:
    
    def __init__(self, numChib=1, ratio_21=1, s_numChib=0, s_ratio_21=0, numBkg = 1, s_numBkg = 0, corr_NB=0, corr_NR=0, name='numChib', q0=0, alpha=0, beta=0, s_q0=0, s_alpha=0, s_beta=0, chiSquare=1):
        self.numChib = numChib
        self.ratio_21 = ratio_21
        self.numBkg = numBkg
        self.s_numChib = s_numChib
        self.s_ratio_21 = s_ratio_21
        self.s_numBkg = s_numBkg
        self.corr_NB = corr_NB
        self.cov_NB = corr_NB*s_numChib*s_numBkg # I pass from correlation to covariance
        self.cov_NR = corr_NR*s_numChib*s_ratio_21 # I pass from correlation to covariance
        self.name = name

        self.numChib1 = numChib/(1+ratio_21)
        self.numChib2 = numChib/(1+1/ratio_21)
        self.s_numChib1 = sqrt( s_numChib**2/(1+ratio_21)**2 + numChib**2/(1+ratio_21)**4*s_ratio_21**2 - 2*numChib/(1+ratio_21)**3*self.cov_NR)
        self.s_numChib2 = sqrt( s_numChib**2/(1+1/ratio_21)**2 + numChib**2/(1+ratio_21)**4*s_ratio_21**2 + 2*numChib/( (1+ratio_21)**2 * (1+1/ratio_21) )*self.cov_NR)
        
        self.q0 = q0
        self.s_q0 = s_q0
        self.alpha = alpha
        self.s_alpha = s_alpha
        self.beta = beta
        self.s_beta = s_beta
        self.chiSquare = chiSquare
        

    def saveToFile(self, fileName = 'numChib.txt'):
        pickle.dump(self, open( fileName, "wb" ))
        
    def loadFromFile(self,fileName = 'numChib.txt'):
        return pickle.load( open( fileName, "rb" ) )

    def calcRatio(self):
        return self.ratio_21

    def calcRatioError(self):
        return self.s_ratio_21

    

    def calcSignificance(self):
        Sig = self.numChib/sqrt(self.numChib + self.numBkg)
        return Sig

    def calcSignificanceError(self):
        S = self.numChib
        s_S = self.s_numChib
        B = self.numBkg
        s_B = self.s_numBkg
        s_SB = self.cov_NB
        
        return sqrt(((S+2*B)**2 * s_S**2 + S**2 * s_B**2 - 2 * S * (S+2*B) * s_SB)/(4 * (S+B)**3))


class Sistematici:
    '''
    Contains the relative sistematics uncertentains coming from various sources
    All the errors are relative (fraction) to value
    '''
    def __init__(self, value=0., statistic=0., ratioEff=0., dscbParam=0., funzParam=0., ptDistro=0., polar_hx=(0.,0.), polar_cs =(0.,0.)):
        self.value = value
        self.statistic = statistic
        self.ratioEff = ratioEff
        self.dscbParam = dscbParam
        self.funzParam = funzParam
        self.ptDistro = ptDistro
        self.polar_hx = polar_hx
        self.polar_cs =polar_cs

    def s_rel_N(self):
        pass

    def s_rel_eff(self):
        pass

    def s_rel_total(self):
        return sqrt( self.ratioEff**2 + self.dscbParam**2 + self.ptDistro**2 + self.funzParam**2)
    
    def saveToFile(self, fileName = 'sistematici.txt'):
        pickle.dump(self, open( fileName, "wb" ))
        
    def loadFromFile(self,fileName = 'sistematici.txt'):
        return pickle.load( open( fileName, "rb" ) )

BF_chib1 = 33.9 #35.
s_BF_chib1 = 2.2 #8.
BF_chib2 = 19.1 #22.
s_BF_chib2 = 1.2 #4.
s_rel_BR_rcs = sqrt( (s_BF_chib1/BF_chib1)**2 + (s_BF_chib2/BF_chib2)**2 )

class RCS:
    '''
    Contains ratio of cross sections
    '''
       

    def __init__(self, numChib, chib1_eff, chib2_eff, sistematici):
        self.numChib = numChib
        self.chib1_eff = chib1_eff
        self.chib2_eff = chib2_eff
        self.sistematici = sistematici
        
    def ratioEfficiency(self):
        return self.chib1_eff.eff/self.chib2_eff.eff

    def s_ratioEfficiency(self):
        return self.ratioEfficiency() * sqrt( (self.chib1_eff.s_eff/self.chib1_eff.eff)**2 +  (self.chib2_eff.s_eff/self.chib2_eff.eff)**2 )

    def ratioN(self):
        return self.numChib.calcRatio()

    def s_ratioN(self):
        return self.numChib.calcRatioError()

    def R(self):
        return self.ratioN() * self.ratioEfficiency()

    def s_R(self):
        return self.R() * self.s_ratioN()/self.ratioN()

    def rcs(self):
        return self.R() * (BF_chib1/BF_chib2)

    def s_stat(self):
        return self.rcs() * self.s_ratioN()/self.ratioN()

    def s_rel_sist(self):
        return self.sistematici.s_rel_total()

    def s_rel_polar_hx(self):
        return self.sistematici.polar_hx

    def s_rel_polar_cs(self):
        return self.sistematici.polar_cs

    def s_BR(self):
        return self.rcs() * s_rel_BR_rcs 


############## CIMITERO CODICE VECCHIO ################

class OldNumChib: # it takes the two numbers N1 and N2 instead of total number and ratio, is now obsolete
    
    def __init__(self, numChib1=1, numChib2=1, s_numChib1=0, s_numChib2=0, numBkg = 1, s_numBkg = 0, corr_12=0, corr_1b=0, corr_2b=0 ,name='numChib'):
        self.numChib1 = numChib1
        self.numChib2 = numChib2
        self.numBkg = numBkg
        self.s_numChib1 = s_numChib1
        self.s_numChib2 = s_numChib2
        self.s_numBkg = s_numBkg
        self.corr_12 = corr_12
        self.corr_1b = corr_1b
        self.corr_2b = corr_2b
        self.cov_12 = corr_12*s_numChib1*s_numChib2 # I pass from correlation to covariance
        self.cov_1b = corr_1b*s_numChib1*s_numBkg
        self.cov_2b = corr_2b*s_numChib2*s_numBkg
        self.name = name

    def saveToFile(self, fileName = 'numChib.txt'):
        pickle.dump(self, open( fileName, "wb" ))
        
    def loadFromFile(self,fileName = 'numChib.txt'):
        return pickle.load(( open( fileName, "rb" ) ))

    def calcRatio(self):
        return self.numChib2/self.numChib1

    def calcRatioError(self):
        return self.calcRatio() * sqrt( (self.s_numChib2/self.numChib2)**2 + (self.s_numChib1/self.numChib1)**2 - 2*self.cov_12/(self.numChib1*self.numChib2))

    def calcSignificance(self):
        Sig = (self.numChib1 + self.numChib2)/sqrt(self.numChib1 + self.numChib2 + self.numBkg)
        return Sig

    def calcSignificanceError(self):
        x = self.numChib1
        y = self.numChib1
        z = self.numBkg
        a = (x + y + 2*z)/(2*(x+y+z)**(3/2))
        b = -(x + y)/(2*(x+y+z)**(3/2))
        return sqrt(a**2 * (self.s_numChib1**2 + self.s_numChib2**2 + self.cov_12) + b*(b*self.s_numBkg**2 + a*(self.cov_1b + self.cov_2b) ) ) 
