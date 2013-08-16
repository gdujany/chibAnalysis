
from ROOT import TMath, TVector,TLorentzVector,TRotation
from math import sqrt,sin,cos

def jpsimuDirections(chiccand,jpsicand, frame='hx'):
    """return two directions: 
       1. direction vector of jpsi in the chic rest frame, wrt to the direction of chic
       2. direction vector of muon in the jpsi rest frame, wrt to the direction of the psi as seen in the chic rest frame"""


    pbeam = 3500
    Mpsi  = 3.097
    Mprot = 0.938

    Ebeam = sqrt(pbeam**2 + Mprot**2)
    
    targ = TLorentzVector(0.,0.,-pbeam,Ebeam)
    beam = TLorentzVector(0.,0., pbeam,Ebeam)


    # chic 4vector in lab frame
    chi = TLorentzVector()
    #chi.SetXYZM(chiccand.px(), chiccand.py(), chiccand.pz(), mass)
    chi.SetXYZM(chiccand.px(), chiccand.py(), chiccand.pz(), chiccand.mass())

    chi_direction = chi.Vect().Unit()
    
    # psi 4vector in lab fram
    psi = TLorentzVector()
    psi.SetXYZM(jpsicand.px(),jpsicand.py(),jpsicand.pz(),jpsicand.mass())
    
    cm_to_chi = -chi.BoostVector()
    chi_to_cm = chi.BoostVector()

    cm_to_psi = -psi.BoostVector()
 
    beam_chi = beam
    beam_chi.Boost(cm_to_chi)         # beam in the chi rest frame

    targ_chi = targ
    targ_chi.Boost(cm_to_chi)         # target in the chi rest frame
 
    beam_direction_chi     = beam_chi.Vect().Unit()
    targ_direction_chi     = targ_chi.Vect().Unit()
    beam_targ_bisec_chi    = (beam_direction_chi - targ_direction_chi).Unit();

    
    psi_chi = psi
    psi_chi.Boost(cm_to_chi)  # psi in the chi rest frame

    psi_direction_chi = psi_chi.Vect().Unit()
  
    # all polarization frames have the same Y axis = the normal to the plane
    # formed by the directions of the colliding hadrons

    Yaxis = ( beam_direction_chi.Cross( targ_direction_chi ) ).Unit()

    # transform(rotation) psi momentum components from polarization axis system
    # to the system with x,y,z axes as in the laboratory

    ChiPolAxis = chi_direction;  # helicity frame
    if frame is 'cs':
        ChiPolAxis =  beam_targ_bisec_chi  

    newZaxis = ChiPolAxis
    newYaxis = Yaxis
    newXaxis = newYaxis.Cross(newZaxis)

    # rotation needed to go to the chi rest frame
    rotation = TRotation()
    rotation.SetToIdentity()
    rotation.RotateAxes(newXaxis, newYaxis, newZaxis);
    rotation.Invert()
    
    psi_chi_rotated = psi_chi.Vect()
                   
    psi_chi_rotated.Transform(rotation); # direction of the psi in the chi rest frame
                                         # relative to direction of chi in the lab


    # now calculate muon direction in the jpsi rest frame wrt chi direction    
    mumass= 0.105

    Yaxis = (ChiPolAxis.Cross(psi_direction_chi)).Unit()

    newZaxis = psi_direction_chi
    newYaxis = Yaxis 
    newXaxis = newYaxis.Cross(newZaxis)

    rotation.SetToIdentity()
    rotation.RotateAxes(newXaxis,newYaxis,newZaxis)
    rotation.Invert()

    #muon in the lab
    lepton = TLorentzVector()
    lepton.SetXYZM(jpsicand.daughter(0).px(),
                   jpsicand.daughter(0).py(),
                   jpsicand.daughter(0).pz(),
                   mumass)

    #muon in the psi rest frame
    lepton_psi = lepton
    lepton_psi.Boost(cm_to_psi)

    lepton_psi_rotated = lepton_psi.Vect()
    lepton_psi_rotated.Transform(rotation)

    return psi_chi_rotated, lepton_psi_rotated


def chooseR(chic_state,helicity=None):
    """ Set values of R and R2 for given helicity, None for unpolarized"""

    R=0
    R2=0
    
    if helicity is None: 
        if chic_state is 1:
            R = 0.6667     
            R2 = 0.0

        if chic_state is 2:
            R= 2./5.
            R2 = 2./5.

    if helicity is 0:        
        R=0.
        R2=0.


    if helicity is 1: 
            R = 1.     
            R2 = 0.
    
    if helicity is 2:
        if chic_state is 1:
            print 'ERROR, chic_1 cannot have helicity=2' 
        R=0.
        R2=1.

    return R,R2


def angDist(jpsidir, leptondir, chic_state, helicity, ):
    """ Arguments:
        jpsidir   : direction of the jpsi in chi rest frame (wrt to chi dir)
        leptondir : direction of mu in the jpsi rest frams (wrt to jpsi dir as
                    seen in the chic rest frame)
        chic_state: 0,1,2
        helicity : 0,1 for chi1, 0,1,2 for chi2, None for unpolarized

        Return probability of finding the state """


    R,R2 = chooseR(chic_state,helicity)

        
    costh_chihe = leptondir.CosTheta()
    phi_chihe   = leptondir.Phi()
    
    cosTH_psi   = jpsidir.CosTheta()
    
    cosTH2_psi = cosTH_psi*cosTH_psi;
    cosTH4_psi = cosTH2_psi*cosTH2_psi;
    costh2_chihe = costh_chihe*costh_chihe;
    sinth2_chihe = 1 - costh2_chihe;
    
    sinTH_psi   = sqrt( 1. -   cosTH2_psi );
    sinth_chihe = sqrt( sinth2_chihe );
    
    sin2TH_psi   = 2.*sinTH_psi*cosTH_psi;
    sin2th_chihe = 2.*sinth_chihe*costh_chihe;
         
    cosphi_chihe = cos( phi_chihe );
    cos2phi_chihe = 2.*cosphi_chihe*cosphi_chihe -1.;

    
    angdistr =0
    PIG = TMath.Pi()
        
    # chic_0 angular distribution
    if chic_state == 0  :

        angdistr = 1. + costh2_chihe
        angdistr *= 3. / ( 64.* PIG*PIG )
        
       


    # chic_1 angular distribution
    if chic_state == 1  :
        
        a2 = 0.
        
        a1 = sqrt( 1. - a2*a2 )   # (a1 taken to be positive)
        
        A0 = sqrt(1./2.) * ( a1 + a2 )
        A1 = sqrt(1./2.) * ( a1 - a2 )
        
        k1 = A1*A1 + 1./2.* R * ( A0*A0 - A1*A1 )
        k2 = ( 1. - 3./2.* R ) * ( A0*A0 - A1*A1 )
        k3 = -A1*A1 + 1./2.* R
        k4 = 1. - 3./2.* R
        k5 = 1./4.* A1*A0 * ( 3.* R - 2. )
        
        angdistr = k1 + k2 * cosTH2_psi + ( k3 + k4 * cosTH2_psi ) * costh2_chihe + k5 * sin2TH_psi * sin2th_chihe * cosphi_chihe 
          
        angdistr *= 9. / ( 64.* PIG*PIG )



    # chic_2 angular distribution
    if  chic_state == 2  :
          
        a3 = 0.
        
        a2 = 0.
        
        a1 = sqrt( 1. - a2*a2 - a3*a3 )  # (a1 taken to be positive)
          
        A0 = sqrt(1./10.)*a1 + sqrt(1./2.)*a2 + sqrt(2./5.)*a3
        A1 = sqrt(3./10.)*a1 + sqrt(1./6.)*a2 - sqrt(8./15.)*a3
        A2 = sqrt(3./5.)*a1  - sqrt(1./3.)*a2 + sqrt(1./15.)*a3
        
        
        
        k1_0 = 1./4.* A0*A0 + 3./8.* A2*A2
        k1_1 = 1./2.* A1*A1 + 1./4.* A2*A2
        k1_2 = 3./8.* A0*A0 + 1./2.* A1*A1 + 1./16.* A2*A2
        k2_0 = -3./2.* A0*A0 + 3.* A1*A1 - 3./4.* A2*A2
        k2_1 = 3./2.* A0*A0 - 3./2.* A1*A1
        k2_2 = -3./4.* A0*A0 + 3./8.* A2*A2
        k3_0 = 9./4.* A0*A0 - 3.* A1*A1 + 3./8.* A2*A2
        k3_1 = -3./2.* A0*A0 + 2.* A1*A1 - 1./4.* A2*A2
        k3_2 = 3./8.* A0*A0 - 1./2.* A1*A1 + 1./16.* A2*A2
        k4_0 = 1./4.* A0*A0 + 3./8.* A2*A2
        k4_1 = -1./2.* A1*A1 + 1./4.* A2*A2
        k4_2 = 3./8.* A0*A0 - 1./2.* A1*A1 + 1./16.* A2*A2
        k5_0 = -3./2.* A0*A0 - 3.*A1*A1 - 3./4.* A2*A2
        k5_1 = 3./2.* A0*A0 + 3./2.* A1*A1
        k5_2 = -3./4.* A0*A0 + 3./8.* A2*A2
        k6_0 = 9./4.* A0*A0 + 3.* A1*A1 + 3./8.* A2*A2
        k6_1 = -3./2.* A0*A0 - 2.*A1*A1 - 1./4.* A2*A2
        k6_2 = 3./8.* A0*A0 + 1./2.* A1*A1 + 1./16.* A2*A2
        k7_0 = -sqrt(6.)/4.* A0*A2
        k7_1 = 0.
        k7_2 = sqrt(6.)/8.* A0*A2
        k8_0 = sqrt(6.)* A0*A2
        k8_1 = -sqrt(6.)/2. *A0*A2
        k8_2 = 0.
        k9_0 = -3.* sqrt(6.)/4.* A0*A2
        k9_1 = sqrt(6.)/2.* A0*A2
        k9_2 = -sqrt(6.)/8.* A0*A2
        k10_0 = sqrt(3.)/4.* A0*A1 + 3.*sqrt(2.)/8.* A1*A2
        k10_1 = -sqrt(3.)/4.* A0*A1
        k10_2 = sqrt(3.)/8.* A0*A1 - 3.*sqrt(2.)/16.* A1*A2
        k11_0 = -3.*sqrt(3.)/4.* A0*A1 - 3.*sqrt(2.)/8.*A1*A2
        k11_1 = sqrt(3.)/2.* A0*A1 + sqrt(2.)/4.* A1*A2
        k11_2 = -sqrt(3.)/8.* A0*A1 - sqrt(2.)/16.* A1*A2
        
        R1 = R
        R0=1.-R1-R2
        
        k1 = R0*k1_0 +R1*k1_1 +R2*k1_2
        k2 = R0*k2_0 +R1*k2_1 +R2*k2_2
        k3 = R0*k3_0 +R1*k3_1 +R2*k3_2
        k4 = R0*k4_0 +R1*k4_1 +R2*k4_2
        k5 = R0*k5_0 +R1*k5_1 +R2*k5_2
        k6 = R0*k6_0 +R1*k6_1 +R2*k6_2
        k7 = R0*k7_0 +R1*k7_1 +R2*k7_2
        k8 = R0*k8_0 +R1*k8_1 +R2*k8_2
        k9 = R0*k9_0 +R1*k9_1 +R2*k9_2
        k10= R0*k10_0+R1*k10_1+R2*k10_2
        k11= R0*k11_0+R1*k11_1+R2*k11_2
        
        angdistr = k1 + k2 * cosTH2_psi + k3 * cosTH4_psi + ( k4 + k5 * cosTH2_psi + k6 * cosTH4_psi ) * costh2_chihe + ( k7 + k8 * cosTH2_psi + k9 * cosTH4_psi ) * sinth2_chihe * cos2phi_chihe + ( k10 + k11 * cosTH2_psi )* sin2TH_psi * sin2th_chihe * cosphi_chihe 
        
        angdistr *= 15. / ( 64.* PIG*PIG )
        

    return angdistr
