# Authors: F. Anders, I. García-Soriano, J. Dolcet (ICCUB)

"""
Galah data for tSNE
"""

import numpy as np
from astropy.io import fits as pyfits
from astropy.table import Table
import scipy

def sqadd(a, b):
    "Add 2 values in quadrature"
    return np.sqrt(a*a + b*b)

class galah(object):
    def __init__(self, tefflogg=True, abundances=True, cannon=True):
        """
        Open a portion of galah data
        """
        hdu = pyfits.open('../data/GALAH_DR2_withSH_feb2020.fits', names=True)
        data = hdu[1].data
        if tefflogg:
            data = data[ (data['teff']>5300) * (data['teff']<6500) *
                         (data['logg']>3) * (data['logg']<5) ]
        if abundances:
            data = data[ (data['flag_o_fe']==0) * (data['flag_na_fe']==0) * 
                         (data['flag_mg_fe']==0) * (data['flag_al_fe']==0) * 
                         (data['flag_si_fe']==0) * (data['flag_k_fe']==0) * 
                         (data['flag_ca_fe']==0) * (data['flag_sc_fe']==0) * 
                         (data['flag_ti_fe']==0) * (data['flag_v_fe']==0) * 
                         (data['flag_cr_fe']==0) * (data['flag_mn_fe']==0) *
                         (data['flag_ni_fe']==0) * (data['flag_cu_fe']==0) * 
                         (data['flag_zn_fe']==0) * (data['flag_y_fe']==0) * 
                         (data['flag_ba_fe']==0)  #* (data['flag_la_fe']==0)
                         ] #(data['e_ba_fe]<1)
        if cannon:
            data = data[ (data['flag_cannon']==0) ]
        self.data = data
        return None

    def create_abundancearray(self, mc=1, age=False,
                           kin=False, feh=True):
        """
        Here, I take only the columns necessary from the lines selected before
        """
        data = self.data
        X        = np.c_[data['fe_h'],data['o_fe'],data['na_fe'],data['mg_fe'],
                         data['al_fe'],data['si_fe'],data['k_fe'],data['ca_fe'],
                         data['sc_fe'],data['ti_fe'],data['v_fe'],data['cr_fe'],
                         data['mn_fe'],data['ni_fe'],data['cu_fe'],data['zn_fe'],
                         data['y_fe'],data['ba_fe'],data['age50']]
        Xerr1    = np.c_[data['e_fe_h'],data['e_o_fe'],data['e_na_fe'],data['e_mg_fe'],
                         data['e_al_fe'],data['e_si_fe'],data['e_k_fe'],data['e_ca_fe'],
                         data['e_sc_fe'],data['e_ti_fe'],data['e_v_fe'],data['e_cr_fe'],
                         data['e_mn_fe'],data['e_ni_fe'],data['e_cu_fe'],data['e_zn_fe'],
                         data['e_y_fe'],data['e_ba_fe'],(data['age84']-data['age16'])*0.5]

        # Take care of the 0.0 uncertainties: forced minimum to 0.03 
        Xerr1[:, :] = np.maximum(Xerr1, 0.03*np.ones(Xerr1.shape))
        Xerr     = np.mean( Xerr1, axis=0)

        if mc > 1:
            Y        = np.zeros(( mc*len(data), 19 ))
            for ii in np.arange(mc):
                Y[ii::mc, :] = scipy.random.normal(loc=X, size=X.shape, scale=Xerr1)                                      
            X = Y
        if not age:
            X = X[:,:-1]; Xerr = Xerr[:-1]
        if kin:
            X    = np.append(X, np.vstack((data['Ulsr'],data['Vlsr'],data['Wlsr'])).T, axis=1)
            Xerr = np.append(Xerr, np.array([10,10,10]).T, axis=0)
        if not feh:
            X = X[:,1:]; Xerr = Xerr[1:]
        # Defining properties
        self.Xerr     =  Xerr
        self.X        =  X
        self.Xnorm    = (X/Xerr[np.newaxis,:])

class apogeedr16_rc(object):
    def __init__(self, path="../data/apogee-rc-DR16.fits", okayflags=True, sample="all", add_astroNN=False):
        # Define useful sub-samples from APOGEE DR16 RC catalogue (Bovy+2014):
        apogee = Table.read(path)
        self.sample = sample
        # First, clean the typical APOGEE flags:
        if okayflags:
            # These are the conditions that all subsets used should fulfil:
            okayflagssnrchi2 = (apogee["SNREV"] > 100) & (apogee["ASPCAP_CHI2"] < 25) & \
                           ("BAD" not in apogee['ASPCAPFLAGS']) & ("NO_ASPCAP" not in apogee['ASPCAPFLAGS']) &\
                           ("TELLURIC" not in apogee['TARGFLAGS']) & ("CLUSTER" not in apogee['TARGFLAGS']) &\
                           ("SERENDIPITOUS" not in apogee['TARGFLAGS']) & ("MASSIVE" not in apogee['TARGFLAGS']) &\
                           ("EMISSION" not in apogee['TARGFLAGS']) & ("BAD" not in apogee['STARFLAGS']) &\
                           ("COMMISSIONING" not in apogee['STARFLAGS']) & ("SUSPECT" not in apogee['STARFLAGS'])
                           #(apogee["GALVR"] >= -15168.6631) & (apogee["GALVT"] >= -20029.2009) & (apogee["GALVZ"] >= -16388.8357) &\
                           #(apogee["GALVR"] <= 14982.2969) & (apogee["GALVT"] <= 20049.5991) & (apogee["GALVZ"] <= 16318.1843)
        # In addition, the APOGEE DR16 element-wise flags should be good. This is how they are enumerated (NEW in DR16!):
        # 0 1  2 3  4  5  6  7 8 9 10 11 12  13  14 15 16 17 18 19 20 | 21 22 23 24 25
        # C CI N O Na Mg Al Si P S K  Ca Ti TiII V  Cr Mn Fe Co Ni Cu | Ge Rb Ce Nd Yb
        # Based on a quick look at the column statistics (so that we don't loose to many stars), 
        # we use all elements up to Cu:
        elemflags = np.sum(apogee['ELEMFLAG'][:, :20], axis=1) == 0

        # Now join all conditions and return what's left.
        self.data = apogee[ okayflagssnrchi2 & elemflags]
        print(len(self.data), "stars in your APOGEE DR16 sample", sample)
        if add_astroNN:
            # Add astroNN data - Read table:
            astronn_rc = Table.read("../data/apogee-rc-DR16_astroNN.fits")
            self.nndata= astronn_rc[okayflagssnrchi2 & elemflags]
        return
    
    def get_ndimspace(self, cn=True, age=False, norm="stdev"):
        """
        Cut out missing data and prepare t-SNE input array
        Optional:
            age: Bool  - include age in the analysis
            cn:  Bool  - include [C/Fe] & [N/Fe], default: True
            norm: str  - normalisation method, default: "stdev"
        """
        # For giants, everything up to Cu is okay
        X        = np.c_[self.data['X_M'][:, :20], self.data['M_H']]
        Xerr     = np.mean(np.c_[self.data['X_M_ERR'][:, :20], self.data['M_H_ERR']], axis=0)
        if not cn:
            X = X[:,2:]; Xerr = Xerr[2:]
        if norm == "hogg2016":
            # Normalise everything by the typical range as in Hogg+2016: 
            Xnorm    = (X/Xerr[np.newaxis,:])
        elif norm == "stdev":
            # Normalise everything by the typical range defined by the std deviation: 
            Xnorm    = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
        elif norm == None:
            # Do not normalise at all:
            self.Xnorm    = X
        else:
            raise ValueError("Please pass a valid 'norm' keyword.")
        self.X, self.Xerr, self.Xnorm = X, Xerr, Xnorm
        return

    def get_ndimspace_H(self, cn=True, age=False, norm="stdev"):
        """
        Cut out missing data and prepare t-SNE input array
        Optional:
            age: Bool  - include age in the analysis
            cn:  Bool  - include [C/H] & [N/H], default: True
            norm: str  - normalisation method, default: "stdev"
        """
        # For giants, everything up to Cu is okay
        X        = self.data['X_H'][:, :20]
        Xerr     = np.mean(self.data['X_H_ERR'][:, :20], axis=0)
        if not cn:
            X = X[:,2:]; Xerr = Xerr[2:]
        if norm == "hogg2016":
            # Normalise everything by the typical range as in Hogg+2016: 
            Xnorm    = (X/Xerr[np.newaxis,:])
        elif norm == "stdev":
            # Normalise everything by the typical range defined by the std deviation: 
            Xnorm    = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
        elif norm == None:
            # Do not normalise at all:
            self.Xnorm    = X
        else:
            raise ValueError("Please pass a valid 'norm' keyword.")
        self.X, self.Xerr, self.Xnorm = X, Xerr, Xnorm
        return

    def get_ndimspace_H_astroNN(self, cn=True, age=False, norm="stdev"):
        """
        Cut out missing data and prepare t-SNE input array
        Optional:
            age: Bool  - include age in the analysis
            cn:  Bool  - include [C/H] & [N/H], default: True
            norm: str  - normalisation method, default: "stdev"
        """
        # For giants, everything up to Cu is okay
        X        = np.c_[self.nndata['C_H'], self.nndata['CI_H'], self.nndata['N_H'], self.nndata['O_H'], 
                         self.nndata['NA_H'], self.nndata['MG_H'], self.nndata['AL_H'], self.nndata['SI_H'], 
                         self.nndata['P_H'], self.nndata['S_H'], self.nndata['K_H'], self.nndata['CA_H'], 
                         self.nndata['TI_H'], self.nndata['TIII_H'], self.nndata['V_H'], self.nndata['CR_H'], 
                         self.nndata['MN_H'], self.nndata['FE_H'], self.nndata['CO_H'], self.nndata['NI_H'] 
                         ]
        Xerr     = np.mean(np.c_[self.nndata['C_H_ERR'], self.nndata['CI_H_ERR'], self.nndata['N_H_ERR'], self.nndata['O_H_ERR'], 
                         self.nndata['NA_H_ERR'], self.nndata['MG_H_ERR'], self.nndata['AL_H_ERR'], self.nndata['SI_H_ERR'], 
                         self.nndata['P_H_ERR'], self.nndata['S_H_ERR'], self.nndata['K_H_ERR'], self.nndata['CA_H_ERR'], 
                         self.nndata['TI_H_ERR'], self.nndata['TIII_H_ERR'], self.nndata['V_H_ERR'], self.nndata['CR_H_ERR'], 
                         self.nndata['MN_H_ERR'], self.nndata['FE_H_ERR'], self.nndata['CO_H_ERR'], self.nndata['NI_H_ERR'] 
                         ], axis=0)
        if not cn:
            X = X[:,2:]; Xerr = Xerr[2:]
        if norm == "hogg2016":
            # Normalise everything by the typical range as in Hogg+2016: 
            Xnorm    = (X/Xerr[np.newaxis,:])
        elif norm == "stdev":
            # Normalise everything by the typical range defined by the std deviation: 
            Xnorm    = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
        elif norm == None:
            # Do not normalise at all:
            self.Xnorm    = X
        else:
            raise ValueError("Please pass a valid 'norm' keyword.")
        self.X, self.Xerr, self.Xnorm = X, Xerr, Xnorm
        return

    def get_umap_tsne_colours(self, p=None, lr=None, nn=None, md=None, metric="euclidean", version="_H"):
        """
        Get the umap & t-SNE results for the optimal hyperparameters,
        and also the arrays used to colour the maps.
        """
        #Read umap/t-SNE results table
        results = Table.read("../data/dimred_results/apogee_rc"+version+"_dimred_hyperparametertest.fits")
        # Get the relevant columns
        self.Xp = results["X_PCA"]
        self.Yp = results["Y_PCA"]
        self.Xt = results["X_tSNE_"+metric+"_p"+str(p) + "_lr"+str(lr)]
        self.Yt = results["Y_tSNE_"+metric+"_p"+str(p) + "_lr"+str(lr)]
        self.Xu = results["X_umap_"+metric+"_nn"+str(nn) + "_md"+str(md)]
        self.Yu = results["Y_umap_"+metric+"_nn"+str(nn) + "_md"+str(md)]
        # And now the colours
        data  = self.data
        self.colors   = [data['FE_H'], data['MG_FE'], data['MG_FE']-data['SI_FE'], data['MG_FE']-(data['FE_H']+data['O_FE']),
                        data['N_FE']-data['O_FE'],data['C_FE']-data['N_FE'],data['AL_FE']-data['MG_FE'], data['P_FE']-(data['SI_FE']),
                        data['K_FE']-(data['S_FE']),data['S_FE']-(data['SI_FE']),data['MN_FE']-data['CR_FE'],data['CU_FE']-data['NI_FE'],
                        data['TEFF'],data['RC_GALR'], data['RC_GALZ'],np.log10(sqadd(data['VSCATTER'],data['VERR'])),np.log10(data['SNREV']),
                        data['GALVR'],data['GALVT'],data['GALVZ']]
        self.titles   = [r'$\rm [Fe/H]$', r'$\rm [Mg/Fe]$', r'$\rm [Mg/Si]$', r'$\rm [Mg/O]$',
                      r'$\rm [N/O]$', r'$\rm [C/N]$', r'$\rm [Al/Mg]$', r'$\rm [P/Si]$',
                      r'$\rm [K/S]$', r'$\rm [S/Si]$', r'$\rm [Mn/Cr]$', r'$\rm [Cu/Ni]$',
                      r'$T_{\rm eff}$', r'$R$', r'$Z$', r'log $\sigma_{RV}$',#r'log $S/N$',
                      r'SNREV', r'$v_R$', r'$v_T$', r'$v_Z$']
        return

    def get_umap_subsets(self, nn=100, md=0.1, **kwargs):
        """
        Get the names and indices of the t-sne-defined subsets
        """
        # First get umap results:
        results = Table.read("../data/dimred_results/apogee_rc_dimred_hyperparametertest.fits")
        self.Xu = results["X_umap_euclidean_nn"+str(nn) + "_md"+str(md)]
        self.Yu = results["Y_umap_euclidean_nn"+str(nn) + "_md"+str(md)]
        
        # Now run HDBSCAN to define the subsets
        import hdbscan
        clusterer = hdbscan.HDBSCAN(**kwargs)
        clusterer.fit( np.vstack((self.Xu, self.Yu)).T )
        self.classcol = clusterer.labels_
        self.classprob= clusterer.probabilities_
        self.subsets  = np.unique(clusterer.labels_)
        #self.classcol= np.char.rstrip(self.data["tsne_class_teffcut40"],b' ')#.decode('utf8').strip()
        #self.subsets = ["thin", "thick1", "thick2", "thick3", "thick4",
        #           "mpthin", "mpthintrans", "smr", "t4trans", "youngthin",
        #           "debris1", "debris2", "debris3", "debris4", "debris5", 
        #           "smr2", "t2trans1", "highTi","lowMg","highAlMg?"]
        self.names   = ["", "", "", "",
                   "", "", "Transition group", "", "",
                   "Young local disc", "", "", "[s/Fe]-enhanced", "", "", r"", "Debris candidate", 
                   r"Extreme-Ti star", r"Low-[Mg/Fe] star", "High-[Al/Mg] star"]
        self.Xcoords = [10, 11, 4.5, -12,  18, -31, 22, 26,-22.5, -14, -2, -25]
        self.Ycoords = [5.5,.5,  -2, -4,   6,  0,   1.5, -.5, -7, -2, -6, 14]
        self.fsize   = [20 , 16,  12, 12,  15,  13, 11, 11, 11, 11, 11, 11]
        self.sym = ["o", "v", "^", ">", "<", "s", "o", "*", "<", "o",
                    "h", "d", "D", "v", "p", "*", "D", "p", "s", "8"]
        self.al  = [.6, .8, .8, .8, .8, .8, .8,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        self.lw  = [0,.5,.5,.5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5]
        self.size= [7,12,12,12,12,15,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]
        self.col = ["grey", "m", "hotpink", "crimson", "r",
                    "g", "brown", "orange", "gold", "k",
                    "yellow", 
                    "gold", "lime", "k", "royalblue"]

