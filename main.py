
# Begin main function
from importlib import reload
from matplotlib.pyplot import figure
from numpy import linspace,logspace,log10
import CRUM as c

reload(c.main)
crm=c.main.CRUM()

ne=1e13

if 1==1:
    l=  [   [5e-2,  1],
            [5e-3,  1],
            [5e-3,  1.5],
            [5e-3,  2],
            [5e-3,  3],
            [5e-4,  3],
            [5e-4,  4],
            [5e-4,  5],
            [5e-5,  5],
            [5e-5,  10],
            [5e-5,  15],
            [5e-5,  30],
        ]
else:
    l=[[5e-3,1]]





''' TIME-EVOLUTION ON PLASMA BACKGROUND'''
if 1==1:

    for i in l:
        ''' GREENLAND VS FULL CRM '''
        if 1==0:
            fig1=crm.plot_nt(i[0],i[1],ne,Nq=True,labelapp=' full')
            crm.plot_nt(i[0],i[1],ne,gl=True,linestyle='-.',linewidth=3,fig=fig1,labelapp=' gl',pretitle='Gl (--) vs Full CRM (-), ', savename='t-dep-box_gl_v_full_te{}_t{}'.format(i[1],i[0]).replace('.','_'))
            ''' UEDGE VS CRUM '''
        if 1==0:
            ue=c.main.CRUM('/input/UE.dat')
            nuclei=False
            fig2=crm.plot_nt(i[0],i[1],ne,gl=True,labelapp=' CRUM',linewidth=3,nuclei=nuclei)
            ue.plot_nt(i[0],i[1],ne,gl=True,linestyle='-.',fig=fig2,labelapp=' UE',pretitle='UEDGE (--) vs. CRUM (-), ', savename='t-dep-box_ue_v_crum_te{}_t{}'.format(i[1],i[0]).replace('.','_'),nuclei=nuclei)
            #print('T={} eV, tau_m={:.2E} s, tau_m,UE={:.2E} s'.format(i[1],1/abs(crm.crm.gl_crm(i[1],ne,matrices=True)[2][1]),1/abs(ue.crm.gl_crm(i[1],ne,matrices=True)[2][1])))

            ''' MAR VS NO MAR '''
        if 1==0:
            noMAR=c.main.CRUM('/input/CRUM_no2vl2_no2vl3.dat')
            noH2plus=c.main.CRUM('/input/CRUM_no2vl2.dat')
            noHminus=c.main.CRUM('/input/CRUM_no2vl3.dat')
            fig3=noMAR.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no MAR', linewidth=3,linestyle='-')
            noH2plus.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no H2+', linewidth=3,fig=fig3,linestyle='-.')
            noHminus.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no H-', linewidth=3,fig=fig3,linestyle='--')
            ue.plot_nt(i[0],i[1],ne,gl=True,fig=fig3,labelapp=' UE',pretitle='UEDGE (--) vs. CRUM (-), ', savename='t-dep-box_MAR_te{}_t{}'.format(i[1],i[0]).replace('.','_'),linestyle=':')
            
            ''' MAR VS MAD VS MAI '''
        if 1==1:
            noMAD=c.main.CRUM('/input/CRUM_no_MAD.dat')
            noMAR=c.main.CRUM('/input/CRUM_no_MAR.dat')
            noMAI=c.main.CRUM('/input/CRUM_no_MAI.dat')
            noMA=c.main.CRUM('/input/CRUM_no_MA.dat')
            nuclei=False
            fig3=crm.plot_nt(i[0],i[1],ne,gl=True,labelapp=' all on', linewidth=3,linestyle='-',nuclei=nuclei)
            #fig3=noMAR.plot_nt(i[0],i[1],ne,gl=True,labelapp=' all on', linewidth=3,linestyle='-',nuclei=nuclei)
            noMAR.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no MAR', linewidth=3,fig=fig3,linestyle='-.',nuclei=nuclei)
            noMAD.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no MAD', linewidth=3,fig=fig3,linestyle=':',nuclei=nuclei)
            noMAI.plot_nt(i[0],i[1],ne,gl=True,linestyle='--',fig=fig3,labelapp=' no MAI',pretitle='MAR MAD MAI, ', savename='t-dep-box_ma_te{}_t{}'.format(i[1],i[0]).replace('.','_'),nuclei=nuclei, color=['b','b'])
            #noMAI.plot_nt(i[0],i[1],ne,gl=True,labelapp=' no MAI', linewidth=3,fig=fig3,linestyle=':',nuclei=nuclei)
            #noMA.plot_nt(i[0],i[1],ne,gl=True,linestyle='-',fig=fig2,labelapp=' all off',pretitle='MAR MAD MAI, ', savename='t-dep-box_ma_te{}_t{}'.format(i[1],i[0]).replace('.','_'),nuclei=nuclei, color=['b','b'])

    ''' Ti!=Te, ni!=ne and E!=0.1 eV comparison '''
if 1==0:
    
    #for i in range(len(l)):
    #    l[i][0]=l[i][0]*3    
    #ne=1e11

    ''' E!=0.1 eV '''
    if 1==1:
        for i in l:
            gl=False
            fig1=crm.plot_nt(i[0],i[1],ne,E=0.1,gl=gl,linestyle='-',linewidth=3,labelapp=' E=0.1',Nq=True)
            crm.plot_nt(i[0],i[1],ne,E=0.5,gl=gl,linestyle='--',linewidth=3,labelapp=' E=0.5',fig=fig1,Nq=True)
            crm.plot_nt(i[0],i[1],ne,E=1,gl=gl,linestyle='-.',linewidth=3,labelapp=' E=1',fig=fig1,Nq=True)
            crm.plot_nt(i[0],i[1],ne,E=5,gl=gl,linestyle='-.',linewidth=3,fig=fig1,labelapp=' E=5',pretitle='Mol. E comparison, ', savename='t-dep-box_E_te{}_t{}'.format(i[1],i[0]).replace('.','_'),Nq=True)
        
    ''' Ti!=Te '''
    if 1==1:
        for i in l[:4]:
            gl=True
            d=0.2
            fig1=crm.plot_nt(i[0],i[1],ne,Ti=1*i[1],gl=gl,linestyle='-',linewidth=3,labelapp=' Ti=Te',Nq=True)
            crm.plot_nt(i[0],i[1],ne,gl=gl,Ti=(1-d)*i[1],linestyle='--',linewidth=3,labelapp=' Ti={}*Te'.format(1-d),fig=fig1,Nq=True)
            crm.plot_nt(i[0],i[1],ne,gl=gl,Ti=(1+d)*i[1],linestyle='-.',linewidth=3,fig=fig1,labelapp=' Ti={}*Te'.format(1+d),pretitle='Ti!=Te comparison, ', savename='t-dep-box_Ti_te{}_t{}'.format(i[1],i[0]).replace('.','_'),Nq=True)
        
    ''' ni!=ne '''
    if 1==1:
        for i in l[:4]:
            gl=False
            d=0.2
            fig1=crm.plot_nt(i[0],i[1],ne,ni=1*ne,gl=gl,linestyle='-',linewidth=3,labelapp=' ni=Tn',Nq=True)
            crm.plot_nt(i[0],i[1],ne,gl=gl,ni=(1-d)*ne,linestyle='--',linewidth=3,labelapp=' ni={}*ne'.format(1-d),fig=fig1,Nq=True)
            crm.plot_nt(i[0],i[1],ne,gl=gl,ni=(1+d)*ne,linestyle='-.',linewidth=3,fig=fig1,labelapp=' ni={}*ne'.format(1+d),pretitle='ni!=ne comparison, ', savename='t-dep-box_Ti_te{}_t{}'.format(i[1],i[0]).replace('.','_'),Nq=True)
        
        



    ''' OPTIMAL CRM '''
if 1==0:
    print('TODO')


    ''' REACTION RATES '''
if 1==0:

    # Te-dependence
    if 1==1:
        ne=1e13
        t=[]
        for i in range(6,40):#29):
            t.append(10**(-1.2+i/10))
        for i in [0,1,2,4]:
            fig1=crm.plot_uerate('output/uerates.dat',t,ne,labelapp=' CRUM',idx=i) 
            crm.plot_uerate('output/uerates_old.dat',t,ne,fig=fig1,labelapp=' UE',idx=i,color='r',linestyle='-.',savename='Te-dep_ue_v_crum_rates_idx{}_ne{:1.1E}'.format(i,ne).replace('.','_')) 



    # ne-dependence
    if 1==1:
        te=[0.5,1,2,5,10]
        n=[]
        for i in range(0,15):#29):
            n.append(10**(10+i/2))
        for t in te:
            for i in [0,1,2,4]:
                fig1=crm.plot_uerate('output/uerates.dat',t,n,labelapp=' CRUM',idx=i) 
                crm.plot_uerate('output/uerates_old.dat',t,n,fig=fig1,labelapp=' UE',idx=i,color='r',linestyle='-.',ylim=(1e-5,1e10),savename='ne-dep_ue_v_crum_rates_idx{}_Te{}'.format(i,t).replace('.','_')) 



    ''' OPTIMIZE CRM '''
if 1==0:
    for k in [7e-2,5e-2,3.45e-2,1e-3,4e-4]:
        for j in [1e11,1e12,1e13,1e14]:
            for i in [1,1.5,2,10,15]:
                print('')
                print('T={:.2E}, n={:.2E}, kappa={:.2E}'.format(i,j,k))
                print('++++++++++++++++++++++++++++++++')
                crm.crm.generate_CRM(i,j,k,epsilon=5)
    '''
    crm.crm.generate_CRM(2,1e13,7e-2)
    '''

