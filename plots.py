from plot_helper import *

from matplotlib.legend_handler import HandlerTuple
# import matplotlib.ticker as ticker
import pandas as pd
from ast import literal_eval

qtp.settings.auto_tidyup = False
fontsize=14
font = {"family" : "serif",
        "serif" : ["Computer Modern Serif"],
        'size'   : fontsize}

plt.rc('xtick',labelsize=fontsize)
plt.rc('ytick',labelsize=fontsize)
plt.rcParams['axes.axisbelow'] = True

print("start:")
# PLOT=int(input())
PLOT=8


if PLOT==2:
    chain=[4.      ,   4.    ,     4.      ,   4.       ,  4.   ,      3.56,
            2.19 ,2.74]
    pointsToPlot = 2*len(chain)
    dimHS=[60,60]
    ax2,last_state=plotEvo(chain, pointsToPlot,dimHS, 1/26,[0,7e-4])
    print(qtp.fidelity(last_state.ptrace(0),qtp.fock(dimHS[0],2))**2)

    # os.chdir(str(Path(__file__).parent/"tempQU"))
    # savefig("1_k26t16FockSimul.pdf")

    plt.show()
    

    
if PLOT==3:
    colors=generate_cmap(3)
    ts,data_k26_E4=np.transpose(qtp.qload("plot31"))
    data_k13_E4=np.transpose(qtp.qload("plot32"))[1]
    data_kFlex_E4=np.transpose(qtp.qload("plot33"))[1]
    plt.plot(ts,data_kFlex_E4,'-o',label=r'(a) optimal',color=colors[0])
    plt.plot(ts,data_k26_E4,'--o',label=r'(b) $g_0=\omega_\mathrm{m}/26$',color=colors[1])
    plt.plot(ts,data_k13_E4,':o',label=r'(c) $g_0=\omega_\mathrm{m}/13$',color=colors[2])
    # plot (df_k13_E3['t'],df_k13_E3['fid'],':o',label='(d)')
    plt.xlabel(r"evolution time ($\cal T$)",fontsize=fontsize-2)
    plt.ylabel("fidelity",fontsize=fontsize-2)
    plt.legend()
    # plt.show()
        
    plt.savefig("1_lowFid_fid-t.pdf", bbox_inches='tight')
            
if PLOT==4:
    ks = [26, 25, 24.5, 23.5, 22.5, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13]
    ks=[1/k for k in ks]
    
    ts=[16-i for i in range(15)]
    plt.scatter(ts,ks)
    # plot (df_k13_E3['t'],df_k13_E3['fid'],':o',label='(d)')
    # plt.plot(ts,[ts[j]*ks[j]**2 for j in range(len(ts))])
    # exit()
    plt.xlabel(r"evolution time ($\cal T$)",fontsize=fontsize-2)
    plt.ylabel(r"coupling strength $g_0/\omega_\mathrm{m}$",fontsize=fontsize-2)
    plt.ylim([0.035,0.08])
    plt.show()
    # plt.savefig("1_lowFid_k-t.pdf", bbox_inches='tight')


if PLOT==5:
    k=1/26
    chain=[4., 3.69291531, 3.51161495, 3.686096 , 2.60729845, 3.19989213,
           3.13977697,1.51847303,2.00439795,3.9995126]
    finPhase=(sum([2*x**2 for x in chain])*4*np.pi*k**2/3)
    pointsToPlot = len(chain)
    dimHS=[60,30]
    ax2,last_state=plotEvo(chain, pointsToPlot,dimHS, 1/26,[0,1.8*1e-4])

    # maxangle=0
    # maxfid=0
    # for angle in np.arange(0,2*np.pi+0.01,0.01):
    #     fid=qtp.fidelity(last_state.ptrace(0),(qtp.fock(dimHS[0],2)+np.exp(1j*angle)*qtp.fock(dimHS[0],0)).unit())**2
    #     if maxfid<fid:
    #         maxangle=angle
    #         maxfid=fid
    # print(maxangle)
    # print(maxfid)
    # print(qtp.fidelity(last_state.ptrace(0),np.exp(-1j*1.86)*qtp.fock(dimHS[0],2)+qtp.fock(dimHS[0],0)).unit())
    print(qtp.fidelity(last_state.ptrace(0),(-np.exp(1j*finPhase)*qtp.fock(dimHS[0],2)+qtp.fock(dimHS[0],0)).unit()))

    # plt.savefig("1_k26t16SuperPosSimul.pdf")

    # plt.show()
    

if PLOT==7:
    df=pd.read_csv('optimizedDrivingData.csv',converters={"eList":literal_eval },index_col=0)
    # input(fidnum(df.values.tolist()[100]))
    # res=qtp.parallel_map(fidnum,df.values.tolist())
    # input(res)
    # df['fidnum']=res
    # df.to_csv('optimizedDriv.csv')

    # eMaxList=[2,3,4,5,6,8,10]
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,1]},figsize=(8.0, 6.0))
    ax1, ax2 = axs
    eMaxList=[2,3,4,5,6]
    neMax=len(eMaxList)
    colors=generate_cmap(neMax)
    
    for j in range(neMax):
        g013=df[np.logical_and(df['eMax']==eMaxList[j], df['g0']>1/26+0.01)]
        g026=df[np.logical_and(df['eMax']==eMaxList[j], df['g0']<1/13-0.01)]
        g0s=[g026,g013]
        for i in range(2):
            axs[i].plot(g0s[i]['eLen'],g0s[i]['fid'],color=colors[j],linestyle='solid',label=r"$E_\mathrm{max}=$"+str(eMaxList[j]), marker='.')
            axs[i].plot(g0s[i]['eLen'],g0s[i]['fidnum'],color=colors[j],linestyle='dashed', marker='.')

    title=[r"(a) $g_0=\omega_\mathrm{m}/26$", r"(b) $g_0=\omega_\mathrm{m}/13$"]
    for i in range(2):
        axs[i].set_ylim([0.3,1.05])
        axs[i].set_title(title[i],fontsize=fontsize)
        axs[i].title.set_position([.05,.75])
        axs[i].set_ylabel(r"fidelity",fontsize=fontsize)
    axs[0].legend(fontsize=fontsize-2,loc="upper left")
    # plt.show()
    plt.tight_layout()
    axs[1].set_xlabel(r"evolution time $(\cal T)$",fontsize=fontsize)
    
    # plt.show()    
    plt.savefig("1_fidVsT.pdf", bbox_inches='tight')

if PLOT==8:
    os.chdir(Path(__file__).parent/"plot_data"/"5")
    # print(os.getcwd())
    # exit()
    # alphas=["(\u2170)",
    #         "(\u2171)",
    #         "(\u2172)"]
    labels=[r'$t_\mathrm{f}=16\mathcal{T}$',r'$t_\mathrm{f}=10\mathcal{T}$',r'$t_\mathrm{f}=5\mathcal{T}$']
    dimHS=[60,60]
    dataFN=["{}plot6{}".format(dimHS,i) for i in range(3)]
    data=[qtp.qload(fn) for fn in dataFN]
    data=[list(map(list, zip(*d))) for d in data]
    data=[d[-3:] for d in data]
    # input([[data[2][1][10*i],data[2][2][10*i]] for i in [1,2,5]])
    # input(data[1])
    # print([d[-2:] for d in data])
    # exit()
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,1]})
    colors=generate_cmap(len(data)+1)
    for i in range(2):
        for j in range(3):
            axs[i].plot(data[j][0],data[j][i+1],label=labels[j],color=colors[j])
        axs[i].plot([0,data[0][0][-1]],[1.,1.],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.95,0.95],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.9,0.9],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.85,0.85],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.8,0.8],linestyle='dotted',color='grey')
        axs[i].plot([1,1],[0.45,1],linestyle='dotted',color='grey')
        axs[i].title.set_position([.1,1.05])
        axs[i].set_ylim([0.6,1.01])
        axs[i].set_xlim([0,5])
        axs[i].set_yticks([0.6,0.8,0.9,1.0])
        axs[i].set_yticklabels([0.6,0.8,0.9,1.0])

    axs[1].legend(loc=3,fontsize=fontsize-2)
    
    axs[0].set_title(r"(a) Fidelity with $|2\rangle$ ($F_\mathrm{l}$)",fontsize=fontsize,loc='left')
    axs[1].set_title(r"(b) Fidelity with $\rho_\mathrm{th}(\bar{n}_\mathrm{th}=0)$ ($F_\mathrm{i}$)",fontsize=fontsize,loc='left')
    # axs[0].set_ylabel("fidelity",fontsize=fontsize)
    # axs[1].set_ylabel("fidelity",fontsize=fontsize)
    # axs[0].title.set_position([.15,.75])
    # axs[1].title.set_position([.2,.75])
    fig.tight_layout(pad=1.5)
    
    # ax[0].legend(loc=6)
    # ax[1].legend(loc=6)
    
    # axs[1].set_xlabel(r"mean thermal phonon numer $\bar{n}_\mathrm{th}$",fontsize=fontsize-2)
    

    plt.xlabel(r"mean thermal phonon numer ($\bar{n}_\mathrm{th}$)",fontsize=fontsize-2)

    plt.savefig("1_thermalInit.pdf")
# 
    # plt.show()

if PLOT==9:
    os.chdir(Path(__file__).parent/"plot_data"/"5")
    # print(os.getcwd())
    # exit()
    # alphas=["(\u2170)",
    #         "(\u2171)",
    #         "(\u2172)"]
    labels=[r'$t_\mathrm{f}=16\mathcal{T}$',r'$t_\mathrm{f}=10\mathcal{T}$',r'$t_\mathrm{f}=5\mathcal{T}$']
    
    dimHS=[60,30]
    dataFN=["{}plot7{}".format(dimHS,i) for i in range(3)]
    data=[qtp.qload(fn) for fn in dataFN]
    data=[list(map(list, zip(*d))) for d in data]
    data=[d[-3:] for d in data]
    for d in data:
        d[0]=[dd[0] for dd in d[0]]
    # data[0]=[d[0] for d in data[0]]
    # print([d[-2:] for d in data])
    # exit()
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,1]})
    colors=generate_cmap(len(data)+1)

    for i in range(2):
        for j in range(3):
            axs[i].plot(data[j][0],data[j][i+1],label=labels[j],color=colors[j])
        axs[i].plot([0,data[0][0][-1]],[1.,1.],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.95,0.95],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.9,0.9],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.85,0.85],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.8,0.8],linestyle='dotted',color='grey')
        axs[i].plot([1e-4,1e-4],[0.6,1],linestyle='dotted',color='grey')
        axs[i].plot([1e-3,1e-3],[0.6,1],linestyle='dotted',color='grey')
        axs[i].title.set_position([.1,1.05])
        axs[i].set_ylim([0.6,1.01])
        axs[i].set_xlim([1e-5,2e-3])
        axs[i].set_yticks([0.6,0.8,0.9,1.0])
        axs[i].set_yticklabels([0.6,0.8,0.9,1.0])
        axs[i].set_xscale("log")
        
        axs[i].set_xticks([1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3])
        axs[i].set_xticklabels([r'$0.1$',r'$0.2$',r'$0.5$',
                           r'$1$',r'$2$',r'$5$',
                           r'$10$',r'$20$'])
        # axs[i].set_xticks([j*1e-5 for j in range(1,10)]+[j*1e-4 for j in range(1,10)]+[1e-3,2e-3])
        # axs[i].set_xticklabels([r'$10^{-5}$',r'$2\times10^{-5}$',None,None,r'$5\times10^{-5}$',None,None,None,None,
        #                         r'$10^{-4}$',r'$2\times10^{-4}$',None,None,r'$5\times10^{-4}$',None,None,None,None,
        #                         r'$10^{-3}$',r'$2\times10^{-3}$'])
        
        # axs[i].get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        # axs[i].tick_params(axis='x', which='minor')
        # axs[i].xaxis.set_minor_formatter(ticker.FormatStrFormatter("%.1f"))
    
    axs[1].legend(loc=3,fontsize=fontsize)
    
    axs[0].set_title(r"(a) Fidelity with $|2\rangle$ ($F_\mathrm{l}$)",fontsize=fontsize,loc='left')
    axs[1].set_title(r"(b) Fidelity with close-system final state ($F_\mathrm{i}$)",fontsize=fontsize,loc='left')
    
    
    # ax[0].legend(loc=6)
    # ax[1].legend(loc=6)
    
    
    plt.xlabel(r"optical decay rate ($10^4\kappa/\omega_\mathrm{m}$)",fontsize=fontsize-2)
    fig.tight_layout(pad=1.5)

    plt.savefig("1_opticaldecay.pdf")

    # plt.show()

if PLOT==10:
    os.chdir(Path(__file__).parent/"plot_data"/"4")
    # print(os.getcwd())
    # exit()
    # alphas=["(\u2170)",
    #         "(\u2171)",
    #         "(\u2172)"]
    labels=[r'$t_\mathrm{f}=16\mathcal{T}$',r'$t_\mathrm{f}=10\mathcal{T}$',r'$t_\mathrm{f}=5\mathcal{T}$']
    
    # dimHS=[60,60]
    dimHS=[60,25]
    dataFN=["{}plot8{}".format(dimHS,i) for i in range(9)]
    data=[qtp.qload(fn) for fn in dataFN]
    data=[list(map(list, zip(*d))) for d in data]
    data=[d[-3:] for d in data]
    for j in range(len(data)):
        data[j][0]=[dd[1] for dd in data[j][0]]
        if j >=6:
            # encounter serious truncation error for the last few data sets.
            data[j][0]=data[j][0][:15]
            data[j][1]=data[j][1][:15]
            data[j][2]=data[j][2][:15]
    # data[0]=[d[0] for d in data[0]]
    # print([d[-2:] for d in data])
    # exit()
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,1]})
    colors=generate_cmap(4)
    plots=[]
    linestyles=['solid','dashed','dotted']
    for i in range(2):
        for j in range(3):
            p1=axs[i].plot(data[3*j][0],data[3*j][i+1],color=colors[0],linestyle=linestyles[j])
            p2=axs[i].plot(data[3*j+1][0],data[3*j+1][i+1],color=colors[1],linestyle=linestyles[j])
            p3=axs[i].plot(data[3*j+2][0],data[3*j+2][i+1],color=colors[2],linestyle=linestyles[j])
            plots+=[p1,p2,p3]
        axs[i].plot([0,data[0][0][-1]],[1.,1.],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.95,0.95],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.9,0.9],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.85,0.85],linestyle='dotted',color='grey')
        axs[i].plot([0,data[0][0][-1]],[0.8,0.8],linestyle='dotted',color='grey')
        axs[i].plot([1e-4,1e-4],[0.6,1],linestyle='dotted',color='grey')
        axs[i].plot([1e-3,1e-3],[0.6,1],linestyle='dotted',color='grey')
        axs[i].title.set_position([.1,1.05])
        axs[i].set_ylim([0.6,1.01])
        axs[i].set_xlim([1e-5,2e-3])
        axs[i].set_yticks([0.6,0.8,0.9,1.0])
        axs[i].set_yticklabels([0.6,0.8,0.9,1.0])
        axs[i].set_xscale("log")
        
        axs[i].set_xticks([1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3])
        axs[i].set_xticklabels([r'$0.1$',r'$0.2$',r'$0.5$',
                           r'$1$',r'$2$',r'$5$',
                           r'$10$',r'$20$'])
    
    tmp=[tuple([plots[3*i+j][0] for i in range(3)]) for j in range(3)]
    # input(tmp)
    l = axs[1].legend(tmp, labels,
                      handler_map={tuple: HandlerTuple(ndivide=None)},loc=3,fontsize=fontsize)
    
    # axs[1].legend(loc=3,fontsize=fontsize)
    
    axs[0].set_title(r"(a) Fidelity with $|2\rangle$ ($F_\mathrm{l}$)",fontsize=fontsize,loc='left')
    axs[1].set_title(r"(b) Fidelity with close-system final state ($F_\mathrm{i}$)",fontsize=fontsize,loc='left')
    
    
    # ax[0].legend(loc=6)
    # ax[1].legend(loc=6)
    
    
        # savefig("6_thermalInit.pdf")

    plt.xlabel(r"mechanical decay rate ($10^4\gamma/\omega_\mathrm{m}$)",fontsize=fontsize-2)
    fig.tight_layout(pad=1.5)
    
    plt.show()

    # plt.savefig("1_mechanicaldecay.pdf")
# 