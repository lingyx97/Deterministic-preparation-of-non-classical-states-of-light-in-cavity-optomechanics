import qutip as qtp
import numpy as np
from timeit import default_timer as timer

from pathlib import Path
from generateOps import geOps
from parseEvoString import NumericStringParser

# https://journals-aps-org.iclibezp1.cc.ic.ac.uk/prapplied/pdf/10.1103/PhysRevApplied.3.054009




def op_dist(A,B=0):
    if B==0:
        return np.sum(np.abs(A.full()))
    if A.isket:
        A=qtp.ket2dm(A)
    if B.isket:
        B=qtp.ket2dm(B)
        
    res=np.sum(np.abs(A.full()-B.full()))
    if np.isnan(res):
        print("resulting matrix contains nan")
        return -1
    return res
    

def diss_superop(diss_ops,state_dm,t):
    # benchmark=310
    # memoryoverflow for [60,60] 2.39PB required
    if len(diss_ops)==0:
        return state_dm
    ld=0
    for op in diss_ops:
        ld+=t*qtp.lindblad_dissipator(op,op)
    return qtp.vector_to_operator(ld.expm()*qtp.operator_to_vector(state_dm))

def diss_exp(diss_ops,state_dm,t,tol=1e-20):
    # benchmark=0.12
    if len(diss_ops)==0:
        return state_dm
    res=0
    count=1
    corr=state_dm
    while (not res) or op_dist(corr,0)>tol:
        if count==20:
            print("Lindbladian probably not converging. err:{}".format(op_dist(corr,0)))
        if count>100:
            print("count:{},err:{}".format(count,op_dist(corr,0)))
            
        if not res:
            res=state_dm #initialize res when res is not defined
            
        tmp=0
        for op in diss_ops:
            n=op.dag()*op
            # input(n.diag()[:10])
            tmp+=t*(op*corr*(op.dag())-(n*corr+corr*n)/2)
        corr=tmp/count
        res+=corr
        # print(res.ptrace(1).diag())
        # input(res.tr())
        count+=1
    return res/(np.abs(res.tr()))


def hamHalfInterv(g0,ops):
    ad,a,bd,b=ops
    evo1 = 1j*np.pi*g0**2*(ad*a)**2
    evo2 = -2*g0*ad*a*(b-bd)
    evo3 = -np.pi*1j*bd*b
    return evo1.expm()*evo2.expm()*evo3.expm()


def findEvo(params):
    data,dimHS,args=params
    # eta,g0,psi=args
    t0=timer()
    nsp = NumericStringParser(dimHS, args)
    res=(nsp.eval(data)).expm()
    print("{} finished in {}s".format(params[-1],int(timer()-t0)))
    if np.isnan(np.sum(res.full())):
        print("contains nan. args={}, res={}".format(args,res))
    return res


def evo_trot(dimHS, g0,etaList,loss=[],evoOrder=5,initOccu=True,recal=0,parallel=1):
    tFin=len(etaList)
    ad, a, bd, b, idMat, initS = geOps(dimHS, initOccu)

    loss_ops=0
    if len(loss)==1:
        loss_ops=[np.sqrt(loss[0])*a]
    elif len(loss)==3:
        kap,gam,nbath=loss
        if gam==0:
            loss_ops=[np.sqrt(kap)*a]
        elif nbath:
            tmp_coef=g0*np.sqrt(gam)*np.sqrt(4/np.log(1+1/nbath)+(2*nbath+1))
            
            loss_ops_undriv=[np.sqrt(gam*(nbath+1))*(b-g0*ad*a),
                    np.sqrt(gam*nbath)*(bd-g0*ad*a),tmp_coef*ad*a]
            loss_ops_driv=[np.sqrt(gam*(nbath+1))*b,np.sqrt(gam*nbath)*bd]
        # loss_ops=[np.sqrt(kap)*a,np.sqrt(gam*(nbath+1))*b,
        #           np.sqrt(gam*nbath)*bd]
    else:
        loss_ops=[]
        
    if initS.isket:
        initS=qtp.ket2dm(initS)
    # print(str(Path(__file__).resolve()))
    # print(str(Path().resolve()))
    data = str(np.genfromtxt(str(Path().resolve().parent /
                            'common/data/evolution{}PY.txt'.format(int(evoOrder))), dtype='str', delimiter="\n"))
    
    params=[]
    psiNow = 0
    for i in range(tFin*4):
        etaPos = int(i/4)
        psiNow += 4/3*(g0*etaList[etaPos])**2*np.pi
        params.append((data,dimHS,(etaList[etaPos],g0,psiNow+i*np.pi)))

    savedFiles = Path.cwd().glob("*.qu")
    savedFilesList = [f.stem for f in savedFiles]

    fn_evoDriv="dim{},eta{},k{},loss{}EvolutionOps".format(dimHS,round(sum(etaList),5),g0,loss)
    evoDriv=0
    if (not recal) and (fn_evoDriv in savedFilesList):
        evoDriv=qtp.qload(fn_evoDriv)
        print("evolution operators loaded")
    elif parallel:
        t0=timer()
        try:
            print("generating evos")
            evoDriv=qtp.parallel_map(findEvo,params,num_cpus=parallel)
            print("evos generated. Time:{}s".format(int(timer()-t0)))
            qtp.qsave(evoDriv,fn_evoDriv)
        except:
            print("memory overflow, parallel generation aborted.")
            evoDriv=0
    evoUndriv=hamHalfInterv(g0,(ad,a,bd,b))
    # evoUndriv2=hamHalfInterv(g0,(ad,a,-bd,-b))
    
    t0=timer()
    for tPos in range(4*tFin):
        evoD=0
        if evoDriv:
            evoD=evoDriv[tPos]
        else:
            # fn_tmp="dim {},eta {},EvolutionOps{}".format(dimHS,round(sum(etaList),tFin),tPos)
            # if fn_tmp in savedFilesList:
            #     evoD=qtp.qload(fn_tmp)
            # else:
            evoD=findEvo(params[tPos])
            # qtp.qsave(evoD,fn_tmp)
        initS=evoD*initS*(evoD.dag())
        if tPos%4==0:
            print(initS.ptrace(0).diag()[:5])
        if len(loss):
            if loss_ops:
                initS=diss_exp(loss_ops,initS,2*np.pi)
            else:
                eta,g0,psi=params[tPos][-1]
                tmp_coef=g0**2*eta**2*gam/2*(4/np.log(1+1/nbath)+2*nbath+1)
                loss_ops_driv_tmp=loss_ops_driv+[np.sqrt(tmp_coef)*(a*np.exp(-1j*psi)+ad*np.exp(1j*psi))]
                initS=diss_exp(loss_ops_driv_tmp,initS,2*np.pi)
        # initS=diss_superop(loss_ops,initS,2*np.pi)
        if tPos%2==1:
            initS=evoUndriv*initS*(evoUndriv.dag())
            if len(loss):
                if loss_ops:
                    initS=diss_exp(loss_ops,initS,np.pi)
                else:
                    initS=diss_exp(loss_ops_undriv,initS,np.pi)

            
        print(initS.tr())
        curT=timer()-t0
        print("Progress: {}%, time spent: {}m{}s.".format(100*(tPos+1)/(4*tFin),int(curT/60),int(curT%60)))
    qtp.qsave(initS,"dim{},eta{},k{},loss{}FinS".format(dimHS,round(sum(etaList),5),g0,loss))
    return initS        

if __name__=="__main__":
    import os
    os.chdir(str(Path(__file__).parent.parent/"tempQU"))
    res=evo_trot([40,20],1/16,
                   [4., 4., 4., 3.563877, 3.28142102],recal=1,parallel=8,evoOrder=3)
    print(qtp.fidelity(res.ptrace(0),qtp.basis(40,2)))
    exit()
    # from lindblad_numeric import *
    # import os
    # from pathlib import Path
    # os.chdir(str(Path(__file__).parent.parent/"tempQU"))
    
    # qtp.settings.auto_tidyup = False
    # etaList=[4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291]
    # g0=0.0384615384615384
    # # dimHS = [60, 31]
    # dimHS2 = [30, 10]
    # res1=evo_trot(dimHS2,g0,etaList,recal=1,evoOrder=4,loss=[0,0,0],parallel=0)
    # # res2=evo_stepwisePsi2(dimHS, g0, etaList, len(etaList))[-1][-1][-1]
    # print(res1.ptrace(0).diag())
    # print(res1.ptrace(1).diag())
    # # print(res2.ptrace(0).diag()[:5])
    # # print(qtp.fidelity(res2.ptrace(0),qtp.fock(60,2)))
    # print(qtp.fidelity(res1.ptrace(0),qtp.fock(30,2)))
    
    # ad, a, bd, b, idMat, initS = geOps([10,10], True)
    # data = str(np.genfromtxt(str(Path(__file__).parent /
    #                         'common/data/evolution5PY.txt'), dtype='str', delimiter="\n"))
    # # input(data)
    # nsp = NumericStringParser([10,10], [4,1/20,0.15])
    # tmp=nsp.eval(data)
    # # print(tmp)
    # # print(tmp.dag()+tmp)
    # print(np.sum((tmp.dag()+tmp).full()))
    