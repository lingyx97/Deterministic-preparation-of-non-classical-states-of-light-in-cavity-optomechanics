import numpy as np
import qutip as qtp
from pathlib import Path
from bisect import bisect
import os
from generateOps import geOps
from smoothStep import generateSmoStepData as gSSD
import numbers

qtp.settings.auto_tidyup = False

###########################################################################################


def HI_coef(ttt, args):
        tfin,tIntervEachEta,eta,psiList,wc=args.values()
        if ttt>tfin:
            return 0
        tPos = int(ttt/np.pi)
        if tPos % 5 == 4:
            tType = 0
            tPos = 0
        else:
            tType = 1
            tPos = 2*int(tPos/5)+int((tPos % 5)/2)
        etaPos = min(int(tPos/(4*tIntervEachEta)), len(eta)-1)
        # print(str(ttt/(2*np.pi))+"___"+str(tPos),end="\t\t\t\r")
        return eta[etaPos]*tType*np.exp(1j*wc*ttt)*np.exp(-1j*psiList[tPos])*np.cos(2*ttt)

def HI_coef_dag(ttt, args):
            return np.conjugate(HI_coef(ttt, args))

def hamiltonian(params,ops):
        ad,a,bd,b=ops
        wc,k=params
        H0 = bd*b
        H1 = (-k*ad*a*(bd+b)+wc*ad*a)
        Hdriv = -2j*a

        HI = [Hdriv, HI_coef]
        HI2 = [Hdriv.dag(), HI_coef_dag]

        return [H0+H1, HI, HI2]


def evo_stepwisePsi1(dimHS, k, eta, t, wc=2, step_trans_params=[0, 0], thermalOccu=0, loss=0,show_diags=True,recal=0):
    # only one mesolve for all driv intervals
    step_trans_dur, step_trans_reso = step_trans_params
    
    etaLen=0
    try:
        eta[0]
        etaLen = len(eta)
    except TypeError:
        eta = [eta]
        etaLen = 1

    tIntervEachEta = int(t/etaLen)
    if(tIntervEachEta != t/etaLen):
        print("tIntervEachEta!=t/etaLen")

    ad, a, bd, b, idMat, initS1 = geOps(dimHS, True)
    if not thermalOccu:
        initS = initS1
    else:
        initS = qtp.tensor(qtp.fock_dm(
            dimHS[0], 0), qtp.thermal_dm(dimHS[1], thermalOccu))

    step_trans_dur *= (2*np.pi)


    psiCorr = [0]
    if etaLen == 1:
        psiCorr = [4/3*(k*eta[0])**2*np.pi*i
                   for i in range(etaLen*4*tIntervEachEta)]
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
    else:
        psiNow = 0
        for i in range(etaLen*4*tIntervEachEta):
            etaPos = int(i/(4*tIntervEachEta))
            psiNow += 4/3*(k*eta[etaPos])**2*np.pi
            psiCorr.append(psiNow)
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
        print("smoothing the step func", end='\r')
        if step_trans_params != [0, 0]:
            psiRefPos = np.arange(-step_trans_dur/2, step_trans_dur /
                                  2+step_trans_reso+1e-6, step_trans_reso)
            transPsis = [gSSD(psiRefPos, 0, step_trans_dur, [psiList[i], psiList[i+1]], label=i/(len(psiList)-1))
                         for i in range(len(psiList)-1)]
        else:
            transPsis = [0, 0]
        print("smoothing completed", end="\r")
        try:
            qtp.qsave(transPsis, 'temp')
        except MemoryError:
            print("saving failed (memory full)")

    def hamiltonian():

        def psiSmooth(tPos, t):
            tRenorm = t % (2*np.pi)
            if tRenorm < step_trans_dur/2 and tPos != 0:
                pos = bisect(psiRefPos, tRenorm)
                return transPsis[tPos-1][pos]
            if tRenorm > 2*np.pi-step_trans_dur/2 and tPos != len(psiList)-1:
                pos = bisect(psiRefPos, tRenorm-2*np.pi)
                return transPsis[tPos][pos]
            return psiList[tPos]

        def HI_coef(ttt, args):
            if ttt>tfin:
                return 0
            tPos = int(ttt/np.pi)
            if tPos % 5 == 4:
                tType = 0
                tPos = 0
            else:
                tType = 1
                tPos = 2*int(tPos/5)+int((tPos % 5)/2)
            etaPos = min(int(tPos/(4*tIntervEachEta)), len(eta)-1)
            # print(str(ttt/(2*np.pi))+"___"+str(tPos),end="\t\t\t\r")

            if step_trans_dur == 0:
                return eta[etaPos]*tType*np.exp(1j*wc*ttt)*np.exp(-1j*psiList[tPos])*np.cos(2*ttt)
                
            return eta[etaPos]*tType*np.exp(1j*wc*ttt)*np.exp(-1j*psiSmooth(tPos, ttt))*np.cos(2*ttt)

        def HI_coef_dag(ttt, args):
            return np.conjugate(HI_coef(ttt, args))

        H0 = bd*b
        H1 = (-k*ad*a*(bd+b)+wc*ad*a)
        Hdriv = -2j*a

        HI = [Hdriv, HI_coef]
        HI2 = [Hdriv.dag(), HI_coef_dag]

        return [H0+H1, HI, HI2]

    trange = np.arange(0, t*5*2*np.pi,np.pi/10)
    tfin = trange[-1]
    tol = 1e-10

    def expects(tnow, state):
        diags0=0
        diags1=0
        if show_diags:
            diags0 = state.ptrace(0).diag()
            diags1 = state.ptrace(1).diag()
        if tnow == tfin:
            return [diags0, diags1, state]
        return [diags0, diags1]

    result = 0

    savedFiles = Path.cwd().glob("*.qu")
    savedFilesList = [f.stem for f in savedFiles]
    fn = str((dimHS, k, sum(eta), t))
    if thermalOccu:
        fn = str((dimHS, k, sum(eta), t, thermalOccu))
    if loss:
        fn = str((dimHS, k, sum(eta), t, loss))
    if (fn in savedFilesList) and (not recal):
        try:
            res=qtp.qload(fn)
            print(fn+" data loaded")
            return res
        except FileNotFoundError:
            print("load failed")
            print(fn)
            print(fn in savedFilesList)

    
    if loss:
        if len(loss)!=3:
            print("wrong loss param:")
            print(loss)
            exit()
        kap, gam, nbath=loss
        tmp_coef=0
        if nbath:
            tmp_coef=np.sqrt((4*gam*k**2)/np.log(1+1/nbath))
        loss_ops=[np.sqrt(kap)*a,np.sqrt(gam*(nbath+1))*(b-k*ad*a),
                  np.sqrt(gam*nbath)*(bd-k*ad*a), tmp_coef*ad*a]
        result = qtp.mesolve(hamiltonian(), initS, trange, loss_ops, expects,
                            options=qtp.Options(atol=tol,
                                                rtol=tol,
                                                nsteps=2**31-1,
                                                norm_tol=tol))#, progress_bar=True)
    else:
        try:
            result = qtp.mesolve(hamiltonian(), initS, trange, [], expects,
                                options=qtp.Options(atol=tol,
                                                    rtol=tol,
                                                    nsteps=2**31-1,
                                                    norm_tol=tol))#, progress_bar=True)
        except TypeError:
            print(k)
            print(eta)
            print(dimHS)
            print(tfin)
            exit()
        
    result = [[t/(10*np.pi) for t in trange], result.expect]
    try:
        qtp.qsave(result, fn)
    except FileNotFoundError:
        print(fn)
        print(os.getcwd())
        print("save failed")

    return result


def evo_stepwisePsi2(dimHS, k, eta, t, wc=2, step_trans_params=[0, 0], 
                     thermalOccu=0, loss=0,show_diags=True,recal=0,progress_bar=None,save=1):
    # mesolve for each interval
    step_trans_dur, step_trans_reso = step_trans_params
    
    etaLen=0
    try:
        eta[0]
        etaLen = len(eta)
    except TypeError:
        eta = [eta]
        etaLen = 1

    tIntervEachEta = int(t/etaLen)
    if(tIntervEachEta != t/etaLen):
        print("tIntervEachEta!=t/etaLen")

    ad, a, bd, b, idMat, initS1 = geOps(dimHS, True)
    if not thermalOccu:
        initS = initS1
    else:
        initS = qtp.tensor(qtp.fock_dm(
            dimHS[0], 0), qtp.thermal_dm(dimHS[1], thermalOccu))

    step_trans_dur *= (2*np.pi)


    psiCorr = [0]
    if etaLen == 1:
        psiCorr = [4/3*(k*eta[0])**2*np.pi*i
                   for i in range(etaLen*4*tIntervEachEta)]
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
    else:
        psiNow = 0
        for i in range(etaLen*4*tIntervEachEta):
            etaPos = int(i/(4*tIntervEachEta))
            psiNow += 4/3*(k*eta[etaPos])**2*np.pi
            psiCorr.append(psiNow)
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
        print("smoothing the step func", end='\r')
        if step_trans_params != [0, 0]:
            psiRefPos = np.arange(-step_trans_dur/2, step_trans_dur /
                                  2+step_trans_reso+1e-6, step_trans_reso)
            transPsis = [gSSD(psiRefPos, 0, step_trans_dur, [psiList[i], psiList[i+1]], label=i/(len(psiList)-1))
                         for i in range(len(psiList)-1)]
        else:
            transPsis = [0, 0]
        print("smoothing completed", end="\r")

    # def hamiltonian():

    #     def HI_coef(ttt, args):
    #         if ttt>tfin:
    #             return 0
    #         tPos = int(ttt/np.pi)
    #         if tPos % 5 == 4:
    #             tType = 0
    #             tPos = 0
    #         else:
    #             tType = 1
    #             tPos = 2*int(tPos/5)+int((tPos % 5)/2)
    #         etaPos = min(int(tPos/(4*tIntervEachEta)), len(eta)-1)
    #         # print(str(ttt/(2*np.pi))+"___"+str(tPos),end="\t\t\t\r")
    #         return eta[etaPos]*tType*np.exp(1j*wc*ttt)*np.exp(-1j*psiList[tPos])*np.cos(2*ttt)
                

    #     def HI_coef_dag(ttt, args):
    #         return np.conjugate(HI_coef(ttt, args))

    #     H0 = bd*b
    #     H1 = (-k*ad*a*(bd+b)+wc*ad*a)
    #     Hdriv = -2j*a

    #     HI = [Hdriv, HI_coef]
    #     HI2 = [Hdriv.dag(), HI_coef_dag]

    #     return [H0+H1, HI, HI2]

    tol = 1e-10


    result = 0

    savedFiles = Path.cwd().glob("*.qu")
    savedFilesList = [f.stem for f in savedFiles]
    fn = str((dimHS, round(k,5), round(sum(eta),5), t))
    
    tPosLast=-1

    if thermalOccu:
        fn = str((dimHS, round(k,5), round(sum(eta),5), t, thermalOccu))
    if loss:
        fn = str((dimHS, round(k,5), round(sum(eta),5), t, loss))
    if (fn in savedFilesList) and (not recal):
        try:
            res=qtp.qload(fn)
            print(fn+" data loaded")
            if isinstance(res[0], numbers.Number):
                tPosLast, resLast=res
                print("{} not completed, to be resumed".format(fn))
            else:
                return res
        except FileNotFoundError:
            print("load failed")
            print(fn)
            print(fn in savedFilesList)
    else:
        print("{} file not ready, to be calculated".format(fn))

    
    if loss:
        if len(loss)!=3:
            print("wrong loss param:")
            print(loss)
            exit()
        kap, gam, nbath=loss
        tmp_coef=0
        if nbath:
            tmp_coef=np.sqrt((4*gam*k**2)/np.log(1+1/nbath))
        loss_ops=[np.sqrt(kap)*a,np.sqrt(gam*(nbath+1))*(b-k*ad*a),
                  np.sqrt(gam*nbath)*(bd-k*ad*a), tmp_coef*ad*a]
        result=[]
        trange=[]
        for tPos in range(1,int(t*5*2+1)):
            if tPos<tPosLast:
                continue
            if tPos==tPosLast:
                trange,result=resLast
                initS=result[-1][-1]
            tstmp= np.arange((tPos-1)*np.pi, tPos*np.pi,np.pi/10)
            
            tfin=tstmp[-1]

            def expects(tnow, state):
                diags0=0
                diags1=0
                if show_diags:
                    diags0 = state.ptrace(0).diag()
                    diags1 = state.ptrace(1).diag()
                if tnow == tfin:
                    return [diags0, diags1, state]
                return [diags0, diags1]

            
            
            args=[tfin,tIntervEachEta,eta,psiList,wc]
            keys=["tfin","tIntervEachEta","eta","psiList","wc"]
            args=dict(zip(keys,args))
            ham=hamiltonian((wc,k),(ad,a,bd,b))
            restmp = qtp.mesolve(ham, initS, tstmp, loss_ops, expects,args=args,
                                options=qtp.Options(atol=tol,
                                                    rtol=tol,
                                                    nsteps=2**31-1,
                                                    norm_tol=tol), progress_bar=progress_bar).expect#, progress_bar=True)
            initS=restmp[-1][-1]
            result=result+list(restmp)
            trange=trange+list(tstmp)
            if save:
                qtp.qsave([tPos,[trange,result]],fn)
            print(tPos/int(t*5*2),end="\r")
            
    else:
        try:
            result=[]
            trange=[]
            for tPos in range(1,int(t*5*2+1)):
                if tPos<tPosLast:
                    continue
                if tPos==tPosLast:
                    trange,result=resLast
                    initS=result[-1][-1]
                tstmp= np.arange((tPos-1)*np.pi, tPos*np.pi+1e-4,np.pi/10)
                
                tfin=tstmp[-1]

                def expects(tnow, state):
                    diags0=0
                    diags1=0
                    if show_diags:
                        diags0 = state.ptrace(0).diag()
                        diags1 = state.ptrace(1).diag()
                    # print("{}:{}".format(tnow,tfin))
                    if tnow == tfin:
                        return [diags0, diags1, state]
                    return [diags0, diags1]

                
                args=[tfin,tIntervEachEta,eta,psiList,wc]
                keys=["tfin","tIntervEachEta","eta","psiList","wc"]
                args=dict(zip(keys,args))
                ham=hamiltonian((wc,k),(ad,a,bd,b))
                restmp = qtp.mesolve(ham, initS, tstmp, [], expects,args=args,
                                    options=qtp.Options(atol=tol,
                                                        rtol=tol,
                                                        nsteps=2**31-1,
                                                        norm_tol=tol), progress_bar=progress_bar).expect#, progress_bar=True)
                initS=restmp[-1][-1]
                result=result+list(restmp)
                trange=trange+list(tstmp)
                if save:
                    qtp.qsave([tPos,[trange,result]],fn)
                print(tPos/int(t*5*2),end="\r")
                
        except TypeError:
        # except FileExistsError:
            print("error:")
            print(k)
            print(eta)
            print(dimHS)
            print(tfin)
            exit()
        
    result = [[t/(10*np.pi) for t in trange], result]
    if save:
        try:
            qtp.qsave(result, fn)
        except FileNotFoundError:
            print(fn)
            print(os.getcwd())
            print("save failed")

    return result


def evo_stepwisePsi_MC(dimHS, k, eta, t, wc=2, step_trans_params=[0, 0], 
                     thermalOccu=0, loss=0,show_diags=True,recal=0,progress_bar=None,save=1):
    # mesolve for each interval
    step_trans_dur, step_trans_reso = step_trans_params
    
    etaLen=0
    try:
        eta[0]
        etaLen = len(eta)
    except TypeError:
        eta = [eta]
        etaLen = 1

    tIntervEachEta = int(t/etaLen)
    if(tIntervEachEta != t/etaLen):
        print("tIntervEachEta!=t/etaLen")

    ad, a, bd, b, idMat, initS = geOps(dimHS, True)

    step_trans_dur *= (2*np.pi)


    psiCorr = [0]
    if etaLen == 1:
        psiCorr = [4/3*(k*eta[0])**2*np.pi*i
                   for i in range(etaLen*4*tIntervEachEta)]
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
    else:
        psiNow = 0
        for i in range(etaLen*4*tIntervEachEta):
            etaPos = int(i/(4*tIntervEachEta))
            psiNow += 4/3*(k*eta[etaPos])**2*np.pi
            psiCorr.append(psiNow)
        psiList = [(psiCorr[i]+i*np.pi) for i in range(len(psiCorr))]
        print("smoothing the step func", end='\r')
        if step_trans_params != [0, 0]:
            psiRefPos = np.arange(-step_trans_dur/2, step_trans_dur /
                                  2+step_trans_reso+1e-6, step_trans_reso)
            transPsis = [gSSD(psiRefPos, 0, step_trans_dur, [psiList[i], psiList[i+1]], label=i/(len(psiList)-1))
                         for i in range(len(psiList)-1)]
        else:
            transPsis = [0, 0]
        print("smoothing completed", end="\r")


    tol = 1e-10


    result = 0

    savedFiles = Path.cwd().glob("*.qu")
    savedFilesList = [f.stem for f in savedFiles]
    fn = "{},{},{}MC".format(dimHS, round(k,5), round(sum(eta),t))
    
    tPosLast=-1

    if loss:
        fn = "{},{},{},{}MC".format(dimHS, round(k,5), round(sum(eta),t),loss)
    if (fn in savedFilesList) and (not recal):
        try:
            res=qtp.qload(fn)
            print(fn+" data loaded")
            if isinstance(res[0], numbers.Number):
                tPosLast, resLast=res
                print("{} not completed, to be resumed".format(fn))
            else:
                return res
        except FileNotFoundError:
            print("load failed")
            print(fn)
            print(fn in savedFilesList)
    else:
        print("{} file not ready, to be calculated".format(fn))

    if loss:
        if len(loss)!=3:
            print("wrong loss param:")
            print(loss)
            exit()
        kap, gam, nbath=loss
        tmp_coef=0
        if nbath:
            tmp_coef=np.sqrt((4*gam*k**2)/np.log(1+1/nbath))
        loss_ops=[np.sqrt(kap)*a,np.sqrt(gam*(nbath+1))*(b-k*ad*a),
                  np.sqrt(gam*nbath)*(bd-k*ad*a), tmp_coef*ad*a]
        result=[]
        trange=[]
        for tPos in range(1,int(t*5*2+1)):
            if tPos<tPosLast:
                continue
            if tPos==tPosLast:
                trange,result=resLast
                initS=result[-1][-1]
            tstmp= np.arange((tPos-1)*np.pi, tPos*np.pi,np.pi/10)
            
            tfin=tstmp[-1]

            
            args=[tfin,tIntervEachEta,eta,psiList,wc]
            keys=["tfin","tIntervEachEta","eta","psiList","wc"]
            args=dict(zip(keys,args))
            ham=hamiltonian((wc,k),(ad,a,bd,b))
            restmp = qtp.mcsolve(ham, initS, tstmp, loss_ops, [],
                                options=qtp.Options(atol=tol,
                                                    rtol=tol,
                                                    nsteps=2**31-1,
                                                    norm_tol=tol), progress_bar=progress_bar,args=args).expect#, progress_bar=True)
            
            initS=np.average(tstmp,axis=0)[-1]
            tstmp=tstmp[-1]
            result.append(initS)
            trange.append(tstmp)
            if save:
                qtp.qsave([tPos,[trange,result]],fn)
            print(tPos/int(t*5*2),end="\r")
            
    else:
        print("no loss ops. Please use mesolve")
        return -1
    
    result = [[t/(10*np.pi) for t in trange], result]
    if save:
        try:
            qtp.qsave(result, fn)
        except FileNotFoundError:
            print(fn)
            print(os.getcwd())
            print("save failed")

    return result


if __name__ == "__main__":
    from timeit import default_timer as timer
    # t0=timer()
    res=evo_stepwisePsi2([60,30], 1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5,save=1)[1][-1][-1]
    
    print(qtp.fidelity(res.ptrace(0),qtp.basis(60,2)))
    # print(timer()-t0)
    # t0=timer()
    # res2=evo_stepwisePsi_MC([30,10], 1/26,
    #                      [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
    #             2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291], 
    #             16,loss=[1e-4,0,0],save=0)[1][-1][-1]
    # print(timer()-t0)
    # print(qtp.fidelity(res,res2))
    # print(res.ptrace(0).diag()[:10])
    # print(res2.ptrace(0).diag()[:10])