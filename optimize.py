from headers import *
import scipy as scp
import pandas as pd
from ast import literal_eval
import random

# from timeit import default_timer as timer
# from ast import literal_eval

qtp.settings.auto_tidyup = False

def fidnum_wrapper(params):
    etaList,tTotal, kk=params
    dimHS=[120,60]
    return fidnum(etaList, tTotal, dimHS, kk)

def fid2_wrapper(params):
    etaList,tTotal, kk=params
    print("{},fid1:{}".format(params,1-evoFullInfid(30,kk,etaList,fins=0)))
    return 1-evoFullInfid([30,10],kk,etaList,fins=0)

def fidnum(etaList, tTotal, dimHS=[60,30], kk=1/50,fins=0):
    results = evo_stepwisePsi2(dimHS, kk, etaList, tTotal, wc=0)[-1]
          
    last_state = results[-1][-1].ptrace(0)
    print("diagnum={}".format(last_state.diag()[:30]))
    
    if fins:
        psirange=np.arange(0,2*np.pi,np.pi/1000)
        maxfid=0
        maxpsi=0
        for psi in psirange:
            curfid=qtp.fidelity(last_state,(np.exp(1j*psi)*fins[1]*qtp.basis(dimHS[0],2)+fins[0]*qtp.basis(dimHS[0],0)).unit())
            if curfid>maxfid:
                maxfid=curfid
                maxpsi=psi
        print("rot fit: psi={}, maxfid={}".format(maxpsi,maxfid))
        return qtp.fidelity(last_state,(fins[1]*np.exp(1j*fins[2])*qtp.basis(dimHS[0],2)+fins[0]*qtp.basis(dimHS[0],0)).unit())

    return qtp.fidelity(last_state,qtp.fock_dm(dimHS[0],2))


def evo1T(g0,E,ops):
    ad,a=ops
    return (-1j*np.pi*g0**2*(4*E**2/3*(ad**2+a**2)-10*(ad*a)**2)).expm()

def comm(x,y):
    return x*y-y*x

def evo1T2(g0,E,ops):
    ad,a,bd,b=ops
    t1=g0**2*(4*E**2/3*(ad**2+a**2)-10*(ad*a)**2)
    t2=g0**3*8*E**2/3*(ad**2-a**2)*(b-bd)
    # t2=-4j*g0**4*E**4*4*np.pi/3*4/3*2*(ad**2-a**2)
    t3=comm(2j*np.pi*g0**2*(ad*a)**2,t1)/2
    return (-1j*np.pi*(t1+t2+t3)).expm()


def evoFull(dimHS,g0,Elist):
    if type(dimHS)==type(1):
        ops=[qtp.create(dimHS),qtp.destroy(dimHS)]
        res=qtp.basis(dimHS,0)
        for E in Elist:
            res=evo1T(g0,E,ops)*res
        return res
    else:
        ops=geOps(dimHS,True)
        res=ops[-1]
        ops=ops[:4]
        for E in Elist:
            res=evo1T2(g0,E,ops)*res
        return res.ptrace(0)

def evoFullInfid(dim,g0,Elist,fins=0):
    dimHS=dim
    if type(dimHS)==type([]):
        dim=dim[0]
    # fins=qtp.basis(dim,2)
    if not fins:
        fins=qtp.basis(dim,2)
    elif fins[2]==None:
        return np.abs(0.5-qtp.ket2dm(evoFull(dimHS,g0,Elist)).diag()[2])+np.abs(0.5-qtp.ket2dm(evoFull(dimHS,g0,Elist)).diag()[0])
    else:
        fins=(qtp.basis(dim,0)*fins[0]+qtp.basis(dim,2)*fins[1]*np.exp(1j*fins[2])).unit()
        
    return 1-qtp.fidelity(evoFull(dimHS,g0,Elist),fins)#+np.mean(Elist)/100
    
# initVals=tuple([random.uniform(0,eMax) for i in range(eLen)])
def rotFins(fins,es,g0):
    return [fins[0],fins[1],fins[2]-sum([4*np.pi*g0**2/3*8*e**2 for e in es])]
def main(params):
    g0,eMax,eLen,dimHS=params
    initVals=tuple([random.random()*eMax for i in range(eLen)])
    initVals=[4.        , 4.        , 4.        , 2.96293001,
       4.        , 2.31117902, 3.25263082, 4.        ,
       4.        , 1.90189397, 2.32084675, 4.        , 2.88453756,
       0.06322757, 3.17709161, 4.        ]
    # initVals=tuple([0 for i in range(eLen)])
    bounds=tuple([(0,eMax) for i in range(eLen)])
    # t0=timer()
    fins=[1,1,np.pi]
    print("start optimizing {}".format(dimHS))
    res=scp.optimize.dual_annealing(lambda x: evoFullInfid(dimHS,
                                                           g0,
                                                           x,
                                                           fins=rotFins(fins,x,g0)),
                                                           bounds=bounds,x0=initVals)
    # print(timer()-t0)

    # print(res.x)
    es=res.x
    output=[g0,eMax,eLen,es,1-evoFullInfid(dimHS,g0,res.x,fins=rotFins(fins,res.x,g0)),
            fidnum(list(res.x),eLen,[60,30],g0,fins=fins)]
    print(params)
    return output


if __name__ == "__main__":
    test=0
    numfid_update=0
    # print(main([1/26,4,16,10]))
    print(fid2_wrapper([[2.0, 2.0, 8.11010579e-05, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0],
                       0,1/13]))
    if test:
        # print(fidnum([6.0, 6.0, 6, 6, 6, 6, 6, 6, 6.0, 6.0],10,[120,30],1/26))
    #     es=[3.99999994, 3.99999994, 3.99999996, 3.99983028, 3.95261874,
    #    3.91088494, 3.90915495, 3.95843091, 3.99999995]
    #     g0=1/26
    #     print(fidnum(es,len(es),[60,30],g0,fins=[1,1]))
        # print(1-evoFullInfid(30,g0,es,fins=rotFins([1,1],es)))
        # print(fidnum([6.0, 6.0, 6, 6, 6, 6, 6, 6, 6.0, 6.0],10,[150,30],1/26))
        # print(1-evoFullInfid([30,10],1/26,[6.0, 6.0, 6, 6, 6, 6, 6, 6, 6.0, 6.0]))
        # print(1-evoFullInfid([50,20],1/26,[6.0, 6.0, 6, 6, 6, 6, 6, 6, 6.0, 6.0]))
        exit()
    if numfid_update:
        df=pd.read_csv('optimizedDriv.csv',converters={"eList":literal_eval },index_col=0)
        data=df[["eList","eLen","g0"]].values
        # res=qtp.parallel_map(fidnum_wrapper,data)
        res=qtp.parallel_map(fid2_wrapper,data)
        df["fid2"]=res
        try:
            df.to_csv("optimizedDriv.csv")
        except Exception as e:
            print("An error occurred:", e)
            qtp.qsave(res,"tmp")

        exit()
    if 0:
        dim=30
        params=[]
        for g0 in [1/13,1/26]:
            for eLen in list(range(2,17)):
                for eMax in [2,3,4,5,6,8,10]:
                    params.append([g0,eMax,eLen,dim])   
        res=qtp.parallel_map(main,params,num_cpus=4)
        g0,eMax,eLen,eList,fid,num_fid=list(map(list, zip(*res)))
        eList=[list(d) for d in eList]
        df=pd.DataFrame({'g0':g0,'eMax':eMax,'eLen':eLen,'eList':eList,'fid':fid,'fidnum':num_fid})
        print(df)
        try:
            df.to_csv("optimizedDriv.csv")
        except Exception as e:
            print("An error occurred:", e)
            qtp.qsave(res,"tmp")

# 
    # main([1/13,4,16])
    # print(evoFullInfid(1/13,[1.76176231 ,0.82617765, 0.64129756, 3.52407489, 2.17848505 ,1.37837895, 0.37760698 ,0.,1.18303533,
    #        2.2981967 ,2.79501525 ,3.74115919, 3.46800745, 2.88499882 ,0.89651358 ,0.13146713]))
    