from qutip import *


def geOps(dimHS,initS=False):
    if type(dimHS)==int:
        HSdimc=HSdimm=dimHS
    else:
        HSdimc,HSdimm=dimHS

    idMatc=qeye(HSdimc)
    idMatm=qeye(HSdimm)

    ad=tensor(create(N=HSdimc),idMatm)
    a=tensor(destroy(N=HSdimc),idMatm)
    bd=tensor(idMatc,create(N=HSdimm))
    b=tensor(idMatc,destroy(N=HSdimm))
    idMat=tensor(idMatc,idMatm)
    if initS==False:
        return [ad,a,bd,b,idMat]
    else:
        try:
            return [ad,a,bd,b,idMat,tensor(thermal_dm(HSdimc,initS[0]),thermal_dm(HSdimm,initS[1]))]
        except:
            return [ad,a,bd,b,idMat,tensor(basis(HSdimc,0),basis(HSdimm,0))]
