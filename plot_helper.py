from headers import *
qtp.settings.auto_tidyup = True
import matplotlib

num_cpus=4


def plotEvo(etaList, tTotal, dimHS=[60, 30], kk=1/50, ax2ylim=[0, 1]):

    trange, results = evo_stepwisePsi2(dimHS, kk, etaList, tTotal, wc=0)

    last_state = results[-1][-1]

    trange = trange[::5000]
    results = results[::5000]

    focks_to_plot = range(7)
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 3]})
    ax1, ax2 = axs

    colors = generate_cmap(int((len(focks_to_plot)+1)/2))

    for i in focks_to_plot[::2]:
        ax1.plot(trange, [np.abs(r[0][i]) for r in results],
                 '-o', label=str(i), color=colors[int(i/2)])
    for i in focks_to_plot[1::2]:
        ax2.plot(trange, [np.abs(r[0][i]) for r in results],
                 '-o', label=str(i), color=colors[int((i-1)/2)])

    ax1.plot([0, 0], [0, 1], linestyle=':', alpha=0.6, color="gray")
    ax2.plot([0, 0], [0, 0.1], linestyle=':', alpha=0.6, color="gray")
    for t in range(len(etaList)):
        xpos = trange[-1]*(t+1)/len(etaList)
        ax1.plot([xpos, xpos], [0, 1], linestyle=':', alpha=0.6, color="gray")
        ax2.plot([xpos, xpos], [0, 0.1], linestyle=':',
                 alpha=0.6, color="gray")
    ax1.plot([0, trange[-1]], [1, 1], linestyle=':', alpha=0.6, color="gray")
    ax1.set_title('(a)', fontsize=14)
    ax2.set_title('(b)', fontsize=14)
    ax1.legend(loc=6)
    ax2.legend(loc=6)

    ax1.set_xlim([0, trange[-1]+.5])
    ax2.set_xlim([0, trange[-1]+.5])
    ax1.set_ylim([0, 1])
    ax2.set_ylim(ax2ylim)
    ax1.set_ylabel("population", fontsize=14)
    ax2.set_ylabel("population", fontsize=14)
    ax2.set_xlabel(r"evolution time $(\cal T)$", fontsize=14)
    ax2.ticklabel_format(axis="y", style="sci",
                         useMathText=True, scilimits=(0, 0))
    ax2.get_yaxis().get_offset_text().set_position((-0.05, 0))
    ax1.title.set_position([0.05, .85])
    ax2.title.set_position([0.05, .75])
    fig.tight_layout(pad=1.5)

    return [ax2, last_state]


def plot3_helper(params):
    t, k, eta = params
    fin_state = evo_stepwisePsi2([60, 30], k, eta, t)[1][-1][-1]
    return [t, qtp.fidelity(qtp.fock_dm(60, 2), fin_state.ptrace(0))]


def plot3_data():
    etas = [[4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291], 
            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.17006597, 4.0, 1.22723673, 1.52223483, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.5495041, 4.0, 1.5971311, 2.36274737, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.71577367, 4.0, 1.63366753, 3.22583102, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.87411802, 4.0, 0.96209439, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.88956397, 3.85112561, 1.61025478, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.40516261, 2.13973876, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0, 3.01725078, 2.33951478, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0, 3.53009253, 0.8734978, 4.0, 4.0], [4.0, 4.0, 4.0, 3.97676322, 0.0, 4.0, 4.0], 
            [4.0, 4.0, 4.0, 1.25045869, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0, 4.0], [4.0, 4.0, 4.0, 4.0], [4.0, 4.0, 4.0], [4.0, 4.0]]
    # [
    # [4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4.,
    #     2.49797375, 3.59489053, 2.08846815, 0., 3.70603786],
    # [4., 4., 4., 4., 4.,        4., 4., 4., 4., 4.,
    #     3.38196694, 4., 0.53912259, 2.01857255, 4.],
    # [4., 4., 4., 4., 4., 4., 4., 4., 4.,
    #     3.8656867, 4., 0.69056364, 3.01814093, 4.],
    # [4., 4., 4., 4., 4.,   4., 4., 4., 4.,
    #     4.,      1.10563985, 3.80210452, 4.],
    # [4., 4., 4., 4., 4.,    4., 4., 4., 3.96051742, 1.94042786,     4., 4.],
    # [4., 4., 4., 4., 4.,    4., 4., 3.69248497, 2.69687054, 4.,   4.],
    # [4., 4., 4., 4., 4., 4.,        3.276679, 3.3070132, 4., 4.],
    # [4., 4., 4., 4., 4., 2.5040263, 4., 4., 4.],
    # [4., 4., 4., 4., 2.57041217, 4., 4., 4.],
    # [4., 4., 4., 4., 3.18409448, 4., 4.],
    # [4., 4., 4., 4., 4., 4.],
    # [4., 4., 4., 4., 4.],
    # [4., 4., 4., 4.],
    # [4., 4., 4.],
    # [4., 4.]
    # ]
    ts = [len(eta) for eta in etas]
    k = 1/26
    params = [(ts[j], k, etas[j]) for j in range(len(ts))]
    res1 = qtp.parallel_map(plot3_helper, params)
    qtp.qsave(res1, "plot31")

    etas = [[0.83093504, 0.0, 0.0, 2.10521188, 3.46665956, 2.77389036, 0.22988597, 0.1306046, 0.86739215, 0.16249172, 1.21165408, 2.85668221, 3.247835, 3.24689537, 2.85721714, 2.3522673], [0.674783, 1.98725721, 2.77445101, 3.37351876, 3.15248317, 2.86872136, 1.69225217, 0.9479075, 0.4195115, 1.69867812, 1.5718756, 2.95950184, 3.08234397, 2.57322173, 0.13098301], [2.32243241, 0.311356, 1.68925154, 3.37609797, 2.79420661, 1.8517636, 2.00062808, 1.25809964, 2.56401952, 2.45209519, 2.76029606, 3.98158015, 3.83377598, 3.02951596], [2.94037727, 2.33163804, 1.78828629, 1.53691764, 0.792807, 0.5236126, 0.8754181, 1.84818364, 3.49894809, 3.56598942, 3.4559403, 2.10826174, 1.40657165], [2.31313424, 3.00783591, 0.0, 0.0, 1.06233684, 0.0, 0.15660283, 2.61055237, 3.44946681, 3.82153135, 3.13792614, 2.59067397], [4.0, 3.78323189, 3.99953547, 3.49068232, 2.72282416, 0.32744249, 0.07638135, 1.85979568, 2.50356829, 2.12706333, 1.42569084], [3.49038164, 0.0, 0.0, 0.0, 2.75168283, 0.88558901, 4.0, 4.0, 4.0, 2.60935021], [0.0, 0.0, 4.0, 4.0, 3.76746468, 2.42546703, 1.14962237, 0.0, 0.000999952145], [4.0, 4.0, 4.0, 2.95518955, 0.97831183, 0.0, 1.0529061, 1.76085499], [4.0, 4.0, 3.76746682, 2.42533462, 1.14979777, 0.0, 0.00462229], [0.0, 3.12325688, 4.0, 4.0, 3.05750842, 1.4690436], [4.0, 4.0, 3.76745115, 2.42540346, 1.14966943],
            [4.0, 4.0, 3.69069554, 2.26722955], [4.0, 4.0, 3.39113346], [4.0, 3.47665854]]
    # [
    #     [3.9999999, 1.69462525, 2.36943397, 0., 0., 0., 2.32540925, 2.45119854,
    #         3.62888092, 3.18525243,       3.48116141, 2.37310122, 0., 0., 0.,       2.1086259],
    #     [3.9999999, 3.34487626, 2.90176591, 0., 1.9236662, 1.98116924, 1.45652025,
    #         0., 2.92582452, 2.9967507, 3.1087518, 2.24876703, 2.21506765, 0., 0.],
    #     [3.9999999, 3.48124632, 2.89161745, 0.91875404, 0.4765506, 0.10041421, 0.,
    #         0.02332447, 2.0193055, 2.86899834, 2.38379075, 2.01947063, 1.29743498, 1.72375327],
    #     [3.9999999, 3.54116719, 2.76939085, 0., 0., 0.34818725, 0., 2.22194408,
    #         2.76459647, 2.32878124, 2.50157588, 2.75529737, 2.56162155],
    #     [3.9999999, 3.40378334, 2.78335668, 2.45494324, 2.19553132, 0.23313043,
    #         0.13510194, 2.092391, 2.00891672, 2.53748046, 2.81992011, 2.64037264],
    #     [3.9999999, 3.91717001, 3.7122676, 2.83006775, 2.55580864,
    #         1.95958183, 0., 0.14030912, 2.23447418, 2.55981859, 1.8119227],
    #     [4., 3.97257957, 4., 3.1996309, 2.4839237,
    #         0., 0., 1.83584526, 2.30775097, 1.69335884],
    #     [4., 4., 3.90880801, 2.98339554, 0.,       0., 0., 1.83714234, 0.],
    #     [4., 4., 4., 2.98911131, 0.81626507, 0., 0.89465581, 1.81979884],
    #     [4., 4., 3.80381825, 2.41662053, 1.01843143, 0., 0.],
    #     [4., 4., 3.80381905, 2.41661554, 1.01837679,       0.],
    #     [4., 4., 3.80381696, 2.41661958, 1.01840755],
    #     [4., 4., 3.76601341, 2.35615187],
    #     [4., 4., 3.57082097],
    #     [4., 3.95031349]
    # ]
    k = 1/13
    params = [(ts[j], k, etas[j]) for j in range(len(ts))]
    res2 = qtp.parallel_map(plot3_helper, params)
    qtp.qsave(res2, "plot32")

    ks = [26, 25, 24.5, 23.5, 22.5, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13]
    ks = [1/k for k in ks]
    etas = [
        [4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4.,
            2.49797375, 3.59489053, 2.08846815, 0., 3.70603786],
        [4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 2.98076997,
            3.18952502, 2.43134345, 1.18002863, 3.31344385],
        [4., 4., 4., 4., 4., 4., 4., 4., 4., 3.68650052,
            2.990641,  2.89324996, 0., 3.63966483],
        [4., 4., 4., 4., 4., 4., 4., 4., 3.93645222,
            2.84478589, 3.17641015, 0., 3.46080055],
        [4., 4., 4., 4., 4., 4., 4., 4., 3.09525006,
            3.23263751, 0.43243956, 3.32943828],
        [4., 4., 4., 4., 4., 4., 4., 3.68606787, 3.37912976, 0.,     3.63974771],
        [4., 4., 4., 4., 4., 4., 3.91014959,
            3.42725666, 1.29097229, 3.47449182],
        [4., 4., 4., 4., 4.,    4., 3.66716892, 1.66618426, 3.37162242],
        [4., 4., 4., 4., 4.,   3.94224786, 2.04301788, 3.29227687],
        [4., 4., 4., 4., 4.,   2.67924787, 3.13455576],
        [4., 4., 4., 4., 3.21014382,        3.09264885],
        [4., 4., 4., 3.563877, 3.28142102],
        [4., 4., 3.7853928, 3.58642551],
        [4.0, 4.0, 3.39113346],
        [4.0, 3.47665854]
    ]
    params = [(ts[j], ks[j], etas[j]) for j in range(len(ts))]
    res3 = qtp.parallel_map(plot3_helper, params)
    qtp.qsave(res3, "plot33")
    return [res1, res2, res3]


def plot678HelperNum(params):
    k, eta, t, state2, dimHS, imperfect = params
    targetStates = [qtp.fock_dm(dimHS[0], 2), state2]
    if type(imperfect) == type([]):
        expects = evo_stepwisePsi2(
            dimHS, k, eta, t, loss=imperfect, show_diags=False)[1]
    else:
        expects = evo_stepwisePsi2(
            dimHS, k, eta, t, thermalOccu=imperfect, show_diags=False)[1]
    finState = (expects[-1][-1]).ptrace(0)
    res = params
    for s in targetStates:
        res.append(qtp.fidelity(s, finState))
    print("{} finished".format([k,eta,t,dimHS,imperfect]))
    print(res[-2:])
    return res

def plot678HelperTrot(params):
    k, eta, t, state2, dimHS, imperfect = params
    print("{} started".format([k, eta, dimHS, imperfect]))
    targetStates = [qtp.fock_dm(dimHS[0], 2), state2]
    if type(imperfect) == type([]):
        fins = evo_trot(dimHS, k,eta,loss=imperfect,evoOrder=5,recal=1)
    else:
        fins = evo_trot(dimHS, k,eta,initOccu=[0,imperfect],evoOrder=5,recal=1)
    finState = fins.ptrace(0)
    res = params
    for s in targetStates:
        res.append(qtp.fidelity(s, finState))
    print("{} finished".format([k, eta, dimHS, imperfect]))
    return res

def plot678base(dimHS):
    res=-1
    try:
        return qtp.qload("plot678base{}".format(dimHS))
    except FileNotFoundError:

        params0 = [[1/26,
                [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                    2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
                16]]
        params0.append([1/21,
                    [4., 4., 4., 4., 4., 4., 3.91014959,
                        3.42725666, 1.29097229, 3.47449182],
                    10])
        params0.append([1/16,
                    [4., 4., 4., 3.563877, 3.28142102],
                    5])
        res= [evo_trot(dimHS, p[0],p[1],evoOrder=5,recal=0).ptrace(0) for p in params0]
        qtp.qsave(res,"plot678base{}".format(dimHS))
    return res

def plot6_data(dimHS=[60,60]):
    
    params = []
    params = [[1/26,
            #    [4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4.,
            #        2.49797375, 3.59489053, 2.08846815, 0., 3.70603786],
            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
             2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16]]
    params.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10])
    params.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5])

    initNList = np.arange(0, 5+1e-10, 0.1)
    targs_magnus=plot678base(dimHS)
    for j in range(len(params)):
        # try:
        #     qtp.qload("{}plot6{}".format(dimHS,j))
        # except FileNotFoundError:
        params[j].append(targs_magnus[j%3])
        full_params = [params[j]+[dimHS, i] for i in initNList]
        res = qtp.parallel_map(plot678HelperTrot, full_params,num_cpus=num_cpus)
        print("{}plot6{}Finished\n".format(dimHS,j))
        qtp.qsave(res, "{}plot6{}".format(dimHS,j))
        print("{}plot6{}saved\n".format(dimHS,j))

    return "plot6 finished"


def plot7_data(dimHS=[60,30]):
    params = []
    params = [[1/26,
               [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
             2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16]]
    params.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10])
    params.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5])

    lossRange = np.logspace(-5, -2.7, 20)
    targs_magnus=plot678base(dimHS)
    for j in range(len(params)):
        # try:
        #     qtp.qload("{}plot7{}".format(dimHS,j))
        # except FileNotFoundError:
        params[j].append(targs_magnus[j%3])
        full_params = [params[j]+[dimHS, [i, 0, 0]] for i in lossRange]
        res = qtp.parallel_map(plot678HelperTrot, full_params,num_cpus=num_cpus)
        print("{}plot7{}Finished\n".format(dimHS,j))
        qtp.qsave(res, "{}plot7{}".format(dimHS,j))
        print("{}plot7{}saved\n".format(dimHS,j))

    return "plot7 finished"


def plot8_data(dimHS=[60,60]):
    

    params = [[1/26,
               [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16, 1]]
    params.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 1])
    params.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 1])
    params.append([1/26,
                   [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                    2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
                   16, 10])
    params.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 10])
    params.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 10])
    params.append([1/26,
                   [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                    2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
                   16, 100])
    params.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 100])
    params.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 100])

    lossRange = np.logspace(-5, -2.7, 20)
    targs_magnus=plot678base(dimHS)
    for j in range(len(params)):
        # try:
        #     qtp.qload("{}plot8{}".format(dimHS,j))            
        #     print("{}plot8{} finished".format(dimHS,j))
        # except FileNotFoundError:
        thermalN=params[j][-1]
        params[j]=params[j][:-1]
        params[j].append(targs_magnus[j%3])
        full_params = [params[j]+[dimHS, [0, i, thermalN]]
                    for i in lossRange]
        res = qtp.parallel_map(plot678HelperTrot, full_params,num_cpus=num_cpus)
        print("{}plot8{}Finished\n".format(dimHS,j))
        qtp.qsave(res, "{}plot8{}".format(dimHS,j))
        print("{}plot8{}saved\n".format(dimHS,j))

    return "plot8 finished"


def plot678_data_sample():
    dimHS = [60, 60]
    params1 = []
    params1 = [[1/26,
            #    [4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4.,
            #        2.49797375, 3.59489053, 2.08846815, 0., 3.70603786],
            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
             2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16]]
    params1.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10])
    params1.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5])

    initNList = [1,2,5]
    full_params1=[]
    for j in range(len(params1)):
        params1[j].append((evo_stepwisePsi2(dimHS, params1[j][0],
                         params1[j][1], params1[j][2])[1][-1][-1]).ptrace(0))
        for i in initNList:
            full_params1.append(params1[j]+[dimHS, i])

    dimHS = [60, 30]
    params2 = []
    params2.append([1/26,
               [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
             2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16])
    params2.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5])
    params2.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10])

    lossRange = [1e-5,1e-4,1e-3]
    full_params2=[]
    for j in range(len(params2)):
        params2[j].append((evo_stepwisePsi2(dimHS, params2[j][0],
                         params2[j][1], params2[j][2])[1][-1][-1]).ptrace(0))
        for i in lossRange:
            full_params2.append(params2[j]+[dimHS, [i, 0, 0]])

    dimHS = [60, 60]
    params3 = [[1/26,
               [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
               16, 1]]
    params3.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 1])
    params3.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 1])
    params3.append([1/26,
                   [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                    2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
                   16, 10])
    params3.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 10])
    params3.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 10])
    params3.append([1/26,
                   [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.97615762, 
                    2.43055998, 3.6873501, 1.95440809, 0.04559224, 3.7363291],
                   16, 100])
    params3.append([1/21,
                   [4., 4., 4., 4., 4., 4., 3.91014959,
                       3.42725666, 1.29097229, 3.47449182],
                   10, 100])
    params3.append([1/16,
                   [4., 4., 4., 3.563877, 3.28142102],
                   5, 100])

    lossRange = [1e-5,1e-3]

    full_params3=[]
    for j in range(len(params3)):
        thermalN=params3[j][-1]
        params3[j]=params3[j][:-1]
        params3[j].append((evo_stepwisePsi2(dimHS, params3[j][0],
                         params3[j][1], params3[j][2])[1][-1][-1]).ptrace(0))
        for i in lossRange:
            full_params3.append(params3[j]+[dimHS, [0, i, thermalN]])
        
    fullparams=full_params1+full_params2+full_params3
    print("calculation started")
    # print([d[-3].ptrace(0).diag()[2] for d in fullparams])
    # exit()
    res = qtp.parallel_map(plot678HelperNum, fullparams,num_cpus=num_cpus)
    qtp.qsave(res, "plot678samples")

    return "sample_finished"

if __name__ == "__main__":
    0
    # print(plot3_data())
    # print(plot6_data([60,25]))
    # print(plot7_data([60,25]))
    # print(plot8_data([60,25]))
    # print(plot6_data([60,60]))
    # print(plot7_data([60,30]))
    # print(plot8_data([60,60]))
    print(plot678_data_sample())
    # dimHS=[20,20]
    # params=[1/26,
    #         [4.,4.,4.,4.,4.,4., 4.,4.,4.,4.,4.,2.49797375 ,3.59489053 ,2.08846815 ,0.,3.70603786],
    #         16]

    # lossRange=np.logspace(-5,-2,50)
    # loss=lossRange[23]
    # print(plot678Helper(params+[(evo_stepwisePsi(dimHS, params[0], params[1], params[2])[1][-1][-1]).ptrace(0),
    #                             dimHS,[loss,0,0]]))
# SimpleMode=True #True if all timestamps are integers
# SAVE=0
# RECALCULATE=0
# METHOD_SIMUL='numeric' #numeric or theoretic
# PLOT1='fidelity' #eigenvalue or fidelity
# TARGET_STATE='fock' #superpos or fock
# SMOOTH_PHASE_MODULATION=[0,1e-4]
# # dimHS=30
# count = 0
# # k=1/22
# # chain =[4.        , 4.        , 4.        , 4.        , 4.        ,
# #        4.        , 3.81604926, 3.57817146, 0.92920533, 3.56379942]#k=1/30
# # chain=[4.,4.,4.,4.,3.98279931 ,3.0684555,2.40084713, 2.02849304] #k=1/25
# # chain=[4.        , 4.        , 4.        , 4.        , 3.95853402,
#     #    2.51821434, 2.90380013] #k=1/25
# # chain=[4.        , 4.        , 4.        , 4.        , 4.        ,
# #        3.56718706, 2.18928661, 2.74160269] #k=1/26
# # chain=[4.      ]#k=1/26 test

# if TARGET_STATE=='fock':
#     k=1/26
#     chain=[4.      ,   4.    ,     4.      ,   4.       ,  4.   ,      3.56,
#             2.19 ,2.74]#k=1/26 better model
#     # chain=[3.56125]*8
#     pointsToPlot = 2*len(chain)
# elif TARGET_STATE=='superpos':
#     k=1/26
#     chain=[4., 3.69291531, 3.51161495, 3.686096 , 2.60729845, 3.19989213,
#            3.13977697,1.51847303,2.00439795,3.9995126] #equal superposition
#     pointsToPlot = len(chain)
# elif TARGET_STATE=='test':
#     k=1/26
#     chain=[2.20000000e+00, 1.79168856e+00, 1.96560239e-04, 1.03257383e+00,
#        1.85834411e+00, 1.50822955e+00, 9.63259222e-01, 1.48973580e+00,
#        1.45596352e+00, 2.45626944e-01, 1.30634134e-03, 1.90094975e+00]
#     pointsToPlot = len(chain)
#     TARGET_STATE_S=[0,1]
# else:
#     input("invalid target")

# focks_to_plot=range(7)

# plotChain(chain, pointsToPlot, k)
