function _torord_rush_larsen_impl!(u_new::AbstractVector{T}, u::AbstractVector{T}, parameters::AbstractVector{T}, celltype, x, t::T, dt::T, spatial_funcs::F) where {T, F}

    # Assign states
    #hL = u[1]
    hL = u[63]
    jca = u[2]
    m = u[3]
    mL = u[4]
    xs1_P = u[5]
    xs1 = u[6]
    a = u[7]
    ap = u[8]
    CaMKt = u[9]
    hLp = u[10]
    h_P = u[11]
    h = u[12]
    hp_P = u[13]
    hp = u[14]
    j = u[15]
    j_P = u[16]
    xs2 = u[17]
    C0 = u[18]
    C1 = u[19]
    O_ = u[20]
    clss = u[21]
    cli = u[22]
    iF = u[23]
    iS = u[24]
    TmBlocked = u[25]
    jp_P = u[26]
    jp = u[27]
    d_P = u[28]
    d = u[29]
    fcaf_P = u[30]
    fcaf = u[31]
    ff_P = u[32]
    ff_ = u[33]
    fcas_P = u[34]
    fcas = u[35]
    fs_P = u[36]
    fs = u[37]
    C2 = u[38]
    I_ = u[39]
    iFp = u[40]
    iSp = u[41]
    Zetas = u[42]
    Zetaw = u[43]
    Ca_TRPN = u[44]
    fcaBPf = u[45]
    fcafp = u[46]
    fBPf = u[47]
    ffp = u[48]
    XS = u[49]
    XW = u[50]
    nca_i = u[51]
    nca_ss = u[52]
    cajsr = u[53]
    cansr = u[54]
    kss = u[55]
    ki = u[56]
    Jrel_np = u[57]
    Jrel_p = u[58]
    cai = u[59]
    nai = u[60]
    cass = u[61]
    nass = u[62]
    #v = u[63]
    v = u[1]
    Jrel_p_P = u[64]
    Jrel_np_P = u[65]

    # Assign parameters
    Aff = parameters[1]
    BSLmax = parameters[2]
    BSRmax = parameters[3]
    CaMKo = parameters[4]
    EKshift = parameters[5]
    F_ = parameters[6]
    Fjunc = parameters[7]
    GClCa = parameters[8]
    GClb = parameters[9]
    GK1_b = parameters[10]
    GKb_b = parameters[11]
    GKr_b = parameters[12]
    GKs_b = parameters[13]
    GNa = parameters[14]
    GNaL_b = parameters[15]
    Gncx_b = parameters[16]
    GpCa = parameters[17]
    Gto_b = parameters[18]
    H = parameters[19]
    ICaL_Multiplier = parameters[20]
    ICaL_fractionSS = parameters[21]
    ICab_Multiplier = parameters[22]
    IClCa_Multiplier = parameters[23]
    IClb_Multiplier = parameters[24]
    IK1_Multiplier = parameters[25]
    IKb_Multiplier = parameters[26]
    IKr_Multiplier_spatial = F !== Nothing ? T(spatial_funcs.IKr_Multiplier(x, t)) : one(T)
    IKr_Multiplier = IKr_Multiplier_spatial * parameters[27]
    IKs_Multiplier = parameters[28]
    INaCa_Multiplier = parameters[29]
    INaCa_fractionSS = parameters[30]
    INaK_Multiplier = parameters[31]
    INaL_Multiplier = parameters[32]
    INa_Multiplier = parameters[33]
    INab_Multiplier = parameters[34]
    IpCa_Multiplier = parameters[35]
    Ito_Multiplier = parameters[36]
    Jrel_Multiplier = parameters[37]
    Jrel_b = parameters[38]
    Jup_Multiplier = parameters[39]
    KdClCa = parameters[40]
    Khp = parameters[41]
    Kki = parameters[42]
    Kko = parameters[43]
    KmBSL = parameters[44]
    KmBSR = parameters[45]
    KmCaAct = parameters[46]
    KmCaM = parameters[47]
    KmCaMK = parameters[48]
    KmCap = parameters[49]
    Kmgatp = parameters[50]
    Kmn_b = parameters[51]
    Knai0_b = parameters[52]
    Knai0_np = parameters[53]
    Knao0 = parameters[54]
    Knap = parameters[55]
    Kxkur = parameters[56]
    L = parameters[57]
    MgADP = parameters[58]
    MgATP = parameters[59]
    PCa_P_b = parameters[60]
    PCa_b = parameters[61]
    PCab = parameters[62]
    PKNa = parameters[63]
    PNab = parameters[64]
    Pnak_b = parameters[65]
    Q10CaL = parameters[66]
    Q10ICaL_a = parameters[67]
    Q10ICaL_ff = parameters[68]
    Q10ICaL_fs = parameters[69]
    Q10K = parameters[70]
    Q10Kb = parameters[71]
    Q10NCX = parameters[72]
    Q10NaK = parameters[73]
    Q10SLCaP = parameters[74]
    Q10SRCaP = parameters[75]
    R = parameters[76]
    T_base = parameters[77]
    T_val = F !== Nothing && hasproperty(spatial_funcs, :T) ? T(spatial_funcs.T(x, t)) : T_base
    TOT_A = parameters[78]
    TRPN_n = parameters[79]
    Tref_b = parameters[80]
    Whole_cell_PP1 = parameters[81]
    aCaMK = parameters[82]
    alpha_1 = parameters[83]
    bCaMK = parameters[84]
    beta_0 = parameters[85]
    beta_1 = parameters[86]
    beta_1_mech = parameters[87]
    bt = parameters[88]
    ca50_b = parameters[89]
    cajsr_half = parameters[90]
    cao = parameters[91]
    celltype_val = F !== Nothing ? T(spatial_funcs.celltype(x, t)) : T(celltype)
    clo = parameters[93]
    cmdnmax_b = parameters[94]
    csqnmax = parameters[95]
    delta = parameters[96]
    dielConstant = parameters[97]
    dr = parameters[98]
    eP = parameters[99]
    fICaLP = parameters[100]
    fIKsP = parameters[101]
    fIKurP = parameters[102]
    fINaKP = parameters[103]
    fINaP = parameters[104]
    fPLBP = parameters[105]
    fRyRP = parameters[106]
    fTnIP = parameters[107]
    gamma = parameters[108]
    gamma_wu = parameters[109]
    i_Stim_Amplitude = F !== Nothing && hasproperty(spatial_funcs, :stim) ? T(spatial_funcs.stim(x, t)) : parameters[110]
    i_Stim_Period = parameters[111]
    i_Stim_PulseDuration = parameters[112]
    i_Stim_Start = parameters[113]
    isHypoxic = F !== Nothing ? T(spatial_funcs.isHypoxic(x, t)) : zero(T)
    celltype_val = F !== Nothing ? T(spatial_funcs.celltype(x, t)) : T(celltype)
    jsrMidpoint = parameters[115]
    k1m = parameters[116]
    k1p = parameters[117]
    k2m = parameters[118]
    k2n = parameters[119]
    k2p = parameters[120]
    k3m = parameters[121]
    k3p = parameters[122]
    k4m = parameters[123]
    k4p = parameters[124]
    k_uw = parameters[125]
    k_ws = parameters[126]
    kasymm = parameters[127]
    kcaoff = parameters[128]
    kcaon = parameters[129]
    kmcmdn = parameters[130]
    kmcsqn = parameters[131]
    kmtrpn_b = parameters[132]
    kna1 = parameters[133]
    kna2 = parameters[134]
    kna3 = parameters[135]
    ko = parameters[136]
    koff = parameters[137]
    ktm_unblock = parameters[138]
    lambda = parameters[139]
    lambda_max = parameters[140]
    lambda_min = parameters[141]
    lambda_rate = parameters[142]
    mode = parameters[143]
    mu_b = parameters[144]
    nK1 = parameters[145]
    nNaCa = parameters[146]
    nTnI = parameters[147]
    nao = parameters[148]
    nperm_b = parameters[149]
    nrel = parameters[150]
    nu_b = parameters[151]
    nup = parameters[152]
    offset = parameters[153]
    pH_base = parameters[154]
    pH = F !== Nothing && hasproperty(spatial_funcs, :pH) ? T(spatial_funcs.pH(x, t)) : pH_base
    perm50 = parameters[155]
    ph_bt = parameters[156]
    phi = parameters[157]
    pkK1 = parameters[158]
    pkNaCa = parameters[159]
    pkTnI = parameters[160]
    pkrel = parameters[161]
    pkup = parameters[162]
    qca = parameters[163]
    qna = parameters[164]
    rad_ = parameters[165]
    tauCa = parameters[166]
    tauK = parameters[167]
    tauNa = parameters[168]
    thL = parameters[169]
    tjca = parameters[170]
    trpnmax = parameters[171]
    vShift = parameters[172]
    wca = parameters[173]
    wfrac = parameters[174]
    wna = parameters[175]
    wnaca = parameters[176]
    zca = parameters[177]
    zcl = parameters[178]
    zk = parameters[179]
    zna = parameters[180]

    # Assign expressions
    GNa_P = T(1.7) * T(11.7802)
    A = (dr / T(0.25)) .* ((T(0.25) * TOT_A) ./ (dr + wfrac .* (T(1) - dr)))
    XSSS = T(0.5) * dr
    XWSS = T(0.5) * (wfrac .* (T(1) - dr))
    Afcaf =
        T(0.3) + T(0.6) ./ (exp((v - T(10.0) / T(1)) / T(10.0)) + T(1.0))
    ah = (
        if (v > T(40.0) * (T(-1)))
            (T(0.0))
        else
            (T(4.43126792958051e-07) * exp((T(-0.147058823529412)) * v))
        end
    )
    aj = (
        if (v > T(40.0) * (T(-1)))
            (T(0.0))
        else
            (
                (
                    (
                        (-(v + T(37.78))) .*
                        (T(25428.0) * exp(T(0.28831) * v) + T(6.948e-06))
                    ) .* exp((T(-0.04391)) * v)
                ) ./ (T(50262745825.954) * exp(T(0.311) * v) + T(1.0))
            )
        end
    )
    bh = (
        if (v > T(40.0) * (T(-1)))
            (
                (T(0.77) * exp(T(0.0900900900900901) * v)) ./
                (T(0.13) * exp(T(0.0900900900900901) * v) + T(0.0497581410839387))
            )
        else
            (T(2.7) * exp(T(0.079) * v) + T(310000.0) * exp(T(0.3485) * v))
        end
    )
    bj = (
        if (v > T(40.0) * (T(-1)))
            (
                (T(0.6) * exp(T(0.157) * v)) ./
                (T(1.0) * exp(T(0.1) * v) + T(0.0407622039783662))
            )
        else
            (
                (T(0.02424) * exp(T(0.12728) * v)) ./
                (T(1.0) * exp(T(0.1378) * v) + T(0.00396086833990426))
            )
        end
    )
    dPss_b = T(1.0323) * exp((T(-1.0553)) * exp((T(-0.081)) * (v + T(9.5))))
    dss = (
        if (v >= T(31.4978))
            (T(1.0))
        else
            (T(1.0763) * exp((T(-1.007)) * exp((T(-0.0829)) * v)))
        end
    )
    fss = T(1.0) ./ (exp((v + T(19.58)) / T(3.696)) + T(1.0))
    fss_P = T(1.0) ./ (exp(((v + T(19.58)) + T(4.0)) / T(3.696)) + T(1.0))
    hLss = T(1.0) ./ (exp((v + T(87.61)) / T(7.488)) + T(1.0))
    hLssp = T(1.0) ./ (exp((v + T(93.81)) / T(7.488)) + T(1.0))
    hss = T(1.0) ./ (exp((v + T(71.55)) / T(7.43)) + T(1.0)) .^ T(2)
    hss_P =
        T(1) ./
        (T(1) * (exp(((v + T(71.55)) + T(5.0)) / T(7.43)) + T(1)) .^ T(2))
    hssp = T(1.0) ./ (exp((v + T(77.55)) / T(7.43)) + T(1.0)) .^ T(2)
    hssp_P =
        T(1) ./ (
            T(1) *
            (exp((((v + T(71.55)) + T(6)) + T(5.0)) / T(7.43)) + T(1)) .^
            T(2)
        )
    jcass = T(1.0) ./ (exp((v + T(18.08)) / T(2.7916)) + T(1.0))
    mLss = T(1.0) ./ (exp((-(v + T(42.85))) / T(5.264)) + T(1.0))
    mss = T(1.0) ./ (exp((-(v + T(56.86))) / T(9.03)) + T(1.0)) .^ T(2)
    taum =
        T(0.06487) * exp(-((v - T(4.823) / T(1)) / T(51.12)) .^ T(2)) +
        T(0.1292) * exp(-((v + T(45.79)) / T(15.54)) .^ T(2))
    tmL =
        T(0.06487) * exp(-((v - T(4.823) / T(1)) / T(51.12)) .^ T(2)) +
        T(0.1292) * exp(-((v + T(45.79)) / T(15.54)) .^ T(2))
    txs1 =
        T(817.3) +
        T(1.0) ./ (
            T(0.0002326) * exp((v + T(48.28)) / T(17.8)) +
            T(0.001292) * exp((-(v + T(210.0))) / T(230.0))
        )
    txs1_P =
        (
            T(817.3) -
            T(1.75) / (
                T(0.001292) * exp((-(T(10) + T(210.0))) / T(230.0)) +
                T(0.0002326) * exp((T(10) + T(48.28)) / T(17.8))
            )
        ) +
        T(2.75) ./ (
            T(0.0002326) * exp((v + T(48.28)) / T(17.8)) +
            T(0.001292) * exp((-(v + T(210.0))) / T(230.0))
        )
    txs2 =
        T(1.0) ./ (
            T(0.01) * exp((v - T(50.0) / T(1)) / T(20.0)) +
            T(0.0193) * exp((-(v + T(66.54))) / T(31.0))
        )
    xkb = T(1.0) ./ (exp((-(v - T(10.8968) / T(1))) / T(23.9871)) + T(1.0))
    xs1ss = T(1.0) ./ (exp((-(v + T(11.6))) / T(8.932)) + T(1.0))
    Afs = T(1.0) - Aff
    Ageo = L .* ((T(2) * T(3.14)) * rad_) + rad_ .* ((T(2) * T(3.14)) * rad_)
    vcell = L .* (rad_ .* ((T(3.14) * T(1000.0)) * rad_))
    AiF =
        T(1.0) ./
        (exp(((EKshift + v) - T(213.6) / T(1)) / T(151.2)) + T(1.0))
    ass =
        T(1.0) ./
        (exp((-((EKshift + v) - T(14.34) / T(1))) / T(14.82)) + T(1.0))
    assp =
        T(1.0) ./
        (exp((-((EKshift + v) - T(24.34) / T(1))) / T(14.82)) + T(1.0))
    dti_develop =
        T(1.354) +
        T(0.0001) ./ (
            exp((-((EKshift + v) - T(12.23) / T(1))) / T(0.2154)) +
            exp(((EKshift + v) - T(167.4) / T(1)) / T(15.89))
        )
    dti_recover =
        T(1.0) -
        T(0.5) ./ (exp(((EKshift + v) + T(70.0)) / T(20.0)) + T(1.0))
    iss = T(1.0) ./ (exp(((EKshift + v) + T(43.94)) / T(5.711)) + T(1.0))
    ta =
        T(1.0515) ./ (
            T(1.0) ./ ((
                T(1.2089) * (
                    exp((-((EKshift + v) - T(18.4099) / T(1))) / T(29.3814)) +
                    T(1.0)
                )
            )) +
            T(3.5) ./ (exp(((EKshift + v) + T(100.0)) / T(29.3814)) + T(1.0))
        )
    tiF_b =
        T(4.562) +
        T(1.0) ./ (
            T(0.3933) * exp((-((EKshift + v) + T(100.0))) / T(100.0)) +
            T(0.08004) * exp(((EKshift + v) + T(50.0)) / T(16.59))
        )
    tiS_b =
        T(23.62) +
        T(1.0) ./ (
            T(0.001416) * exp((-((EKshift + v) + T(96.52))) / T(59.05)) +
            T(1.78e-08) * exp(((EKshift + v) + T(114.1)) / T(8.079))
        )
    KsCa = T(1.0) + T(0.6) ./ ((T(3.8e-05) ./ cai) .^ T(1.4) + T(1.0))
    Bcajsr = T(1.0) ./ ((csqnmax .* kmcsqn) ./ (cajsr + kmcsqn) .^ T(2) + T(1.0))
    Bcass =
        T(1.0) ./ (
            (BSLmax .* KmBSL) ./ (KmBSL + cass) .^ T(2) +
            ((BSRmax .* KmBSR) ./ (KmBSR + cass) .^ T(2) + T(1.0))
        )
    CaMKb = (CaMKo .* (T(1.0) - CaMKt)) ./ (KmCaM ./ cass + T(1.0))
    Qpow = (T_val - T(310)) / T(10)
    ECl = ((R .* T_val) ./ ((F_ .* zcl))) .* log(clo ./ cli)
    vffrt = (F_ .* (F_ .* v)) ./ ((R .* T_val))
    vfrt = (F_ .* v) ./ ((R .* T_val))
    EClss = ((R .* T_val) ./ ((F_ .* zcl))) .* log(clo ./ clss)
    EK = ((R .* T_val) ./ ((F_ .* zk))) .* log(ko ./ ki)
    Ii = (T(0.5) * (T(4.0) * cai + (cli + (ki + nai)))) / T(1000.0)
    EKs = ((R .* T_val) ./ ((F_ .* zk))) .* log((PKNa .* nao + ko) ./ (PKNa .* nai + ki))
    ENa = ((R .* T_val) ./ ((F_ .* zna))) .* log(nao ./ nai)
    GClCa_1 = GClCa .* IClCa_Multiplier
    delta_epi = (
        if (celltype_val == T(1.0))
            (
                T(1.0) -
                T(0.95) ./ (exp(((EKshift + v) + T(70.0)) / T(5.0)) + T(1.0))
            )
        else
            (T(1.0))
        end
    )
    km2n = T(1.0) * jca
    ICaL_Multiplier_1 = (
        (isHypoxic == T(1)) ? (T(0.75) * ICaL_Multiplier) : (ICaL_Multiplier)
    )
    Io = (T(0.5) * (T(4.0) * cao + (clo + (ko + nao)))) / T(1000.0)
    IKs_Multiplier_1 = (
        (isHypoxic == T(1)) ? ((T(22) * IKs_Multiplier) / T(30)) : (IKs_Multiplier)
    )
    INaL_Multiplier_1 = (
        (isHypoxic == T(1)) ? (T(1.5) * INaL_Multiplier) : (INaL_Multiplier)
    )
    INa_Multiplier_1 = (
        (isHypoxic == T(1)) ? (T(0.9) * INa_Multiplier) : (INa_Multiplier)
    )
    Iss = (T(0.5) * (T(4.0) * cass + (clss + (kss + nass)))) / T(1000.0)
    Istim = (
        if (
            i_Stim_PulseDuration >=
            t + (
                (-i_Stim_Period) .* floor((-(i_Stim_Start - t)) ./ i_Stim_Period) -
                i_Stim_Start
            ) && i_Stim_Start <= t
        )
            (i_Stim_Amplitude)
        else
            (T(0.0))
        end
    )
    Jdiff = (-cai + cass) ./ tauCa
    JdiffCl = (-cli + clss) ./ tauNa
    JdiffNa = (-nai + nass) ./ tauNa
    JdiffK = (-ki + kss) ./ tauK
    Jtr = (-cajsr + cansr) / T(60.0)
    Knai0_P = T(0.7) * Knai0_np
    P = eP ./ (((H ./ Khp + T(1.0)) + nai ./ Knap) + ki ./ Kxkur)
    XU = -XS + (-XW + (T(1) - TmBlocked))
    a2 = k2p
    a4 = ((MgATP .* k4p) ./ Kmgatp) ./ (T(1.0) + MgATP ./ Kmgatp)
    btp = T(1.25) * bt
    tau_rel_P_b = (T(0.75) * bt) ./ (T(1.0) + T(0.0123) ./ cajsr)
    tau_rel_b = bt ./ (T(1.0) + T(0.0123) ./ cajsr)
    allo_i = T(1.0) ./ ((KmCaAct ./ cai) .^ T(2) + T(1.0))
    allo_ss = T(1.0) ./ ((KmCaAct ./ cass) .^ T(2) + T(1.0))
    b1 = MgADP .* k1m
    betaCaMKII = bCaMK .* ((T(0.9) * Whole_cell_PP1) / T(0.1371) + T(0.1))
    cmdnmax = ((celltype_val == T(1.0)) ? (T(1.3) * cmdnmax_b) : (cmdnmax_b))
    constA = T(1820000.0) ./ (T_val .* dielConstant) .^ T(1.5)
    fICaL_P = fICaLP
    fINaK_PKA = fINaKP
    fINa_P = fINaP
    fJrel_PKA = fRyRP
    fJup_P = fPLBP
    gamma_rate =
        gamma .* (
            if (
                (Zetas >= T(-1) || Zetas <= T(0) || Zetas > -Zetas - T(1)) && (
                    Zetas > T(0) && Zetas < T(-1) ||
                    (Zetas >= T(-1) || Zetas > T(-1)) && Zetas < T(-1) ||
                    Zetas > T(0)
                )
            )
                (Zetas .* ((Zetas > T(0)) ? (T(1)) : (T(0))))
            else
                ((-Zetas - T(1)) .* ((Zetas < T(-1)) ? (T(1)) : (T(0))))
            end
        )
    gamma_rate_w = gamma_wu .* abs(Zetaw)
    h4_i = (nai ./ kna1) .* (T(1.0) + nai ./ kna2) + T(1.0)
    h4_ss = (nass ./ kna1) .* (T(1.0) + nass ./ kna2) + T(1.0)
    h10_i = (nao ./ kna1) .* (T(1.0) + nao ./ kna2) + (kasymm + T(1.0))
    h10_ss = (nao ./ kna1) .* (T(1.0) + nao ./ kna2) + (kasymm + T(1.0))
    k2_i = kcaoff
    k2_ss = kcaoff
    k5_i = kcaoff
    k5_ss = kcaoff
    kmtrpn = kmtrpn_b .* (T(1.6) * fTnIP) + kmtrpn_b .* (T(1) - fTnIP)
    lambda0 = ((lambda > lambda_max) ? (lambda_max) : (lambda))
    mu = ((mode == T(1)) ? (T(1)) : (mu_b))
    nca = nca_ss
    nperm = ((mode == T(1)) ? (T(2.2)) : (nperm_b))
    nu = ((mode == T(1)) ? (T(1)) : (nu_b))
    pH_Multiplier_IK1 = (T(1.0) * T(1.2453)) ./ (exp10(nK1 .* (-pH + pkK1)) + T(1))
    pH_Multiplier_Jup = (T(1.0) * T(3.3779)) ./ (exp10(nup .* (-pH + pkup)) + T(1))
    pH_Multiplier_NCX =
        (T(1.0) * T(2.474)) ./ (exp10(nNaCa .* (-pH + pkNaCa)) + T(1))
    pH_Multiplier_RyR =
        (T(1.0) * T(1.0897)) ./ (exp10(nrel .* (-pH + pkrel)) + T(1))
    pH_Multiplier_TnICa =
        (exp10(nTnI .* (-pH + pkTnI)) + T(1)) ./
        (exp10(nTnI .* (pkTnI - T(7.2))) + T(1))
    pH_Multiplier_maxTension = ph_bt .* (pH - T(7.2)) + T(1)
    thLp = T(3.0) * thL
    Afcas = T(1.0) - Afcaf
    tauh = T(1.0) ./ (ah + bh)
    tauj = T(1.0) ./ (aj + bj)
    dPss = ((dPss_b <= T(1)) ? (dPss_b) : (T(1)))
    fcass = fss
    fBPss = fss_P
    fcass_P = fss_P
    dhL_dt = (-hL + hLss) ./ thL
    dhL_dt_linearized = T(-1) ./ thL
    #u_new[1] = dhL_dt .* (exp(dhL_dt_linearized .* dt) - T(1)) ./ dhL_dt_linearized + hL
    u_new[63] =
        dhL_dt .* (exp(dhL_dt_linearized .* dt) - T(1)) ./ dhL_dt_linearized + hL
    jss = hss
    jss_P = hss_P
    jssp_P = hssp_P
    djca_dt = (-jca + jcass) ./ tjca
    djca_dt_linearized = T(-1) ./ tjca
    u_new[2] =
        djca_dt .* (exp(djca_dt_linearized .* dt) - T(1)) ./ djca_dt_linearized + jca
    dm_dt = (-m + mss) ./ taum
    dm_dt_linearized = T(-1) ./ taum
    u_new[3] = dm_dt .* (exp(dm_dt_linearized .* dt) - T(1)) ./ dm_dt_linearized + m
    dmL_dt = (-mL + mLss) ./ tmL
    dmL_dt_linearized = T(-1) ./ tmL
    u_new[4] = dmL_dt .* (exp(dmL_dt_linearized .* dt) - T(1)) ./ dmL_dt_linearized + mL
    xs2ss = xs1ss
    dxs1_P_dt = (-xs1_P + xs1ss) ./ txs1_P
    dxs1_P_dt_linearized = T(-1) ./ txs1_P
    u_new[5] =
        dxs1_P_dt .* (exp(dt .* dxs1_P_dt_linearized) - T(1)) ./ dxs1_P_dt_linearized +
        xs1_P
    dxs1_dt = (-xs1 + xs1ss) ./ txs1
    dxs1_dt_linearized = T(-1) ./ txs1
    u_new[6] =
        dxs1_dt .* (exp(dt .* dxs1_dt_linearized) - T(1)) ./ dxs1_dt_linearized + xs1
    f = Aff .* ff_ + Afs .* fs
    fBP = Aff .* fBPf + Afs .* fs_P
    f_P = Aff .* ff_P + Afs .* fs_P
    fp = Aff .* ffp + Afs .* fs
    Acap = T(2) * Ageo
    vjsr = T(0.0048) * vcell
    vmyo = T(0.68) * vcell
    vnsr = T(0.0552) * vcell
    vss = T(0.02) * vcell
    AiS = T(1.0) - AiF
    da_dt = (-a + ass) ./ ta
    da_dt_linearized = T(-1) ./ ta
    u_new[7] = a + da_dt .* (exp(da_dt_linearized .* dt) - T(1)) ./ da_dt_linearized
    dap_dt = (-ap + assp) ./ ta
    dap_dt_linearized = T(-1) ./ ta
    u_new[8] = ap + dap_dt .* (exp(dap_dt_linearized .* dt) - T(1)) ./ dap_dt_linearized
    CaMKa = CaMKb + CaMKt
    ICaL_taua_multiplier = ((T_val != T(310)) ? (Q10ICaL_a .^ Qpow) : (T(1)))
    ICaL_tauff_multiplier = ((T_val != T(310)) ? (Q10ICaL_ff .^ Qpow) : (T(1)))
    ICaL_taufs_multiplier = ((T_val != T(310)) ? (Q10ICaL_fs .^ Qpow) : (T(1)))
    IKb_Multiplier_1 = (
        (T_val != T(310)) ? (IKb_Multiplier .* Q10Kb .^ Qpow) : (IKb_Multiplier)
    )
    IKr_Multiplier_1 = (
        if (T_val != T(310))
            (Q10K .^ Qpow .* (IKr_Multiplier .* IKr_Multiplier))
        else
            (IKr_Multiplier)
        end
    )
    INaK_Multiplier_1 = (
        (T_val != T(310)) ? (INaK_Multiplier .* Q10NaK .^ Qpow) : (INaK_Multiplier)
    )
    IpCa_Multiplier_1 = (
        (T_val != T(310)) ? (IpCa_Multiplier .* Q10SLCaP .^ Qpow) : (IpCa_Multiplier)
    )
    IClb = (GClb .* IClb_Multiplier) .* (-ECl + v)
    INab_INab =
        ((vffrt .* (INab_Multiplier .* PNab)) .* (nai .* exp(vfrt) - nao)) ./
        (exp(vfrt) - T(1.0) / T(1))
    Knao = Knao0 .* exp((vfrt .* (T(1.0) - delta)) / T(3.0))
    alpha = T(0.1161) * exp(T(0.299) * vfrt)
    alpha_2 = T(0.0578) * exp(T(0.971) * vfrt)
    alpha_C2ToI = T(5.2e-05) * exp(T(1.525) * vfrt)
    alpha_i = T(0.2533) * exp(T(0.5953) * vfrt)
    beta_ = T(0.2442) * exp(vfrt .* (T(1.604) * (T(-1))))
    beta_2 = T(0.000349) * exp(vfrt .* (T(1.062) * (T(-1))))
    beta_i = T(0.06525) * exp(vfrt .* (T(0.8209) * (T(-1))))
    hca = exp(qca .* vfrt)
    hna = exp(qna .* vfrt)
    aK1 =
        T(4.094) ./
        (exp(T(0.1217) * ((-EK + v) - T(49.934) / T(1))) + T(1.0))
    bK1 =
        (
            T(15.72) * exp(T(0.0674) * ((-EK + v) - T(3.257) / T(1))) +
            exp(T(0.0618) * ((-EK + v) - T(594.31) / T(1)))
        ) ./ (exp((T(0.1629) * (T(-1))) * ((-EK + v) + T(14.207))) + T(1.0))
    INa_BP = jp_P .* (hp_P .* (m .^ T(3) .* (GNa_P .* (-ENa + v))))
    INa_CaMK = jp .* (hp .* (m .^ T(3) .* (GNa .* (-ENa + v))))
    INa_NP = j .* (h .* (m .^ T(3) .* (GNa .* (-ENa + v))))
    INa_PKA = j_P .* (h_P .* (m .^ T(3) .* (GNa_P .* (-ENa + v))))
    IClCa_junc = ((Fjunc .* GClCa_1) ./ (KdClCa ./ cass + T(1.0))) .* (-EClss + v)
    IClCa_sl =
        ((GClCa_1 .* (T(1.0) - Fjunc)) ./ (KdClCa ./ cai + T(1.0))) .* (-ECl + v)
    tiF = delta_epi .* tiF_b
    tiS = delta_epi .* tiS_b
    ICaL_Multiplier_2 = (
        (T_val != T(310)) ? (ICaL_Multiplier_1 .* Q10CaL .^ Qpow) : (ICaL_Multiplier_1)
    )
    IKs_Multiplier_2 = (
        if (T_val != T(310))
            (Q10K .^ Qpow .* (IKs_Multiplier .* IKs_Multiplier_1))
        else
            (IKs_Multiplier_1)
        end
    )
    GNaL = INaL_Multiplier_1 .* ((celltype_val == T(1.0)) ? (T(0.6) * GNaL_b) : (GNaL_b))
    b3 = (H .* (P .* k3m)) ./ (T(1.0) + MgATP ./ Kmgatp)
    a_relp = (T(0.5) * btp) / T(1.0)
    tau_relBP_b = (T(0.75) * btp) ./ (T(1.0) + T(0.0123) ./ cajsr)
    tau_relp_b = btp ./ (T(1.0) + T(0.0123) ./ cajsr)
    tau_rel_P = ((tau_rel_P_b < T(0.001)) ? (T(0.001)) : (tau_rel_P_b))
    tau_rel = ((tau_rel_b < T(0.001)) ? (T(0.001)) : (tau_rel_b))
    dCaMKt_dt = -CaMKt .* betaCaMKII + (CaMKb .* aCaMK) .* (CaMKb + CaMKt)
    dCaMKt_dt_linearized = CaMKb .* aCaMK - betaCaMKII
    u_new[9] =
        CaMKt + (
            if (abs(dCaMKt_dt_linearized) > T(1.0e-08))
                (
                    dCaMKt_dt .* (exp(dCaMKt_dt_linearized .* dt) - T(1)) ./
                    dCaMKt_dt_linearized
                )
            else
                (dCaMKt_dt .* dt)
            end
        )
    Bcai = T(1.0) ./ ((cmdnmax .* kmcmdn) ./ (cai + kmcmdn) .^ T(2) + T(1.0))
    gamma_cai = exp(
        (T(4.0) * (-constA)) .* (sqrt(Ii) ./ (sqrt(Ii) + T(1.0)) - T(0.3) * Ii)
    )
    gamma_cao = exp(
        (T(4.0) * (-constA)) .* (sqrt(Io) ./ (sqrt(Io) + T(1.0)) - T(0.3) * Io)
    )
    gamma_ki = exp(
        (T(1.0) * (-constA)) .* (sqrt(Ii) ./ (sqrt(Ii) + T(1.0)) - T(0.3) * Ii)
    )
    gamma_ko = exp(
        (T(1.0) * (-constA)) .* (sqrt(Io) ./ (sqrt(Io) + T(1.0)) - T(0.3) * Io)
    )
    gamma_kss = exp(
        (T(1.0) * (-constA)) .* (sqrt(Iss) ./ (sqrt(Iss) + T(1.0)) - T(0.3) * Iss)
    )
    gamma_nai = exp(
        (T(1.0) * (-constA)) .* (sqrt(Ii) ./ (sqrt(Ii) + T(1.0)) - T(0.3) * Ii)
    )
    gamma_nao = exp(
        (T(1.0) * (-constA)) .* (sqrt(Io) ./ (sqrt(Io) + T(1.0)) - T(0.3) * Io)
    )
    gamma_nass = exp(
        (T(1.0) * (-constA)) .* (sqrt(Iss) ./ (sqrt(Iss) + T(1.0)) - T(0.3) * Iss)
    )
    Knai0 = (
        if (fINaK_PKA == T(0))
            (Knai0_b)
        else
            (Knai0_P .* fINaK_PKA + Knai0_np .* (T(1) - fINaK_PKA))
        end
    )
    xb_su_gamma = XS .* gamma_rate
    xb_wu_gamma = XW .* gamma_rate_w
    h5_i = (nai .* nai) ./ ((kna2 .* (h4_i .* kna1)))
    h6_i = T(1.0) ./ h4_i
    h5_ss = (nass .* nass) ./ ((kna2 .* (h4_ss .* kna1)))
    h6_ss = T(1.0) ./ h4_ss
    h11_i = (nao .* nao) ./ ((kna2 .* (h10_i .* kna1)))
    h12_i = T(1.0) ./ h10_i
    h11_ss = (nao .* nao) ./ ((kna2 .* (h10_ss .* kna1)))
    h12_ss = T(1.0) ./ h10_ss
    Lfac_value =
        beta_0 .* (
            (lambda0 + ((lambda0 > lambda_min) ? (lambda_min) : (lambda0))) -
            (lambda_min + T(1))
        ) + T(1)
    k_ws_1 = k_ws .* mu
    ktm_block =
        (T(0.5) * (ktm_unblock .* perm50 .^ nperm)) ./ (-XWSS + (T(0.5) - XSSS))
    tmb_tmp = ((Ca_TRPN < T(1.0e-12)) ? (T(0)) : (Ca_TRPN .^ (-nperm / T(2))))
    k_uw_1 = k_uw .* nu
    IK1_Multiplier_1 = IK1_Multiplier .* pH_Multiplier_IK1
    pH_Multiplier_Ito = pH_Multiplier_IK1
    Jup_Multiplier_1 =
        pH_Multiplier_Jup .*
        ((T_val != T(310)) ? (Jup_Multiplier .* Q10SRCaP .^ Qpow) : (Jup_Multiplier))
    INaCa_Multiplier_1 =
        pH_Multiplier_NCX .*
        ((T_val != T(310)) ? (INaCa_Multiplier .* Q10NCX .^ Qpow) : (INaCa_Multiplier))
    a_rel = ((T(0.5) * bt) .* pH_Multiplier_RyR) / T(1.0)
    ca50 = pH_Multiplier_TnICa .* ((mode == T(1)) ? (T(2.5)) : (ca50_b))
    Tref = pH_Multiplier_maxTension .* ((mode == T(1)) ? (T(40.5)) : (Tref_b))
    dhLp_dt = (-hLp + hLssp) ./ thLp
    dhLp_dt_linearized = T(-1) ./ thLp
    u_new[10] =
        dhLp_dt .* (exp(dhLp_dt_linearized .* dt) - T(1)) ./ dhLp_dt_linearized + hLp
    fca = Afcaf .* fcaf + Afcas .* fcas
    fcaBP = Afcaf .* fcaBPf + Afcas .* fcas_P
    fcap = Afcaf .* fcafp + Afcas .* fcas
    fcap_P = Afcaf .* fcaf_P + Afcas .* fcas_P
    dh_P_dt = (-h_P + hss_P) ./ tauh
    dh_P_dt_linearized = T(-1) ./ tauh
    u_new[11] =
        dh_P_dt .* (exp(dh_P_dt_linearized .* dt) - T(1)) ./ dh_P_dt_linearized + h_P
    dh_dt = (-h + hss) ./ tauh
    dh_dt_linearized = T(-1) ./ tauh
    u_new[12] = dh_dt .* (exp(dh_dt_linearized .* dt) - T(1)) ./ dh_dt_linearized + h
    dhp_P_dt = (-hp_P + hssp_P) ./ tauh
    dhp_P_dt_linearized = T(-1) ./ tauh
    u_new[13] =
        dhp_P_dt .* (exp(dhp_P_dt_linearized .* dt) - T(1)) ./ dhp_P_dt_linearized + hp_P
    dhp_dt = (-hp + hssp) ./ tauh
    dhp_dt_linearized = T(-1) ./ tauh
    u_new[14] =
        dhp_dt .* (exp(dhp_dt_linearized .* dt) - T(1)) ./ dhp_dt_linearized + hp
    taujp = T(1.46) * tauj
    fcaBPss = fcass_P
    dj_dt = (-j + jss) ./ tauj
    dj_dt_linearized = T(-1) ./ tauj
    u_new[15] = dj_dt .* (exp(dj_dt_linearized .* dt) - T(1)) ./ dj_dt_linearized + j
    dj_P_dt = (-j_P + jss_P) ./ tauj
    dj_P_dt_linearized = T(-1) ./ tauj
    u_new[16] =
        dj_P_dt .* (exp(dj_P_dt_linearized .* dt) - T(1)) ./ dj_P_dt_linearized + j_P
    dxs2_dt = (-xs2 + xs2ss) ./ txs2
    dxs2_dt_linearized = T(-1) ./ txs2
    u_new[17] =
        dxs2_dt .* (exp(dt .* dxs2_dt_linearized) - T(1)) ./ dxs2_dt_linearized + xs2
    i = AiF .* iF + AiS .* iS
    ip = AiF .* iFp + AiS .* iSp
    fICaLp = T(1) ./ (T(1) * (T(1) + KmCaMK ./ CaMKa))
    fINaLp = T(1.0) ./ (T(1.0) + KmCaMK ./ CaMKa)
    fINap = T(1.0) ./ (T(1.0) + KmCaMK ./ CaMKa)
    fItop = T(1.0) ./ (T(1.0) + KmCaMK ./ CaMKa)
    fJrelp = T(1.0) ./ (T(1.0) + KmCaMK ./ CaMKa)
    fJupp = T(1.0) ./ (T(1.0) + KmCaMK ./ CaMKa)
    td =
        (
            (offset + T(0.6)) +
            T(1.0) ./ (
                exp((T(0.05) * (T(-1))) * ((v + vShift) + T(6.0))) +
                exp(T(0.09) * ((v + vShift) + T(14.0)))
            )
        ) ./ ICaL_taua_multiplier
    tfcaf =
        (
            T(7.0) +
            T(1.0) ./ (
                T(0.04) * exp((-(v - T(4.0) / T(1))) / T(7.0)) +
                T(0.04) * exp((v - T(4.0) / T(1)) / T(7.0))
            )
        ) ./ ICaL_tauff_multiplier
    tff =
        (
            T(7.0) +
            T(1.0) ./ (
                T(0.0045) * exp((-(v + T(20.0))) / T(10.0)) +
                T(0.0045) * exp((v + T(20.0)) / T(10.0))
            )
        ) ./ ICaL_tauff_multiplier
    kmnCDIincrease = ((T_val < T(310)) ? (T(1.0) ./ ICaL_taufs_multiplier) : (T(1.0)))
    tfcas =
        (
            T(100.0) +
            T(1.0) ./
            (T(0.00012) * exp((-v) / T(3.0)) + T(0.00012) * exp(v / T(7.0)))
        ) ./ ICaL_taufs_multiplier
    tfs =
        (
            T(1000.0) +
            T(1.0) ./ (
                T(3.5e-05) * exp((-(v + T(5.0))) / T(4.0)) +
                T(3.5e-05) * exp((v + T(5.0)) / T(6.0))
            )
        ) ./ ICaL_taufs_multiplier
    GKb = IKb_Multiplier_1 .* ((celltype_val == T(1.0)) ? (T(0.6) * GKb_b) : (GKb_b))
    GKr =
        IKr_Multiplier_1 .* (
            if (celltype_val == T(1.0))
                (T(1.3) * GKr_b)
            else
                (((celltype_val == T(2)) ? (T(0.8) * GKr_b) : (GKr_b)))
            end
        )
    Pnak =
        INaK_Multiplier_1 .* (
            if (celltype_val == T(1.0))
                (T(0.9) * Pnak_b)
            else
                (((celltype_val == T(2)) ? (T(0.7) * Pnak_b) : (Pnak_b)))
            end
        )
    IpCa_IpCa = (cai .* (GpCa .* IpCa_Multiplier_1)) ./ (KmCap + cai)
    a3 =
        (k3p .* (ko ./ Kko) .^ T(2)) ./ (
            ((T(1.0) + ko ./ Kko) .^ T(2) + (T(1.0) + nao ./ Knao) .^ T(3)) -
            T(1.0) / T(1)
        )
    b2 =
        (k2m .* (nao ./ Knao) .^ T(3)) ./ (
            ((T(1.0) + ko ./ Kko) .^ T(2) + (T(1.0) + nao ./ Knao) .^ T(3)) -
            T(1.0) / T(1)
        )
    dC0_dt = -C0 .* alpha + C1 .* beta_
    dC0_dt_linearized = -alpha
    u_new[18] =
        C0 + (
            if (abs(dC0_dt_linearized) > T(1.0e-08))
                (dC0_dt .* (exp(dC0_dt_linearized .* dt) - T(1)) ./ dC0_dt_linearized)
            else
                (dC0_dt .* dt)
            end
        )
    dC1_dt = (-C1) .* (alpha_1 + beta_) + (C0 .* alpha + C2 .* beta_1)
    dC1_dt_linearized = -alpha_1 - beta_
    u_new[19] =
        C1 + (
            if (abs(dC1_dt_linearized) > T(1.0e-08))
                (dC1_dt .* (exp(dC1_dt_linearized .* dt) - T(1)) ./ dC1_dt_linearized)
            else
                (dC1_dt .* dt)
            end
        )
    beta_ItoC2 = (alpha_C2ToI .* (beta_2 .* beta_i)) ./ ((alpha_2 .* alpha_i))
    dO__dt = (-O_) .* (alpha_i + beta_2) + (C2 .* alpha_2 + I_ .* beta_i)
    dO__dt_linearized = -alpha_i - beta_2
    u_new[20] =
        O_ + (
            if (abs(dO__dt_linearized) > T(1.0e-08))
                (dO__dt .* (exp(dO__dt_linearized .* dt) - T(1)) ./ dO__dt_linearized)
            else
                (dO__dt .* dt)
            end
        )
    h1_i = (nai ./ kna3) .* (hna + T(1.0)) + T(1.0)
    h1_ss = (nass ./ kna3) .* (hna + T(1.0)) + T(1.0)
    h7_i = (nao ./ kna3) .* (T(1.0) + T(1.0) ./ hna) + T(1.0)
    h7_ss = (nao ./ kna3) .* (T(1.0) + T(1.0) ./ hna) + T(1.0)
    K1ss = aK1 ./ (aK1 + bK1)
    dclss_dt = -JdiffCl + (Acap .* IClCa_junc) ./ ((F_ .* vss))
    u_new[21] = clss + dclss_dt .* dt
    IClCa = IClCa_junc + IClCa_sl
    dcli_dt = (Acap .* (IClCa_sl + IClb)) ./ ((F_ .* vmyo)) + (JdiffCl .* vss) ./ vmyo
    u_new[22] = cli + dcli_dt .* dt
    tiFp = tiF .* (dti_develop .* dti_recover)
    diF_dt = (-iF + iss) ./ tiF
    diF_dt_linearized = T(-1) ./ tiF
    u_new[23] =
        diF_dt .* (exp(diF_dt_linearized .* dt) - T(1)) ./ diF_dt_linearized + iF
    tiSp = tiS .* (dti_develop .* dti_recover)
    diS_dt = (-iS + iss) ./ tiS
    diS_dt_linearized = T(-1) ./ tiS
    u_new[24] =
        diS_dt .* (exp(diS_dt_linearized .* dt) - T(1)) ./ diS_dt_linearized + iS
    PCa = ICaL_Multiplier_2 .* (
        if (celltype_val == T(1.0))
            (T(1.2) * PCa_b)
        else
            (((celltype_val == T(2)) ? (T(2) * PCa_b) : (PCa_b)))
        end
    )
    PCa_P =
        ICaL_Multiplier_2 .* (
            if (celltype_val == T(1))
                (T(1.2) * PCa_P_b)
            else
                (((celltype_val == T(2)) ? (T(2) * PCa_P_b) : (PCa_P_b)))
            end
        )
    GKs = IKs_Multiplier_2 .* ((celltype_val == T(1.0)) ? (T(1.4) * GKs_b) : (GKs_b))
    a_relBP = T(1.4) * a_relp
    tau_relBP = ((tau_relBP_b < T(0.001)) ? (T(0.001)) : (tau_relBP_b))
    tau_relp = ((tau_relp_b < T(0.001)) ? (T(0.001)) : (tau_relp_b))
    ICab_ICab =
        (
            (vffrt .* (T(4.0) * (ICab_Multiplier .* PCab))) .*
            ((-cao) .* gamma_cao + (cai .* gamma_cai) .* exp(T(2) * vfrt))
        ) ./ (exp(T(2) * vfrt) - T(1.0) / T(1))
    PhiCaL_i =
        (
            (T(4.0) * vffrt) .*
            ((-cao) .* gamma_cao + (cai .* gamma_cai) .* exp(T(2) * vfrt))
        ) ./ (exp(T(2) * vfrt) - T(1.0) / T(1))
    PhiCaL_ss =
        (
            (T(4.0) * vffrt) .*
            (-cao .* gamma_cao + (cass .* gamma_cai) .* exp(T(2) * vfrt))
        ) ./ (exp(T(2) * vfrt) - T(1.0))
    PhiCaK_i =
        (
            (T(1.0) * vffrt) .*
            ((-gamma_ko) .* ko + (gamma_ki .* ki) .* exp(T(1.0) * vfrt))
        ) ./ (exp(T(1.0) * vfrt) - T(1.0) / T(1))
    PhiCaK_ss =
        (
            (T(1.0) * vffrt) .*
            ((-gamma_ko) .* ko + (gamma_kss .* kss) .* exp(T(1.0) * vfrt))
        ) ./ (exp(T(1.0) * vfrt) - T(1.0) / T(1))
    PhiCaNa_i =
        (
            (T(1.0) * vffrt) .*
            ((-gamma_nao) .* nao + (gamma_nai .* nai) .* exp(T(1.0) * vfrt))
        ) ./ (exp(T(1.0) * vfrt) - T(1.0) / T(1))
    PhiCaNa_ss =
        (
            (T(1.0) * vffrt) .*
            ((-gamma_nao) .* nao + (gamma_nass .* nass) .* exp(T(1.0) * vfrt))
        ) ./ (exp(T(1.0) * vfrt) - T(1.0) / T(1))
    Knai = Knai0 .* exp((delta .* vfrt) / T(3.0))
    k6_i = kcaon .* (cai .* h6_i)
    k6_ss = kcaon .* (cass .* h6_ss)
    k1_i = kcaon .* (cao .* h12_i)
    k1_ss = kcaon .* (cao .* h12_ss)
    Lfac = ((Lfac_value < T(0)) ? (T(0)) : (Lfac_value))
    cds = (wfrac .* ((k_ws_1 .* phi) .* (T(1) - dr))) ./ dr
    k_su = wfrac .* (k_ws_1 .* (T(-1) + T(1) ./ (T(1) * dr)))
    xb_ws = XW .* k_ws_1
    dTmBlocked_dt =
        -TmBlocked .* Ca_TRPN .^ (nperm / T(2)) .* ktm_unblock +
        XU .* (ktm_block .* ((tmb_tmp > T(100)) ? (T(100)) : (tmb_tmp)))
    dTmBlocked_dt_linearized = -Ca_TRPN .^ (nperm / T(2)) .* ktm_unblock
    u_new[25] =
        TmBlocked + (
            if (abs(dTmBlocked_dt_linearized) > T(1.0e-08))
                (
                    dTmBlocked_dt .* (exp(dTmBlocked_dt_linearized .* dt) - T(1)) ./
                    dTmBlocked_dt_linearized
                )
            else
                (dTmBlocked_dt .* dt)
            end
        )
    cdw =
        (((k_uw_1 .* phi) .* (T(1) - dr)) .* (T(1) - wfrac)) ./
        ((wfrac .* (T(1) - dr)))
    k_wu = k_uw_1 .* (T(-1) + T(1) ./ (T(1) * wfrac)) - k_ws_1
    xb_uw = XU .* k_uw_1
    GK1 =
        IK1_Multiplier_1 .* (
            if (celltype_val == T(1.0))
                (T(1.2) * GK1_b)
            else
                (((celltype_val == T(2)) ? (T(1.3) * GK1_b) : (GK1_b)))
            end
        )
    Ito_Multiplier_1 = Ito_Multiplier .* pH_Multiplier_Ito
    Jleak = (cansr .* (T(0.0048825) * Jup_Multiplier_1)) / T(15.0)
    Jup_BP_b =
        (cai .* (T(0.005425) * (T(2.75) * Jup_Multiplier_1))) ./
        (cai + T(0.54) * (T(0.00092) - T(0.00017)))
    Jup_P_b =
        (cai .* (T(0.005425) * Jup_Multiplier_1)) ./ (cai + T(0.00092) * T(0.54))
    Jupnp_b =
        (cai .* (T(0.005425) * (T(0.9) * Jup_Multiplier_1))) ./ (cai + T(0.00092))
    Jupp_b =
        (cai .* (T(0.005425) * (T(2.75) * (T(0.9) * Jup_Multiplier_1)))) ./
        ((cai + T(0.00092)) - T(0.00017))
    Gncx =
        INaCa_Multiplier_1 .* (
            if (celltype_val == T(1.0))
                (T(1.1) * Gncx_b)
            else
                (((celltype_val == T(2)) ? (T(1.4) * Gncx_b) : (Gncx_b)))
            end
        )
    a_relP = T(1.4) * a_rel
    ca50_1 =
        beta_1_mech .* ((lambda - T(1) < T(0.2)) ? (lambda - T(1)) : (T(0.2))) +
        ca50
    djp_P_dt = (-jp_P + jssp_P) ./ taujp
    djp_P_dt_linearized = T(-1) ./ taujp
    u_new[26] =
        djp_P_dt .* (exp(djp_P_dt_linearized .* dt) - T(1)) ./ djp_P_dt_linearized + jp_P
    djp_dt = (-jp + jss) ./ taujp
    djp_dt_linearized = T(-1) ./ taujp
    u_new[27] =
        djp_dt .* (exp(djp_dt_linearized .* dt) - T(1)) ./ djp_dt_linearized + jp
    fICaL_BP = fICaL_P .* fICaLp
    INaL_INaL = (mL .* (GNaL .* (-ENa + v))) .* (fINaLp .* hLp + hL .* (T(1.0) - fINaLp))
    fINa_BP = fINa_P .* fINap
    fJrel_BP = fJrel_PKA .* fJrelp
    fJup_BP = fJup_P .* fJupp
    dd_P_dt = (dPss - d_P) ./ td
    dd_P_dt_linearized = T(-1) ./ td
    u_new[28] =
        d_P + dd_P_dt .* (exp(dd_P_dt_linearized .* dt) - T(1)) ./ dd_P_dt_linearized
    dd_dt = (-d + dss) ./ td
    dd_dt_linearized = T(-1) ./ td
    u_new[29] = d + dd_dt .* (exp(dd_dt_linearized .* dt) - T(1)) ./ dd_dt_linearized
    tfcafp = T(2.5) * tfcaf
    dfcaf_P_dt = (-fcaf_P + fcass_P) ./ tfcaf
    dfcaf_P_dt_linearized = T(-1) ./ tfcaf
    u_new[30] =
        dfcaf_P_dt .* (exp(dfcaf_P_dt_linearized .* dt) - T(1)) ./
        dfcaf_P_dt_linearized + fcaf_P
    dfcaf_dt = (-fcaf + fcass) ./ tfcaf
    dfcaf_dt_linearized = T(-1) ./ tfcaf
    u_new[31] =
        dfcaf_dt .* (exp(dfcaf_dt_linearized .* dt) - T(1)) ./ dfcaf_dt_linearized + fcaf
    tffp = T(2.5) * tff
    dff_P_dt = (-ff_P + fss_P) ./ tff
    dff_P_dt_linearized = T(-1) ./ tff
    u_new[32] =
        dff_P_dt .* (exp(dff_P_dt_linearized .* dt) - T(1)) ./ dff_P_dt_linearized + ff_P
    dff__dt = (-ff_ + fss) ./ tff
    dff__dt_linearized = T(-1) ./ tff
    u_new[33] =
        dff__dt .* (exp(dff__dt_linearized .* dt) - T(1)) ./ dff__dt_linearized + ff_
    Kmn = Kmn_b .* kmnCDIincrease
    dfcas_P_dt = (-fcas_P + fcass_P) ./ tfcas
    dfcas_P_dt_linearized = T(-1) ./ tfcas
    u_new[34] =
        dfcas_P_dt .* (exp(dfcas_P_dt_linearized .* dt) - T(1)) ./
        dfcas_P_dt_linearized + fcas_P
    dfcas_dt = (-fcas + fcass) ./ tfcas
    dfcas_dt_linearized = T(-1) ./ tfcas
    u_new[35] =
        dfcas_dt .* (exp(dfcas_dt_linearized .* dt) - T(1)) ./ dfcas_dt_linearized + fcas
    dfs_P_dt = (-fs_P + fss_P) ./ tfs
    dfs_P_dt_linearized = T(-1) ./ tfs
    u_new[36] =
        dfs_P_dt .* (exp(dfs_P_dt_linearized .* dt) - T(1)) ./ dfs_P_dt_linearized + fs_P
    dfs_dt = (-fs + fss) ./ tfs
    dfs_dt_linearized = T(-1) ./ tfs
    u_new[37] =
        dfs_dt .* (exp(dfs_dt_linearized .* dt) - T(1)) ./ dfs_dt_linearized + fs
    GKbNP = GKb
    GKbP = T(1.2) * GKb
    IKr_IKr = (O_ .* (GKr .* (T(0.4472135954999579) * sqrt(ko)))) .* (-EK + v)
    x3 = b1 .* (a3 .* a4) + (a4 .* (b1 .* b2) + (a4 .* (a2 .* a3) + b1 .* (b2 .* b3)))
    dC2_dt =
        (-C2) .* (alpha_C2ToI + (alpha_2 + beta_1)) +
        (I_ .* beta_ItoC2 + (C1 .* alpha_1 + O_ .* beta_2))
    dC2_dt_linearized = -alpha_2 - alpha_C2ToI - beta_1
    u_new[38] =
        C2 + (
            if (abs(dC2_dt_linearized) > T(1.0e-08))
                (dC2_dt .* (exp(dC2_dt_linearized .* dt) - T(1)) ./ dC2_dt_linearized)
            else
                (dC2_dt .* dt)
            end
        )
    dI__dt = (-I_) .* (beta_ItoC2 + beta_i) + (C2 .* alpha_C2ToI + O_ .* alpha_i)
    dI__dt_linearized = -beta_ItoC2 - beta_i
    u_new[39] =
        I_ + (
            if (abs(dI__dt_linearized) > T(1.0e-08))
                (dI__dt .* (exp(dI__dt_linearized .* dt) - T(1)) ./ dI__dt_linearized)
            else
                (dI__dt .* dt)
            end
        )
    h2_i = (hna .* nai) ./ ((h1_i .* kna3))
    h3_i = T(1.0) ./ h1_i
    h2_ss = (hna .* nass) ./ ((h1_ss .* kna3))
    h3_ss = T(1.0) ./ h1_ss
    h8_i = nao ./ ((h7_i .* (hna .* kna3)))
    h9_i = T(1.0) ./ h7_i
    h8_ss = nao ./ ((h7_ss .* (hna .* kna3)))
    h9_ss = T(1.0) ./ h7_ss
    diFp_dt = (-iFp + iss) ./ tiFp
    diFp_dt_linearized = T(-1) ./ tiFp
    u_new[40] =
        diFp_dt .* (exp(diFp_dt_linearized .* dt) - T(1)) ./ diFp_dt_linearized + iFp
    diSp_dt = (-iSp + iss) ./ tiSp
    diSp_dt_linearized = T(-1) ./ tiSp
    u_new[41] =
        diSp_dt .* (exp(diSp_dt_linearized .* dt) - T(1)) ./ diSp_dt_linearized + iSp
    PCaK = T(0.0003574) * PCa
    PCaNa = T(0.00125) * PCa
    PCap = T(1.1) * PCa
    PCaK_P = T(0.0003574) * PCa_P
    PCaNa_P = T(0.00125) * PCa_P
    GKs_P = T(50) * GKs
    IKs_NP = (xs2 .* (xs1 .* (GKs .* KsCa))) .* (-EKs + v)
    ICaL_i_BP =
        (d_P .* (PCa_P .* PhiCaL_i)) .*
        (fBP .* (T(1.0) - nca_i) + nca_i .* (fcaBP .* jca))
    ICaL_i_NP =
        (d .* (PCa .* PhiCaL_i)) .* (f .* (T(1.0) - nca_i) + nca_i .* (fca .* jca))
    ICaL_i_PKA =
        (d_P .* (PCa_P .* PhiCaL_i)) .*
        (f_P .* (T(1.0) - nca_i) + nca_i .* (fcap_P .* jca))
    ICaL_ss_BP =
        (d_P .* (PCa_P .* PhiCaL_ss)) .* (fBP .* (T(1.0) - nca) + nca .* (fcaBP .* jca))
    ICaL_ss_NP = (d .* (PCa .* PhiCaL_ss)) .* (f .* (T(1.0) - nca) + nca .* (fca .* jca))
    ICaL_ss_PKA =
        (d_P .* (PCa_P .* PhiCaL_ss)) .* (f_P .* (T(1.0) - nca) + nca .* (fcap_P .* jca))
    a1 =
        (k1p .* (nai ./ Knai) .^ T(3)) ./ (
            ((T(1.0) + ki ./ Kki) .^ T(2) + (T(1.0) + nai ./ Knai) .^ T(3)) -
            T(1.0) / T(1)
        )
    b4 =
        (k4m .* (ki ./ Kki) .^ T(2)) ./ (
            ((T(1.0) + ki ./ Kki) .^ T(2) + (T(1.0) + nai ./ Knai) .^ T(3)) -
            T(1.0) / T(1)
        )
    Ta = (Lfac .* (Tref ./ dr)) .* (XS .* (Zetas + T(1)) + XW .* Zetaw)
    dZetas_dt = A .* lambda_rate - Zetas .* cds
    dZetas_dt_linearized = -cds
    u_new[42] =
        Zetas + (
            if (abs(dZetas_dt_linearized) > T(1.0e-08))
                (
                    dZetas_dt .* (exp(dZetas_dt_linearized .* dt) - T(1)) ./
                    dZetas_dt_linearized
                )
            else
                (dZetas_dt .* dt)
            end
        )
    xb_su = XS .* k_su
    dZetaw_dt = A .* lambda_rate - Zetaw .* cdw
    dZetaw_dt_linearized = -cdw
    u_new[43] =
        Zetaw + (
            if (abs(dZetaw_dt_linearized) > T(1.0e-08))
                (
                    dZetaw_dt .* (exp(dZetaw_dt_linearized .* dt) - T(1)) ./
                    dZetaw_dt_linearized
                )
            else
                (dZetaw_dt .* dt)
            end
        )
    xb_wu = XW .* k_wu
    IK1_IK1 = (K1ss .* (GK1 .* (T(0.4472135954999579) * sqrt(ko)))) .* (-EK + v)
    Gto =
        Ito_Multiplier_1 .* ((celltype_val in (T(1), T(2))) ? (T(2) * Gto_b) : (Gto_b))
    Jup_BP = ((celltype_val == T(1)) ? (T(1.3) * Jup_BP_b) : (Jup_BP_b))
    Jup_P = ((celltype_val == T(1)) ? (T(1.3) * Jup_P_b) : (Jup_P_b))
    Jupnp = ((celltype_val == T(1)) ? (T(1.3) * Jupnp_b) : (Jupnp_b))
    Jupp = ((celltype_val == T(1)) ? (T(1.3) * Jupp_b) : (Jupp_b))
    dCa_TRPN_dt =
        koff .* (-Ca_TRPN + ((T(1000) * cai) ./ ca50_1) .^ TRPN_n .* (T(1) - Ca_TRPN))
    dCa_TRPN_dt_linearized = koff .* (-((T(1000) * cai) ./ ca50_1) .^ TRPN_n - T(1))
    u_new[44] =
        Ca_TRPN + (
            if (abs(dCa_TRPN_dt_linearized) > T(1.0e-08))
                (
                    dCa_TRPN_dt .* (exp(dCa_TRPN_dt_linearized .* dt) - T(1)) ./
                    dCa_TRPN_dt_linearized
                )
            else
                (dCa_TRPN_dt .* dt)
            end
        )
    fICaL_CaMKonly = -fICaL_BP + fICaLp
    fICaL_PKAonly = -fICaL_BP + fICaL_P
    fINa_CaMKonly = -fINa_BP + fINap
    fINa_PKAonly = -fINa_BP + fINa_P
    fJrel_CaMKonly = -fJrel_BP + fJrelp
    fJrel_PKAonly = -fJrel_BP + fJrel_PKA
    fJup_CaMKonly = -fJup_BP + fJupp
    fJup_PKAonly = -fJup_BP + fJup_P
    dfcaBPf_dt = (-fcaBPf + fcaBPss) ./ tfcafp
    dfcaBPf_dt_linearized = T(-1) ./ tfcafp
    u_new[45] =
        dfcaBPf_dt .* (exp(dfcaBPf_dt_linearized .* dt) - T(1)) ./
        dfcaBPf_dt_linearized + fcaBPf
    dfcafp_dt = (-fcafp + fcass) ./ tfcafp
    dfcafp_dt_linearized = T(-1) ./ tfcafp
    u_new[46] =
        dfcafp_dt .* (exp(dfcafp_dt_linearized .* dt) - T(1)) ./ dfcafp_dt_linearized +
        fcafp
    dfBPf_dt = (-fBPf + fBPss) ./ tffp
    dfBPf_dt_linearized = T(-1) ./ tffp
    u_new[47] =
        dfBPf_dt .* (exp(dfBPf_dt_linearized .* dt) - T(1)) ./ dfBPf_dt_linearized + fBPf
    dffp_dt = (-ffp + fss) ./ tffp
    dffp_dt_linearized = T(-1) ./ tffp
    u_new[48] =
        dffp_dt .* (exp(dffp_dt_linearized .* dt) - T(1)) ./ dffp_dt_linearized + ffp
    anca_i = T(1.0) ./ (k2n ./ km2n + (Kmn ./ cai + T(1.0)) .^ T(4))
    anca_ss = T(1.0) ./ (k2n ./ km2n + (Kmn ./ cass + T(1.0)) .^ T(4))
    IKb_NP = (GKbNP .* xkb) .* (-EK + v)
    IKb_P = (GKbP .* xkb) .* (-EK + v)
    k4pp_i = h2_i .* wnaca
    k7_i = wna .* (h2_i .* h5_i)
    k4p_i = (h3_i .* wca) ./ hca
    k4pp_ss = h2_ss .* wnaca
    k7_ss = wna .* (h2_ss .* h5_ss)
    k4p_ss = (h3_ss .* wca) ./ hca
    k3pp_i = h8_i .* wnaca
    k8_i = wna .* (h11_i .* h8_i)
    k3p_i = h9_i .* wca
    k3pp_ss = h8_ss .* wnaca
    k8_ss = wna .* (h11_ss .* h8_ss)
    k3p_ss = h9_ss .* wca
    ICaK_i_NP =
        (d .* (PCaK .* PhiCaK_i)) .* (f .* (T(1.0) - nca_i) + nca_i .* (fca .* jca))
    ICaK_ss_NP =
        (d .* (PCaK .* PhiCaK_ss)) .* (f .* (T(1.0) - nca) + nca .* (fca .* jca))
    ICaNa_i_NP =
        (d .* (PCaNa .* PhiCaNa_i)) .* (f .* (T(1.0) - nca_i) + nca_i .* (fca .* jca))
    ICaNa_ss_NP =
        (d .* (PCaNa .* PhiCaNa_ss)) .* (f .* (T(1.0) - nca) + nca .* (fca .* jca))
    ICaL_i_CaMK =
        (d .* (PCap .* PhiCaL_i)) .* (fp .* (T(1.0) - nca_i) + nca_i .* (fcap .* jca))
    ICaL_ss_CaMK =
        (d .* (PCap .* PhiCaL_ss)) .* (fp .* (T(1.0) - nca) + nca .* (fcap .* jca))
    PCaKp = T(0.0003574) * PCap
    PCaNap = T(0.00125) * PCap
    ICaK_i_BP =
        (d_P .* (PCaK_P .* PhiCaK_i)) .*
        (fBP .* (T(1.0) - nca_i) + nca_i .* (fcaBP .* jca))
    ICaK_i_PKA =
        (d_P .* (PCaK_P .* PhiCaK_i)) .*
        (f_P .* (T(1.0) - nca_i) + nca_i .* (fcap_P .* jca))
    ICaK_ss_BP =
        (d_P .* (PCaK_P .* PhiCaK_ss)) .* (fBP .* (T(1.0) - nca) + nca .* (fcaBP .* jca))
    ICaK_ss_PKA =
        (d_P .* (PCaK_P .* PhiCaK_ss)) .*
        (f_P .* (T(1.0) - nca) + nca .* (fcap_P .* jca))
    ICaNa_i_BP =
        (d_P .* (PCaNa_P .* PhiCaNa_i)) .*
        (fBP .* (T(1.0) - nca_i) + nca_i .* (fcaBP .* jca))
    ICaNa_i_PKA =
        (d_P .* (PCaNa_P .* PhiCaNa_i)) .*
        (f_P .* (T(1.0) - nca_i) + nca_i .* (fcap_P .* jca))
    ICaNa_ss_BP =
        (d_P .* (PCaNa_P .* PhiCaNa_ss)) .*
        (fBP .* (T(1.0) - nca) + nca .* (fcaBP .* jca))
    ICaNa_ss_PKA =
        (d_P .* (PCaNa_P .* PhiCaNa_ss)) .*
        (f_P .* (T(1.0) - nca) + nca .* (fcap_P .* jca))
    IKs_P = (xs2 .* (xs1_P .* (GKs_P .* KsCa))) .* (-EKs + v)
    x1 = a2 .* (a1 .* b3) + (b3 .* (a2 .* b4) + (a2 .* (a1 .* a4) + b3 .* (b2 .* b4)))
    x2 = b4 .* (a2 .* a3) + (b4 .* (a3 .* b1) + (a3 .* (a1 .* a2) + b4 .* (b1 .* b2)))
    x4 = a1 .* (b2 .* b3) + (a1 .* (a4 .* b2) + (a1 .* (a3 .* a4) + b2 .* (b3 .* b4)))
    dXS_dt = -xb_su_gamma + (-xb_su + xb_ws)
    u_new[49] = XS + dXS_dt .* dt
    dXW_dt = -xb_wu_gamma + (-xb_ws + (xb_uw - xb_wu))
    u_new[50] = XW + dXW_dt .* dt
    Ito_Ito = (Gto .* (-EK + v)) .* (i .* (a .* (T(1.0) - fItop)) + ip .* (ap .* fItop))
    INa_INa =
        INa_Multiplier_1 .* (
            INa_BP .* fINa_BP + (
                INa_PKA .* fINa_PKAonly + (
                    INa_CaMK .* fINa_CaMKonly +
                    INa_NP .* (-fINa_BP + (-fINa_PKAonly + (T(1) - fINa_CaMKonly)))
                )
            )
        )
    Jrel =
        (Jrel_Multiplier .* Jrel_b) .* (
            Jrel_p_P .* fJrel_BP + (
                Jrel_np_P .* fJrel_PKAonly + (
                    Jrel_np .*
                    (-fJrel_BP + (-fJrel_PKAonly + (T(1.0) - fJrel_CaMKonly))) +
                    Jrel_p .* fJrel_CaMKonly
                )
            )
        )
    Jup =
        -Jleak + (
            Jup_BP .* fJup_BP + (
                Jup_P .* fJup_PKAonly + (
                    Jupnp .* (-fJup_BP + (-fJup_PKAonly + (T(1.0) - fJup_CaMKonly))) +
                    Jupp .* fJup_CaMKonly
                )
            )
        )
    dnca_i_dt = anca_i .* k2n - km2n .* nca_i
    dnca_i_dt_linearized = -km2n
    u_new[51] =
        nca_i + (
            if (abs(dnca_i_dt_linearized) > T(1.0e-08))
                (
                    dnca_i_dt .* (exp(dnca_i_dt_linearized .* dt) - T(1)) ./
                    dnca_i_dt_linearized
                )
            else
                (dnca_i_dt .* dt)
            end
        )
    dnca_ss_dt = anca_ss .* k2n - km2n .* nca_ss
    dnca_ss_dt_linearized = -km2n
    u_new[52] =
        nca_ss + (
            if (abs(dnca_ss_dt_linearized) > T(1.0e-08))
                (
                    dnca_ss_dt .* (exp(dnca_ss_dt_linearized .* dt) - T(1)) ./
                    dnca_ss_dt_linearized
                )
            else
                (dnca_ss_dt .* dt)
            end
        )
    IKb_IKb = IKb_NP .* (T(1) - fIKurP) + IKb_P .* fIKurP
    k4_i = k4p_i + k4pp_i
    k4_ss = k4p_ss + k4pp_ss
    k3_i = k3p_i + k3pp_i
    k3_ss = k3p_ss + k3pp_ss
    ICaL_i_b =
        ICaL_i_BP .* fICaL_BP + (
            ICaL_i_PKA .* fICaL_PKAonly + (
                ICaL_i_CaMK .* fICaL_CaMKonly +
                ICaL_i_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    ICaL_ss_b =
        ICaL_ss_BP .* fICaL_BP + (
            ICaL_ss_PKA .* fICaL_PKAonly + (
                ICaL_ss_CaMK .* fICaL_CaMKonly +
                ICaL_ss_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    ICaK_i_CaMK =
        (d .* (PCaKp .* PhiCaK_i)) .* (fp .* (T(1.0) - nca_i) + nca_i .* (fcap .* jca))
    ICaK_ss_CaMK =
        (d .* (PCaKp .* PhiCaK_ss)) .* (fp .* (T(1.0) - nca) + nca .* (fcap .* jca))
    ICaNa_i_CaMK =
        (d .* (PCaNap .* PhiCaNa_i)) .* (fp .* (T(1.0) - nca_i) + nca_i .* (fcap .* jca))
    ICaNa_ss_CaMK =
        (d .* (PCaNap .* PhiCaNa_ss)) .* (fp .* (T(1.0) - nca) + nca .* (fcap .* jca))
    IKs_IKs = IKs_NP .* (T(1) - fIKsP) + IKs_P .* fIKsP
    E1_ = x1 ./ (x4 + (x3 + (x1 + x2)))
    E2 = x2 ./ (x4 + (x3 + (x1 + x2)))
    E3 = x3 ./ (x4 + (x3 + (x1 + x2)))
    E4 = x4 ./ (x4 + (x3 + (x1 + x2)))
    dcajsr_dt = Bcajsr .* (-Jrel + Jtr)
    u_new[53] = cajsr + dcajsr_dt .* dt
    dcansr_dt = Jup - Jtr .* vjsr ./ vnsr
    u_new[54] = cansr + dcansr_dt .* dt
    x2_i = (k1_i .* k7_i) .* (k4_i + k5_i) + (k4_i .* k6_i) .* (k1_i + k8_i)
    x2_ss = (k1_ss .* k7_ss) .* (k4_ss + k5_ss) + (k4_ss .* k6_ss) .* (k1_ss + k8_ss)
    x1_i = (k2_i .* k4_i) .* (k6_i + k7_i) + (k5_i .* k7_i) .* (k2_i + k3_i)
    x3_i = (k1_i .* k3_i) .* (k6_i + k7_i) + (k6_i .* k8_i) .* (k2_i + k3_i)
    x4_i = (k2_i .* k8_i) .* (k4_i + k5_i) + (k3_i .* k5_i) .* (k1_i + k8_i)
    x1_ss = (k2_ss .* k4_ss) .* (k6_ss + k7_ss) + (k5_ss .* k7_ss) .* (k2_ss + k3_ss)
    x3_ss = (k1_ss .* k3_ss) .* (k6_ss + k7_ss) + (k6_ss .* k8_ss) .* (k2_ss + k3_ss)
    x4_ss = (k2_ss .* k8_ss) .* (k4_ss + k5_ss) + (k3_ss .* k5_ss) .* (k1_ss + k8_ss)
    ICaL_i = ICaL_i_b .* (T(1) - ICaL_fractionSS)
    ICaL_ss = ICaL_fractionSS .* ICaL_ss_b
    ICaK_i_b =
        ICaK_i_BP .* fICaL_BP + (
            ICaK_i_PKA .* fICaL_PKAonly + (
                ICaK_i_CaMK .* fICaL_CaMKonly +
                ICaK_i_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    ICaK_ss_b =
        ICaK_ss_BP .* fICaL_BP + (
            ICaK_ss_PKA .* fICaL_PKAonly + (
                ICaK_ss_CaMK .* fICaL_CaMKonly +
                ICaK_ss_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    ICaNa_i_b =
        ICaNa_i_BP .* fICaL_BP + (
            ICaNa_i_PKA .* fICaL_PKAonly + (
                ICaNa_i_CaMK .* fICaL_CaMKonly +
                ICaNa_i_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    ICaNa_ss_b =
        ICaNa_ss_BP .* fICaL_BP + (
            ICaNa_ss_PKA .* fICaL_PKAonly + (
                ICaNa_ss_CaMK .* fICaL_CaMKonly +
                ICaNa_ss_NP .* (-fICaL_BP + (-fICaL_PKAonly + (T(1) - fICaL_CaMKonly)))
            )
        )
    JnakNa = T(3.0) * (E1_ .* a3 - E2 .* b3)
    JnakK = T(2) * ((-E3) .* a1 + E4 .* b1)
    E1_i = x1_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E2_i = x2_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E3_i = x3_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E4_i = x4_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E1_ss = x1_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E2_ss = x2_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E3_ss = x3_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E4_ss = x4_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    ICaL_ICaL = ICaL_i + ICaL_ss
    Jrel_inf_b =
        ((ICaL_ss .* (-a_rel)) / T(1.0)) ./
        ((cajsr_half ./ cajsr) .^ T(8.0) + T(1.0))
    Jrel_infp_b =
        ((ICaL_ss .* (-a_relp)) / T(1.0)) ./
        ((cajsr_half ./ cajsr) .^ T(8.0) + T(1.0))
    ICaK_i = ICaK_i_b .* (T(1) - ICaL_fractionSS)
    ICaK_ss = ICaK_ss_b .* ICaL_fractionSS
    ICaNa_i = ICaNa_i_b .* (T(1) - ICaL_fractionSS)
    ICaNa_ss = ICaL_fractionSS .* ICaNa_ss_b
    INaK_INaK = Pnak .* (JnakK .* zk + JnakNa .* zna)
    JncxCa_i = (-E1_i) .* k1_i + E2_i .* k2_i
    JncxNa_i =
        (-E2_i) .* k3pp_i + (E3_i .* k4pp_i + T(3.0) * ((-E1_i) .* k8_i + E4_i .* k7_i))
    JncxCa_ss = (-E1_ss) .* k1_ss + E2_ss .* k2_ss
    JncxNa_ss =
        (-E2_ss) .* k3pp_ss +
        (E3_ss .* k4pp_ss + T(3.0) * ((-E1_ss) .* k8_ss + E4_ss .* k7_ss))
    Jrel_infBP_b =
        ((-ICaL_ICaL) .* a_relBP) ./ ((jsrMidpoint ./ cajsr) .^ T(8.0) + T(1.0))
    Jrel_inf_P_b =
        ((-ICaL_ICaL) .* a_relP) ./ ((jsrMidpoint ./ cajsr) .^ T(8.0) + T(1.0))
    Jrel_inf = ((celltype_val == T(2)) ? (T(1.7) * Jrel_inf_b) : (Jrel_inf_b))
    Jrel_infp = ((celltype_val == T(2)) ? (T(1.7) * Jrel_infp_b) : (Jrel_infp_b))
    ICaK = ICaK_i + ICaK_ss
    dkss_dt = -JdiffK + (Acap .* (-ICaK_ss)) ./ ((F_ .* vss))
    u_new[55] = dkss_dt .* dt + kss
    ICaNa = ICaNa_i + ICaNa_ss
    dki_dt =
        (
            Acap .* (
                -(
                    ICaK_i + (
                        (T(-2)) * INaK_INaK +
                        (Istim + (IKb_IKb + (IK1_IK1 + (IKs_IKs + (IKr_IKr + Ito_Ito)))))
                    )
                )
            )
        ) ./ ((F_ .* vmyo)) + (JdiffK .* vss) ./ vmyo
    u_new[56] = dki_dt .* dt + ki
    INaCa_i =
        (allo_i .* (Gncx .* (T(1.0) - INaCa_fractionSS))) .*
        (JncxCa_i .* zca + JncxNa_i .* zna)
    INaCa_ss =
        (allo_ss .* (Gncx .* INaCa_fractionSS)) .* (JncxCa_ss .* zca + JncxNa_ss .* zna)
    Jrel_infBP = ((celltype_val == T(2)) ? (T(1.7) * Jrel_infBP_b) : (Jrel_infBP_b))
    Jrel_inf_P = ((celltype_val == T(2)) ? (T(1.7) * Jrel_inf_P_b) : (Jrel_inf_P_b))
    dJrel_np_dt = (Jrel_inf - Jrel_np) ./ tau_rel
    dJrel_np_dt_linearized = T(-1) ./ tau_rel
    u_new[57] =
        Jrel_np +
        dJrel_np_dt .* (exp(dJrel_np_dt_linearized .* dt) - T(1)) ./
        dJrel_np_dt_linearized
    dJrel_p_dt = (Jrel_infp - Jrel_p) ./ tau_relp
    dJrel_p_dt_linearized = T(-1) ./ tau_relp
    u_new[58] =
        Jrel_p +
        dJrel_p_dt .* (exp(dJrel_p_dt_linearized .* dt) - T(1)) ./ dJrel_p_dt_linearized
    dcai_dt =
        Bcai .* (
            -dCa_TRPN_dt .* trpnmax + (
                (
                    (Acap .* (-(T(-2) * INaCa_i + (ICab_ICab + (ICaL_i + IpCa_IpCa))))) ./
                    (((T(2) * F_) .* vmyo)) - Jup .* vnsr ./ vmyo
                ) + (Jdiff .* vss) ./ vmyo
            )
        )
    u_new[59] = cai + dcai_dt .* dt
    dnai_dt =
        (
            Acap .* (
                -(
                    INab_INab + (
                        T(3.0) * INaK_INaK +
                        (ICaNa_i + (T(3.0) * INaCa_i + (INaL_INaL + INa_INa)))
                    )
                )
            )
        ) ./ ((F_ .* vmyo)) + (JdiffNa .* vss) ./ vmyo
    u_new[60] = dnai_dt .* dt + nai
    dcass_dt =
        Bcass .* (
            -Jdiff + (
                (Acap .* (-(ICaL_ss - T(2) * INaCa_ss))) ./ (((T(2) * F_) .* vss)) +
                (Jrel .* vjsr) ./ vss
            )
        )
    u_new[61] = cass + dcass_dt .* dt
    dnass_dt = -JdiffNa + (Acap .* (-(ICaNa_ss + T(3.0) * INaCa_ss))) ./ ((F_ .* vss))
    u_new[62] = dnass_dt .* dt + nass
    dv_dt = -(
        INa_INa +
        INaL_INaL +
        Ito_Ito +
        ICaL_ICaL +
        ICaNa +
        ICaK +
        IKr_IKr +
        IKs_IKs +
        IK1_IK1 +
        INaCa_i +
        INaCa_ss +
        INaK_INaK +
        INab_INab +
        IKb_IKb +
        IpCa_IpCa +
        ICab_ICab +
        IClCa +
        IClb +
        Istim
    )
    #u_new[63] = dt .* dv_dt + v
    u_new[1] = dt .* dv_dt + v
    dJrel_p_P_dt = (Jrel_infBP - Jrel_p_P) ./ tau_relBP
    dJrel_p_P_dt_linearized = T(-1) ./ tau_relBP
    u_new[64] =
        Jrel_p_P +
        dJrel_p_P_dt .* (exp(dJrel_p_P_dt_linearized .* dt) - T(1)) ./
        dJrel_p_P_dt_linearized
    dJrel_np_P_dt = (Jrel_inf_P - Jrel_np_P) ./ tau_rel_P
    dJrel_np_P_dt_linearized = T(-1) ./ tau_rel_P
    u_new[65] =
        Jrel_np_P +
        dJrel_np_P_dt .* (exp(dJrel_np_P_dt_linearized .* dt) - T(1)) ./
        dJrel_np_P_dt_linearized
    return nothing
end
