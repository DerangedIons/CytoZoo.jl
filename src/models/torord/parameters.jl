const TORORD_NUM_PARAMS = 181

const TORORD_PARAMETER_NAMES = (
    :Aff,               # 1
    :BSLmax,            # 2
    :BSRmax,            # 3
    :CaMKo,             # 4
    :EKshift,           # 5
    :F,                 # 6
    :Fjunc,             # 7
    :GClCa,             # 8
    :GClb,              # 9
    :GK1_b,             # 10
    :GKb_b,             # 11
    :GKr_b,             # 12
    :GKs_b,             # 13
    :GNa,               # 14
    :GNaL_b,            # 15
    :Gncx_b,            # 16
    :GpCa,              # 17
    :Gto_b,             # 18
    :H,                 # 19
    :ICaL_Multiplier,   # 20
    :ICaL_fractionSS,   # 21
    :ICab_Multiplier,   # 22
    :IClCa_Multiplier,  # 23
    :IClb_Multiplier,   # 24
    :IK1_Multiplier,    # 25
    :IKb_Multiplier,    # 26
    :IKr_Multiplier,    # 27
    :IKs_Multiplier,    # 28
    :INaCa_Multiplier,  # 29
    :INaCa_fractionSS,  # 30
    :INaK_Multiplier,   # 31
    :INaL_Multiplier,   # 32
    :INa_Multiplier,    # 33
    :INab_Multiplier,   # 34
    :IpCa_Multiplier,   # 35
    :Ito_Multiplier,    # 36
    :Jrel_Multiplier,   # 37
    :Jrel_b,            # 38
    :Jup_Multiplier,    # 39
    :KdClCa,            # 40
    :Khp,               # 41
    :Kki,               # 42
    :Kko,               # 43
    :KmBSL,             # 44
    :KmBSR,             # 45
    :KmCaAct,           # 46
    :KmCaM,             # 47
    :KmCaMK,            # 48
    :KmCap,             # 49
    :Kmgatp,            # 50
    :Kmn_b,             # 51
    :Knai0_b,           # 52
    :Knai0_np,          # 53
    :Knao0,             # 54
    :Knap,              # 55
    :Kxkur,             # 56
    :L,                 # 57
    :MgADP,             # 58
    :MgATP,             # 59
    :PCa_P_b,           # 60
    :PCa_b,             # 61
    :PCab,              # 62
    :PKNa,              # 63
    :PNab,              # 64
    :Pnak_b,            # 65
    :Q10CaL,            # 66
    :Q10ICaL_a,         # 67
    :Q10ICaL_ff,        # 68
    :Q10ICaL_fs,        # 69
    :Q10K,              # 70
    :Q10Kb,             # 71
    :Q10NCX,            # 72
    :Q10NaK,            # 73
    :Q10SLCaP,          # 74
    :Q10SRCaP,          # 75
    :R,                 # 76
    :T,                 # 77
    :TOT_A,             # 78
    :TRPN_n,            # 79
    :Tref_b,            # 80
    :Whole_cell_PP1,    # 81
    :aCaMK,             # 82
    :alpha_1,           # 83
    :bCaMK,             # 84
    :beta_0,            # 85
    :beta_1,            # 86
    :beta_1_mech,       # 87
    :bt,                # 88
    :ca50_b,            # 89
    :cajsr_half,        # 90
    :cao,               # 91
    :celltype,          # 92
    :clo,               # 93
    :cmdnmax_b,         # 94
    :csqnmax,           # 95
    :delta,             # 96
    :dielConstant,      # 97
    :dr,                # 98
    :eP,                # 99
    :fICaLP,            # 100
    :fIKsP,             # 101
    :fIKurP,            # 102
    :fINaKP,            # 103
    :fINaP,             # 104
    :fPLBP,             # 105
    :fRyRP,             # 106
    :fTnIP,             # 107
    :gamma,             # 108
    :gamma_wu,          # 109
    :i_Stim_Amplitude,  # 110
    :i_Stim_Period,     # 111
    :i_Stim_PulseDuration, # 112
    :i_Stim_Start,      # 113
    :isHypoxic,         # 114
    :jsrMidpoint,       # 115
    :k1m,               # 116
    :k1p,               # 117
    :k2m,               # 118
    :k2n,               # 119
    :k2p,               # 120
    :k3m,               # 121
    :k3p,               # 122
    :k4m,               # 123
    :k4p,               # 124
    :k_uw,              # 125
    :k_ws,              # 126
    :kasymm,            # 127
    :kcaoff,            # 128
    :kcaon,             # 129
    :kmcmdn,            # 130
    :kmcsqn,            # 131
    :kmtrpn_b,          # 132
    :kna1,              # 133
    :kna2,              # 134
    :kna3,              # 135
    :ko,                # 136
    :koff,              # 137
    :ktm_unblock,       # 138
    :lambda,            # 139
    :lambda_max,        # 140
    :lambda_min,        # 141
    :lambda_rate,       # 142
    :mode,              # 143
    :mu_b,              # 144
    :nK1,               # 145
    :nNaCa,             # 146
    :nTnI,              # 147
    :nao,               # 148
    :nperm_b,           # 149
    :nrel,              # 150
    :nu_b,              # 151
    :nup,               # 152
    :offset,            # 153
    :pH,                # 154
    :perm50,            # 155
    :ph_bt,             # 156
    :phi,               # 157
    :pkK1,              # 158
    :pkNaCa,            # 159
    :pkTnI,             # 160
    :pkrel,             # 161
    :pkup,              # 162
    :qca,               # 163
    :qna,               # 164
    :rad_,              # 165
    :tauCa,             # 166
    :tauK,              # 167
    :tauNa,             # 168
    :thL,               # 169
    :tjca,              # 170
    :trpnmax,           # 171
    :vShift,            # 172
    :wca,               # 173
    :wfrac,             # 174
    :wna,               # 175
    :wnaca,             # 176
    :zca,               # 177
    :zcl,               # 178
    :zk,                # 179
    :zna,               # 180
    :IKr_Multiplier_hetero, # 181
)

const TORORD_PARAM_INDEX = Dict{Symbol,Int}(
    name => i for (i, name) in enumerate(TORORD_PARAMETER_NAMES)
)

function _torord_init_parameters!(p)
    p[1] = 0.6                # Aff
    p[2] = 1.124              # BSLmax
    p[3] = 0.047              # BSRmax
    p[4] = 0.05               # CaMKo
    p[5] = 0.0                # EKshift
    p[6] = 96485.0            # F
    p[7] = 1.0                # Fjunc
    p[8] = 0.2843             # GClCa
    p[9] = 0.00198            # GClb
    p[10] = 0.6992            # GK1_b
    p[11] = 0.0189            # GKb_b
    p[12] = 0.0321            # GKr_b
    p[13] = 0.0011            # GKs_b
    p[14] = 11.7802           # GNa
    p[15] = 0.0279            # GNaL_b
    p[16] = 0.0034            # Gncx_b
    p[17] = 0.0005            # GpCa
    p[18] = 0.16              # Gto_b
    p[19] = 1.0e-07           # H
    p[20] = 1.0               # ICaL_Multiplier
    p[21] = 0.8               # ICaL_fractionSS
    p[22] = 1.0               # ICab_Multiplier
    p[23] = 1.0               # IClCa_Multiplier
    p[24] = 1.0               # IClb_Multiplier
    p[25] = 1.0               # IK1_Multiplier
    p[26] = 1.0               # IKb_Multiplier
    p[27] = 1.0               # IKr_Multiplier
    p[28] = 1.0               # IKs_Multiplier
    p[29] = 1.0               # INaCa_Multiplier
    p[30] = 0.35              # INaCa_fractionSS
    p[31] = 1.0               # INaK_Multiplier
    p[32] = 1.0               # INaL_Multiplier
    p[33] = 1.0               # INa_Multiplier
    p[34] = 1.0               # INab_Multiplier
    p[35] = 1.0               # IpCa_Multiplier
    p[36] = 1.0               # Ito_Multiplier
    p[37] = 1.0               # Jrel_Multiplier
    p[38] = 1.5378            # Jrel_b
    p[39] = 1.0               # Jup_Multiplier
    p[40] = 0.1               # KdClCa
    p[41] = 1.698e-07         # Khp
    p[42] = 0.5               # Kki
    p[43] = 0.3582            # Kko
    p[44] = 0.0087            # KmBSL
    p[45] = 0.00087           # KmBSR
    p[46] = 0.00015           # KmCaAct
    p[47] = 0.0015            # KmCaM
    p[48] = 0.15              # KmCaMK
    p[49] = 0.0005            # KmCap
    p[50] = 1.698e-07         # Kmgatp
    p[51] = 0.002             # Kmn_b
    p[52] = 9.073             # Knai0_b
    p[53] = 9.073             # Knai0_np
    p[54] = 27.78             # Knao0
    p[55] = 224.0             # Knap
    p[56] = 292.0             # Kxkur
    p[57] = 0.01              # L
    p[58] = 0.05              # MgADP
    p[59] = 9.8               # MgATP
    p[60] = 8.3757e-05 * 2.0  # PCa_P_b
    p[61] = 8.3757e-05        # PCa_b
    p[62] = 5.9194e-08        # PCab
    p[63] = 0.01833           # PKNa
    p[64] = 1.9239e-09        # PNab
    p[65] = 15.4509           # Pnak_b
    p[66] = 1.8               # Q10CaL
    p[67] = 2.5               # Q10ICaL_a
    p[68] = 3.0               # Q10ICaL_ff
    p[69] = 3.0               # Q10ICaL_fs
    p[70] = 3.0               # Q10K
    p[71] = 1.5               # Q10Kb
    p[72] = 3.0               # Q10NCX
    p[73] = 1.63              # Q10NaK
    p[74] = 2.35              # Q10SLCaP
    p[75] = 2.6               # Q10SRCaP
    p[76] = 8314.0            # R
    p[77] = 310.0             # T
    p[78] = 25.0              # TOT_A
    p[79] = 1.65              # TRPN_n
    p[80] = 80.0              # Tref_b
    p[81] = 0.1371            # Whole_cell_PP1
    p[82] = 0.05              # aCaMK
    p[83] = 0.154375          # alpha_1
    p[84] = 0.00068           # bCaMK
    p[85] = 2.3               # beta_0
    p[86] = 0.1911            # beta_1
    p[87] = -2.4              # beta_1_mech
    p[88] = 4.75              # bt
    p[89] = 0.7645            # ca50_b
    p[90] = 1.7               # cajsr_half
    p[91] = 1.8               # cao
    p[92] = 2.0               # celltype
    p[93] = 150.0             # clo
    p[94] = 0.05              # cmdnmax_b
    p[95] = 10.0              # csqnmax
    p[96] = -0.155            # delta
    p[97] = 74.0              # dielConstant
    p[98] = 0.25              # dr
    p[99] = 4.2               # eP
    p[100] = 1.0              # fICaLP
    p[101] = 1.0              # fIKsP
    p[102] = 1.0              # fIKurP
    p[103] = 1.0              # fINaKP
    p[104] = 1.0              # fINaP
    p[105] = 1.0              # fPLBP
    p[106] = 1.0              # fRyRP
    p[107] = 1.0              # fTnIP
    p[108] = 0.0085           # gamma
    p[109] = 0.615            # gamma_wu
    p[110] = -53.0            # i_Stim_Amplitude
    p[111] = 1000.0           # i_Stim_Period
    p[112] = 1.0              # i_Stim_PulseDuration
    p[113] = 0.0              # i_Stim_Start
    p[114] = 1.0              # isHypoxic
    p[115] = 1.7              # jsrMidpoint
    p[116] = 182.4            # k1m
    p[117] = 949.5            # k1p
    p[118] = 39.4             # k2m
    p[119] = 500.0            # k2n
    p[120] = 687.2            # k2p
    p[121] = 79300.0          # k3m
    p[122] = 1899.0           # k3p
    p[123] = 40.0             # k4m
    p[124] = 639.0            # k4p
    p[125] = 0.026            # k_uw
    p[126] = 0.004            # k_ws
    p[127] = 12.5             # kasymm
    p[128] = 5000.0           # kcaoff
    p[129] = 1500000.0        # kcaon
    p[130] = 0.00238          # kmcmdn
    p[131] = 0.8              # kmcsqn
    p[132] = 0.0005           # kmtrpn_b
    p[133] = 15.0             # kna1
    p[134] = 5.0              # kna2
    p[135] = 88.12            # kna3
    p[136] = 5.0              # ko
    p[137] = 0.07854          # koff
    p[138] = 0.02626          # ktm_unblock
    p[139] = 1.0              # lambda
    p[140] = 1.2              # lambda_max
    p[141] = 0.87             # lambda_min
    p[142] = 0.0              # lambda_rate
    p[143] = 0.0              # mode
    p[144] = 3.94046          # mu_b
    p[145] = 1.795            # nK1
    p[146] = 0.991            # nNaCa
    p[147] = 1.65             # nTnI
    p[148] = 140.0            # nao
    p[149] = 2.036            # nperm_b
    p[150] = 1.87             # nrel
    p[151] = 10.15996         # nu_b
    p[152] = 1.14             # nup
    p[153] = 0.0              # offset
    p[154] = 7.2              # pH
    p[155] = 0.35             # perm50
    p[156] = 0.621            # ph_bt
    p[157] = 2.23             # phi
    p[158] = 6.86             # pkK1
    p[159] = 7.37             # pkNaCa
    p[160] = 6.79             # pkTnI
    p[161] = 6.64             # pkrel
    p[162] = 7.53             # pkup
    p[163] = 0.167            # qca
    p[164] = 0.5224           # qna
    p[165] = 0.0011           # rad_
    p[166] = 0.2              # tauCa
    p[167] = 2.0              # tauK
    p[168] = 2.0              # tauNa
    p[169] = 200.0            # thL
    p[170] = 72.5             # tjca
    p[171] = 0.07             # trpnmax
    p[172] = 0.0              # vShift
    p[173] = 60000.0          # wca
    p[174] = 0.5              # wfrac
    p[175] = 60000.0          # wna
    p[176] = 5000.0           # wnaca
    p[177] = 2.0              # zca
    p[178] = -1.0             # zcl
    p[179] = 1.0              # zk
    p[180] = 1.0              # zna
    p[181] = 1.0              # IKr_Multiplier_hetero
    return nothing
end
