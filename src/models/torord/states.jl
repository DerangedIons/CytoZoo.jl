const TORORD_NUM_STATES = 65

const TORORD_STATE_NAMES = (
    :v,          # 1 (transmembrane potential)
    :jca,        # 2
    :m,          # 3
    :mL,         # 4
    :xs1_P,      # 5
    :xs1,        # 6
    :a,          # 7
    :ap,         # 8
    :CaMKt,      # 9
    :hLp,        # 10
    :h_P,        # 11
    :h,          # 12
    :hp_P,       # 13
    :hp,         # 14
    :j,          # 15
    :j_P,        # 16
    :xs2,        # 17
    :C0,         # 18
    :C1,         # 19
    :O_,         # 20
    :clss,       # 21
    :cli,        # 22
    :iF,         # 23
    :iS,         # 24
    :TmBlocked,  # 25
    :jp_P,       # 26
    :jp,         # 27
    :d_P,        # 28
    :d,          # 29
    :fcaf_P,     # 30
    :fcaf,       # 31
    :ff_P,       # 32
    :ff_,        # 33
    :fcas_P,     # 34
    :fcas,       # 35
    :fs_P,       # 36
    :fs,         # 37
    :C2,         # 38
    :I_,         # 39
    :iFp,        # 40
    :iSp,        # 41
    :Zetas,      # 42
    :Zetaw,      # 43
    :Ca_TRPN,    # 44
    :fcaBPf,     # 45
    :fcafp,      # 46
    :fBPf,       # 47
    :ffp,        # 48
    :XS,         # 49
    :XW,         # 50
    :nca_i,      # 51
    :nca_ss,     # 52
    :cajsr,      # 53
    :cansr,      # 54
    :kss,        # 55
    :ki,         # 56
    :Jrel_np,    # 57
    :Jrel_p,     # 58
    :cai,        # 59
    :nai,        # 60
    :cass,       # 61
    :nass,       # 62
    :hL,         # 63
    :Jrel_p_P,   # 64
    :Jrel_np_P,  # 65
)

const TORORD_STATE_INDEX = Dict{Symbol,Int}(
    name => i for (i, name) in enumerate(TORORD_STATE_NAMES)
)

function _torord_init_state_values!(u)
    u[1] = -91.33918          # v (transmembrane potential)
    u[2] = 0.9999743          # jca
    u[3] = 0.0004619565       # m
    u[4] = 9.987709e-05       # mL
    u[5] = 0.0                # xs1_P
    u[6] = 0.2642293          # xs1
    u[7] = 0.0007994042       # a
    u[8] = 0.0004072777       # ap
    u[9] = 0.0201882          # CaMKt
    u[10] = 0.3339899         # hLp
    u[11] = 0.8562483         # h_P
    u[12] = 0.8739077         # h
    u[13] = 0.8560362         # hp_P
    u[14] = 0.7478972         # hp
    u[15] = 0.8737841         # j
    u[16] = 0.545264167       # j_P
    u[17] = 0.0001327348      # xs2
    u[18] = 0.9983451         # C0
    u[19] = 0.000708646       # C1
    u[20] = 0.0003445733      # O_
    u[21] = 48.91274          # clss
    u[22] = 48.91277          # cli
    u[23] = 0.9997514         # iF
    u[24] = 0.5702538         # iS
    u[25] = 1.0               # TmBlocked
    u[26] = 0.0               # jp_P
    u[27] = 0.8735375         # jp
    u[28] = 1.0               # d_P
    u[29] = -8.334604e-30     # d
    u[30] = 1.0               # fcaf_P
    u[31] = 1.0               # fcaf
    u[32] = 0.9334            # ff_P
    u[33] = 1.0               # ff_
    u[34] = 1.0               # fcas_P
    u[35] = 0.999754          # fcas
    u[36] = 1.0               # fs_P
    u[37] = 0.9183587         # fs
    u[38] = 0.0005910047      # C2
    u[39] = 1.06423e-05       # I_
    u[40] = 0.9997514         # iFp
    u[41] = 0.6351927         # iSp
    u[42] = 0.0               # Zetas
    u[43] = 0.0               # Zetaw
    u[44] = 0.0               # Ca_TRPN
    u[45] = 1.0               # fcaBPf
    u[46] = 1.0               # fcafp
    u[47] = 1.0               # fBPf
    u[48] = 1.0               # ffp
    u[49] = 0.0               # XS
    u[50] = 0.0               # XW
    u[51] = 0.001257861       # nca_i
    u[52] = 0.000533652       # nca_ss
    u[53] = 2.016415          # cajsr
    u[54] = 2.012225          # cansr
    u[55] = 156.713           # kss
    u[56] = 156.7131          # ki
    u[57] = -1.300486e-21     # Jrel_np
    u[58] = -7.610714e-20     # Jrel_p
    u[59] = 8.297576e-05      # cai
    u[60] = 15.94867          # nai
    u[61] = 6.642187e-05      # cass
    u[62] = 15.94922          # nass
    u[63] = 0.5986118         # hL
    u[64] = 8.236e-06         # Jrel_p_P
    u[65] = 0.0               # Jrel_np_P
    return nothing
end
