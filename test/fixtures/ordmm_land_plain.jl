
const NUM_STATES = 48;
const NUM_PARAMS = 139;
const NUM_MONITORED = 309;
# Parameter index
function parameter_index(name::String)

    if name == "Aff"
        return 1


    elseif name == "Ahf"
        return 2


    elseif name == "BSLmax"
        return 3


    elseif name == "BSRmax"
        return 4


    elseif name == "Beta0"
        return 5


    elseif name == "Beta1"
        return 6


    elseif name == "CaMKo"
        return 7


    elseif name == "Esac_ns"
        return 8


    elseif name == "F"
        return 9


    elseif name == "GKb"
        return 10


    elseif name == "GNa"
        return 11


    elseif name == "Gncx"
        return 12


    elseif name == "GpCa"
        return 13


    elseif name == "Gsac_k"
        return 14


    elseif name == "Gsac_ns"
        return 15


    elseif name == "Gto"
        return 16


    elseif name == "H"
        return 17


    elseif name == "Khp"
        return 18


    elseif name == "Kki"
        return 19


    elseif name == "Kko"
        return 20


    elseif name == "KmBSL"
        return 21


    elseif name == "KmBSR"
        return 22


    elseif name == "KmCaAct"
        return 23


    elseif name == "KmCaM"
        return 24


    elseif name == "KmCaMK"
        return 25


    elseif name == "Kmgatp"
        return 26


    elseif name == "Kmn"
        return 27


    elseif name == "Knai0"
        return 28


    elseif name == "Knao0"
        return 29


    elseif name == "Knap"
        return 30


    elseif name == "Kxkur"
        return 31


    elseif name == "L"
        return 32


    elseif name == "MgADP"
        return 33


    elseif name == "MgATP"
        return 34


    elseif name == "PCab"
        return 35


    elseif name == "PKNa"
        return 36


    elseif name == "PNab"
        return 37


    elseif name == "Pnak"
        return 38


    elseif name == "R"
        return 39


    elseif name == "T"
        return 40


    elseif name == "Tot_A"
        return 41


    elseif name == "Tref"
        return 42


    elseif name == "Trpn50"
        return 43


    elseif name == "aCaMK"
        return 44


    elseif name == "amp"
        return 45


    elseif name == "bCaMK"
        return 46


    elseif name == "bt"
        return 47


    elseif name == "calib"
        return 48


    elseif name == "cao"
        return 49


    elseif name == "cat50_ref"
        return 50


    elseif name == "celltype"
        return 51


    elseif name == "cmdnmax"
        return 52


    elseif name == "csqnmax"
        return 53


    elseif name == "dLambda"
        return 54


    elseif name == "delta"
        return 55


    elseif name == "delta_epi"
        return 56


    elseif name == "duration"
        return 57


    elseif name == "eP"
        return 58


    elseif name == "emcoupling"
        return 59


    elseif name == "etal"
        return 60


    elseif name == "etas"
        return 61


    elseif name == "gammas"
        return 62


    elseif name == "gammaw"
        return 63


    elseif name == "isacs"
        return 64


    elseif name == "k1m"
        return 65


    elseif name == "k1p"
        return 66


    elseif name == "k2m"
        return 67


    elseif name == "k2n"
        return 68


    elseif name == "k2p"
        return 69


    elseif name == "k3m"
        return 70


    elseif name == "k3p"
        return 71


    elseif name == "k4m"
        return 72


    elseif name == "k4p"
        return 73


    elseif name == "kasymm"
        return 74


    elseif name == "kcaoff"
        return 75


    elseif name == "kcaon"
        return 76


    elseif name == "kmcmdn"
        return 77


    elseif name == "kmcsqn"
        return 78


    elseif name == "kmtrpn"
        return 79


    elseif name == "kna1"
        return 80


    elseif name == "kna2"
        return 81


    elseif name == "kna3"
        return 82


    elseif name == "ko"
        return 83


    elseif name == "ktrpn"
        return 84


    elseif name == "ku"
        return 85


    elseif name == "kuw"
        return 86


    elseif name == "kws"
        return 87


    elseif name == "lambda_max"
        return 88


    elseif name == "lmbda"
        return 89


    elseif name == "mode"
        return 90


    elseif name == "nao"
        return 91


    elseif name == "ntm"
        return 92


    elseif name == "ntrpn"
        return 93


    elseif name == "p_a"
        return 94


    elseif name == "p_b"
        return 95


    elseif name == "p_k"
        return 96


    elseif name == "phi"
        return 97


    elseif name == "qca"
        return 98


    elseif name == "qna"
        return 99


    elseif name == "rad"
        return 100


    elseif name == "rs"
        return 101


    elseif name == "rw"
        return 102


    elseif name == "scale_HF_CaMKa"
        return 103


    elseif name == "scale_HF_GK1"
        return 104


    elseif name == "scale_HF_GNaL"
        return 105


    elseif name == "scale_HF_Gncx"
        return 106


    elseif name == "scale_HF_Gto"
        return 107


    elseif name == "scale_HF_Jleak"
        return 108


    elseif name == "scale_HF_Jrel_inf"
        return 109


    elseif name == "scale_HF_Jup"
        return 110


    elseif name == "scale_HF_Pnak"
        return 111


    elseif name == "scale_HF_cat50_ref"
        return 112


    elseif name == "scale_HF_thL"
        return 113


    elseif name == "scale_ICaL"
        return 114


    elseif name == "scale_IK1"
        return 115


    elseif name == "scale_IKr"
        return 116


    elseif name == "scale_IKs"
        return 117


    elseif name == "scale_INaL"
        return 118


    elseif name == "scale_drug_ICaL"
        return 119


    elseif name == "scale_drug_ICab"
        return 120


    elseif name == "scale_drug_IK1"
        return 121


    elseif name == "scale_drug_IKb"
        return 122


    elseif name == "scale_drug_IKr"
        return 123


    elseif name == "scale_drug_IKs"
        return 124


    elseif name == "scale_drug_INa"
        return 125


    elseif name == "scale_drug_INaL"
        return 126


    elseif name == "scale_drug_INab"
        return 127


    elseif name == "scale_drug_IpCa"
        return 128


    elseif name == "scale_drug_Isack"
        return 129


    elseif name == "scale_drug_Isacns"
        return 130


    elseif name == "scale_drug_Ito"
        return 131


    elseif name == "thL"
        return 132


    elseif name == "tjca"
        return 133


    elseif name == "trpnmax"
        return 134


    elseif name == "wca"
        return 135


    elseif name == "wna"
        return 136


    elseif name == "wnaca"
        return 137


    elseif name == "zca"
        return 138


    elseif name == "zk"
        return 139

    end
    return -1
end
# State index
function state_index(name::String)

    if name == "hL"
        return 1


    elseif name == "a"
        return 2


    elseif name == "ap"
        return 3


    elseif name == "d"
        return 4


    elseif name == "ff"
        return 5


    elseif name == "fs"
        return 6


    elseif name == "hf"
        return 7


    elseif name == "hs"
        return 8


    elseif name == "m"
        return 9


    elseif name == "xrf"
        return 10


    elseif name == "xrs"
        return 11


    elseif name == "xs1"
        return 12


    elseif name == "CaMKt"
        return 13


    elseif name == "xk1"
        return 14


    elseif name == "Zetaw"
        return 15


    elseif name == "XS"
        return 16


    elseif name == "XW"
        return 17


    elseif name == "TmB"
        return 18


    elseif name == "hLp"
        return 19


    elseif name == "iF"
        return 20


    elseif name == "iS"
        return 21


    elseif name == "fcaf"
        return 22


    elseif name == "fcas"
        return 23


    elseif name == "jca"
        return 24


    elseif name == "j"
        return 25


    elseif name == "fcafp"
        return 26


    elseif name == "ffp"
        return 27


    elseif name == "hsp"
        return 28


    elseif name == "jp"
        return 29


    elseif name == "mL"
        return 30


    elseif name == "xs2"
        return 31


    elseif name == "Zetas"
        return 32


    elseif name == "nca"
        return 33


    elseif name == "CaTrpn"
        return 34


    elseif name == "iFp"
        return 35


    elseif name == "iSp"
        return 36


    elseif name == "cajsr"
        return 37


    elseif name == "cansr"
        return 38


    elseif name == "kss"
        return 39


    elseif name == "Cd"
        return 40


    elseif name == "Jrelnp"
        return 41


    elseif name == "Jrelp"
        return 42


    elseif name == "ki"
        return 43


    elseif name == "cass"
        return 44


    elseif name == "nass"
        return 45


    elseif name == "cai"
        return 46


    elseif name == "nai"
        return 47


    elseif name == "v"
        return 48

    end
    return -1
end
# Monitor index
function monitor_index(name::String)

    if name == "zna"
        return 1


    elseif name == "Isac_P_k"
        return 2


    elseif name == "Isac_P_ns"
        return 3


    elseif name == "Afcaf"
        return 4


    elseif name == "AiF"
        return 5


    elseif name == "Axrf"
        return 6


    elseif name == "ass"
        return 7


    elseif name == "assp"
        return 8


    elseif name == "dss"
        return 9


    elseif name == "dti_develop"
        return 10


    elseif name == "dti_recover"
        return 11


    elseif name == "fss"
        return 12


    elseif name == "hLss"
        return 13


    elseif name == "hLssp"
        return 14


    elseif name == "hss"
        return 15


    elseif name == "hssp"
        return 16


    elseif name == "iss"
        return 17


    elseif name == "mLss"
        return 18


    elseif name == "mss"
        return 19


    elseif name == "rkr"
        return 20


    elseif name == "ta"
        return 21


    elseif name == "td"
        return 22


    elseif name == "tfcaf"
        return 23


    elseif name == "tfcas"
        return 24


    elseif name == "tff"
        return 25


    elseif name == "tfs"
        return 26


    elseif name == "thf"
        return 27


    elseif name == "ths"
        return 28


    elseif name == "tj"
        return 29


    elseif name == "tm"
        return 30


    elseif name == "txk1"
        return 31


    elseif name == "txrf"
        return 32


    elseif name == "txrs"
        return 33


    elseif name == "txs1"
        return 34


    elseif name == "txs2"
        return 35


    elseif name == "xkb"
        return 36


    elseif name == "xrss"
        return 37


    elseif name == "xs1ss"
        return 38


    elseif name == "Afs"
        return 39


    elseif name == "Ageo"
        return 40


    elseif name == "vcell"
        return 41


    elseif name == "Ahs"
        return 42


    elseif name == "Aw"
        return 43


    elseif name == "Jupnp"
        return 44


    elseif name == "Jupp"
        return 45


    elseif name == "KsCa"
        return 46


    elseif name == "Bcai"
        return 47


    elseif name == "Bcajsr"
        return 48


    elseif name == "Bcass"
        return 49


    elseif name == "Jdiff"
        return 50


    elseif name == "CaMKb"
        return 51


    elseif name == "CaTrpn_max"
        return 52


    elseif name == "vffrt"
        return 53


    elseif name == "vfrt"
        return 54


    elseif name == "EK"
        return 55


    elseif name == "rk1"
        return 56


    elseif name == "xk1ss"
        return 57


    elseif name == "EKs"
        return 58


    elseif name == "ENa"
        return 59


    elseif name == "GK1"
        return 60


    elseif name == "GKr"
        return 61


    elseif name == "GKs"
        return 62


    elseif name == "GNaL"
        return 63


    elseif name == "km2n"
        return 64


    elseif name == "IpCa"
        return 65


    elseif name == "Istim"
        return 66


    elseif name == "JdiffK"
        return 67


    elseif name == "JdiffNa"
        return 68


    elseif name == "Jtr"
        return 69


    elseif name == "Jleak"
        return 70


    elseif name == "Knai"
        return 71


    elseif name == "Knao"
        return 72


    elseif name == "P"
        return 73


    elseif name == "PCa"
        return 74


    elseif name == "XS_max"
        return 75


    elseif name == "XW_max"
        return 76


    elseif name == "XU"
        return 77


    elseif name == "a2"
        return 78


    elseif name == "a4"
        return 79


    elseif name == "a_rel"
        return 80


    elseif name == "btp"
        return 81


    elseif name == "tau_rel_tmp"
        return 82


    elseif name == "allo_i"
        return 83


    elseif name == "allo_ss"
        return 84


    elseif name == "b1"
        return 85


    elseif name == "ksu"
        return 86


    elseif name == "cs"
        return 87


    elseif name == "cw"
        return 88


    elseif name == "kwu"
        return 89


    elseif name == "gammasu"
        return 90


    elseif name == "gammawu"
        return 91


    elseif name == "h10"
        return 92


    elseif name == "h10_i"
        return 93


    elseif name == "h4"
        return 94


    elseif name == "h4_i"
        return 95


    elseif name == "hca"
        return 96


    elseif name == "hna"
        return 97


    elseif name == "k2"
        return 98


    elseif name == "k2_i"
        return 99


    elseif name == "k5"
        return 100


    elseif name == "k5_i"
        return 101


    elseif name == "kb"
        return 102


    elseif name == "lambda_min12"
        return 103


    elseif name == "thLp"
        return 104


    elseif name == "tiF"
        return 105


    elseif name == "tiS"
        return 106


    elseif name == "Afcas"
        return 107


    elseif name == "AiS"
        return 108


    elseif name == "Axrs"
        return 109


    elseif name == "fcass"
        return 110


    elseif name == "dhL_dt"
        return 111


    elseif name == "jss"
        return 112


    elseif name == "da_dt"
        return 113


    elseif name == "dap_dt"
        return 114


    elseif name == "dd_dt"
        return 115


    elseif name == "tfcafp"
        return 116


    elseif name == "tffp"
        return 117


    elseif name == "dff_dt"
        return 118


    elseif name == "dfs_dt"
        return 119


    elseif name == "dhf_dt"
        return 120


    elseif name == "thsp"
        return 121


    elseif name == "dhs_dt"
        return 122


    elseif name == "tjp"
        return 123


    elseif name == "tmL"
        return 124


    elseif name == "dm_dt"
        return 125


    elseif name == "dxrf_dt"
        return 126


    elseif name == "dxrs_dt"
        return 127


    elseif name == "xs2ss"
        return 128


    elseif name == "dxs1_dt"
        return 129


    elseif name == "f"
        return 130


    elseif name == "fp"
        return 131


    elseif name == "Acap"
        return 132


    elseif name == "vjsr"
        return 133


    elseif name == "vmyo"
        return 134


    elseif name == "vnsr"
        return 135


    elseif name == "vss"
        return 136


    elseif name == "h"
        return 137


    elseif name == "hp"
        return 138


    elseif name == "As"
        return 139


    elseif name == "CaMKa"
        return 140


    elseif name == "dCaMKt_dt"
        return 141


    elseif name == "ICab"
        return 142


    elseif name == "INab"
        return 143


    elseif name == "PhiCaK"
        return 144


    elseif name == "PhiCaL"
        return 145


    elseif name == "PhiCaNa"
        return 146


    elseif name == "IKb"
        return 147


    elseif name == "dxk1_dt"
        return 148


    elseif name == "IK1"
        return 149


    elseif name == "IKs"
        return 150


    elseif name == "anca"
        return 151


    elseif name == "a1"
        return 152


    elseif name == "b4"
        return 153


    elseif name == "a3"
        return 154


    elseif name == "b2"
        return 155


    elseif name == "b3"
        return 156


    elseif name == "PCaK"
        return 157


    elseif name == "PCaNa"
        return 158


    elseif name == "PCap"
        return 159


    elseif name == "a_relp"
        return 160


    elseif name == "tau_relp_tmp"
        return 161


    elseif name == "tau_rel"
        return 162


    elseif name == "dZetaw_dt"
        return 163


    elseif name == "dXS_dt"
        return 164


    elseif name == "dXW_dt"
        return 165


    elseif name == "h11"
        return 166


    elseif name == "h12"
        return 167


    elseif name == "h11_i"
        return 168


    elseif name == "h12_i"
        return 169


    elseif name == "h5"
        return 170


    elseif name == "h6"
        return 171


    elseif name == "h5_i"
        return 172


    elseif name == "h6_i"
        return 173


    elseif name == "h1"
        return 174


    elseif name == "h1_i"
        return 175


    elseif name == "h7"
        return 176


    elseif name == "h7_i"
        return 177


    elseif name == "dTmB_dt"
        return 178


    elseif name == "C"
        return 179


    elseif name == "cat50"
        return 180


    elseif name == "lambda_min087"
        return 181


    elseif name == "dhLp_dt"
        return 182


    elseif name == "tiFp"
        return 183


    elseif name == "diF_dt"
        return 184


    elseif name == "tiSp"
        return 185


    elseif name == "diS_dt"
        return 186


    elseif name == "fca"
        return 187


    elseif name == "fcap"
        return 188


    elseif name == "i"
        return 189


    elseif name == "ip"
        return 190


    elseif name == "xr"
        return 191


    elseif name == "dfcaf_dt"
        return 192


    elseif name == "dfcas_dt"
        return 193


    elseif name == "djca_dt"
        return 194


    elseif name == "dj_dt"
        return 195


    elseif name == "dfcafp_dt"
        return 196


    elseif name == "dffp_dt"
        return 197


    elseif name == "dhsp_dt"
        return 198


    elseif name == "djp_dt"
        return 199


    elseif name == "dmL_dt"
        return 200


    elseif name == "dxs2_dt"
        return 201


    elseif name == "dZetas_dt"
        return 202


    elseif name == "fICaLp"
        return 203


    elseif name == "fINaLp"
        return 204


    elseif name == "fINap"
        return 205


    elseif name == "fItop"
        return 206


    elseif name == "fJrelp"
        return 207


    elseif name == "fJupp"
        return 208


    elseif name == "dnca_dt"
        return 209


    elseif name == "x2"
        return 210


    elseif name == "x1"
        return 211


    elseif name == "x3"
        return 212


    elseif name == "x4"
        return 213


    elseif name == "PCaKp"
        return 214


    elseif name == "PCaNap"
        return 215


    elseif name == "tau_relp"
        return 216


    elseif name == "k1"
        return 217


    elseif name == "k1_i"
        return 218


    elseif name == "k6"
        return 219


    elseif name == "k6_i"
        return 220


    elseif name == "h2"
        return 221


    elseif name == "h3"
        return 222


    elseif name == "h2_i"
        return 223


    elseif name == "h3_i"
        return 224


    elseif name == "h8"
        return 225


    elseif name == "h9"
        return 226


    elseif name == "h8_i"
        return 227


    elseif name == "h9_i"
        return 228


    elseif name == "F1"
        return 229


    elseif name == "dCd"
        return 230


    elseif name == "dCaTrpn_dt"
        return 231


    elseif name == "h_lambda_prima"
        return 232


    elseif name == "diFp_dt"
        return 233


    elseif name == "diSp_dt"
        return 234


    elseif name == "IKr"
        return 235


    elseif name == "ICaL"
        return 236


    elseif name == "INaL"
        return 237


    elseif name == "INa"
        return 238


    elseif name == "Ito"
        return 239


    elseif name == "Jrel"
        return 240


    elseif name == "Jup"
        return 241


    elseif name == "E1"
        return 242


    elseif name == "E2"
        return 243


    elseif name == "E3"
        return 244


    elseif name == "E4"
        return 245


    elseif name == "ICaK"
        return 246


    elseif name == "ICaNa"
        return 247


    elseif name == "k4pp"
        return 248


    elseif name == "k7"
        return 249


    elseif name == "k4p_ss"
        return 250


    elseif name == "k4pp_i"
        return 251


    elseif name == "k7_i"
        return 252


    elseif name == "k4p_i"
        return 253


    elseif name == "k3pp"
        return 254


    elseif name == "k8"
        return 255


    elseif name == "k3p_ss"
        return 256


    elseif name == "k3pp_i"
        return 257


    elseif name == "k8_i"
        return 258


    elseif name == "k3p_i"
        return 259


    elseif name == "eta"
        return 260


    elseif name == "J_TRPN"
        return 261


    elseif name == "h_lambda"
        return 262


    elseif name == "Jrel_inf"
        return 263


    elseif name == "Jrel_infp"
        return 264


    elseif name == "dcajsr_dt"
        return 265


    elseif name == "dcansr_dt"
        return 266


    elseif name == "JnakNa"
        return 267


    elseif name == "JnakK"
        return 268


    elseif name == "dkss_dt"
        return 269


    elseif name == "k4"
        return 270


    elseif name == "k4_i"
        return 271


    elseif name == "k3"
        return 272


    elseif name == "k3_i"
        return 273


    elseif name == "Fd"
        return 274


    elseif name == "dCd_dt"
        return 275


    elseif name == "Ta"
        return 276


    elseif name == "dJrelnp_dt"
        return 277


    elseif name == "dJrelp_dt"
        return 278


    elseif name == "INaK"
        return 279


    elseif name == "x2_ss"
        return 280


    elseif name == "x2_i"
        return 281


    elseif name == "x1_ss"
        return 282


    elseif name == "x3_ss"
        return 283


    elseif name == "x4_ss"
        return 284


    elseif name == "x1_i"
        return 285


    elseif name == "x3_i"
        return 286


    elseif name == "x4_i"
        return 287


    elseif name == "Tp"
        return 288


    elseif name == "dki_dt"
        return 289


    elseif name == "E1_ss"
        return 290


    elseif name == "E2_ss"
        return 291


    elseif name == "E3_ss"
        return 292


    elseif name == "E4_ss"
        return 293


    elseif name == "E1_i"
        return 294


    elseif name == "E2_i"
        return 295


    elseif name == "E3_i"
        return 296


    elseif name == "E4_i"
        return 297


    elseif name == "Ttot"
        return 298


    elseif name == "JncxCa_ss"
        return 299


    elseif name == "JncxNa_ss"
        return 300


    elseif name == "JncxCa_i"
        return 301


    elseif name == "JncxNa_i"
        return 302


    elseif name == "INaCa_ss"
        return 303


    elseif name == "INaCa_i"
        return 304


    elseif name == "dcass_dt"
        return 305


    elseif name == "dnass_dt"
        return 306


    elseif name == "dcai_dt"
        return 307


    elseif name == "dnai_dt"
        return 308


    elseif name == "dv_dt"
        return 309

    end
    return -1
end


function init_parameter_values!(parameters)
    #=
    Aff=0.6, Ahf=0.99, BSLmax=1.124, BSRmax=0.047, Beta0=2.3, Beta1=-2.4, CaMKo=0.05, Esac_ns=-10, F=96485.0, GKb=0.003, GNa=31, Gncx=0.0008, GpCa=0.0005, Gsac_k=(0.2882 * 800) / 210, Gsac_ns=0.006, Gto=0.02, H=1e-07, Khp=1.698e-07, Kki=0.5, Kko=0.3582, KmBSL=0.0087, KmBSR=0.00087, KmCaAct=0.00015, KmCaM=0.0015, KmCaMK=0.15, Kmgatp=1.698e-07, Kmn=0.002, Knai0=9.073, Knao0=27.78, Knap=224.0, Kxkur=292.0, L=0.01, MgADP=0.05, MgATP=9.8, PCab=2.5e-08, PKNa=0.01833, PNab=3.75e-10, Pnak=30, R=8314.0, T=310.0, Tot_A=25, Tref=120, Trpn50=0.35, aCaMK=0.05, amp=-80.0, bCaMK=0.00068, bt=4.75, calib=1, cao=1.8, cat50_ref=0.805, celltype=0, cmdnmax=0.05, csqnmax=10.0, dLambda=0, delta=-0.155, delta_epi=1.0, duration=0.5, eP=4.2, emcoupling=1, etal=200, etas=20, gammas=0.0085, gammaw=0.615, isacs=0, k1m=182.4, k1p=949.5, k2m=39.4, k2n=1000.0, k2p=687.2, k3m=79300.0, k3p=1899.0, k4m=40.0, k4p=639.0, kasymm=12.5, kcaoff=5000.0, kcaon=1500000.0, kmcmdn=0.00238, kmcsqn=0.8, kmtrpn=0.0005, kna1=15.0, kna2=5.0, kna3=88.12, ko=5.4, ktrpn=0.1, ku=0.04, kuw=0.182, kws=0.012, lambda_max=1.1, lmbda=1, mode=1, nao=140.0, ntm=2.4, ntrpn=2, p_a=2.1, p_b=9.1, p_k=7, phi=2.23, qca=0.167, qna=0.5224, rad=0.0011, rs=0.25, rw=0.5, scale_HF_CaMKa=1.0, scale_HF_GK1=1.0, scale_HF_GNaL=1.0, scale_HF_Gncx=1.0, scale_HF_Gto=1.0, scale_HF_Jleak=1.0, scale_HF_Jrel_inf=1.0, scale_HF_Jup=1.0, scale_HF_Pnak=1.0, scale_HF_cat50_ref=1.0, scale_HF_thL=1.0, scale_ICaL=1.018, scale_IK1=1.414, scale_IKr=1.119, scale_IKs=1.648, scale_INaL=2.274, scale_drug_ICaL=1.0, scale_drug_ICab=1.0, scale_drug_IK1=1.0, scale_drug_IKb=1.0, scale_drug_IKr=1.0, scale_drug_IKs=1.0, scale_drug_INa=1.0, scale_drug_INaL=1.0, scale_drug_INab=1.0, scale_drug_IpCa=1.0, scale_drug_Isack=1.0, scale_drug_Isacns=1.0, scale_drug_Ito=1.0, thL=200.0, tjca=75.0, trpnmax=0.07, wca=60000.0, wna=60000.0, wnaca=5000.0, zca=2.0, zk=1.0
    =#
    parameters[1] = 0.6
    parameters[2] = 0.99
    parameters[3] = 1.124
    parameters[4] = 0.047
    parameters[5] = 2.3
    parameters[6] = -2.4
    parameters[7] = 0.05
    parameters[8] = -10
    parameters[9] = 96485.0
    parameters[10] = 0.003
    parameters[11] = 31
    parameters[12] = 0.0008
    parameters[13] = 0.0005
    parameters[14] = (0.2882 * 800) / 210
    parameters[15] = 0.006
    parameters[16] = 0.02
    parameters[17] = 1e-07
    parameters[18] = 1.698e-07
    parameters[19] = 0.5
    parameters[20] = 0.3582
    parameters[21] = 0.0087
    parameters[22] = 0.00087
    parameters[23] = 0.00015
    parameters[24] = 0.0015
    parameters[25] = 0.15
    parameters[26] = 1.698e-07
    parameters[27] = 0.002
    parameters[28] = 9.073
    parameters[29] = 27.78
    parameters[30] = 224.0
    parameters[31] = 292.0
    parameters[32] = 0.01
    parameters[33] = 0.05
    parameters[34] = 9.8
    parameters[35] = 2.5e-08
    parameters[36] = 0.01833
    parameters[37] = 3.75e-10
    parameters[38] = 30
    parameters[39] = 8314.0
    parameters[40] = 310.0
    parameters[41] = 25
    parameters[42] = 120
    parameters[43] = 0.35
    parameters[44] = 0.05
    parameters[45] = -80.0
    parameters[46] = 0.00068
    parameters[47] = 4.75
    parameters[48] = 1
    parameters[49] = 1.8
    parameters[50] = 0.805
    parameters[51] = 0
    parameters[52] = 0.05
    parameters[53] = 10.0
    parameters[54] = 0
    parameters[55] = -0.155
    parameters[56] = 1.0
    parameters[57] = 0.5
    parameters[58] = 4.2
    parameters[59] = 1
    parameters[60] = 200
    parameters[61] = 20
    parameters[62] = 0.0085
    parameters[63] = 0.615
    parameters[64] = 0
    parameters[65] = 182.4
    parameters[66] = 949.5
    parameters[67] = 39.4
    parameters[68] = 1000.0
    parameters[69] = 687.2
    parameters[70] = 79300.0
    parameters[71] = 1899.0
    parameters[72] = 40.0
    parameters[73] = 639.0
    parameters[74] = 12.5
    parameters[75] = 5000.0
    parameters[76] = 1500000.0
    parameters[77] = 0.00238
    parameters[78] = 0.8
    parameters[79] = 0.0005
    parameters[80] = 15.0
    parameters[81] = 5.0
    parameters[82] = 88.12
    parameters[83] = 5.4
    parameters[84] = 0.1
    parameters[85] = 0.04
    parameters[86] = 0.182
    parameters[87] = 0.012
    parameters[88] = 1.1
    parameters[89] = 1
    parameters[90] = 1
    parameters[91] = 140.0
    parameters[92] = 2.4
    parameters[93] = 2
    parameters[94] = 2.1
    parameters[95] = 9.1
    parameters[96] = 7
    parameters[97] = 2.23
    parameters[98] = 0.167
    parameters[99] = 0.5224
    parameters[100] = 0.0011
    parameters[101] = 0.25
    parameters[102] = 0.5
    parameters[103] = 1.0
    parameters[104] = 1.0
    parameters[105] = 1.0
    parameters[106] = 1.0
    parameters[107] = 1.0
    parameters[108] = 1.0
    parameters[109] = 1.0
    parameters[110] = 1.0
    parameters[111] = 1.0
    parameters[112] = 1.0
    parameters[113] = 1.0
    parameters[114] = 1.018
    parameters[115] = 1.414
    parameters[116] = 1.119
    parameters[117] = 1.648
    parameters[118] = 2.274
    parameters[119] = 1.0
    parameters[120] = 1.0
    parameters[121] = 1.0
    parameters[122] = 1.0
    parameters[123] = 1.0
    parameters[124] = 1.0
    parameters[125] = 1.0
    parameters[126] = 1.0
    parameters[127] = 1.0
    parameters[128] = 1.0
    parameters[129] = 1.0
    parameters[130] = 1.0
    parameters[131] = 1.0
    parameters[132] = 200.0
    parameters[133] = 75.0
    parameters[134] = 0.07
    parameters[135] = 60000.0
    parameters[136] = 60000.0
    parameters[137] = 5000.0
    parameters[138] = 2.0
    parameters[139] = 1.0
end


function init_state_values!(states)
    #=
    hL=1, a=0, ap=0, d=0, ff=1, fs=1, hf=1, hs=1, m=0, xrf=0, xrs=0, xs1=0, CaMKt=0, xk1=1, Zetaw=0, XS=0, XW=0, TmB=1, hLp=1, iF=1, iS=1, fcaf=1, fcas=1, jca=1, j=1, fcafp=1, ffp=1, hsp=1, jp=1, mL=0, xs2=0, Zetas=0, nca=0, CaTrpn=0, iFp=1, iSp=1, cajsr=1.2, cansr=1.2, kss=145, Cd=0, Jrelnp=0, Jrelp=0, ki=145, cass=0.0001, nass=7, cai=0.0001, nai=7, v=-87
    =#
    states[1] = 1
    states[2] = 0
    states[3] = 0
    states[4] = 0
    states[5] = 1
    states[6] = 1
    states[7] = 1
    states[8] = 1
    states[9] = 0
    states[10] = 0
    states[11] = 0
    states[12] = 0
    states[13] = 0
    states[14] = 1
    states[15] = 0
    states[16] = 0
    states[17] = 0
    states[18] = 1
    states[19] = 1
    states[20] = 1
    states[21] = 1
    states[22] = 1
    states[23] = 1
    states[24] = 1
    states[25] = 1
    states[26] = 1
    states[27] = 1
    states[28] = 1
    states[29] = 1
    states[30] = 0
    states[31] = 0
    states[32] = 0
    states[33] = 0
    states[34] = 0
    states[35] = 1
    states[36] = 1
    states[37] = 1.2
    states[38] = 1.2
    states[39] = 145
    states[40] = 0
    states[41] = 0
    states[42] = 0
    states[43] = 145
    states[44] = 0.0001
    states[45] = 7
    states[46] = 0.0001
    states[47] = 7
    states[48] = -87
end


function rhs(t, states, parameters, values)

    # Assign states
    hL = states[1]
    a = states[2]
    ap = states[3]
    d = states[4]
    ff = states[5]
    fs = states[6]
    hf = states[7]
    hs = states[8]
    m = states[9]
    xrf = states[10]
    xrs = states[11]
    xs1 = states[12]
    CaMKt = states[13]
    xk1 = states[14]
    Zetaw = states[15]
    XS = states[16]
    XW = states[17]
    TmB = states[18]
    hLp = states[19]
    iF = states[20]
    iS = states[21]
    fcaf = states[22]
    fcas = states[23]
    jca = states[24]
    j = states[25]
    fcafp = states[26]
    ffp = states[27]
    hsp = states[28]
    jp = states[29]
    mL = states[30]
    xs2 = states[31]
    Zetas = states[32]
    nca = states[33]
    CaTrpn = states[34]
    iFp = states[35]
    iSp = states[36]
    cajsr = states[37]
    cansr = states[38]
    kss = states[39]
    Cd = states[40]
    Jrelnp = states[41]
    Jrelp = states[42]
    ki = states[43]
    cass = states[44]
    nass = states[45]
    cai = states[46]
    nai = states[47]
    v = states[48]

    # Assign parameters
    Aff = parameters[1]
    Ahf = parameters[2]
    BSLmax = parameters[3]
    BSRmax = parameters[4]
    Beta0 = parameters[5]
    Beta1 = parameters[6]
    CaMKo = parameters[7]
    Esac_ns = parameters[8]
    F = parameters[9]
    GKb = parameters[10]
    GNa = parameters[11]
    Gncx = parameters[12]
    GpCa = parameters[13]
    Gsac_k = parameters[14]
    Gsac_ns = parameters[15]
    Gto = parameters[16]
    H = parameters[17]
    Khp = parameters[18]
    Kki = parameters[19]
    Kko = parameters[20]
    KmBSL = parameters[21]
    KmBSR = parameters[22]
    KmCaAct = parameters[23]
    KmCaM = parameters[24]
    KmCaMK = parameters[25]
    Kmgatp = parameters[26]
    Kmn = parameters[27]
    Knai0 = parameters[28]
    Knao0 = parameters[29]
    Knap = parameters[30]
    Kxkur = parameters[31]
    L = parameters[32]
    MgADP = parameters[33]
    MgATP = parameters[34]
    PCab = parameters[35]
    PKNa = parameters[36]
    PNab = parameters[37]
    Pnak = parameters[38]
    R = parameters[39]
    T = parameters[40]
    Tot_A = parameters[41]
    Tref = parameters[42]
    Trpn50 = parameters[43]
    aCaMK = parameters[44]
    amp = parameters[45]
    bCaMK = parameters[46]
    bt = parameters[47]
    calib = parameters[48]
    cao = parameters[49]
    cat50_ref = parameters[50]
    celltype = parameters[51]
    cmdnmax = parameters[52]
    csqnmax = parameters[53]
    dLambda = parameters[54]
    delta = parameters[55]
    delta_epi = parameters[56]
    duration = parameters[57]
    eP = parameters[58]
    emcoupling = parameters[59]
    etal = parameters[60]
    etas = parameters[61]
    gammas = parameters[62]
    gammaw = parameters[63]
    isacs = parameters[64]
    k1m = parameters[65]
    k1p = parameters[66]
    k2m = parameters[67]
    k2n = parameters[68]
    k2p = parameters[69]
    k3m = parameters[70]
    k3p = parameters[71]
    k4m = parameters[72]
    k4p = parameters[73]
    kasymm = parameters[74]
    kcaoff = parameters[75]
    kcaon = parameters[76]
    kmcmdn = parameters[77]
    kmcsqn = parameters[78]
    kmtrpn = parameters[79]
    kna1 = parameters[80]
    kna2 = parameters[81]
    kna3 = parameters[82]
    ko = parameters[83]
    ktrpn = parameters[84]
    ku = parameters[85]
    kuw = parameters[86]
    kws = parameters[87]
    lambda_max = parameters[88]
    lmbda = parameters[89]
    mode = parameters[90]
    nao = parameters[91]
    ntm = parameters[92]
    ntrpn = parameters[93]
    p_a = parameters[94]
    p_b = parameters[95]
    p_k = parameters[96]
    phi = parameters[97]
    qca = parameters[98]
    qna = parameters[99]
    rad = parameters[100]
    rs = parameters[101]
    rw = parameters[102]
    scale_HF_CaMKa = parameters[103]
    scale_HF_GK1 = parameters[104]
    scale_HF_GNaL = parameters[105]
    scale_HF_Gncx = parameters[106]
    scale_HF_Gto = parameters[107]
    scale_HF_Jleak = parameters[108]
    scale_HF_Jrel_inf = parameters[109]
    scale_HF_Jup = parameters[110]
    scale_HF_Pnak = parameters[111]
    scale_HF_cat50_ref = parameters[112]
    scale_HF_thL = parameters[113]
    scale_ICaL = parameters[114]
    scale_IK1 = parameters[115]
    scale_IKr = parameters[116]
    scale_IKs = parameters[117]
    scale_INaL = parameters[118]
    scale_drug_ICaL = parameters[119]
    scale_drug_ICab = parameters[120]
    scale_drug_IK1 = parameters[121]
    scale_drug_IKb = parameters[122]
    scale_drug_IKr = parameters[123]
    scale_drug_IKs = parameters[124]
    scale_drug_INa = parameters[125]
    scale_drug_INaL = parameters[126]
    scale_drug_INab = parameters[127]
    scale_drug_IpCa = parameters[128]
    scale_drug_Isack = parameters[129]
    scale_drug_Isacns = parameters[130]
    scale_drug_Ito = parameters[131]
    thL = parameters[132]
    tjca = parameters[133]
    trpnmax = parameters[134]
    wca = parameters[135]
    wna = parameters[136]
    wnaca = parameters[137]
    zca = parameters[138]
    zk = parameters[139]

    # Assign expressions
    zna = 1.0
    Isac_P_k = 0
    Isac_P_ns = 0
    Afcaf = 0.3 + 0.6 ./ (exp((v - 10.0) / 10.0) + 1.0)
    AiF = 1.0 ./ (exp((v - 213.6) / 151.2) + 1.0)
    Axrf = 1.0 ./ (exp((v + 54.81) / 38.21) + 1.0)
    ass = 1.0 ./ (exp((-(v - 14.34)) / 14.82) + 1.0)
    assp = 1.0 ./ (exp((-(v - 24.34)) / 14.82) + 1.0)
    dss = 1.0 ./ (exp((-(v + 3.94)) / 4.23) + 1.0)
    dti_develop = 1.354 + 0.0001 ./ (exp((-(v - 12.23)) / 0.2154) + exp((v - 167.4) / 15.89))
    dti_recover = 1.0 - 0.5 ./ (exp((v + 70.0) / 20.0) + 1.0)
    fss = 1.0 ./ (exp((v + 19.58) / 3.696) + 1.0)
    hLss = 1.0 ./ (exp((v + 87.61) / 7.488) + 1.0)
    hLssp = 1.0 ./ (exp((v + 93.81) / 7.488) + 1.0)
    hss = 1.0 ./ (exp((v + 78.5) / 6.22) + 1)
    hssp = 1.0 ./ (exp(((v + 78.5) + 6.2) / 6.22) + 1)
    iss = 1.0 ./ (exp((v + 43.94) / 5.711) + 1.0)
    mLss = 1.0 ./ (exp((-(v + 42.85)) / 5.264) + 1.0)
    mss = 1.0 ./ (exp((-((v + 39.57) + 9.4)) / 7.5) + 1.0)
    rkr = (1.0 * (1.0 ./ (exp((v + 55.0) / 75.0) + 1.0))) ./ (exp((v - 10.0) / 30.0) + 1.0)
    ta = 1.0515 ./ (1.0 ./ ((1.2089 * (exp((-(v - 18.4099)) / 29.3814) + 1.0))) + 3.5 ./ (exp((v + 100.0) / 29.3814) + 1.0))
    td = 0.6 + 1.0 ./ (exp((-0.05) * (v + 6.0)) + exp(0.09 * (v + 14.0)))
    tfcaf = 7.0 + 1.0 ./ (0.04 * exp((-(v - 4.0)) / 7.0) + 0.04 * exp((v - 4.0) / 7.0))
    tfcas = 100.0 + 1.0 ./ (0.00012 * exp((-v) / 3.0) + 0.00012 * exp(v / 7.0))
    tff = 7.0 + 1.0 ./ (0.0045 * exp((-(v + 20.0)) / 10.0) + 0.0045 * exp((v + 20.0) / 10.0))
    tfs = 1000.0 + 1.0 ./ (3.5e-05 * exp((-(v + 5.0)) / 4.0) + 3.5e-05 * exp((v + 5.0) / 6.0))
    thf = 1.0 ./ (6.149 * exp((v + 0.5096) / 20.27) + 1.432e-05 * exp((-(v + 1.196)) / 6.285))
    ths = 1.0 ./ (0.009794 * exp((-(v + 17.95)) / 28.05) + 0.3343 * exp((v + 5.73) / 56.66))
    tj = 2.038 + 1.0 ./ (0.3052 * exp((v + 0.9941) / 38.45) + 0.02136 * exp((-(v + 100.6)) / 8.281))
    tm = 1.0 ./ (6.765 * exp((v + 11.64) / 34.77) + 8.552 * exp((-(v + 77.42)) / 5.955))
    txk1 = 122.2 ./ (exp((-(v + 127.2)) / 20.36) + exp((v + 236.8) / 69.33))
    txrf = 12.98 + 1.0 ./ (4.123e-05 * exp((-(v - 47.78)) / 20.38) + 0.3652 * exp((v - 31.66) / 3.869))
    txrs = 1.865 + 1.0 ./ (1.128e-05 * exp((-(v - 29.74)) / 25.94) + 0.06629 * exp((v - 34.7) / 7.355))
    txs1 = 817.3 + 1.0 ./ (0.0002326 * exp((v + 48.28) / 17.8) + 0.001292 * exp((-(v + 210.0)) / 230.0))
    txs2 = 1.0 ./ (0.01 * exp((v - 50.0) / 20.0) + 0.0193 * exp((-(v + 66.54)) / 31.0))
    xkb = 1.0 ./ (exp((-(v - 14.48)) / 18.34) + 1.0)
    xrss = 1.0 ./ (exp((-(v + 8.337)) / 6.789) + 1.0)
    xs1ss = 1.0 ./ (exp((-(v + 11.6)) / 8.932) + 1.0)
    Afs = 1.0 - Aff
    Ageo = L .* ((2 * 3.14) * rad) + rad .* ((2 * 3.14) * rad)
    vcell = L .* (rad .* ((3.14 * 1000) * rad))
    Ahs = 1.0 - Ahf
    Aw = (Tot_A .* rs) ./ (rs + rw .* (1 - rs))
    Jupnp = (0.004375 * cai) ./ (cai + 0.00092)
    Jupp = ((0.004375 * 2.75) * cai) ./ ((cai + 0.00092) - 0.00017)
    KsCa = 1.0 + 0.6 ./ ((3.8e-05 ./ cai) .^ 1.4 + 1.0)
    Bcai = 1.0 ./ ((cmdnmax .* kmcmdn) ./ (cai + kmcmdn) .^ 2.0 + 1.0)
    Bcajsr = 1.0 ./ ((csqnmax .* kmcsqn) ./ (cajsr + kmcsqn) .^ 2.0 + 1.0)
    Bcass = 1.0 ./ ((BSLmax .* KmBSL) ./ (KmBSL + cass) .^ 2.0 + ((BSRmax .* KmBSR) ./ (KmBSR + cass) .^ 2.0 + 1.0))
    Jdiff = (-cai + cass) / 0.2
    CaMKb = (CaMKo .* (1.0 - CaMKt)) ./ (KmCaM ./ cass + 1.0)
    CaTrpn_max = ((CaTrpn > 0) ? (CaTrpn) : (0))
    vffrt = (F .* (F .* v)) ./ ((R .* T))
    vfrt = (F .* v) ./ ((R .* T))
    EK = ((R .* T) ./ F) .* log(ko ./ ki)
    rk1 = 1.0 ./ (exp((-2.6 * ko + (v + 105.8)) / 9.493) + 1.0)
    xk1ss = 1.0 ./ (exp((-((2.5538 * ko + v) + 144.59)) ./ (1.5692 * ko + 3.8115)) + 1.0)
    EKs = ((R .* T) ./ F) .* log((PKNa .* nao + ko) ./ (PKNa .* nai + ki))
    ENa = ((R .* T) ./ F) .* log(nao ./ nai)
    GK1 = scale_HF_GK1 .* ((0.1908 * scale_IK1) .* scale_drug_IK1)
    GKr = (0.046 * scale_IKr) .* scale_drug_IKr
    GKs = (0.0034 * scale_IKs) .* scale_drug_IKs
    GNaL = scale_HF_GNaL .* ((0.0075 * scale_INaL) .* scale_drug_INaL)
    km2n = 1.0 * jca
    IpCa = (cai .* (GpCa .* scale_drug_IpCa)) ./ (cai + 0.0005)
    Istim = ((duration >= t) ? (amp) : (0))
    JdiffK = (-ki + kss) / 2.0
    JdiffNa = (-nai + nass) / 2.0
    Jtr = (-cajsr + cansr) / 100.0
    Jleak = ((0.0039375 * cansr) .* scale_HF_Jleak) / 15.0
    Knai = Knai0 .* exp((F .* (delta .* v)) ./ (((3.0 * R) .* T)))
    Knao = Knao0 .* exp((F .* (v .* (1.0 - delta))) ./ (((3.0 * R) .* T)))
    P = eP ./ (((H ./ Khp + 1.0) + nai ./ Knap) + ki ./ Kxkur)
    PCa = (0.0001 * scale_ICaL) .* scale_drug_ICaL
    XS_max = ((XS > 0) ? (XS) : (0))
    XW_max = ((XW > 0) ? (XW) : (0))
    XU = -XW + (-XS + (1 - TmB))
    a2 = k2p
    a4 = ((MgATP .* k4p) ./ Kmgatp) ./ (1.0 + MgATP ./ Kmgatp)
    a_rel = 0.5 * bt
    btp = 1.25 * bt
    tau_rel_tmp = bt ./ (1.0 + 0.0123 ./ cajsr)
    allo_i = 1.0 ./ ((KmCaAct ./ cai) .^ 2.0 + 1.0)
    allo_ss = 1.0 ./ ((KmCaAct ./ cass) .^ 2.0 + 1.0)
    b1 = MgADP .* k1m
    ksu = (kws .* rw) .* (-1 / 1 + 1 ./ (1 * rs))
    cs = ((kws .* phi) .* (rw .* (1 - rs))) ./ rs
    cw = ((kuw .* phi) .* ((1 - rs) .* (1 - rw))) ./ ((rw .* (1 - rs)))
    kwu = kuw .* (-1 / 1 + 1 ./ (1 * rw)) - kws
    gammasu = gammas .* (((Zetas >= -1 || Zetas <= 0 || Zetas > -Zetas - 1) && (Zetas > 0 && Zetas < -1 || (Zetas >= -1 || Zetas > -1) && Zetas < -1 || Zetas > 0)) ? (Zetas .* ((Zetas > 0) ? (1) : (0))) : ((-Zetas - 1 / 1) .* ((Zetas < -1 / 1) ? (1) : (0))))
    gammawu = gammaw .* abs(Zetaw)
    h10 = (nao ./ kna1) .* (1 + nao ./ kna2) + (kasymm + 1.0)
    h10_i = (nao ./ kna1) .* (1.0 + nao ./ kna2) + (kasymm + 1.0)
    h4 = (nass ./ kna1) .* (1 + nass ./ kna2) + 1.0
    h4_i = (nai ./ kna1) .* (1 + nai ./ kna2) + 1.0
    hca = exp((F .* (qca .* v)) ./ ((R .* T)))
    hna = exp((F .* (qna .* v)) ./ ((R .* T)))
    k2 = kcaoff
    k2_i = kcaoff
    k5 = kcaoff
    k5_i = kcaoff
    kb = (Trpn50 .^ ntm .* ku) ./ (-rw .* (1 - rs) + (1 - rs))
    lambda_min12 = ((lmbda < 1.2) ? (lmbda) : (1.2))
    thLp = scale_HF_thL .* (3.0 * thL)
    tiF = delta_epi .* (1 ./ (1 * (0.3933 * exp((-(v + 100.0)) / 100.0) + 0.08004 * exp((v + 50.0) / 16.59)))) + 4.562
    tiS = delta_epi .* (1 ./ (1 * (0.001416 * exp((-(v + 96.52)) / 59.05) + 1.78e-08 * exp((v + 114.1) / 8.079)))) + 23.62
    Afcas = 1.0 - Afcaf
    AiS = 1.0 - AiF
    Axrs = 1.0 - Axrf
    fcass = fss
    dhL_dt = (-hL + hLss) ./ ((scale_HF_thL .* thL))
    values[1] = dhL_dt
    jss = hss
    da_dt = (-a + ass) ./ ta
    values[2] = da_dt
    dap_dt = (-ap + assp) ./ ta
    values[3] = dap_dt
    dd_dt = (-d + dss) ./ td
    values[4] = dd_dt
    tfcafp = 2.5 * tfcaf
    tffp = 2.5 * tff
    dff_dt = (-ff + fss) ./ tff
    values[5] = dff_dt
    dfs_dt = (-fs + fss) ./ tfs
    values[6] = dfs_dt
    dhf_dt = (-hf + hss) ./ thf
    values[7] = dhf_dt
    thsp = 3.0 * ths
    dhs_dt = (-hs + hss) ./ ths
    values[8] = dhs_dt
    tjp = 1.46 * tj
    tmL = tm
    dm_dt = (-m + mss) ./ tm
    values[9] = dm_dt
    dxrf_dt = (-xrf + xrss) ./ txrf
    values[10] = dxrf_dt
    dxrs_dt = (-xrs + xrss) ./ txrs
    values[11] = dxrs_dt
    xs2ss = xs1ss
    dxs1_dt = (-xs1 + xs1ss) ./ txs1
    values[12] = dxs1_dt
    f = Aff .* ff + Afs .* fs
    fp = Aff .* ffp + Afs .* fs
    Acap = 2 * Ageo
    vjsr = 0.0048 * vcell
    vmyo = 0.68 * vcell
    vnsr = 0.0552 * vcell
    vss = 0.02 * vcell
    h = Ahf .* hf + Ahs .* hs
    hp = Ahf .* hf + Ahs .* hsp
    As = Aw
    CaMKa = scale_HF_CaMKa .* (CaMKb + CaMKt)
    dCaMKt_dt = -CaMKt .* bCaMK + (CaMKb .* aCaMK) .* (CaMKb + CaMKt)
    values[13] = dCaMKt_dt
    ICab = ((abs(v) < 0.019306977288832503) ? (0.6666666666666667 * F .^ 3 .* PCab .* cai .* scale_drug_ICab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) - 0.2273333333333334 * F .^ 3 .* PCab .* cao .* scale_drug_ICab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 2.0 * F .^ 2 .* PCab .* cai .* scale_drug_ICab .* v ./ (R .* T) + 0.682 * F .^ 2 .* PCab .* cao .* scale_drug_ICab .* v ./ (R .* T) + 2.0 * F .* PCab .* cai .* scale_drug_ICab - 0.682 * F .* PCab .* cao .* scale_drug_ICab) : ((((4.0 * (PCab .* scale_drug_ICab)) .* ((F .* (F .* v)) ./ ((R .* T)))) .* (cai .* exp(2.0 * ((F .* v) ./ ((R .* T)))) - 0.341 * cao)) ./ (exp(2.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    INab = ((abs(v) < 0.05179474679231207) ? (0.08333333333333334 * F .^ 3 .* PNab .* nai .* scale_drug_INab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) - 0.08333333333333334 * F .^ 3 .* PNab .* nao .* scale_drug_INab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.5 * F .^ 2 .* PNab .* nai .* scale_drug_INab .* v ./ (R .* T) + 0.5 * F .^ 2 .* PNab .* nao .* scale_drug_INab .* v ./ (R .* T) + F .* PNab .* nai .* scale_drug_INab - F .* PNab .* nao .* scale_drug_INab) : ((((PNab .* scale_drug_INab) .* ((F .* (F .* v)) ./ ((R .* T)))) .* (nai .* exp((F .* v) ./ ((R .* T))) - nao)) ./ (exp((F .* v) ./ ((R .* T))) - 1.0)))
    PhiCaK = ((abs(v) < 0.05179474679231207) ? (-0.0625 * F .^ 3 .* ko .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.0625 * F .^ 3 .* kss .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.375 * F .^ 2 .* ko .* v ./ (R .* T) + 0.375 * F .^ 2 .* kss .* v ./ (R .* T) - 0.75 * F .* ko + 0.75 * F .* kss) : (((1.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.75 * ko + (0.75 * kss) .* exp(1.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(1.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    PhiCaL = ((abs(v) < 0.019306977288832503) ? (-0.2273333333333334 * F .^ 3 .* cao .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.6666666666666667 * F .^ 3 .* cass .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.682 * F .^ 2 .* cao .* v ./ (R .* T) + 2.0 * F .^ 2 .* cass .* v ./ (R .* T) - 0.682 * F .* cao + 2.0 * F .* cass) : (((4.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.341 * cao + cass .* exp(2.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(2.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    PhiCaNa = ((abs(v) < 0.05179474679231207) ? (-0.0625 * F .^ 3 .* nao .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.0625 * F .^ 3 .* nass .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.375 * F .^ 2 .* nao .* v ./ (R .* T) + 0.375 * F .^ 2 .* nass .* v ./ (R .* T) - 0.75 * F .* nao + 0.75 * F .* nass) : (((1.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.75 * nao + (0.75 * nass) .* exp(1.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(1.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    IKb = (xkb .* (GKb .* scale_drug_IKb)) .* (-EK + v)
    dxk1_dt = (-xk1 + xk1ss) ./ txk1
    values[14] = dxk1_dt
    IK1 = (xk1 .* (rk1 .* (GK1 .* sqrt(ko)))) .* (-EK + v)
    IKs = (xs2 .* (xs1 .* (GKs .* KsCa))) .* (-EKs + v)
    anca = 1.0 ./ (k2n ./ km2n + (Kmn ./ cass + 1.0) .^ 4.0)
    a1 = (k1p .* (nai ./ Knai) .^ 3.0) ./ (((1.0 + ki ./ Kki) .^ 2.0 + (1.0 + nai ./ Knai) .^ 3.0) - 1.0)
    b4 = (k4m .* (ki ./ Kki) .^ 2.0) ./ (((1.0 + ki ./ Kki) .^ 2.0 + (1.0 + nai ./ Knai) .^ 3.0) - 1.0)
    a3 = (k3p .* (ko ./ Kko) .^ 2.0) ./ (((1.0 + ko ./ Kko) .^ 2.0 + (1.0 + nao ./ Knao) .^ 3.0) - 1.0)
    b2 = (k2m .* (nao ./ Knao) .^ 3.0) ./ (((1.0 + ko ./ Kko) .^ 2.0 + (1.0 + nao ./ Knao) .^ 3.0) - 1.0)
    b3 = (H .* (P .* k3m)) ./ (1.0 + MgATP ./ Kmgatp)
    PCaK = 0.0003574 * PCa
    PCaNa = 0.00125 * PCa
    PCap = 1.1 * PCa
    a_relp = 0.5 * btp
    tau_relp_tmp = btp ./ (1.0 + 0.0123 ./ cajsr)
    tau_rel = ((tau_rel_tmp < 0.001) ? (0.001) : (tau_rel_tmp))
    dZetaw_dt = Aw .* dLambda - Zetaw .* cw
    values[15] = dZetaw_dt
    dXS_dt = -XS .* gammasu + (-XS .* ksu + XW .* kws)
    values[16] = dXS_dt
    dXW_dt = -XW .* gammawu + (-XW .* kws + (XU .* kuw - XW .* kwu))
    values[17] = dXW_dt
    h11 = (nao .* nao) ./ ((kna2 .* (h10 .* kna1)))
    h12 = 1.0 ./ h10
    h11_i = (nao .* nao) ./ ((kna2 .* (h10_i .* kna1)))
    h12_i = 1.0 ./ h10_i
    h5 = (nass .* nass) ./ ((kna2 .* (h4 .* kna1)))
    h6 = 1.0 ./ h4
    h5_i = (nai .* nai) ./ ((kna2 .* (h4_i .* kna1)))
    h6_i = 1.0 ./ h4_i
    h1 = (nass ./ kna3) .* (hna + 1) + 1
    h1_i = (nai ./ kna3) .* (hna + 1) + 1
    h7 = (nao ./ kna3) .* (1.0 + 1.0 ./ hna) + 1.0
    h7_i = (nao ./ kna3) .* (1.0 + 1.0 ./ hna) + 1.0
    dTmB_dt = -TmB .* CaTrpn .^ (ntm / 2) .* ku + XU .* (kb .* ((CaTrpn .^ ((-ntm) / 2) < 100) ? (CaTrpn .^ ((-ntm) / 2)) : (100)))
    values[18] = dTmB_dt
    C = lambda_min12 - 1 / 1
    cat50 = scale_HF_cat50_ref .* (Beta1 .* (lambda_min12 - 1 / 1) + cat50_ref)
    lambda_min087 = ((lambda_min12 < 0.87) ? (lambda_min12) : (0.87))
    dhLp_dt = (-hLp + hLssp) ./ thLp
    values[19] = dhLp_dt
    tiFp = tiF .* (dti_develop .* dti_recover)
    diF_dt = (-iF + iss) ./ tiF
    values[20] = diF_dt
    tiSp = tiS .* (dti_develop .* dti_recover)
    diS_dt = (-iS + iss) ./ tiS
    values[21] = diS_dt
    fca = Afcaf .* fcaf + Afcas .* fcas
    fcap = Afcaf .* fcafp + Afcas .* fcas
    i = AiF .* iF + AiS .* iS
    ip = AiF .* iFp + AiS .* iSp
    xr = Axrf .* xrf + Axrs .* xrs
    dfcaf_dt = (-fcaf + fcass) ./ tfcaf
    values[22] = dfcaf_dt
    dfcas_dt = (-fcas + fcass) ./ tfcas
    values[23] = dfcas_dt
    djca_dt = (fcass - jca) ./ tjca
    values[24] = djca_dt
    dj_dt = (-j + jss) ./ tj
    values[25] = dj_dt
    dfcafp_dt = (-fcafp + fcass) ./ tfcafp
    values[26] = dfcafp_dt
    dffp_dt = (-ffp + fss) ./ tffp
    values[27] = dffp_dt
    dhsp_dt = (-hsp + hssp) ./ thsp
    values[28] = dhsp_dt
    djp_dt = (-jp + jss) ./ tjp
    values[29] = djp_dt
    dmL_dt = (-mL + mLss) ./ tmL
    values[30] = dmL_dt
    dxs2_dt = (-xs2 + xs2ss) ./ txs2
    values[31] = dxs2_dt
    dZetas_dt = As .* dLambda - Zetas .* cs
    values[32] = dZetas_dt
    fICaLp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    fINaLp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    fINap = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    fItop = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    fJrelp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    fJupp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    dnca_dt = anca .* k2n - km2n .* nca
    values[33] = dnca_dt
    x2 = b4 .* (a2 .* a3) + (b4 .* (a3 .* b1) + (a3 .* (a1 .* a2) + b4 .* (b1 .* b2)))
    x1 = a2 .* (a1 .* b3) + (b3 .* (a2 .* b4) + (a2 .* (a1 .* a4) + b3 .* (b2 .* b4)))
    x3 = b1 .* (a3 .* a4) + (a4 .* (b1 .* b2) + (a4 .* (a2 .* a3) + b1 .* (b2 .* b3)))
    x4 = a1 .* (b2 .* b3) + (a1 .* (a4 .* b2) + (a1 .* (a3 .* a4) + b2 .* (b3 .* b4)))
    PCaKp = 0.0003574 * PCap
    PCaNap = 0.00125 * PCap
    tau_relp = ((tau_relp_tmp < 0.001) ? (0.001) : (tau_relp_tmp))
    k1 = kcaon .* (cao .* h12)
    k1_i = kcaon .* (cao .* h12_i)
    k6 = kcaon .* (cass .* h6)
    k6_i = kcaon .* (cai .* h6_i)
    h2 = (hna .* nass) ./ ((h1 .* kna3))
    h3 = 1.0 ./ h1
    h2_i = (hna .* nai) ./ ((h1_i .* kna3))
    h3_i = 1.0 ./ h1_i
    h8 = nao ./ ((h7 .* (hna .* kna3)))
    h9 = 1.0 ./ h7
    h8_i = nao ./ ((h7_i .* (hna .* kna3)))
    h9_i = 1.0 ./ h7_i
    F1 = exp(C .* p_b) - 1 / 1
    dCd = C - Cd
    dCaTrpn_dt = ktrpn .* (-CaTrpn + ((1000 * cai) ./ cat50) .^ ntrpn .* (1 - CaTrpn))
    values[34] = dCaTrpn_dt
    h_lambda_prima = Beta0 .* ((lambda_min087 + lambda_min12) - 1.87) + 1
    diFp_dt = (-iFp + iss) ./ tiFp
    values[35] = diFp_dt
    diSp_dt = (-iSp + iss) ./ tiSp
    values[36] = diSp_dt
    IKr = (rkr .* (xr .* (GKr .* (0.4303314829119352 * sqrt(ko))))) .* (-EK + v)
    ICaL = (d .* (PhiCaL .* (PCa .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaL .* (PCap .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    INaL = (mL .* (GNaL .* (-ENa + v))) .* (fINaLp .* hLp + hL .* (1.0 - fINaLp))
    INa = (m .^ 3.0 .* ((GNa .* scale_drug_INa) .* (-ENa + v))) .* (j .* (h .* (1.0 - fINap)) + jp .* (fINap .* hp))
    Ito = ((scale_HF_Gto .* (Gto .* scale_drug_Ito)) .* (-EK + v)) .* (i .* (a .* (1.0 - fItop)) + ip .* (ap .* fItop))
    Jrel = Jrelnp .* (1.0 - fJrelp) + Jrelp .* fJrelp
    Jup = -Jleak + (Jupnp .* (1.0 - fJupp) + scale_HF_Jup .* (Jupp .* fJupp))
    E1 = x1 ./ (x4 + (x3 + (x1 + x2)))
    E2 = x2 ./ (x4 + (x3 + (x1 + x2)))
    E3 = x3 ./ (x4 + (x3 + (x1 + x2)))
    E4 = x4 ./ (x4 + (x3 + (x1 + x2)))
    ICaK = (d .* (PhiCaK .* (PCaK .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaK .* (PCaKp .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    ICaNa = (d .* (PhiCaNa .* (PCaNa .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaNa .* (PCaNap .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    k4pp = h2 .* wnaca
    k7 = wna .* (h2 .* h5)
    k4p_ss = (h3 .* wca) ./ hca
    k4pp_i = h2_i .* wnaca
    k7_i = wna .* (h2_i .* h5_i)
    k4p_i = (h3_i .* wca) ./ hca
    k3pp = h8 .* wnaca
    k8 = wna .* (h11 .* h8)
    k3p_ss = h9 .* wca
    k3pp_i = h8_i .* wnaca
    k8_i = wna .* (h11_i .* h8_i)
    k3p_i = h9_i .* wca
    eta = ((dCd < 0) ? (etas) : (etal))
    J_TRPN = dCaTrpn_dt .* trpnmax
    h_lambda = ((h_lambda_prima > 0) ? (h_lambda_prima) : (0))
    Jrel_inf = ((-ICaL) .* a_rel) ./ (((1.5 * scale_HF_Jrel_inf) ./ cajsr) .^ 8.0 + 1.0)
    Jrel_infp = ((-ICaL) .* a_relp) ./ (((1.5 * scale_HF_Jrel_inf) ./ cajsr) .^ 8.0 + 1.0)
    dcajsr_dt = Bcajsr .* (-Jrel + Jtr)
    values[37] = dcajsr_dt
    dcansr_dt = Jup - Jtr .* vjsr ./ vnsr
    values[38] = dcansr_dt
    JnakNa = 3.0 * (E1 .* a3 - E2 .* b3)
    JnakK = 2.0 * (-E3 .* a1 + E4 .* b1)
    dkss_dt = -JdiffK + (Acap .* (-ICaK)) ./ ((F .* vss))
    values[39] = dkss_dt
    k4 = k4p_ss + k4pp
    k4_i = k4p_i + k4pp_i
    k3 = k3p_ss + k3pp
    k3_i = k3p_i + k3pp_i
    Fd = dCd .* eta
    dCd_dt = (p_k .* (C - Cd)) ./ eta
    values[40] = dCd_dt
    Ta = (h_lambda .* (Tref ./ rs)) .* (XS .* (Zetas + 1) + XW .* Zetaw)
    dJrelnp_dt = (Jrel_inf - Jrelnp) ./ tau_rel
    values[41] = dJrelnp_dt
    dJrelp_dt = (Jrel_infp - Jrelp) ./ tau_relp
    values[42] = dJrelp_dt
    INaK = (Pnak .* scale_HF_Pnak) .* (JnakK .* zk + JnakNa .* zna)
    x2_ss = (k1 .* k7) .* (k4 + k5) + (k4 .* k6) .* (k1 + k8)
    x2_i = (k1_i .* k7_i) .* (k4_i + k5_i) + (k4_i .* k6_i) .* (k1_i + k8_i)
    x1_ss = (k2 .* k4) .* (k6 + k7) + (k5 .* k7) .* (k2 + k3)
    x3_ss = (k1 .* k3) .* (k6 + k7) + (k6 .* k8) .* (k2 + k3)
    x4_ss = (k2 .* k8) .* (k4 + k5) + (k3 .* k5) .* (k1 + k8)
    x1_i = (k2_i .* k4_i) .* (k6_i + k7_i) + (k5_i .* k7_i) .* (k2_i + k3_i)
    x3_i = (k1_i .* k3_i) .* (k6_i + k7_i) + (k6_i .* k8_i) .* (k2_i + k3_i)
    x4_i = (k2_i .* k8_i) .* (k4_i + k5_i) + (k3_i .* k5_i) .* (k1_i + k8_i)
    Tp = p_a .* (F1 + Fd)
    dki_dt = (Acap .* (-(-2.0 * INaK + (Istim + (Isac_P_ns / 3 + (Isac_P_k + (IKb + (IK1 + (IKs + (IKr + Ito)))))))))) ./ ((F .* vmyo)) + (JdiffK .* vss) ./ vmyo
    values[43] = dki_dt
    E1_ss = x1_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E2_ss = x2_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E3_ss = x3_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E4_ss = x4_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    E1_i = x1_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E2_i = x2_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E3_i = x3_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    E4_i = x4_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    Ttot = Ta + Tp
    JncxCa_ss = -E1_ss .* k1 + E2_ss .* k2
    JncxNa_ss = -E2_ss .* k3pp + (E3_ss .* k4pp + 3.0 * (-E1_ss .* k8 + E4_ss .* k7))
    JncxCa_i = -E1_i .* k1_i + E2_i .* k2_i
    JncxNa_i = -E2_i .* k3pp_i + (E3_i .* k4pp_i + 3.0 * (-E1_i .* k8_i + E4_i .* k7_i))
    INaCa_ss = (allo_ss .* ((0.2 * Gncx) .* scale_HF_Gncx)) .* (JncxCa_ss .* zca + JncxNa_ss .* zna)
    INaCa_i = (allo_i .* ((0.8 * Gncx) .* scale_HF_Gncx)) .* (JncxCa_i .* zca + JncxNa_i .* zna)
    dcass_dt = Bcass .* (-Jdiff + ((Acap .* (-(ICaL - 2.0 * INaCa_ss))) ./ (((2.0 * F) .* vss)) + (Jrel .* vjsr) ./ vss))
    values[44] = dcass_dt
    dnass_dt = -JdiffNa + (Acap .* (-(ICaNa + 3.0 * INaCa_ss))) ./ ((F .* vss))
    values[45] = dnass_dt
    dcai_dt = Bcai .* (-J_TRPN + (((Acap .* (-(Isac_P_ns / 3 + (-2.0 * INaCa_i + (ICab + IpCa))))) ./ (((2.0 * F) .* vmyo)) - Jup .* vnsr ./ vmyo) + (Jdiff .* vss) ./ vmyo))
    values[46] = dcai_dt
    dnai_dt = (Acap .* (-(Isac_P_ns / 3 + (INab + (3.0 * INaK + (3.0 * INaCa_i + (INa + INaL))))))) ./ ((F .* vmyo)) + (JdiffNa .* vss) ./ vmyo
    values[47] = dnai_dt
    dv_dt = -(Isac_P_k + (Isac_P_ns + (Istim + (ICab + (IpCa + (IKb + (INab + (INaK + (INaCa_ss + (INaCa_i + (IK1 + (IKs + (IKr + (ICaK + (ICaNa + (ICaL + (Ito + (INa + INaL))))))))))))))))))
    values[48] = dv_dt
end


function monitor_values(t, states, parameters, values)

    # Assign states
    hL = states[1]
    a = states[2]
    ap = states[3]
    d = states[4]
    ff = states[5]
    fs = states[6]
    hf = states[7]
    hs = states[8]
    m = states[9]
    xrf = states[10]
    xrs = states[11]
    xs1 = states[12]
    CaMKt = states[13]
    xk1 = states[14]
    Zetaw = states[15]
    XS = states[16]
    XW = states[17]
    TmB = states[18]
    hLp = states[19]
    iF = states[20]
    iS = states[21]
    fcaf = states[22]
    fcas = states[23]
    jca = states[24]
    j = states[25]
    fcafp = states[26]
    ffp = states[27]
    hsp = states[28]
    jp = states[29]
    mL = states[30]
    xs2 = states[31]
    Zetas = states[32]
    nca = states[33]
    CaTrpn = states[34]
    iFp = states[35]
    iSp = states[36]
    cajsr = states[37]
    cansr = states[38]
    kss = states[39]
    Cd = states[40]
    Jrelnp = states[41]
    Jrelp = states[42]
    ki = states[43]
    cass = states[44]
    nass = states[45]
    cai = states[46]
    nai = states[47]
    v = states[48]

    # Assign parameters
    Aff = parameters[1]
    Ahf = parameters[2]
    BSLmax = parameters[3]
    BSRmax = parameters[4]
    Beta0 = parameters[5]
    Beta1 = parameters[6]
    CaMKo = parameters[7]
    Esac_ns = parameters[8]
    F = parameters[9]
    GKb = parameters[10]
    GNa = parameters[11]
    Gncx = parameters[12]
    GpCa = parameters[13]
    Gsac_k = parameters[14]
    Gsac_ns = parameters[15]
    Gto = parameters[16]
    H = parameters[17]
    Khp = parameters[18]
    Kki = parameters[19]
    Kko = parameters[20]
    KmBSL = parameters[21]
    KmBSR = parameters[22]
    KmCaAct = parameters[23]
    KmCaM = parameters[24]
    KmCaMK = parameters[25]
    Kmgatp = parameters[26]
    Kmn = parameters[27]
    Knai0 = parameters[28]
    Knao0 = parameters[29]
    Knap = parameters[30]
    Kxkur = parameters[31]
    L = parameters[32]
    MgADP = parameters[33]
    MgATP = parameters[34]
    PCab = parameters[35]
    PKNa = parameters[36]
    PNab = parameters[37]
    Pnak = parameters[38]
    R = parameters[39]
    T = parameters[40]
    Tot_A = parameters[41]
    Tref = parameters[42]
    Trpn50 = parameters[43]
    aCaMK = parameters[44]
    amp = parameters[45]
    bCaMK = parameters[46]
    bt = parameters[47]
    calib = parameters[48]
    cao = parameters[49]
    cat50_ref = parameters[50]
    celltype = parameters[51]
    cmdnmax = parameters[52]
    csqnmax = parameters[53]
    dLambda = parameters[54]
    delta = parameters[55]
    delta_epi = parameters[56]
    duration = parameters[57]
    eP = parameters[58]
    emcoupling = parameters[59]
    etal = parameters[60]
    etas = parameters[61]
    gammas = parameters[62]
    gammaw = parameters[63]
    isacs = parameters[64]
    k1m = parameters[65]
    k1p = parameters[66]
    k2m = parameters[67]
    k2n = parameters[68]
    k2p = parameters[69]
    k3m = parameters[70]
    k3p = parameters[71]
    k4m = parameters[72]
    k4p = parameters[73]
    kasymm = parameters[74]
    kcaoff = parameters[75]
    kcaon = parameters[76]
    kmcmdn = parameters[77]
    kmcsqn = parameters[78]
    kmtrpn = parameters[79]
    kna1 = parameters[80]
    kna2 = parameters[81]
    kna3 = parameters[82]
    ko = parameters[83]
    ktrpn = parameters[84]
    ku = parameters[85]
    kuw = parameters[86]
    kws = parameters[87]
    lambda_max = parameters[88]
    lmbda = parameters[89]
    mode = parameters[90]
    nao = parameters[91]
    ntm = parameters[92]
    ntrpn = parameters[93]
    p_a = parameters[94]
    p_b = parameters[95]
    p_k = parameters[96]
    phi = parameters[97]
    qca = parameters[98]
    qna = parameters[99]
    rad = parameters[100]
    rs = parameters[101]
    rw = parameters[102]
    scale_HF_CaMKa = parameters[103]
    scale_HF_GK1 = parameters[104]
    scale_HF_GNaL = parameters[105]
    scale_HF_Gncx = parameters[106]
    scale_HF_Gto = parameters[107]
    scale_HF_Jleak = parameters[108]
    scale_HF_Jrel_inf = parameters[109]
    scale_HF_Jup = parameters[110]
    scale_HF_Pnak = parameters[111]
    scale_HF_cat50_ref = parameters[112]
    scale_HF_thL = parameters[113]
    scale_ICaL = parameters[114]
    scale_IK1 = parameters[115]
    scale_IKr = parameters[116]
    scale_IKs = parameters[117]
    scale_INaL = parameters[118]
    scale_drug_ICaL = parameters[119]
    scale_drug_ICab = parameters[120]
    scale_drug_IK1 = parameters[121]
    scale_drug_IKb = parameters[122]
    scale_drug_IKr = parameters[123]
    scale_drug_IKs = parameters[124]
    scale_drug_INa = parameters[125]
    scale_drug_INaL = parameters[126]
    scale_drug_INab = parameters[127]
    scale_drug_IpCa = parameters[128]
    scale_drug_Isack = parameters[129]
    scale_drug_Isacns = parameters[130]
    scale_drug_Ito = parameters[131]
    thL = parameters[132]
    tjca = parameters[133]
    trpnmax = parameters[134]
    wca = parameters[135]
    wna = parameters[136]
    wnaca = parameters[137]
    zca = parameters[138]
    zk = parameters[139]

    # Assign expressions
    zna = 1.0
    values[1] = zna
    Isac_P_k = 0
    values[2] = Isac_P_k
    Isac_P_ns = 0
    values[3] = Isac_P_ns
    Afcaf = 0.3 + 0.6 ./ (exp((v - 10.0) / 10.0) + 1.0)
    values[4] = Afcaf
    AiF = 1.0 ./ (exp((v - 213.6) / 151.2) + 1.0)
    values[5] = AiF
    Axrf = 1.0 ./ (exp((v + 54.81) / 38.21) + 1.0)
    values[6] = Axrf
    ass = 1.0 ./ (exp((-(v - 14.34)) / 14.82) + 1.0)
    values[7] = ass
    assp = 1.0 ./ (exp((-(v - 24.34)) / 14.82) + 1.0)
    values[8] = assp
    dss = 1.0 ./ (exp((-(v + 3.94)) / 4.23) + 1.0)
    values[9] = dss
    dti_develop = 1.354 + 0.0001 ./ (exp((-(v - 12.23)) / 0.2154) + exp((v - 167.4) / 15.89))
    values[10] = dti_develop
    dti_recover = 1.0 - 0.5 ./ (exp((v + 70.0) / 20.0) + 1.0)
    values[11] = dti_recover
    fss = 1.0 ./ (exp((v + 19.58) / 3.696) + 1.0)
    values[12] = fss
    hLss = 1.0 ./ (exp((v + 87.61) / 7.488) + 1.0)
    values[13] = hLss
    hLssp = 1.0 ./ (exp((v + 93.81) / 7.488) + 1.0)
    values[14] = hLssp
    hss = 1.0 ./ (exp((v + 78.5) / 6.22) + 1)
    values[15] = hss
    hssp = 1.0 ./ (exp(((v + 78.5) + 6.2) / 6.22) + 1)
    values[16] = hssp
    iss = 1.0 ./ (exp((v + 43.94) / 5.711) + 1.0)
    values[17] = iss
    mLss = 1.0 ./ (exp((-(v + 42.85)) / 5.264) + 1.0)
    values[18] = mLss
    mss = 1.0 ./ (exp((-((v + 39.57) + 9.4)) / 7.5) + 1.0)
    values[19] = mss
    rkr = (1.0 * (1.0 ./ (exp((v + 55.0) / 75.0) + 1.0))) ./ (exp((v - 10.0) / 30.0) + 1.0)
    values[20] = rkr
    ta = 1.0515 ./ (1.0 ./ ((1.2089 * (exp((-(v - 18.4099)) / 29.3814) + 1.0))) + 3.5 ./ (exp((v + 100.0) / 29.3814) + 1.0))
    values[21] = ta
    td = 0.6 + 1.0 ./ (exp((-0.05) * (v + 6.0)) + exp(0.09 * (v + 14.0)))
    values[22] = td
    tfcaf = 7.0 + 1.0 ./ (0.04 * exp((-(v - 4.0)) / 7.0) + 0.04 * exp((v - 4.0) / 7.0))
    values[23] = tfcaf
    tfcas = 100.0 + 1.0 ./ (0.00012 * exp((-v) / 3.0) + 0.00012 * exp(v / 7.0))
    values[24] = tfcas
    tff = 7.0 + 1.0 ./ (0.0045 * exp((-(v + 20.0)) / 10.0) + 0.0045 * exp((v + 20.0) / 10.0))
    values[25] = tff
    tfs = 1000.0 + 1.0 ./ (3.5e-05 * exp((-(v + 5.0)) / 4.0) + 3.5e-05 * exp((v + 5.0) / 6.0))
    values[26] = tfs
    thf = 1.0 ./ (6.149 * exp((v + 0.5096) / 20.27) + 1.432e-05 * exp((-(v + 1.196)) / 6.285))
    values[27] = thf
    ths = 1.0 ./ (0.009794 * exp((-(v + 17.95)) / 28.05) + 0.3343 * exp((v + 5.73) / 56.66))
    values[28] = ths
    tj = 2.038 + 1.0 ./ (0.3052 * exp((v + 0.9941) / 38.45) + 0.02136 * exp((-(v + 100.6)) / 8.281))
    values[29] = tj
    tm = 1.0 ./ (6.765 * exp((v + 11.64) / 34.77) + 8.552 * exp((-(v + 77.42)) / 5.955))
    values[30] = tm
    txk1 = 122.2 ./ (exp((-(v + 127.2)) / 20.36) + exp((v + 236.8) / 69.33))
    values[31] = txk1
    txrf = 12.98 + 1.0 ./ (4.123e-05 * exp((-(v - 47.78)) / 20.38) + 0.3652 * exp((v - 31.66) / 3.869))
    values[32] = txrf
    txrs = 1.865 + 1.0 ./ (1.128e-05 * exp((-(v - 29.74)) / 25.94) + 0.06629 * exp((v - 34.7) / 7.355))
    values[33] = txrs
    txs1 = 817.3 + 1.0 ./ (0.0002326 * exp((v + 48.28) / 17.8) + 0.001292 * exp((-(v + 210.0)) / 230.0))
    values[34] = txs1
    txs2 = 1.0 ./ (0.01 * exp((v - 50.0) / 20.0) + 0.0193 * exp((-(v + 66.54)) / 31.0))
    values[35] = txs2
    xkb = 1.0 ./ (exp((-(v - 14.48)) / 18.34) + 1.0)
    values[36] = xkb
    xrss = 1.0 ./ (exp((-(v + 8.337)) / 6.789) + 1.0)
    values[37] = xrss
    xs1ss = 1.0 ./ (exp((-(v + 11.6)) / 8.932) + 1.0)
    values[38] = xs1ss
    Afs = 1.0 - Aff
    values[39] = Afs
    Ageo = L .* ((2 * 3.14) * rad) + rad .* ((2 * 3.14) * rad)
    values[40] = Ageo
    vcell = L .* (rad .* ((3.14 * 1000) * rad))
    values[41] = vcell
    Ahs = 1.0 - Ahf
    values[42] = Ahs
    Aw = (Tot_A .* rs) ./ (rs + rw .* (1 - rs))
    values[43] = Aw
    Jupnp = (0.004375 * cai) ./ (cai + 0.00092)
    values[44] = Jupnp
    Jupp = ((0.004375 * 2.75) * cai) ./ ((cai + 0.00092) - 0.00017)
    values[45] = Jupp
    KsCa = 1.0 + 0.6 ./ ((3.8e-05 ./ cai) .^ 1.4 + 1.0)
    values[46] = KsCa
    Bcai = 1.0 ./ ((cmdnmax .* kmcmdn) ./ (cai + kmcmdn) .^ 2.0 + 1.0)
    values[47] = Bcai
    Bcajsr = 1.0 ./ ((csqnmax .* kmcsqn) ./ (cajsr + kmcsqn) .^ 2.0 + 1.0)
    values[48] = Bcajsr
    Bcass = 1.0 ./ ((BSLmax .* KmBSL) ./ (KmBSL + cass) .^ 2.0 + ((BSRmax .* KmBSR) ./ (KmBSR + cass) .^ 2.0 + 1.0))
    values[49] = Bcass
    Jdiff = (-cai + cass) / 0.2
    values[50] = Jdiff
    CaMKb = (CaMKo .* (1.0 - CaMKt)) ./ (KmCaM ./ cass + 1.0)
    values[51] = CaMKb
    CaTrpn_max = ((CaTrpn > 0) ? (CaTrpn) : (0))
    values[52] = CaTrpn_max
    vffrt = (F .* (F .* v)) ./ ((R .* T))
    values[53] = vffrt
    vfrt = (F .* v) ./ ((R .* T))
    values[54] = vfrt
    EK = ((R .* T) ./ F) .* log(ko ./ ki)
    values[55] = EK
    rk1 = 1.0 ./ (exp((-2.6 * ko + (v + 105.8)) / 9.493) + 1.0)
    values[56] = rk1
    xk1ss = 1.0 ./ (exp((-((2.5538 * ko + v) + 144.59)) ./ (1.5692 * ko + 3.8115)) + 1.0)
    values[57] = xk1ss
    EKs = ((R .* T) ./ F) .* log((PKNa .* nao + ko) ./ (PKNa .* nai + ki))
    values[58] = EKs
    ENa = ((R .* T) ./ F) .* log(nao ./ nai)
    values[59] = ENa
    GK1 = scale_HF_GK1 .* ((0.1908 * scale_IK1) .* scale_drug_IK1)
    values[60] = GK1
    GKr = (0.046 * scale_IKr) .* scale_drug_IKr
    values[61] = GKr
    GKs = (0.0034 * scale_IKs) .* scale_drug_IKs
    values[62] = GKs
    GNaL = scale_HF_GNaL .* ((0.0075 * scale_INaL) .* scale_drug_INaL)
    values[63] = GNaL
    km2n = 1.0 * jca
    values[64] = km2n
    IpCa = (cai .* (GpCa .* scale_drug_IpCa)) ./ (cai + 0.0005)
    values[65] = IpCa
    Istim = ((duration >= t) ? (amp) : (0))
    values[66] = Istim
    JdiffK = (-ki + kss) / 2.0
    values[67] = JdiffK
    JdiffNa = (-nai + nass) / 2.0
    values[68] = JdiffNa
    Jtr = (-cajsr + cansr) / 100.0
    values[69] = Jtr
    Jleak = ((0.0039375 * cansr) .* scale_HF_Jleak) / 15.0
    values[70] = Jleak
    Knai = Knai0 .* exp((F .* (delta .* v)) ./ (((3.0 * R) .* T)))
    values[71] = Knai
    Knao = Knao0 .* exp((F .* (v .* (1.0 - delta))) ./ (((3.0 * R) .* T)))
    values[72] = Knao
    P = eP ./ (((H ./ Khp + 1.0) + nai ./ Knap) + ki ./ Kxkur)
    values[73] = P
    PCa = (0.0001 * scale_ICaL) .* scale_drug_ICaL
    values[74] = PCa
    XS_max = ((XS > 0) ? (XS) : (0))
    values[75] = XS_max
    XW_max = ((XW > 0) ? (XW) : (0))
    values[76] = XW_max
    XU = -XW + (-XS + (1 - TmB))
    values[77] = XU
    a2 = k2p
    values[78] = a2
    a4 = ((MgATP .* k4p) ./ Kmgatp) ./ (1.0 + MgATP ./ Kmgatp)
    values[79] = a4
    a_rel = 0.5 * bt
    values[80] = a_rel
    btp = 1.25 * bt
    values[81] = btp
    tau_rel_tmp = bt ./ (1.0 + 0.0123 ./ cajsr)
    values[82] = tau_rel_tmp
    allo_i = 1.0 ./ ((KmCaAct ./ cai) .^ 2.0 + 1.0)
    values[83] = allo_i
    allo_ss = 1.0 ./ ((KmCaAct ./ cass) .^ 2.0 + 1.0)
    values[84] = allo_ss
    b1 = MgADP .* k1m
    values[85] = b1
    ksu = (kws .* rw) .* (-1 / 1 + 1 ./ (1 * rs))
    values[86] = ksu
    cs = ((kws .* phi) .* (rw .* (1 - rs))) ./ rs
    values[87] = cs
    cw = ((kuw .* phi) .* ((1 - rs) .* (1 - rw))) ./ ((rw .* (1 - rs)))
    values[88] = cw
    kwu = kuw .* (-1 / 1 + 1 ./ (1 * rw)) - kws
    values[89] = kwu
    gammasu = gammas .* (((Zetas >= -1 || Zetas <= 0 || Zetas > -Zetas - 1) && (Zetas > 0 && Zetas < -1 || (Zetas >= -1 || Zetas > -1) && Zetas < -1 || Zetas > 0)) ? (Zetas .* ((Zetas > 0) ? (1) : (0))) : ((-Zetas - 1 / 1) .* ((Zetas < -1 / 1) ? (1) : (0))))
    values[90] = gammasu
    gammawu = gammaw .* abs(Zetaw)
    values[91] = gammawu
    h10 = (nao ./ kna1) .* (1 + nao ./ kna2) + (kasymm + 1.0)
    values[92] = h10
    h10_i = (nao ./ kna1) .* (1.0 + nao ./ kna2) + (kasymm + 1.0)
    values[93] = h10_i
    h4 = (nass ./ kna1) .* (1 + nass ./ kna2) + 1.0
    values[94] = h4
    h4_i = (nai ./ kna1) .* (1 + nai ./ kna2) + 1.0
    values[95] = h4_i
    hca = exp((F .* (qca .* v)) ./ ((R .* T)))
    values[96] = hca
    hna = exp((F .* (qna .* v)) ./ ((R .* T)))
    values[97] = hna
    k2 = kcaoff
    values[98] = k2
    k2_i = kcaoff
    values[99] = k2_i
    k5 = kcaoff
    values[100] = k5
    k5_i = kcaoff
    values[101] = k5_i
    kb = (Trpn50 .^ ntm .* ku) ./ (-rw .* (1 - rs) + (1 - rs))
    values[102] = kb
    lambda_min12 = ((lmbda < 1.2) ? (lmbda) : (1.2))
    values[103] = lambda_min12
    thLp = scale_HF_thL .* (3.0 * thL)
    values[104] = thLp
    tiF = delta_epi .* (1 ./ (1 * (0.3933 * exp((-(v + 100.0)) / 100.0) + 0.08004 * exp((v + 50.0) / 16.59)))) + 4.562
    values[105] = tiF
    tiS = delta_epi .* (1 ./ (1 * (0.001416 * exp((-(v + 96.52)) / 59.05) + 1.78e-08 * exp((v + 114.1) / 8.079)))) + 23.62
    values[106] = tiS
    Afcas = 1.0 - Afcaf
    values[107] = Afcas
    AiS = 1.0 - AiF
    values[108] = AiS
    Axrs = 1.0 - Axrf
    values[109] = Axrs
    fcass = fss
    values[110] = fcass
    dhL_dt = (-hL + hLss) ./ ((scale_HF_thL .* thL))
    values[111] = dhL_dt
    jss = hss
    values[112] = jss
    da_dt = (-a + ass) ./ ta
    values[113] = da_dt
    dap_dt = (-ap + assp) ./ ta
    values[114] = dap_dt
    dd_dt = (-d + dss) ./ td
    values[115] = dd_dt
    tfcafp = 2.5 * tfcaf
    values[116] = tfcafp
    tffp = 2.5 * tff
    values[117] = tffp
    dff_dt = (-ff + fss) ./ tff
    values[118] = dff_dt
    dfs_dt = (-fs + fss) ./ tfs
    values[119] = dfs_dt
    dhf_dt = (-hf + hss) ./ thf
    values[120] = dhf_dt
    thsp = 3.0 * ths
    values[121] = thsp
    dhs_dt = (-hs + hss) ./ ths
    values[122] = dhs_dt
    tjp = 1.46 * tj
    values[123] = tjp
    tmL = tm
    values[124] = tmL
    dm_dt = (-m + mss) ./ tm
    values[125] = dm_dt
    dxrf_dt = (-xrf + xrss) ./ txrf
    values[126] = dxrf_dt
    dxrs_dt = (-xrs + xrss) ./ txrs
    values[127] = dxrs_dt
    xs2ss = xs1ss
    values[128] = xs2ss
    dxs1_dt = (-xs1 + xs1ss) ./ txs1
    values[129] = dxs1_dt
    f = Aff .* ff + Afs .* fs
    values[130] = f
    fp = Aff .* ffp + Afs .* fs
    values[131] = fp
    Acap = 2 * Ageo
    values[132] = Acap
    vjsr = 0.0048 * vcell
    values[133] = vjsr
    vmyo = 0.68 * vcell
    values[134] = vmyo
    vnsr = 0.0552 * vcell
    values[135] = vnsr
    vss = 0.02 * vcell
    values[136] = vss
    h = Ahf .* hf + Ahs .* hs
    values[137] = h
    hp = Ahf .* hf + Ahs .* hsp
    values[138] = hp
    As = Aw
    values[139] = As
    CaMKa = scale_HF_CaMKa .* (CaMKb + CaMKt)
    values[140] = CaMKa
    dCaMKt_dt = -CaMKt .* bCaMK + (CaMKb .* aCaMK) .* (CaMKb + CaMKt)
    values[141] = dCaMKt_dt
    ICab = ((abs(v) < 0.019306977288832503) ? (0.6666666666666667 * F .^ 3 .* PCab .* cai .* scale_drug_ICab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) - 0.2273333333333334 * F .^ 3 .* PCab .* cao .* scale_drug_ICab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 2.0 * F .^ 2 .* PCab .* cai .* scale_drug_ICab .* v ./ (R .* T) + 0.682 * F .^ 2 .* PCab .* cao .* scale_drug_ICab .* v ./ (R .* T) + 2.0 * F .* PCab .* cai .* scale_drug_ICab - 0.682 * F .* PCab .* cao .* scale_drug_ICab) : ((((4.0 * (PCab .* scale_drug_ICab)) .* ((F .* (F .* v)) ./ ((R .* T)))) .* (cai .* exp(2.0 * ((F .* v) ./ ((R .* T)))) - 0.341 * cao)) ./ (exp(2.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    values[142] = ICab
    INab = ((abs(v) < 0.05179474679231207) ? (0.08333333333333334 * F .^ 3 .* PNab .* nai .* scale_drug_INab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) - 0.08333333333333334 * F .^ 3 .* PNab .* nao .* scale_drug_INab .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.5 * F .^ 2 .* PNab .* nai .* scale_drug_INab .* v ./ (R .* T) + 0.5 * F .^ 2 .* PNab .* nao .* scale_drug_INab .* v ./ (R .* T) + F .* PNab .* nai .* scale_drug_INab - F .* PNab .* nao .* scale_drug_INab) : ((((PNab .* scale_drug_INab) .* ((F .* (F .* v)) ./ ((R .* T)))) .* (nai .* exp((F .* v) ./ ((R .* T))) - nao)) ./ (exp((F .* v) ./ ((R .* T))) - 1.0)))
    values[143] = INab
    PhiCaK = ((abs(v) < 0.05179474679231207) ? (-0.0625 * F .^ 3 .* ko .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.0625 * F .^ 3 .* kss .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.375 * F .^ 2 .* ko .* v ./ (R .* T) + 0.375 * F .^ 2 .* kss .* v ./ (R .* T) - 0.75 * F .* ko + 0.75 * F .* kss) : (((1.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.75 * ko + (0.75 * kss) .* exp(1.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(1.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    values[144] = PhiCaK
    PhiCaL = ((abs(v) < 0.019306977288832503) ? (-0.2273333333333334 * F .^ 3 .* cao .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.6666666666666667 * F .^ 3 .* cass .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.682 * F .^ 2 .* cao .* v ./ (R .* T) + 2.0 * F .^ 2 .* cass .* v ./ (R .* T) - 0.682 * F .* cao + 2.0 * F .* cass) : (((4.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.341 * cao + cass .* exp(2.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(2.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    values[145] = PhiCaL
    PhiCaNa = ((abs(v) < 0.05179474679231207) ? (-0.0625 * F .^ 3 .* nao .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.0625 * F .^ 3 .* nass .* v .^ 2 ./ (R .^ 2 .* T .^ 2) + 0.375 * F .^ 2 .* nao .* v ./ (R .* T) + 0.375 * F .^ 2 .* nass .* v ./ (R .* T) - 0.75 * F .* nao + 0.75 * F .* nass) : (((1.0 * ((F .* (F .* v)) ./ ((R .* T)))) .* (-0.75 * nao + (0.75 * nass) .* exp(1.0 * ((F .* v) ./ ((R .* T)))))) ./ (exp(1.0 * ((F .* v) ./ ((R .* T)))) - 1.0)))
    values[146] = PhiCaNa
    IKb = (xkb .* (GKb .* scale_drug_IKb)) .* (-EK + v)
    values[147] = IKb
    dxk1_dt = (-xk1 + xk1ss) ./ txk1
    values[148] = dxk1_dt
    IK1 = (xk1 .* (rk1 .* (GK1 .* sqrt(ko)))) .* (-EK + v)
    values[149] = IK1
    IKs = (xs2 .* (xs1 .* (GKs .* KsCa))) .* (-EKs + v)
    values[150] = IKs
    anca = 1.0 ./ (k2n ./ km2n + (Kmn ./ cass + 1.0) .^ 4.0)
    values[151] = anca
    a1 = (k1p .* (nai ./ Knai) .^ 3.0) ./ (((1.0 + ki ./ Kki) .^ 2.0 + (1.0 + nai ./ Knai) .^ 3.0) - 1.0)
    values[152] = a1
    b4 = (k4m .* (ki ./ Kki) .^ 2.0) ./ (((1.0 + ki ./ Kki) .^ 2.0 + (1.0 + nai ./ Knai) .^ 3.0) - 1.0)
    values[153] = b4
    a3 = (k3p .* (ko ./ Kko) .^ 2.0) ./ (((1.0 + ko ./ Kko) .^ 2.0 + (1.0 + nao ./ Knao) .^ 3.0) - 1.0)
    values[154] = a3
    b2 = (k2m .* (nao ./ Knao) .^ 3.0) ./ (((1.0 + ko ./ Kko) .^ 2.0 + (1.0 + nao ./ Knao) .^ 3.0) - 1.0)
    values[155] = b2
    b3 = (H .* (P .* k3m)) ./ (1.0 + MgATP ./ Kmgatp)
    values[156] = b3
    PCaK = 0.0003574 * PCa
    values[157] = PCaK
    PCaNa = 0.00125 * PCa
    values[158] = PCaNa
    PCap = 1.1 * PCa
    values[159] = PCap
    a_relp = 0.5 * btp
    values[160] = a_relp
    tau_relp_tmp = btp ./ (1.0 + 0.0123 ./ cajsr)
    values[161] = tau_relp_tmp
    tau_rel = ((tau_rel_tmp < 0.001) ? (0.001) : (tau_rel_tmp))
    values[162] = tau_rel
    dZetaw_dt = Aw .* dLambda - Zetaw .* cw
    values[163] = dZetaw_dt
    dXS_dt = -XS .* gammasu + (-XS .* ksu + XW .* kws)
    values[164] = dXS_dt
    dXW_dt = -XW .* gammawu + (-XW .* kws + (XU .* kuw - XW .* kwu))
    values[165] = dXW_dt
    h11 = (nao .* nao) ./ ((kna2 .* (h10 .* kna1)))
    values[166] = h11
    h12 = 1.0 ./ h10
    values[167] = h12
    h11_i = (nao .* nao) ./ ((kna2 .* (h10_i .* kna1)))
    values[168] = h11_i
    h12_i = 1.0 ./ h10_i
    values[169] = h12_i
    h5 = (nass .* nass) ./ ((kna2 .* (h4 .* kna1)))
    values[170] = h5
    h6 = 1.0 ./ h4
    values[171] = h6
    h5_i = (nai .* nai) ./ ((kna2 .* (h4_i .* kna1)))
    values[172] = h5_i
    h6_i = 1.0 ./ h4_i
    values[173] = h6_i
    h1 = (nass ./ kna3) .* (hna + 1) + 1
    values[174] = h1
    h1_i = (nai ./ kna3) .* (hna + 1) + 1
    values[175] = h1_i
    h7 = (nao ./ kna3) .* (1.0 + 1.0 ./ hna) + 1.0
    values[176] = h7
    h7_i = (nao ./ kna3) .* (1.0 + 1.0 ./ hna) + 1.0
    values[177] = h7_i
    dTmB_dt = -TmB .* CaTrpn .^ (ntm / 2) .* ku + XU .* (kb .* ((CaTrpn .^ ((-ntm) / 2) < 100) ? (CaTrpn .^ ((-ntm) / 2)) : (100)))
    values[178] = dTmB_dt
    C = lambda_min12 - 1 / 1
    values[179] = C
    cat50 = scale_HF_cat50_ref .* (Beta1 .* (lambda_min12 - 1 / 1) + cat50_ref)
    values[180] = cat50
    lambda_min087 = ((lambda_min12 < 0.87) ? (lambda_min12) : (0.87))
    values[181] = lambda_min087
    dhLp_dt = (-hLp + hLssp) ./ thLp
    values[182] = dhLp_dt
    tiFp = tiF .* (dti_develop .* dti_recover)
    values[183] = tiFp
    diF_dt = (-iF + iss) ./ tiF
    values[184] = diF_dt
    tiSp = tiS .* (dti_develop .* dti_recover)
    values[185] = tiSp
    diS_dt = (-iS + iss) ./ tiS
    values[186] = diS_dt
    fca = Afcaf .* fcaf + Afcas .* fcas
    values[187] = fca
    fcap = Afcaf .* fcafp + Afcas .* fcas
    values[188] = fcap
    i = AiF .* iF + AiS .* iS
    values[189] = i
    ip = AiF .* iFp + AiS .* iSp
    values[190] = ip
    xr = Axrf .* xrf + Axrs .* xrs
    values[191] = xr
    dfcaf_dt = (-fcaf + fcass) ./ tfcaf
    values[192] = dfcaf_dt
    dfcas_dt = (-fcas + fcass) ./ tfcas
    values[193] = dfcas_dt
    djca_dt = (fcass - jca) ./ tjca
    values[194] = djca_dt
    dj_dt = (-j + jss) ./ tj
    values[195] = dj_dt
    dfcafp_dt = (-fcafp + fcass) ./ tfcafp
    values[196] = dfcafp_dt
    dffp_dt = (-ffp + fss) ./ tffp
    values[197] = dffp_dt
    dhsp_dt = (-hsp + hssp) ./ thsp
    values[198] = dhsp_dt
    djp_dt = (-jp + jss) ./ tjp
    values[199] = djp_dt
    dmL_dt = (-mL + mLss) ./ tmL
    values[200] = dmL_dt
    dxs2_dt = (-xs2 + xs2ss) ./ txs2
    values[201] = dxs2_dt
    dZetas_dt = As .* dLambda - Zetas .* cs
    values[202] = dZetas_dt
    fICaLp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[203] = fICaLp
    fINaLp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[204] = fINaLp
    fINap = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[205] = fINap
    fItop = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[206] = fItop
    fJrelp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[207] = fJrelp
    fJupp = 1.0 ./ (1.0 + KmCaMK ./ CaMKa)
    values[208] = fJupp
    dnca_dt = anca .* k2n - km2n .* nca
    values[209] = dnca_dt
    x2 = b4 .* (a2 .* a3) + (b4 .* (a3 .* b1) + (a3 .* (a1 .* a2) + b4 .* (b1 .* b2)))
    values[210] = x2
    x1 = a2 .* (a1 .* b3) + (b3 .* (a2 .* b4) + (a2 .* (a1 .* a4) + b3 .* (b2 .* b4)))
    values[211] = x1
    x3 = b1 .* (a3 .* a4) + (a4 .* (b1 .* b2) + (a4 .* (a2 .* a3) + b1 .* (b2 .* b3)))
    values[212] = x3
    x4 = a1 .* (b2 .* b3) + (a1 .* (a4 .* b2) + (a1 .* (a3 .* a4) + b2 .* (b3 .* b4)))
    values[213] = x4
    PCaKp = 0.0003574 * PCap
    values[214] = PCaKp
    PCaNap = 0.00125 * PCap
    values[215] = PCaNap
    tau_relp = ((tau_relp_tmp < 0.001) ? (0.001) : (tau_relp_tmp))
    values[216] = tau_relp
    k1 = kcaon .* (cao .* h12)
    values[217] = k1
    k1_i = kcaon .* (cao .* h12_i)
    values[218] = k1_i
    k6 = kcaon .* (cass .* h6)
    values[219] = k6
    k6_i = kcaon .* (cai .* h6_i)
    values[220] = k6_i
    h2 = (hna .* nass) ./ ((h1 .* kna3))
    values[221] = h2
    h3 = 1.0 ./ h1
    values[222] = h3
    h2_i = (hna .* nai) ./ ((h1_i .* kna3))
    values[223] = h2_i
    h3_i = 1.0 ./ h1_i
    values[224] = h3_i
    h8 = nao ./ ((h7 .* (hna .* kna3)))
    values[225] = h8
    h9 = 1.0 ./ h7
    values[226] = h9
    h8_i = nao ./ ((h7_i .* (hna .* kna3)))
    values[227] = h8_i
    h9_i = 1.0 ./ h7_i
    values[228] = h9_i
    F1 = exp(C .* p_b) - 1 / 1
    values[229] = F1
    dCd = C - Cd
    values[230] = dCd
    dCaTrpn_dt = ktrpn .* (-CaTrpn + ((1000 * cai) ./ cat50) .^ ntrpn .* (1 - CaTrpn))
    values[231] = dCaTrpn_dt
    h_lambda_prima = Beta0 .* ((lambda_min087 + lambda_min12) - 1.87) + 1
    values[232] = h_lambda_prima
    diFp_dt = (-iFp + iss) ./ tiFp
    values[233] = diFp_dt
    diSp_dt = (-iSp + iss) ./ tiSp
    values[234] = diSp_dt
    IKr = (rkr .* (xr .* (GKr .* (0.4303314829119352 * sqrt(ko))))) .* (-EK + v)
    values[235] = IKr
    ICaL = (d .* (PhiCaL .* (PCa .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaL .* (PCap .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    values[236] = ICaL
    INaL = (mL .* (GNaL .* (-ENa + v))) .* (fINaLp .* hLp + hL .* (1.0 - fINaLp))
    values[237] = INaL
    INa = (m .^ 3.0 .* ((GNa .* scale_drug_INa) .* (-ENa + v))) .* (j .* (h .* (1.0 - fINap)) + jp .* (fINap .* hp))
    values[238] = INa
    Ito = ((scale_HF_Gto .* (Gto .* scale_drug_Ito)) .* (-EK + v)) .* (i .* (a .* (1.0 - fItop)) + ip .* (ap .* fItop))
    values[239] = Ito
    Jrel = Jrelnp .* (1.0 - fJrelp) + Jrelp .* fJrelp
    values[240] = Jrel
    Jup = -Jleak + (Jupnp .* (1.0 - fJupp) + scale_HF_Jup .* (Jupp .* fJupp))
    values[241] = Jup
    E1 = x1 ./ (x4 + (x3 + (x1 + x2)))
    values[242] = E1
    E2 = x2 ./ (x4 + (x3 + (x1 + x2)))
    values[243] = E2
    E3 = x3 ./ (x4 + (x3 + (x1 + x2)))
    values[244] = E3
    E4 = x4 ./ (x4 + (x3 + (x1 + x2)))
    values[245] = E4
    ICaK = (d .* (PhiCaK .* (PCaK .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaK .* (PCaKp .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    values[246] = ICaK
    ICaNa = (d .* (PhiCaNa .* (PCaNa .* (1.0 - fICaLp)))) .* (f .* (1.0 - nca) + nca .* (fca .* jca)) + (d .* (PhiCaNa .* (PCaNap .* fICaLp))) .* (fp .* (1.0 - nca) + nca .* (fcap .* jca))
    values[247] = ICaNa
    k4pp = h2 .* wnaca
    values[248] = k4pp
    k7 = wna .* (h2 .* h5)
    values[249] = k7
    k4p_ss = (h3 .* wca) ./ hca
    values[250] = k4p_ss
    k4pp_i = h2_i .* wnaca
    values[251] = k4pp_i
    k7_i = wna .* (h2_i .* h5_i)
    values[252] = k7_i
    k4p_i = (h3_i .* wca) ./ hca
    values[253] = k4p_i
    k3pp = h8 .* wnaca
    values[254] = k3pp
    k8 = wna .* (h11 .* h8)
    values[255] = k8
    k3p_ss = h9 .* wca
    values[256] = k3p_ss
    k3pp_i = h8_i .* wnaca
    values[257] = k3pp_i
    k8_i = wna .* (h11_i .* h8_i)
    values[258] = k8_i
    k3p_i = h9_i .* wca
    values[259] = k3p_i
    eta = ((dCd < 0) ? (etas) : (etal))
    values[260] = eta
    J_TRPN = dCaTrpn_dt .* trpnmax
    values[261] = J_TRPN
    h_lambda = ((h_lambda_prima > 0) ? (h_lambda_prima) : (0))
    values[262] = h_lambda
    Jrel_inf = ((-ICaL) .* a_rel) ./ (((1.5 * scale_HF_Jrel_inf) ./ cajsr) .^ 8.0 + 1.0)
    values[263] = Jrel_inf
    Jrel_infp = ((-ICaL) .* a_relp) ./ (((1.5 * scale_HF_Jrel_inf) ./ cajsr) .^ 8.0 + 1.0)
    values[264] = Jrel_infp
    dcajsr_dt = Bcajsr .* (-Jrel + Jtr)
    values[265] = dcajsr_dt
    dcansr_dt = Jup - Jtr .* vjsr ./ vnsr
    values[266] = dcansr_dt
    JnakNa = 3.0 * (E1 .* a3 - E2 .* b3)
    values[267] = JnakNa
    JnakK = 2.0 * (-E3 .* a1 + E4 .* b1)
    values[268] = JnakK
    dkss_dt = -JdiffK + (Acap .* (-ICaK)) ./ ((F .* vss))
    values[269] = dkss_dt
    k4 = k4p_ss + k4pp
    values[270] = k4
    k4_i = k4p_i + k4pp_i
    values[271] = k4_i
    k3 = k3p_ss + k3pp
    values[272] = k3
    k3_i = k3p_i + k3pp_i
    values[273] = k3_i
    Fd = dCd .* eta
    values[274] = Fd
    dCd_dt = (p_k .* (C - Cd)) ./ eta
    values[275] = dCd_dt
    Ta = (h_lambda .* (Tref ./ rs)) .* (XS .* (Zetas + 1) + XW .* Zetaw)
    values[276] = Ta
    dJrelnp_dt = (Jrel_inf - Jrelnp) ./ tau_rel
    values[277] = dJrelnp_dt
    dJrelp_dt = (Jrel_infp - Jrelp) ./ tau_relp
    values[278] = dJrelp_dt
    INaK = (Pnak .* scale_HF_Pnak) .* (JnakK .* zk + JnakNa .* zna)
    values[279] = INaK
    x2_ss = (k1 .* k7) .* (k4 + k5) + (k4 .* k6) .* (k1 + k8)
    values[280] = x2_ss
    x2_i = (k1_i .* k7_i) .* (k4_i + k5_i) + (k4_i .* k6_i) .* (k1_i + k8_i)
    values[281] = x2_i
    x1_ss = (k2 .* k4) .* (k6 + k7) + (k5 .* k7) .* (k2 + k3)
    values[282] = x1_ss
    x3_ss = (k1 .* k3) .* (k6 + k7) + (k6 .* k8) .* (k2 + k3)
    values[283] = x3_ss
    x4_ss = (k2 .* k8) .* (k4 + k5) + (k3 .* k5) .* (k1 + k8)
    values[284] = x4_ss
    x1_i = (k2_i .* k4_i) .* (k6_i + k7_i) + (k5_i .* k7_i) .* (k2_i + k3_i)
    values[285] = x1_i
    x3_i = (k1_i .* k3_i) .* (k6_i + k7_i) + (k6_i .* k8_i) .* (k2_i + k3_i)
    values[286] = x3_i
    x4_i = (k2_i .* k8_i) .* (k4_i + k5_i) + (k3_i .* k5_i) .* (k1_i + k8_i)
    values[287] = x4_i
    Tp = p_a .* (F1 + Fd)
    values[288] = Tp
    dki_dt = (Acap .* (-(-2.0 * INaK + (Istim + (Isac_P_ns / 3 + (Isac_P_k + (IKb + (IK1 + (IKs + (IKr + Ito)))))))))) ./ ((F .* vmyo)) + (JdiffK .* vss) ./ vmyo
    values[289] = dki_dt
    E1_ss = x1_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    values[290] = E1_ss
    E2_ss = x2_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    values[291] = E2_ss
    E3_ss = x3_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    values[292] = E3_ss
    E4_ss = x4_ss ./ (x4_ss + (x3_ss + (x1_ss + x2_ss)))
    values[293] = E4_ss
    E1_i = x1_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    values[294] = E1_i
    E2_i = x2_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    values[295] = E2_i
    E3_i = x3_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    values[296] = E3_i
    E4_i = x4_i ./ (x4_i + (x3_i + (x1_i + x2_i)))
    values[297] = E4_i
    Ttot = Ta + Tp
    values[298] = Ttot
    JncxCa_ss = -E1_ss .* k1 + E2_ss .* k2
    values[299] = JncxCa_ss
    JncxNa_ss = -E2_ss .* k3pp + (E3_ss .* k4pp + 3.0 * (-E1_ss .* k8 + E4_ss .* k7))
    values[300] = JncxNa_ss
    JncxCa_i = -E1_i .* k1_i + E2_i .* k2_i
    values[301] = JncxCa_i
    JncxNa_i = -E2_i .* k3pp_i + (E3_i .* k4pp_i + 3.0 * (-E1_i .* k8_i + E4_i .* k7_i))
    values[302] = JncxNa_i
    INaCa_ss = (allo_ss .* ((0.2 * Gncx) .* scale_HF_Gncx)) .* (JncxCa_ss .* zca + JncxNa_ss .* zna)
    values[303] = INaCa_ss
    INaCa_i = (allo_i .* ((0.8 * Gncx) .* scale_HF_Gncx)) .* (JncxCa_i .* zca + JncxNa_i .* zna)
    values[304] = INaCa_i
    dcass_dt = Bcass .* (-Jdiff + ((Acap .* (-(ICaL - 2.0 * INaCa_ss))) ./ (((2.0 * F) .* vss)) + (Jrel .* vjsr) ./ vss))
    values[305] = dcass_dt
    dnass_dt = -JdiffNa + (Acap .* (-(ICaNa + 3.0 * INaCa_ss))) ./ ((F .* vss))
    values[306] = dnass_dt
    dcai_dt = Bcai .* (-J_TRPN + (((Acap .* (-(Isac_P_ns / 3 + (-2.0 * INaCa_i + (ICab + IpCa))))) ./ (((2.0 * F) .* vmyo)) - Jup .* vnsr ./ vmyo) + (Jdiff .* vss) ./ vmyo))
    values[307] = dcai_dt
    dnai_dt = (Acap .* (-(Isac_P_ns / 3 + (INab + (3.0 * INaK + (3.0 * INaCa_i + (INa + INaL))))))) ./ ((F .* vmyo)) + (JdiffNa .* vss) ./ vmyo
    values[308] = dnai_dt
    dv_dt = -(Isac_P_k + (Isac_P_ns + (Istim + (ICab + (IpCa + (IKb + (INab + (INaK + (INaCa_ss + (INaCa_i + (IK1 + (IKs + (IKr + (ICaK + (ICaNa + (ICaL + (Ito + (INa + INaL))))))))))))))))))
    values[309] = dv_dt
end

