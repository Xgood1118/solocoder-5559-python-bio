import re

PROSITE_PATTERNS = {
    'PS00001': {'name': 'ATP_GTP_A', 'pattern': '[AG]-x(4)-G-K-[ST]', 'description': 'ATP/GTP-binding site motif A (P-loop)'},
    'PS00002': {'name': 'ATP_GTP_B', 'pattern': '[DE]-x(2)-[LIVMFY]-x(3)-R', 'description': 'ATP/GTP-binding site motif B'},
    'PS00003': {'name': 'ZINC_FINGER_C2H2', 'pattern': 'C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H', 'description': 'Zinc finger C2H2 type'},
    'PS00004': {'name': 'ZINC_FINGER_C4', 'pattern': 'C-x(2)-C-x(13)-C-x(2)-C-x(14,52)-C-x(2)-C-x(10,16)-C-x(2)-C', 'description': 'Zinc finger C4 type (two domains)'},
    'PS00005': {'name': 'ZINC_FINGER_C6', 'pattern': 'C-x(2)-C-x(6)-C-x(6)-C-x(2)-C-x(6)-C', 'description': 'Zinc finger C6 type (GAL4 family)'},
    'PS00006': {'name': 'ZINC_RING', 'pattern': 'C-x-H-[LIVMFY]-x-C-x-[FY]-x(2)-C-x(2)-C-x-[STAGN]-x(5,15)-[P]-x-[LIVMFY]-x-C-x-[LIVMFY]-C-x(2)-C', 'description': 'Zinc ring finger'},
    'PS00007': {'name': 'EGF_CA', 'pattern': 'C-x(3,5)-S-x(6,12)-C-x(2)-G-x(2,3)-C-x(2,7)-C', 'description': 'Calcium-binding EGF domain'},
    'PS00008': {'name': 'EGF', 'pattern': 'C-x-C-x(5,7)-[FY]-x-C-x-[GADENQKRSTVIL]-x-[FYHV]-x-C', 'description': 'EGF-like domain signature 1'},
    'PS00009': {'name': 'RGD', 'pattern': 'R-G-D', 'description': 'RGD cell attachment sequence'},
    'PS00010': {'name': 'RNP_1', 'pattern': '[RK]-G-{AE}-[FY]-[GA]-[FY]-x-[FY]', 'description': 'RNP-1, RNA-binding domain signature'},
    'PS00011': {'name': 'RNP_2', 'pattern': '[LIVM]-[FY]-{IV}-V-x-N-x-[LIVM]', 'description': 'RNP-2, RNA-binding domain signature'},
    'PS00012': {'name': 'HOMEOBOX', 'pattern': '[LIVMFE]-P-[STAN]-x(3)-[LIVM]-E-[LIVM]-x(3)-[YW]-[FYW]-[JNRQEK]-N-[RK]', 'description': 'Homeobox domain signature'},
    'PS00013': {'name': 'HMG_BOX', 'pattern': '[KR]-x(5,13)-[APV]-x(2)-[FYL]-x(2)-[FYW]-R-x-M-[SA]-[EK]-[KR]', 'description': 'HMG box signature'},
    'PS00014': {'name': 'BZIP_BASIC', 'pattern': '[DEQKRNTNSADRKE]-x(14,25)-[KRAQECGRKQEKDKER]-x-[RKQSGDEQKRN]-x(2,4)-[STA]-x-[STAG]-R', 'description': 'Basic region leucine zipper (bZIP) signature'},
    'PS00015': {'name': 'HELIX_TURN_HELIX', 'pattern': '[IVFYDENSTAVKR]-[EDKRTASGVYLI]-[SENDAKRQTI]-x-[LIYVMFAE]-[APLSIVYCTN]-x(2)-[LIVMFWYSTA]-x-[GRKASTNQEL]-[VILMFWYA]', 'description': 'Helix-turn-helix DNA-binding domain signature'},
    'PS00016': {'name': 'BHLH', 'pattern': 'x-[ENQKR]-x(2)-[LIVT]-[LIVMFA]-x(3,7)-[DE]-x-R-x-R-x(2)-K-[LIVMFY]-x(2)-[LIVMFY]-x-[LIVMFY]', 'description': 'Helix-loop-helix DNA-binding domain signature'},
    'PS00017': {'name': 'ATP_BINDING', 'pattern': '[AG]-x(4)-G-K-[ST]', 'description': 'ATP-binding region signature (ABC transporters)'},
    'PS00018': {'name': 'GLY_RICH', 'pattern': 'G-x-G-x-G-x-[GSAE]-x-G', 'description': 'Glycine-rich RNA-binding protein signature'},
    'PS00019': {'name': 'ASPARTYL_PROTEASE', 'pattern': '[LIVMFGAC]-[LIVMFYWSTAGN]-D-T-[GSTAV]-[LIVMFYW]-[SGC]', 'description': 'Eukaryotic aspartyl protease active site'},
    'PS00020': {'name': 'SERINE_PROTEASE_TRYPSIN', 'pattern': '[ST]-x-[DE]-x(10,14)-G-D-S-G-G-P', 'description': 'Serine proteases, trypsin family histidine active site'},
    'PS00021': {'name': 'SERINE_PROTEASE_SUBTILASE', 'pattern': '[GSTANE]-[FYHDS]-G-[TIV]-[S]-[STAG]-[VA]-P', 'description': 'Serine proteases, subtilase family serine active site'},
    'PS00022': {'name': 'CYSTEINE_PROTEASE_C1', 'pattern': '[LIVMFC]-x(2)-[STAV]-x(2,4)-[LIVMFYWDS]-Q-x(2)-C-G-x-[STAN]-[SGA]', 'description': 'Cysteine proteases C1 (papain family) active site'},
    'PS00023': {'name': 'METALLO_PROTEASE_M1', 'pattern': '[VILMFYW]-H-E-[LIVMFY]-H-A-[LIVMFYG]', 'description': 'Metalloproteases M1 family (aminopeptidases) HEXXH zinc-binding region'},
    'PS00024': {'name': 'PROTEIN_KINASE', 'pattern': '[LIVM]-[RK]-x(3)-G-T-x(2)-[FY]-x-[LIVM]', 'description': 'Eukaryotic protein kinases ATP-binding region signature'},
    'PS00025': {'name': 'PROTEIN_KINASE_ATP', 'pattern': '[G]-G-x-G-[FY]-x-[VILM]-{0,30}-[A]-x-K', 'description': 'Protein kinases catalytic domain signature'},
    'PS00026': {'name': 'PROTEIN_PHOSPHATASE', 'pattern': '[LIVMFY]-[LIVMFWSTAG]-G-x-[LIVMFC]-H-[GAS]-[GAS]-D-x-[LIVM]-x(2,3)-[GN]', 'description': 'Serine/threonine protein phosphatases signature'},
    'PS00027': {'name': 'TYROSINE_PHOSPHATASE', 'pattern': '[LIVMF]-(x)-H-C-[SAG]-x-G-[GSAT]-G-R', 'description': 'Protein tyrosine phosphatases active site signature'},
    'PS00028': {'name': 'G_PROTEIN', 'pattern': '[GK]-x(4)-G-K-[ST]-[LIVMFYW]-[LIVMFYW]', 'description': 'G-proteins G1 (GXXXXGKS/T) motif'},
    'PS00029': {'name': 'SH2', 'pattern': 'G-[LIVMFY]-x(2)-[LIVFYW]-x(2,3)-[LIVMFW]-[LIVMFY]-S-x-W-F-x-G-[LEA]-R', 'description': 'SH2 domain signature'},
    'PS00030': {'name': 'SH3', 'pattern': 'x-[LIVMFY]-x(7,9)-[FYW]-x(1,2)-Y-[ND]-x-[AGP]-[FYW]', 'description': 'SH3 domain signature'},
    'PS00031': {'name': 'PH_DOMAIN', 'pattern': '[RK]-x(2)-[KR]-[LIVM]-x(7,12)-[LIVMFYW]-x(2)-W-x(3)-[LIVM]-x-K', 'description': 'Pleckstrin homology (PH) domain signature'},
    'PS00032': {'name': 'WW_DOMAIN', 'pattern': 'W-x-[FY]-x(2)-[FY]-x-H-x-[LIVM]-[LIVMT]-x(5)-G-[WY]-x(3,8)-W-x-P', 'description': 'WW domain signature'},
    'PS00033': {'name': 'LIM_DOMAIN', 'pattern': 'C-x(2)-C-x(16,23)-H-x(2)-[CH]-x(2)-C-x(2)-C-x(16,21)-C-x(2,4)-[CDH]', 'description': 'LIM domain signature'},
    'PS00034': {'name': 'EF_HAND', 'pattern': 'D-x-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]', 'description': 'EF-hand calcium-binding domain'},
    'PS00035': {'name': 'HEXOKINASE', 'pattern': '[LIVM]-[LIVMFAC]-x(2)-P-T-[LIVM]-x(2)-G-L', 'description': 'Hexokinase signature 1'},
    'PS00036': {'name': 'CYTOCHROME_C', 'pattern': 'C-x-[CH]-{4,30}-C-[LH]-x(2)-C', 'description': 'Cytochrome c family heme-binding site signature'},
    'PS00037': {'name': 'CYTOCHROME_B5', 'pattern': '[HY]-x(2)-[PG]-[HY]-x(2)-[SA]', 'description': 'Cytochrome b5 heme-binding domain signature'},
    'PS00038': {'name': 'FLAVOPROTEIN', 'pattern': '[AG]-x-[GSTACV]-x(2)-G-[LIVM]-x-[SADG]-x(7,11)-x-[LIVMFYW]-x-E-[LIVMFYW]', 'description': 'Flavoproteins FAD/FAD-binding domain signature 1'},
    'PS00039': {'name': 'CU_ZN_SOD', 'pattern': '[GDNEQH]-[FYWP]-H-[IVFYW]-H-[STNQ]-x(2)-[GSTAPIMVQH]-x(2)-G', 'description': 'Copper/zinc superoxide dismutase signature'},
    'PS00040': {'name': 'FE_MN_SOD', 'pattern': 'D-x-[WE]-[WFY]-H-A-[WYF]', 'description': 'Iron/manganese superoxide dismutases signature'},
    'PS00041': {'name': 'CATALASE', 'pattern': '[RKHPSA][LIVM]-x(2)-[LIVM]-R-[LIVMGA]-x(2,3)-[GSAT]-x-[HRKD]-x(2,6)-P', 'description': 'Catalase heme-binding signature'},
    'PS00042': {'name': 'PEROXIDASE', 'pattern': '[LIVM]-[LIVM]-x-[LIVMFYW]-[STAG]-x-R-{0,15}-H-x(2)-[LIVM]-x(3,5)-[DENQK]-[LIVM]-[DE]', 'description': 'Plant heme-dependent peroxidases active site signature'},
    'PS00043': {'name': 'THIOREDOXIN', 'pattern': '[LIVMF]-[LIVM]-W-C-[GAPC]-[PAS]-[CRA]-[KRQSN]', 'description': 'Thioredoxin family active site signature'},
    'PS00044': {'name': 'GLUTAREDOXIN', 'pattern': '[FY]-x(2)-[ST]-C-[PAGV]-[YWF]-C', 'description': 'Glutaredoxin family active site signature'},
    'PS00045': {'name': 'FERREDOXIN', 'pattern': 'C-x-[AIV]-x(2)-C-x(2,3)-C-x(3)-[LIVMCPA]-x-P', 'description': 'Ferredoxin iron-sulfur binding region signature'},
    'PS00046': {'name': 'RUBREDOXIN', 'pattern': 'C-x(2)-C-x(29,32)-C-x(2)-C', 'description': 'Rubredoxin iron-binding region signature'},
    'PS00047': {'name': 'TREFORIAL', 'pattern': 'C-x(3,14)-C-x(3,6)-C-x(4,11)-C-x(3,6)-C-x(4,16)-C', 'description': 'Trefoil (P-type) domain signature'},
    'PS00048': {'name': 'EGF_CA', 'pattern': 'C-x(3,5)-S-x(6,12)-C-x(2)-G-x(2,3)-C-x(2,7)-C', 'description': 'Calcium-binding EGF domain'},
    'PS00049': {'name': 'ANAPHYLATOXIN', 'pattern': 'C-x(28,36)-C-x(6,10)-C-x(2)-C-x(7,16)-C', 'description': 'Anaphylatoxin-like domain signature'},
    'PS00050': {'name': 'SRC_HOMOLOGY_3', 'pattern': 'x-[LIVMFY]-x(7,9)-[FYW]-x(1,2)-Y-[ND]-x-[AGP]-[FYW]', 'description': 'SRC homology 3 (SH3) domain signature'},
}


def prosite_to_regex(pattern):
    result = ''
    i = 0
    p = pattern.upper()
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    while i < len(p):
        c = p[i]
        if c == 'X':
            result += '.'
            i += 1
        elif c == '[':
            j = i + 1
            while j < len(p) and p[j] != ']':
                j += 1
            result += p[i:j+1]
            i = j + 1
        elif c == '{':
            j = i + 1
            while j < len(p) and p[j] != '}':
                j += 1
            chars = p[i+1:j]
            result += '[^' + chars + ']'
            i = j + 1
        elif c == '(':
            j = i + 1
            while j < len(p) and p[j] != ')':
                j += 1
            nums = p[i+1:j].split(',')
            if len(nums) == 1:
                result += '{' + nums[0] + '}'
            else:
                result += '{' + nums[0] + ',' + nums[1] + '}'
            i = j + 1
        elif c == '-':
            i += 1
        elif c in amino_acids:
            result += c
            i += 1
        else:
            result += c
            i += 1
    return result
