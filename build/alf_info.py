import numpy as np
alf_info={}
alf_info['name']='sting_cplx'          # must match prep/sting_cplx.py
alf_info['nsubs']=[7]                  # ONE coupled site: substituent N switches on in LIGA and LIGB together
alf_info['impcons']='custom'           # use G_imp/ (theta_bias_ind8.5, fnex=8.5); without this ALF silently ignores G_imp
alf_info['nblocks']=np.sum(alf_info['nsubs'])
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']=''
alf_info['temp']=298.15
# all 7 ligands are neutral (protonated acid, matching the CHARMM-GUI QAD) -> no alf_info['q'] / PME charge-change correction
