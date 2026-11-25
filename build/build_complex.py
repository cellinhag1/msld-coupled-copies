##
## pyCHARMM system build for the STING COMPLEX leg, coupled two-copy MSLD.
##
## System: CHARMM-GUI solvated 6UKM dimer (complex_prep/step3_pbcsetup.psf/crd,
## loops rebuilt by the user). The two crystal QAD copies (segid HETA, HETB) are
## deleted and rebuilt as alchemical hybrid ligands LIGA (in HETA's pose) and
## LIGB (in HETB's pose), each = one shared core + 7 chain substituents.
##
## COUPLING: one alchemical site, nsubs = [7]. Substituent N of copy A and
## substituent N of copy B are put into ONE selection -> ONE BLOCK -> ONE
## lambda. So both copies always carry the same substituent. ALF/BLOCK see an
## ordinary 1-site, 7-substituent system (nblocks = 7).
##   * A-subN <-> B-subM (N != M) is excluded by adex, as it must be: those two
##     never coexist in a real molecule pair.
##   * A-subN <-> B-subN interact inside block N, scaled by lambda_N (not
##     lambda_N^2). Exact at the physical endpoints (lambda = 0 or 1).
##
## Exec'd by msld_flat.py / msld_prod.py from inside run<i>/ (prep -> ../prep).
## Expects nsubs, nblocks, bias, temp from variables*.py.
##

import numpy as np
import pycharmm
import pycharmm.generate as gen
import pycharmm.read as read
import pycharmm.select as select
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.settings as settings

##############################################
# Global parameters

sysname = 'sting_cplx'
box = 93.0            # CHARMM-GUI step3_pbcsetup.str (cubic, centred at origin)
pmegrid = 96          # smallest FFT-friendly (2/3/5) even integer >= box
builddir = 'prep'
fnex = 8.5            # matches G_imp (theta_bias_ind8.5) and the solvent coupled test

nsites = len(nsubs)
assert nsites == 1 and nsubs[0] == 7, "expected one coupled site with 7 substituents"

copies = [('LIGA', 'a'), ('LIGB', 'b')]
pert = {0: {'subs': [str(i + 1) for i in range(nsubs[0])], 'resid': '1', 'type': 'ligand'}}

##############################################
# Topology / parameters: CHARMM-GUI toppar (CGenFF 5.0, matches the ligand
# CGenFF runs), then the hybrid ligand core + all substituent patches.

settings.set_bomb_level(-2)
pycharmm.lingo.charmm_script('stream ' + builddir + '/toppar_complex.str')
read.rtf(builddir + '/core.rtf', append=True)
read.prm(builddir + '/full_ligand.prm', append=True, flex=True)
for sub in pert[0]['subs']:
    read.rtf(builddir + '/site1_sub' + sub + '_pres.rtf', append=True)

##############################################
# Environment FIRST (plain, non-append read of the full CHARMM-GUI PSF keeps
# all its patches: termini, disulfides, protonation), then drop crystal QAD.

read.psf_card(builddir + '/step3_pbcsetup.psf')
read.coor_card(builddir + '/step3_pbcsetup.crd')
pycharmm.lingo.charmm_script('dele atom sele segid HETA .or. segid HETB end')

# The protein/CGenFF RTFs switch on autogen after EVERY patch, which would
# regenerate angles/dihedrals for the whole ~76k-atom system 14 times.
# Turn that off; autogen runs once below, on the two ligands only.
pycharmm.lingo.charmm_script('autogen nopatch')

# Symmetric protonation: CHARMM-GUI applied ASPP to PROA 205 only. Both chains'
# Asp205 have identical surroundings (helix N-cap, backbone-NH contacts only),
# so make PROA 205 a standard charged ASP like PROB 205, and remove one bulk
# chloride (CLA resid 54, ~45 A from protein) to keep the box neutral.
read.rtf(builddir + '/aspd.rtf', append=True)
pycharmm.lingo.charmm_script('patch ASPD PROA 205 setup')
pycharmm.lingo.charmm_script('dele atom sele segid IONS .and. resid 54 .and. resname CLA end')
pycharmm.lingo.charmm_script('set envq ?cgtot')
envq = float(pycharmm.lingo.get_charmm_variable('ENVQ'))
assert abs(envq) < 1e-4, 'environment net charge after ASP205 fix = {}'.format(envq)
nhd2 = sum(pycharmm.SelectAtoms().by_res_and_type('PROA PROB', '205', 'HD2').get_selection())
assert nhd2 == 0, 'Asp205 HD2 still present in {} chain(s)'.format(nhd2)
print('ASP205 CHECK: both chains charged ASP, environment net charge {:.5f}'.format(envq))

# CHARMM-GUI's toppar streams restore their own bomlev on exit, so set it again:
# the core alone carries a non-integer charge (-1.288) until the patches below
# restore integer charge, which trips a level-0 PSFSUM abort otherwise.
settings.set_bomb_level(-2)

##############################################
# Hybrid ligands: core + 7 substituent patches, per copy, each copy with its
# own aligned coordinates (copy A from HETA's pose, copy B from HETB's).

for segid, tag in copies:
    core_pdb = builddir + '/core' + tag + '.pdb'      # PDB segid column = LIGA / LIGB
    read.sequence_pdb(core_pdb)
    gen.new_segment(segid, setup_ic=True)              # pycharmm default: no terminal patches
    read.pdb(core_pdb, resid=True)

for sub in pert[0]['subs']:
    for segid, tag in copies:
        pycharmm.lingo.charmm_script('patch p1_{} {} {} setup'.format(sub, segid, pert[0]['resid']))
        read.pdb(builddir + '/site1_sub' + sub + '_frag' + tag + '.pdb', resid=True)

# fill any coordinate the PDBs did not provide (should be none -- checked below)
pycharmm.lingo.charmm_script('''
ic generate sele segid LIGA .or. segid LIGB end
ic param
ic build
''')

# lone pairs: lpsites.inp is empty for this panel; streamed once per copy so a
# future halogenated analog works (it references @ligseg/@resnum)
for segid, tag in copies:
    pycharmm.lingo.charmm_script('set ligseg {}\nset resnum {}'.format(segid, pert[0]['resid']))
    pycharmm.lingo.charmm_script('stream ' + builddir + '/lpsites.inp')

pycharmm.lingo.charmm_script('auto angle dihe sele segid LIGA .or. segid LIGB end')

##############################################
# COUPLED selections: substituent N = copy-A atoms OR copy-B atoms

pert[0]['select'] = []
for sub in pert[0]['subs']:
    selname = 'site1sub' + sub
    pert[0]['select'].append(selname)
    names = [l.split()[1].upper() for l in open(builddir + '/site1_sub' + sub + '_pres.rtf') if l[0:4] == 'ATOM']
    atoms = ' '.join(names)
    sel = None
    for segid, tag in copies:
        s = pycharmm.SelectAtoms().by_res_and_type(segid, pert[0]['resid'], atoms)
        assert sum(s.get_selection()) == len(names), "{} sub{}: selected {} of {} atoms".format(
            segid, sub, sum(s.get_selection()), len(names))
        sel = s if sel is None else (sel | s)
    select.store_selection(selname, sel)

# remove bonded terms between different substituents (autogen fabricates
# angles/dihedrals across the shared core); no bonds exist between copies
settings.set_bomb_level(-1)
for i in range(nsubs[0]):
    for j in range(i + 1, nsubs[0]):
        pycharmm.lingo.charmm_script('dele connectivity sele {} show end sele {} show end'.format(
            pert[0]['select'][i], pert[0]['select'][j]))

##############################################
# Hard stop on any uninitialized coordinate (silent 9999.0 placeholders)

pycharmm.lingo.charmm_script('define junk sele .not. init end\nset nuninit ?nsel')
nuninit = int(float(pycharmm.lingo.get_charmm_variable('NUNINIT')))
if nuninit > 0:
    pycharmm.lingo.charmm_script('print coor sele .not. init end')
    raise RuntimeError('{} atoms have uninitialized coordinates'.format(nuninit))
print('BUILD CHECK: 0 uninitialized atoms')
settings.set_bomb_level(0)

##############################################
# Crystal / images: cutoff box/2 (ME1 guide: small fixed cutoffs caused
# water depletion at the box edge)

crystal.define_cubic(box)
crystal.build(box / 2.0)
for seg in ('PROA', 'PROB', 'LIGA', 'LIGB'):
    image.setup_segment(0.0, 0.0, 0.0, seg)
for res in ('TIP3', 'SOD', 'CLA'):
    image.setup_residue(0.0, 0.0, 0.0, res)

pycharmm.lingo.charmm_script('coor copy comp')

##############################################
# BLOCK setup (same generic logic as the solvent coupled test; the coupling
# is already in the selections above)

blockplusone = nblocks + 1
block_init = '''
!! BLOCK setup
BLOCK {}
   clear
END
BLOCK {}
'''.format(blockplusone, blockplusone)

block_call = ''
pert[0]['buffer'] = []
ii = 2
for sub in range(nsubs[0]):
    block_call += 'Call {} sele {} show end\n'.format(ii, pert[0]['select'][sub])
    pert[0]['buffer'].append(ii)
    ii += 1

block_parm = '''
qldm theta
lang temp {}
soft w14
pmel ex

ldin 1 1.0  0.0  5.0  0.0  5.0'''.format(temp)

block_ldin = ''
sitestr = ''
for sub in range(nsubs[0]):
    tmplmb = 1.0 - (0.01 * (nsubs[0] - 1)) if sub == 0 else 0.01
    block_ldin += 'ldin {} {:.4f} 0.0 5.0 {} 5.0\n'.format(pert[0]['buffer'][sub], tmplmb, bias['b'][0, sub])
    sitestr += '1  '

block_adex = ''
for i in range(nsubs[0]):
    for j in range(i + 1, nsubs[0]):
        block_adex += 'adex {} {}\n'.format(pert[0]['buffer'][i], pert[0]['buffer'][j])

block_msld = '''
rmla bond thet impr
msld 0  {} fnex {}
msma
'''.format(sitestr, fnex)

block_varb = 'ldbi {}\n'.format((5 * nblocks * (nblocks - 1)) // 2)
ibias = 0
for i in range(nsubs[0]):
    for j in range(i + 1, nsubs[0]):
        bi, bj = pert[0]['buffer'][i], pert[0]['buffer'][j]
        block_varb += 'ldbv {} {} {}  6  0.00 {} 0\n'.format(ibias + 1, bi, bj, -bias['c'][i, j])
        block_varb += 'ldbv {} {} {} 10 -5.56 {} 0\n'.format(ibias + 2, bi, bj, -bias['x'][i, j])
        block_varb += 'ldbv {} {} {}  8 0.017 {} 0\n'.format(ibias + 3, bi, bj, -bias['s'][i, j])
        block_varb += 'ldbv {} {} {} 10 -5.56 {} 0\n'.format(ibias + 4, bj, bi, -bias['x'][j, i])
        block_varb += 'ldbv {} {} {}  8 0.017 {} 0\n'.format(ibias + 5, bj, bi, -bias['s'][j, i])
        ibias += 5
block_varb += 'END\n'

pycharmm.lingo.charmm_script('\n'.join([block_init, block_call, block_parm, block_ldin,
                                        block_adex, block_msld, block_varb]))
