import numpy as np
import h5py as h5

from molmod.units import kelvin, bar, angstrom, kjmol, liter, amu
from molmod.constants import avogadro

import yaff
from yaff.pes.eos import PREOS
from yaff.sampling.io import MCHDF5Writer
from yaff.sampling.mc import GCMC
from yaff.sampling.mcutils import MCScreenLog
from yaff.system import System
from yaff import log
log.set_level(log.medium)
import datetime

from gemmi import cif
import getpass

def simulate():
    T = 298.0*kelvin
    # Input files
    fn_guest = 'CO2.chk'
    fn_host = 'MIL53.chk'
    fn_pars = ['pars.txt']
    host = System.from_file(fn_host).supercell(1,1,1)
    # Description of allowed MC moves and their corresponding probabilities
    mc_moves =  {'insertion':1.0, 'deletion':1.0,
                 'translation':1.0, 'rotation':1.0}
    # Construct equation of state to link pressure, fugacity and chemical potential
    eos = PREOS.from_name('carbondioxide')
    # Loop over pressures to construct isotherm
    pressures = np.array([0.1,0.5,1.0,3.0,5.0,10.0])*bar
    #pressures = np.array([0.1])*bar
    fugacities = []
    mus = []
    steps = []
    uptake = np.zeros(pressures.shape)
    for iP, P in enumerate(pressures):
        fugacity = eos.calculate_fugacity(T,P)
        mu = eos.calculate_mu(T,P)
        # Screen logger
        screenlog = MCScreenLog(step=1000)
        # HDF5 logger
        fh5 = h5.File('trajectory_%d.h5'%iP,'w')
        hdf5writer = MCHDF5Writer(fh5, step=1000)
        # Setup the GCMC calculation, this generates a lot of output so when
        # force fields are generated, so we silence the logging for a while.
        log.set_level(log.silent)
        gcmc = GCMC.from_files(fn_guest, fn_pars, host=host,
            rcut=12.0*angstrom, tr=None, tailcorrections=True, hooks=[screenlog, hdf5writer],
            reci_ei='ewald_interaction', nguests=30)
        log.set_level(log.medium)
        # Set the external conditions
        fugacities.append(fugacity)
        mus.append(mu)
        gcmc.set_external_conditions(T, fugacity)
        # Run MC simulation
        gcmc.run(10000, mc_moves=mc_moves)
        uptake[iP] = gcmc.Nmean
        steps.append(gcmc.counter)
        fh5.close()

    # initialise

    d = cif.Document()
    d.add_new_block('yaff2aif')

    block = d.sole_block()
    block.set_pair('_audit_aif_version', '6acf6ef')

    #label metadata

    block.set_pair('_exptl_operator',  getpass.getuser())
    block.set_pair('_simltn_date', datetime.datetime.now().isoformat())
    block.set_pair('_simltn_code', 'yaff'+yaff.__version__)

    block.set_pair('_exptl_method', 'GCMC')
    block.set_pair('_exptl_isotherm_type', 'absolute')

    block.set_pair('_exptl_adsorptive', 'CO2')
    block.set_pair('_exptl_temperature', str(T))

    block.set_pair('_simltn_forcefield_adsorptive', fn_guest+'_'+'_'.join(fn_pars))
    block.set_pair('_simltn_forcefield_adsorbent', fn_host+'_'+'_'.join(fn_pars))

    block.set_pair('_adsnt_material_id', 'MIL-53(Al)')
    #record mass to infer simulation size
    block.set_pair('_adsnt_sample_mass', '%.5E' % np.sum(host.masses/amu))


    block.set_pair('_units_temperature', 'K')
    block.set_pair('_units_energy', 'kJ/mol')
    block.set_pair('_units_loading','molecules')
    block.set_pair('_units_pressure','bar')
    block.set_pair('_units_mass','amu')

    #format adsorption

    pressures_bar = pressures/bar
    fugacities_bar = np.array(fugacities)/bar
    mus_kjmol = np.array(mus)/kjmol
    loop_ads = block.init_loop('_adsorp_', ['pressure', 'amount', 'fugacity', 'chemicalpotential', 'nsteps'])
    loop_ads.set_all_values([
        ['%.5E' % item for item in pressures_bar],
        ['%.5E' % item for item in uptake],
        ['%.5E' % item for item in fugacities_bar],
        ['%.5E' % item for item in mus_kjmol],
        ['%.3E' % item for item in steps]
    ])

    d.write_file('test.aif')






    
    np.save('results.npy', np.array([pressures,uptake]).T)
    


if __name__=='__main__':
    simulate()