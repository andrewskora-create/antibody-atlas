"""Local terminal-only model completion. Never substitutes residues or fills internal gaps."""
import json
import math
import sys
from pathlib import Path


def repair(payload):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from Bio.SeqUtils import seq1, seq3
    from .propermab import _terminal_identity
    fixer=PDBFixer(filename=payload['input_path'])
    chains=list(fixer.topology.chains())
    if {c.id for c in chains} != {'H','L'}:
        raise ValueError('Terminal repair requires exactly H and L chains.')
    fixer.missingResidues={}
    omitted={}
    for chain in chains:
        residues=list(chain.residues())
        observed=''.join(seq1(r.name) for r in residues)
        expected=payload['sequences'][chain.id]
        match=_terminal_identity(expected,observed,chain.id)
        omitted[chain.id]=match['omitted_positions']
        first=match['input_positions'][0]-1
        last=match['input_positions'][-1]
        if first:fixer.missingResidues[(chain.index,0)]=[seq3(a).upper() for a in expected[:first]]
        if last<len(expected):fixer.missingResidues[(chain.index,len(residues))]=[seq3(a).upper() for a in expected[last:]]
    if not any(omitted.values()):raise ValueError('No terminal residues require reconstruction.')
    # PDBFixer positions/minimizes new atoms while retaining original atom coordinates.
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=1729)
    # Refine the added peptide with a bonded force field, not only PDBFixer's
    # overlap-removal potential. Preserve the original heavy-atom model.
    import openmm as mm
    from openmm import unit
    from openmm.app import Modeller, ForceField, NoCutoff, element
    forcefield=ForceField('charmm36.xml','charmm36/tip3p-pme-b.xml')
    modeller=Modeller(fixer.topology,fixer.positions)
    modeller.addHydrogens(forcefield,pH=payload.get('ph',7.4))
    system=forcefield.createSystem(modeller.topology,nonbondedMethod=NoCutoff,constraints=None,removeCMMotion=False)
    for chain in modeller.topology.chains():
        for pos,residue in enumerate(chain.residues(),1):
            if pos not in omitted[chain.id]:
                for atom in residue.atoms():
                    if atom.element!=element.hydrogen:system.setParticleMass(atom.index,0)
    integrator=mm.VerletIntegrator(.001*unit.picoseconds)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'2'})
    context.setPositions(modeller.positions)
    mm.LocalEnergyMinimizer.minimize(context,maxIterations=500)
    fixer.topology=modeller.topology
    fixer.positions=context.getState(getPositions=True).getPositions()
    xyz=fixer.positions.value_in_unit(__import__('openmm').unit.nanometer)
    joins=[]
    for chain in fixer.topology.chains():
        residues=list(chain.residues())
        observed=''.join(seq1(r.name) for r in residues)
        if observed!=payload['sequences'][chain.id]:raise ValueError('Repaired model does not match the full input.')
        # Validate peptide connectivity/geometry at every new residue junction.
        for i in range(len(residues)-1):
            if i+1 not in omitted[chain.id] and i+2 not in omitted[chain.id]:continue
            atoms1={a.name:a for a in residues[i].atoms()};atoms2={a.name:a for a in residues[i+1].atoms()}
            if 'C' not in atoms1 or 'N' not in atoms2:raise ValueError('Rebuilt peptide junction is incomplete.')
            distance=math.dist(xyz[atoms1['C'].index],xyz[atoms2['N'].index])*10
            if not 1.20<=distance<=1.45:raise ValueError('Rebuilt peptide bond failed the geometry check.')
            joins.append({'chain':chain.id,'input_positions':[i+1,i+2],'peptide_bond_angstrom':round(distance,3)})
        for i,residue in enumerate(residues):residue.id=str(i+1);residue.insertionCode=''
    with open(payload['output_path'],'w') as handle:PDBFile.writeFile(fixer.topology,fixer.positions,handle,keepIds=True)
    return {'omitted_positions':omitted,'junctions':joins,'method':'Local PDBFixer terminal reconstruction, seed 1729; CHARMM minimization of rebuilt heavy atoms and hydrogens with original heavy atoms fixed; exact full-input identity and peptide-junction checks.'}


if __name__=='__main__':
    payload=json.load(sys.stdin)
    def offline(event,args):
        if event in {'socket.connect','socket.getaddrinfo','socket.connect_ex'}:raise RuntimeError('Network disabled for local terminal repair.')
    sys.addaudithook(offline)
    try:
        import contextlib
        with open('/dev/null','w') as sink,contextlib.redirect_stdout(sink),contextlib.redirect_stderr(sink):
            result=repair(payload)
        print(json.dumps(result))
    except Exception as exc:
        print(json.dumps({'error':type(exc).__name__+': terminal reconstruction or geometry validation failed.'}))
        sys.exit(1)
