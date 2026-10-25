"""Mutations address exact coordinates and H bonds require valid atom chemistry."""
from io import StringIO

import numpy as np
import pytest

app = pytest.importorskip('openmm.app')
from openmm import Vec3, unit
PDBFixer = pytest.importorskip('pdbfixer').PDBFixer
from app.hbonds import find_hydrogen_bonds
from app.repack import _mutate_exact, RepackError


def duplicate_number_fixer():
    topology = app.Topology()
    chain = topology.addChain('H')
    positions = []
    for j, insertion in enumerate(['', 'A']):
        residue = topology.addResidue('ALA', chain, id='112', insertionCode=insertion)
        for i, name in enumerate(['N', 'CA', 'C', 'O']):
            element = app.element.nitrogen if name == 'N' else app.element.oxygen if name == 'O' else app.element.carbon
            topology.addAtom(name, element, residue)
            positions.append(Vec3(j + i / 10, 0, 0))
    handle = StringIO()
    app.PDBFile.writeFile(topology, positions * unit.nanometer, handle, keepIds=True)
    return PDBFixer(pdbfile=StringIO(handle.getvalue()))


@pytest.mark.parametrize('target,expected', [('112', ['GLY', 'ALA']), ('112A', ['ALA', 'GLY'])])
def test_exact_insertion_identity_survives_pdbfixer(target, expected):
    fixer = duplicate_number_fixer()
    _mutate_exact(fixer, 'H', target, 'ALA', 'GLY')
    residues = list(fixer.topology.residues())
    assert [r.name for r in residues] == expected
    assert [(r.id, r.insertionCode.strip()) for r in residues] == [('112', ''), ('112', 'A')]
    assert [r.chain.id for r in residues] == ['H', 'H']


def test_duplicate_exact_identity_fails_closed():
    fixer = duplicate_number_fixer()
    list(fixer.topology.residues())[1].insertionCode = ''
    with pytest.raises(RepackError, match='unique coordinate identity'):
        _mutate_exact(fixer, 'H', '112', 'ALA', 'GLY')


def model_bonds(resname, atomname, *, protonated=False, other_ring_h=False):
    topology = app.Topology()
    chain = topology.addChain('H')
    donor_residue = topology.addResidue('ALA', chain, id='1')
    donor = topology.addAtom('N', app.element.nitrogen, donor_residue)
    hydrogen = topology.addAtom('H', app.element.hydrogen, donor_residue)
    topology.addBond(donor, hydrogen)
    target = topology.addResidue(resname, chain, id='10')
    element = app.element.oxygen if atomname.startswith('O') else app.element.sulfur if atomname.startswith('S') else app.element.nitrogen
    acceptor = topology.addAtom(atomname, element, target)
    coords = [[0., 0, 0], [1., 0, 0], [2.8, 0, 0]]
    if protonated:
        h = topology.addAtom('HX', app.element.hydrogen, target)
        topology.addBond(acceptor, h)
        coords.append([2.8, 1., 0])
    if other_ring_h:
        other = topology.addAtom('NE2' if atomname == 'ND1' else 'ND1', app.element.nitrogen, target)
        h = topology.addAtom('HR', app.element.hydrogen, target)
        topology.addBond(other, h)
        coords.extend([[7., 0, 0], [8., 0, 0]])
    return [b for b in find_hydrogen_bonds(topology, np.array(coords)) if b.donor_resid == '1']


@pytest.mark.parametrize('resname,atomname', [('LYS','NZ'), ('ARG','NE'), ('ARG','NH1'),
                                              ('ASN','ND2'), ('GLN','NE2'), ('TRP','NE1'),
                                              ('ALA','N'), ('PRO','N'), ('CYX','SG')])
def test_nonacceptor_nitrogens_and_disulfide_are_excluded(resname, atomname):
    assert model_bonds(resname, atomname) == []


@pytest.mark.parametrize('resname,atomname', [('ALA','O'), ('ASP','OD1'), ('GLU','OE1'),
                                              ('ASN','OD1'), ('GLN','OE1'), ('SER','OG'), ('MET','SD')])
def test_valid_acceptors_still_produce_geometry(resname, atomname):
    assert len(model_bonds(resname, atomname)) == 1


def test_histidine_depends_on_explicit_tautomer():
    assert model_bonds('HIS', 'ND1') == []
    assert len(model_bonds('HIS', 'ND1', other_ring_h=True)) == 1
    assert model_bonds('HIS', 'ND1', protonated=True, other_ring_h=True) == []
    assert len(model_bonds('HIE', 'ND1')) == 1
    assert model_bonds('HIP', 'ND1') == []


def test_protonated_carboxyl_oxygen_not_acceptor():
    assert model_bonds('ASH', 'OD2', protonated=True) == []
    assert len(model_bonds('SER', 'OG', protonated=True)) == 1
