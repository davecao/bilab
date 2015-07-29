# -*- coding: utf-8 -*-

from bilab.structure.measure import buildDistMatrix

__all__ = ['density1d']

def density1d(mol, cutoff=10.0, save2file=False, selection='protein and name CA'):
    """
    Args:
        mol (class)        : bilab.structure.AtomGroup
        cutoff (float)     : cutoff of distance for a residue 
                             to generate a neighbor list
        selection (string) : the atom name for selection, 
                            default is CA of a residue.

    Returns:
        list of probilities of residues.

    Note: AtomGroup.select much faster than findNeighbors
    """
    atoms = mol.select(selection)
    hv = atoms.getHierView()
    chains = list(hv)
    # iterate over residues
    #for i, residue in enumerate(hv.iterResidues()):
    #    print residue
    result = {}
    mol_title = 'UNNAMED'
    if save2file:
        if atoms.getTitle() != 'Unknown':
            mol_title = atoms.getTitle()

    for ch in chains:
        ch_id = ch.getChid()
        distMat = bilab.structure.measure.buildDistMatrix(ch)
        distMat = np.asarray(distMat)
        # normalization by the maximum of distance
        n_distMat = 1.0*distMat/distMat.max()
        # Student distribution degenerates difference of distances
        t_distMat = np.reciprocal(1+distMat*distMat)
        dens = []
        result[ch_id] = {'org_dist':distMat, 'no_distMat':n_distMat,'t_dist':t_distMat, 'density':dens}
        if save2file:
            base_fname = mol_title + '_' + ch_id
            np.savetxt( base_fname + '_distMat.txt', distMat , delimiter=',', fmt="%8.3f")
            np.savetxt( base_fname + '_ndistMat.txt', n_distMat, delimiter=',', fmt="%8.3f")
            np.savetxt( base_fname + '_tdistMat.txt', t_distMat, delimiter=',', fmt="%8.3f")

    return result
