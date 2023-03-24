import plotly.graph_objects as go
import numpy as np
import os
import re
from typing import List
from data import protein
from data import residue_constants
from scipy.spatial.transform import Rotation

CA_IDX = residue_constants.atom_order['CA']

def create_scatter(pos_3d,
                   mode='markers',
                   marker_size=None,
                   name=None,
                   opacity=None,
                   color=None,
                   colorscale=None,
                   ):
    """Creates Scatter3D objects for use in plotly.
    Args:
        pos_3d: [N, 3] array containing N points with
            euclidean coordinates.
        mode: How to display points.
            Use 'markers' for scatter.
            Use 'lines' for lines connecting consecutive points.
            Use 'lines+markers' for scatter and lines.
        marker_size: Size of markers.
        name: Label of plotting layer to be displayed in legend.
        opacity: Transparency of points.
    """
    x, y, z = np.split(pos_3d, 3, axis=-1)
    args_dict = {
        'x': x[:, 0],
        'y': y[:, 0],
        'z': z[:, 0],
        'mode': mode,
        'marker': {}
    }
    if marker_size:
        args_dict['marker']['size'] = marker_size
    if name:
        args_dict['name'] = name
    if opacity:
        args_dict['opacity'] = opacity
    if color:
        args_dict['marker']['color'] = color
    if colorscale:
        args_dict['marker']['colorscale'] = colorscale
    return go.Scatter3d(**args_dict)

def write_prot_to_pdb(
        prot_pos: np.ndarray, file_path: str, overwrite=False, no_indexing=False):
    if overwrite:
        max_existing_idx = 0
    else:
        file_dir = os.path.dirname(file_path)
        file_name = os.path.basename(file_path).strip('.pdb')
        existing_files = [x for x in os.listdir(file_dir) if file_name in x]
        max_existing_idx = max([
            int(re.findall(r'_(\d+).pdb', x)[0]) for x in existing_files if re.findall(r'_(\d+).pdb', x)
            if re.findall(r'_(\d+).pdb', x)] + [0])
    if no_indexing:
        save_path = file_path.strip('.pdb') + '.pdb'
    else:
        save_path = file_path.strip('.pdb') + f'_{max_existing_idx+1}.pdb'
    with open(save_path, 'w') as f:
        if prot_pos.ndim == 3:
            for t, bb_pos in enumerate(prot_pos):
                bb_prot = create_bb_prot(bb_pos)
                pdb_prot = protein.to_pdb(bb_prot, model=t + 1, add_end=False)
                f.write(pdb_prot)
        elif prot_pos.ndim == 2:
            bb_prot = create_bb_prot(prot_pos)
            pdb_prot = protein.to_pdb(bb_prot, model=1, add_end=False)
            f.write(pdb_prot)
        else:
            raise ValueError(f'Invalid positions shape {prot_pos.shape}')
        f.write('END')
    return save_path

def create_bb_prot(bb_pos: np.ndarray):
    assert bb_pos.ndim == 2
    assert bb_pos.shape[1] == 3
    n = bb_pos.shape[0]
    imputed_atom_pos = np.zeros([n, 37, 3])
    imputed_atom_pos[:, CA_IDX] = bb_pos
    imputed_atom_mask = np.zeros([n, 37])
    imputed_atom_mask[:, CA_IDX] = 1.0
    residue_index = np.arange(n)
    chain_index = np.zeros(n)
    b_factors = np.zeros([n, 37])
    aatype = np.zeros(n, dtype=np.int)
    return protein.Protein(
        atom_positions=imputed_atom_pos,
        atom_mask=imputed_atom_mask,
        aatype=aatype,
        residue_index=residue_index,
        chain_index=chain_index,
        b_factors=b_factors)

