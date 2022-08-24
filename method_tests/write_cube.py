import numpy as np

def write_cube(filepath:str, vox_grid:np.ndarray, coords:np.ndarray, elements:list[str]):
    with open(filepath, 'w') as cube:
        # write header stuff
        cube.write('title\n')
        cube.write('title2\n')
        cube.write('num species >> origin xyz\n')
        cube.write('voxel/lattice info1 \n')
        cube.write('voxel/lattice infor \n')
        cube.write('voxel/lattice info3 \n')

        for element, [x, y, z] in zip(elements, coords):
            cube.write(f'{element} 0.0 {x} {y} {z}\n')
        
        # write voxel grid to cube file from tensor
        nvox = np.prod(vox_grid.shape)
        cube_vox_lines = nvox//6
        leftover_vox = nvox%6

        vox_flat = vox_grid.flatten()

        for k in range(0,cube_vox_lines,6):
            cube.write(f'{vox_flat[k]} {vox_flat[k+1]} {vox_flat[k+2]} {vox_flat[k+3]} {vox_flat[k+4]} {vox_flat[k+5]}\n')

        cube.write(f' '.join(map(str, vox_flat[-leftover_vox:])))

