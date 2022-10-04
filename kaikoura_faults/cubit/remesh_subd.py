import os
import glob
import subprocess
import meshio
from fault_mesh_building.fault_mesh.utilities.cubit import make_journal_file_surface

indir="/home/UOCNT/cpe88/PycharmProjects/cfm_leapfrog/kaikoura_faults/Subduction"
outdir=indir

#make .jou files

for file in glob.iglob(f'{indir}/subduction_raoul_special10k_inside_Boundary.obj'):
    #for file in glob.iglob(f'{indir}/1. Alpine combined.obj'):
    fault = "hikurangi"
    print(fault)

    # convert to .stl before using with meshio
    mesh=meshio.read(file)
    inmesh = os.path.join(indir,fault+".stl")
    meshio.write(inmesh,mesh)


    outmesh=os.path.join(outdir,fault+"_remeshed.stl")
    outjou=os.path.join(outdir,fault+".jou")
    make_journal_file_surface(inmesh,outmesh,outjou,mesh_size=2000,min_mesh_size=1500)
    subprocess.run(['/home/UOCNT/cpe88/programs/Coreform-Cubit-2022.4/bin/coreform_cubit', '-nographics', '-nojournal', outjou])


