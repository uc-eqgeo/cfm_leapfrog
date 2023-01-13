import os
import glob
import subprocess
import meshio
from fault_mesh.utilities.cubit import make_journal_file_surface

indir=r"Z:\Penney\leapfrog_cfm\kaikoura_faults\new_version_19_12_22\leapfrog_meshes"
outdir=indir
print(indir,os.path.exists(indir))
#make .jou files

for file in glob.iglob(f'{indir}/2. Hope combined.obj'):
#for file in glob.iglob(f'{indir}/*.obj'):
    print(file)
    #for file in glob.iglob(f'{indir}/1. Alpine combined.obj'):
    fault = file.split('.')[1].replace(" ", "")
    print(fault)

    # convert to .stl before using with meshio
    mesh=meshio.read(file)
    inmesh = os.path.join(indir,fault+".stl")
    meshio.write(inmesh,mesh)


    outmesh=os.path.join(outdir,fault+"_remeshed.stl")
    outjou=os.path.join(outdir,fault+".jou")
    make_journal_file_surface(inmesh,outmesh,outjou,mesh_size=2000,min_mesh_size=1500)

    subprocess.run(["coreform_cubit.exe", '-nographics', '-nojournal', outjou],executable=r"C:\Program Files\Coreform Cubit 2022.11\bin\coreform_cubit.exe")

#
