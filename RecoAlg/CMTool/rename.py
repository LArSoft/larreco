import os

hh = [x for x in os.listdir('.') if x.endswith('.h')]

cxx = [x for x in os.listdir('.') if x.endswith('.cxx')]

files = hh + cxx

for f in files:

    contents = open(f,'r').read()

    contents = contents.replace('RecoUtilException','CRUException')
    contents = contents.replace('ClusterParamsAlgNew','ClusterParamsAlg')
    contents = contents.replace(".hh\"",".h\"")
    contents = contents.replace('\"ClusterParamsAlg.h\"','\"RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h\"')
    contents = contents.replace('\"CRUException.h\"','\"RecoAlg/ClusterRecoUtil/CRUException.h\"')
    contents = contents.replace('\"GeometryUtilities.h\"','\"Utilities/GeometryUtilities.h\"')
    contents = contents.replace("larutil","util")

    newfile = open(f,'w')
    newfile.write(contents)
    newfile.close()
    
