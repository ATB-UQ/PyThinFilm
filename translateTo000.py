import pmx
import sys
model = pmx.Model(sys.argv[1])

model.translate(map(lambda x:-x, model.com(vector_only=True)))
model.write("cbp000.pdb", "CBP molecule at 0,0,0", 0)
