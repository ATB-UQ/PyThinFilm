import __main__
from random import randint

__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import pymol
from pymol import cmd
from pymol.cgo import LINEWIDTH,  BEGIN, LINES, COLOR, VERTEX, END

DPI = 600
ANTIALIAS = 2
RENDER_PROPERLY = True
def draw_bounding_box(selection="(all)", padding=0.0, linewidth=1.0, r=1.0, g=0.0, b=0.0, box_lengths=None):
        """                                                                  
        DESCRIPTION                                                          
                Given selection, draw the bounding box around it.          

        USAGE:
                drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]

        PARAMETERS:
                selection,              the selection to enboxen.  :-)
                                        defaults to (all)
   
                padding,                defaults to 0

                linewidth,              width of box lines
                                        defaults to 2.0

                r,                      red color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                               

                g,                      green color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                                 

                b,                      blue color component, valid range is [0.0, 1.0]
                                        defaults to 1.0                                

        RETURNS
                string, the name of the CGO box

        NOTES
                * This function creates a randomly named CGO box that minimally spans the protein. The
                user can specify the width of the lines, the padding and also the color.                            
        """                                                                                                    

        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection) if box_lengths is None else ([0,0,0], box_lengths)

        print "Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)

        minX = minX - float(padding)
        minY = minY - float(padding)
        minZ = minZ - float(padding)
        maxX = maxX + float(padding)
        maxY = maxY + float(padding)
        maxZ = maxZ + float(padding)

        if padding != 0:
            print "Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ)

        boundingBox = [
                LINEWIDTH, float(linewidth),

                BEGIN, LINES,
                COLOR, float(r), float(g), float(b),

                VERTEX, minX, minY, minZ,       #1
                VERTEX, minX, minY, maxZ,       #2

                VERTEX, minX, maxY, minZ,       #3
                VERTEX, minX, maxY, maxZ,       #4

                VERTEX, maxX, minY, minZ,       #5
                VERTEX, maxX, minY, maxZ,       #6

                VERTEX, maxX, maxY, minZ,       #7
                VERTEX, maxX, maxY, maxZ,       #8


                VERTEX, minX, minY, minZ,       #1
                VERTEX, maxX, minY, minZ,       #5

                VERTEX, minX, maxY, minZ,       #3
                VERTEX, maxX, maxY, minZ,       #7

                VERTEX, minX, maxY, maxZ,       #4
                VERTEX, maxX, maxY, maxZ,       #8

                VERTEX, minX, minY, maxZ,       #2
                VERTEX, maxX, minY, maxZ,       #6


                VERTEX, minX, minY, minZ,       #1
                VERTEX, minX, maxY, minZ,       #3

                VERTEX, maxX, minY, minZ,       #5
                VERTEX, maxX, maxY, minZ,       #7

                VERTEX, minX, minY, maxZ,       #2
                VERTEX, minX, maxY, maxZ,       #4

                VERTEX, maxX, minY, maxZ,       #6
                VERTEX, maxX, maxY, maxZ,       #8

                END
        ]

        boxName = "box_" + str(randint(0,10000))
        while boxName in cmd.get_names():
                boxName = "box_" + str(randint(0,10000))

        cmd.load_cgo(boundingBox,boxName)
        return boxName

def save_png(model_name):
    print "saving png"
    filename="{0}.png".format(model_name)
    if RENDER_PROPERLY:
        cmd.set("antialias", ANTIALIAS)
        cmd.png(filename, ray=1, width=2400, height=1800, dpi=DPI)
    else:
        cmd.png(filename)
    # Get out!

def close():
    cmd.quit()

def load_model(model_name):
    cmd.delete("all")
    print "loading model"
    cmd.load("{0}.pdb".format(model_name))
    cmd.set("depth_cue",0)
    cmd.set("dash_gap", 0)
    cmd.set("dash_radius", 0.8)
    cmd.set("dash_color", "blue")
    #cmd.bg_color("white")
    cmd.color("green", "name C*")
    cmd.set("ray_opaque_background", 0)

def init():
    pymol.finish_launching()

def draw_connectivities(model_name, points, edges, cutoff_distance, box_lengths):

    load_model(model_name)
    draw_bounding_box(box_lengths=box_lengths)
    print "creating com dummy atoms"
    for i, pos in enumerate(points):
        pseudo_atom_name = "com_{0}".format(i)
        cmd.pseudoatom(pseudo_atom_name, pos=[p*10 for p in pos])
        cmd.hide("nonbonded", pseudo_atom_name)
    print "adding com distances"
    for i, e in enumerate(edges):
        distance_name = "distance_{0}".format(i)
        cmd.distance(distance_name, "com_{0}".format(e[0]), "com_{0}".format(e[1]))
        cmd.hide("labels", distance_name)
    cmd.show("sticks", model_name)
    cmd.set("stick_transparency", 1)
    cmd.set_view( (\
     0.459229559,    0.209184960,   -0.863334477,\
    -0.888164818,    0.090163931,   -0.450589985,\
    -0.0,    0.973708749,    0.227197826,\
    -0.000013376,   -0.000241756, -409.681274414,\
    46.355365753,   37.701179504,   50.704589844,\
   197.865356445,  621.506164551,   20.000000000 ))
    save_png("images/" + model_name + "_{0}".format(cutoff_distance))


if __name__=="__main__":
    # Load Structures

    cmd.load("whole_n15.pdb")
    draw_connectivities("whole_n15", [[10,5,5], [10,5,6]], [[0,1]])