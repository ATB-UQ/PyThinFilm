import __main__
from random import randint
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import pymol
from pymol import cmd
from pymol.cgo import LINEWIDTH,  BEGIN, LINES, COLOR, VERTEX, END

DPI = 300
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
        maxZ = 110

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
    cmd.set_view( (\
     0.822059155,    0.141536161,   -0.551531196,\
     0.011512177,    0.964284837,    0.264617831,\
     0.569285274,   -0.223881051,    0.791070879,\
     0.000000000,    0.000000000, -486.167602539,\
    47.096778870,   37.321388245,   54.275535583,\
   219.350936890,  752.984191895,  -20.000000000 ))
    filename="{0}.png".format(model_name)

    if RENDER_PROPERLY:
        cmd.set("antialias", ANTIALIAS)
        cmd.png(filename, ray=1, width=1200, height=900, dpi=DPI)
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
    cmd.bg_color("white")
    cmd.color("green", "name C*")

def init():
    pymol.finish_launching()

def draw_connectivities(model_name, points, edges, cutoff_distance, box_lengths):

    load_model(model_name)
    draw_bounding_box()
    print "creating com dummy atoms"
    for i, pos in enumerate(points):
        cmd.pseudoatom("com_{0}".format(i), pos=[p*10 for p in pos])
    print "adding com distances"
    for i, e in enumerate(edges):
        distance_name = "distance_{0}".format(i)
        cmd.distance(distance_name, "com_{0}".format(e[0]), "com_{0}".format(e[1]))
        cmd.hide("labels", distance_name)
    cmd.show("sticks", model_name)
    cmd.set("stick_transparency", 1)
    save_png("images/" + model_name + "_{0}".format(cutoff_distance))


if __name__=="__main__":
    # Load Structures

    cmd.load("whole_n15.pdb")
    draw_connectivities("whole_n15", [[10,5,5], [10,5,6]], [[0,1]])