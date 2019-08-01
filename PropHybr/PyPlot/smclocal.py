"""
##  The smclocal function plots a local region of a given smc grid.
##  cel is the full cell array in shape(nc, 5).
##                                JGLi04Mar2019
"""

def smclocal(cel, verts, ncels, colrs, config, Arctic=False, 
             mdlname='SMC', buoys='', psfile='output.ps', 
             paprorn='portrait', paprtyp='a3'):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readcell import readtext
    from steromap import steromap
    from colrboxy import colrboxy
    from scale_depth import scale_depth

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]
    if( Arctic ): ncabgm=config[4]

##  Maximum mapping radius.
    radius=rdpols[0]
    pangle=rdpols[1]
    plon  =rdpols[2]
    plat  =rdpols[3]

##  Define depth to color index conversion parameters and marks.
##  Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=131)

##  Cell color is decided by its depth value.
    ndeps = ncstr + np.rint( (cstar-np.log10(cel[:,4]))*factr ).astype(np.int)
    nc = cel.shape[0]
    if( Arctic ):
        na=int(ncabgm[1])
        nb=int(ncabgm[2])
        jm=int(ncabgm[3])

##  Workout output file format from its extension
#   for i in np.arange(len(psfile))[::-1]:
#       if( psfile[i] == '.' ):
#           psfmt = psfile[i+1:]
#           break
#   print " Output file format will be "+psfmt
##  Python plt.savefig will do this automatically.

##  Use selected cells to draw the plot.
    fig=plt.figure(figsize=sztpxy[0:2])
    ax=fig.add_subplot(1,1,1)
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    ax.set(**xprop)
    ax.set(**yprop)
    ax.set_aspect('equal')
    ax.set_autoscale_on(False)
    ax.set_axis_off()
    plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

##  Create color array for this plot
    pcface = []
    pcedge = []
    for i in ncels:
##  Mark the arctical map-east reference region by first ring  
##  of its boundary cells at jm-3
        if( Arctic and cel[i,1] == jm-24 ):
            pcface.append(colrs(246))
            pcedge.append(colrs(246))
        elif( Arctic and cel[i,1] == jm ):
            pcface.append(colrs(168))
            pcedge.append(colrs(168))
        else:
            pcface.append(colrs(255))
            pcedge.append(colrs(ndeps[i]))

##  Create PolyCollection from selected verts and define edge and face color.
    polynorth = PolyCollection(verts)
    polynorth.set_facecolor(pcface)
    polynorth.set_edgecolor(pcedge)
    polynorth.set_linewidth( 0.2 )

##  Draw the selected cells as colored polygons.
    ax.add_collection(polynorth)

##  Draw colorbar inside plot.
    xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks, nclrm=nclrm)
    ax.add_collection(clrply)
    dkx=clrbxy[2]; dky=clrbxy[3]
    for i in range(len(depth)):
        m = marks[i]
        if( dkx < dky ):    
            plt.text(xkeys[0]+1.15*dkx, ykeys[m], str(depth[i]),
                verticalalignment='center', fontsize=11, color='b' )
        else:
            plt.text(xkeys[m], ykeys[0]+1.15*dky, str(depth[i]),
                horizontalalignment='center', fontsize=11, color='b' )

    if( dkx < dky ):    
        plt.text(xkeys[0]+1.9*dkx, ykeys[marks[2]], 'Depth m',
                 rotation=-90,verticalalignment='center', fontsize=15, color='k' )
    else:
        plt.text(xkeys[marks[2]], ykeys[0]+1.9*dky, 'Depth m',
                 rotation=0,horizontalalignment='center', fontsize=15, color='k' )

##  Put cell information inside plot
    tpx=sztpxy[2] 
    tpy=sztpxy[3]
    dpy=-0.6 

    plt.text(tpx, tpy+dpy*1, mdlname+' Grid',  
             horizontalalignment='center', fontsize=15, color='k' )
    plt.text(tpx, tpy+dpy*2, 'NC='+str(nc), 
             horizontalalignment='center', fontsize=13, color='r' )
    if( Arctic ):
        plt.text(tpx, tpy+dpy*3, 'NA='+str(na), 
             horizontalalignment='center', fontsize=13, color='r' )
        plt.text(tpx, tpy+dpy*4, 'NB='+str(nb), 
             horizontalalignment='center', fontsize=13, color='r' )

#;  Overlay buoy sits on grid map if buoy file is provided.
    if( len(buoys) > 3 ):
        hdr, buoyll = readtext(buoys)
        nmbu=int(hdr[0])
        buoyids=buoyll[:,0].astype(str)
        buoylat=buoyll[:,1].astype(np.float)
        buoylon=buoyll[:,2].astype(np.float)

#; Convert slat slon to elat elon with given new pole
        elat,elon,sxc,syc = steromap(buoylat,buoylon,plat,plon,Pangl=pangle)

#; Mark buoy position on map
        print (' Selected buoys in this plot:')
        for i in range(nmbu):
            if( (elat[i] >= 25.0) and (rngsxy[0] < sxc[i] < rngsxy[1])
                              and (rngsxy[2] < syc[i] < rngsxy[3]) ):
                print (' {:6} {:8.3f} {:8.3f}'.format( buoyids[i], buoylat[i], buoylon[i] ))
                txtsz=int( abs(np.sin(elat[i]*d2rad)*12.0) )
                plt.text(sxc[i], syc[i], 'r',  fontsize=txtsz,
                     horizontalalignment='center', color='r' )
                plt.text(sxc[i], syc[i], '.',  fontsize=txtsz*2,
                     horizontalalignment='center', color='r' )

##  Save plot as ps file
    print (" Save the smc grid local plot ... " )
    plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
                orientation=paprorn,papertype=paprtyp)

##  End of smclocal plot function. ##

