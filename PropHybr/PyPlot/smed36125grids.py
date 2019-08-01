"""
##
##  It reads cell files and uses own projection to draw 
##  the SMC grid. Projected polygons are collected into 
##  vert variables for grid and subsequent swh plots. 
##
#;  First created: For 4 views     J G Li   26 Nov 2008
#;  Modified for Color grids       J G Li   26 Aug 2009
#;  Adapted for 25km SMC grids     J G Li    3 Feb 2010
#;  Updated G25SMC-12c grids.      J G Li   25 May 2011
#;  Sterographic projection nR.    J G Li   30 Jun 2011
#;  Adapted for G50SMC grid.       J G Li   19 Aug 2011
#;  Extended to fill the Arctic.   J G Li    5 Oct 2011
#;  Rectify polar cell position.   J G Li   25 Oct 2011
#;  Simplify with readcell and steromap.  JGLi12Nov2014
##
##  Converted into a Python function.     JGLi05Dec2018
##  Save ELat/Lon and sx/yc in file.      JGLi11Dec2018
##  Add color map and draw color bar.     JGLi14Dec2018
##  Adapted for SMC36125 grid plot.       JGLi03Jan2019
##  Import resource and set stacksize.    JGLi07Jan2019
##  Use polycollections for two plots.    JGLi30Jan2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Updated to Python 3.6.8 GCC 7.3.0.    JGLi02Apr2019
##  Adapted for UK12H rotated SMC grid.   JGLi20Jun2019
##  Adapted for SMC Med36125 grid.        JGLi08Jul2019
##  Modified to run on xcs machine.       JGLi16Jul2019
##
"""

if( 1 > 0 ):
##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from datetime import datetime
    from readcell import readcell   
    from steromap import steromap
    from rgbcolor import rgbcolor

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Read global and Arctic part cells. 
    DatGMC='/data/d02/frjl/PropHybr/DatGMC/'
    Wrkdir='/data/d02/frjl/PropHybr/PyPlot/'
    MyCode='/data/d02/frjl/PropHybr/PyPlot/'
    Cel_file = DatGMC+'Med36125Cel0.dat'

    headrs, cel = readcell( [Cel_file] ) 
    nc = int( headrs[0].split()[0] )
    print (' Med36125 total cel number = %d' % nc )

##  Size-1 cell increments
#   dxlon=0.01350
#   dylat=0.01350
#   dxlon=0.703125
#   dylat=0.468750
#   dxlon=0.087890625
#   dylat=0.058593750
    dxlon=0.0439453125
    dylat=0.029296875

##  Lon and lat grid numbers + 1  for cell faces
#   nx=1459
#   ny=1345
#   nx=513
#   ny=385
#   nx=4097
#   ny=3073
    nx=8193
    ny=6145

##  Rotated grid starting i, j locations
#   ZLon=-10.88950 - dxlon*0.5
#   ZLat= -7.29420 - dylat*0.5
    imnrgn = np.min( cel[:,0] )
    ishft = -(int(imnrgn)/16 - 1)*16
    jequt = 3072
    print ("ishft, jequt =", ishft, jequt)
    ZLon= -0.0
    ZLat= -0.0
    print ("ZLon, ZLat =", ZLon, ZLat)

##  Half the latitude grid to set j=0 on the Equator.
    xlon=(np.arange(nx)-ishft)*dxlon + ZLon
    ylat=(np.arange(ny)-jequt)*dylat + ZLat

##  Adjust cel i,j by ishft,jequt for easy indexing in xlon, ylat
    cel[:,0]=cel[:,0]+ishft
    cel[:,1]=cel[:,1]+jequt

##  Use own color map and defined depth colors 
    colrfile = MyCode+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Maximum mapping radius.
    radius=10.0

    print (" Draw SMC Med36125 grid ... ")

##  Mediterranean regional plot
    pangle= 11.50 
    plon= 14.8
    plat= 38.0
    clrbxy=[  3.5,  5.0,   8.0,  0.8]
    sztpxy=[ 15.0,  9.0,  -9.0, -3.0]
    rngsxy=[-15.0, 15.0,  -9.0,  9.0]
    papror='landscape'

    print( " Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Initial verts and ncels variable for polycollections.
    nvrts = []
    ncels = []
##  Loop over all cells. 
    for i in range(nc): 
        xc=[cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0]]
        yc=[cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1]]
        slat=ylat[yc]
        slon=xlon[xc]

##  Convert slat slon to elat elon with given new pole
        elat,elon,sxc,syc = steromap(slat,slon,plat,plon,Pangl=pangle,Onecl=True)

        if( (elat[0] >= 0.0) and (rngsxy[0] < sxc[0] < rngsxy[1])
                             and (rngsxy[2] < syc[0] < rngsxy[3]) ):
                nvrts.append(list(zip(sxc,syc)))
                ncels.append(i)

##  Set plot size and limits and message out anchor point.
    rdpols=[radius, pangle, plon, plat]
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy])
    pzfile=DatGMC+'VrtsMed325.npz'

##  Store selected north and south verts and cell numbers for swh plots.
##  Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

##  These variables could be loaded back by
#   vrtcls = np.load(DatGMC+'VrtsUK12H.npz')
#   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
##

##  Draw your selected grid plot.
    psfile=Wrkdir+'sMed36125grd.ps' 

    from smclocal import smclocal
    smclocal( cel, nvrts,ncels,colrs,config, 
              mdlname='Med36125', psfile=psfile,
              paprorn=papror)

    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of smed36126grids program ##

