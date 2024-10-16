#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.2 patchlevel 6 
#    	last modified Sep 2009
#    	System: Linux 2.6.32-55-generic
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004, 2007 - 2009
#    	Thomas Williams, Colin Kelley and many others
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from http://www.gnuplot.info/faq/
#    
#    	Send bug reports and suggestions to <http://sourceforge.net/projects/gnuplot>
#    
unset clip points
set clip one
unset clip two
set bar 1.000000
set xdata
set ydata
set zdata
set x2data
set y2data
set boxwidth
set style fill  empty border
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set angles radians
set grid nopolar
set grid xtics mxtics ytics mytics noztics nomztics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set key title ""
set key off
unset label
unset arrow
unset style line
unset style arrow
unset logscale
set logscale x 10
set logscale y 10
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1  
set samples 1000, 1000
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size 1,1
set origin 0,0
set style data lines
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics 10
set mytics 10
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics autofreq
set ytics autofreq
set ztics autofreq
set nox2tics
set noy2tics
set cbtics autofreq
set title "Empirical Roofline Graph (Results.snellius.genoa.02/Run.003)" 
set timestamp bottom 
set timestamp "" 
set rrange [ * : * ] noreverse nowriteback  # (currently [0.00000:10.0000] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set xlabel "FLOPs / Byte" 
set x2label "" 
set xrange [1.000000e-02 : 1.000000e+02] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set ylabel "GFLOPs / sec" 
set y2label "" 
set yrange [1.000000e+00 : *] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zlabel "" 
set zrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set cblabel "" 
set cbrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "C"
set pm3d explicit at s
set pm3d scansautomatic
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set loadpath 
set fit noerrorvariables

set term postscript solid color rounded
set output "Results.snellius.genoa.02/Run.003/roofline.ps"

# Plotting goes after this...
set label '76.3 GFLOPs/sec (FP64 Maximum)' at 2.0000000e+00,9.1560000e+01 left textcolor rgb '#000080'
set label 'L1 - 609.8 GB/s' at 2.8351540e-02,2.0919411e+01 left rotate by 45 textcolor rgb '#800000'
set label 'L2 - 460.6 GB/s' at 3.9471136e-02,1.5026107e+01 left rotate by 45 textcolor rgb '#800000'
set label 'L3 - 353.5 GB/s' at 4.5055597e-02,1.3163681e+01 left rotate by 45 textcolor rgb '#800000'
set label 'L4 - 241.8 GB/s' at 5.4475428e-02,1.0887432e+01 left rotate by 45 textcolor rgb '#800000'
set label 'L5 - 145.9 GB/s' at 7.0121902e-02,8.4580923e+00 left rotate by 45 textcolor rgb '#800000'
set label 'L6 - 82.6 GB/s' at 9.3221975e-02,6.3622072e+00 left rotate by 45 textcolor rgb '#800000'
set label 'DRAM - 44.3 GB/s' at 1.1224131e-01,4.1074753e+00 left rotate by 45 textcolor rgb '#800000'
plot \
     (x <= 1.2512299e-01 ? 6.0980000e+02 * x : 1/0) lc 1 lw 2,\
     (x <= 1.6564271e-01 ? 4.6063000e+02 * x : 1/0) lc 1 lw 2,\
     (x <= 2.1582937e-01 ? 3.5352000e+02 * x : 1/0) lc 1 lw 2,\
     (x <= 3.1551090e-01 ? 2.4183000e+02 * x : 1/0) lc 1 lw 2,\
     (x <= 5.2278177e-01 ? 1.4595000e+02 * x : 1/0) lc 1 lw 2,\
     (x <= 9.2395253e-01 ? 8.2580000e+01 * x : 1/0) lc 1 lw 2,\
     (x <= 1.7231256e+00 ? 4.4280000e+01 * x : 1/0) lc 1 lw 2,\
     (x >= 1.2512299e-01 ? 7.6300000e+01 : 1/0) lc 3 lw 2
