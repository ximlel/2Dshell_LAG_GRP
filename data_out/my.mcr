#!MC 1120
## Set up Export file type and file name.
$!EXPORTSETUP EXPORTFORMAT = MPEG4
$!EXPORTSETUP EXPORTFNAME = "/home/leixin/Desktop/RESULT/timeseries.mp4"
$!EXPORTSETUP FFMPEGQSCALE = 1
$!EXPORTSETUP QUALITY = 100
$!EXPORTSETUP SUPERSAMPLEFACTOR = 16
$!EXPORTSETUP ANIMATIONSPEED = 10
$!EXPORTSTART
EXPORTREGION = CURRENTFRAME
$!LOOP 229
$!GLOBALTIME SolutionTime = (0.001*|LOOP|)
$!GLOBALCONTOUR COLORMAPFILTER
{
COLORMAPDISTRIBUTION = CONTINUOUS
CONTINUOUSCOLOR
{
CMIN = |MINC|
CMAX = |MAXC|
}
}
$!EXPORTNEXTFRAME
$!ENDLOOP
$!EXPORTFINISH