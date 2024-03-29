Data1 = "/Users/hiroki/Work/DCPAM/hs94/last100days/Temp.nc@Temp,time=^-1"


options = 
['--var',     #                            | set the variable name and slicing parameters.
'--vf',      #                             | use variable whose name is included in the filename.
#'--pry',     #                             | start pry after drawing the first figure
#'--auto_pry',# <args>                      | automatically execute the given scripts like as in pry after drawing the first figure.
'--list',    #                             | execute gplist instead of any other options.
'--silent',  #                             | show only calculated values. do not show any messages.
'--parallel',#                             | use parallel anyway.
'--help',    #
###
### visualize options ###
###
'--clrmap 5 ',  # <1- or filename>            | set colormap to draw tone/contour. filename of DCL colormap is acceptable.
'--itr 30 ',     # <1-4,5-7,10-15,20-23,30-33> | set axis scale. default is 1.
                #                             | 1 : linear scale for x/y axis
                #                             | 2 : linear scale for x , log scale for y axis
                #                             | 3 : log scale for x , linear scale for y axis
                #                             | 4 : log scale for x/y axis
#'--similar', # <simfac,vxoff,vyoff>        | (for 5<=itr<=7) set similarity parameters which are fed in DCL.grssim.
'--map_axis 180,0,23.4 --itr 30',# <uxc,uyc,rot>               | (for 10<=itr<=33) set mapping parameters which are fed in DCL.umpcnt.
'--map_radius 60 --itr 30', # <radius>                 | (for itr>=20) set clipping radius (degree) around the tangential point. Deafault=90.
#'--xcoord',  # <xcoord>                    | name of x-coordinate (use for associate coordinates)
#'--ycoord',  # <ycoord>                    | name of y-coordinate (use for associate coordinates)
'--title Title',   # <title>                     | set title of figure
'--notitle', #                             | do not print title. equivalent to --title ""
'--aspect 1',  # <aspect>                    | set aspect ratio of Viewport. default is 2.0.
#'--anim', '--animate', # <dim>             | plot animation along <dim>. <dim> must be name of dimension.
'--noannotate', #                          | not draw annotations.
#'--alternate', '--Ga', #                   | enable to backing store.
#'--nowait', '--Gw', #                      | not wait for any actions if animate
#'--smooth', '--Gaw', #                     | equal to --anlternate && --nowait
#'--reverse', '--Gr', #                     | plot animation reversible if animate
'--exch', #                                | exchange(transpose) x/y axis.
'--map coast_world --itr 10', # <map_type>               | plot map. itr number must be set. this option is neglect if itr number is 1-4.
                   #                          | abailable map type is coast_world, border_world, plate_world, state_usa, coast_japan, pref_japan
#'--time_ax',    # <nil|false|h|ymd>        | specify type of calendar-type time axis:
                   #                          |  nil   (=> auto slection)
                   #                          | false (=> do not use the time axis even if
                   #                          |           the units of the axis is a time one with since field)
                   #                          | "h"   (=> like nil, but always use the hour-resolving datetime_ax method
                   #                          |           in dclext_datetime_ax.rb)
                   #                          | "ymd" (=> like "h" but for y-m-d type using DCL.uc[xy]acl)
'--sldiv y,2,2',      # <<y|t>,m,n>              | split the drawing window into multiple panel.
                   #                          | argument should be provied in DCL.sldiv form.
                   #                          | e.g., --sldiv y,2,2 make 2x2 window.
]
=begin
'--scatter',    # [xtitle,ytitle]          | make scatter plot with provided 2 gturls.
                   #                          | 1st and 2nd gturls will be used for x- and y-axes, respectively.
                   #                          | titles of x- and y- axes can be specified by arguments xtitle,ytitle (optional)
                   #                          | range of each axis can be specified by option --range [xmin:xmax,ymin:ymax] format.
'--color_scatter', # [xtitle,ytitle]       | make color scatter plot with provided 3 gturls.
                      #                       | 1st and 2nd gturls are used for x- and y-axes, respectively.
                      #                       | 3rd gturl is used for color of marks.
                      #                       | titles of x- and y- axes can be specified by arguments xtitle,ytitle (optional)
                      #                       | range of each axis and color can be specified by option --range [xmin:xmax,ymin:ymax,zmin:zmax] format
'--cot',           #                       | contour over tone. draw tone with 1st gturl and overplot contour with 2nd gturl.
'--histogram',     # [any]                 | draw 1D histogram.
                      #                       | options "--int -n" and "--range xmin:xmax[,ymin:ymax]" can be used to set bins number and range.
                      #                       | options "--exch", "--title", and "--overplot n" also can be used.
                      #                       | if any argument is given, histogram are shown in percentage (%) not in numbers.
'--histogram2D',   # [any]                 | draw 2D histogram with two gtulrs.
                      #                       | if any argument is given, histogram are shown in percentage (%) not in numbers.
                      #                       | options "--int -n" and "--range xmin:xmax[,ymin:ymax]" can be used to set bins number and range.
                      #                       | options "--sint [-]n" and "--srange xmin:xmax" can be used to set tone interval[number] and range.
                      #                       | options "--exch" and "--title" also can be used.
'--rmap',          # <dim>                 | draw 2D map of correaltion between 1st and 2nd gturls along <dim> axis.
                      #                       | <dim> can be axis name (string) or dimension number (integer), but only one
                      #                       | axis is acceptable.
'--linearline',    # <arg>                 | draw linear line expressed bay <arg>.
                      #                       | <arg> should be "x=<value>", "y=<value>", "x=y", or "y=x", where <value> is integer or float.
                      #                       | multiple lines can be drawn by setting multiple args separated by comma, such as "x=0,y=0".
'--land',          #                       | mask data with land info, https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc.
                      #                       | you need to download landsea.nc and place it on correct dir.
'--ocean',         #                       | mask data with ocean info, https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc.
                      #                       | you need to download landsea.nc and place it on correct dir.
'--line',          #                       | make line plot forced. (about first 1D)
'--mark',          #                       | make mark plot forced. (about first 1D)
'--index',         # <index_num>           | set DCL line index, which set the color/thickness of the line primitive. please see DCL documents.
'--type',          # <type_num>            | set line type.
'--overplot',      # <num>                 | set number of lines on each figure with color.
                      #                       | use with --nocolor if you want to distinguish by line type.
'--nocolor',       #                       | draw without using color.
'--right_axis',    #                       | use with "--overplot 2" and draw 2nd y-axis in right side for 2nd gturl.
                      #                       | range of 2nd y-axis can be specified by --range left_min1:left_max,right_min:right_max.
'--top_axis',      #                       | 使い方を忘れた...
'--overplot_rm',   # <span>                | overplot the running mean with thick line. running mean span is required.
'--overplot_stddev',#                      | overplot the mean +/- stddev (like errorbars). this must be used with "--mean" option.
  ### tone or cont option ###
'--nocont',        #                       | make tone plot, without contour.
'--noshade',       #                       | make contour plot, without tone.
'--log_int',       #                       | use log interval in contour/tone.
'--range',         # <min:max>             | set min/max value for contour/tone/line/mark plot. min or max must be set.
'--crange',        # <min:max>             | set min/max value for contour plot. this is more dominant than --range
'--srange',        # <min:max>             | set min/max value for tone plot. this is more dominant than --interval/int
'--interval', '--int',  # <num>            | set interval value for contour/tone plot. set the number of lines if you set negative value.
'--cint',          # <num>                 | set interval value for contour plot. this is more dominant than --interval/int
'--sint',          # <num>                 | set interval value for tone plot. this is more dominant than --interval/int.
'--levels',        # <val1,val2,val3,...>  | set values of contour/tone levels.
'--clevels',       # <val1,val2,val3,...>  | set values of contour levels.
'--slevels',       # <val1,val2,val3,...>  | set values of tone levels.
'--patterns',      # <pattern1,pattern2,..>| set each patterns for tone plot.
'--tone',          # <a|e|f|b|c>           | set tone subroutine:
                      #                       | a (=> tone routine is selected automatically depending on the datasize)
                      #                       | e (=> DCL.uetone is used)
                      #                       | f (=> DCL.uetonf is used)
                      #                       | b (=> DCL.uetonb is used)
                      #                       | c (=> DCL.uetonc is used)
'--udsfmt',        # <strings>             | change contour label format. see UDCNTR/DCL manual for the format.
'--nocolorbar',    #                       | do not draw color bar.
'--nozero',        #                       | do not draw zero contour.
'--nodraw',        #                       | do not draw any figures.
'--file',       # <png|eps|pdf>[,filename] | raw in file with given format. available formats are png, eps, and pdf.
'--crop',       #                          | crop (trim) output png or pdf.
'--title_array',# <string,string,...>      | set multiple titles given as --title_array hoge,foo,bar; this gives the title "hoge" to
                   #                          | the 1st figure, "foo" to the 2nd one, and "bar" to the 3rd one.
'--miscindex',  # <num>                    | set line index (color and width) of title, label, frame, etc.
'--index_array',# <index1,index2,...>      | set array of index used for 1st, 2nd,... lines/mark.
'--type_array', # <index1,index2,...>      | set array of type index used for 1st, 2nd,... lines/mark.
'--subtitle',   # <string>                 | set subtitle, which locates just below the title and is common for multple figures.
'--irange',     # <min:max>                | same as --range but with +/- infinity for both boundaries. use with --int is recomended.
'--xrange',     # <min:max>                | set the range of x-axis for line plot and 2D plot.
'--yrange',     # <min:max>                | set the range of y-axis for line plot and 2D plot.
'--clr_range',  # <min:max>                | set color range. maximum range is 10:99. default is 15:94.
'--wm',         #                          | alias for --itr 10 --map coast_world --nocont.
'--size',       # <num>                    | set the size factor of x-window (default: 1.25)
'--fullcolor',  # <1|2>                    | draw fullcolor fig. --fullcolor 1 draws with 1 var.
                   #                          | --fullcolor 2 draws with 2 vars; 1st gturl is used for color (hue) and
                   #                          | 2nd one is used for brightness (value).
                   #                          | option of --range xmin:xmax,ymin:ymax can be used.
'--step',       #                          | make 1D line plot to step shape.
'--zerocenter', # [max]                    | set y-axis center in line plot figure to zero.
'--vector',     #                          | draw a vector plot.
                   #                          | When 2 gturls are given, 1st one and 2nd one are used for x- and y-components of the vectors.
                   #                          | When 3 or 4 gturls are given, 2nd one and 3rd one are used for x- and y-components of the vectors,
                   #                          | and tone of 1st one and counter of 4th one are overlaid.
                   #                          | If the 1st dim is lon or lat and the 2nd dim is not lon or lat, vectors are scaled by the geometry aspect ratio.
'--vfact',      # <factor>                 | use with "--vector" to change the size of the vector by <factor>.
'--vint',       # <int>                    | use with "--vector" to change the gird-interval to draw vectors.
'--vkeep',      #                          | use with "--vector" to keep the size of unit vector as that of the first plot.
'--nolegend',   #                          | do not draw legend of line plot.
'--textbox',    # <text>,<num>,[l|c]       | draw textbox, legends, and/or color bar in frame #<num>.
'--fact',       # <fact>                   | factor for charcter size
'--xmean',      #                          | draw line plot of mean along x-axis on the right hand side of the tone/contour plot.
'--ymean',      #                          | draw line plot of mean along y-axis below the tone/contour plot.
'--maskshading',# <gturl,maskval,ipat>     | draw a mask for values less than or equal to maskval (optional; defalult is 0)
                   #                          | with given gturl (optional; default is same as previously drawn GPhys)
                   #                          | by ipat pattern number (optional; defalult is 1515).
'--nolabels',   #                          | do not draw labels.
'--fill',       #                          | fill regions between plotted line and zero (or min or max) line. if histogram, fill the boxes.
'--lwf',        # <float>                  | set line with factor for pdf drawing. defalut is 999, which looks like DISP drawing.
'--axis_options',# <var1=val1,var2=val2>,..| set some axis options of USPACK (i.e., xoff/yoff, xfac/yfac, dxt/dyt, dxl,dyl)
'--bothsides',   #                         | draw bothsides of the planet when itr = 30
'--xint',  # <xlabelint[:division]>        | set x-axis' label-tickmark interval and number of number of division between two label-tickmarks
              #                               | e.g., --xint 10:5 draws label interval with 10 and sub-tickmark interval of 2 (=10/5).
'--yint',  # <ylabelint[:division]>        | same as --xint but for y-axis.
'--panelfit', # [margin]                   | force the figure to fit into the panel. margin can be given by ratio [0-1].
###
### analysis options ###
###
'--mean',   # <dim>                        | mean along axis <dim>.
'--stddev', # <dim>                        | standard deviation along axis <dim>.
'--eddy',   # <dim>                        | deviation from mean along axis <dim>.
'--diff',   #                              | calculate 1st gturl - 2nd gturl
'--divide', #                              | calculate 1st gturl / 2nd gturl
'--add',    #                              | calculate the sum of multiple variables.
'--operation',  # <math_func>              | operation of the specified math function on the data.
                   #                          | <math_func> should be a math function with one argument such as log10, sqrt, sin, etc.
'--running_mean','--rm', # <dim>:<span>[:<skip>][,<dim>:<span>[:<skip>]]
                            #                 | calculate the running mean along <dim> with a window of <span>.
                            #                 | multiple dims are not supported.
                            #                 | if integer <skip> is given, the date thinning is performed to the running mean.
                            #                 | for example, --rm t,12,12 to monthly data gives annual mean series (one value per year).
'--global_mean','--gm',  #                 | calculate the global mean considering the weight along latitude.
                            #                 | dims of longitude and latitude must be the 1st and 2nd dims, respectivly.
'--ensemble_mean','--em', # [num]          | calculate the ensemble mean.
                             #                | you must provide two or more gturls whose size and dimension are same.
                             #                | if <num> is given, ensemble mean will be calculated for each <num> gturls.
                             #                | if "concat" is given for <num>, provided gturls are concatenated along new dimension "member".
'--normalize',           # [ |mean|float]  | normalize data by initial (when no arg.) or mean or given value when data is one-dimension.
                            #                 | normalization for other data is not implemented yet.
'--lowpass','--lp',      # <dim>,<wn>      | apply lowpass filter to <dim>-th dimension with wavenumber <wn> using fft.
                            #                 | i.e., this cut off the wavenumbers greater than <wn>.
                            #                 | <dim> must be integer. dimname is not supported yet.
'--highpass','--hp',     # <dim>,<wn>      | apply highpass filter to <dim>-th dimension with wavenumber <wn> using fft.
                            #                 | i.e., this cut off the wavenumbers less than <wn>.
                            #                 | <dim> must be integer 0, 1, 2, or 3. dimnames and dim > 3 are not supported yet.
'--power_spectra','--ps',# <dim>[,any]     | calculate power spectra using fft along <dim>-th dimension.
                            #                 | by default, wavenumber is used for axis.
                            #                 | if any argument is given after comma, wavelength or period is used for axis.
'--derivative',          # <dim>           | calculate the derivative alog axis <dim>. GPhys.threepoint_O2nd_deriv is used.
'--srot',                #                 | calculate rotation on sphere with 1st gturl and 2nd gturl. radius can be set by --radius.
'--sdiv',                #                 | calculate divergence on sphere with 1st gturl and 2nd gturl. radius can be set by --radius.
'--sht',                 #                 | use spherical_harmonics_next module for srot, sdiv.
'--HKE',                 #                 | calculate horizontal Kinetic energy with 1st gturl and 2nd gturl.
'--integrate', # <dims>                    | calculate integration along <dims>.
'--sum',       # <dims>                    | calculate sum along <dims>.
'--addingup',  # <dim>[,<dim>,...]         | adding up values along <dim>; i.e., val[i] =  val[0..i].sum
'--max',       # <dims>                    | max values along axis <dims>.
'--min',       # <dims>                    | min values along axis <dims>.
'--ES',        # [total|rot|div][,vor_div_given]
                  #                           | calculate Horizontal Kinetic Energy Spectra from U given as 1st gturl
                  #                           | and V given as 2nd gturl. output total, rot, and/or div component.
                  #                           | if "vor_div_given" is given, 1st and 2nd gturls must be Vor and Div, respectivly,
                  #                           | and they will be used for calculation (faster).
'--wnf_analysis',# <sym|asym|bg|sym/bg|asym/bg|log_sym/bg|log_asym/bg>
                    #                         | calculate time-space spectra following Wheeler and Kiladis (1999)
                    #                         | name of component to be output has to be given.
'--mvo',  # <mathematical operation>       | mathematical operation for multiple gturls. For mathematical operation,
             #                                | 1st gturl must be indicated by "x", 2nd by "y", 3rd by "z", 4th by "u",
             #                                | 5th by "v", 6th by "w". This cannot treat vars more than 6.
             #                                | Also, this cannot be used repeatedly.
'--mvo_only',   #                          | use with "--mvo" to execute "multi variable operation" only.





###
### data/constants arrange options ###
###
'--cut',# dimname=pos1[:pos2[:thinning_intv]][,dimname=...]: | apply cut, slice, and/or thinning after mathmatical operation.
'--radius',  # <raidus>                    | set planetary radius.
'--Omega',   # <Omega>                     | set planetary rotation rate [s-1].
'--rename',  # <varname>                   | rename the output variable name by <varname>. if --long_name is not used, "long_name" is also renamed.
'--long_name', # <long_name>               | rename the output variable "long_name" by <long_name>.
'--regrid',  # <latgrid_num>               | interpolate to the Gausian latitude and cooresponding longitude.
                #                             | number of lat-grid points is required and available numbers are
                #                             | 8, 16, 32, 64, 94, 96, 120, 128, 160, 240, 256, 320, 480, 640. gglat.rb required [TODO: 自前計算にする]
'--GLAT',
'--read_csv',  # <num>                     | read csv file. the first column is used for the axis and the other columns are used for data values.
                  #                           | <num> specifies how many data columns you want to read.
'--composite',  # <dim>,<span>,<skip>      | if <skip> is not given or zero, calculate mean with interval <span>.
                   #                          | i.e., --composite time,12 calculates mean of each month from monthly mean data.
                   #                          | if <skip> is given, calculate composite mean for dimension <dim>
                   #                          | with length <span> with interval of <skip>.
                   #                          | i.e., --composite time,3,9 calculates mean for time=0,1,2,12,13,14,24,...
                   #                          | this can be used to calculate seasonal mean such as DJF from mounthly mean data.
'--split_axis',# <dim>,<span>,[newaxisname]| split the <dim> axis with <span>.
                  #                           | for example, --split_axis time,12 applied to monthly data splits the time axis
                  #                           | to year and month axis. name of the new axis can be given by [newnaxisname]
'--calendar360',  #                        | set the calendar type to "360_day calendar". normally, this will be set automatically.
'--calendar365',  #                        | set the calendar type to "365_day calendar". normally, this will be set automatically.
'--flip_data',    # <dim[,dim,..]>         | flip (turn over) the data and axis of <dim>th dimension. <dim> must by given by integer(s).
'--interpolate_pole', #                    | extend the latitude axis to have grid points on each pole (-90/+90 deg), and
                         #                    | interpolate the pole values by mean of values on surrounding grid points.
                         #                    | this assumes the data has lon and lat in its first and second dimension, respectivly.
'--set_missing_value',# [float]            | replace the missing values to given value and validate then. default: 0.0.
'--ignoreunit',       #                    | ignore unit
'--set_assoc_coords', # <gturl>            | set <gturl> as an associate coordinate.
'--sig2p', # <gturl of surface pressure>   | transform sigma coordinate to pressure coordinate. gturl of surface pressure must be given.
'--sig2z', # <gturl of Temperature>        | transform sigma coordinate to z (not log_p) coordinate. gturl of Temperature must be given.
'--interpolate', # <args>                  |   interpolate after operations.
                    #                         | <args> can be axis_name=90 or axis_name=[0,90,180,270] or axis_name=0:360:4 for example.
                    #                         | 1) axis_name=90 style gives interpolation at the location and rank shrinks.
                    #                         | 2) axis_name=[0,90,180,270] style gives interpolation for multiple location and use then as an axis.
                    #                         | 3) axis_name=0:360:4 style similar to (2) but 0:360:4 means from 0 to 360 with separation of (360-0)/4.
                    #                         | multiple interpolation is not supported yet.
'--extrapolation',   #                     | enable extrapolation, when using --interpolate option.
'--unit',        # <unit>                  | set the operated variable's unit.
'--invalidation',# <num or range>          | invalidate the given value or range.
'--join',        #                         | join muliple gturls along with an axis whose values are continual across the gturls and are not duplicate.
'--extend_backward_along', # <dim>,<len>[,<value>] | extend backward the object along <dim>-axis for <len> numbers with <value>.
'--thinning', # <num>                      | data thinning by picking up data with <num> interval
'--nothinning', #                          | prevent auto data thinning.
'--read_fits',  #                          | read fits format image. automatically added if file extension is .fits. (RubyFits required)
'--offset',     # <float>                  | offset the data by <float>.
'--cyclic_extension',# <dim:num[,dim:num]> | extend the data along <dim> with <num> number cyclically. e.g., --cyclic_extension x:4,y:4.
'--axis_units_name', # <axdim:unit[:axname:operation]>[,axdim:unit:axname:operation]
                        #                     | change unit (and name) of axis. multiple axes can be given as
                        #                     | --axis_units_name ax0:deg:lon,ax2:km:height:/1000  for example.
                        #                     | mathematical operation can be applied to the axis values.
'--replace_axis',#<axN:val1,val2,val3,...> | replace the value of the axis <axN> (N is 0, 1, 2,...) with given Array [val1, val2, val3, val4,...]
'--zonal_shift', #<speed at equator in m/s>| shift entire field zonally as it is advected by the solid body rotation of <speed at equator in m/s>.
###
### output options ###
###
'--nc',    # [outfilename]                 | output operated GPhys object as NetCDF file. output file name can be specified.
              #                               | Default name is "out.nc".
              #                               | if used with "--anim", GPhys object for last shown figure is outputtedd.
'--nc4',   # [outfilename]                 | output operated GPhys object as NetCDF ver 4 with compression.
              #                               | may take longer time to write out.
              #                               | output file name can be specified. Default name is "out.nc".
'--nc4a',  # [outfilename]                 | output operated gphys continually with along time dim.
'--csv',   # [outfilename]                 | output operated GPhys object (1 dim only) as CSV file. output file name can be specified.
              #                               | Default name is "out.csv".
              #                               | if used with "--anim", GPhys object for last shown figure is outputtedd.
'--produce', # <something>                 | produce <something>. <something> can be "solar_zenith_angle"
'--output_assoc_coords',  #                | use with --nc4 or --nc to output associate coordinate as netcdf file.
'--gif',   #                               | output gif animation. this must be used with --anim and --file png options.
'--mp4',   # [number of frames per second] | output mp4 animation. frame rate can be given. this must be used with --anim and --file png options.
'--ntype', # <type>                  
=end         


num = options.length
num.times{|i|
  n = format('%03d', i)
  cmd = "gpv #{Data1} --file png,#{n} #{options[i]} --title '#{options[i]}'"
  system(cmd)
}