#!/usr/bin/env ruby
##################################################
=begin rdoc
= NAME
gpv - quick viewer and manipulater for the values of a variable specified by a gtool4-type URL.

gpv is based on and inspired by gpview, which is developed by Dr Shin-ichi Takehiro and other dcmodel developers.

= AUTHOR
Hiroki Kashimura (hiroki@gfd-dennou.org)

= USAGE
  gpv hoge.nc@hoge

  see http://www.gfd-dennou.org/member/hiroki/homepage/main007.html#gpv for details.
=end

#################################################

GPV_DIR = "/Users/hiroki/Dropbox/RubyScripts/gpv_dir/"

require GPV_DIR+"gpv_config.rb"
require GPV_DIR+"gpv_analysis.rb"
require GPV_DIR+"gpv_arrange.rb"
require GPV_DIR+"gpv_utils.rb"
require GPV_DIR+"gpv_visualize.rb"
require GPV_DIR+"gpv_gphysmod.rb"
require GPV_DIR+"gpv_spherical_harmonics_next.rb"



class GPV
#↑↑See this source code for available options↑↑
def __set_options(opts=nil)
  if (opts.class == String) then
    prev_argv = ARGV.clone
    while (ARGV[0])
      ARGV.shift
    end
    ARGV.concat(opts.split(" "))
  end
  ## options for DCL
  check_dclopts
  ## parse options
  parser = GetoptLong.new
  parser.set_options(




###
### global options ###
###
  ['--var',     #                             | set the variable name and slicing parameters.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--vf',      #                             | use variable whose name is included in the filename.
    GetoptLong::NO_ARGUMENT],
  ['--pry',     #                             | start pry after drawing the first figure
    GetoptLong::NO_ARGUMENT],
  ['--auto_pry',# <args>                      | automatically execute the given scripts like as in pry after drawing the first figure.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--list',    #                             | execute gplist instead of any other options.
    GetoptLong::NO_ARGUMENT],
  ['--silent',  #                             | show only calculated values. do not show any messages.
    GetoptLong::NO_ARGUMENT],
  ['--parallel',#                             | use parallel anyway.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--help',    #
    GetoptLong::NO_ARGUMENT],



###
### visualize options ###
###
  ['--clrmap',  # <1- or filename>            | set colormap to draw tone/contour. filename of DCL colormap is acceptable.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--itr',     # <1-4,5-7,10-15,20-23,30-33> | set axis scale. default is 1.
                #                             | 1 : linear scale for x/y axis
                #                             | 2 : linear scale for x , log scale for y axis
                #                             | 3 : log scale for x , linear scale for y axis
                #                             | 4 : log scale for x/y axis
    GetoptLong::REQUIRED_ARGUMENT],
  ['--similar', # <simfac,vxoff,vyoff>        | (for 5<=itr<=7) set similarity parameters which are fed in DCL.grssim.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--map_axis',# <uxc,uyc,rot>               | (for 10<=itr<=33) set mapping parameters which are fed in DCL.umpcnt.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--map_radius', # <radius>                 | (for itr>=20) set clipping radius (degree) around the tangential point. Deafault=90.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--xcoord',  # <xcoord>                    | name of x-coordinate (use for associate coordinates)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--ycoord',  # <ycoord>                    | name of y-coordinate (use for associate coordinates)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--title',   # <title>                     | set title of figure
    GetoptLong::REQUIRED_ARGUMENT],
  ['--notitle', #                             | do not print title. equivalent to --title ""
    GetoptLong::NO_ARGUMENT],
  ['--aspect',  # <aspect>                    | set aspect ratio of Viewport. default is 2.0.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--anim', '--animate', # <dim>             | plot animation along <dim>. <dim> must be name of dimension.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--noannotate', #                          | not draw annotations.
    GetoptLong::NO_ARGUMENT],
  ['--alternate', '--Ga', #                   | enable to backing store.
    GetoptLong::NO_ARGUMENT],
  ['--nowait', '--Gw', #                      | not wait for any actions if animate
    GetoptLong::NO_ARGUMENT],
  ['--smooth', '--Gaw', #                     | equal to --anlternate && --nowait
    GetoptLong::NO_ARGUMENT],
  ['--reverse', '--Gr', #                     | plot animation reversible if animate
    GetoptLong::NO_ARGUMENT],
  ['--exch', #                                | exchange(transpose) x/y axis.
    GetoptLong::NO_ARGUMENT],
  ['--map', '--m', # <map_type>               | plot map. itr number must be set. this option is neglect if itr number is 1-4.
                   #                          | abailable map type is coast_world, border_world, plate_world, state_usa, coast_japan, pref_japan
    GetoptLong::REQUIRED_ARGUMENT],
  ['--time_ax',    # <nil|false|h|ymd>        | specify type of calendar-type time axis:
                   #                          |  nil   (=> auto slection)
                   #                          | false (=> do not use the time axis even if
                   #                          |           the units of the axis is a time one with since field)
                   #                          | "h"   (=> like nil, but always use the hour-resolving datetime_ax method
                   #                          |           in dclext_datetime_ax.rb)
                   #                          | "ymd" (=> like "h" but for y-m-d type using DCL.uc[xy]acl)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--sldiv',      # <<y|t>,m,n>              | split the drawing window into multiple panel.
                   #                          | argument should be provied in DCL.sldiv form.
                   #                          | e.g., --sldiv y,2,2 make 2x2 window.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--scatter',    # [xtitle,ytitle]          | make scatter plot with provided 2 gturls.
                   #                          | 1st and 2nd gturls will be used for x- and y-axes, respectively.
                   #                          | titles of x- and y- axes can be specified by arguments xtitle,ytitle (optional)
                   #                          | range of each axis can be specified by option --range [xmin:xmax,ymin:ymax] format.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--color_scatter', # [xtitle,ytitle]       | make color scatter plot with provided 3 gturls.
                      #                       | 1st and 2nd gturls are used for x- and y-axes, respectively.
                      #                       | 3rd gturl is used for color of marks.
                      #                       | titles of x- and y- axes can be specified by arguments xtitle,ytitle (optional)
                      #                       | range of each axis and color can be specified by option --range [xmin:xmax,ymin:ymax,zmin:zmax] format
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--cot',           #                       | contour over tone. draw tone with 1st gturl and overplot contour with 2nd gturl.
    GetoptLong::NO_ARGUMENT],
  ['--histogram',     # [any]                 | draw 1D histogram.
                      #                       | options "--int -n" and "--range xmin:xmax[,ymin:ymax]" can be used to set bins number and range.
                      #                       | options "--exch", "--title", and "--overplot n" also can be used.
                      #                       | if any argument is given, histogram are shown in percentage (%) not in numbers.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--histogram2D',   # [any]                 | draw 2D histogram with two gtulrs.
                      #                       | if any argument is given, histogram are shown in percentage (%) not in numbers.
                      #                       | options "--int -n" and "--range xmin:xmax[,ymin:ymax]" can be used to set bins number and range.
                      #                       | options "--sint [-]n" and "--srange xmin:xmax" can be used to set tone interval[number] and range.
                      #                       | options "--exch" and "--title" also can be used.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--rmap',          # <dim>                 | draw 2D map of correaltion between 1st and 2nd gturls along <dim> axis.
                      #                       | <dim> can be axis name (string) or dimension number (integer), but only one
                      #                       | axis is acceptable.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--linearline',    # <arg>                 | draw linear line expressed bay <arg>.
                      #                       | <arg> should be "x=<value>", "y=<value>", "x=y", or "y=x", where <value> is integer or float.
                      #                       | multiple lines can be drawn by setting multiple args separated by comma, such as "x=0,y=0".
    GetoptLong::REQUIRED_ARGUMENT],
  ['--land',          #                       | mask data with land info, https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc.
                      #                       | you need to download landsea.nc and place it on correct dir.
    GetoptLong::NO_ARGUMENT],
  ['--ocean',         #                       | mask data with ocean info, https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc.
                      #                       | you need to download landsea.nc and place it on correct dir.
    GetoptLong::NO_ARGUMENT],
  ['--line',          #                       | make line plot forced. (about first 1D)
    GetoptLong::NO_ARGUMENT],
  ['--mark',          #                       | make mark plot forced. (about first 1D)
    GetoptLong::NO_ARGUMENT],
  ['--index',         # <index_num>           | set DCL line index, which set the color/thickness of the line primitive. please see DCL documents.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--type',          # <type_num>            | set line type.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--overplot',      # <num>                 | set number of lines on each figure with color.
                      #                       | use with --nocolor if you want to distinguish by line type.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--nocolor',       #                       | draw without using color.
    GetoptLong::NO_ARGUMENT],
  ['--right_axis',    #                       | use with "--overplot 2" and draw 2nd y-axis in right side for 2nd gturl.
                      #                       | range of 2nd y-axis can be specified by --range left_min1:left_max,right_min:right_max.
    GetoptLong::NO_ARGUMENT],
  ['--top_axis',      #                       | 使い方を忘れた...
    GetoptLong::REQUIRED_ARGUMENT],
  ['--overplot_rm',   # <span>                | overplot the running mean with thick line. running mean span is required.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--overplot_stddev',#                      | overplot the mean +/- stddev (like errorbars). this must be used with "--mean" option.
    GetoptLong::NO_ARGUMENT],
  ### tone or cont option ###
  ['--nocont',        #                       | make tone plot, without contour.
    GetoptLong::NO_ARGUMENT],
  ['--noshade',       #                       | make contour plot, without tone.
    GetoptLong::NO_ARGUMENT],
  ['--log_int',       #                       | use log interval in contour/tone.
    GetoptLong::NO_ARGUMENT],
  ['--range',         # <min:max>             | set min/max value for contour/tone/line/mark plot. min or max must be set.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--crange',        # <min:max>             | set min/max value for contour plot. this is more dominant than --range
    GetoptLong::REQUIRED_ARGUMENT],
  ['--srange',        # <min:max>             | set min/max value for tone plot. this is more dominant than --interval/int
    GetoptLong::REQUIRED_ARGUMENT],
  ['--interval', '--int',  # <num>            | set interval value for contour/tone plot. set the number of lines if you set negative value.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--cint',          # <num>                 | set interval value for contour plot. this is more dominant than --interval/int
    GetoptLong::REQUIRED_ARGUMENT],
  ['--sint',          # <num>                 | set interval value for tone plot. this is more dominant than --interval/int.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--levels',        # <val1,val2,val3,...>  | set values of contour/tone levels.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--clevels',       # <val1,val2,val3,...>  | set values of contour levels.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--slevels',       # <val1,val2,val3,...>  | set values of tone levels.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--patterns',      # <pattern1,pattern2,..>| set each patterns for tone plot.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--tone',          # <a|e|f|b|c>           | set tone subroutine:
                      #                       | a (=> tone routine is selected automatically depending on the datasize)
                      #                       | e (=> DCL.uetone is used)
                      #                       | f (=> DCL.uetonf is used)
                      #                       | b (=> DCL.uetonb is used)
                      #                       | c (=> DCL.uetonc is used)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--udsfmt',        # <strings>             | change contour label format. see UDCNTR/DCL manual for the format.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--nocolorbar',    #                       | do not draw color bar.
    GetoptLong::NO_ARGUMENT],
  ['--nozero',        #                       | do not draw zero contour.
    GetoptLong::NO_ARGUMENT],
  ['--nodraw',        #                       | do not draw any figures.
    GetoptLong::NO_ARGUMENT],
  ['--file',       # <png|eps|pdf>[,filename] | raw in file with given format. available formats are png, eps, and pdf.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--crop',       #                          | crop (trim) output png or pdf.
    GetoptLong::NO_ARGUMENT],
  ['--title_array',# <string,string,...>      | set multiple titles given as --title_array hoge,foo,bar; this gives the title "hoge" to
                   #                          | the 1st figure, "foo" to the 2nd one, and "bar" to the 3rd one.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--miscindex',  # <num>                    | set line index (color and width) of title, label, frame, etc.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--index_array',# <index1,index2,...>      | set array of index used for 1st, 2nd,... lines/mark.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--type_array', # <index1,index2,...>      | set array of type index used for 1st, 2nd,... lines/mark.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--subtitle',   # <string>                 | set subtitle, which locates just below the title and is common for multple figures.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--irange',     # <min:max>                | same as --range but with +/- infinity for both boundaries. use with --int is recomended.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--xrange',     # <min:max>                | set the range of x-axis for line plot and 2D plot.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--yrange',     # <min:max>                | set the range of y-axis for line plot and 2D plot.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--clr_range',  # <min:max>                | set color range. maximum range is 10:99. default is 15:94.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--wm',         #                          | alias for --itr 10 --map coast_world --nocont.
    GetoptLong::NO_ARGUMENT],
  ['--size',       # <num>                    | set the size factor of x-window (default: 1.25)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--fullcolor',  # <1|2>                    | draw fullcolor fig. --fullcolor 1 draws with 1 var.
                   #                          | --fullcolor 2 draws with 2 vars; 1st gturl is used for color (hue) and
                   #                          | 2nd one is used for brightness (value).
                   #                          | option of --range xmin:xmax,ymin:ymax can be used.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--step',       #                          | make 1D line plot to step shape.
    GetoptLong::NO_ARGUMENT],
  ['--zerocenter', # [max]                    | set y-axis center in line plot figure to zero.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--vector',     #                          | draw a vector plot.
                   #                          | When 2 gturls are given, 1st one and 2nd one are used for x- and y-components of the vectors.
                   #                          | When 3 or 4 gturls are given, 2nd one and 3rd one are used for x- and y-components of the vectors,
                   #                          | and tone of 1st one and counter of 4th one are overlaid.
                   #                          | If the 1st dim is lon or lat and the 2nd dim is not lon or lat, vectors are scaled by the geometry aspect ratio.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--vfact',      # <factor>                 | use with "--vector" to change the size of the vector by <factor>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--vint',       # <int>                    | use with "--vector" to change the gird-interval to draw vectors.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--vkeep',      #                          | use with "--vector" to keep the size of unit vector as that of the first plot.
    GetoptLong::NO_ARGUMENT],
  ['--nolegend',   #                          | do not draw legend of line plot.
    GetoptLong::NO_ARGUMENT],
  ['--textbox',    # <text>,<num>,[l|c]       | draw textbox, legends, and/or color bar in frame #<num>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--fact',       # <fact>                   | factor for charcter size
    GetoptLong::REQUIRED_ARGUMENT],
  ['--xmean',      #                          | draw line plot of mean along x-axis on the right hand side of the tone/contour plot.
    GetoptLong::NO_ARGUMENT],
  ['--ymean',      #                          | draw line plot of mean along y-axis below the tone/contour plot.
    GetoptLong::NO_ARGUMENT],
  ['--maskshading',# <gturl,maskval,ipat>     | draw a mask for values less than or equal to maskval (optional; defalult is 0)
                   #                          | with given gturl (optional; default is same as previously drawn GPhys)
                   #                          | by ipat pattern number (optional; defalult is 1515).
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--nolabels',   #                          | do not draw labels.
    GetoptLong::NO_ARGUMENT],
  ['--fill',       #                          | fill regions between plotted line and zero (or min or max) line. if histogram, fill the boxes.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--lwf',        # <float>                  | set line with factor for pdf drawing. defalut is 999, which looks like DISP drawing.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--axis_options',# <var1=val1,var2=val2>,..| set some axis options of USPACK (i.e., xoff/yoff, xfac/yfac, dxt/dyt, dxl,dyl)
    GetoptLong::REQUIRED_ARGUMENT],
  ['--bothsides',   #                         | draw bothsides of the planet when itr = 30
    GetoptLong::NO_ARGUMENT],
  ['--xint',  # <xlabelint[:division]>        | set x-axis' label-tickmark interval and number of number of division between two label-tickmarks
              #                               | e.g., --xint 10:5 draws label interval with 10 and sub-tickmark interval of 2 (=10/5).
    GetoptLong::REQUIRED_ARGUMENT],
  ['--yint',  # <ylabelint[:division]>        | same as --xint but for y-axis.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--panelfit', # [margin]                   | force the figure to fit into the panel. margin can be given by ratio [0-1].
    GetoptLong::OPTIONAL_ARGUMENT],
###
### analysis options ###
###
  ['--mean',   # <dim>                        | mean along axis <dim>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--stddev', # <dim>                        | standard deviation along axis <dim>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--eddy',   # <dim>                        | deviation from mean along axis <dim>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--diff',   #                              | calculate 1st gturl - 2nd gturl
    GetoptLong::NO_ARGUMENT],
  ['--divide', #                              | calculate 1st gturl / 2nd gturl
    GetoptLong::NO_ARGUMENT],
  ['--add',    #                              | calculate the sum of multiple variables.
    GetoptLong::NO_ARGUMENT],
  ['--operation',  # <math_func>              | operation of the specified math function on the data.
                   #                          | <math_func> should be a math function with one argument such as log10, sqrt, sin, etc.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--running_mean','--rm', # <dim>:<span>[:<skip>][,<dim>:<span>[:<skip>]]
                            #                 | calculate the running mean along <dim> with a window of <span>.
                            #                 | multiple dims are not supported.
                            #                 | if integer <skip> is given, the date thinning is performed to the running mean.
                            #                 | for example, --rm t,12,12 to monthly data gives annual mean series (one value per year).
    GetoptLong::REQUIRED_ARGUMENT],
  ['--global_mean','--gm',  #                 | calculate the global mean considering the weight along latitude.
                            #                 | dims of longitude and latitude must be the 1st and 2nd dims, respectivly.
    GetoptLong::NO_ARGUMENT],
  ['--ensemble_mean','--em', # [num]          | calculate the ensemble mean.
                             #                | you must provide two or more gturls whose size and dimension are same.
                             #                | if <num> is given, ensemble mean will be calculated for each <num> gturls.
                             #                | if "concat" is given for <num>, provided gturls are concatenated along new dimension "member".
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--normalize',           # [ |mean|float]  | normalize data by initial (when no arg.) or mean or given value when data is one-dimension.
                            #                 | normalization for other data is not implemented yet.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--lowpass','--lp',      # <dim>,<wn>      | apply lowpass filter to <dim>-th dimension with wavenumber <wn> using fft.
                            #                 | i.e., this cut off the wavenumbers greater than <wn>.
                            #                 | <dim> must be integer. dimname is not supported yet.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--highpass','--hp',     # <dim>,<wn>      | apply highpass filter to <dim>-th dimension with wavenumber <wn> using fft.
                            #                 | i.e., this cut off the wavenumbers less than <wn>.
                            #                 | <dim> must be integer 0, 1, 2, or 3. dimnames and dim > 3 are not supported yet.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--power_spectra','--ps',# <dim>[,any]     | calculate power spectra using fft along <dim>-th dimension.
                            #                 | by default, wavenumber is used for axis.
                            #                 | if any argument is given after comma, wavelength or period is used for axis.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--derivative',          # <dim>           | calculate the derivative alog axis <dim>. GPhys.threepoint_O2nd_deriv is used.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--srot',                #                 | calculate rotation on sphere with 1st gturl and 2nd gturl. radius can be set by --radius.
    GetoptLong::NO_ARGUMENT],
  ['--sdiv',                #                 | calculate divergence on sphere with 1st gturl and 2nd gturl. radius can be set by --radius.
    GetoptLong::NO_ARGUMENT],
  ['--sht',                 #                 | use spherical_harmonics_next module for srot, sdiv.
    GetoptLong::NO_ARGUMENT],
  ['--HKE',                 #                 | calculate horizontal Kinetic energy with 1st gturl and 2nd gturl.
    GetoptLong::NO_ARGUMENT],
  ['--integrate', # <dims>                    | calculate integration along <dims>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--sum',       # <dims>                    | calculate sum along <dims>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--addingup',  # <dim>[,<dim>,...]         | adding up values along <dim>; i.e., val[i] =  val[0..i].sum
    GetoptLong::REQUIRED_ARGUMENT],
  ['--max',       # <dims>                    | max values along axis <dims>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--min',       # <dims>                    | min values along axis <dims>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--ES',        # [total|rot|div][,vor_div_given]
                  #                           | calculate Horizontal Kinetic Energy Spectra from U given as 1st gturl
                  #                           | and V given as 2nd gturl. output total, rot, and/or div component.
                  #                           | if "vor_div_given" is given, 1st and 2nd gturls must be Vor and Div, respectivly,
                  #                           | and they will be used for calculation (faster).
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--wnf_analysis',# <sym|asym|bg|sym/bg|asym/bg|log_sym/bg|log_asym/bg>
                    #                         | calculate time-space spectra following Wheeler and Kiladis (1999)
                    #                         | name of component to be output has to be given.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--mvo',  # <mathematical operation>       | mathematical operation for multiple gturls. For mathematical operation,
             #                                | 1st gturl must be indicated by "x", 2nd by "y", 3rd by "z", 4th by "u",
             #                                | 5th by "v", 6th by "w". This cannot treat vars more than 6.
             #                                | Also, this cannot be used repeatedly.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--mvo_only',   #                          | use with "--mvo" to execute "multi variable operation" only.
    GetoptLong::NO_ARGUMENT],





###
### data/constants arrange options ###
###
  ['--cut',# dimname=pos1[:pos2[:thinning_intv]][,dimname=...]: | apply cut, slice, and/or thinning after mathmatical operation.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--radius',  # <raidus>                    | set planetary radius.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--Omega',   # <Omega>                     | set planetary rotation rate [s-1].
    GetoptLong::REQUIRED_ARGUMENT],
  ['--rename',  # <varname>                   | rename the output variable name by <varname>. if --long_name is not used, "long_name" is also renamed.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--long_name', # <long_name>               | rename the output variable "long_name" by <long_name>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--regrid',  # <latgrid_num>               | interpolate to the Gausian latitude and cooresponding longitude.
                #                             | number of lat-grid points is required and available numbers are
                #                             | 8, 16, 32, 64, 94, 96, 120, 128, 160, 240, 256, 320, 480, 640. gglat.rb required [TODO: 自前計算にする]
    GetoptLong::REQUIRED_ARGUMENT],
  ['--GLAT',
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--read_csv',  # <num>                     | read csv file. the first column is used for the axis and the other columns are used for data values.
                  #                           | <num> specifies how many data columns you want to read.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--composite',  # <dim>,<span>,<skip>      | if <skip> is not given or zero, calculate mean with interval <span>.
                   #                          | i.e., --composite time,12 calculates mean of each month from monthly mean data.
                   #                          | if <skip> is given, calculate composite mean for dimension <dim>
                   #                          | with length <span> with interval of <skip>.
                   #                          | i.e., --composite time,3,9 calculates mean for time=0,1,2,12,13,14,24,...
                   #                          | this can be used to calculate seasonal mean such as DJF from mounthly mean data.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--split_axis',# <dim>,<span>,[newaxisname]| split the <dim> axis with <span>.
                  #                           | for example, --split_axis time,12 applied to monthly data splits the time axis
                  #                           | to year and month axis. name of the new axis can be given by [newnaxisname]
    GetoptLong::REQUIRED_ARGUMENT],
  ['--calendar360',  #                        | set the calendar type to "360_day calendar". normally, this will be set automatically.
    GetoptLong::NO_ARGUMENT],
  ['--calendar365',  #                        | set the calendar type to "365_day calendar". normally, this will be set automatically.
    GetoptLong::NO_ARGUMENT],
  ['--flip_data',    # <dim[,dim,..]>         | flip (turn over) the data and axis of <dim>th dimension. <dim> must by given by integer(s).
    GetoptLong::REQUIRED_ARGUMENT],
  ['--interpolate_pole', #                    | extend the latitude axis to have grid points on each pole (-90/+90 deg), and
                         #                    | interpolate the pole values by mean of values on surrounding grid points.
                         #                    | this assumes the data has lon and lat in its first and second dimension, respectivly.
    GetoptLong::NO_ARGUMENT],
  ['--set_missing_value',# [float]            | replace the missing values to given value and validate then. default: 0.0.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--ignoreunit',       #                    | ignore unit
    GetoptLong::NO_ARGUMENT],
  ['--set_assoc_coords', # <gturl>            | set <gturl> as an associate coordinate.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--sig2p', # <gturl of surface pressure>   | transform sigma coordinate to pressure coordinate. gturl of surface pressure must be given.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--sig2z', # <gturl of Temperature>        | transform sigma coordinate to z (not log_p) coordinate. gturl of Temperature must be given.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--interpolate', # <args>                  |   interpolate after operations.
                    #                         | <args> can be axis_name=90 or axis_name=[0,90,180,270] or axis_name=0:360:4 for example.
                    #                         | 1) axis_name=90 style gives interpolation at the location and rank shrinks.
                    #                         | 2) axis_name=[0,90,180,270] style gives interpolation for multiple location and use then as an axis.
                    #                         | 3) axis_name=0:360:4 style similar to (2) but 0:360:4 means from 0 to 360 with separation of (360-0)/4.
                    #                         | multiple interpolation is not supported yet.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--extrapolation',   #                     | enable extrapolation, when using --interpolate option.
    GetoptLong::NO_ARGUMENT],
  ['--unit',        # <unit>                  | set the operated variable's unit.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--invalidation',# <num or range>          | invalidate the given value or range.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--join',        #                         | join muliple gturls along with an axis whose values are continual across the gturls and are not duplicate.
    GetoptLong::NO_ARGUMENT],
  ['--extend_backward_along', # <dim>,<len>[,<value>] | extend backward the object along <dim>-axis for <len> numbers with <value>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--thinning', # <num>                      | data thinning by picking up data with <num> interval
    GetoptLong::REQUIRED_ARGUMENT],
  ['--nothinning', #                          | prevent auto data thinning.
    GetoptLong::NO_ARGUMENT],
  ['--read_fits',  #                          | read fits format image. automatically added if file extension is .fits. (RubyFits required)
    GetoptLong::NO_ARGUMENT],
  ['--offset',     # <float>                  | offset the data by <float>.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--cyclic_extension',# <dim:num[,dim:num]> | extend the data along <dim> with <num> number cyclically. e.g., --cyclic_extension x:4,y:4.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--axis_units_name', # <axdim:unit[:axname:operation]>[,axdim:unit:axname:operation]
                        #                     | change unit (and name) of axis. multiple axes can be given as
                        #                     | --axis_units_name ax0:deg:lon,ax2:km:height:/1000  for example.
                        #                     | mathematical operation can be applied to the axis values.
    GetoptLong::REQUIRED_ARGUMENT],
  ['--replace_axis',#<axN:val1,val2,val3,...> | replace the value of the axis <axN> (N is 0, 1, 2,...) with given Array [val1, val2, val3, val4,...]
    GetoptLong::REQUIRED_ARGUMENT],
  ['--zonal_shift', #<speed at equator in m/s>| shift entire field zonally as it is advected by the solid body rotation of <speed at equator in m/s>.
    GetoptLong::REQUIRED_ARGUMENT],







###
### output options ###
###
  ['--nc',    # [outfilename]                 | output operated GPhys object as NetCDF file. output file name can be specified.
              #                               | Default name is "out.nc".
              #                               | if used with "--anim", GPhys object for last shown figure is outputtedd.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--nc4',   # [outfilename]                 | output operated GPhys object as NetCDF ver 4 with compression.
              #                               | may take longer time to write out.
              #                               | output file name can be specified. Default name is "out.nc".
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--nc4a',  # [outfilename]                 | output operated gphys continually with along time dim.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--csv',   # [outfilename]                 | output operated GPhys object (1 dim only) as CSV file. output file name can be specified.
              #                               | Default name is "out.csv".
              #                               | if used with "--anim", GPhys object for last shown figure is outputtedd.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--produce', # <something>                 | produce <something>. <something> can be "solar_zenith_angle"
    GetoptLong::REQUIRED_ARGUMENT],
  ['--output_assoc_coords',  #                | use with --nc4 or --nc to output associate coordinate as netcdf file.
    GetoptLong::NO_ARGUMENT],
  ['--gif',   #                               | output gif animation. this must be used with --anim and --file png options.
    GetoptLong::NO_ARGUMENT],
  ['--mp4',   # [number of frames per second] | output mp4 animation. frame rate can be given. this must be used with --anim and --file png options.
    GetoptLong::OPTIONAL_ARGUMENT],
  ['--ntype', # <type>                        | convert type of values. <type> should be "sfloat", "float", etc.
    GetoptLong::REQUIRED_ARGUMENT],

  #                   ['--dlon',                     GetoptLong::NO_ARGUMENT],
  #                   ['--dlat',                     GetoptLong::NO_ARGUMENT],
  #                   ['--version',                  GetoptLong::NO_ARGUMENT]  # to be defined

                     )
  @USED_OPTIONS = []
  begin
    parser.each_option do |name, arg|
      eval "@OPT_#{name.sub(/^--/, '').gsub(/-/, '_')} = '#{arg}'"  # strage option value to @OPT_val
      @USED_OPTIONS << "OPT_#{name.sub(/^--/, '').gsub(/-/, '_')}"
    end
  rescue
    help
    raise
  end
  ARGV.concat(prev_argv) if prev_argv
end








# gp: input GPhys or GPhps Array, argv: string of gpv optiones,
# window: flag for DCL window (0: open & close, 1: close only,
# 2: open only, 3: don't open or close)
def main(gp = nil, argv=nil, window=0)
  if (gp or argv.class == String) then
    until (ARGV.length == 0) do
      ARGV.shift
    end
    ARGV.concat(argv.split(" ")) if argv
  end
  @flag_window_open = window


#####################################################
###++++++           Main Routine            ++++++###
## store command line and current directory
@commandline = $0.split("/")[-1] + " " + ARGV.join(" ")
@current_dir = Dir.pwd.to_s + "/"
@comment = @commandline + " (in #{@current_dir} @ #{Time.now})"

# set options
__set_options

# store source file names to @sources
@sources = []
ARGV.each{|s| @sources << s.split("@")[0] if File.exist?(s.split("@")[0]) }

## Print out help message
if (@OPT_help)
  help
  exit(1)
end

## exec gplist
if (@OPT_list) then
  ARGV.each{|v|
    v = v.split("@")[0];  system ("gplist #{v}")
    begin
      gtime = GPhys::IO.open(v,"time")
      cal = gtime.get_att("calendar")
      sdate = UNumeric.new(gtime.val[0],gtime.units.to_s).to_datetime(0.0,cal)
      edate = UNumeric.new(gtime.val[-1],gtime.units.to_s).to_datetime(0.0,cal)
      print "   time is from #{sdate} \n"
      print "             to #{edate} in calendar-type '#{cal}'\n"
    rescue
    end
    vars_and_axes = GPhys::IO.var_names(v); vars = GPhys::IO.var_names_except_coordinates(v)
    vars_and_axes.each{|a|
      unless (vars.include?(a)) then
        axis = GPhys::IO.open(v,a).val
        direction = axis[0]<axis[-1] ? "POSITIVE" : "NEGATIVE"
        print "   #{a}-axis is packed from #{axis[0]} to #{axis[-1]} (#{direction} direction). \n"
      end

    }

  }
  exit
end



## alias options
if (@OPT_wm)
  @OPT_itr = 10 unless @OPT_itr
  @OPT_map = "coast_world"  unless @OPT_map
  @OPT_nocont = true unless (@OPT_cint or @OPT_crange)
end

@OPT_title = "" if @OPT_notitle



## set radius (defalut is 6370 km)
if (@OPT_radius)
  GAnalysis::Planet.radius=UNumeric.new(@OPT_radius.to_f, "m")
  @Radius = @OPT_radius.to_f
  print "Planetary radius is set to #{@Radius} m.\n"
else
  @Radius = GAnalysis::Planet.radius
end

## set planetary rotation rate (defalut is 7.27220521664304e-05 s-1)
if (@OPT_Omega)
  GAnalysis::Planet.omega=UNumeric.new(@OPT_Omega.to_f, "s-1")
  @Omega = @OPT_Omega.to_f
  print "Omega, plaetary rotation rate, is set to #{@Omega} s-1.\n"
else
  @Omega = GAnalysis::Planet.omega
end

## require spherical harmonics library if needed
if (@OPT_sph)
  require "spherical_harmonics"
  include SphericalHarmonics
end


## set some figure option
DCL::swlset('lwait', false) if (@OPT_nowait    || @OPT_Gw || @OPT_smooth || @OPT_Gaw)
                                           # set wait or nowait
DCL::swlset('lalt',  true)  if (@OPT_alternate || @OPT_Ga || @OPT_smooth || @OPT_Gaw)
                                           # set backing store option
if (@OPT_noannotate || @OPT_panelfit)
  @annotate = false
else
  @annotate = true
end

@Overplot_max = ( @OPT_overplot.to_i || 1 )
@Overplot = 1

if (@OPT_nocolorbar)
  @colorbar = false
else
  @colorbar = true
end

if (@OPT_tone)
  case @OPT_tone
  when "a"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = true,false,false,false,false
  when "e"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = false,false,false,false,false
  when "f"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = false,true,false,false,false
  when "b"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = false,false,true,false,false
  when "c"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = false,false,false,true,false
  when "full"
    @auto, @tonf, @tonb, @tonc, @tone_fullcolor = false,false,false,false,true
  else
    raise "The value of option --tone should be 'a','e','f','b','c', or 'full'."
  end
else
  @auto, @tonf, @tonb, @tonc, @tone_fullcolor = true,false,false,false,false
end

if (@OPT_log_int)
  @log_int = true
else
  @log_int = false
end


## decide VIEWPORT
##@VIEWPORT = set_vpsize( VIEWPORT, (@OPT_aspect||2.0) )

## tune the size of axis parameters.
#DCL.uzfact(0.7)

## set the format of contour labels.
udsfmt = (@OPT_udsfmt||DCL.udqfmt())
DCL.udsfmt(udsfmt)

## draw figure
loopdim   = ( @OPT_animate || @OPT_anim )
loopdim = loopdim.to_i if loopdim.to_i.to_s == loopdim
loopdim = loopdim[-1].to_i if loopdim =~ /#{AX}[0-9]/

### set colormap
#binding.pry
if (@OPT_clrmap) then
  if ( File.exist?(@OPT_clrmap) ) then
    DCL.swpset("CLRMAP", @OPT_clrmap)
  else
    DCL.sgscmn(@OPT_clrmap||63)
  end
else
  DCL.sgscmn(63)
end

### set title_array
@title_array = @OPT_title_array.split(",") if @OPT_title_array

### set index_array
if (@OPT_index_array) then
  @index_array = @OPT_index_array.split(",")
  @index_array.map!{|s|
    n, t = s.split("x")
    if t then Array.new(t.to_i).fill(n.to_i)
    else n.to_i end
    }
  @index_array.flatten!
elsif (@OPT_index.to_i >= 10) then # if index >= 10 is set, use same index for all lines.
  @OPT_index = @OPT_index.to_i
elsif (@OPT_nocolor)
  @OPT_index = 2
else #default
  @index_array = [10,20,40,400,650,300,930,12,22,42,402,652,302,932]
  @index_array.map!{|i| i + (@OPT_index.to_i||2)}
end

### set type array
if (@OPT_type_array) then
  @type_array = @OPT_type_array.split(",")
  @type_array.map!{|s|
    n, t = s.split("x")
    if t then Array.new(t.to_i).fill(n.to_i)
    else n.to_i end
    }
  @type_array.flatten!

else
  @type_array = [1,3,4,2,5]
end
### set alphabet array
@alphabet = ["a"]; 25.times{|i| @alphabet << @alphabet[-1].next}


#@flag_window_open = false
@flag_margin = true if (@flag_window_open == 0 or @flag_window_open == 2)

#settings for text box

if (@OPT_textbox.class == String) then
  @OPT_textbox.split(",").each{|s|
    if (s == s.to_i.to_s) then
      @textbox_fnum = s.to_i
    elsif (s == "l" or s == "c") then
      @textbox_lc = true
    else
      @textbox = s
    end
  }
  @textbox = "" unless @textbox
  @textbox_fnum = (ARGV.size)/($Overplot_max) unless @textbox_fnum
end



gary = Array.new

@flag_read_csv = 1


## set Land/Ocean mask
if (@OPT_land or @OPT_ocean or @OPT_color_scatter)
  landsea = GPhys::IO.open("/Users/hiroki/Dropbox/RubyScripts/landsea.nc","LSMASK").copy #海陸マスクを読み込む 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
  lsgrid = landsea.grid.copy

  land_mask = NArrayMiss.to_nam(landsea.val.ne(0)) # 0=Ocean, 1=Ocean以外にする
  land_mask.invalidation(landsea.val.eq(0))
  lary = VArray.new(land_mask,{"long_name"=>"land mask"} ,"lmask")
  lmask = GPhys.new(lsgrid, lary)

  ocean_mask = NArrayMiss.to_nam(landsea.val.eq(0)) # 0=Ocean, 1=Ocean以外にする
  ocean_mask.invalidation(landsea.val.ne(0))
  oary = VArray.new(ocean_mask,{"long_name"=>"ocean mask"} ,"omask")
  omask = GPhys.new(lsgrid, oary)
end

# set ARGV if gp is given
if (gp.class == GPhys) then
  ARGV.unshift(gp)
elsif (gp.class == Array) then
  gp.reverse_each{|g| ARGV.unshift(g) }
end


# preparation for set_assoc_coords
if (@OPT_set_assoc_coords) then
  @assoc_coords =  open_gturl(@OPT_set_assoc_coords)
end

# preparation for sig2p
if (@OPT_sig2p) then
  @PS = open_gturl(@OPT_sig2p)
  @PS = @PS[true,true,0,false] if @PS.shape[2] == 1 # 長さ１の鉛直軸を持ってる場合は除去する
end

# preparation for sig2z
if (@OPT_sig2z) then
  @T = open_gturl(@OPT_sig2z) # 温度 T が必要
end

# enable extrapolation or not
if (@OPT_extrapolation) then
  GPhys.extrapolation=true
end



#
#---------------------------------------------------------------------
#      MAIN LOOP
#---------------------------------------------------------------------
## open netcdf variables
while ARGV[0] do

  gturl = ARGV[0]

  if (@OPT_read_fits ) then
    gp = fits2gphys(gturl)
    ARGV.shift
  else
    gp, gturl = open_gturl_wildcard(gturl)
    next unless gp # gp == false なら、つぎのgturlを検索
    ARGV.shift
  end

  ## check whether gp has time axis and its calendar type
  time_axis = gp.axnames.find{|i| ["time","TIME","Time","t","T"].include?(i)}
  if (time_axis) then
    calendar = gp.axis(time_axis).to_gphys.get_att("calendar")
    case calendar
    when "360_day", "uniform30day" then
      @OPT_calendar360 = true
    when "all_leap", "366_day" then #not implemented yet.
    when "noleap", "365_day" then
      @OPT_calendar365 = true
    when "gregorian", "standard", "", nil then
    end
  end

  # replace axis values
  if (@OPT_replace_axis) then
    axnum, newaxis = @OPT_replace_axis.split(":")
    num = axnum.sub("ax","").to_i
    axary = newaxis.split(",")
    axary.map!{|a| a.to_f}
    unit = gp.coordinate(num).units.to_s
    axname = gp.axis(num).name
    gp.axis(num).set_pos(VArray.new(NArray.to_na(axary), {"units"=>unit}, axname))
    print "#{axname}'s value was replaced by #{axary}"
  end


  # force axis unit
  if (@OPT_axis_units_name) then
    @OPT_axis_units_name.split(",").each{|it|
      axnum, unit, axname, ope = it.split(":")
      num = axnum.sub("ax","").to_i
      unit = gp.coordinate(num).units.to_s if unit == "" # 単位未指定時は変更しない
      axname = gp.axis(num).name unless axname # 軸名未指定時は変更しない
      print "axis #{gp.axis(num).name} [#{gp.coord(num).units}] was changed to "
      typecode = gp.coordinate(num).typecode
      axis_val = gp.coordinate(num).val.to_type(typecode)
      axis_val = (eval "axis_val #{ope}") if ope
      gp.axis(num).set_pos(VArray.new(axis_val, {"units"=>unit}, axname))
#      gp.axis(num).set_pos(VArray.new(gp.coordinate(num).val, {"units"=>unit}, axname))
      print "#{axname} [#{unit}] \n"
    }
  end

  # cyclic_extension
  if (@OPT_cyclic_extension) then
    ary = @OPT_cyclic_extension.split(",")
    ary.each{|a|
      dim, num = a.split(":")
      num = 1 if num == nil
      gp = cyclic_extension(gp,dim,num.to_i).copy
    }
  end


  # check whether gp's size and if large turn on the regrid option.
  unless (@OPT_regrid || @OPT_nothinning) then
    thres = 180
    if (@OPT_thinning) then
      gp = data_thinning(gp, @OPT_thinning.to_i)
    elsif (gp.shape.length >= 2) then
      @OPT_nocont = true if (gp.first2D.shape.max > thres*2 && !@OPT_noshade)
      #gp = data_thinning(gp, thres)
    end
  end


  ## flip data and axis of <dim> th dimension
  ## standard deviation along any axis (preparation)
  if (@OPT_flip_data) then
    dims_flip = set_dims(@OPT_flip_data)
    dims_flip.each{|dim|
      gp = reverse(gp,dim)
    }
  end


  ## extend gphys object along dim-axis
  if (@OPT_extend_backward_along) then
    dim, len, val = @OPT_extend_backward_along.split(",")
    gp = extend_backward_along(gp, dim, len.to_i, val.to_f)
  end


  # interpolate the value on the poles (i.e., lat = -90/+90 deg by the mean of values on surrounding grid points)
  gp = interpolate_pole(gp) if @OPT_interpolate_pole

  gp.set_att("units","") if @OPT_ignoreunit


  # preparation to use spherical_harmonics_next module
  if (@OPT_sht or @OPT_ES) then
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, true)
    nmax = gp.coordinate(lon_dim).length/3
    gw, pmn, dpmn = SphericalHarmonics::sh_init(gp.coordinate(lon_dim).val, gp.coordinate(lat_dim).val, nmax, gp)
  end

  if (@OPT_zonal_shift) then
    gp = zonal_shift_on_sphere(gp, @OPT_zonal_shift)
  end

  if (@OPT_wnf_analysis) then
    gp = wnf_analysis(gp)
  end



  ## for case of calculating difference of two gphys object
  if (@OPT_diff or @OPT_srot or @OPT_sdiv or @OPT_divide or @OPT_HKE or @OPT_ES)
    raise "--diff/rot/div option must be used with even numbers (2, 4, 6,...) of gturls" if ARGV[0] == nil
    prev_gturl = gturl; prev_gp = gp
    gturl = ARGV[0]
    gp, gturl = open_gturl_wildcard(gturl)
    unless (gp)
      gturl = ARGV[0]
      gp, gturl = open_gturl_wildcard(gturl)
    end

    gp.set_att("units","") if @OPT_ignoreunit


    if (@OPT_diff)
      print "  Taking difference: #{prev_gturl} - #{gturl}\n" unless @OPT_silent
      gturl = prev_gturl+" - "+gturl
      gp = prev_gp - gp
    elsif (@OPT_divide)
      print "  Taking division: #{prev_gturl} / #{gturl}\n" unless @OPT_silent
      gturl = prev_gturl+" / "+gturl
      gp = prev_gp / gp
    elsif (@OPT_srot)
      print "  calculating rot(#{prev_gturl},#{gturl}) on sphere with radius #{@Radius}\n" unless @OPT_silent
      gturl = "rot("+prev_gturl+","+gturl+") on sphere with radius #{@Radius}"
      gp = GAnalysis::Planet.rot_s(prev_gp,gp) unless @OPT_sht
      gp = SphericalHarmonics::sh_vorticity(prev_gp,gp,@Radius) if @OPT_sht
    elsif (@OPT_sdiv)
      print "  calculating div(#{prev_gturl},#{gturl}) on sphere with radius #{@Radius}\n" unless @OPT_silent
      gturl = "div("+prev_gturl+","+gturl+") on sphere with radius #{@Radius}"
      gp = GAnalysis::Planet.div_s(prev_gp,gp) unless @OPT_sht
      gp = SphericalHarmonics::sh_divergence(prev_gp,gp,@Radius) if @OPT_sht
    elsif (@OPT_HKE)
      print "  Calculating Horizontal Kinetic Energy: [(#{prev_gturl})^2 + (#{gturl})^2]/2\n" unless @OPT_silent
      gturl = "[("+prev_gturl+")^2 + ("+gturl+")^2]/2"
      gp = (prev_gp*prev_gp + gp*gp)*0.5
      gp.rename("HKE")
      gp.set_att("long_name","Horizontal Kinetic Energy")

    elsif (@OPT_ES)
      comp_ary = @OPT_ES.split(",")
      if (comp_ary.include?("vor_div_given")) then
        flag_vor_div_given = true
        comp_ary.delete("vor_div_given")
      else
        flag_vor_div_given = false
      end
      comp_ary << "total" if comp_ary.size == 0


      print "  Calculating Energy Spectra \n" unless @OPT_silent
      gturl = "Energy spectra from " +prev_gturl + " and " + gturl
      if (gp.rank > 2) then # use parallel
        prev_gp = prev_gp.copy; gp = gp.copy
        np = (ENV['OMP_NUM_THREADS']|| 4).to_i
        dname = gp.axnames[-1]
        print "Operation is processed along #{dname} dim with #{np} proccesses. This may take a while...\n"
        gpa = Parallel.map(gp.coordinate(-1).val.to_a, :in_processes=>np, :progress=>"progress"){|i|
    # #        gp, gturl = open_gturl_wildcard(gturl, true)
    #       tmp_gp = reopen(gp).cut_rank_conserving(dims_remained[-1]=>i)
    #       if @OPT_set_assoc_coords then
    #         tmp_assoc_coords = reopen(@assoc_coords).cut_rank_conserving(dims_remained[-1]=>i)
    #         tmp_gp.set_assoc_coords(tmp_assoc_coords)
    #       end
          SphericalHarmonics::sh_energyspectrum(prev_gp.cut_rank_conserving(dname=>i), gp.cut_rank_conserving(dname=>i) , @Radius, comp_ary, flag_vor_div_given)
        }
        gp = gpa.transpose.map{|ga| GPhys.join(ga) }
      else
        gp = SphericalHarmonics::sh_energyspectrum(prev_gp, gp, @Radius, comp_ary, flag_vor_div_given)
      end
      if (gp.size == 1) then
        gp = gp[0]
      else
        ARGV << gp[1..-1]; ARGV.flatten!
        gp = gp[0]
      end
      @OPT_ES = nil
    end
    ARGV.shift

  end


#  gp = GAnalysis::Planet.grad_sx(gp) if @OPT_dlon
#  gp = GAnalysis::Planet.grad_sy_cosphifact(gp,1) if @OPT_dlat

  ## for case of calculating sum of gphys objects
  if (@OPT_add or @OPT_ensemble_mean)
    raise "--add/em option must be used with gturls more than one" if ARGV[0] == nil
    print "  Summing up: #{gturl} \n" unless @OPT_silent
    count = 1;
    gp = regrid(gp) if @OPT_regrid
    gpary = [gp] if @OPT_ensemble_mean == "concat"
    while (ARGV[0]) do
      prev_gturl = gturl; prev_gp = gp
      gturl = ARGV[0]
      gp, gturl = open_gturl_wildcard(gturl)
      gp = regrid(gp) if @OPT_regrid
      gpary << gp if @OPT_ensemble_mean == "concat"
      gp = prev_gp + gp
      print "  adding up:   + #{gturl}\n" unless @OPT_silent
      gturl = prev_gturl+" + "+gturl
      ARGV.shift
      count += 1
      break if count == @OPT_ensemble_mean.to_i
    end
    if (@OPT_ensemble_mean == "concat") then
      gp = GPhys.concat(gpary, NArray.int(count).indgen!, "member", {"unit"=>"1", "long_name"=>"ensemble member"}).copy
      print "  Ensemble members are concatenated.\n" unless @OPT_silent
    elsif @OPT_ensemble_mean then
      gp.assoc_coords=nil if (gp.assoc_coords) # assoc_coordsがあると演算がうまくいかない
      gp = gp/(count.to_f)
      print "  Taking the ensemble mean.\n" unless @OPT_silent
    end
  end

  ## for case of appling mathematical operation to multiple variables
  if (@OPT_mvo) then
    if (@OPT_ntype == "sfloat") then typecode = 4 else typecode = nil end
    unused_gturl = []
    if @OPT_mvo.include?("ax0")
      ax0 = gp.axis(0).to_gphys.val
      ax0 = ax0.to_type(@OPT_ntype) if @OPT_ntype
    end
    if @OPT_mvo.include?("ax1") then
      ax1_tmp = gp.axis(1).to_gphys.val; ax1 = NArray.new((typecode || ax1_tmp.typecode), *gp.first2D.shape)
      ax1.shape[0].times{|i| ax1[i,true] = ax1_tmp}
      ax1 = VArray.new(ax1, {"units"=>gp.coord(1).units.to_s}, gp.axis(1).name)
    end
    if @OPT_mvo.include?("ax2") then
      ax2_tmp = gp.axis(2).to_gphys.val; ax2 = NArray.new((typecode || ax2_tmp.typecode), *gp.shape[0..2])
      ax2.shape[0].times{|i|
        ax2.shape[1].times{|j| ax2[i,j,true] = ax2_tmp}
      }
      ax2 = VArray.new(ax2, {"units"=>gp.coord(2).units.to_s}, gp.axis(2).name)
    end
    if @OPT_mvo.include?("ax3") then
      ax3_tmp = gp.axis(3).to_gphys.val; ax3 = NArray.new((typecode || ax3_tmp.typecode), *gp.shape[0..3])
      ax3.shape[0].times{|i|
        ax3.shape[1].times{|j|
          ax3.shape[2].times{|k| ax3[i,j,k,true] = ax3_tmp}
        }
      }
      ax3 = VArray.new(ax3, {"units"=>gp.coord(3).units.to_s}, gp.axis(3).name)
    end

    case gp.ntype
    when "sfloat" then  nb = 4
    when "float" then nb = 8
    else nb = 4
    end

    if ( (gp.length*nb > 2.2E9 || @OPT_parallel) && NArrayType == "standard" ) then # for GPhys ogject larger than 2GB
      np = (ENV['OMP_NUM_THREADS']|| 4).to_i
      da_name = gp.coordinate(-1).name
      if (@OPT_parallel && @OPT_parallel != "") then
        da_name = @OPT_parallel # 分割処理する軸を陽に指定する
      end
      da_val_ary = gp.coordinate(da_name).val.to_a
      print "Operation is processed along #{da_name} dim with #{np} proccesses. This may take a while...\n"
      xx = gp ; print " x : #{gturl}\n"
      unless (ARGV.empty?) then
        yy, zz, uu, vv, ww = ["y","z","u","v","w"].map{|i|
          if (@OPT_mvo.include?(i) && !ARGV.empty?) then
            gp, gturl = open_gturl_wildcard(ARGV[0])
            gp, gturl = open_gturl_wildcard(ARGV[0]) if gp == false
            print " #{i} : #{gturl}\n"
            ARGV.shift
            gp
          else
            unused_gturl << ARGV[0]
            ARGV.shift; nil
          end
        }
      end

      if (xx.data.file.class == NetCDF) then # NetCDF4 対策
        if (xx.data.file.format >= 3 ) then
          xx, yy, zz, uu, vv, ww = [xx, yy, zz, uu, vv, ww].map{|tt|
            if (tt) then
              tt = (0..(da_val_ary.length-1)).to_a.map{|i|
                if tt.axnames.include?(da_name) then
                  tt.cut_rank_conserving(da_name=>da_val_ary[i]).copy
                else
                  tt.copy unless tt.axnames.include?(da_name)
                end
              }
            end
          }
          @flag_on_memory = true
        end
      end

      print "operation: #{@OPT_mvo} \n"
      gpa = Parallel.map((0..(da_val_ary.length-1)).to_a, :in_processes=>np, :progress=>"progress"){|i|
        if (@flag_on_memory) then
          x = xx[i]
        else
          x = reopen(xx).cut_rank_conserving(da_name=>i)
        end
        y, z, u, v, w = [yy, zz, uu, vv, ww].map{|tt|
          if (tt) then
            if (@flag_on_memory) then
              t = tt[i]
            else
              t = reopen(tt).cut_rank_conserving(da_name=>da_val_ary[i]) if tt.axnames.include?(da_name)
              t = tt unless tt.axnames.include?(da_name)
            end
          else
            t = nil
          end
          t
        }
        eval "x = #{@OPT_mvo}"; x
      }
      x = xx; y = yy; z = zz; u = uu; v = vv; w = ww
      x, y, z, u, v, w = [x,y,z,u,v,w].map{|t|
        if (t.class == Array) then
          t[0]
        else
          t
        end
      }
#      @flag_on_memory = nil; @OPT_parallel = false
      if (@OPT_nc4a) then
        @flag_mvo_gpa = true; @OPT_parallel = true
      else
        print "combining #{gpa.size} gphys objects... this may take very long time."
        gp = GPhys.join(gpa)
      end


    else #single proccess
      x = gp ; print " x : #{gturl}\n"
      unless (ARGV.empty?) then
        y, z, u, v, w = ["y","z","u","v","w"].map{|i|
          if (@OPT_mvo.include?(i) && !ARGV.empty?) then
            if (ARGV[0].class == String) then
              gp, gturl = open_gturl_wildcard(ARGV[0])
              gp, gturl = open_gturl_wildcard(ARGV[0]) if gp == false
            elsif (ARGV[0].class == GPhys) then
              gp = ARGV[0]; gturl = gp.name
            end
            print " #{i} : #{gturl}\n"
            ARGV.shift
            gp
          else
            unused_gturl << ARGV[0]
            ARGV.shift; nil
          end
        }
      end
#      z = z.copy; z.rename("theta")
      print "operation: #{@OPT_mvo} \n"
#      binding.pry
      eval "gp = #{@OPT_mvo}"
    end
    varname = @OPT_mvo.clone
    varname.gsub!(/exp/, "EXP")
    varname.gsub!(/ax0/, x.axnames[0]) if ax0
    varname.gsub!(/ax1/, x.axnames[1]) if ax1
    varname.gsub!(/ax2/, x.axnames[2]) if ax2
    varname.gsub!(/ax3/, x.axnames[3]) if ax3

    varname.gsub!(/x/, x.name) if x
    varname.gsub!(/y/, y.name) if y
    varname.gsub!(/z/, z.name) if z
    varname.gsub!(/u/, u.name) if u
    varname.gsub!(/v/, v.name) if v
    varname.gsub!(/w/, w.name) if w
    unless (gp.class == GPhys) then
      p gp; abort
    else
      begin
        gp.rename(varname)
        gp.set_att("long_name", varname)
      rescue

      end
      gturl = "operation: #{varname}"
    end

    # put back the unused gturl to ARGV and set @OPT_mvo to false
    unused_gturl.each{|g|  ARGV << g if g }
    @OPT_mvo = false
  end

  ## for case of calculating sum of gphys objects
  if (@OPT_join)
    raise "--add/em option must be used with gturls more than one" if ARGV[0] == nil
    print "  Joining: #{gturl} \n" unless @OPT_silent
    count = 1
    while (ARGV[0]) do
      prev_gturl = gturl; prev_gp = gp
      gturl = ARGV[0]
      gp, gturl = open_gturl_wildcard(gturl)
      gp = GPhys.join([prev_gp, gp])
      print "  Joining:   + #{gturl}\n" unless @OPT_silent
      gturl = prev_gturl+" + "+gturl
      ARGV.shift
      count += 1
    end
    # 軸の値に重複がないか調べる。あればエラーを返す。
    gp.axnames.each{|a|
      axis = gp.coord(a).val.to_a
      if (axis.length != axis.uniq.length) then
        raise "ERROR: joined GPhys object has dupulicate value(s) in axis '#{a}''"
      end
    }
  end




  ###---------------------------------------------

  ## apply composite
  gp = composite(gp) if (@OPT_composite)

  ## apply split axis
  gp = split_axis(gp) if (@OPT_split_axis)

  ## apply lowpass filter
  gp = lowpass(gp) if (@OPT_lowpass)

  ## apply highpass filter
  gp = highpass(gp) if (@OPT_highpass)

  ## interpolate to gausian latitude and correspond longitude
  ## --regrid requires number of latitude grid poins.
  ## --GLAT can automatically set the number of grid points.
  gp = gglat(gp) if (@OPT_GLAT)
  gp = regrid(gp) if (@OPT_regrid)


  ## apply Land/Ocean mask
  if (@OPT_land or @OPT_ocean) # or @OPT_color_scatter)
    gp_lon = gp.axis(0).to_gphys.val; gp_lat = gp.axis(1).to_gphys.val
    gp.assoc_coords=nil if (gp.assoc_coords) # assoc_coordsがあると演算がうまくいかない
    if (@OPT_land)
      gp_lm = lmask.cut(gp_lon,gp_lat)
      gp = gp*gp_lm
    elsif (@OPT_ocean)
      gp_om = omask.cut(gp_lon,gp_lat)
      gp = gp*gp_om
    end
  end

  ## power spectra
  # gp = powerspectra(gp) if (@OPT_power_spectra)

  ## derivative along any axis (preparation)
  dims_deriv = set_dims(@OPT_derivative) if @OPT_derivative

  ## mean along any axis (preparation)
  dims_mean = set_dims(@OPT_mean) if (@OPT_mean)

  ## running mean along any axis (preparation)
  if (@OPT_running_mean)
    rm_ary = @OPT_running_mean.split(",")
    rm_ary.map!{|a|
      rm_dim, rm_span, rm_skip = a.split(/\s*:\s*/)
      if rm_dim.to_i.to_s  == rm_dim then
        rm_dim = rm_dim.to_i
      end
      rm_span = rm_span.to_i
      rm_skip = rm_skip.to_i
      [rm_dim, rm_span, rm_skip]
    }
  end

  ## standard deviation along any axis (preparation)
  dims_stddev = set_dims(@OPT_stddev) if (@OPT_stddev)

  ## max along any axis (preparation)
  dims_max = set_dims(@OPT_max) if (@OPT_max)

  ## min along any axis (preparation)
  dims_min = set_dims(@OPT_min) if (@OPT_min)

  ## deviation from mean along any axis (preparation)
  dims_eddy = set_dims(@OPT_eddy) if (@OPT_eddy)

  ## integration along any axis (preparation)
  dims_integrate = set_dims(@OPT_integrate) if (@OPT_integrate)

  ## sum along any axis (preparation)
  dims_sum = set_dims(@OPT_sum) if (@OPT_sum)

  ## adding up along any axis (preparation)
  dims_addingup = set_dims(@OPT_addingup) if (@OPT_addingup)

  GGraph.margin_info($0, gturl, nil, nil, 0.0, 0.0, nil, 0.0) if @annotate && @flag_margin # draw margin infomation
  DCL.slsttl(gturl, 'b',  1.0, -1.0, 0.008, 2) if @annotate
  @flag_margin = false
  @kind_of_fig = nil

  # set missing value to specific value.
  if (@OPT_set_missing_value) then
    gp = gp.copy
    gp.replace_val(gp.val.set_missing_value(@OPT_set_missing_value.to_f).all_valid)
  end


  ope_string = ""
  if (@OPT_power_spectra) then
    ope_string = ope_string << "->power_spectra(#{@OPT_power_spectra})"
  end

  @USED_OPTIONS.each{|opt|
    case opt
    when "OPT_derivative" then
      dims_deriv.each{|dim| ope_string = ope_string << "->derivative(#{dim})" }
    when "OPT_running_mean" then
      rm_ary.each{|rm| ope_string = ope_string << "->running_mean(#{rm})" }
    when "OPT_mean" then
      dims_mean.each{|dim| ope_string = ope_string << "->mean(#{dim})" }
    when "OPT_global_mean" then
      ope_string = ope_string << "->global_mean"
    when "OPT_stddev" then
      dims_stddev.each{|dim| ope_string = ope_string << "->stddev(#{dim})" }
    when "OPT_max" then
      dims_max.each{|dim| ope_string = ope_string << "->max(#{dim})" }
    when "OPT_min" then
      dims_min.each{|dim| ope_string = ope_string << "->min(#{dim})" }
    when "OPT_eddy" then
      dims_eddy.each{|dim| ope_string = ope_string << "->eddy(#{dim})"}
    when "OPT_integrate" then
      dims_integrate.each{|dim| ope_string = ope_string << "->integrate(#{dim})"}
    when "OPT_sum" then
      dims_sum.each{|dim| ope_string = ope_string << "->sum(#{dim})"}
    when "OPT_addingup" then
      dims_addingup.each{|dim| ope_string = ope_string << "->addingup(#{dim})"}
    else
    end
  }
  print "Operations:#{ope_string} in this order.\n" if (ope_string != "" && !@OPT_silent )

#------------------------------------------------------------------------
  proc = Proc.new do |g,i|

      ## power spectra
      g = powerspectra(g) if (@OPT_power_spectra)


      ## simple operation with given order
      @USED_OPTIONS.each{|opt|

        ## derivative
        if (opt == "OPT_derivative")
        dims_deriv.each{|dim|
          begin
            lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(g, true)
          rescue
          end
          dim = g.dim_index(dim) if dim.class == String
          if (dim == lon_dim && @OPT_sht) then # SHTモジュールを使う
            g = SphericalHarmonics::sh_dlon(g)/@Radius # d/dlamda
          elsif (dim == lat_dim && @OPT_sht) then # SHTモジュールを使う
            g = SphericalHarmonics::sh_dlat(g)/@Radius # d/dmu = cos^-1(phi)*d/dphi
            g = SphericalHarmonics::multiply_cos(g) # cosをかけて d/dphi にする
          else
            g = g.threepoint_O2nd_deriv(dim)
          end
        }
        end


        #### running mean
        if (opt == "OPT_running_mean")
          # g = g.running_mean(rm_dim,rm_span,true)
          # GPhys付属のメソッド running_mean(dim,span,BC,minimal length)
          # BCは 10:SIMPLE, 11:CYCLIC, 12:TRIM
          # minmal length は、移動平均区間に有効なデータがその数以下のときに欠損にする。
          # 以下の設定で、従来と同じ振る舞いになる。
          rm_ary.each{|rm|
            rm_dim = rm[0]; rm_span=rm[1]; rm_skip=rm[2]
            g = g.running_mean(rm_dim,rm_span,10,rm_span)
            unless rm_skip == 0
              thinning = Hash.new; thinning[rm_dim] = {(rm_span/2..(-rm_span/2))=>rm_skip}
              g = g[thinning]
            end
          }
        end

        ## mean along any axis
        if (opt == "OPT_mean")
            @stddev = g.stddev(*dims_mean) if @OPT_overplot_stddev
            g = g.mean(*dims_mean)
          dims_mean.each{|dim|
            # @stddev = @stddev.stddev(dim) if @OPT_overplot_stddev
            # g = g.average(dim)
            # g = g.mean(dim)
          }
        end

        #### global mean
        if (opt == "OPT_global_mean")
          if (true) then
            lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(g, true)
            g.axis(lat_dim).set_pos(VArray.new(NMath::sin(g.axis(lat_dim).pos.val*PI/180.0), nil, "sin_lat") )
            g = g.mean(lon_dim).average("sin_lat") # global mean
            # g = GAnalysis::Planet.ave_s(g)
          else # 経度がない場合
            lat_dim = 0
            g.axis(lat_dim).set_pos(VArray.new(NMath::sin(g.axis(lat_dim).pos.val*PI/180.0), nil, "sin_lat") )
            g = g.average(0) # global mean
          end

       end

        ## standard deviation along any axis
        if (opt == "OPT_stddev")
          dims_stddev.each{|dim|
            g = g.stddev(dim)
          }
        end

        ## max values along any axis
        if (opt == "OPT_max")
          dims_max.each{|dim|
            g = g.max(dim)
          }
        end

        ## min values along any axis
        if (opt == "OPT_min")
          dims_min.each{|dim|
            g = g.min(dim)
          }
        end

        ## deviation from mean along any axis
        if (opt == "OPT_eddy")
          dims_eddy.each{|dim|
            g = g.eddy(dim)
          }
        end

        ## integration along any axis
        if (opt == "OPT_integrate")
          dims_integrate.each{|dim|
            g = g.integrate(dim)
          }
        end

        ## sum along any axis
        if (opt == "OPT_sum")
          dims_sum.each{|dim|
            g = g.sum(dim)
          }
        end

        if (opt == "OPT_addingup")
          dims_addingup.each{|dim|
            g = addingup(g,dim)
          }
        end
      } ## simple operations end



      ## operation of a mathematical function
      if (@OPT_operation)
        p "g = g.#{@OPT_operation}" unless @OPT_silent
        eval "g = g.#{@OPT_operation}"
      end

      if (@OPT_normalize)
        case g.rank
        when 1
          if (@OPT_normalize == "")
            g = g/g[0].val
            g.set_att("long_name",g.get_att("long_name")+" normalized by init" )
            g.units="1"
          elsif (@OPT_normalize == "mean")
            g = g/g.val.mean
            g.set_att("long_name",g.get_att("long_name")+" normalized by mean" )
            g.units="1"
          elsif (@OPT_normalize == "diff")
            g = g - g[5].val
            g.set_att("long_name",g.get_att("long_name")+" diff from init" )
            g.units="1"
          elsif (@OPT_normalize == "absmax")
            g = g/g.val.abs.max
            g.set_att("long_name",g.get_att("long_name")+" diff from init" )
            g.units="1"
          else
            norm = @OPT_normalize.to_f
            g = g/norm
            g.set_att("long_name",g.get_att("long_name")+" normalized by "+norm.to_s)
            g.units="1"
          end
        when 2


        else
          norm = @OPT_normalize.to_f
          g = g/norm
          g.set_att("long_name",g.get_att("long_name")+" normalized by "+norm.to_s)
          g.units="1"

        end
      end


      ### interpolate after operation
      if (@OPT_interpolate) then
        if (@OPT_sig2p) then # for large file
          dim = (GAnalysis::SigmaCoord::find_sigma_d(g) || 2) # if cannot find sigma axis, assume 3rd axis.
          sig = g.axis(dim).to_gphys
          press = GAnalysis::SigmaCoord::sig_ps2p(@PS.cut(g.axnames[g.shape.index(1)]=>i), sig, dim)
          press.units=("Pa"); press.set_att("positive","down")
          g.set_assoc_coords([press])
          # binding.pry
        end

        axes = g.axnames; assoc_coords = (g.assoccoordnames || [] )
        if (@OPT_interpolate.include?("=")) then
          cutaxis, val = @OPT_interpolate.split("=")
          if val.include?(":") then
            vmin, vmax, step = val.split(":")
            vmin = vmin.to_f; vmax = vmax.to_f
            step = 30 unless step # use default steps if not given
            step = step.to_f
            vwidth = (vmax - vmin)/step
            nval = NArray.float(step.to_i).indgen!*vwidth + vmin if NArrayType == "standard"
            nval = NumRu::NArray.float(step.to_i).indgen!*vwidth + vmin if NArrayType == "bigmem"
          elsif(val == "auto")
            # set latter
          else
            val = (eval val)
            if val.class == Array then
              nval = NArray.to_na(val).to_type("sfloat")
            else
              nval = NArray[val].to_type("sfloat")
            end
          end

          if axes.include?(cutaxis) then
            va = g.coordinate(cutaxis)
            va = VArray.new(nval,{"units"=>va.units.to_s},cutaxis)
            g = g.interpolate(va)
            g = g.cut(cutaxis=>val) unless (val.class == Array or val.class == String)
          elsif assoc_coords.include?(cutaxis) then
            va = g.assoc_coords.coord(cutaxis) # assoc_coordsをVArray化
            ax_dim = 0
            if (va.rank >= 2) then # va[true, 0, 0].stddev.valが最大となる軸が対応する軸だと判断する。
              max = 0
              va.rank.times{|i|
                cut_array = va.shape_current.clone.fill(0) ; cut_array[i]=true
                tmp = va[*cut_array].stddev
                if (tmp) then tmp = tmp.val else tmp = 0 end
                if (max < tmp) then
                  max = tmp; ax_dim = i
                end
              }
            end
            if (val == "auto") then
              tmp = NArray.int(va.rank).indgen.to_a if NArrayType == "standard"
              tmp = NumRu::NArray.int(va.rank).indgen.to_a if NArrayType == "bigmem"
              tmp.delete(ax_dim)
              aval = rgn(va.mean(*tmp).val).to_a.uniq # 重複を取り除く
              nval = NArray.to_na(aval) if NArrayType == "standard"
              nval = NumRu::NArray.to_na(aval) if NArrayType == "bigmem"
            end
            va = VArray.new(nval,{"units"=>va.units.to_s,"long_name"=>va.long_name},cutaxis)
            g = g.interpolate(axes[ax_dim]=>va)
            g = g.cut(cutaxis=>val) unless (val.class == Array or val.class == String)
          else
            raise "axis #{cutaxis} is not contained."
          end
        end

      end

      ### add offset
      g = g +  @OPT_offset.to_f if (@OPT_offset)

      ### cut, slice, thinning after operation
      if (@OPT_cut)
        g = cut_slice_thinning(g,@OPT_cut)
      end

      if (@OPT_produce) then
        case @OPT_produce
        when "solar_zenith_angle"
          g = solar_zenith_angle(g)
        when "reflectance"
          g = reflectance(g)
        end


      end

      next g # return this object

  end

  case gp.ntype
  when "sfloat"
    nb = 4
  when "float"
    nb = 8
  else
    nb = 4
  end

  if (gp.get_att("missing_value")) then
    fact = 1.0 # 欠損値を含んでいると、内部で2GB越えしやすい
  else
    fact = 1.0
  end

  if ((gp.length*nb*fact > 2.2E9 || @OPT_parallel )&& !loopdim && NArrayType == "standard") then # for GPhys ogject larger than 2GB
    dims_used = [dims_addingup, dims_deriv, dims_eddy, dims_flip, dims_integrate,
                 dims_max, dims_mean, dims_min, dims_stddev, dims_sum].flatten.compact
    dims_remained = gp.axnames - dims_used
    dims_remained.delete_at(-1) if (dims_remained.length > 1 && gp.coordinate(dims_remained[-1]).length > 100) # 分割数が超えると結合に時間が掛かるので。
    if (@OPT_parallel && dims_remained.include?(@OPT_parallel)) # 陽に指定する
      dims_remained = [@OPT_parallel]
    end
    # multi-proccesses, many memory use.
    np = (ENV['OMP_NUM_THREADS']|| 4).to_i
    print "Required object size is #{(gp.length*nb*fact/1.0E9).round(2)}GB, which is larger than 2GB.\n" unless @OPT_parallel
    print "Operation is processed along #{dims_remained[-1]} dim with #{np} proccesses. This may take a while......\n"
    remaind_dim_val_array = gp.coordinate(dims_remained[-1]).val.to_a

    unless (@OPT_mvo_only or @flag_mvo_gpa) then

      if (gp.data.file.class == NetCDF) then
        if (gp.data.file.format >= 3) then # NetCDF4 対策
          gpary = (0..(remaind_dim_val_array.length-1)).to_a.map{|i|
            gp.cut_rank_conserving(dims_remained[-1]=>remaind_dim_val_array[i]).copy
          }
          @flag_on_memory = true
        end
      end

      gpa = Parallel.map( (0..(remaind_dim_val_array.length-1)).to_a, :in_processes=>np, :progress=>"progress"){|i|
  #        gp, gturl = open_gturl_wildcard(gturl, true)
        if (@flag_mvo_gpa or @flag_on_memory) then # mvo オペレーションした場合は、ファイルから再読み込みしない。
          tmp_gp = gpary[i]
        else
          tmp_gp = reopen(gp).cut_rank_conserving(dims_remained[-1]=>remaind_dim_val_array[i])
        end
        # force axis unit
        if (@OPT_axis_units_name) then
          @OPT_axis_units_name.split(",").each{|it|
            axnum, unit, axname = it.split(":")
            num = axnum.sub("ax","").to_i
            unit = tmp_gp.coordinate(num).units.to_s if unit == "" # 単位未指定時は変更しない
            axname = tmp_gp.axis(num).name unless axname # 軸名未指定時は変更しない
            typecode = tmp_gp.coordinate(num).typecode
            tmp_gp.axis(num).set_pos(VArray.new(tmp_gp.coordinate(num).val.to_type(typecode), {"units"=>unit}, axname))
          }
        end

        if @OPT_set_assoc_coords then
          tmp_assoc_coords = reopen(@assoc_coords).cut_rank_conserving(dims_remained[-1]=>remaind_dim_val_array[i])
          tmp_gp.set_assoc_coords(tmp_assoc_coords)
        end
        proc.call(tmp_gp, remaind_dim_val_array[i])
      }
    end

    begin
      unless (@OPT_nc4a) then
        print "Joining #{gpa.length} GPhys objects to create a single GPhys object of #{gpa[0].length*gpa.length*nb/1.0E9} GB...\n"
        gp = GPhys.join(gpa)
        gary << visualize_and_output(gp)
      end
    rescue
      print "\n  ***  Error Occurred during visualize and output **** \n"
      print "       Operated data size might be too large.\n"
      loop do
        print "       Do you want to output operated data 'gp' as 'tmp.nc' [Y] or start pry [p] or just abort [n]?\n"
        key = $stdin.gets.chomp
        if (key == "Y" || key == "") then
          @OPT_nc4a = "tmp.nc"
          break
        elsif (key == "p") then
          binding.pry
        else
          raise "Abort"
        end
      end
    end

    if (@OPT_nc4a) then
      if (@OPT_nc4a == ""); outfilename = "out.nc"; else; outfilename = @OPT_nc4a; end
      print "large data is being written as #{@OPT_nc4a} ..."
      gpa.size.times{|t|
        if (@OPT_rename || @OPT_long_name || @OPT_unit) then
          #gpa[t] = gpa[t].copy if (gpa[t].data.file != nil || gpa[t].data.file.class == NArray)
          gpa[t].rename(@OPT_rename) if @OPT_rename
          gpa[t].set_att("long_name", @OPT_long_name) if @OPT_long_name
          gpa[t].set_att("long_name", @OPT_rename) if (!@OPT_long_name && @OPT_rename)
          gpa[t].units=@OPT_unit if @OPT_unit
        end

        if (@OPT_axis_units_name) then
          @OPT_axis_units_name.split(",").each{|it|
            axnum, unit, axname = it.split(":")
            num = axnum.sub("ax","").to_i
            unit = gpa[t].coordinate(num).units.to_s if unit == "" # 単位未指定時は変更しない
            axname = gpa[t].axis(num).name unless axname # 軸名未指定時は変更しない
            typecode = gpa[t].coordinate(num).typecode
            gpa[t].axis(num).set_pos(VArray.new(gpa[t].coordinate(num).val.to_type(typecode), {"units"=>unit}, axname))
          }
        end



        auto_write(outfilename, gpa[t].to_type("sfloat"), t, false) # 第４引数は圧縮の有無（圧縮すると読み込みが遅い）
      }
      print " finished\n"
      exit
    end


  elsif loopdim && !(@OPT_scatter or @OPT_color_scatter or @OPT_histogram2D or @OPT_cot or @OPT_rmap or @OPT_fullcolor == "2" or @OPT_vector) then         # animation
    each_along_dims(gp, loopdim){|gp_subset|
      gp_subset = proc.call(gp_subset, 0)
      gary << visualize_and_output(gp_subset)
    }
  else
    gp = proc.call(gp, 0)
    gary << visualize_and_output(gp)
  end


  if (@OPT_textbox) then
    draw_text(@textbox, gary) if (DCL.sgpget("nframe") % (@textbox_fnum) == @textbox_fnum-1) && (@Overplot == 1)
  end


end # end of main loop for given gturl.


if (@OPT_scatter) then
  draw_multi_vars(gary, "scatter")
end

if (@OPT_color_scatter) then
#  gary << landsea.cut(gp_lon,gp_lat)
  draw_multi_vars(gary, "color_scatter")
end

if (@OPT_cot) then
  draw_multi_vars(gary, "contour_over_tone")
end

if (@OPT_histogram2D) then
  draw_multi_vars(gary, "histogram2D")
end

if (@OPT_rmap) then
  draw_multi_vars(gary, "rmap")
end

if (@OPT_fullcolor == "2") then
  draw_multi_vars(gary, "fullcolor2")
end

if (@OPT_vector) then
  draw_multi_vars(gary, "vector")
end


if (@OPT_pry or @OPT_auto_pry) then
  require "pry"
  # ## open work station and setup draw parameters
  # if (@flag_window_open == false)
  #   DCL.swpset('lalt',false) if @OPT_pry
  #   DCL.gropn(@OPT_wsn||1)
  #   draw_setup
  #   @flag_window_open = true
  # end

  #include GGraph
  g = gary[0] if gary.length == 1
  binding.pry if @OPT_pry
  p @OPT_auto_pry
  eval "#{@OPT_auto_pry}" if @OPT_auto_pry
end



# dcl window close
if (@flag_window_open == 1)
  DCL.grcls
  DCL.sgpset("nframe", 0) # フレームの位置を戻す
  # DCL.sgpset("npage", 0) # ページの初期化は不要？

end

if (@outncfile) then
  @outncfile.put_att("command",@comment)
  @outncfile.close
end

if (@OPT_nc4) then
  @OPT_nc4 = "out.nc" if @OPT_nc4 == ""
  write_netcdf4(gary, @OPT_nc4)
end


if (@OPT_file == "eps" && @flag_window_open == 1) #add comment
  fname = DCL.swpget("fname").strip
  DCL.sgpget("npage").times{|i|
    pnum = sprintf("%04d",(i+1).to_s)
    system ("#{CMD_SSED} -i -e '2i %%Comment: #{@comment}' #{fname}_#{pnum}.eps") if CMD_SSED
  }
end

if (@OPT_file == "png" && @flag_window_open == 1) # add comment & crop
  # gif アニメを作成する
  if (@OPT_gif) then
    system ("#{CMD_CONVERT} -delay 50 -loop 0 dcl_????.png dcl.gif") if CMD_CONVERT
  elsif (@OPT_mp4) then
    frate = 2; frate = @OPT_mp4.to_i if (@OPT_mp4 == @OPT_mp4.to_i.to_s)
    system ("#{CMD_FFMPEG} -framerate #{frate} -i dcl_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p dcl.mp4") if CMD_FFMPEG
  else
    # ImageMagic を使って、png に実行コマンドをメタデータとして記述する。
    # exiftool dcl_001.png | grep Comment などで参照可能
    fname = DCL.swpget("fname").strip
    DCL.sgpget("npage").times{|i|
      pnum = sprintf("%04d",(i+1).to_s)
      system ("#{CMD_MOGRIFY}  -comment '#{@comment}' #{fname}_#{pnum}.png") if CMD_MOGRIFY
      system ("#{CMD_MOGRIFY}  -trim +repage #{fname}_#{pnum}.png") if (@OPT_crop && CMD_MOGRIFY)
    }
  end
end

if (@OPT_file == "pdf" && @flag_window_open == 1) # rotate to landscape & add comment & crop
  fname = DCL.swpget("fname").strip
#  system ("/opt/local/bin/pdftk #{fname}.pdf cat 1-endeast output #{fname}-tmp.pdf")
  system ("#{CMD_CPDF} #{fname}.pdf -rotate 90 -o #{fname}-tmp.pdf") if CMD_CPDF
  if (@OPT_crop) then
    system ("#{CMD_PDFCROP} #{fname}-tmp.pdf #{fname}.pdf") if CMD_PDFCROP
    system ("#{CMD_RM}  #{fname}-tmp.pdf" )
  else
    system ("#{CMD_MV} #{fname}-tmp.pdf #{fname}.pdf")
  end
  system ("#{CMD_EXIFTOOL} -Title='#{@comment}' #{fname}.pdf -overwrite_original") if CMD_EXIFTOOL
end

return gary

end

end




# 実行文
if $0 == __FILE__


if ENV['PLANET'] then
  if (ENV['PLANET'].downcase == "venus") then
    GAnalysis::Met::Kappa=0.1914 # for Venus AFES
    GAnalysis::Met::P00=UNumeric[9.2E6,"Pa"]
    GAnalysis::Met::Cp=UNumeric[1.0E3,"J.kg-1.K-1"]
    GAnalysis::Met::R=UNumeric[191.4,"J.kg-1.K-1"] # Gas constants for dry air
    GAnalysis::Met::set_g(UNumeric[8.87,"m.s-2"])
    GAnalysis::Planet::radius=(UNumeric[6.05E6, "m"])
    GAnalysis::Planet::omega=(UNumeric[3.0E-7, "s-1"])

    print "Constants for Venus are used.\n"
    print "Cp    = ", GAnalysis::Met::Cp, "\n"
    print "Cpv   = ", GAnalysis::Met::Cpv, " (Earth)\n"
    print "Kappa = ", GAnalysis::Met::Kappa, "\n"
    print "Lat0  = ", GAnalysis::Met::Lat0, " (Earth)\n"
    print "P00   = ", GAnalysis::Met::P00, "\n"
    print "R     = ", GAnalysis::Met::R, "\n"
    print "Rv    = ", GAnalysis::Met::Rv, " (Earth)\n"
    print "T0    = ", GAnalysis::Met::T0, " (Earth)\n"

    print "g     = ", GAnalysis::Met.g , "\n"
    print "radius= ", GAnalysis::Planet.radius , "\n"
    print "omega = ", GAnalysis::Planet.omega , "\n"
  elsif (ENV['PLANET'].downcase == "mars") then
    GAnalysis::Met::P00=UNumeric[6.1E2,"Pa"]
    GAnalysis::Met::Cp=UNumeric[843.9,"J.kg-1.K-1"] #
    GAnalysis::Met::R=UNumeric[191.0,"J.kg-1.K-1"] # Gas constants for dry air
    GAnalysis::Met::Kappa=GAnalysis::Met::R/GAnalysis::Met::Cp
    GAnalysis::Met::set_g(UNumeric[3.72,"m.s-2"])
    GAnalysis::Planet::radius=(UNumeric[3.389E6, "m"])
    GAnalysis::Planet::omega=(UNumeric[7.09E-5, "s-1"])

    print "Constants for Mars are used.\n"
    print "Cp    = ", GAnalysis::Met::Cp, "\n"
    print "Cpv   = ", GAnalysis::Met::Cpv, " (Earth)\n"
    print "Kappa = ", GAnalysis::Met::Kappa, "\n"
    print "Lat0  = ", GAnalysis::Met::Lat0, " (Earth)\n"
    print "P00   = ", GAnalysis::Met::P00, "\n"
    print "R     = ", GAnalysis::Met::R, "\n"
    print "Rv    = ", GAnalysis::Met::Rv, " (Earth)\n"
    print "T0    = ", GAnalysis::Met::T0, " (Earth)\n"

    print "g     = ", GAnalysis::Met.g , "\n"
    print "radius= ", GAnalysis::Planet.radius , "\n"
    print "omega = ", GAnalysis::Planet.omega , "\n"

  else
    print "Constants for Earth are used.\n"
    print "Cp    = ", GAnalysis::Met::Cp, "\n"
    print "Cpv   = ", GAnalysis::Met::Cpv, "\n"
    print "Kappa = ", GAnalysis::Met::Kappa, "\n"
    print "Lat0  = ", GAnalysis::Met::Lat0, "\n"
    print "P00   = ", GAnalysis::Met::P00, "\n"
    print "R     = ", GAnalysis::Met::R, "\n"
    print "Rv    = ", GAnalysis::Met::Rv, "\n"
    print "T0    = ", GAnalysis::Met::T0, "\n"

    print "g     = ", GAnalysis::Met.g , "\n"
    print "radius= ", GAnalysis::Planet::radius , "\n"
    print "omega = ", GAnalysis::Planet::omega , "\n"


  end
end

GPV.new.main



end


=begin
物理量の計算・変換メモ

=座標変換
* sigma座標 → p座標
  * 要：地表気圧 PS
  gpv Hoge.nc@Hoge --sig2p PS.nc@PS --interpolate p=auto --nc4a  Hoge_p.nc --nodraw

* 気圧座標 p
  * 要：気圧 p
  gpv Hoge.nc@Hoge --set_assoc_coords prs.nc@prs --interpolate prs=auto --nodraw --nc4 Hoge_p.nc

* 高度座標　Z
  * 要：ジオポテンシャルハイト Z
  gpv Hoge.nc@Hoge --set_assoc_coords Z.nc@Z --interpolate Z=auto --nodraw --nc4 Hoge_z.nc

  * 要：温度 T (sigma → z)
  gpv Hoge.nc@Hoge --sig2z T.nc@T --interpolate Z=auto --nc4  Hoge_z.nc --nodraw


=物理量計算

* ジオポテンシャルハイト Z
  * 要：温度 T, sigma座標
  gpv Hoge.nc@Hoge --sig2z T.nc@T --nc4 Z.nc --nodraw  --output_assoc_coords

* 気圧 p
  * 要：地表気圧 PS, sigma座標
  gpv Hoge.nc@Hoge --sig2p PS.ctl@PS --nc4 p.nc --nodraw  --output_assoc_coords

* 密度 RHO
  * 要：温度 T, 気圧座標 p, 気体定数 R
  gpv T.nc@T --mvo "1.0/(x*GAnalysis::Met::R)*ax2" --nc4 RHO.nc --rename "RHO" --unit "kg m-3" --nodraw

  * 要：温度 T, 気圧 p, 気体定数 R
  gpv T.nc@T p.nc@p --mvo "1.0/(x*GAnalysis::Met::R)*y" --nc4 RHO.nc --rename "RHO" --unit "kg m-3" --nodraw

* 鉛直流 w
  * 要：鉛直p速度 OMG, 密度 RHO, 重力加速度 g
  gpv {OMG.nc@OMG,RHO.nc@RHO} --mvo "-x/(y*GAnalysis::Met.g)" --unit "m.s-1" --rename "w" --long_name "vertical velocity" --nodraw --nc4 w.nc

* 絶対角運動量 AngMom
  * 要：東西風速 u (２つ目の次元が 緯度[deg])
  gpv U.nc@U --mvo "(x + @Radius*@Omega*(ax1*D2R).cos)*@Radius*(ax1*D2R).cos" --nodraw --rename AngMom --long_name "Angular momentum" --nc4 AngMom.nc

* 温位 Theta
  * 要：温度 T, 気圧 p
  gpv T.nc@T p.nc@p  --mvo "GAnalysis::Met::temp2theta(x,y)" --nc4 Theta.nc --rename Theta --nodraw --unit "K"

* 温位面上の PV
  * 要：温位座標での 東西風速 U, 南北風速 V, 気圧 p
  gpv U.nc@U V.nc@V --mvo "GAnalysis::Planet.absvor_s(x,y)" --nc AbsVor.nc --nodraw --rename AbsVor --unit "s-1" --long_name "absolute vorticity
  gpv p.nc --derivative theta --nc dp_dtheta.nc --nodraw
  gpv AbsVor.nc@AbsVor dp_dtheta.nc@dp_dtheta --mvo "-x/y*GAnalysis::Met::g" --nc PV.nc --mvo_only --rename "PV" --nodraw


=NetCDFフォーマット変換
  gpv hoge.ctl --nc4a hoge.nc --long_name "hogehoge" --axis_units_name ax2:m:height --unit "hoge"

  * SCALE-GM出力の変換
gpv u.ctl --nc4a nc/u.nc --long_name "zonal velocity" --axis_units_name ax2:m:z --unit "m.s-1"
gpv v.ctl --nc4a nc/v.nc --long_name "meridional velocity" --axis_units_name ax2:m:z --unit "m.s-1"
gpv w.ctl --nc4a nc/w.nc --long_name "vertical velocity" --axis_units_name ax2:m:z --unit "m.s-1"
gpv prs.ctl --nc4a nc/p.nc --long_name "pressure" --rename p --axis_units_name ax2:m:z --unit "Pa"
gpv t.ctl --nc4a nc/T.nc --long_name "temperature" --rename T --axis_units_name ax2:m:z --unit "K"
gpv ps.ctl --nc nc/ps.nc --long_name "surface pressure" --rename ps --axis_units_name ax2:m:z --unit "Pa" --nodraw


=シークエンス








=end
