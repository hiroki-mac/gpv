#!/usr/bin/env ruby
##################################################
=begin rdoc
= NAME
gpv - quick viewer and manipulater for the values of a variable specified by a gtool4-type URL.

gpv is based on and inspired by gpview, which is developed by Dr. Shin-ichi Takehiro and other dcmodel developers.

= AUTHOR
Hiroki Kashimura (hiroki@gfd-dennou.org)

= USAGE
  gpv hoge.nc@hoge

  see https://github.com/hiroki-mac/gpv
=end

#################################################
if File.ftype($0) == 'link' then
  link = File.readlink($0)
  path = File.expand_path(link, File.dirname($0))
else
  path = File.expand_path($0)
end
GPV_DIR = File.dirname(path) + "/"

require GPV_DIR+"gpv_config.rb"
require GPV_DIR+"gpv_options.rb"
require GPV_DIR+"gpv_analysis.rb"
require GPV_DIR+"gpv_arrange.rb"
require GPV_DIR+"gpv_utils.rb"
require GPV_DIR+"gpv_visualize.rb"
require GPV_DIR+"gpv_gphysmod.rb"
require GPV_DIR+"gpv_spherical_harmonics_next.rb"
require GPV_DIR+"gpv_lagrange.rb"
require GPV_DIR+"gpv_sequence.rb"


#require "profile"



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
  parser.set_options(*OPTIONS_for_GetoptLong)
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
  if (@OPT_help != "") then
    OPTIONS.each{|opt|
      if (opt.include?("--"+@OPT_help)) then
        print "\n= OPTION\n"
        opt.each{|v|
          print "\n   " if (v == opt[-1])
          print v
          print " "
          }
        print "\n"
        exit(1)
      end
      }
      print "\n ERROR: option --#{@OPT_help} was not found.\n\n"
  end
  print "\n= OPTIONS\n"
  print "Available OPTIONS are follows:\n"
  OPTIONS.each{|opt| opt.each{|v|
    print "#{v}, " if (v.to_s.start_with?("--"))
  }}
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

## exec ncatted
if (@OPT_edit_ncatt) then
  var, att, val = @OPT_edit_ncatt.split(":")
  flag_delete = false
  if (val == nil) then
    flag_delete = true
  elsif (val == val.to_i.to_s) then
    type = "l" # long
  elsif (val == val.to_f.to_s) then
    type = "f" # float
  else
    type = "c" # character
  end

  ARGV.each{|v|
    file = v.split("@")[0]
    if (flag_delete == false) then
      system ("ncatted -a #{att},#{var},c,#{type},#{val} #{file}") # create att if there is not
      system ("ncatted -a #{att},#{var},m,#{type},#{val} #{file}") # modify att
      print "att:#{att} of var:#{var} in #{file} was created/modified to #{val}.\n"
    else
      system ("ncatted -a #{att},#{var},d,, #{file}")                # delete att
      print "att:#{att} of var:#{var} in #{file} was deleted.\n"

    end
  }
  exit
end

## exec sequence scripts
if (@OPT_sequence) then
  sequence = @OPT_sequence.split(",")
  sequence.each{|key|
    print "sequence \"#{key}\" is processed.\n"
    begin
      eval "GPVSequence.#{key}"
    rescue
      raise "sequence key:#{key} is wrong...\n"
    end
    print "sequence \"#{key}\" has ended.\n"
  }
  exit
end






## alias options
if (@OPT_wm) then
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
    DCL.sgscmn(@OPT_clrmap.to_i||63)
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
if (@OPT_land or @OPT_ocean)# or @OPT_color_scatter)
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
    if @OPT_replace_axis.include?("@")
      axnum, newaxis_gturl = @OPT_replace_axis.split(":")
      num = axnum.sub("ax","").to_i
      newaxis_gp, newaxis_gturl = open_gturl_wildcard(newaxis_gturl)
      gp.axis(num).set_pos(VArray.new(newaxis_gp.val, {"units"=>newaxis_gp.units.to_s}, newaxis_gp.name))
      print "axis #{num} was replaced by values of #{newaxis_gturl}"
    else
      axnum, newaxis = @OPT_replace_axis.split(":")
      num = axnum.sub("ax","").to_i
      axary = newaxis.split(",")
      axary.map!{|a| a.to_f}
      unit = gp.coordinate(num).units.to_s
      axname = gp.axis(num).name
      gp.axis(num).set_pos(VArray.new(NArray.to_na(axary), {"units"=>unit}, axname))
      print "#{axname}'s value was replaced by #{axary}"
    end
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
      axis_val = (eval "axis_val#{ope}") if ope
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
    thres = 720
    if (@OPT_thinning) then
      gp = data_thinning(gp, @OPT_thinning.to_i)
    elsif (gp.shape.length >= 2 && gp.first2D.shape.max > thres*2) then
      @OPT_nocont = true if (!@OPT_noshade)

      case gp.ntype
      when "sfloat" then  nb = 4
      when "float" then nb = 8
      else nb = 4
      end

      if (gp.length*nb > 2.2E9) then
        np = (ENV['OMP_NUM_THREADS']|| 4).to_i
        dname = gp.axnames[-1]
#        print "Operation is processed along #{dname} dim with #{np} proccesses. This may take a while...\n"
        gpa = Parallel.map(gp.coordinate(-1).val.to_a, :in_processes=>1, :progress=>"progress"){|i|
#          gp_subset = reopen(gp).cut_rank_conserving(dname=>i)
          gp_subset = gp.cut_rank_conserving(dname=>i)
          gp_subset = data_thinning(gp_subset, thres)
        }
        gp = GPhys.join(gpa)

      else
        gp = data_thinning(gp, thres)
      end
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


  ## calculate ground altitude from z-coordinate and topography data and set it as associate coordinate
  if (@OPT_set_GZ_ac) then
    topo = open_gturl_wildcard(@OPT_set_GZ_ac)[0]
    dim = gp.axnames.find_index("lev")
    dim = 2 unless dim
    gp.set_assoc_coords([calc_ground_altitude(gp,topo)])
  end

  # interpolate the value on the poles (i.e., lat = -90/+90 deg by the mean of values on surrounding grid points)
  gp = interpolate_pole(gp) if @OPT_interpolate_pole

  gp.set_att("units","") if @OPT_ignoreunit


  # preparation to use spherical_harmonics_next module
  if (@OPT_sht or @OPT_ES or @OPT_uvcomp) then
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, true)
    nmax = gp.coordinate(lon_dim).length/3
    if (@OPT_ES == "output_pmn") then
      lat_na = gp.coordinate(lat_dim).val
      lat_na.length.times{|j|
        gpmn, gdpmn = SphericalHarmonics::output_pmn(gp.coordinate(lon_dim).val,lat_na[j..j],nmax,gp)
        auto_write("pmn.nc",  [gpmn, gdpmn],  j)
        p j
      }
      p "pmn and dpmn were written in pmn.nc"
      exit
    end
    gw, pmn, dpmn = SphericalHarmonics::sh_init(gp.coordinate(lon_dim).val, gp.coordinate(lat_dim).val, nmax, gp, nil)
  end

  if (@OPT_zonal_shift) then
    gp = zonal_shift_on_sphere(gp, @OPT_zonal_shift)
  end

  if (@OPT_wnf_analysis) then
    gp = wnf_analysis(gp)
  end



  ## for case of calculating difference of two gphys object
  if (@OPT_diff or @OPT_srot or @OPT_sdiv or @OPT_divide or @OPT_HKE or @OPT_ES or @OPT_uvcomp)
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
    elsif (@OPT_uvcomp)
      print "  Calculating #{@OPT_uvcomp}."
      comp_ary = @OPT_uvcomp.split(",")
      gpa = SphericalHarmonics::sh_uvcomps(prev_gp, gp, @Radius, comp_ary)
      if (gpa.length > 1) then
        ARGV << gpa[1..-1]; ARGV.flatten!
      end
      gp = gpa[0]
      @OPT_uvcomp = nil

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
  #####################################################
  ## for case of perfoming Lagrange tracer advection ##
  #####################################################
  if (@OPT_particle_advection) then
    u = gp
    v = open_gturl(ARGV[0])
    w = open_gturl(ARGV[1]) if ARGV[1]

    # defaul values
    initloc = []
    t_div = 24
    np = 1
    proc = 0
    out_div = 1

    @OPT_particle_advection.split(" ").each{|a|
      if (a.include?("=")) then
        case a.split("=")[0]
        when "tdiv"
          t_div = a.split("=")[1].to_i
        when "outdiv"
          out_div = a.split("=")[1].to_i
        when "parallel"
          np = a.split("=")[1].to_i
        when "proc"
          proc = a.split("=")[1].to_i
        end
      elsif (a.include?("@")) then

      else
        x_loc, y_loc, z_loc = a.split(",")
        if (x_loc.include?(":")) then x0, x1, xn = x_loc.split(":"); xn = 10 unless xn
        else x0 = x_loc; x1 = x_loc; xn = 1 end
        if (y_loc.include?(":")) then y0, y1, yn = y_loc.split(":"); yn = 10 unless yn
        else y0 = y_loc; y1 = y_loc; yn = 1 end
        if (z_loc.include?(":")) then z0, z1, zn = z_loc.split(":"); zn = 10 unless zn
        else z0 = z_loc; z1 = z_loc; zn = 1 end
        x0 = x0.to_f; x1 = x1.to_f; xn = xn.to_i ; dx = (x1 - x0)/[xn-1,1].max
        y0 = y0.to_f; y1 = y1.to_f; yn = yn.to_i ; dy = (y1 - y0)/[yn-1,1].max
        z0 = z0.to_f; z1 = z1.to_f; zn = zn.to_i ; dz = (z1 - z0)/[zn-1,1].max
        zn.times{|k| yn.times{|j| xn.times{|i|
              initloc << [x0+dx*i, y0+dy*j, z0+dz*k]
        } } }
      end
    }
    p_lon, p_lat, p_z = set_particles(u,initloc) # 粒子オブジェクト（初期位置）生成
    print "(#{p_lon.val[0]}, #{p_lat.val[0]}, #{p_z.val[0]}) t = 0\n"

    # 粒子群を分割する
    slen = p_lon.length/np
    p_loc = []
    np.times{|n|
      range = (n*slen)..((n+1)*slen-1)
      range = (n*slen)..(-1) if (np-1 == n)
      p_loc << [p_lon[range,true], p_lat[range,true], p_z[range,true]]
    }

    if (u.rank > 2 && u.axnames[-1].include?('t') ) then # 時間軸（最後の軸と仮定）があるかどうかの判定
      t_axis = u.coord(-1)
      t_len = t_axis.length
      flag_t_exist = true
    else
      t_len = 100 # デフォルト値
      flag_t_exist = false
    end

    # 流速場をメモリ上に展開（2GB制限あり）
    begin
      u = u.copy; v = v.copy; w = w.copy
      print "u, v, w were copied on memory\n"
    rescue
      print "Sizes of u, v, w exceed 2GB, so that flow variables were NOT copied on memory\n"
      filename = u.axis(0).inspect.split("'")[3]
      if (filename) then
        if filename.include?(".nc") then
          dl = `ncdump -sh #{filename} | grep DeflateLevel`
          if (!dl.empty?) then
            print dl, "\n"
            print "CAUTION: Deflated NetCDF4 files require many computational costs.\n"
            print "         Converting files to NetCDF3 is recomended.\n"
          end
        end
      end
    end


    if (p_z) then # 3次元移流
      auto_write("particles_#{proc}.nc", p_loc[proc], 0, false) # 初期位置の記録
      p_lon_s = p_loc[proc][0]; p_lat_s = p_loc[proc][1];  p_z_s = p_loc[proc][2]
      LagrangeSphere.init_3D(u)
    else # 水平2次元
      auto_write("particles_#{proc}.nc", [p_lon,p_lat], 0, false) # 初期位置の記録
      LagrangeSphere.init(u)
    end
    id0 = p_lon_s.coord(0).val[0]
    id1 = p_lon_s.coord(0).val[-1]
    # 流速データの時間間隔の時刻ループ
    (t_len-1).times{|t|
      if (flag_t_exist) then
        data_dt = (t_axis[t+1] - t_axis[t]).to_type(5)
        dt = data_dt/t_div
      else
        data_dt = UNumeric[1.0,"day"]
        dt = data_dt/t_div
        u_now = u; v_now = v
      end

      # 流速データの時間間隔を t_div 分割して計算する。
      t_div.times{|tn|
        if (p_z) then # 3次元移流
          if (flag_t_exist) then
            u_now = u[true,true,true,t..(t+1)] # (u[true,true,true,t]*(t_div-n) + u[true,true,true,t+1]*n)/t_div
            v_now = v[true,true,true,t..(t+1)] # (v[true,true,true,t]*(t_div-n) + v[true,true,true,t+1]*n)/t_div
            w_now = w[true,true,true,t..(t+1)] # (w[true,true,true,t]*(t_div-n) + w[true,true,true,t+1]*n)/t_div
          end
          p_lon_s, p_lat_s, p_z_s = LagrangeSphere.particle_advection_3D(p_lon_s,p_lat_s,p_z_s,u_now,v_now,w_now,dt,0,t_div,tn)
        else # 水平2次元移流
          if (flag_t_exist) then
            u_now = (u[true,true,t]*(t_div-tn) + u[true,true,t+1]*tn)/t_div
            v_now = (v[true,true,t]*(t_div-tn) + v[true,true,t+1]*tn)/t_div
          end
          p_lon, p_lat = LagrangeSphere.particle_advection_2D(p_lon,p_lat,u_now,v_now,dt)
        end
        # auto_write("particles.nc", [p_lon,p_lat], t+1, false)
        auto_write("particles_#{proc}.nc", [p_lon_s,p_lat_s,p_z_s], t+1, false) if (tn+1)%(t_div/out_div) == 0
#        print "(#{p_lon_s.val[0]}, #{p_lat_s.val[0]}, #{p_z_s.val[0]}) t = #{t+1}, p = #{proc} (id:#{id0}-#{id1})\n"

      }
      if (p_z) then
#        auto_write("particles_#{proc}.nc", [p_lon_s,p_lat_s,p_z_s], t+1, false)
        print "(#{p_lon_s.val[0]}, #{p_lat_s.val[0]}, #{p_z_s.val[0]}) t = #{t+1}, p = #{proc} (id:#{id0}-#{id1})\n"
      else
        auto_write("particles.nc", [p_lon,p_lat], t+1, false)
        print "(#{p_lon.val[0]}, #{p_lat.val[0]}) t = #{t+1}\n"
      end
    }
    exit!
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
          p "data is being copied on memory"
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
          x = reopen(xx).cut_rank_conserving(da_name=>da_val_ary[i])
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
        @flag_mvo_gpa = true; @OPT_parallel = true unless @OPT_parallel
      else
        print "combining #{gpa.size} gphys objects... this may take very long time."
        gp = GPV.join(gpa)
        @OPT_parallel = nil
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

  ## apply lowpass filter (FFT)
  gp = lowpass(gp) if (@OPT_lowpass && !@OPT_sht)
  ## apply highpass filter (FFT)
  gp = highpass(gp) if (@OPT_highpass && !@OPT_sht)

  ## apply lowpass filter (SphericalHarmonics)
  gp = lowpass_sht(gp) if (@OPT_lowpass && !@OPT_highpass && @OPT_sht)
  ## apply lowpass filter (SphericalHarmonics)
  gp = highpass_sht(gp) if (@OPT_highpass && !@OPT_lowpass && @OPT_sht)
  ## apply bandpass (highpass and lowpass) filter (SphericalHarmonics)
  gp = bandpass_sht(gp) if (@OPT_highpass && @OPT_lowpass && @OPT_sht)

  ## apply 2D median filter
  gp = median_filter2D(gp,@OPT_median_filter2D.to_i) if @OPT_median_filter2D


  ## interpolate to gausian latitude and correspond longitude
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
            num = 1
            dims_mean.each{|d| num = num * g.coord(d).length }
            if (g.ntype == "sfloat")
              if (num > 1E6) then
                raise "ERROR: Data number for mean is too large and loss of trailing digits occours. Try with --ntype float."
              elsif (num > 1E5) then
                print "WARNING: Data number for mean is large and prescision may be very low.\n"
                print "         Consider to use --ntype float.\n"
              end
            end
            @stddev = g.stddev(*dims_mean) if @OPT_overplot_stddev

            if (@OPT_rank_conserving)
              tmp = g.cut_rank_conserving(dims_mean[0]=>0)
              tmp.replace_val(g.mean(dims_mean[0]).val.reshape(*tmp.shape))
              tmp.coordinate(dims_mean[0]).replace_val(NArray.to_na([g.coordinate(dims_mean[0]).val.mean]))
              g = tmp
            else
              g = g.mean(*dims_mean)
            # dims_mean.each{|dim|
              # @stddev = @stddev.stddev(dim) if @OPT_overplot_stddev
              # g = g.average(dim)
              # g = g.mean(dim)
            # }
            end
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

          cutaxis, targetaxis = cutaxis.split("[")
          targetaxis.sub!("]","") if (targetaxis)

          if axes.include?(cutaxis) then
            va = g.coordinate(cutaxis)
            va = VArray.new(nval,{"units"=>va.units.to_s},cutaxis)
            g = g.interpolate(va)
            g = g.cut(cutaxis=>val) unless (val.class == Array or val.class == String)
          elsif assoc_coords.include?(cutaxis) then
            va = g.assoc_coords.coord(cutaxis) # assoc_coordsをVArray化
            ax_dim = 0
            if (targetaxis) then
              axes.size.times{|i| ax_dim = i if axes[i] == targetaxis }
            elsif (va.rank >= 2) then # va[true, 0, 0].stddev.valが最大となる軸が対応する軸だと判断する。
              max = 0
              va.rank.times{|i|
                cut_array = va.shape_current.clone.fill(0) ; cut_array[i]=true
                tmp = va[*cut_array].stddev
                if (tmp) then tmp = tmp.val else tmp = 0 end
                if (max < tmp) then
                  max = tmp; ax_dim = i
                end
              }
              print "NOTE: interpolation with #{cutaxis} is performed along #{axes[i]}. if this is NOT correct, please specify the correct dim by assoc_coord[dim] (ex. p[z])\n"
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
#            binding.pry
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

  end # end of proc
#---------------------------------------------------
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
    if (@OPT_parallel && @OPT_parallel != "" && @OPT_parallel != true) # 陽に指定する
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
        gp = GPV.join(gpa)
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

        if (@OPT_output_assoc_coords) then
          if (@OPT_sig2p) then
            a = gpa[t].axnames.map{|ax|
              if (@PS.axnames.include?(ax))
                if (gpa[t].coord(ax).length == @PS.coord(ax).length); true; else; t; end
              else
                nil
              end
            }
            a.delete(nil)
            g = sig2p(gpa[t],@PS.cut_rank_conserving(*a))
            g = g.assoc_coord_gphys(g.assoccoordnames[0])
          else
            g = gpa[t].assoc_coord_gphys(gpa[t].assoccoordnames[0])
          end
          auto_write(outfilename, g.to_type("sfloat"), t, false) # 第４引数は圧縮の有無（圧縮すると読み込みが遅い）
        else
          auto_write(outfilename, gpa[t].to_type("sfloat"), t, false) # 第４引数は圧縮の有無（圧縮すると読み込みが遅い）
        end
      }
      print " finished\n"
      exit
    end


  elsif loopdim && !(@OPT_scatter or @OPT_color_scatter or @OPT_histogram2D or @OPT_cot or @OPT_rmap or @OPT_fullcolor == "2" or @OPT_vector) then         # animation
    if (@OPT_anim_div) then
      gp_ld = gp.coord(loopdim); gp_ld_val = gp_ld.val; ln = gp_ld_val.length
      t_div = @OPT_anim_div.to_i
      (ln-1).times{|n|
        t_div.times{|m|
          # 線形補間 (あまり綺麗に見えない)
#          gp_subset = (gp.cut_rank_conserving(loopdim=>gp_ld_val[n])*(t_div-m) + gp.cut_rank_conserving(loopdim=>gp_ld_val[n+1])*m )/t_div
          gp_subset = gp.cut_rank_conserving(loopdim=>gp_ld_val[n])
          ld_now = (gp_ld[n]*(t_div-m) + gp_ld[n+1]*m)/t_div
          gp_subset.axis(loopdim).set_pos(ld_now)
          gp_subset = proc.call(gp_subset.cut(loopdim=>ld_now.val[0]), 0)
          gary << visualize_and_output(gp_subset)
          plot_particles(gp_subset) if @OPT_plot_particles
        }
      }
    else
      @LOOP_COUNTER = 0
      each_along_dims(gp, loopdim){|gp_subset|
        gp_subset = proc.call(gp_subset, 0)
        gary << visualize_and_output(gp_subset)
        plot_particles(gp_subset) if @OPT_plot_particles
        @LOOP_COUNTER += 1
      }
    end
  elsif (@OPT_zoomin) then
    x_dist, y_dist, zoom_type = @OPT_zoomin.split(",")
    x_axis, x_distrange = x_dist.split("="); y_axis, y_distrange = y_dist.split("=")
    # ズーム先の範囲
    x_min,  x_max = x_distrange.split(":") ; y_min,  y_max = y_distrange.split(":")
    x_min = x_min.to_f; x_max = x_max.to_f ; y_min = y_min.to_f; y_max = y_max.to_f
    # 初期範囲
    ix_min = gp.coord(x_axis).min.to_f     ; ix_max = gp.coord(x_axis).max.to_f
    iy_min = gp.coord(y_axis).val.min      ; iy_max = gp.coord(y_axis).val.max
    divnum = 300 - 1 # 分割数
    if (zoom_type == "linear") then
      # 線形にズームイン
      dx_min = (ix_min - x_min)/divnum ; dx_max = (ix_max - x_max)/divnum
      dy_min = (iy_min - y_min)/divnum ; dy_max = (iy_max - y_max)/divnum
      (divnum+1).times{ |n|
        @OPT_xrange = (ix_min - dx_min*n).to_s + ":" + (ix_max - dx_max*n).to_s
        @OPT_yrange = (iy_min - dy_min*n).to_s + ":" + (iy_max - dy_max*n).to_s
        visualize_and_output(gp)
      }
    else
      # 比でズームイン
      x_center = (x_min + x_max)*0.5 ; y_center = (y_min + y_max)*0.5
      r_xmax = ((x_max - x_center)/(ix_max - x_center))**(1.0/divnum)
      r_ymax = ((y_max - y_center)/(iy_max - y_center))**(1.0/divnum)
      r_xmin = ((x_min - x_center)/(ix_min - x_center))**(1.0/divnum)
      r_ymin = ((y_min - y_center)/(iy_min - y_center))**(1.0/divnum)
      (divnum+1).times{ |n|
        @OPT_xrange = (x_center-(x_center-ix_min)*r_xmin**n).to_s + ":" + (x_center+(ix_max-x_center)*r_xmax**n  ).to_s
        @OPT_yrange = (y_center-(y_center-iy_min)*r_ymin**n).to_s + ":" + (y_center+(iy_max-y_center)*r_ymax**n  ).to_s
        visualize_and_output(gp)
      }
    end
  else # 1変数、通常時
    gp = proc.call(gp, 0)
    gary << visualize_and_output(gp)
    plot_particles(gp) if @OPT_plot_particles
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

if (@OPT_particle_advection) then
  u = gary[0]; v = gary[1]
  w = gary[2] if gary[2]
  t_div   = 24 # defaul value
  initloc = []

  @OPT_particle_advection.split(" ").each{|a|
    if (a.include?("=")) then
      t_div = a.split("=")[1].to_i
    elsif (a.include?("@")) then

    else
      x_loc, y_loc, z_loc = a.split(",")
      if (x_loc.include?(":")) then x0, x1, xn = x_loc.split(":"); xn = 10 unless xn
      else x0 = x_loc; x1 = x_loc; xn = 1 end
      if (y_loc.include?(":")) then y0, y1, yn = y_loc.split(":"); yn = 10 unless yn
      else y0 = y_loc; y1 = y_loc; yn = 1 end
      if (z_loc.include?(":")) then z0, z1, zn = z_loc.split(":"); zn = 10 unless zn
      else z0 = z_loc; z1 = z_loc; zn = 1 end
      x0 = x0.to_f; x1 = x1.to_f; xn = xn.to_i ; dx = (x1 - x0)/[xn-1,1].max
      y0 = y0.to_f; y1 = y1.to_f; yn = yn.to_i ; dy = (y1 - y0)/[yn-1,1].max
      z0 = z0.to_f; z1 = z1.to_f; zn = zn.to_i ; dz = (z1 - z0)/[zn-1,1].max
      zn.times{|k| yn.times{|j| xn.times{|i|
            initloc << [x0+dx*i, y0+dy*j, z0+dz*k]
      } } }
    end
  }


  p_lon, p_lat, p_z = set_particles(u,initloc) # 粒子オブジェクト（初期位置）生成
  print "(#{p_lon.val[0]}, #{p_lat.val[0]}, #{p_z.val[0]}) t = 0\n"

  if (p_z) then # 3次元移流
    auto_write("particles.nc", [p_lon,p_lat,p_z], 0, false) # 初期位置の記録
    LagrangeSphere.init_3D(u)
  else # 水平2次元
    auto_write("particles.nc", [p_lon,p_lat], 0, false) # 初期位置の記録
    LagrangeSphere.init(u)
  end

  if (u.rank > 2 && u.axnames[-1].include?('t') ) then # 時間軸（最後の軸と仮定）があるかどうかの判定
    t_axis = u.coord(-1)
    t_len = t_axis.length
    flag_t_exist = true
  else
    t_len = 100 # デフォルト値
    flag_t_exist = false
  end
  # 以下 parallel 用
  # results = Parallel.map([0..161,162..323,324..485,486..-1], :in_processes=>4){|i| #このブロックが並列処理される
  #   p_lon_s = p_lon[i,0]  # サブセットを取り出す
  #   p_lat_s = p_lat[i,0]  # サブセットを取り出す
  #   p_z_s = p_z[i,0]  # サブセットを取り出す

  (t_len-1).times{|t|
    if (flag_t_exist) then
      data_dt = (t_axis[t+1] - t_axis[t]).to_type(5)
      dt = data_dt/t_div
    else
      data_dt = UNumeric[1.0,"day"]
      dt = data_dt/t_div
      u_now = u; v_now = v
    end

    t_div.times{|n|
      if (p_z) then # 3次元移流
        if (flag_t_exist) then
          u_now = (u[true,true,true,t]*(t_div-n) + u[true,true,true,t+1]*n)/t_div
          v_now = (v[true,true,true,t]*(t_div-n) + v[true,true,true,t+1]*n)/t_div
          w_now = (w[true,true,true,t]*(t_div-n) + w[true,true,true,t+1]*n)/t_div
        end
        if (false) then # parallel用 ← 並列処理にしても速くならなかった
          p_lon_s, p_lat_s, p_z_s = LagrangeSphere.particle_advection_3D(p_lon_s,p_lat_s,p_z_s,u_now,v_now,w_now,dt)
        else
          p_lon, p_lat, p_z = LagrangeSphere.particle_advection_3D(p_lon,p_lat,p_z,u_now,v_now,w_now,dt)
        end
      else # 水平2次元移流
        if (flag_t_exist) then
          u_now = (u[true,true,t]*(t_div-n) + u[true,true,t+1]*n)/t_div
          v_now = (v[true,true,t]*(t_div-n) + v[true,true,t+1]*n)/t_div
        end
        p_lon, p_lat = LagrangeSphere.particle_advection_2D(p_lon,p_lat,u_now,v_now,dt)
      end
      # auto_write("particles.nc", [p_lon,p_lat], t+1, false)
    }
    if (p_z) then
      auto_write("particles.nc", [p_lon,p_lat,p_z], t+1, false)
      print "(#{p_lon.val[0]}, #{p_lat.val[0]}, #{p_z.val[0]}) t = #{t+1}\n"
    else
      auto_write("particles.nc", [p_lon,p_lat], t+1, false)
      print "(#{p_lon.val[0]}, #{p_lat.val[0]}) t = #{t+1}\n"
    end

  }
  # 以下 parallel 用
  # [p_lon_s,p_lat_s, p_z_s]
  # }
  # results = results.transpose
  # p_lon = GPhys.join(results[0])
  # p_lat = GPhys.join(results[1])
  # p_z = GPhys.join(results[2])


end



# dcl window close
if (@flag_window_open == 1)
  DCL.grcls
  DCL.sgpset("nframe", 0) # フレームの位置を戻す
  # DCL.sgpset("npage", 0) # ページの初期化は不要？

end

if (@outncfile) then
  @outncfile.put_att("command",@comment)
  @outncfile.put_att("title", "see command")
  @outncfile.put_att("comment", "see command")
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

* 地表温位 theta_s
  * 要：地表温度 Ts, 地表気圧 ps
  gpv Ts.nc@Ts ps.nc@ps  --mvo "GAnalysis::Met::temp2theta(x,y)" --nc4 Theta_s.nc --rename Theta_s --nodraw --unit "K"

* 温位面上の PV
  * 要：温位座標での 東西風速 U, 南北風速 V, 気圧 p
  gpv U.nc@U V.nc@V --mvo "GAnalysis::Planet.absvor_s(x,y)" --nc AbsVor.nc --nodraw --rename AbsVor --unit "s-1" --long_name "absolute vorticity
  gpv p.nc --derivative theta --nc dp_dtheta.nc --nodraw
  gpv AbsVor.nc@AbsVor dp_dtheta.nc@dp_dtheta --mvo "-x/y*GAnalysis::Met::g" --nc PV.nc --mvo_only --rename "PV" --nodraw

* 地表面応力 tau_s
  * 要：温度 t、気圧 prs、地表面温度 ts、東西風 u 、南北風 v in z-座標
  gpv {t.nc@t,prs.nc@prs} --mvo "GAnalysis::Met::temp2theta(x,y)" --nc Theta.nc --rename Theta --nodraw --unit "K" --ntype sfloat --nodraw

  gpv ts.nc@ts,time=^-1 ps.nc@ps,time=^-1 --mvo "GAnalysis::Met::temp2theta(x,y)" --nc Theta_s.nc --rename Theta_s --nodraw --unit K --ntype sfloat --nodraw

  gpv Theta.nc@Theta,lev=0 Theta_s.nc@Theta_s,lev=0 u.nc@u,lev=0,time=^-1 v.nc@v,lev=0,time=^-1 rho.nc@rho,lev=0,time=^-1 --mvo "surface_stress(x,y,z,u,v)" --rename tau_s --nc tau_s.nc --nodraw



=NetCDFフォーマット変換
  gpv hoge.ctl --nc4a hoge.nc --long_name "hogehoge" --axis_units_name ax2:m:height --unit "hoge"

  * SCALE-GM出力の変換
  gpv u.ctl --nc4a nc/u.nc --long_name "zonal velocity" --axis_units_name ax2:m:z --unit "m.s-1"
  gpv v.ctl --nc4a nc/v.nc --long_name "meridional velocity" --axis_units_name ax2:m:z --unit "m.s-1"
  gpv w.ctl --nc4a nc/w.nc --long_name "vertical velocity" --axis_units_name ax2:m:z --unit "m.s-1"
  gpv prs.ctl --nc4a nc/p.nc --long_name "pressure" --rename p --axis_units_name ax2:m:z --unit "Pa"
  gpv t.ctl --nc4a nc/T.nc --long_name "temperature" --rename T --axis_units_name ax2:m:z --unit "K"
  gpv ps.ctl --nc nc/ps.nc --long_name "surface pressure" --rename ps --axis_units_name ax2:m:z --unit "Pa" --nodraw









=end
