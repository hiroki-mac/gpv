#--
# =DESCRIPTION
#  Utility methods for gpv.
# =AUTHOR
#  Hiroki Kashimura
#++
class GPV
  def each_along_dims(gphys, loopdim)

    raise ArgumentError,"1st argument must be an GPhys." if !gphys.is_a?(GPhys)
    if loopdim.is_a?(String)
      dimname = loopdim
    elsif
      if loopdim < 0
        dimname = gphys.coord(gphys.rank + loopdim).name
      else
        dimname = gphys.coord(loopdim).name
      end
    else
      raise ArgumentError,"loopdims must consist of Integer and/or String"
    end

    loopdim_na = gphys.coord(dimname).val                      # get coord ary
    loopdim_na = loopdim_na[-1..0] if @OPT_reverse || @OPT_Gr  # reverse
    loopdim_na.each { |x|
      yield( gphys.cut(dimname=>x) )
    }
  end

  def rgn(na) # DCLの関数 RGN** を使って、入力したNArray配列に近いキリのいい数字の配列を返す。
    out = na.collect{|v|
      le = DCL.rgnle(v); ge = DCL.rgnge(v)
      (v - le < ge - v ? le : ge) # より近い方のキリのいい数字を採用する
    }
    return out
  end

  def ymd2time(unit, ymd)
    return ymd.to_f if (ymd.to_i.to_s == ymd or ymd.to_f.to_s == ymd)
    if (unit.include?("since")) then
      unit, sdate = unit.split("since")
      if (unit.include?("day") or unit.include?("month") or unit.include?("year"))
        sdate = Date.parse(sdate)
        ymd += "m1" unless ymd.include?("m")
        ymd += "d1" unless ymd.include?("d")
        ymd = ymd.sub("y","").sub("m","-").sub("d","-")
        vdate = Date.parse(ymd)
        value = date_diff(vdate, sdate, "year" ) if unit.include?("year")
        value = date_diff(vdate, sdate, "month") if unit.include?("month")
        value = date_diff(vdate, sdate, "day"  ) if unit.include?("day")
        return value
      elsif (unit.include?("sec") or unit.include?("min") or unit.include?("hour"))
        sdate = DateTime.parse(sdate)
        ymd += "m1" unless ymd.include?("m")
        ymd += "d1" unless ymd.include?("d")
        ymd = ymd.sub("y","").sub("m","-").sub("d","-")
        vdate = Date.parse(ymd)
        value = date_diff(vdate, sdate, "hour")    if unit.include?("hour")
        value = date_diff(vdate, sdate, "min" )    if unit.include?("min")
        value = date_diff(vdate, sdate, "sec" )    if unit.include?("sec")
        return value
      else
        raise "unsupported time unit\n"
      end
    else
      return ymd.to_f
    end
  end

  def ymd2date(subset)
    # treat y????m?d? like expression.
    if (subset.match(/y([0-9]+)m([0-9]+)d([0-9]+)/)) then
      ymd = subset.scan(/y([0-9]+)m([0-9]+)d([0-9]+)/)
      ymd.each{|tmp|
        subset["y"+tmp[0]+"m"+tmp[1]+"d"+tmp[2]]='Date.parse("'+tmp.join('-')+'")'
      }
    elsif (subset.match(/y([0-9]+)m([0-9]+)/)) then
      ym = subset.scan(/y([0-9]+)m([0-9]+)/)
      ym.length.times{|i|
        date_parse = 'Date.parse("'+ym[i].join('-')+'-01'+'")' unless i == 1
        unless (@OPT_calendar360) then
          date_parse = 'Date.parse("'+ym[i].join('-')+'-28'+'")' if i == 1 && ym[i][1].to_i == 2
          date_parse = 'Date.parse("'+ym[i].join('-')+'-30'+'")' if i == 1 && [4,6,9,11].include?(ym[i][1].to_i)
          date_parse = 'Date.parse("'+ym[i].join('-')+'-31'+'")' if i == 1 && [1,3,5,7,8,10,12].include?(ym[i][1].to_i)
        else
          date_parse = 'Date.parse("'+ym[i].join('-')+'-30'+'")' if i == 1
        end
        subset["y"+ym[i][0]+"m"+ym[i][1]] = date_parse
      }
    elsif (subset.match(/y([0-9]+)/)) then
      y = subset.scan(/y([0-9]+)/)
      y.length.times{|i|
        subset["y"+y[i][0]]='Date.parse("'+y[i][0]+'-01-01'+'")' unless i == 1
        unless (@OPT_calendar360) then
          subset["y"+y[i][0]]='Date.parse("'+y[i][0]+'-12-31'+'")' if i == 1
        else
          subset["y"+y[i][0]]='Date.parse("'+y[i][0]+'-12-30'+'")' if i == 1
        end
      }
    end
    return subset
  end

  def counter(string1, init=0.0, delta=1.0, string2=nil)
    @counter = 0 if @counter == nil
    now = init + delta*@counter
    string = string1 + sprintf("%2.2f",now).to_s + string2

    # for time day and hour
      # binding.pry
      day = (now/24).floor; hour = now%24
      string = string1 + sprintf("%2.0f",day).to_s + " Eday " + sprintf("%2.0f",hour).to_s + " h"

    title(string, 0.5)
    @counter += 1
  end

  def get_calendar(gp, dim)
    gp_axis = gp.axis(dim).to_gphys
    @Calendar = gp_axis.get_att("calendar")
    @Units = gp_axis.get_att("units")
  end

  def set_calendar(gp,dim)
    gp_axis = gp.axis(dim).to_gphys
    cal = gp_axis.get_att("calendar")
    units = gp_axis.get_att("units")

    cal = "360_day" if cal == "uniform30day"

    if (@Calendar == cal && @Units == units) then
      return gp
    elsif (!cal) then
      return gp
    else
      time_array = gp_axis.val
      new_time_array = time_array.map{|i|
        datetime = UNumeric.new(i,units).to_datetime(0.0, cal)
        UNumeric.from_date(datetime,@Units,@Calendar).to_f
        }
      gp.axis(dim).set_pos(VArray.new(new_time_array, {"units"=>@Units}, @Units) )
      return gp
    end
  end

  def date_diff(date1, date2, unit="day")
    if (@OPT_calendar360) then
      day1 = (date1.year-1)*360.0 + (date1.month-1)*30.0 + date1.day
      day2 = (date2.year-1)*360.0 + (date2.month-1)*30.0 + date2.day
      diff = day1 - day2 # diff in days
    elsif (@OPT_calendar365) then
      day1 = (date1.year-1)*365.0 + date1.yday.to_f
      day2 = (date2.year-1)*365.0 + date2.yday.to_f
      day1 = day1 - 1.0 if (date1.leap? && date1.month >= 3) #閏年で3月以降なら2月29日分を除く
      day2 = day2 - 1.0 if (date2.leap? && date2.month >= 3)
      diff = day1 - day2 # diff in days
    else
      diff = (date1 - date2).to_f # diff in days
    end
    case unit
    when "day","days","d"
      return diff
    when "month","months","m"
      raise "not supported yet"
    when "year","years","y"
      raise "not supported yet"
    when "hour","hours","h"
      return diff*24.0
    when "min","minit","minits"
      return diff*24*60.0
    when "sec","second","seconds"
      return diff*24*60.0
    end
      raise "not supported yet"
  end

  # Float()で変換できれば数値、例外発生したら違う
  # cf. http://taro-tnk.hatenablog.com/entry/2012/12/17/001552
  def float_string?(str)
    Float(str)
    true
  rescue ArgumentError
    false
  end

  def common_chars(str_ary)
    flag_scom = true; flag_ecom = true; i = 0; j = -1
    while (flag_scom or flag_ecom)
      flag_scom = str_ary.all?{|s| s[i] == str_ary[0][i]}
      i = i + 1 if flag_scom
      flag_ecom = str_ary.all?{|s| s[j] == str_ary[0][j]}
      j = j - 1 if flag_ecom
    end
    scom = (i == 0) ? nil : str_ary[0][0..(i-1)]
    ecom = (j ==-1) ? nil : str_ary[0][(j+1)..-1]
    return [scom, ecom]
  end

  def common_blocks(str_ary)
    spstr_ary = str_ary.map{|s| s.split("_")}
    length_ary = spstr_ary.map{|s| s.length}
    lmin = length_ary.min
    com_block = []
    lmin.times{|i|
      com_block << i if spstr_ary.all?{|a| a[i]==spstr_ary[0][i]}
    }
    com_block.map!{|i| spstr_ary[0][i]}
    return com_block
  end

  # TODO 複数ファイル・変数への対応
  def fits2gphys(gturl)
  #  require "RubyFits"
    file, num = gturl.split("@")
    f=FitsFile.new(file)
    hdu=f.hdu(num.to_i)
    imageArray=hdu.getAsArray()
    na = NArray.to_na(imageArray).transpose(1,0) # NArray化
    xaxis = Axis.new
    xaxis.pos = VArray.new(NArray.int(hdu.getXSize).indgen!,nil,"x")
    yaxis = Axis.new
    yaxis.pos = VArray.new(NArray.int(hdu.getYSize).indgen!,nil,"y")
    grid = Grid.new(xaxis, yaxis)
    nam = NArrayMiss.to_nam(na, na.eq(na)).set_missing_value(-9.99e33) if na[0].class == Float
    nam = NArrayMiss.to_nam(na, na.eq(na)).set_missing_value(-999)     if na[0].class == Fixnum
    #nam = na
    varray = VArray.new(nam,{"long_name"=>hdu.getName,"units"=>hdu.getHeader(17).toString}, hdu.getName)
    gp = GPhys.new(grid,varray)

    #  binding.pry

    return gp

  end

  # TODO: 複数ファイル・変数への対応
  def csv2gphys(gturl)
    csvfile, vcol = gturl.split("@")
    array = CSV.read(csvfile) # csvの読み込み
    if (vcol.to_i > 0) then # 開く列を数字で指定している場合。
      #ヘッダーの有無を調べ、ある場合は変数名等に使う。ヘッダーは1行目の場合のみ対応。
      unless (float_string?(array[0][0])) then # 数値の文字列じゃない場合はヘッダーだと判断
        axname = array[0][0]; varname = array[0][vcol.to_i]; longname = array[0][vcol.to_i]
        array = array.drop(1)
      else
        axname = "x"; varname = "v"+vcol.to_s; longname ="value"+vcol.to_s
      end
    else # 開く列を変数名で指定している場合は、必ずヘッダーがあるはず。
      array[0].length.times{|i|
        if (array[0][i].strip == vcol) then vcol = i; break; end
      }
      if (vcol.to_i <= 0) then
        raise "Including var names are #{array[0][1..-1].join(", ")}. You can also use natural numbers to specify the column number."
      end
      axname = array[0][0]; varname = array[0][vcol]; longname = array[0][vcol]
      array = array.drop(1)
    end
    array.map!{|j| j.map{|i| i.to_f}} # 要素の数値化
    narray = NArray.to_na(array).transpose(1,0) # NArray化
    xaxis = Axis.new
    xaxis.pos = VArray.new(narray[true,0],nil,axname)
    grid = Grid.new(xaxis)
    varray = VArray.new(narray[true,vcol.to_i],{"long_name"=>longname,"units"=>""},varname)
    gp = GPhys.new(grid,varray)
    return gp
  end


  def output_csv(g)
    if (@OPT_csv == ""); outfilename = "out.csv"; else; outfilename = @OPT_csv; end
    axis0 = g.axis(0).to_gphys # get first axis as gphys object
    axis0_units = axis0.units.to_s # get units of the axis as string
    time_axis = []

    if axis0_units.include?("since") then
      cal_unit, start_date = axis0_units.split("since")
      if (cal_unit.include?("year"))
        start_year = Date.parse(start_date)
        axis0.val.each{|y|
          time_axis << start_year.next_year(y)
        }
      elsif (cal_unit.include?("month"))
        start_month = Date.parse(start_date)
        axis0.val.each{|m|
          time_axis << start_month.next_month(m)
        }
      elsif (cal_unit.include?("day"))
        start_day = Date.parse(start_date)
        axis0.val.each{|d|
          time_axis << start_day.next_day(d)
        }
      elsif (cal_unit.include?("hour"))
        start_hour = DateTime.parse(start_date)
        axis0.val.each{|h|
          time_axis << DateTime.parse(start_hour.time + h*3600)
        }
      elsif (cal_unit.include?("min"))
        start_min = DateTime.parse(start_date)
        axis0.val.each{|m|
          time_axis << DateTime.parse(start_min.time + m*60)
        }

      elsif (cal_unit.include?("sec"))
        start_sec = DateTime.parse(start_date)
        axis0.val.each{|s|
          time_axis << DateTime.parse(start_sec.time + s)
        }
      else
        raise "check the unit of time axis"
      end
      ary0 = time_axis
    else
      ary0 = axis0.val
    end

    ary1 = g.val
    ary1_mask = ary1.valid? if ary1.class == NArrayMiss
    require "csv"
    CSV.open(outfilename, "wb") do |csv|
      ary1.length.times{|i|
        if (ary1.class == NArrayMiss) then
          next unless ary1_mask[i]
        end
        csv << [ary0[i].to_s,ary1[i].to_s]
        }
    end
    print "Operated GPhys object is written as #{outfilename}.\n" unless @OPT_silent
  end


  def auto_write(outfilename, gary, index, compress=false) # 出力ファイル名, Array of GPhys, 時刻のindex, 圧縮の可否
    NetCDF.creation_format=(NetCDF::NC_NETCDF4 | NetCDF::NC_CLASSIC_MODEL) # if compress

    gary = [gary] if gary.class == GPhys
    udname = gary[0].axnames[-1] # UNLIMITED次元の名前
    udval = gary[0].coordinate(udname).val

    if (index == 0) then # 初期設定
      outncfile = NetCDF.create(outfilename)
      gary[0].axnames.each{|n|
        a = gary[0].coordinate(n)
        if (n == udname) then
          outncfile.def_dim(n, 0) # 次元定義 [UNLIMITED次元]
        else
          outncfile.def_dim(n, a.length) # 次元定義
        end
        outncfile.def_var(n, a.ntype, [n]) # 次元変数定義
        a.att_names.each{|m| outncfile.var(n).put_att(m, a.get_att(m), nil) } # 次元変数属性の設定
      }
      gary.each{|g|
        outncfile.def_var(g.name, g.ntype, g.axnames) # 変数定義
        outncfile.var(g.name).deflate(1, true) if compress # set commpression level and shuffle.
        g.att_names.each{|n| outncfile.var(g.name).put_att(n, (g.get_att(n)||""), nil) } # 変数属性の設定
      }

      # set global atts if used in gpv.rb
      if (@sources) then
        if (@sources.length == 1 && @sources[0].include?(".nc")) then # copy the global atts of the source file.
          old_nc = NetCDF.open(@sources[0], "r")
          old_nc.each_att{|at| outncfile.put_att(at.name, at.get) }
        end
        outncfile.put_att("sources", @sources.join(", "))
        outncfile.put_att("command", @commandline )
        outncfile.put_att("working_directory", @current_dir )
        outncfile.put_att("process_date", "#{Time.now}" )
      end

      outncfile.enddef # 定義モード終了
      gary[0].axnames[0..-2].each{|n|
        outncfile.var(n).put(gary[0].coordinate(n).val) # 次元変数の値代入
      }
      outncfile.close # 一旦閉じる
    end

    outncfile = NetCDF.open(outfilename, "a") # 追記モードで開く
    st = outncfile.dim(udname).length
    et = outncfile.dim(udname).length + (udval.size - 1)
    outncfile.var(udname).put(udval,  "start"=>[st], "end"=>[et]) # UNLIMITED次元の値代入

    gary.each{|g|
      sta = Array.new(g.rank-1,0) << st
      eta = Array.new(g.rank-1,-1) << et
      outncfile.var(g.name).put(g.val, "start"=>sta,"end"=>eta) # 変数の値代入
    }
    outncfile.close # 処理の途中でも、別プロセスで可視化できるように毎回閉じる。
  end

  def reopen(gp) # for using in Parallel.map
    file = gp.data.file
    if (!file || file.class == NArray) then # gp is on memory
      return gp
    else # gp is on file
      filepath = file.path
      varname = gp.name
      mapping = gp.data.marshal_dump[1]
      if (mapping) then
        slicer = gp.data.marshal_dump[1].slicer
        return GPhys::IO::open(filepath,varname)[*slicer]
      else
        return GPhys::IO::open(filepath,varname)
      end
    end
  end

  # write out gphys object to a single netcdf ver 4 file.
  # gary is an Array of GPhys, ncfilename is name of output netcdf file.
  # All GPhys objects must have same shape.
  def write_netcdf4(gary, ncfilename, deflate=0, shuffle=true)
    NetCDF.creation_format=(NetCDF::NC_NETCDF4 | NetCDF::NC_CLASSIC_MODEL)
    nc = NetCDF.create(ncfilename, false, false)
    # define dims and axes
    if (@OPT_output_assoc_coords) then
      tmp_ary = []
      gary[0].assoccoordnames.each{|an|
        tmp_ary << gary[0].assoc_coord_gphys(an)
      }
      gary = tmp_ary
    end
  #  binding.pry
    gary.each{|g|
      g.rank.times{|i|
        unless (nc.dim_names.include?(g.coordinate(i).name)) then
          ga = g.coordinate(i)
          nc.def_dim(ga.name, ga.length)
  #        nc.def_var(ga.name, ga.ntype, [nc.dims[-1]])
          nc.def_var(ga.name, "sfloat", [nc.dims[-1]])
          # set axes att
          ga.att_names.each{|an| nc.var(ga.name).put_att(an, ga.get_att(an))}
        end
      }
    }
    # set global atts
    if (@sources.length == 1 && @sources[0].include?(".nc")) then # copy the global atts of the source file.
      old_nc = NetCDF.open(@sources[0], "r")
      old_nc.each_att{|at| nc.put_att(at.name, at.get) }
    end
    nc.put_att("sources", @sources.join(", "))
    nc.put_att("command", @commandline )
    nc.put_att("working_directory", @current_dir )
    nc.put_att("process_date", "#{Time.now}" )

  #  binding.pry
    # define vars
    gary.each{|g|
      g = g.to_type(@OPT_ntype) if @OPT_ntype
      g = g.copy if ( (g.data.file != nil || g.data.file.class == NArray) )# && g.rank != 0 && g.assoccoordnames == nil)
      nc.def_var(g.name, g.ntype, nc.dims(g.axnames))
      nc.var(g.name).deflate(deflate, shuffle) # set commpression level and shuffle.
      if (g.val.class == NArrayMiss) then # to set missing values to output netcdf files.
        if (g.val.get_mask.eq(0).sum > 0) then # skip if there are no missing value points.
          val = g.val
          rmiss = (val.to_na*val.get_mask.eq(0).to_f).max
          rmiss = (val.to_na*val.get_mask.eq(0).to_f).min if rmiss == 0.0 # get missing value
          g.set_att("missing_value", [rmiss])
        end
      end
  #    g.del_att("positive") if (g.get_att("positive") && g.rank != 0)
      # set var att
      g.att_names.each{|an| nc.var(g.name).put_att(an, g.get_att(an))  }
    }
    nc.enddef # change to data models

    # write axes
    nc.each_dim{|d| gary.each{|g|
      gr = g.grid.copy # to avoid error when axis name has changed.
      if (gr.axnames.include?(d.name)) then
        nc.var(d.name).put(gr.coord(d.name).val); break
      end
    }  }

    # write vars
    gary.each{|g|
      nc.var(g.name).put_with_miss(g.val) if g.val.class == NArrayMiss
      nc.var(g.name).put(g.val) if g.val.class == NArray
      print "Operated GPhys object #{g.name} is written in #{ncfilename}.\n" unless @OPT_silent
    }
    nc.close
  end

  def open_gturl_wildcard(gturl, print=false)
  ###   Match wildcard expression for variables if "*" and/or "?" is included.
  ###   Matched variables are added to ARGV array.
  ###   Return gphys object and gturl or false; when false "next" should be executed after this function.
    if (gturl.class == String) then
      @Var_array = [] if @Var_array == nil
      if (@OPT_var =~ /{(.*)}/) then
        @Var_array = $1.split(",")
        @OPT_var = @OPT_var.sub("{#{$1}}","*")
        gturl = gturl + '@' + @OPT_var
      else
        gturl = gturl+'@'+@OPT_var if (@OPT_var && !gturl.include?('@'))
        gturl = gturl+'@*' unless gturl.include?('@')
      end

      # For CSV format
      if (gturl.include?(".csv@")) then
        csvfile, vcol = gturl.split("@")
        if (vcol == "*") then
          array = CSV.read(csvfile); maxcol =  array[0].length - 1
          ARGV.shift
          maxcol.downto(1){|n|
            if (@Var_array==[] or @Var_array.include?(n.to_s) or @Var_array.include?(array[0][n])) then
              ARGV.unshift(csvfile+"@"+array[0][n]) unless float_string?(array[0][n])
              ARGV.unshift(csvfile+"@"+n.to_s) if float_string?(array[0][n])
            end
          }
          return false, gturl
        else
          gp = csv2gphys(gturl)
          return gp, gturl
        end
      end

      # For file format supported by GPhys.
      if (gturl =~ /(,\w.*)/) # convert y????m??d?? expression to Date.parser
        gturl = gturl.sub($1,ymd2date($1))
      end
      file, var, slice, cut_slice, thinning  = parse_gturl(gturl) #GPhys::IO.parse_gturl(gturl)
      file = file[0] if (file.class == Array)
      if (var.class == Array) then
        ARGV.shift
        var.reverse_each{|v|
          if (@Var_array == [])
            next if (@OPT_vf && !file.include?(v))
          else
            next if (!@Var_array.include?(v))
          end
          gturl_splitted =  gturl.split("@")
          gtmpary = gturl_splitted[1].split(",")
          # matched variables are converted to gturl and added to ARGV array.
          ARGV.unshift( gturl_splitted[0] + "@" + gturl_splitted[1].sub(gtmpary[0],v) )
        }
        #next
        return false, gturl
      end
      ###---------------------------------------------
      gturl = find_axisnames(gturl) if gturl.include?(","+AX)
      begin
  #      p "Date.parse(\"0005-02-01\")".to_i
        p gturl
        gp = open_gturl(gturl)
        print "  Reading #{gturl}\n" unless (@OPT_silent || print)

      rescue
        raise "given multiple files cannot be open in one gphys object by #{gturl}. give the variable name.\n"
      end

      gp = invalidation(gp) if @OPT_invalidation
      gp = sig2p(gp) if @OPT_sig2p # change sigma coord to pressure coord
      gp = sig2z(gp) if @OPT_sig2z # change sigma coord to z coord

      gp.set_assoc_coords([@assoc_coords]) if @OPT_set_assoc_coords

      return gp, gturl
    elsif (gturl.class == GPhys) then
      gp = gturl
      gturl = gp.name + "(" + gp.axnames.join(", ") + ")"

      gp = invalidation(gp) if @OPT_invalidation
      gp = sig2p(gp) if @OPT_sig2p # change sigma coord to pressure coord
      gp = sig2z(gp) if @OPT_sig2z # change sigma coord to pressure coord
      gp.set_assoc_coords([@assoc_coords]) if @OPT_set_assoc_coords

      return gp, gturl
    else
      raise "unsupported format."
    end
  end

  def __split_range(range)

    if /(.*):(.*)/ =~ range
      if $1 == ""
        min = nil
      else
        min = $1.to_f
      end
      if $2 == ""
        max = nil
      else
        max = $2.to_f
      end
    elsif range == nil
      min = max = nil
    else
      raise "invalid range: variable subset specification error. split range with ':'\n\n"
    end

    return min, max
  end

  def find_axisnames(gturl)
    file, var, slice, cut_slice, thinning  = parse_gturl(gturl) # GPhys::IO.parse_gturl(gturl)
    axes = GPhys::IO.open(file, var).axnames
    axes.length.times{|i|
  	  gturl = gturl.sub(AX+(i).to_s, axes[i]) if gturl.include?(","+AX+(i).to_s)
  	  }
  	gturl
  end

  ## dimension preparation
  def set_dims(opt_dims)
    dims = (opt_dims).split(/\s*,\s*/)
    dims = dims.map{|dim|
    if (dim.to_i.to_s == dim) then 	dim.to_i
    elsif (dim =~ /#{AX}[0-9]/) then dim[-1].to_i
    else dim
    end
    }
    return dims
  end

  # 粒子群をGPhysオブジェクトで定義する
  # 軸は粒子IDと時間
  # x座標、y座標、z座標を変数とする3つのGPhysを出力する
  # gpは初期時刻（時間軸）を取り出すための流速のGPhys
  # initlocは初期位置を入れたArray [[x0,y0,z0], [x1,y1,z1], ...]
  def set_particles(gp,initlocary)
    if (initlocary[0].size == 3) then flag_3D = true else  flag_3D = false  end
    x_axis = gp.coordinate(0)
    y_axis = gp.coordinate(1)
    z_axis = gp.coordinate(2) if flag_3D

    id_num = initlocary.size
    init_loc_na = NArray.to_na(initlocary)

    x_va = VArray.new(init_loc_na[0,true].reshape(id_num,1), {"long_name"=>x_axis.long_name , "units"=> x_axis.units.to_s}, x_axis.name)
    y_va = VArray.new(init_loc_na[1,true].reshape(id_num,1), {"long_name"=>y_axis.long_name , "units"=> y_axis.units.to_s}, y_axis.name)
    z_va = VArray.new(init_loc_na[2,true].reshape(id_num,1), {"long_name"=>z_axis.long_name , "units"=> z_axis.units.to_s}, z_axis.name) if flag_3D

    id_axis = Axis.new
    id_axis.pos = VArray.new(NArray.int(id_num).indgen,{"long_name"=>"particle ID", "units"=> ""}, "id")

    t_axis = Axis.new
    if (gp.rank > 2 && gp.axnames[-1].include?('t') ) then # 時間軸（最後の軸と仮定）があるかどうかの判定
      t_gp = gp.coord(-1)
      t_axis.pos = VArray.new(t_gp[0].val, {"long_name"=>t_gp.long_name, "units"=>t_gp.units.to_s}, t_gp.name )
    else
      t_axis.pos = VArray.new(NArray.sfloat(1).indgen, {"long_name"=>"time", "units"=>"s"}, "time" )
    end

    grid = Grid.new(id_axis, t_axis)

    print "#{id_axis.length} particles were set.\n"

    return GPhys.new(grid,x_va), GPhys.new(grid,y_va), GPhys.new(grid,z_va) if flag_3D
    return GPhys.new(grid,x_va), GPhys.new(grid,y_va), nil

  end

  def self.join(ary)
    if (ary.size <= 110) then 
      return GPhys.join(ary)
    else 
      n = 0; tmp_ary = []
      while (n*100 < ary.size) 
        range = (n*100)..([(n+1)*100-1, ary.size-1].min)
        tmp_ary << ary[range]; n = n + 1
      end
      np = (ENV['OMP_NUM_THREADS']|| 4).to_i
      tmp_ary = Parallel.map(tmp_ary, :in_processes=>np, :progress=>"progress"){|i| GPhys.join(i) }
      return GPhys.join(tmp_ary)
    end
  end

  # def fit_rank(gp, gpref)
  #   gprank = gp.rank
  #   if (gp.class == GPhys) then
  #
  #   elsif (gp.class == NArrayMiss) then
  #
  #   elsif (gp.class == NArray) then
  #     when gprank
  #     case 1
  #       gp.shape[0]
  #     case 2
  #     end
  #   end
  # end

end
