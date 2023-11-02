#--
# data arrange methods for gpv.
#++

class GPV

  # Interpolate to Gaussian latitudes.
  def gglat(gp)
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp)
    nlat = SphericalHarmonics.gauss_lat(@OPT_GLAT.to_i)
    vlat = VArray.new(nlat,{"units"=>"degree_north"},gp.axnames[lat_dim])
    gp = gp.interpolate(lat_dim=> vlat)
    if (lon_dim) then
      nlon = NArray.float(vlat.length*2)
      (vlat.length*2).times{|i| nlon[i] = 360.0/(vlat.length*2)*i }
      vlon = VArray.new(nlon, {"units"=>"degree_east"},gp.axnames[lon_dim])
      gp = gp.cyclic_ext(lon_dim).interpolate(lon_dim=> vlon)
    end
    return gp
  end

  # Regrid to given horizontal resolition with Gaussian latitudes
  def regrid(gp)
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp)
    nlat = SphericalHarmonics.gauss_lat(@OPT_regrid.to_i)
    vlat = VArray.new(nlat,{"units"=>"degree_north"},gp.axnames[lat_dim])
    nlon = NArray.float(vlat.length*2)
    (vlat.length*2).times{|i| nlon[i] = 360.0/(vlat.length*2)*i }
    vlon = VArray.new(nlon, {"units"=>"degree_east"},gp.axnames[lon_dim])
    gp = interpolate_pole(gp) if @OPT_interpolate_pole
    gp = gp.interpolate(lat_dim=> vlat)
    if (gp.axis(lon_dim).to_gphys.val[0] == 0.0) then
      gp = gp.cyclic_ext(lon_dim).interpolate(lon_dim => vlon)
    else
      lon_axis = gp.axis(lon_dim)
      v0 = lon_axis.to_gphys.val[0]; v1 = lon_axis.to_gphys.val[-1]
      gp = gp.cut(lon_axis.name=> -v0..v1).interpolate(lon_dim => vlon)
    end
  end

  # split axis
  def split_axis(gp)
    sa_dim, sa_span, sa_axisname = (@OPT_split_axis).split(/\s*,\s*/)
    if sa_dim.to_i.to_s == sa_dim then
      sa_dim = sa_dim.to_i
    end
    sa_span = sa_span.to_i

    tmp = Array.new
    newaxis = Array.new
    sa_span.times{|i|
      tmp << gp[{sa_dim=>{(i..-1)=>sa_span}}]
      newaxis << (i+1).to_f
    }
    sa_axisname = "newaxis" unless sa_axisname
    gp = GPhys.concat(tmp,newaxis,sa_axisname)
    return gp
  end

  # 次元 dim を反転させる（データ格納の南北方向が逆の場合などに使う）
  def reverse(gp,dim)
    raise "dim must be given by integer: 0,1,2,..." if dim.class != Fixnum
    # 値を反転させる
    val_nam = gp.val
    if (val_nam.class == NArrayMiss) then
      new_mask = val_nam.get_mask.reverse(dim)
      new_na = val_nam.to_na.reverse(dim)
      new_nam = NArrayMiss.to_nam(new_na, new_mask)
    else
      new_nam = val_nam.reverse(dim)
    end

    new_gp = gp.copy
    new_gp.replace_val(new_nam)

    # 軸を反転させる
    new_pos = new_gp.axis(dim).pos
    new_axis_na = new_gp.axis(dim).to_gphys.val.reverse(0)
    new_pos.replace_val(new_axis_na)
    new_axis = Axis.new
    new_axis.pos = new_pos
    new_grid = new_gp.grid.change_axis(dim,new_axis)
    new_gp = GPhys.new(new_grid, new_gp.data)

    return new_gp
  end

  # 緯度座標軸に極の上 -/+90 deg の格子点を付け加える。データの値は元々の最南/最北点の経度平均値とする。
  ### これだとベクトル量には使ってはいけない。スカラー量にのみ使える。
  def interpolate_pole(gp)
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp)
    unless (gp.axis(lat_dim).to_gphys.val[0].abs == 90.0) then
      packed_from_south = (gp.axis(lat_dim).to_gphys.val[0] < 0) ? true : false #緯度の格納順を調べる
      # 緯度軸を拡張する
      new_axis_va = gp.axis(lat_dim).to_gphys.data
      if (packed_from_south) then
        new_axis_ar = new_axis_va.val.to_a.unshift(-90.0).push(90.0)
      else
        new_axis_ar = new_axis_va.val.to_a.unshift(90.0).push(-90.0)
      end

      new_axis_na = NArray.to_na(new_axis_ar)
      new_axis_va = VArray.new(new_axis_na,{"units"=>"degree_north"},gp.axnames[lat_dim])
      # interpolateメソッドでgphysオブジェクトを拡張する
      gp = gp.interpolate(lat_dim=> new_axis_va)
      # interpolate メソッドでは極で欠損値になるので、改めて値を設定する
      val = gp.val
      sp_val = val[true, 1, false].mean(0) # lon, latの次元の位置を仮定している
      np_val = val[true,-2, false].mean(0)
      gp.axis(lon_dim).length.times{|i|
        val[i, 0, false] = sp_val
        val[i,-1, false] = np_val
      }
      gp.replace_val(val)
      return gp
    else
      return gp
    end
  end

  # Extend the data backward.
  def extend_backward_along(gp,dim,len,value=0.0)
    dim = gp.axnames[dim.to_i]  if (dim.to_i.to_s == dim)
    ax_val = gp.coord(dim).copy.val
    ext_gp = gp.cut(dim=>ax_val[0]..ax_val[len-1])
    dif = ax_val[len] - ax_val[0]
    ext_gp.axis(dim).set_pos(ext_gp.coord(dim)-dif)
    ext_gp.replace_val(ext_gp.val*0.0+value)
    return GPhys.join([ext_gp,gp])
  end

  # Extend the data cyclically.
  def cyclic_extension(gp,dim,n=1)
    dim = gp.axnames[dim.to_i]  if (dim.to_i.to_s == dim || dim.class == Fixnum)
    ndim = gp.axnames.index(dim)
    gaxis = gp.coordinate(dim)
    dx = (gaxis[1..-1] - gaxis[0..-2]).val.mean
    if (n==1) then # 長さ 1 のVArray/GPhysは仕様が複雑なので、以下のように回避する
      gal = gaxis[0..1]; gar = gaxis[-2..-1]
      sl = gp.shape.clone.fill(true); sl[ndim]=(-2..-1)
      sr = gp.shape.clone.fill(true); sr[ndim]=(0..1)
      sc = gp.shape.clone.fill(true); sc[ndim]=(1..-2)
      gple = gp[*sl];  gpre = gp[*sr]
      gple.axis(ndim).set_pos(gal - dx*2)
      gpre.axis(ndim).set_pos(gar + dx*2)
      return GPhys.join([gple,gp,gpre])[*sc]
    else
      gal = gaxis[0..(n-1)]; gar = gaxis[-n..-1]
      sl = gp.shape.clone.fill(true); sl[ndim]=(-n..-1)
      sr = gp.shape.clone.fill(true); sr[ndim]=(0..(n-1))
      gple = gp[*sl];  gpre = gp[*sr]
      gple.axis(ndim).set_pos(gal - dx*n)
      gpre.axis(ndim).set_pos(gar + dx*n)
      return GPhys.join([gple,gp,gpre])
    end
  end

  # Thinning the data
  def data_thinning(gp, num)
    att = gp.data.attr_copy # copy attributes
    step = gp.first2D.shape.max/num
    return gp if (step == 0)
    shape = gp.shape; index = []
    shape.size.times{|i|
      if (i > 1) then index << true
      elsif (shape[i] < shape[0]/step/2) then index << true
      elsif (shape[i]/step > 1) then
        index << (NArray.int(shape[i]/step).indgen!*step).to_a
      else index << true
      end
    }
    case index.size
    when 1
      val = gp.val[index[0]]; grid = gp.grid[index[0]]
    when 2
      val = gp.val[index[0],index[1]]; grid = gp.grid[index[0],index[1]]
    when 3
      val = gp.val[index[0],index[1],index[2]]; grid = gp.grid[index[0],index[1],index[2]]
    else
      val = gp.val[index[0],index[1],index[2],index[3],false]
      grid = gp.grid[index[0],index[1],index[2],index[3],false]
    end
    print "NOTE: Input data is thinned to #{grid.shape} (from #{gp.shape}).\n" unless @OPT_silent
    print "      Use --nothinning to prevent data thinning.\n" unless @OPT_silent
    return GPhys.new(grid, VArray.new(val, att, gp.name))
  end

  def invalidation(gp)
    rmiss = -999.0
    miss_val = (eval @OPT_invalidation)
    val = gp.val
    val = NArrayMiss.to_nam(val) if (val.class == NArray)
    if (miss_val.class == Float or miss_val.class == Fixnum) then
      mask = val.eq(miss_val).eq(0) #note mask must be 0 or 1; other value causes strange value
    elsif (miss_val.class == Range) then
      mask = val.ge(miss_val.min)*val.le(miss_val.max).eq(0)
    else
      raise "invalidation value must be given by Float (1.0) or Fixnum (1) or Range (0..1)"
    end
    val.set_mask(mask)
    new_gp = gp.copy ; new_gp.replace_val(val)
  #  new_gp.set_att("missing_value",[miss_val])
    return new_gp
  end

  def cut_slice_thinning(gp, var)
    slice = Hash.new
    cut_slice = Hash.new
    thinning = Hash.new
    var_descr = var.split(/,/)
    var_descr.each do |s|
      if /(.*)=(.*)/ =~ s
        dimname = $1
        subset = $2
        subset = ymd2date(subset)      # convert y????m??d?? expression to Date.parser
        case subset
        when /\^(.*):(.*):(.*)/
          slice[dimname] = (eval $1)..(eval $2)
          thinning[dimname] = {(0..-1) => $3.to_i}
        when /\^(.*):(.*)/
          slice[dimname] = (eval $1)..(eval $2)
        when /\^(.*)/
          slice[dimname] = (eval $1)
        when /(.*):(.*):(.*)/
          cut_slice[dimname] = (eval $1)..(eval $2)
          thinning[dimname] = {(0..-1) => $3.to_i}
        when /(.*):(.*)/
          cut_slice[dimname] = (eval $1)..(eval $2)
        else
          cut_slice[dimname] = (eval subset)
        end
        else
          raise "invalid URL: variable subset specification error\n\n" +
          "URL format: " + GTURLfmt
        end
      end
      slice = nil if slice.length == 0
      cut_slice = nil if cut_slice.length == 0
      thinning = nil if thinning.length == 0

      gp = gp[slice] if slice
      gp = gp.cut(cut_slice) if cut_slice
      gp = gp[thinning] if thinning
      return gp
  end

  def make_bnd_grid(grid, num, ztype=nil) # ztype should be "height", "theta", "sigma", or "pressure"
    z = grid.axis(num) # z軸を取り出す
    # z軸が、高度かどうか判定
    if (ztype == nil) then
      if (["m", "km"].include?(z.pos.units.to_s)) then ztype = "height"
      elsif (["K"].include?(z.pos.units.to_s)) then ztype = "theta"
      else ztype = "sigma" end
    end
    z.pos.replace_val(z.pos.val.log) if (ztype == "sigma" or ztype == "pressure") # σ 軸 or p 軸
    new_z = Axis.new(true,false,z.name+"_bnd") # cell typeの軸を作成
    new_z.set_cell_guess_bounds(z.pos).set_pos_to_bounds
    bound_val = new_z.cell_bounds.val
    new_z.pos.name=(z.name + "_bnd")

    # 間隔が狭すぎる or 間隔の変化が大きいと(?)失敗する。その場合は、単純内挿（平均）で置き換える。
    if (bound_val != bound_val.sort && bound_val[-1..0] != bound_val.sort) then 
      z_val = z.to_gphys.val
      bound_val[1..-2] = (z_val[0..-2] + z_val[1..-1])/2.0
    end

    # 上端・下端の値の修正
    if (ztype == "height") then
      bound_val[0]  = [0.0, bound_val[0]  ].max
    elsif (ztype == "sigma")
      bound_val = bound_val.exp
      bound_val[0]  = [1.0, bound_val[0]  ].min
      bound_val[-1] = [0.0, bound_val[-1] ].max
    elsif (ztype == "pressure")
      bound_val = bound_val.exp
    end
    new_z.cell_bounds.replace_val(bound_val)
    new_grid = grid.change_axis(num, new_z)
    return new_grid
  end


  def take_first_valid(gp, dim, rev=false)
    ary = Array.new(gp.rank).fill(true)
    ary[dim] = 0

    new_gp = gp[*ary]*0.0
    new_val = new_gp.val.to_na

    val = gp[*ary].val
    if (val.class != NArrayMiss) then
      raise "Given GPhys object has no missing value attribut (is not NArrayMiss)"
    end
    cum_mask = gp[*ary].val.get_mask*0 # 初期値 全部ゼロ

    num = gp.shape[dim]
    num.times{|k|
      if (rev == false) then
        ary[dim] = k
      else
        ary[dim] = (num - 1) - k
      end
      val = gp[*ary].val
      cum_mask = cum_mask + val.get_mask
      mask = cum_mask.eq(1) # 1以外はゼロのマスク
      new_val = new_val + (val.to_na)*mask
      if (cum_mask.eq(0).max == 0) then
        break
      end
    }
    new_gp.replace_val(NArrayMiss.to_nam(new_val,cum_mask))
    return new_gp
  end


  def calc_ground_altitude(gp,topo,dim=2)
    z = gp.coord(dim)
    gz_na = NArray.float(*(gp.shape))
    ranks = [*0..(gp.rank-1)] # 転置用のindex配列
    ranks.delete_at(dim)
    ranks.unshift(dim)
    gz_na = gz_na.transpose(*ranks) + z.val
    ranks = [*0..(gp.rank-1)]
    ranks.delete_at(0)
    ranks.insert(dim,0)
    gz_na = gz_na.transpose(*ranks) - topo.val.to_na
    gz = GPhys.new(gp.grid.copy,VArray.new(gz_na,{"long_name"=>"ground altitude","units"=>"m"},"GZ"))
    return gz
  end




end
