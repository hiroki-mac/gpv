#--
# =DESCRIPTION
#  Modification of GPhys and GGraph methods, and NArray.
# =AUTHOR 
#  Hiroki Kashimura
#++

class Float
  def floor # Infinity にも floor メソッドを追加する。
    if (self.abs == Float::INFINITY) then
      return self
    else
      return super()
    end
  end
end

class NArray

  def to_na
    self
  end

  def log10; NMath::log10(self);  end
  def log;   NMath::log(self) ;   end
  def exp;   NMath::exp(self) ;   end
  def sqrt;  NMath::sqrt(self) ;  end
  def log2;  NMath::log2(self) ;  end
  def sin;   NMath::sin(self) ;   end
  def cos;   NMath::cos(self) ;   end
  def tan;   NMath::tan(self) ;   end
  def sinh;  NMath::sinh(self) ;  end
  def cosh;  NMath::cosh(self) ;  end
  def tanh;  NMath::tanh(self) ;  end
  def asin;  NMath::asin(self) ;  end
  def acos;  NMath::acos(self) ;  end
  def atan;  NMath::atan(self) ;  end
  def asinh; NMath::asinh(self) ; end
  def acosh; NMath::acosh(self) ; end
  def atanh; NMath::atanh(self) ; end
  def cot;   NMath::cot(self) ;   end
  def csc;   NMath::csc(self) ;   end
  def sec;   NMath::sec(self) ;   end
  def coth;  NMath::coth(self) ;  end
  def csch;  NMath::csch(self) ;  end
  def sech;  NMath::sech(self) ;  end
  def acot;  NMath::acot(self) ;  end
  def acsc;  NMath::acsc(self) ;  end
  def asec;  NMath::asec(self) ;  end
  def acoth; NMath::acoth(self) ; end
  def acsch; NMath::acsch(self) ; end
  def asech; NMath::asech(self) ; end
end


class GPhys
  # Join multiple GPhys objects (not need for any pre-ordering).
  #
  # ARGUMENT
  # * gpnarray [Array (or 1D NArray) of GPhys]
  # 
  def GPhys.join(gpary, ignore_overlap=false)

  #< initialization with the first GPhys object >
    gp = gpary[0]
    rank = gp.rank
    gpstore = MDStorage.new(rank)
    gpstore[ *Array.new(rank, 0) ] = gp     # first element
    x0s = (0...rank).collect{|d|
      pos = gp.axis(d).pos
      x0 = UNumeric[ pos.val[0].round(12), pos.units ] # round(12)を追加した
      [ x0 ]   # first values of each coordinate
    }

    #< scan the coordiantes of the remaining GPhys objects >
    for k in 1...gpary.length
      gp = gpary[k]
      idx = Array.new
      for d in 0...rank
        pos = gp.axis(d).pos
        x0 = UNumeric[ pos.val[0].round(12), pos.units ] # round(12)を追加した
        i = x0s[d].index(x0)
        if i.nil?
          x0s[d].push(x0)
          i = x0s[d].length-1
        end
        idx.push(i)
      end
      gpstore[*idx] = gp
    end

    if !ignore_overlap && gpstore.count_non_nil != gpary.length
      raise(ArgumentError,"Cannot uniquely locate one or more objects; some overlap in the grids?")
    end

    gpnary = gpstore.to_na

    #< Sort along dimensions to join >
    gpnary = __sort_gpnary(gpnary)

    #< Join! >
    self.join_md_nocheck(gpnary)
  end
end



class GPV

  include GGraph

  #class_function
def GGraph::annotate(str_ary, noff=nil)
  lclip = DCL.sgpget('lclip')
  DCL.sgpset('lclip',nil)
  lnum = 0
  str_ary.each{ |str|lnum += 1 }
    charsize = 0.7 * DCL.uzpget('rsizec1')
  dvx = 0.01
  dvy = charsize*1.5
  raise TypeError,"Array expected" if ! str_ary.is_a?(Array)
  vxmin,vxmax,vymin,vymax = DCL.sgqvpt
  #vx = 0.70
  #vy = 0.0045 + (lnum-1)*dvy
  vx = [0.5 + 0.3*(vxmax-vxmin), 0.7].min
  vy = [vymin-0.1, 0].max + (lnum-1)*dvy

  if @@nannot > 1 then # ３つ目以降の annotation は図の右上に書く
    vx = vxmax + 0.01
    vy = vymax - charsize/2
  end
  str_ary.each{|str|
    DCL::sgtxzv(vx,vy,str,charsize,0,-1,DCL.sgiget('index'))
    vy -= dvy
    @@nannot += 1
  }
  DCL.sgpset('lclip',lclip)
  nil
end







def parse_gturl(gturl)
  if /(.*)@(.*)/ =~ gturl
    file = $1
    var = $2
  elsif /(.*)\/(.*)/ =~ gturl
    file = $1
    var = $2
  else
    raise "invalid URL: '[@|/]' between path & variable is not found\n\n" +
    "URL format: " + GTURLfmt
  end
  if /[\*\?]/ =~ file ## match file names if wildcard expression included.
    raise "\n Any files did not match the given expression: \""+file+\
    "\" \n" if Dir[file].empty?
    file = Dir[file].sort
  elsif (!File.exist?(file))
    raise "File specified by gturl \"#{gturl}\" was not found"
  end
  if /,/ =~ var
    slice = Hash.new
    cut_slice = Hash.new
    thinning = Hash.new
    var_descr = var.split(/,/)
    var = var_descr.shift
    var_descr.each do |s|
      if /(.*)=(.*)/ =~ s
        dimname = $1
        subset = $2
        case subset
        when /\^(.*):(.*):(.*)/
          slice[dimname] = (eval $1)..(eval $2)
          thinning[dimname] = {(0..-1) => $3.to_i}
        when /\^(.*):(.*)/
          slice[dimname] = (eval $1)..(eval $2)
        when /\^(.*)/
          slice[dimname] = $1.to_i
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
  else
    slice = nil
    cut_slice = nil
    thinning = nil
  end
  if /[\*\?]/ =~ var ## match var names if wildcard expression included.
    var_reg = var.gsub("*",".*").gsub("?",".") # convert to regular exp.
    case file
    when String
      vars = GPhys::IO.var_names_except_coordinates(file)
    when Array
      vars_t = file.collect{|f| GPhys::IO.var_names_except_coordinates(f)}
      vars = vars_t.flatten.uniq
    end
    vars_matched = vars.select{|v| v.match(/#{var_reg}/)}
    if (vars_matched.empty?)
      raise "\n Any variables in \"#{file}\" did not match the given
      expression: \"#{var}\".\n Included variable(s): #{vars.join(", ")}."
    end
    # put the variable name by String if the number of matched variable is
    # one, otherwise by Array of Strings.
    var = vars_matched.length == 1 ? vars_matched[0] : vars_matched
  end

  if (slice) then 
    axnames = GPhys::IO.open(file,var).axnames
    slice.each_key{|k| 
      unless (axnames.include?(k)) then 
        p "NOTE: #{k}-axis does not exist in #{file}; ignored."
        slice.delete(k) 
      end
    }
    slice = nil if (slice.length == 0)
  end
  if (cut_slice) then 
    axnames = GPhys::IO.open(file,var).axnames
    cut_slice.each_key{|k| 
      unless (axnames.include?(k)) then 
        p "NOTE: #{k}-axis does not exist in #{file}; ignored."
        cut_slice.delete(k)
      end
    }
    cut_slice = nil if (cut_slice.length == 0)
  end


  [file, var, slice, cut_slice, thinning]
end   # def parse_gturl

def open_gturl(gturl)
  file, var, slice, cut_slice, thinning = parse_gturl(gturl)
  if var.is_a?(Array)
    raise "This method treats a gturl of a single GPhys object. " +
    "Use open_multi_gturl to treat a gturl of multiple objects."
  end
  gp = GPhys::IO.open(file,var)
  if @OPT_rank_conserving then
    slice = slice.map{|k, v| [k, gp.coord(k).val[v]] }.to_h if slice
    gp = gp.cut_rank_conserving(slice) if slice
    gp = gp.cut_rank_conserving(cut_slice) if cut_slice
  else
    gp = gp[slice] if slice
    gp = gp.cut(cut_slice) if cut_slice
  end
  if (thinning &&  @sources.to_s.include?(".ctl")) then
    gp = gp.copy ; gp = gp[thinning]
  end
  gp
end   # def open_gturl




# File ../../lib/numru/ganalysis/covariance.rb, line 8
def covariance(gphys0, gphys1, *dims)
  unless GPhys===gphys0 && GPhys===gphys1
    raise "gphys0 and gphys1 must be GPhys"
  end
  unless gphys0.shape == gphys1.shape
    raise "gphys0 and gphys1 must have the same shape"
  end
  units = gphys0.units*gphys1.units
  if dims.length == 0
    dims = Array.new
    gphys0.rank.times{|i| dims.push i }
  else
    dims = dims.map{|dim| gphys0.dim_index(dim) }
  end
  val0 = gphys0.val
  val1 = gphys1.val
  if val0.is_a?(NArrayMiss)
    if val1.is_a?(NArrayMiss)
      mask = val0.get_mask * val1.get_mask
      ndiv = mask.to_type(NArray::LINT).accum(*dims)
      val0 = val0.set_mask(mask)
      val1 = val1.set_mask(mask)
    else
      ndiv = val0.get_mask.to_type(NArray::LINT).accum(*dims)
      val1 = NArrayMiss.to_nam(val1,val0.get_mask)
    end
  elsif val1.is_a?(NArrayMiss)
    ndiv = val1.get_mask.to_type(NArray::LINT).accum(*dims)
    val0 = NArrayMiss.to_nam(val0,val1.get_mask)
  else
    ndiv = 1
    gphys0.shape.each_with_index{|s,i|
      ndiv *= s if dims.include?(i)
    }
  end
  val0 -= val0.accum(*dims).div!(ndiv)
  val1 -= val1.accum(*dims).div!(ndiv)
  nary = val0.mul_add(val1,*dims)
  if Float === nary
    ndiv = ndiv[0] if ndiv.is_a?(NArray)
    nary /= (ndiv-1)
    return UNumeric.new(nary, units), ndiv
  else
    nary = nary/(ndiv.sum(*dims)-1)
    vary = VArray.new(nary,
                      {"long_name"=>"covariance","units"=>units.to_s},
                      "covariance")
    new_grid = gphys0.grid.delete_axes(dims, "covariance").copy
    return GPhys.new(new_grid,vary), ndiv
  end
end

def corelation(gphys0, gphys1, *dims)
  val0 = gphys0.val
  val1 = gphys1.val
  if val0.is_a?(NArrayMiss)
    mask = val0.get_mask
  else
    mask = NArray.byte(*(val0.shape)).fill!(1)
  end
  if val1.is_a?(NArrayMiss)
    mask2 = val1.get_mask
  else
    mask2 = NArray.byte(*(val1.shape)).fill!(1)
  end
  mask.mul!(mask2)
  val0 = NArrayMiss.to_nam(val0) unless val0.is_a?(NArrayMiss)
  val1 = NArrayMiss.to_nam(val1) unless val1.is_a?(NArrayMiss)
  val0 = val0.set_mask(mask)
  val1 = val1.set_mask(mask2)
  p val0,val1 if @DEBUG
  gphys0 = gphys0.copy.replace_val(val0)
  gphys1 = gphys1.copy.replace_val(val1)

  covar, ndiv = covariance(gphys0, gphys1,*dims)
#  return covar/(gphys0.stddev(*dims)*gphys1.stddev(*dims)), mask.to_type(NArray::LINT).sum(*dims)
  cor = covar/(gphys0.stddev(*dims)*gphys1.stddev(*dims))
  cor.rename("correlation")
  cor.set_att("long_name","R(#{gphys0.name},#{gphys1.name}) along "+gphys0.axnames[*dims])
  return cor
end

end


class GPhys
  # 鉛直軸の種類と軸番号を判別する。
  def self.vertical_axis_kind?(gp)
    if (GAnalysis::Met.find_prs_d(gp)) then
      zd = GAnalysis::Met.find_prs_d(gp)
      type = "pressure"
    elsif (GAnalysis::SigmaCoord.find_sigma_d(gp)) then
      zd = GAnalysis::SigmaCoord.find_sigma_d(gp)
      type = "sigma"
    else
      zd_name = gp.axnames.select{|name| ["lev","level","altitude","z","height"].include?(name.downcase)}
      if (zd_name.empty?) then
        zd = nil; type = nil
      else
        zd = gp.axnames.find_index(zd_name[0])
        type = "height"
      end
    end
    return [zd, type]
  end

  # dn軸を最初の軸にする。
  def bring_axis_to_first(dn)
    dn = self.axnames.find_index(dn) if (dn.class == String)
    return self if (dn == 0)
    dims = []
    self.rank.times{|n| dims << n} # [0,1,2,3,..]
    dims.delete(dn)
    return self.transpose(dn,*dims)
  end

  # 最初の軸を dn番目にする。
  def bring_firstaxis_to(dn)
    return self if dn == 0
    dims = []
    self.rank.times{|n|
      if (n < dn) then dims << n+1 elsif (n == dn) then dims << 0 else dims << n end
    }
    return self.transpose(*dims)
  end

  # 鉛直微分。鉛直軸を自動判別して対応する。
  def dz(temp=nil)
    zd, type = GPhys.vertical_axis_kind?(self)
    g = GAnalysis::Met.g # 重力加速度
    r = GAnalysis::Met::R # 乾燥大気の気体定数
    case type
    when "height"
      return self.threepoint_O2nd_deriv(zd)
    when "pressure" # d/dz = - g*p/(R*temp) d/dp
      raise "Temperature is required for conversion; give as .dz(temp)." if (temp == nil)
      prs = self.coord(zd)
      return ((-g)/(r*temp.bring_axis_to_first(zd))*prs).bring_firstaxis_to(zd) * self.threepoint_O2nd_deriv(zd)
    when "sigma" # dA/dz = - g*sig/(R*temp)*dA/dsig
      sig = self.coord(zd)
      raise "temperature is required for conversion; give as .dz(temp)." if (temp == nil)
      return ((-g)/(r*temp.bring_axis_to_first(zd))*sig).bring_firstaxis_to(zd) * self.threepoint_O2nd_deriv(zd)
    when nil
      raise "vertical axis was not found."
    end
  end


end
