#--
# modification of GPhys and GGraph methods. and NArray
#++

class NArray
  def to_na
    self
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
  [file, var, slice, cut_slice, thinning]
end   # def parse_gturl

def open_gturl(gturl)
  file, var, slice, cut_slice, thinning = parse_gturl(gturl)
  if var.is_a?(Array)
    raise "This method treats a gturl of a single GPhys object. " +
    "Use open_multi_gturl to treat a gturl of multiple objects."
  end
  gp = GPhys::IO.open(file,var)
  gp = gp[slice] if slice
  gp = gp.cut(cut_slice) if cut_slice
  gp = gp[thinning] if thinning
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
