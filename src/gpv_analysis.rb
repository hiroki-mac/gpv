#--
# data analysis methods for gpv.
#++



class GPV

  # Lowpass filter using FFT.
    def lowpass(gp)
    lp_dim = @OPT_lowpass.split(",")[0].to_i
    lp_wn = (@OPT_lowpass.split(",")[1].to_i)+1
    sp = gp.fft(false,lp_dim)
    mask = sp.val.real.fill(1.0)
    mask[lp_wn..-lp_wn,false] = 0.0 if lp_dim == 0
    mask[true, lp_wn..-lp_wn,false] = 0.0 if lp_dim == 1
    mask[true, true, lp_wn..-lp_wn,false] = 0.0 if lp_dim == 2
    mask[true, true, true, lp_wn..-lp_wn,false] = 0.0 if lp_dim == 3
    gp = (sp*mask).fft(true,lp_dim).real
    # rad => deg for lon
    # gp.axis(0).set_pos(VArray.new(gp.coordinate(0).val*R2D, {"units"=>"deg"}, "lon"))
    return gp
  end

  # Highpass filter using FFT.
  def highpass(gp)
    hp_dim = @OPT_highpass.split(",")[0].to_i
    hp_wn = @OPT_highpass.split(",")[1].to_i
    sp = gp.fft(false,hp_dim)
    mask = sp.val.real.fill(0.0)
    mask[hp_wn..-hp_wn,false] = 1.0 if hp_dim == 0
    mask[true, hp_wn..-hp_wn,false] = 1.0 if hp_dim == 1
    mask[true, true, hp_wn..-hp_wn,false] = 1.0 if hp_dim == 2
    mask[true, true, true, hp_wn..-hp_wn,false] = 1.0 if hp_dim == 3
    gp = (sp*mask).fft(true,hp_dim).real
    return gp
  end

  # Lowpass filter using Spherical harmonics transformation.
    def lowpass_sht(gp)
    lp_wn = @OPT_lowpass.to_i + 1
    sp_mn = SphericalHarmonics.sh_trans(gp) #正変換
    mask = sp_mn.val.real.fill(1.0)
    mask[lp_wn..-lp_wn,true,false] = 0.0
    mask[true,lp_wn..-1,false] = 0.0
    gp = SphericalHarmonics.sh_invtrans(sp_mn*mask)
    # rad => deg for lon
    # gp.axis(0).set_pos(VArray.new(gp.coordinate(0).val*R2D, {"units"=>"deg"}, "lon"))
    return gp
  end
  # Highpass filter using Spherical harmonics transformation.
  def highpass_sht(gp)
    hp_wn = @OPT_highpass.to_i
    sp_mn = SphericalHarmonics.sh_trans(gp) #正変換
    mask = sp_mn.val.real.fill(0.0)
    mask[hp_wn..-hp_wn,false] = 1.0
    mask[true, hp_wn..-1,false] = 1.0
    gp = SphericalHarmonics.sh_invtrans(sp_mn*mask)
    return gp
  end
  # Bandpass filter using Spherical harmonics transformation.
  def bandpass_sht(gp)
    lp_wn = @OPT_lowpass.to_i + 1
    hp_wn = @OPT_highpass.to_i
    sp_mn = SphericalHarmonics.sh_trans(gp) #正変換
    mask = sp_mn.val.real.fill(0.0)
    mask[hp_wn..-hp_wn,false] = 1.0
    mask[true, hp_wn..-1,false] = 1.0
    mask[lp_wn..-lp_wn,true,false] = 0.0
    mask[true,lp_wn..-1,false] = 0.0
    gp = SphericalHarmonics.sh_invtrans(sp_mn*mask)
    return gp
  end


  # Composite mean.
  def composite(gp)
    comp_dim, comp_span, comp_skip = (@OPT_composite).split(/\s*,\s*/)

    if comp_dim.to_i.to_s  == comp_dim then
      comp_dim = comp_dim.to_i
    end

    comp_span = comp_span.to_i
    comp_skip = comp_skip.to_i

    if (comp_skip != 0)
      tmp = 0.0
      (gp.axis(comp_dim).length/(comp_span+comp_skip)).times{|i|
        sp = i*(comp_span+comp_skip)
        tmp = tmp + gp[{comp_dim=>sp..(sp+comp_span-1)}].mean(comp_dim)
      }
      gp = tmp
    else
      tmp = 0.0
      (gp.axis(comp_dim).length/(comp_span)).times{|i|
        sp = i*comp_span
        tmp = tmp + gp[{comp_dim=>sp..(sp+comp_span-1)}]
      }
      gp.axis(comp_dim).length/(comp_span)
      gp = tmp / (gp.axis(comp_dim).length/(comp_span))
    end
    return gp
  end

  # Powerspectra using FFT.
  def powerspectra(gp)
    ps_dim = @OPT_power_spectra.split(",")[0].to_i
    flag_wl = @OPT_power_spectra.split(",")[1]
    axunit = gp.axis(ps_dim).to_gphys.units.to_s
    axname = gp.axnames[ps_dim]
    gp = (gp.fft(false,ps_dim).abs**2).rawspect2powerspect(ps_dim)
    wn = gp.shape[ps_dim]
    ax = gp.axis(ps_dim).to_gphys.val
    gp = gp.cut(gp.axnames[ps_dim]=>ax[1]..ax[wn/2]) # 波数 0 は出さないようにしている
    if (flag_wl)
      axunit = axunit.split("since")[0] if (axunit.include?("since"))
      axat = Attribute.new
      axat["units"]=axunit
      if (axname.include?("lon") or axname.include?("lat") or axname=="x" or axname=="y")
        gp.axis(ps_dim).set_pos(VArray.new(2.0*PI/ax[1..wn/2]*180/PI, axat, "wavelength - "+axname))
      elsif (axname.include?("time") or axname=="t")
        gp.axis(ps_dim).set_pos(VArray.new(2.0*PI/ax[1..wn/2], axat, "period - "+axname))
      end
    end
    print "Operation: power spectra along #{axname} is calculated.\n"
    return gp
  end

  # Convert sigma coordinate ot pressure coordinate.
  def sig2p(gp,ps=@PS)
    dim = (GAnalysis::SigmaCoord::find_sigma_d(gp) || 2) # if cannot find sigma axis, assume 3rd axis.
    sig = gp.axis(dim).to_gphys
    if (sig.length*ps.length*8>2.2E9 && NArrayType == "standard") then # 2GB超のオブジェクト生成を避ける
      @OPT_parallel = true; @OPT_sig2p = true
      return gp
    else
      @OPT_sig2p = false unless @OPT_nc4a
      press = GAnalysis::SigmaCoord::sig_ps2p(ps, sig, dim)
      press.units=("Pa"); press.set_att("positive","down")
      gp.set_assoc_coords([press])
  #    gp = gp.assoc_coord_gphys("p")
      return gp
    end
  end

  # Convert sigma coordinate to height coordinate.
  def sig2z(gp, gt=nil)
    print "CAUTION: sig2z needs data from lowest layer for calculation.\n"
    @T = gt if gt
    dim = (GAnalysis::SigmaCoord::find_sigma_d(gp) || 2) # if cannot find sigma axis, assume 3rd axis.
    sig = gp.axis(dim).to_gphys.val
    if (@T.length*8>2.2E9 && NArrayType == "standard") then # 2GB超のオブジェクト生成を避ける
      @OPT_parallel = true; @OPT_sig2z = true
      return gp
    else
      @OPT_sig2z = false
      t_na = @T.val
      z_na = t_na*0.0 # z の NArray を用意
      fact = (-GAnalysis::Met::R/GAnalysis::Met.g).to_f
      if (dim == 0) # 鉛直軸 が 1次元目の場合
        z_na[0,false] = 2.0*fact*(sig[0]-1.0)/(sig[0]+1.0)*t_na[0,false] # z_1
        for k in 1..(sig.length-1) # z_2...K
          z_na[k,false] = z_na[true,true,k-1,false] + fact*(sig[k]-sig[k-1])/(sig[k]+sig[k-1]) \
                          *(t_na[k,false]+t_na[k-1,false])
        end
      elsif (dim == 1) # 鉛直軸 が 2次元目の場合
        z_na[true,0,false] = 2.0*fact*(sig[0]-1.0)/(sig[0]+1.0)*t_na[true,0,false] # z_1
        for k in 1..(sig.length-1) # z_2...K
          z_na[true,k,false] = z_na[true,k-1,false] + fact*(sig[k]-sig[k-1])/(sig[k]+sig[k-1]) \
                                      *(t_na[true,k,false]+t_na[true,k-1,false])
        end

      elsif (dim == 2) # 鉛直軸 が 3次元目の場合
        z_na[true,true,0,false] = 2.0*fact*(sig[0]-1.0)/(sig[0]+1.0)*t_na[true,true,0,false] # z_1
        for k in 1..(sig.length-1) # z_2...K
          z_na[true,true,k,false] = z_na[true,true,k-1,false] \
                                    + fact*(sig[k]-sig[k-1])/(sig[k]+sig[k-1]) \
                                      *(t_na[true,true,k,false]+t_na[true,true,k-1,false])
        end
      else
        raise "vertical axis must 1st, 2nd, or 3rd axis\n"
      end


      z = GPhys.new(gp.grid.copy,VArray.new(z_na,{"long_name"=>"height","units"=>"m"},"Z"))
      gp.set_assoc_coords([z])
      return gp
    end
  end

  # Adding up along the given dimension.
  def addingup(gp, dim)
    gpaxes = gp.axnames
    if (dim.class == String) then
      d = nil; gpaxes.length.times{|i| d = i if gpaxes[i] == dim}; dimname = dim
      raise "Error::dimension of #{dim} for addingup was not found." if d == nil
    else d = dim; dimname = gpaxes[d] end
    newdata = gp.val*0.0
    len = gp.shape[d]
    gpval = gp.val
    sum = gpval.sum(d)*0.0
    if (d == 0) then
      len.times{|i| sum = sum + gpval[i,false];  newdata[i,false] = sum }
    elsif (d == 1) then
      len.times{|i| sum = sum + gpval[true,i,false];  newdata[true,i,false] = sum }
    elsif (d == 2) then
      len.times{|i| sum = sum + gpval[true,true,i,false];  newdata[true,true,i,false] = sum }
    elsif (d == 3) then
      len.times{|i| sum = sum + gpval[true,true,true,i,false];  newdata[true,true,true,i,false] = sum }
    elsif (d == 4) then
      len.times{|i| sum = sum + gpval[true,true,true,true,i,false];  newdata[true,true,true,true,i,false] = sum }
    elsif (d == 5) then
      len.times{|i| sum = sum + gpval[true,true,true,true,true,i,false];  newdata[true,true,true,true,true,i,false] = sum }
    else raise "rank greater than 6 is not supported"
    end
    varray = VArray.new(newdata,nil,gp.name+" addingup along "+dimname)
    return GPhys.new(gp.grid, varray)
  end

  # z座標でのPVを計算する。引数はz座標上の東西風、南北風、温位、密度。
  # 惑星半径と自転角速度は適切に設定する必要がある。
  # 鉛直軸は第３軸と想定。
  def pv_on_z(u,v,theta,rho)
    z = 2
    uz = u.threepoint_O2nd_deriv(z); vz = v.threepoint_O2nd_deriv(z)
    thetaz = theta.threepoint_O2nd_deriv(z)
    if (@OPT_sht) # 球面調和関数変換を使う
      absvor = SphericalHarmonics::sh_vorticity(u,v,@Radius) + SphericalHarmonics::sinlat*2.0*@Omega.to_f
      thetax = SphericalHarmonics::sh_dlon(theta)/(SphericalHarmonics::coslat*@Radius.to_f)
      thetay = SphericalHarmonics::sh_dlat(theta)*(SphericalHarmonics::coslat)/@Radius.to_f # (注) sh_dlat は d/dmu のこと。
    else
      absvor = GAnalysis::Planet.absvor_s(u,v) # 絶対渦度
      thetax, thetay = GAnalysis::Planet.grad_s(theta)
    end
    return (absvor*thetaz - vz*thetax + uz*thetay )/rho
  #  return (absvor*thetaz + uz*thetay)/rho
  end

  # z座標での質量流線関数を求める。
  # \Psi = cos(\phi) \sum (\bar{ \rho * v })dz を計算して、質量流線関数を求める。
  # -1/cos(\phi) * ∂ψ/∂z = \bar{\rho * v}, 1/(a*cos(\phi)) * ∂ψ/∂φ = \bar{\rho w} で定義されている。
  # 半グリッドずれた(cell boundary)で定義されたPsiを返す
  def mstrmfunc_on_z(v, rho)
    print "CAUTION: mstrmfunc_on_z needs data from lowest layer for calculation.\n"
    vrm = (v*rho).mean(0)
    vrm = vrm*(vrm.axis(0).to_gphys*D2R).cos # \bar{rho*v}cos(\phi)
    vrm_na = vrm.val.to_na
    bnd_grid = make_bnd_grid(vrm.grid, 1)
    psi = NArray.float(*bnd_grid.shape) # psi の NArray を用意
    z = bnd_grid.coord(1).val       # sigma の NArray
    km = z.length-1

    psi[true,0,false] = 0.0 
    for k in 0..(km-1) # 積み上げ
      psi[true,k+1,false] = psi[true,k,false] + (vrm_na[true,k,false])*(z[k+1]-z[k])
    end

    return GPhys.new(bnd_grid,VArray.new(-psi,{"long_name"=>"mass streamfunction","units"=>"kg.m-1.s-1"},"msf"))
  end

  # σ座標での質量流線関数を求める。
  # \Psi = cos(\phi) \sum (\bar{ v })dsigma を計算して、質量流線関数を求める。
  # -1/cos(\phi) * ∂ψ/∂sigma = \bar{v}, 1/(a*cos(\phi)) * ∂ψ/∂φ = \bar{sig_dot} で定義されている。
  # 半グリッドずれた(cell boundary)で定義されたPsiを返す
  def mstrmfunc_on_sigma(v)
    print "CAUTION: mstrmfunc_on_sigma needs data from lowest layer for calculation.\n"
    v = v.mean(0)
    v = (v*(v.axis(0).to_gphys*D2R).cos)  # \bar{v}cos(\phi)
    v_na = v.val.to_na
    bnd_grid = make_bnd_grid(v.grid, 1)
    psi = NArray.float(*bnd_grid.shape) # psi の NArray を用意
    sigma = bnd_grid.coord(1).val       # sigma の NArray
    km = sigma.length-1
    
    psi[true,0,false] = 0.0 # 最下層
    for k in 0..(km-1) # 積み上げ
      psi[true,k+1,false] = psi[true,k,false] + (v_na[true,k,false])*(sigma[k+1]-sigma[k])
    end

    return GPhys.new(bnd_grid,VArray.new(psi,{"long_name"=>"mass streamfunction","units"=>"m.s-1"},"msf"))
  end


  # 時空間スペクトル（とりあえずGCM出力データのみに対応）
  # 入力は（lon, lat, time）の3次元データ。ただし lat は赤道域（±15 deg）とする。
  def wnf_analysis(gp)
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, true)
    # detrend
    gp = gp.detrend(-1)
    # calculate symmetric and anti-symmetric components
    gp_sym  = (gp + gp[true,-1..0,true])*0.5
    gp_asym = (gp - gp[true,-1..0,true])*0.5
    # calculate raw spectra
    GPhys::fft_ignore_missing(true)
    nt = gp.shape[-1]-1 # 時間方向のデータ数
    t = 0
    unit = gp.coord(-1).units.to_s
    unit, since = unit.split("since") if unit.include?("since")
    if (unit.downcase.include?("min")) then 
      day_in_unit = 1440
    elsif (unit.downcase.include?("sec")) then 
      day_in_unit = 86400
    elsif (unit.downcase.include?("day")) then
      day_in_unit = 1
    else
      raise "unit: #{gp.coord(-1).units.to_s} is unsupported. "
    end
    t_val = gp.coord(-1).val
#    span = ((t_val[-1]-t_val[0])/day_in_unit).to_i # 1解析期間 (days)
    span = 1 # 1解析期間 (days)
    ndpd = day_in_unit/(t_val[1]-t_val[0]) # 1 day あたりのデータ数
    while (((t+1)*span*ndpd)-1 < nt)
      trange = (t*span*ndpd)..[((t+1)*span*ndpd - 1), nt].min
      gp_sym_sp  = gp_sym[true,true,trange].detrend(-1).cos_taper(-1).fft(false,0,2)
      gp_asym_sp = gp_asym[true,true,trange].detrend(-1).cos_taper(-1).fft(false,0,2)
      gp_sym_pw, gp_asym_pw = [gp_sym_sp, gp_asym_sp].map{|g|
        # change axes
        z0=g.axis(0).pos
        g.axis(0).set_pos(-z0); g.coord(0).name = 'k'; g.coord(0).set_att('long_name', 'k-wavenumber')
        g.coord(0).del_att('standard_name'); g.coord(0).del_att('axis'); g.coord(0).del_att('gt_calc_weight')
        z1=g.axis(-1).pos.convert_units("day-1")/(2.0*Math::PI)
        g.axis(-1).set_pos(z1); g.coord(-1).name = 'freq'; g.coord(-1).set_att('long_name', 'frequency')
        g.coord(-1).del_att('standard_name'); g.coord(-1).del_att('axis'); g.coord(-1).del_att('delta_t')
        # calculate powerspectra
        power=(g.abs**2).mean(1).rawspect2powerspect(0,-1).spect_zero_centering(0).spect_one_sided(-1)* GPhys::COS_TAPER_SP_FACTOR
      }
      if (t == 0) then
        gp_sym_pw_sum = gp_sym_pw.copy; gp_asym_pw_sum = gp_asym_pw.copy
      else
        gp_sym_pw_sum = gp_sym_pw_sum + gp_sym_pw; gp_asym_pw_sum = gp_asym_pw_sum + gp_asym_pw
      end
      t = t + 1
    end
    gp_sym_pw = gp_sym_pw_sum/t ; gp_asym_pw = gp_asym_pw_sum/t

    # estimate background spectra by smoothing
    tmp = (gp_sym_pw + gp_asym_pw).val
    10.times{ # smoothing in frequency
      tmp0 = (2.0*tmp[true,0] + tmp[true,1])/3.0; tmp1 = (tmp[true,-2] + 2.0*tmp[true,-1])/3.0
      tmp[true,1..-2] = (tmp[true,0..-3] + 2.0*tmp[true,1..-2] + tmp[true,2..-1])/4.0
      tmp[true,0] = tmp0; tmp[true,-1] = tmp1
    }
    40.times{ # smoothing in k
      tmp0 = (2*tmp[0,true] + tmp[1,true])/3.0; tmp1 = (tmp[-2,true] + 2*tmp[-1,true])/3.0
      tmp[1..-2,true] = (tmp[0..-3,true] + 2*tmp[1..-2,true] + tmp[2..-1,true])/4.0
      tmp[0,true] = tmp0; tmp[-1,true] = tmp1
    }
    gp_bg_pw = gp_sym_pw.copy; gp_bg_pw.replace_val(tmp)

    return gp_sym_pw.set_att('long_name', 'sym-spectra of '+gp_sym_pw.name ) if @OPT_wnf_analysis == "sym"
    return gp_asym_pw.set_att('long_name', 'asym-spectra of '+gp_asym_pw.name ) if @OPT_wnf_analysis == "asym"
    return gp_bg_pw.set_att('long_name', 'log10   bg-spectra of '+gp_bg_pw.name ).log10 if @OPT_wnf_analysis == "bg"
    return (gp_sym_pw/gp_bg_pw).set_att('long_name', 'sym/bg-spectra of '+gp_sym_pw.name ) if @OPT_wnf_analysis == "sym/bg"
    return (gp_asym_pw/gp_bg_pw).set_att('long_name', 'asym/bg-spectra of '+gp_asym_pw.name ) if @OPT_wnf_analysis == "asym/bg"
    return ((gp_sym_pw/gp_bg_pw).log10).set_att('long_name', 'log10(sym/bg-spectra) of '+gp_sym_pw.name ) if @OPT_wnf_analysis == "log_sym/bg"
    return ((gp_asym_pw/gp_bg_pw).log10).set_att('long_name', 'log10(asym/bg-spectra) of '+gp_asym_pw.name ) if @OPT_wnf_analysis == "log_asym/bg"
    raise "argument for wnf_analysis must be \"sym, asym, bg, sym/bg, asym/bg, log_sym/bg, or log_asym/bg \" \n"
  end

  # EP-flux
  def epflux(gp_u,gp_v,gp_omega,gp_t)
    print "Caution: EP-flux calculation requires *** p-coordinate *** for vertical axis. \n"
    require "numru/gphys/ep_flux"
    GPhys::EP_Flux.radius=GAnalysis::Planet.radius
    GPhys::EP_Flux.rot_period=2.0*PI/GAnalysis::Planet.omega
    GPhys::EP_Flux.g_forces=GAnalysis::Met.g
    GPhys::EP_Flux.p00=GAnalysis::Met::P00
    GPhys::EP_Flux.gas_const=GAnalysis::Met::R
    GPhys::EP_Flux.cp=GAnalysis::Met::Cp
    GPhys::EP_Flux.scale_height=UNumeric.new(16000.0,  "m") # set scale_height if needed

    ofile = NetCDF.create("epflux.nc")

    print "         Scale height of #{GPhys::EP_Flux.scale_height} is used.\n"

    GPhys::IO.each_along_dims_write([gp_u, gp_v, gp_omega, gp_t], ofile, -1){|u, v, omega, t|
      epflx_y, epflx_z, v_rmean, w_rmean, gp_lat, gp_z, u_mean, theta_mean,
            uv_dash, vt_dash, uw_dash, dtheta_dz = ary = GPhys::EP_Flux::ep_full_sphere(u, v, omega, t, false) # false if potential temperature

      epflx_div = GPhys::EP_Flux::div_sphere(epflx_y, epflx_z)
      strm_rmean = GPhys::EP_Flux::strm_rmean(v_rmean)
#      binding.pry
      [ epflx_y, epflx_z, v_rmean, w_rmean, strm_rmean, epflx_div, u_mean, theta_mean, uv_dash, vt_dash, uw_dash, dtheta_dz ]
    }

    ofile.close
    abort
  end

  # Generalized EP flux derived from Vallis (2006) texts. See Kashimura(Yamamoto)'s Master Thesis'
  def g_epflux(gp_u,gp_v,gp_omega,gp_t)
    print "Caution: EP-flux calculation requires *** p-coordinate *** for vertical axis. \n"
    require "numru/gphys/ep_flux"
    GPhys::EP_Flux.radius=GAnalysis::Planet.radius
    GPhys::EP_Flux.rot_period=2.0*PI/GAnalysis::Planet.omega
    GPhys::EP_Flux.g_forces=GAnalysis::Met.g
    GPhys::EP_Flux.p00=GAnalysis::Met::P00
    GPhys::EP_Flux.gas_const=GAnalysis::Met::R
    GPhys::EP_Flux.cp=GAnalysis::Met::Cp
    # t_s =
    # GPhys::EP_Flux.scale_height=UNumeric.new((GPhys::EP_Flux.gas_const*t_s/GPhys::EP_Flux.g_forces),  "m") # set scale_height if needed
    print "         Scale height of #{GPhys::EP_Flux.scale_height} is used.\n"
    print "         [Note] You should multiply by a*rho_s = a*p_s/(gH) = #{@Radius*GPhys::EP_Flux.p00/(GPhys::EP_Flux.g_forces*GPhys::EP_Flux.scale_height)} to convert to Andrews et al.'s (1987) form.\n"
    scale_height = GPhys::EP_Flux.scale_height
    ofile = NetCDF.create("g_epflux.nc")
    flag_temp_or_theta = true # false if potential temperature

    GPhys::IO.each_along_dims_write([gp_u, gp_v, gp_omega, gp_t], ofile, -1){|u, v, omega, t|
      xyzdims = [0,1,2]
      ax_lon = u.axis(xyzdims[0]) # Axis of longitude
      ax_lat = u.axis(xyzdims[1]) # Axis of latitude
      ax_z =   u.axis(xyzdims[2]) # Axis of vertical
      lon_nm, lat_nm, z_nm = ax_lon.pos.name, ax_lat.pos.name, ax_z.pos.name
      gp_lon, gp_lat, gp_z = GPhys::EP_Flux.make_gphys(ax_lon, ax_lat, ax_z)

      ## convert axes
      gp_z   = GPhys::EP_Flux.to_z_if_pressure(gp_z)     # P => z=-H*log(P/P00) (units-based)
      gp_lon = GPhys::EP_Flux.to_rad_if_deg(gp_lon)    # deg => rad (unit convesion)
      gp_lat = GPhys::EP_Flux.to_rad_if_deg(gp_lat)    # deg => rad (unit convesion)
      w      = GPhys::EP_Flux.to_w_if_omega(omega, gp_z)  # dP/dt => dz/dt (units-based)
      t      = GPhys::EP_Flux.to_theta_if_temperature(t, gp_z, flag_temp_or_theta)
                    # temperature => potential temperature (if flag is true)
  #    gp_lon.units=("radian"); gp_lat.units=("radian")

      ## replace grid (without duplicating data)
      grid = u.grid_copy
      old_grid = u.grid_copy                 # saved to use in outputs
      grid.axis(lon_nm).pos = gp_lon.data       # in radian
      grid.axis(lat_nm).pos = gp_lat.data       # in radian
      grid.axis(z_nm).pos = gp_z.data           # log-p height
      u = GPhys.new(grid, u.data)
      v = GPhys.new(grid, v.data)
      w = GPhys.new(grid, w.data)
      t = GPhys.new(grid, t.data)
      ## get each term
      #  needed in F_y and F_z
      uv_dash, vt_dash, uw_dash = GPhys::EP_Flux.eddy_products(u, v, w, t, lon_nm)
      wt_dash = (t.eddy(lon_nm)*w.eddy(lon_nm)).mean(lon_nm)
      wt_dash.data.set_att("long_name","W'T'"); wt_dash.data.rename!("wt_dash")
      theta_mean = t.mean(lon_nm)
      dtheta_dz = GPhys::EP_Flux.deriv(theta_mean, z_nm)
      # dtheta_dz * dtheta_dz.gt(0) # to eliminate negative value
  #    binding.pry

      cos_lat = gp_lat.cos
      a_cos_lat = @Radius * cos_lat
      a_cos_lat.data.rename!('a_cos_lat')
      a_cos_lat.data.set_att('long_name', 'radius * cos_lat')
      GPhys::EP_Flux.remove_0_at_poles(a_cos_lat)
      #  needed in F_y only
      u_mean = u.mean(lon_nm)
      du_dz  = GPhys::EP_Flux.deriv(u_mean, z_nm)
      #  needed in F_z only
      f_cor = 2 * @Omega * gp_lat.sin
      f_cor.data.rename!('f_cor')
      f_cor.data.set_att('long_name', 'Coriolis parameter')
      dtheta_dphi = GPhys::EP_Flux.deriv( theta_mean, lat_nm)
      ducos_dphi = GPhys::EP_Flux.deriv( u_mean * cos_lat, lat_nm)
      avort = (-ducos_dphi/a_cos_lat) + f_cor        # -- absolute vorticity
      avort.data.units = "s-1"
      avort.data.rename!('avort')
      avort.data.set_att('long_name', 'zonal mean absolute vorticity')

      # calculate psi (see Eq. C.45 in Kashimura(Yamamoto)'s Master Thesis)
      # [note] psi can be approximated to vt_dash/dtheta_dz in typical Earth cases.
      psi = (vt_dash*dtheta_dz - wt_dash*dtheta_dphi/@Radius)/((dtheta_dphi/@Radius)**2 + dtheta_dz**2)

      ## F_y, F_z
      sigma = (-gp_z/scale_height).exp
      epflx_y = ( - uv_dash + du_dz*psi ) * cos_lat * sigma
      epflx_z = ( - uw_dash + avort*psi ) * cos_lat * sigma
      epflx_y.data.name = "epflx_y"; epflx_z.data.name = "epflx_z"
      epflx_y.data.set_att("long_name", "Generalized EP flux y component")
      epflx_z.data.set_att("long_name", "Generalized EP flux z component")

      ## v_rmean, w_rmean
      z_nm = gp_z.data.name    # change z_nm from pressure to z
      v_mean = v.mean(lon_nm); w_mean = w.mean(lon_nm)
      v_rmean = ( v_mean - GPhys::EP_Flux.deriv( (psi*sigma), z_nm )/sigma )
      w_rmean = ( w_mean + GPhys::EP_Flux.deriv( (psi*cos_lat), lat_nm )/a_cos_lat )
      v_rmean.data.name = "v_rmean"; w_rmean.data.name = "w_rmean"
      v_rmean.data.set_att("long_name", "residual zonal mean V")
      w_rmean.data.set_att("long_name", "residual zonal mean W")

      ## convert with past grid
    	gp_ary = [] # grid convertes gphyss into
    	grid_xmean = old_grid.delete_axes(lon_nm)
    	[epflx_y, epflx_z, v_rmean, w_rmean, gp_lat, gp_z, u_mean, theta_mean, uv_dash, vt_dash, wt_dash, uw_dash, dtheta_dz, dtheta_dphi/@Radius].each {|g|
    	  if grid_xmean.shape.size != g.shape.size
    	    gp_ary << g
    	  else
    	    gp_ary << GPhys.new(grid_xmean, g.data) #back to the original grid
    	  end
      }

      gp_ary[5] = GPhys::EP_Flux::div_sphere(gp_ary[0], gp_ary[1])
      gp_ary[4] = GPhys::EP_Flux::strm_rmean(gp_ary[2])

  #    [ epflx_y, epflx_z, v_rmean, w_rmean, strm_rmean, epflx_div, u_mean, theta_mean, uv_dash, vt_dash, uw_dash, dtheta_dz ]
  #   binding.pry
      gp_ary
    }

    ofile.close
    abort
  end

  # Divergence on meridional plane
  def div_on_meridional_plane(v, w)
    require "numru/gphys/ep_flux"
    GPhys::EP_Flux.radius=GAnalysis::Planet.radius
    GPhys::EP_Flux.rot_period=2.0*PI/GAnalysis::Planet.omega
    GPhys::EP_Flux.g_forces=GAnalysis::Met.g
    GPhys::EP_Flux.p00=GAnalysis::Met::P00
    GPhys::EP_Flux.gas_const=GAnalysis::Met::R
    GPhys::EP_Flux.cp=GAnalysis::Met::Cp

    div = -GPhys::EP_Flux::div_sphere(v, w)
    div.rename("div_#{v.name}_#{w.name}")
    div.set_att("long_name", "divergence of (#{v.name}, #{w.name})")
    return div
  end

  # Gradient along lat direction on sphere.
  def grad_sy(gp,lat=0)
    lam, phi, lond, latd  = GAnalysis::Planet.get_lambda_phi(gp,false)
    cos_phi = phi.cos
    ys = gp.cderiv(latd,GAnalysis::Planet.latbc(phi),phi) / GAnalysis::Planet.radius
    return ys
  end

  def solar_zenith_angle(gp) # solar_zenith_angle at local noon
    gt = gp.axis("time").to_gphys; vt = gt.val;  day_array = []
    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, false)
    glat = gp.axis(lat_dim).to_gphys; vlat = glat.val*PI/180.0

    ndays = NArray.float(vlat.length, vt.length)
    hpara = NArrayMiss.float(vlat.length, vt.length).fill!(1.0)
    cos_theta0 = NArray.float(vlat.length, vt.length)

    vt.length.times{|i|
      un = UNumeric.new(vt[i],gt.units.to_s)
      ndays[true,i] = un.to_datetime.yday + un.to_datetime.hour/24.0
    }
    gamma = (ndays - 1.0)*2.0*PI/365.0
    gamma[0,5..-1]
    # Solar declination (太陽赤緯)
    delta = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2.0*gamma) + 0.000907*sin(2.0*gamma) - 0.002697*cos(3.0*gamma) + 0.001480*sin(3.0*gamma)

    # 日積算量の計算には必要↓
    cosH = -tan(delta)*tan(vlat)
    # 極夜・白夜に対応する
    mask_polar_night = cosH.ge(1.0); mask_polar_day = cosH.le(-1.0)
    cosH = cosH*(mask_polar_night.eq(0.0)*mask_polar_day.eq(0.0)) +
           mask_polar_night.to_f*1.0 + mask_polar_day.to_f*(-1.0)
    hpara = acos(cosH)
    cos_theta0 = (hpara*sin(vlat)*sin(delta) + cos(vlat)*cos(delta)*sin(hpara))/PI
    theta0 = acos(cos_theta0)*180.0/PI

    # 南中時（時角=0)の天頂角で、その日の天頂角を代表させる。
    cos_theta0 = (sin(vlat)*sin(delta)+cos(vlat)*cos(delta))
    theta0 = acos(cos_theta0)*180.0/PI
    theta0 = NArrayMiss.to_nam(theta0).set_mask(mask_polar_night.eq(0))

    gtheta0 = gp*0.0 + theta0
    gtheta0.name=("solar_zenith_angle")
    gtheta0.set_att("long_name", "solar zenith angle at local noon")
    gtheta0.set_att("units","deg")

    return gtheta0
  end

  def solar_zenith_angle_planet(gp) # 公転軌道は円軌道を仮定
    # unit_time = min の場合
    unit_time_in_day  = 1440.0
    days_in_year      = 365*2
    unit_time_in_year = unit_time_in_day*days_in_year

    # 赤道傾斜角
    axial_tilt = 25.9*D2R
    # 春分の日付 (単位付き)
    equinox = UNumeric.new(0.0, "minutes since 0000-01-01 00:00:00")

    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, false)
    lat = gp.coord(lat_dim).val*D2R
    lon = gp.coord(lon_dim).val*D2R
    time = gp.coord("time").convert_units(equinox.units).val

    # 赤経(?) [time]
    ra = time.mod(unit_time_in_year - equinox)/unit_time_in_year*2.0*PI
    # lon=0での時角 [time]
    ha0 = time.mod(unit_time_in_day)/unit_time_in_day*2.0*PI-PI
    # sin(赤緯) [time]
    sin_delta = ra.sin * sin(axial_tilt)
    # cos(太陽天頂角)
    cos_Z = NArray.float(lon.length, lat.length, time.length)

    puts "===NOTE==============================="
    puts "unit_time_in_day = #{unit_time_in_day}"
    puts "days_in_year = #{days_in_year}"
    puts "axial_tilt = #{axial_tilt*R2D}"
    puts "equinox = #{equinox}"
    puts "===NOTE==============================="


    time.length.times{|t|
      ha = lon + ha0[t]
      lat.length.times{|j|
        cos_Z[true,j,t] = ha.cos * cos(lat[j]) * sqrt(1.0 - sin_delta[t]*sin_delta[t]) + sin(lat[j])*sin_delta[t]
      }
    }
    # 日積算量の計算には必要↓
    # cosH = -tan(delta)*tan(vlat)
    # 極夜・白夜に対応する
    # mask_polar_night = cosH.ge(1.0); mask_polar_day = cosH.le(-1.0)
    # cosH = cosH*(mask_polar_night.eq(0.0)*mask_polar_day.eq(0.0)) +
    #        mask_polar_night.to_f*1.0 + mask_polar_day.to_f*(-1.0)
    # hpara = acos(cosH)
    # cos_theta0 = (hpara*sin(vlat)*sin(delta) + cos(vlat)*cos(delta)*sin(hpara))/PI
    # theta0 = acos(cos_theta0)*180.0/PI

    # 南中時（時角=0)の天頂角で、その日の天頂角を代表させる。
    # cos_theta0 = (sin(vlat)*sin(delta)+cos(vlat)*cos(delta))
    # theta0 = acos(cos_theta0)*180.0/PI
    # theta0 = NArrayMiss.to_nam(theta0).set_mask(mask_polar_night.eq(0))

    gtheta0 = gp*0.0 + cos_Z.acos*R2D
    gtheta0.name=("solar_zenith_angle")
    gtheta0.set_att("long_name", "solar zenith angle")
    gtheta0.set_att("units","deg")

    return gtheta0
  end


  def reflectance(gp) # 浅野正二「大気放射学の基礎」式 6.23
    # 入力 gp は 太陽天頂角のGPhysオブジェクト
    mask = gp.val.get_mask
    mu0 = cos(gp.val.to_na*PI/180.0)
    # p mu0
    # mu0 = cos(80.0*PI/180.0)
    omega = 0.99999; g = 0.5; tau = 0.05
    # Eddington近似
    gamma1 = 0.25*(7-omega*(4+3*g)); gamma2 = -0.25*(1-omega*(4-3*g))
    gamma3 = 0.25*(2-3*g*mu0); gamma4 = 1.0 - gamma3
    alpha1 = gamma4*gamma1 + gamma3*gamma2; alpha2 = gamma3*gamma1 + gamma4*gamma2
    kappa = (gamma1*gamma1 - gamma2*gamma2)**0.5

    r = omega/(1.0 - kappa*kappa*mu0*mu0)
    r = r/((kappa + gamma1)*exp(kappa*tau) + (kappa - gamma1)*exp(-kappa*tau))
    r = r*( (1.0 - kappa*mu0)*(alpha2 + kappa*gamma3)*exp(kappa*tau) - (1.0+kappa*mu0)*(alpha2 - kappa*gamma3)*exp(-kappa*tau) - 2.0*kappa*(gamma3 - alpha2*mu0)*exp(-tau/mu0))
    return gp*0 + r
  end

  # 剛体回転（角速度：angvel）分、東西にシフトさせる。最後の軸は時間。
  # 経度軸の単位は deg, 等間隔。
  # シフトは最近傍のグリッドデータを使う（補間・内挿はしない）
  def zonal_shift_on_sphere(gp, u)
    if (u.include?("@")) then
      u = open_gturl(u).mean(0).val
      angvel = u/@Radius
    else
      u = u.to_f
      angvel = u/@Radius.to_f
    end


    lon_dim, lat_dim = GAnalysis::Planet.find_lon_lat_dims(gp, true)
    time_ax = gp.axis(-1).to_gphys
    lon = gp.axis(lon_dim).to_gphys.val
    del_lon = lon[1] - lon[0]
    lat = gp.axis(lat_dim).to_gphys.val
    new_val = gp.val

    # 時間軸の単位の秒数
    if (time_ax.units.to_s.include?("day")) then
      unit_time = 24.0*60.0*60.0
    elsif (time_ax.units.to_s.include?("hour")) then
      unit_time = 60.0*60.0
    elsif (time_ax.units.to_s.include?("min")) then
      unit_time = 60.0
    elsif (time_ax.units.to_s.include?("sec")) then
      unit_time = 1.0
    else
      unit_time = 1.0
    end

    time_ax.length.times{|t|
      if (angvel.class == Float) then
        delta_deg = (time_ax.val[t] - time_ax.val[0])*unit_time*angvel*R2D % 360
        delta_i = (delta_deg/del_lon).round % lon.length
        if (gp.rank == 4) then
          # gp に鉛直軸がある（x,y,z,t）場合
          new_val[0..-(delta_i+1),true,true,t] = gp[(delta_i)..-1,true,true,t].val
          new_val[-(delta_i)..-1,true,true,t]  = gp[0..(delta_i-1),true,true,t].val
        else
          # gp に鉛直軸がない（x,y,t）場合
          new_val[0..-(delta_i+1),true,t] = gp[(delta_i)..-1,true,t].val
          new_val[-(delta_i)..-1,true,t]  = gp[0..(delta_i-1),true,t].val
        end
      elsif (angvel.rank == 1) then # 剛体回転
        delta_deg = (time_ax.val[t] - time_ax.val[0])*unit_time*(angvel[0..t].mean)*R2D % 360
        delta_i = (delta_deg/del_lon).round % lon.length
        # gp に鉛直軸がない（x,y,t）場合
        new_val[0..-(delta_i+1),true,t] = gp[(delta_i)..-1,true,t].val
        new_val[-(delta_i)..-1,true,t]  = gp[0..(delta_i-1),true,t].val
      else # 緯度ごとに角速度を変える
        lat.length.times{|j|
          delta_deg = (time_ax.val[t] - time_ax.val[0])*unit_time*(angvel[j,0..t].mean/cos(lat[j]*D2R))*R2D % 360
          delta_i = (delta_deg/del_lon).round % lon.length
          # gp に鉛直軸がない（x,y,t）場合
          new_val[0..-(delta_i+1),j,t] = gp[(delta_i)..-1,j,t].val
          new_val[-(delta_i)..-1,j,t]  = gp[0..(delta_i-1),j,t].val
        }
      end
    }
    ngp = gp.copy
    return ngp.replace_val(new_val)
  end

  def bulkRicherdsonNum(theta, theta_s, u, v, z=nil)
    g = GAnalysis::Met.g
    if (z.nil?) then
      theta.lost_axes.each{|la|
        axname, valunit = la.split("=")
        if (["lev","level","z","height"].include?(axname.downcase)) then
          z = UNumeric[valunit.to_f, "m"]
        end
      }
    elsif (z.class == Float)
      z = UNumeric[z, "m"]
    end
    return (theta - theta_s)*g*z/(theta_s*(u*u+v*v))
  end

  def surface_stress(theta, theta_s, u, v, rho, z=nil)
    # constant parameters (Louis, 1979)
    c_star = 7.4; b = 9.4; k = 0.4; z0 = UNumeric[0.01, "m"]

    if (z.nil?) then
      theta.lost_axes.each{|la|
        axname, valunit = la.split("=")
        if (["lev","level","z","height"].include?(axname.downcase)) then
          z = UNumeric[valunit.to_f, "m"]
        end
      }
    elsif (z.class == Float) then 
      z = UNumeric[z, "m"]
    end

    a2 = k*k/((z/z0).log)**2
    c  = c_star*a2*b*((z/z0).sqrt)
    rib = bulkRicherdsonNum(theta, theta_s, u, v, z) # BulkRicherdsonNum

    f_unstable   = 1.0 - (b*rib)/(1.0 + c*rib.abs.sqrt)
    f_stable = 1.0/(1.0+0.5*b*rib)**2

    theta_z = (theta - theta_s)/z # 安定度の計算はこれでよいのか？

    tau_s = rho*a2*(u*u+v*v)*( f_stable*(theta_z.gt(0)) + f_unstable*(theta_z.le(0)) )

    return tau_s
  end

  def surface_stress_L82(theta, theta_s, u, v, rho, z=nil)
    # constant parameters (Louis, 1982)
    k = 0.4; z0 = UNumeric[0.01, "m"] # 粗度長（本当は計算で求める）

    if (z.nil?) then
      theta.lost_axes.each{|la|
        axname, valunit = la.split("=")
        if (["lev","level","z","height"].include?(axname.downcase)) then
          z = UNumeric[valunit.to_f, "m"]
        end
      }
    elsif (z.class == Float) then 
      z = UNumeric[z, "m"]
    end

    a2 = k*k/((z/z0).log)**2
    rib = bulkRicherdsonNum(theta, theta_s, u, v, z) # BulkRicherdsonNum

    f_unstable   = 1.0 - (10.0*rib)/(1.0 + 75.0*a2*(((rib.abs)*(z/z0)).sqrt))
    f_stable = 1.0/(1.0+10.0*rib/((1.0+5.0*rib).sqrt) )

    theta_z = (theta - theta_s)/z # 安定度の計算はこれでよいのか？

    tau_s = rho*a2*(u*u+v*v)*( f_stable*(rib.gt(0)) + f_unstable*(rib.le(0)) )

    return tau_s
  end

  # ブラントバイサラ振動数を計算する｜N^2= g(∂T/∂z+ g/Cp)/T
  # temp: 温度
  def brunt_vaisala_N2(temp)
    g = GAnalysis::Met.g
    cp = GAnalysis::Met::Cp
    n2 = g*(temp.dz(temp) + g/cp)/temp    
    return n2
  end

  # 大気安定度 Γ - Γd を計算する
  def stability(temp)
    g = GAnalysis::Met.g
    cp = GAnalysis::Met::Cp
    s = temp.dz(temp) + g/cp
    return s
  end

  # 密度と気圧から静力学的な気圧を求める。z座標。
  def static_pressure(rho,prs)
    g = GAnalysis::Met.g
    zd, type = GPhys.vertical_axis_kind?(rho)
    if (type == "height") then 
      z = rho.coord(zd); kmax = z.length - 1  
      rho = rho.bring_axis_to_first(zd)
      pre = rho.copy
      pre.set_att("units","Pa"); pre.set_att("long_name", "static pressure"); pre.rename("pre")
      prs = prs.bring_axis_to_first(zd)
      # pre[kmax, false] = 0.0 
      pre[kmax, false] = prs[kmax,false].val.mean # 平均値を使うのがよい？
      (0..(kmax-1)).reverse_each{|k|
        pre[k,false] = pre[k+1,false] + (rho[k+1,false]+rho[k,false])*0.5*g*(z[k+1]-z[k])     
      }
      return pre.bring_firstaxis_to(zd)
      # ps を使って、下から計算
      # ps = ps.bring_axis_to_first(zd) 
      # pre[0,false] = ps[0,false] - rho[0,false]*g*z[0] 
      # # pre[0,false] = ps[0,false] - (rho[0,false]*z[1] - rho[1,false]*z[0])/(z[1]-z[0])*g*z[0] # rhoを外挿する場合
      # (1..kmax).each{|k|
      #   pre[k,false] = pre[k-1,false] - (rho[k,false]+rho[k-1,false])*0.5*g*(z[k]-z[k-1])     
      # }
      # return pre.bring_firstaxis_to(zd)
    else 
      raise "vertical axis must be z (height)."
    end
  end

  def median_filter2D(gp,n) 
    raise "num median filter2D must be odd number!!" if n.even?
    val = gp.val; new_val = NArrayMiss.sfloat(*(val.shape))
    im, jm = val.shape
    (n/2..(jm-n/2-1)).each{|j|
      (n/2..(im-n/2-1)).each{|i|
        new_val[i,j] = val[(i-n/2)..(i+n/2), (j-n/2)..(j+n/2)].median 
      }
    }
    return gp.replace_val(new_val)
  end




end
