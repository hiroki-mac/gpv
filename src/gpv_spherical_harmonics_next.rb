=begin
=DESCRIPTION
 ノルムを２に正規化した 球面調和関数ライブラリ
 依存ライブラリ: NArray, ruby-gsl
=AUTHOR
 Hiroki Kashimura

提供モジュール関数:
p_m_nT
deg2mu
gauss_w
sh_trans
sh_invtrans
lg_trans
lg_invtrans
multiply_cos
multiply_sin
sh_init
sh_dlon
sh_dlat
lg_dlat
sh_vorticity
sh_divergence
sh_energyspectrum
=end

require "narray"
require "gsl"
require "complex"
require "objspace"

class Complex
  def round(n)
    return Complex(self.real.round(n), self.imag.round(n))
  end
end
# class NArray
#   def round(n)
#     return self.collect{|v| v.round(n)}
#   end
# end



module SphericalHarmonics
  include NMath
  module_function
  if (RUBY_VERSION < "2.4" ) then
    Integer = Fixnum
  end

  # 緯度（Deg）をサイン緯度（μ）に変換
  # lat : 1D NArray
  def deg2mu(lat)
  	mu = NMath.sin(lat*PI/180.0)
  end

  def mu2deg(mu)
    lat = NMath.asin(mu)*180.0/PI
  end

  # gslを使って、ルジャンドル陪関数 P^m_n を求める
  # (三角切断；P^m_nがない要素は欠損値とする)
  def p_m_nT(x, mmax)
    mmax = mmax + 1
    len = x.length
    case @precision
    when "double" then
      pmn = NArray.float(len, mmax+2, mmax+2).fill!(-9.99E36) #微分計算のために n 方向に１つ余分に作っておく。m方向の余分はダミー(ゼロ)。
    when "single" then
      pmn = NArray.sfloat(len, mmax+2, mmax+2).fill!(-9.99E36) #微分計算のために n 方向に１つ余分に作っておく
    end

#    for m in 0..mmax
    np = (ENV['OMP_NUM_THREADS']|| 4).to_i
    pn_ary = Parallel.map(0..mmax, in_progress:np, :progress=>"progress"){|m|
      fact = sqrt(4.0*PI)*(-1.0)**m # GSLで計算されるP^m_nを2で正規化したものに変換する係数

      # 全てのルジャンドル陪関数をGSLで計算する方法（切断波数が高いと低速）
      # for l in m..mmax+1
      #   for j in 0..len-1
      #     # pmn[j, m, m..mmax+1] = (GSL::Sf::legendre_sphPlm_array(mmax+1, m, x[j]).to_na)*fact
      #     pmn[j, m, l] =  (GSL::Sf::legendre_sphPlm(l, m, x[j]))*fact
      #   end
      # end

      # 一部のルジャンドル陪関数をGSLで計算し、残りは３項漸化式を用いて計算する方法（高速）
      for j in 0..len-1
        # pmn[j, m, m..mmax+1] = (GSL::Sf::legendre_sphPlm_array(mmax+1, m, x[j]).to_na)*fact
        l = m  ; pmn[j, m, l] =  (GSL::Sf::legendre_sphPlm(l, m, x[j]))*fact
        l = m+1; pmn[j, m, l] =  (GSL::Sf::legendre_sphPlm(l, m, x[j]))*fact
      end
      # ３項漸化式（石岡 2004, p.120）
      for l in (m+1)..(mmax)
        pmn[true,m,l+1] = sqrt( (2.0*l+1.0)*(2.0*l+3.0)/((l+1.0)*(l+1.0)-m*m) )*x*pmn[true,m,l] \
                      - sqrt( (2.0*l+3.0)/(2.0*l-1.0) * (l*l-m*m)/((l+1.0)*(l+1.0)-m*m) )*pmn[true,m,l-1]
      end

      pmn[true,m,true]
    }


    for m in 0..mmax
      pmn[true,m,true] = pn_ary[m]
    end

#    end
    pmn[true,0,0] = 1.0
    pmn[true,mmax+1,mmax+1] = 0.0 # ダミー

    # ルジャンドル陪関数の微分の計算

    case @precision
    when "double" then
      dpmn = NArray.float(len, mmax+1, mmax+1).fill!(-9.99E36)
      emn = NArray.float(mmax+1, mmax+2)    # 係数を計算しておく
    when "single" then
      dpmn = NArray.sfloat(len, mmax+1, mmax+1).fill!(-9.99E36)
      emn = NArray.sfloat(mmax+1, mmax+2)    # 係数を計算しておく
    end

    dpn_ary = Parallel.map(0..mmax, in_progress:np, :progress=>"progress"){|m|

#    for m in 0..mmax
      for n in (m+1)..(mmax+1)
	       emn[m, n] = sqrt( (n*n - m*m)/((2.0*n+1.0)*(2.0*n-1.0)) )
      end
#    end

    x_tmp = 1.0/(1.0 - x*x)
#    for m in 0..mmax
      dpmn[true, m, m] = ( - m*emn[m,m+1]*pmn[true,m,m+1]) *x_tmp # m = n の場合
      for n in (m+1)..mmax
        dpmn[true, m, n] = ((n+1.0)*emn[m,n]*pmn[true,m,n-1] - n*emn[m,n+1]*pmn[true,m,n+1]) *x_tmp
      end
      dpmn[true,m,true]
#    end
    } 
    for m in 0..mmax
      dpmn[true,m,true] = dpn_ary[m]
    end


#    binding.pry
#    return pmn, dpmn, emn
    return pmn[true,0..-2,0..-2], dpmn[true,0..-1,0..-1], emn[0..-2,0..-2]
  end

  # サインガウス緯度（sin(θ_k)）とガウス重み（w_k）を求める
  # lat: Integer or 1D NArray
  def gauss_w(lat)
    if (lat.class == Integer) then
      jmax = lat
    else
      jmax = lat.size
    end

    case @precision
    when "single" then
      mu = NArray.sfloat(jmax)
      gw = NArray.sfloat(jmax)
      pjmax_m1 = NArray.sfloat(jmax)
      pjmax_p1 = NArray.sfloat(jmax)
      pjmax    = NArray.sfloat(jmax)
    else
      mu = NArray.float(jmax)
      gw = NArray.float(jmax)
      pjmax_m1 = NArray.float(jmax)
      pjmax_p1 = NArray.float(jmax)
      pjmax    = NArray.float(jmax)
    end

    # ニュートン法でガウス緯度を計算する。（参考：石岡 2004）
    # 初期推定値。muが与えられているときはそれを初期推定値とする。
    if (lat.class == Integer) then
      for j in 0..jmax-1
        mu[j] = cos(PI*(jmax-j-0.25)/(jmax+0.5))
      end
    else
      mu = deg2mu(lat)
    end

    # 必要なルジャンドル多項式を求める
    50.times{|t|
      for j in 0..jmax-1
        pjmax_m1[j] = GSL::Sf::legendre_Pl(jmax-1, mu[j])*sqrt(2.0*jmax-1.0)
        pjmax_p1[j] = GSL::Sf::legendre_Pl(jmax+1, mu[j])*sqrt(2.0*jmax+3.0)
        pjmax[j]    = GSL::Sf::legendre_Pl(jmax, mu[j])*sqrt(2.0*jmax+1.0)
      end
      n = jmax
      dpjmax = ( (n+1.0)*n/sqrt((2.0*n+1.0)*(2.0*n-1.0))*pjmax_m1 - \
                n*(n+1.0)/sqrt((2.0*n+3.0)*(2.0*n+1.0))*pjmax_p1 )/(1.0-mu*mu)
      # 次の推定値
      mu_next = mu - pjmax/dpjmax
      # 収束判定
      if ((mu - mu_next).abs.max < 1.0E-15) then
        mu = mu_next ; print "mu has calculated by #{t+1} times iteration.\n"; break
      else
        mu = mu_next ; print "mu calculation was not converged. \n" if (t == 49)
      end
    }

    # ガウス重みを計算する
    for j in 0..jmax-1
      pjmax_m1[j] = GSL::Sf::legendre_Pl(jmax-1, mu[j])*sqrt(2.0*jmax-1.0)
      pjmax_p1[j] = GSL::Sf::legendre_Pl(jmax+1, mu[j])*sqrt(2.0*jmax+3.0)
    end
    n = jmax
    dpjmax = ( (n+1.0)*n/sqrt((2.0*n+1.0)*(2.0*n-1.0))*pjmax_m1 - \
              n*(n+1.0)/sqrt((2.0*n+3.0)*(2.0*n+1.0))*pjmax_p1 )/(1.0-mu*mu)
    gw = 2.0*sqrt((2.0*n-1.0)*(2.0*n+1.0))/(n*pjmax_m1*dpjmax)

#    binding.pry
    10.times{ gw = gw*2.0/gw.to_a.sum } # test ガウス重みの精度が上がる？

#    binding.pry

    return mu, gw
  end





  # 球面調和正変換
  # 引数：グリッドデータ、P_n^m、ガウス重み
  # グリッドデータgijは1次元目がlon、2次元目がlatでなければならない。
  # 第4引数に正整数を与えると、その数(only)の全波数のみ計算する。

  def sh_trans(gij, dlon=false, only=false, nmax=@pmn[0,0,true].size-2 )
  shapes = gij.shape
  imax = shapes[0]
  jmax = shapes[1]

  case @precision
  when "double" then
    if (shapes.size == 2) then
  #	  smn = NArrayMiss.complex(imax, nmax+1).fill!(-9.99E36).all_invalid
  	  smn = NArray.complex(imax, nmax+1)#.fill!(-9.99E36)
    elsif (shapes.size == 3) then
  # 	  smn = NArrayMiss.complex(imax, nmax+1, shapes[2]).fill!(-9.99E36).all_invalid
   	  smn = NArray.complex(imax, nmax+1, shapes[2])#.fill!(-9.99E36)
    elsif (shapes.size == 4) then
  # 	  smn = NArrayMiss.complex(imax, nmax+1, shapes[2], shapes[3]).fill!(-9.99E36).all_invalid
   	  smn = NArray.complex(imax, nmax+1, shapes[2], shapes[3])#.fill!(-9.99E36)
    else
    		raise "number of input grid data dimension must be less than 5"
    end
  when "single" then
    if (shapes.size == 2) then
  #	  smn = NArrayMiss.scomplex(imax, nmax+1).fill!(-9.99E36).all_invalid
  	  smn = NArray.scomplex(imax, nmax+1)#.fill!(-9.99E36)
    elsif (shapes.size == 3) then
  # 	  smn = NArrayMiss.scomplex(imax, nmax+1, shapes[2]).fill!(-9.99E36).all_invalid
   	  smn = NArray.scomplex(imax, nmax+1, shapes[2])#.fill!(-9.99E36)
    elsif (shapes.size == 4) then
  # 	  smn = NArrayMiss.scomplex(imax, nmax+1, shapes[2], shapes[3]).fill!(-9.99E36).all_invalid
   	  smn = NArray.scomplex(imax, nmax+1, shapes[2], shapes[3])#.fill!(-9.99E36)
    else
    		raise "number of input grid data dimension must be less than 5"
    end
  end


  # FFT
  if (shapes.size == 4) then #サイズの大きな配列に対応するため。
    case @precision
    when "double" then
      gmj = NArray.complex(shapes[0], shapes[1], shapes[2], shapes[3])
    when "single" then
      gmj = NArray.scomplex(shapes[0], shapes[1], shapes[2], shapes[3])
    end
    shapes[3].times{|t|
    	gmj[true, true, true, t] = gij[true, true, true, t].fft(false, 0).val
    }
  else
    gmj = gij.fft(false, 0).val
  end

  gjm = gmj.transpose(1,0) if gmj.rank == 2
  gjm = gmj.transpose(1,0,2..-1) if gmj.rank >= 3

  # ルジェンドル正変換
  if (@pmn.class == NArrayMiss) then
 	  pmn_na = @pmn.to_na
  else
	  pmn_na = @pmn
  end
	for n in 0..nmax # 三角切断
    # if (@flag_pmn_file) then 
    #   pmn_na = @gpmn[0..n,n,true].val.transpose(1,0) # @gpmnのindex は m, n, j の順番なので入れ替える
    #   smn[0..n,n, false] = ( (gjm[true, 0..n, false].mul_add(pmn_na.mul!(@gw[true]), 0 )).mul!(0.5) ) #.round(13) # test
    #   smn[-1..-n,n,false] = smn[1..n,n,false].conj if (n != 0)      
    if (only == false) then
      smn[0..n,n, false] = ( (gjm[true, 0..n, false].mul_add(pmn_na[true, 0..n, n].mul!(@gw[true]), 0 )).mul!(0.5) ) #.round(13) # test
  	  smn[-1..-n,n,false] = smn[1..n,n,false].conj if (n != 0)
    elsif (n == only) then # 特定の全波数 n だけ変換する場合
  	  smn[0..n,n, false] = (gjm[true, 0..n,false].mul_add(pmn_na[true, 0..n, n].mul!(@gw[true]), 0 )).mul!(0.5)
  	  smn[-1..-n,n,false] = smn[1..n,n,false].conj if (n != 0)
    end
	end

	# 経度微分をする場合
	if (dlon) then
    dlon = 1 if (dlon.class != Integer)
		dlon.times{
        smn[0..nmax] = smn[0..nmax].mul!(@complex_im)
	      smn[-1..-(nmax), false] = smn[1..(nmax), false].conj # m < 0 部分の値
		}
	end

	# GPhysオブジェクト化
	# 軸の設定
	old_grid = gij.grid_copy #入力GPhysオブジェクトのグリッド情報を取得
	xaxis = Axis.new
	zwn = VArray.new(NArray.int(imax).indgen!,{"long_name"=>"zonal wave number","units"=>"wave number"},"zwn")
	xaxis.pos = zwn
	yaxis = Axis.new
	twn = VArray.new(NArray.int(nmax+1).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"twn")
	yaxis.pos = twn
	new_grid = old_grid.change_axis(0, xaxis).change_axis(1, yaxis) #第0, 1軸を置き換える

	# 変数の設定
	sname = gij.name
	lname = gij.long_name
	unit = gij.units.to_s

	va = VArray.new(smn, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_smn = GPhys.new(new_grid, va)
  return gphys_smn
  end

# 球面調和逆変換
  # 引数：スペクトルデータ、P_n^m、経度(Deg) 、緯度(Deg)
  # 第8引数に正整数を与えると、その数(only)の全波数のみ計算する。
  def sh_invtrans(smn, dlon=false, dlat=false, only=false, nmax=@pmn[0,0,true].size-2)
  jmax = @lat.size
  imax = @lon.size
#  mmax = @pmn[0,true,0].size
  shapes = smn.shape
  case @precision
  when "double" then
    if (shapes.size == 2) then
  #	  gij = NArrayMiss.float(imax, jmax).fill!(0.0)
  #	  gmj = NArrayMiss.complex(imax, jmax).fill!(0.0)
  	  gij = NArray.float(imax, jmax)#.fill!(0.0)
  	  gmj = NArray.complex(imax, jmax)#.fill!(0.0)
    elsif (shapes.size == 3) then
  #	  gij = NArrayMiss.float(imax, jmax, shapes[2]).fill!(0.0)
  #	  gmj = NArrayMiss.complex(imax, jmax, shapes[2]).fill!(0.0)
  	  gij = NArray.float(imax, jmax, shapes[2])#.fill!(0.0)
  	  gmj = NArray.complex(imax, jmax, shapes[2])#.fill!(0.0)
    elsif (shapes.size == 4) then
  #	  gij = NArrayMiss.float(imax, jmax, shapes[2], shapes[3]).fill!(0.0)
  #	  gmj = NArrayMiss.complex(imax, jmax, shapes[2], shapes[3]).fill!(0.0)
  	  gij = NArray.float(imax, jmax, shapes[2], shapes[3])#.fill!(0.0)
  	  gmj = NArray.complex(imax, jmax, shapes[2], shapes[3])#.fill!(0.0)
    else
    		raise "number of input grid data dimension must be less than 5"
    end
  when "single" then
    if (shapes.size == 2) then
  #	  gij = NArrayMiss.sfloat(imax, jmax).fill!(0.0)
  #	  gmj = NArrayMiss.scomplex(imax, jmax).fill!(0.0)
  	  gij = NArray.sfloat(imax, jmax)#.fill!(0.0)
  	  gmj = NArray.scomplex(imax, jmax)#.fill!(0.0)
    elsif (shapes.size == 3) then
  #	  gij = NArrayMiss.sfloat(imax, jmax, shapes[2]).fill!(0.0)
  #	  gmj = NArrayMiss.scomplex(imax, jmax, shapes[2]).fill!(0.0)
  	  gij = NArray.sfloat(imax, jmax, shapes[2])#.fill!(0.0)
  	  gmj = NArray.scomplex(imax, jmax, shapes[2])#.fill!(0.0)
    elsif (shapes.size == 4) then
  #	  gij = NArrayMiss.sfloat(imax, jmax, shapes[2], shapes[3]).fill!(0.0)
  #	  gmj = NArrayMiss.scomplex(imax, jmax, shapes[2], shapes[3]).fill!(0.0)
  	  gij = NArray.sfloat(imax, jmax, shapes[2], shapes[3])#.fill!(0.0)
  	  gmj = NArray.scomplex(imax, jmax, shapes[2], shapes[3])#.fill!(0.0)
    else
    		raise "number of input grid data dimension must be less than 5"
    end
  end
  smn_tmp = smn.val

  # ルジェンドル逆変換
  if (dlat == false ) then
    if (@pmn.class == NArrayMiss) then
  	 	pmn_na = @pmn.to_na
  	else
  		pmn_na = @pmn
  	end
	  for j in 0..jmax-1 #三角切断
		  if (only == false) then
			  gmj[0..nmax, j, false] =  pmn_na[j, 0..nmax, 0..nmax].mul_add(smn_tmp[0..nmax, 0..nmax, false], 1)
        gmj[-1..-(nmax), j, false] = gmj[1..(nmax), j, false].conj # m < 0 部分の値
		  else
        gmj[0..only, j, false] =  pmn_na[j, 0..only, only].mul!(smn_tmp[0..only, only, false])
        gmj[-1..-only, j, false] = gmj[1..only, j, false].conj # m < 0 部分の値
		  end
		end

	elsif ( dlat ) then #緯度微分を作用させたルジャンドル逆変換

	 	if (@pmn.class == NArrayMiss) then
	 	  dpmn_na = @dpmn.to_na
	  else
		  dpmn_na = @dpmn
    end
		for j in 0..jmax-1 #三角切断
		  if (only == false) then
			  gmj[0..nmax, j, false] =  dpmn_na[j, 0..nmax, 0..nmax].mul_add( smn_tmp[0..nmax, 0..nmax, false], 1)
			  gmj[-1..-(nmax),j, false] = gmj[1..(nmax),j, false].conj # m < 0 部分の値
		  else #if (m <= only) then # 特定の全波数 n だけ変換する場合
        gmj[0..only, j, false] =  dpmn_na[j, 0..only, only].mul!(smn_tmp[0..only, only, false])
			  gmj[-1..-only,j, false] = gmj[1..only,j, false].conj # m < 0 部分の値
		  end
		end

	end

	# 経度微分をする場合
	if (dlon) then
    dlon = 1 if (dlon.class != Integer)
			dlon.times{
  			gmj[0..(nmax),false] = gmj[0..(nmax),false].mul!(@complex_im) # I は虚数単位
        gmj[-1..-(nmax), false] = gmj[1..(nmax), false].conj # m < 0 部分の値
			}
#      binding.pry
	end

	# 逆FFT
	if (gmj.class == NArrayMiss) then
		gij = FFTW3.fft(gmj.to_na, 1, 0).real
	else
		gij = FFTW3.fft(gmj, 1, 0).real
	end

#  print gij[true,0,0].mean,", ", gmj[0,0,0],", ",  smn.name, smn.units, "\n"

  case @precision
  when "double" then
  when "single" then
    gij = gij.to_type("sfloat") #FFTW3は倍精度計算なので、単精度に戻す
  end




	#GPhysオブジェクト化
	old_grid = smn.grid_copy
	xaxis = Axis.new
	vlon = VArray.new(@lon,{"longitude"=>" ","units"=>"deg"},"lon")
	xaxis.pos = vlon
	yaxis = Axis.new
	vlat = VArray.new(@lat,{"latitude"=>" ","units"=>"deg"},"lat")
	yaxis.pos = vlat
	new_grid = old_grid.change_axis(0, xaxis).change_axis(1, yaxis) #第0, 1軸を置き換える

	# 変数の設定
	sname = smn.name
	lname = smn.long_name
	unit = smn.units.to_s

  va = VArray.new(gij, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_gmj= GPhys.new(new_grid, va)
  return gphys_gmj
  end




  # ルジャンドル正変換
  # 引数：グリッドデータ、P_n^0、ガウス重み
  #グリッドデータgjは1次元目が lat でなければならない。

  def lg_trans(gj, nmax=@pmn[0,0,true].size-2)
  shapes = gj.shape
  jmax = shapes[0]
  pn = @pmn[true, 0, true]
  dpn = @dpmn[true, 0, true]

  case @precision
  when "double" then
    if (shapes.size == 1) then
  	  sn = NArray.float(nmax+1).fill!(-9.99E36)
    elsif (shapes.size == 2) then
   	  sn = NArray.float(nmax+1, shapes[1]).fill!(-9.99E36)
    elsif (shapes.size == 3) then
   	  sn = NArray.float(nmax+1, shapes[1], shapes[2]).fill!(-9.99E36)
    else
    		raise "number of input grid data dimension must be less than 4"
    end
  when "single" then
    if (shapes.size == 1) then
  	  sn = NArray.sfloat(nmax+1).fill!(-9.99E36)
    elsif (shapes.size == 2) then
   	  sn = NArray.sfloat(nmax+1, shapes[1]).fill!(-9.99E36)
    elsif (shapes.size == 3) then
   	  sn = NArray.sfloat(nmax+1, shapes[1], shapes[2]).fill!(-9.99E36)
    else
    		raise "number of input grid data dimension must be less than 4"
    end
  end


  gj_na = gj.val
  gj_na = gj_na.to_na if (gj_na.class == NArrayMiss)

  # ルジェンドル正変換
 	if (pn.class == NArrayMiss) then
	 	pn_na = pn.to_na
	else
		pn_na = pn
	end
	for n in 0..nmax
	  sn[n, false] = (0.5* ( gj_na[true, false].mul_add(pn_na[true, n].mul!(@gw[true]), 0 )	)  ) # .round(13) # test
	end

	# GPhysオブジェクト化
	# 軸の設定
	old_grid = gj.grid_copy #入力GPhysオブジェクトのグリッド情報を取得
	yaxis = Axis.new
	twn = VArray.new(NArray.int(nmax+1).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"twn")
	yaxis.pos = twn
	new_grid = old_grid.change_axis(0, yaxis) #第 0 軸を置き換える

	# 変数の設定
	sname = gj.name
	lname = gj.long_name
	unit = gj.units.to_s

	va = VArray.new(sn, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_sn = GPhys.new(new_grid, va)

  end


# ルジャンドル逆変換
  # 引数：スペクトルデータ、P_n^0、緯度(Deg)
  def lg_invtrans(sn, dlat=false, nmax=@pmn[0,0,true].size-2)
  jmax = @lat.size
  pn = @pmn[true, 0, true]
  dpn = @dpmn[true, 0, true]

  shapes = sn.shape
  case @precision
  when "double" then
    if (shapes.size == 1) then
      gj = NArray.float(jmax)#.fill!(0.0)
    elsif (shapes.size == 2) then
      gj = NArray.float(jmax, shapes[1])#.fill!(0.0)
    elsif (shapes.size == 3) then
      gj = NArray.float(jmax, shapes[1], shapes[2])#.fill!(0.0)
    else
    		raise "number of input grid data dimension must be less than 4"
    end
  when "single" then
    if (shapes.size == 1) then
      gj = NArray.sfloat(jmax)#.fill!(0.0)
    elsif (shapes.size == 2) then
      gj = NArray.sfloat(jmax, shapes[1])#.fill!(0.0)
    elsif (shapes.size == 3) then
      gj = NArray.sfloat(jmax, shapes[1], shapes[2])#.fill!(0.0)
    else
    		raise "number of input grid data dimension must be less than 4"
    end
  end
  sn_tmp = sn.val



  # ルジェンドル逆変換
  if (dlat == false ) then
   	if (pn.class == NArrayMiss) then
	 	pn_na = pn.to_na
	else
		pn_na = pn
	end

	for j in 0..jmax-1 #三角切断
	  gj[j, false] =  pn_na[j,0..nmax].mul_add(sn_tmp[0..nmax, false], 0)
	end

	elsif ( dlat ) then #緯度微分を作用させたルジャンドル逆変換
	 	if (pn.class == NArrayMiss) then
		 	dpn_na = dpn.to_na
		else
			dpn_na = dpn
	  end

		for j in 0..jmax-1 #三角切断
		  gj[j, false] =  dpn_na[j,0..nmax].mul_add( sn_tmp[ 0..nmax, false], 0)
		end
	end

	#GPhysオブジェクト化
	old_grid = sn.grid_copy
	yaxis = Axis.new
	vlat = VArray.new(@lat,{"latitude"=>" ","units"=>"deg"},"lat")
	yaxis.pos = vlat
	new_grid = old_grid.change_axis(0, yaxis) #第 0軸を置き換える

	# 変数の設定
	sname = sn.name
	lname = sn.long_name
	unit = sn.units.to_s

  va = VArray.new(gj, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_gj= GPhys.new(new_grid, va)

  end





# cos(lat)をかける
def multiply_cos(gphys, num=1)
  # モジュール変数を利用する方法（高速）
  if (@coslat) then
    gshape = gphys.shape
    return gphys*@coslat**num if (gshape == @coslat.shape)
    case @coslat.rank
    when 2
      return gphys*@coslat[0,true]**num if (gshape == @coslat[0,true].shape)
      return gphys*@coslat[true,0]**num if (gshape == @coslat[true,0].shape)
    when 3
      return gphys*@coslat[0,true,true]**num if (gshape == @coslat[0,true,true].shape)
      return gphys*@coslat[true,0,true]**num if (gshape == @coslat[true,0,true].shape)
      return gphys*@coslat[true,true,0]**num if (gshape == @coslat[true,true,0].shape)
      return gphys*@coslat[0,0,true]**num if (gshape == @coslat[0,0,true].shape)
      return gphys*@coslat[0,true,0]**num if (gshape == @coslat[0,true,0].shape)
      return gphys*@coslat[true,0,0]**num if (gshape == @coslat[true,0,0].shape)
    when 4..10
      return gphys*@coslat[0,true,true,true,false]**num if (gshape == @coslat[0,true,true,true,false].shape)
      return gphys*@coslat[true,0,true,true,false]**num if (gshape == @coslat[true,0,true,true,false].shape)
      return gphys*@coslat[true,true,0,true,false]**num if (gshape == @coslat[true,true,0,true,false].shape)
      return gphys*@coslat[true,true,true,0,false]**num if (gshape == @coslat[true,true,true,0,false].shape)
      return gphys*@coslat[0,0,true,true,false]**num if (gshape == @coslat[0,0,true,true,false].shape)
      return gphys*@coslat[true,0,0,true,false]**num if (gshape == @coslat[true,0,0,true,false].shape)
      return gphys*@coslat[true,true,0,0,false]**num if (gshape == @coslat[true,true,0,0,false].shape)
      return gphys*@coslat[0,true,true,0,false]**num if (gshape == @coslat[0,true,true,0,false].shape)
      return gphys*@coslat[0,0,0,true,false]**num if (gshape == @coslat[0,0,0,true,false].shape)
      return gphys*@coslat[true,0,0,0,false]**num if (gshape == @coslat[true,0,0,0,false].shape)
      return gphys*@coslat[0,true,0,0,false]**num if (gshape == @coslat[0,true,0,0,false].shape)
      return gphys*@coslat[0,0,true,0,false]**num if (gshape == @coslat[0,0,true,0,false].shape)
    end
  end

  # 上記が使えない場合（低速）
	# 緯度軸を見つける
	gresult = gphys.copy
    axarray = gresult.axnames
   latpos = 999
	for ln in ["lat", "Lat", "LAT", "latitude", "Latitude", "LATITUDE", "y"]
		if (axarray.index(ln)) then
			lat =  gresult.axis(ln).to_gphys.val
			latpos = axarray.index(ln)
			break
		end
	end
	if (latpos == 999) then
				raise "It seems that the GPhys object does not have lat axis."
	end
	# cos(lat) 作成
	coslat = NMath.cos(lat*PI/180.0)
	if (latpos == 0) then
#		lat.size.times{|j| gresult[j, false] = gphys[j, false]*coslat[j]**num	}
    gresult = gphys*coslat**num
	elsif (latpos == 1) then
		lat.size.times{|j| gresult[true, j, false] = gphys[true, j, false]*coslat[j]**num }
	elsif (latpos == 2) then
		lat.size.times{|j|
			gresult[true, true, j, false] = gphys[true, true, j, false]*coslat[j]**num
		}
	elsif (latpos == 3) then
		lat.size.times{|j|
			gresult[true, true, true, j, false] = gphys[true, true, true, j, false]*coslat[j]**num
		}
	else
		raise "lat should be within 4th dimension"
	end
	return gresult
end








# sin(lat)をかける
def multiply_sin(gphys, num=1)
  # モジュール変数を利用する方法（高速）
  if (@sinlat) then
    gshape = gphys.shape
    return gphys*@sinlat**num if (gshape == @sinlat.shape)
    case @sinlat.rank
    when 2
      return gphys*@sinlat[0,true]**num if (gshape == @sinlat[0,true].shape)
      return gphys*@sinlat[true,0]**num if (gshape == @sinlat[true,0].shape)
    when 3
      return gphys*@sinlat[0,true,true]**num if (gshape == @sinlat[0,true,true].shape)
      return gphys*@sinlat[true,0,true]**num if (gshape == @sinlat[true,0,true].shape)
      return gphys*@sinlat[true,true,0]**num if (gshape == @sinlat[true,true,0].shape)
      return gphys*@sinlat[0,0,true]**num if (gshape == @sinlat[0,0,true].shape)
      return gphys*@sinlat[0,true,0]**num if (gshape == @sinlat[0,true,0].shape)
      return gphys*@sinlat[true,0,0]**num if (gshape == @sinlat[true,0,0].shape)
    when 4..10
      return gphys*@sinlat[0,true,true,true,false]**num if (gshape == @sinlat[0,true,true,true,false].shape)
      return gphys*@sinlat[true,0,true,true,false]**num if (gshape == @sinlat[true,0,true,true,false].shape)
      return gphys*@sinlat[true,true,0,true,false]**num if (gshape == @sinlat[true,true,0,true,false].shape)
      return gphys*@sinlat[true,true,true,0,false]**num if (gshape == @sinlat[true,true,true,0,false].shape)
      return gphys*@sinlat[0,0,true,true,false]**num if (gshape == @sinlat[0,0,true,true,false].shape)
      return gphys*@sinlat[true,0,0,true,false]**num if (gshape == @sinlat[true,0,0,true,false].shape)
      return gphys*@sinlat[true,true,0,0,false]**num if (gshape == @sinlat[true,true,0,0,false].shape)
      return gphys*@sinlat[0,true,true,0,false]**num if (gshape == @sinlat[0,true,true,0,false].shape)
      return gphys*@sinlat[0,0,0,true,false]**num if (gshape == @sinlat[0,0,0,true,false].shape)
      return gphys*@sinlat[true,0,0,0,false]**num if (gshape == @sinlat[true,0,0,0,false].shape)
      return gphys*@sinlat[0,true,0,0,false]**num if (gshape == @sinlat[0,true,0,0,false].shape)
      return gphys*@sinlat[0,0,true,0,false]**num if (gshape == @sinlat[0,0,true,0,false].shape)
    end
  end
  # 上記が使えない場合（低速）
	# 緯度軸を見つける
	gresult = gphys.copy
    axarray = gresult.axnames
   latpos = 999
	for ln in ["lat", "Lat", "LAT", "latitude", "Latitude", "LATITUDE", "y"]
		if (axarray.index(ln)) then
			lat =  gresult.axis(ln).to_gphys.val
			latpos = axarray.index(ln)
			break
		end
	end
	if (latpos == 999) then
				raise "It seems that the GPhys object does not have lat axis."
	end

	# sin(lat) 作成
	sinlat = NMath.sin(lat*PI/180.0)
	if (latpos == 0) then
#		lat.size.times{|j| gresult[j, false] = gphys[j, false]*sinlat[j]**num }
    gresult = gphys*sinlat**num
	elsif (latpos == 1) then
		lat.size.times{|j|
			gresult[true, j, false] = gphys[true, j, false]*sinlat[j]**num
#			gresult[true, j, false] = gphys[true, j, false]*0.0 + sinlat[j]**num
		}
	elsif (latpos == 2) then
		lat.size.times{|j|
			gresult[true, true, j, false] = gphys[true, true, j, false]*sinlat[j]**num
		}
	elsif (latpos == 3) then
		lat.size.times{|j|
			gresult[true, true, true, j, false] = gphys[true, true, true, j, false]*sinlat[j]**num
		}
	else
		raise "lat should be within 4th dimension"
	end

	return gresult
end



############################################
## 以下は、上の関数を呼び出しやすいようにしたラッパー##
############################################

# 初期化：Pmnとガウス重みを求める
# 引数：経度(Deg), 緯度(Deg), [切断波数：与えなければ緯度格子の数をもとに自動で計算する]
#      変換するGPhysサンプル（次元を参照する）、変換精度（single or double）
# ガウス緯度を内部計算する場合は、第２引数に、緯度方向の格子点数を与える。
# また、第２引数に何も与えなければ、lon/2を緯度方向の格子点数とする。
def sh_init(lon, lat = lon.size/2, mmax = ((lat.size*2-1)/3), gp=nil, precision=nil, pmn_file=nil)
  print "Truncation wavenumber is #{mmax}.\n"

   # 倍精度/単精度の選択
  if (precision) then
    @precision = precision
  elsif (gp) then
    case gp.ntype
    when "float" then
      @precision = "double"
    when "sfloat" then
      @precision = "single"
    end
  else
    @precision = "double" #デフォルトの精度
  end

  if (lat.class == Integer) then # ガウス緯度を内部計算
    @mu, @gw = gauss_w(lat) # ガウス緯度・重み
    @lon = lon;  @lat = mu2deg(@mu)
  else # 与えられた、緯度を用いる。ただし、サイン緯度・重みの計算には精確なガウス緯度を計算する。
    @lon = lon;  @lat = lat
    @mu, @gw = gauss_w(lat) #ガウス重み
  end
  print "Gaussian weight and latitude generated.\n"
  if (pmn_file == nil) then 
    @pmn, @dpmn, @emn = p_m_nT(@mu, mmax) #三角切断 P^m_n, dP^m_n
    @flag_pmn_file = false

    print "P^m_n, dP^m_n generated.\n"
  else 
    print "P^m_n, dP^m_n will be loaded from file: #{pmn_file}"
    @pmn, @dpmn, @emn = p_m_nT(@mu[0..0], mmax) #三角切断 P^m_n, dP^m_n
    @gpmn  = GPhys::IO.open_gturl("pmn_gl10.nc@pmn")
    @gdpmn = GPhys::IO.open_gturl("pmn_gl10.nc@dpmn")
    @flag_pmn_file = true
  end

  # モジュール変数の作成
  if (gp) then
    @coslat = gp.val*0.0
    lat.length.times{|j| @coslat[true,j,false] = cos(lat[j]*PI/180.0)}
    @sinlat = gp.val*0.0
    lat.length.times{|j| @sinlat[true,j,false] = sin(lat[j]*PI/180.0)}
    #     @Cp05 = gp.val*0.0 + 0.5 #定数配列 0.5
    #     @Cm05 = -@Cp05 # 定数配列
  end

  case @precision
  when "double" then
    @complex_im = NArray.complex(mmax+1)
  when "single" then
    @complex_im = NArray.scomplex(mmax+1)
  end
  (mmax+1).times{|m| @complex_im[m] = Complex::I*m }

#  binding.pry

  return @gw, @pmn, @dpmn
end

def coslat()
  return @coslat
end

def sinlat()
  return @sinlat
end
      
# 任意の点数のガウス緯度を求める。これは初期化を呼ばなくても使える。
# lat: Integer or NArray
def gauss_lat(lat) 
  mu, wg = gauss_w(lat)
  return mu2deg(mu)
end
      

def sh_trans_multi(gijarray, dlon=false, only=false)
  #p ObjectSpace.memsize_of(gparray)
  num = gijarray.length
  smn = sh_trans(GPhys.concat(gijarray, NArray.sfloat(num).indgen!, "tmp_axis"), dlon, only)
  smnarray = []
  num.times{|n|
    smnarray << smn.cut("tmp_axis"=>n)
  }
  return smnarray
end

def sh_invtrans_multi(smnarray, dlon=false, dlat=false, only=false)
  #p ObjectSpace.memsize_of(gparray)
  num = smnarray.length
  gij = sh_invtrans(GPhys.concat(smnarray, NArray.sfloat(num).indgen!, "tmp_axis"), dlon, dlat, only)
  gijarray = []
  num.times{|n|
    gijarray << gij.cut("tmp_axis"=>n)
  }
  return gijarray
end





# 経度微分d/dlon（グリッド→スペクトルー(経度微分)→グリッド）
# 引数：GPhysオブジェクト、Pmn、ガウス重み、[経度微分の階数]
def sh_dlon(gphys, dn=1)
  if (gphys.class == GPhys) then
  	#グリッド→スペクトル
  	smn = sh_trans(gphys)
#    smn[true,-1,false]=0.0 # test
  	#スペクトル→グリッド（経度微分：d/dlon）
  	return sh_invtrans(smn,dn)
  elsif (gphys.class == Array) then
    smn = sh_invtrans_multi(gphys)
    return sh_invtrans_multi(smn,dn)
  else
    raise "arg must be GPhys or Array of GPhyses which are same shape)"
  end
end


# 緯度微分d/dmu（グリッド→スペクトルー(緯度微分)→グリッド）
# 引数：GPhysオブジェクト
def sh_dlat(gphys)
  if (gphys.class == GPhys) then
  	#グリッド→スペクトル
  	smn = sh_trans( gphys, false, false, @pmn[0,0,true].size-1)
  	#スペクトル→グリッド（緯度微分：d/dmu）
    return sh_invtrans(smn, false, true, false, @pmn[0,0,true].size-1)


#   スペクトル→ cos^2(lat)*d/dmu の スペクトル
    #グリッド/cos^2(lat)→スペクトル/cos^2(lat)
    smn = sh_trans( gphys/(@coslat*@coslat), false, false, @pmn[0,0,true].size-1)
    #スペクトル/cos^2(lat)→スペクトル（緯度微分：d/dmu）
    d_smn = sp_dmu_cos2(smn)
    #スペクトル（緯度微分：d/dmu）→グリッド（緯度微分：d/dmu）
    d_gij = sh_invtrans(d_smn, false, false, false, @pmn[0,0,true].size-2)
    return d_gij

#  	return d_gij/(@coslat*@coslat)

  elsif (gphys.class == Array) then
    smn = sh_trans_multi(gphys)
    return sh_invtrans_multi(smn, false, true)
  else
    raise "arg must be GPhys or Array of GPhyses which are same shape)"
  end

end

# スペクトルデータから、その経度微分のスペクトルを求める。
def sp_dlon(gphys_smn)
  d_gphys_smn = gphys_smn.copy
  nmax = @complex_im.shape[0] - 1
  smn = gphys_smn.val
	smn[0..nmax,false] = smn[0..nmax,false].mul!(@complex_im) # I は虚数単位
  smn[-1..-nmax, false] = smn[1..nmax, false].conj # m < 0 部分の値
  d_gphys_smn.replace_val(smn)
  return d_gphys_smn
end

# スペクトルデータから、そのcos^2(lat)*mu微分 (mu = sin(lat)) のスペクトルを求める。
def sp_cos2_dmu(gphys_smn)
  smn = gphys_smn.val

  imax = smn.shape[0]
  nmax = @pmn[0,0,true].size-2 ## smn.shape[1] - 1 # 切断波数

#  sval = NArray.complex(imax,nmax+2,1)
  if smn.rank == 2 then
    sval = NArray.complex(imax,nmax+2)
  else
    shape_ary = [imax,nmax+2] + smn.shape[2..-1]
    sval = NArray.complex(*shape_ary)
  end


#binding.pry
  sval[0..nmax,0,false] = (2)*smn[0..nmax,1,false]*@emn[0..nmax,1]
  for n in 1..nmax-1
    sval[0..nmax,n,false] = (n+2)*smn[0..nmax,n+1,false]*@emn[0..nmax,n+1] - (n-1)*smn[0..nmax,n-1,false]*@emn[0..nmax,n]
  end
  if (nmax == smn.shape[1]-1) then # smn[true,nmax+1] がない場合
    sval[0..nmax,nmax,false] = -(nmax-1)*smn[0..nmax,nmax-1,false]*@emn[0..nmax,nmax]
  else # smn[true,nmax+1] がある場合
#    binding.pry
    sval[0..nmax,nmax,false] = (nmax+2)*smn[0..nmax,nmax+1,false]*@emn[0..nmax,nmax+1]-(nmax-1)*smn[0..nmax,nmax-1,false]*@emn[0..nmax,nmax]
  end
  sval[0..nmax,nmax+1,false] = -(nmax)*smn[0..nmax,nmax,false]*@emn[0..nmax,nmax+1]

  sval[-1..-nmax,true,false] = sval[1..nmax,true,false].conj


  # GPhysオブジェクト化
	# 軸の設定
	old_grid = gphys_smn.grid_copy #入力GPhysオブジェクトのグリッド情報を取得
	yaxis = Axis.new
	twn = VArray.new(NArray.int(nmax+2).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"twn")
	yaxis.pos = twn
	new_grid = old_grid.change_axis(1, yaxis) #第0, 1軸を置き換える

	# 変数の設定
	sname = gphys_smn.name
	lname = gphys_smn.long_name
	unit = gphys_smn.units.to_s

	va = VArray.new(sval, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_smn_new = GPhys.new(new_grid, va)
  return gphys_smn_new
end

# スペクトルデータから、d(g*cos^2(lat))/dmu  (mu = sin(lat)) のスペクトルを求める。
def sp_dmu_cos2(gphys_smn)
  smn = gphys_smn.val

  imax = smn.shape[0]
  nmax = @pmn[0,0,true].size-2 ## smn.shape[1] - 1 # 切断波数
  if smn.rank == 2 then
    sval = NArray.complex(imax,nmax+1)
  else
    shape_ary = [imax,nmax+1] + smn.shape[2..-1]
    sval = NArray.complex(*shape_ary)
  end

  sval[0..nmax,0,false] = 0.0
  for n in 1..nmax-1
    sval[0..nmax,n,false] = -(n+1)*smn[0..nmax,n-1,false]*@emn[0..nmax,n] + n*smn[0..nmax,n+1,false]*@emn[0..nmax,n+1]
  end
  if (nmax == smn.shape[1]-1) then # smn[true,nmax+1] がない場合
    sval[0..nmax,nmax,false] = -(nmax+1)*smn[0..nmax,nmax-1,false]*@emn[0..nmax,nmax]
  else # smn[true,nmax+1] がある場合
#    binding.pry
    sval[0..nmax,nmax,false] = -(nmax+1)*smn[0..nmax,nmax-1,false]*@emn[0..nmax,nmax] + nmax*smn[0..nmax,nmax+1,false]*@emn[0..nmax,nmax+1]
  end
#  sval[0..nmax,nmax+1,false] = -(nmax+2)*smn[0..nmax,nmax,false]*@emn[0..nmax,nmax+1]

  sval[-1..-nmax,true,false] = sval[1..nmax,true,false].conj


  # GPhysオブジェクト化
	# 軸の設定
	old_grid = gphys_smn.grid_copy #入力GPhysオブジェクトのグリッド情報を取得
	yaxis = Axis.new
	twn = VArray.new(NArray.int(nmax+1).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"twn")
	yaxis.pos = twn
	new_grid = old_grid.change_axis(1, yaxis) #第0, 1軸を置き換える

	# 変数の設定
	sname = gphys_smn.name
	lname = gphys_smn.long_name
	unit = gphys_smn.units.to_s

	va = VArray.new(sval, {"long_name"=>lname, "units"=>unit}, sname)
	gphys_smn_new = GPhys.new(new_grid, va)
  return gphys_smn_new
end



# 経度微分d/dlonしたものと緯度微分d/dmしたものを同時に出力（d^2/dmudlon ではない）
# 入力が [x,y,z] のとき、出力は [dxdlon, dydlon, dzdlon,dxdlat, dydlat, dzdlat] となる
def sh_dlon_and_dlat(gphys)
  if (gphys.class == GPhys) then
  	#グリッド→スペクトル
  	smn = sh_trans(gphys)
  	#スペクトル→グリッド [緯度微分：d/dmu, 緯度微分：d/dmu]
  	return [sh_invtrans(smn, true, false), sh_invtrans(smn, false, true)]
  elsif (gphys.class == Array) then
    smn = sh_trans_multi(gphys)
    return [sh_invtrans_multi(smn, true, false), sh_invtrans_multi(smn, false, true)].flatten
  else
    raise "arg must be GPhys or Array of GPhyses which are same shape)"
  end

end



# 緯度微分d/dmu（グリッド→スペクトルー(緯度微分)→グリッド）
# 引数：GPhysオブジェクト、P0n、ガウス重み
def lg_dlat(gphys)
	#グリッド→スペクトル
	sn = lg_trans( gphys)
	#スペクトル→グリッド（緯度微分：d/dmu）
	gresult = lg_invtrans(sn, true)
	return gresult
end


def sh_laplacian(gphys)
  nmax = @pmn[0,0,true].size-2
  if (gphys.class == GPhys) then
    s = sh_trans(gphys)
    (nmax+1).times{|n|
      s[true,n,false] = -s[true,n,false]*n*(n+1)
    }
#    s[true,-1,false] = 0.0 # test
    return sh_invtrans(s)
  elsif (gphys.class == Array) then
    s = sh_trans_multi(gphys)
    s.length.times{|i|
      (nmax+1).times{|n|
        s[i][true,n,false] = -s[i][true,n,false]*n*(n+1)
      }
    }
    return sh_invtrans_multi(s)
  else
    raise "arg must be GPhys or Array of GPhyses which are same shape)"
  end
end

def sp_laplacian(s)
  sout = s.copy
  nmax = @pmn[0,0,true].size-2
  (nmax+1).times{|n|
    sout[true,n,false] = -s[true,n,false]*n*(n+1)
  }
  return sout
end

def sh_invlaplacian(gphys)
  nmax = @pmn[0,0,true].size-2
  if (gphys.class == GPhys) then
    s = sh_trans(gphys)
    (1..nmax).each{|n|
      s[true,n,false] = -s[true,n,false]/(n*(n+1))
    }
    return sh_invtrans(s)
  elsif (gphys.class == Array) then
    s = sh_trans_multi(gphys)
    s.length.times{|i|
      (1..nmax).each{|n|
        s[i][true,n,false] = -s[i][true,n,false]/(n*(n+1))
      }
    }
    return sh_invtrans_multi(s)
  else
    raise "arg must be GPhys or Array of GPhyses which are same shape)"
  end
end

def sp_invlaplacian(s)
  nmax = @pmn[0,0,true].size-2
  sout = s.copy
  (1..nmax).each{|n|
    sout[true,n,false] = -s[true,n,false]/(n*(n+1))
  }
  return sout
end



# u, v から鉛直渦度を計算（グリッド → グリッド）
def sh_vorticity(u, v, a=1.0, only=false)
	# u にcoslatをかける
	ucoslat = multiply_cos(u)
	vcoslat = multiply_cos(v)
	#グリッド→スペクトル
	su = sh_trans( ucoslat, false, only)
	sv = sh_trans( vcoslat, false, only)

	#スペクトル→グリッド（経度微分：d/dlon）
	v_dlon = sh_invtrans(sv, 1, false, only)
	#スペクトル→グリッド（緯度微分：d/dmu）
	u_dlat = sh_invtrans(su, false, 1, only)

	gresult = multiply_cos(v_dlon, -2)/a - u_dlat/a
	gresult.rename("Vor")
	gresult.set_att("long_name" , "Vorticity")
	gresult.set_att("units" , "s-1")
	return gresult
end

# u, v から水平発散を計算（グリッド → グリッド）
def sh_divergence(u, v, a=1.0, only=false)
	# v にcoslatをかける
	vcoslat = multiply_cos(v, 1)
	ucoslat = multiply_cos(u, 1)

  #微分をどのタイミングで行うかによって、結果が微妙に変わるのだが、なぜ？

	#グリッド→スペクトル
	 su = sh_trans( ucoslat, false, only)
	 sv = sh_trans( vcoslat, false, only)
	# #スペクトル→グリッド（経度微分：d/dlon）
	 u_dlon = sh_invtrans(su, 1, false, only)
	# #スペクトル→グリッド（緯度微分：d/dmu）
	 v_dlat = sh_invtrans(sv, false, 1, only)

	gresult = multiply_cos(u_dlon, -2)/a + v_dlat/a

	gresult.rename("Div")
	gresult.set_att("long_name", "Divergence")
	gresult.set_att("units" ,"s-1")

	return gresult
end

def sh_vor_and_div(u, v, a=1.0, only=false)
	# u にcoslatをかける
	ucoslat = multiply_cos(u)
	vcoslat = multiply_cos(v)

	#グリッド→スペクトル
  su, sv = sh_trans_multi([ucoslat, vcoslat], false, false, only)
	#スペクトル→グリッド（経度微分：d/dlon）
  u_dlon, v_dlon = sh_invtrans_multi([su,sv],1,false,only)
	#スペクトル→グリッド（緯度微分：d/dmu）
  u_dlat, v_dlat = sh_invtrans_multi([su,sv],false,1,only)

	gvor = multiply_cos(v_dlon, -2)/a - u_dlat/a
	gvor.rename("Vor")
	gvor.set_att("long_name" , "Vorticity")
	gvor.set_att("units" , "s-1")

  gdiv = multiply_cos(u_dlon, -2)/a + v_dlat/a
	gdiv.rename("Div")
	gdiv.set_att("long_name", "Divergence")
	gdiv.set_att("units" ,"s-1")

	return [gvor, gdiv]

end

# u, v から鉛直渦度を計算（グリッド → スペクトル）
def sp_vorticity(u, v, a=1.0, only=false)
  # 切断波数
  nmax=@pmn[0,0,true].size-2
	# u, v にcos^-1(lat)をかける
	ucoslat_1 = multiply_cos(u, -1)
	vcoslat_1 = multiply_cos(v, -1)
	#グリッド→スペクトル
	su = sh_trans( ucoslat_1,  false, only, nmax+1)
	sv = sh_trans( vcoslat_1, false, only)
	#スペクトル→スペクトル（経度微分：d/dlon）
	sv_dlon = sp_dlon(sv)
	#スペクトル→スペクトル（緯度微分：d(g*cos^2(lat))/dmu）
	su_dlat = sp_dmu_cos2(su)

  sresult = (sv_dlon - su_dlat)/a
  sresult.rename("Vor")
  sresult.set_att("long_name", "Vorticity")
  sresult.set_att("units" , "s-1")

	# return sh_invtrans(sresult)
  return sresult
end

# u, v から水平発散を計算（グリッド → スペクトル）
def sp_divergence(u, v, a=1.0, only=false)
  # 切断波数
  nmax=@pmn[0,0,true].size-2
	# u, v にcos^-1(lat)をかける
	ucoslat_1 = multiply_cos(u, -1)
	vcoslat_1 = multiply_cos(v, -1)
	#グリッド→スペクトル
	su = sh_trans( ucoslat_1,  false, only)
	sv = sh_trans( vcoslat_1, false, only, nmax+1)
	#スペクトル→スペクトル（経度微分：d/dlon）
	su_dlon = sp_dlon(su)
	#スペクトル→スペクトル（緯度微分：d(g*cos^2(lat))/dmu）
	sv_dlat = sp_dmu_cos2(sv)

  sresult = (su_dlon + sv_dlat)/a
  sresult.rename("Div")
  sresult.set_att("long_name", "Divergence")
  sresult.set_att("units" , "s-1")

	# return sh_invtrans(sresult,false,false,false,nmax)
  return sresult
end


# 水平風(u,v)の回転成分・発散成分を求める
def sh_uvcomps(u, v, a=1.0, comp_ary=["u_rot"])
	nmax = @pmn[0,0,true].size - 2

  # 渦度・発散のスペクトル → 逆ラプラシアン（渦度→流線関数、発散→速度ポテンシャル）
  if (comp_ary.include?("u_rot") or comp_ary.include?("v_rot")) then 
    psi_s = sp_invlaplacian( sp_vorticity(u,v,a) )
  end
  if (comp_ary.include?("u_div") or comp_ary.include?("v_div")) then 
    chi_s = sp_invlaplacian( sp_divergence(u,v,a) )
  end

  # 水平風の回転成分・発散成分
  gphys_es_ary = []
  comp_ary.each{|comp|
    case comp
    when "u_rot"
      gp = -sh_invtrans(sp_cos2_dmu(psi_s))/(@coslat)*a
      gp.rename("u_rot")
      gp.set_att("long_name","u (rot. comp.)")
    when "u_div"
      gp = sh_invtrans(sp_dlon(chi_s))/(@coslat)*a
      gp.rename("u_div")
      gp.set_att("long_name","u (div. comp.)")
    when "v_rot"
      gp = sh_invtrans(sp_dlon(psi_s))/(@coslat)*a
      gp.rename("v_rot")
      gp.set_att("long_name","v (rot. comp.)")
    when "v_div"
      gp = sh_invtrans(sp_cos2_dmu(chi_s))/(@coslat)*a
      gp.rename("v_div")
      gp.set_att("long_name","v (div. comp.)")
    end
    gphys_es_ary << gp
  }
  return gphys_es_ary
end 


# NOTE
# 係数 0.25, 0.125は２倍しなければならないはず。
#
def sh_energyspectrum(u, v, a=1.0, comp_ary=["total"], vor_div_given=false)
# 参照：Burgess, Erler, and Shepherd (2013) The Troposphere-to-Stratosphere Transition in Kinetic Energy Spectra and Nonlinear Spectral Fluxes as Seen in ECMWF Analyses, JAS, vol. 70, pp. 669-687
	nmax = @pmn[0,0,true].size - 2

  if (vor_div_given == false) then
    print "Note: U & V must be given.\n"
    flag_vor = false; flag_div = false
    comp_ary.each{|comp|
      if comp.include?("total")
        flag_vor = true; flag_div = true
      elsif comp.include?("rot")
        flag_vor = true
      elsif comp.include?("div")
        flag_div = true
      else
        raise "the 4th argument must be Array of 'total' (default), 'total_u', 'total_v', 'rot', 'rot_u', 'rot_v', 'div', 'div_u', and/or 'div_v'."
      end
    }
    if flag_vor then
      vor = sh_vorticity( u, v, a )
      print "Vorticity caluculated.\n"
  		svor = sh_trans(vor)
      print "Vorticity trasformed.\n"
  		shapes = svor.shape
  		svorval = svor.val
    end
    if flag_div then
      div = sh_divergence( u, v, a )
      print "Divergence calculated.\n"
  		sdiv = sh_trans(div)
      print "Divergence transformed.\n"
  		shapes = sdiv.shape
  		sdivval = sdiv.val
    end
  elsif (vor_div_given == true) then
    print "Note: Vorticity & Divergence must be given.\n"
    vor = u; div = v
    svor = sh_trans(vor)
    print "Vorticity transformed.\n"
    sdiv = sh_trans(div)
    print "Divergence transformed.\n"
    shapes = sdiv.shape
    sdivval = sdiv.val
    svorval = svor.val
  else
    raise "the 5th argument must be true (vor & div are given) or false (u & v are given; default)"
  end


  gphys_es_ary = []
  comp_ary.each{|comp|

    case @precision
    when "double" then
    	if (shapes.size == 2) then
    	  es = NArrayMiss.float(nmax+1).fill!(0.0).all_invalid
    	elsif (shapes.size == 3) then
    	  es = NArrayMiss.float(nmax+1, shapes[2]).fill!(0.0).all_invalid
    	elsif (shapes.size == 4) then
    	  es = NArrayMiss.float(nmax+1, shapes[2], shapes[3]).fill!(0.0).all_invalid
    	else
      		raise "number of input grid data dimension must be less than 5"
    	end
    when "single" then
      if (shapes.size == 2) then
    	  es = NArrayMiss.sfloat(nmax+1).fill!(0.0).all_invalid
    	elsif (shapes.size == 3) then
    	  es = NArrayMiss.sfloat(nmax+1, shapes[2]).fill!(0.0).all_invalid
    	elsif (shapes.size == 4) then
    	  es = NArrayMiss.sfloat(nmax+1, shapes[2], shapes[3]).fill!(0.0).all_invalid
    	else
      		raise "number of input grid data dimension must be less than 5"
    	end
    end

    es_rank = es.rank

  	for n in 1..nmax
  		es[n, false] = 0.0
  		case comp
  		when "rot" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (svorval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false] if es_rank == 1
  		when "rot_u" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (2.0*n*(n+1) - (2.0*n+1)*(m.abs))*(svorval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1
  		when "rot_v" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (m.abs)*(svorval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)*(2.0*n+1.0)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)*(2.0*n+1.0)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1

  		when "div" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (sdivval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false] if es_rank == 1
  		when "div_u" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (m.abs)*(sdivval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)*(2.0*n+1.0)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)*(2.0*n+1.0)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1
  		when "div_v" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (2.0*n*(n+1) - (2.0*n+1)*(m.abs))*(sdivval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1

  		when "total" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (svorval[m, n, false].abs)**2 + (sdivval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.5*(a**2)/(n*(n+1.0))*es[n, false] if es_rank == 1
  		when "total_u" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (2.0*n*(n+1) - (2.0*n+1)*(m.abs))*(svorval[m, n, false].abs)**2 \
  				                                       + (2.0*n+1.0)*(m.abs)*(sdivval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1
  		when "total_v" then
  			for m in -n..n
  				es[n,false] = es[n, false]  + (2.0*n*(n+1) - (2.0*n+1)*(m.abs))*(sdivval[m, n, false].abs)**2 \
  				                                       + (2.0*n+1.0)*(m.abs)*(svorval[m, n, false].abs)**2
  			end
  			es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false].to_na if es_rank > 1
        es[n, false] = 0.25*(a**2)/(n*n*(n+1.0)*(n+1.0))*es[n, false] if es_rank == 1
  		end
  	end

  	#GPhysオブジェクト化
  	case comp
  	when "rot", "total", "rot_u", "total_u", "rot_v", "total_v" then
  		old_grid = vor.grid_copy
  	when "div", "div_u", "div_v" then
  		old_grid = div.grid_copy
  	end

  	yaxis = Axis.new
  	twn = VArray.new(NArray.int(nmax+1).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"twn")
  	yaxis.pos = twn
  	new_grid = old_grid.change_axis(1, yaxis).delete_axes(0) #第1軸を置き換える、第0軸を消去する

  	case comp
  	when "rot" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (rot)", "units"=>"m2 s-2"}, "EngySpectRot")
  	when "rot_u" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (rot_u)", "units"=>"m2 s-2"}, "EngySpectRot_u")
  	when "rot_v" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (rot_v)", "units"=>"m2 s-2"}, "EngySpectRot_v")

  	when "div" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (div)", "units"=>"m2 s-2"}, "EngySpectDiv")
  	when "div_u" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (div_u)", "units"=>"m2 s-2"}, "EngySpectDiv_u")
  	when "div_v" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (div_v)", "units"=>"m2 s-2"}, "EngySpectDiv_v")

  	when "total" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (total)", "units"=>"m2 s-2"}, "EngySpectTotal")
  	when "total_u" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (total_u)", "units"=>"m2 s-2"}, "EngySpectTotal_u")
  	when "total_v" then
  	   va = VArray.new(es, {"long_name"=>"Energy Spectrum (total_v)", "units"=>"m2 s-2"}, "EngySpectTotal_v")
  	end
  	gphys_es= GPhys.new(new_grid, va)
    gphys_es_ary << gphys_es
  }

	return gphys_es_ary
end

def output_pmn(lon, lat = lon.size/2, mmax = ((lat.size*2-1)/3), gp=nil, precision=nil)
  print "Truncation wavenumber is #{mmax}.\n"
   # 倍精度/単精度の選択
  if (precision) then
    @precision = precision
  elsif (gp) then
    case gp.ntype
    when "float" then
      @precision = "double"
    when "sfloat" then
      @precision = "single"
    end
  else
    @precision = "double" #デフォルトの精度
  end

  mu = deg2mu(lat)

  pmn, dpmn, emn = p_m_nT(NArray.to_na(mu), mmax) #三角切断 P^m_n, dP^m_n
  pmn.reshape!(mmax+2,mmax+2,1)
  dpmn.reshape!(mmax+2,mmax+2,1)

  maxis = Axis.new
  m_va = VArray.new(NArray.int(mmax+2).indgen!,{"long_name"=>"zonal wave number","units"=>"wave number"},"m")
  maxis.pos = m_va
  naxis = Axis.new
  n_va = VArray.new(NArray.int(mmax+2).indgen!,{"long_name"=>"total wave number","units"=>"wave number"},"n")
  naxis.pos = n_va
  jaxis = Axis.new
  mu_va = VArray.new(NArray.to_na(mu),{"long_name"=>"sin(lat)","units"=>"1"},"mu")
  jaxis.pos = mu_va
  grid = Grid.new(maxis, naxis, jaxis)
  gpmn = GPhys.new(grid,VArray.new(pmn,{"long_name"=>"assoc legendre polynominal","units"=>"1"},"pmn"))
  gdpmn = GPhys.new(grid,VArray.new(pmn,{"long_name"=>"assoc legendre polynominal deriv","units"=>"1"},"dpmn"))

  return gpmn, gdpmn


end



end
