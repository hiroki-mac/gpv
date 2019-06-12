#--
# lagrangian methods for gpv.
#++
# reference: https://github.com/tenomoto/advection made by Dr. Takeshi Enomoto

require "narray"
require "gsl"

include NMath
module LagrangeSphere
  module_function
  # lon, lat は radian 単位、latは南極から詰める
  # 入力は NArray or Float
  def grid_ext(glon, glat, n=1, keep_index=true) # grid を両側にnずつ拡張する。ただしindexはそのまま。
    lonlen = glon.length
    latlen = glat.length
    glon_ext = NArray.sfloat(lonlen+2*n)
    glat_ext = NArray.sfloat(latlen+2*n)

    glon_ext[0..(lonlen-1)] = glon
    glat_ext[0..(latlen-1)] = glat

    # 等間隔格子を仮定
    dlon = glon[1]-glon[0]
    dlat = glat[1]-glat[0]

    n.times{|m|
      glon_ext[-(m+1)] = glon[0] - dlon*(m+1)
      glon_ext[lonlen+m] = glon[-1] + dlon*(m+1)
      glat_ext[-(m+1)] = glat[0] - dlat*(m+1)
      glat_ext[latlen+m] = glat[-1] + dlat*(m+1)
    }

    if (keep_index) then
      return glon_ext, glat_ext
    else
      glon_ext_o = NArray.sfloat(lonlen+2*n)
      glat_ext_o = NArray.sfloat(latlen+2*n)
      glon_ext_o[n..-1] = glon_ext[0..(-1-n)]
      glon_ext_o[0..(n-1)] = glon_ext[(-n)..(-1)]
      glat_ext_o[n..-1] = glat_ext[0..(-1-n)]
      glat_ext_o[0..(n-1)] = glat_ext[(-n)..(-1)]
      return glon_ext_o, glat_ext_o
    end
  end


  def data_ext(data, n=1, keep_index=true)
    shape = data.shape
    if shape[2] then
      data_ext = NArray.sfloat(shape[0]+2*n, shape[1]+2*n, shape[2])
    else
      data_ext = NArray.sfloat(shape[0]+2*n, shape[1]+2*n)
    end

    data_ext[0..(shape[0]-1), 0..(shape[1]-1),false] = data
    n.times{|m|
      # 経度方向
      data_ext[-(m+1), 0..(shape[1]-1),false] = data[-(m+1),true,false]
      data_ext[shape[0]+m, 0..(shape[1]-1),false] = data[0+m,true,false]
      # 緯度方向 ↓ 以下の拡張は、東西半球が変わることを考慮していないので、不適では？
      data_ext[0..(shape[0]-1), -(m+1),false] = data[true, m+1,false]
      data_ext[0..(shape[0]-1), shape[1]+m,false] = data[true, shape[1]-(m+2),false]
    }
    if (keep_index) then
      return data_ext
    else
      if (shape[2]) then
        data_ext_o = NArray.sfloat(shape[0]+2*n, shape[1]+2*n, shape[2])
      else
        data_ext_o = NArray.sfloat(shape[0]+2*n, shape[1]+2*n)
      end
      data_ext_o[n..-1, n..-1,false] = data_ext[0..(-1-n), 0..(-1-n),false]
      data_ext_o[0..(n-1), n..-1,false] = data_ext[(-n)..(-1), 0..(-1-n),false]
      data_ext_o[n..-1, 0..(n-1),false] = data_ext[0..(-1-n), (-n)..(-1),false]
      return data_ext_o
    end
  end

  def data_ext_halo_only(gp, n=1, vect_comp=true)
    shape = gp.shape
    if (shape[3]) then
      de_x0 = NArray.sfloat(2*n-1, shape[1], shape[2], shape[3]) # 左端
      de_x1 = NArray.sfloat(2*n-1, shape[1], shape[2], shape[3]) # 右端
      de_y0 = NArray.sfloat(shape[0], 2*n-1, shape[2], shape[3]) # 下端
      de_y1 = NArray.sfloat(shape[0], 2*n-1, shape[2], shape[3]) # 上端
    elsif (shape[2]) then
      de_x0 = NArray.sfloat(2*n-1, shape[1], shape[2]) # 左端
      de_x1 = NArray.sfloat(2*n-1, shape[1], shape[2]) # 右端
      de_y0 = NArray.sfloat(shape[0], 2*n-1, shape[2]) # 下端
      de_y1 = NArray.sfloat(shape[0], 2*n-1, shape[2]) # 上端
    else
      de_x0 = NArray.sfloat(2*n-1, shape[1]) # 左端
      de_x1 = NArray.sfloat(2*n-1, shape[1]) # 右端
      de_y0 = NArray.sfloat(shape[0], 2*n-1) # 下端
      de_y1 = NArray.sfloat(shape[0], 2*n-1) # 上端
    end
    (2*n-1).times{|m|
      # 経度方向
      de_x0[m, true,false] = gp[-n+m,true,false].val.to_na
      de_x1[m, true,false] = gp[(shape[0]-3+m)%shape[0], true,false].val.to_na
      # 緯度方向
      if (m < n) then # 東西半球の反転
        de_y0[true,m,false] = gp[true, n-m, false].val.to_na
        tmp = de_y0[0..(shape[0]/2-1),m,false]
        de_y0[0..(shape[0]/2-1),m,false] = -de_y0[(shape[0]/2)..-1,m,false]
        de_y0[(shape[0]/2)..-1,m,false] = -tmp

        de_y1[true,2*n-2-m,false] = gp[true,shape[1]-1-n+m, false].val.to_na
        de_y1[0..(shape[0]/2-1),m,false] = -de_y0[(shape[0]/2)..-1,m,false]
        de_y1[(shape[0]/2)..-1,m,false] = -tmp
      else
        de_y0[true,m,false] = gp[true, m-n, false].val.to_na
        de_y1[true,2*n-2-m,false] = gp[true,shape[1]-1-(m-n),false].val.to_na
      end
    }
    return [de_x0, de_x1, de_y0, de_y1]
  end



  def lonlat2xyz(lon, lat, r=1.0) # 緯度経度座標からデカルト座標に変換
    if (lon.class == NArray) then
      x = NMath.cos(lon)*NMath.cos(lat)*r
      y = NMath.sin(lon)*NMath.cos(lat)*r
      z = NMath.sin(lat)*r
    else
      x = cos(lon)*cos(lat)*r
      y = sin(lon)*cos(lat)*r
      z = sin(lat)*r
    end
    return [x, y, z]
  end


  def uv2xyz(u, v, lon, lat) # 緯度経度座標の流速 u, v からデカルト座標の流速 xd, yd, zd に変換
    if (u.class == NArray) then
      xd = -u*NMath.sin(lon) - v*NMath.cos(lon)*NMath.sin(lat)
      yd = u*NMath.cos(lon) - v*NMath.sin(lon)*NMath.sin(lat)
      zd = v*NMath.cos(lat)
    else
      xd = -u*sin(lon) - v*cos(lon)*sin(lat)
      yd = u*cos(lon) - v*sin(lon)*sin(lat)
      zd = v*cos(lat)
    end
    return [xd, yd, zd]
  end


  def find_midpoint_xyz(x, y, z, xd, yd, zd, dt, r=1.0) #中間点を推定する
    tmp = r*r + (0.5*dt)*(0.5*dt)*(xd*xd+yd*yd+zd*zd) - 2.0*(0.5*dt)*(x*xd+y*yd+z*zd)
    b = r/(tmp**0.5) # 球面上に束縛するための係数
    xm = b*(x - 0.5*dt*xd) # x.sbt!(xd*dt*0.5).mul!(b)
    ym = b*(y - 0.5*dt*yd) # y.sbt!(yd*dt*0.5).mul!(b)
    zm = b*(z - 0.5*dt*zd) # z.sbt!(zd*dt*0.5).mul!(b)
    return [xm, ym, zm] # [x, y, z]
  end


  def find_departurepoint_xyz(x, y, z, xm, ym, zm, r=1.0) #始点を求める
    b = 2.0*(x*xm + y*ym + z*zm)/(r*r) # 球面上に束縛するための係数
    x0 = b*xm - x # xm.mul!(b).sbt!(x)
    y0 = b*ym - y # ym.mul!(b).sbt!(y)
    z0 = b*zm - z # zm.mul!(b).sbt!(z)
    return [x0, y0, z0] # [xm, ym, zm]
  end


  def xyz2lonlat(x, y, z, r=1.0) # デカルト座標から緯度経度座標に変換（0:2PI）
    if (x.class == NArray) then
      lat = NMath.asin(z/r)
      lon = NMath.atan2(y/r,x/r) +PI*2.0 % PI*2.0
      lon = lon + lon.lt(0).to_type("sfloat")*PI*2.0
    else
      lat = asin(z/r)
      lon = atan2(y/r,x/r)+PI*2.0 % PI*2.0
      lon = lon + PI*2.0 if lon < 0
    end
    return [lon, lat]
  end


  def find_stencil(plon, plat) # 指定座標を囲む4つのグリッド点を見つける
    # plon, plat が指定座標
    # 指定座標を囲むグリッドのインデックスは i_na, i_na + 1, j_na, j_na + 1 になる。
    if (plon.class == NArray) then # 1次元 NArray
      len = plon.length
      i_na = NArray.int(len).fill(-999999)
      j_na = NArray.int(len).fill(-999999)
      len.times{|n|
        i_na[n] = GSL::Interp.bsearch(@gslv_lon_ext, plon[n], 0, @gslv_lon_ext.length-1)
        j_na[n] = GSL::Interp.bsearch(@gslv_lat_ext, plat[n], 0, @gslv_lat_ext.length-1)
      }
      return i_na, j_na
    else
      i_na = GSL::Interp.bsearch(@gslv_lon_ext, plon, 0, @gslv_lon_ext.length-1)
      j_na = GSL::Interp.bsearch(@gslv_lat_ext, plat, 0, @gslv_lat_ext.length-1)
      return i_na, j_na
    end
  end

  def find_stencil_1D(pz) # 指定座標を囲む4つのグリッド点を見つける
    # pz が指定座標
    # 指定座標を囲むグリッドのインデックスは k_na, k_na + 1になる。
    if (pz.class == NArray) then # 1次元 NArray
      len = pz.length
      k_na = NArray.int(len).fill(-999999)
      len.times{|n|
        k_na[n] = GSL::Interp.bsearch(@gslv_z, pz[n], 0, @gslv_z.length-1)
      }
      return k_na
    else
      k_na = GSL::Interp.bsearch(@gslv_z, pz, 0, @gslv_z.length-1)
      return k_na
    end
  end



  def interp1D(xa, ya, xt, type="linear")
    ip = GSL::Interp.alloc(type, xa, ya) # 初期化
    return ip.eval(xa, ya, xt)
  end


  def interp2D(glon, glat, data, plon, plat, i, j, type="linear") # 2次元補完（１点のみ）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    case type
    when "linear", "cspline"
      if (1 <= i && i <= glon.length-3) then
        slon = glon[(i-1)..(i+2)]
        da0 = data[(i-1)..(i+2), j-1]
        da1 = data[(i-1)..(i+2), j  ]
        da2 = data[(i-1)..(i+2), j+1]
        da3 = data[(i-1)..(i+2), j+2]
      else
        slon = NArray[ glon[i-1], glon[i], glon[i+1], glon[i+2] ]
        da0 = NArray[ data[i-1,j-1], data[i,j-1], data[i+1,j-1], data[i+2,j-1] ]
        da1 = NArray[ data[i-1,j  ], data[i,j  ], data[i+1,j  ], data[i+2,j  ] ]
        da2 = NArray[ data[i-1,j+1], data[i,j+1], data[i+1,j+1], data[i+2,j+1] ]
        da3 = NArray[ data[i-1,j+2], data[i,j+2], data[i+1,j+2], data[i+2,j+2] ]
      end
      if (1 <= j && j <= glat.length-3) then
       slat = glat[(j-1)..(j+2)]
      else
       slat = NArray[ glat[j-1], glat[j], glat[j+1], glat[j+2] ]
      end

      #puts slon.to_a, slat.to_a, plon, plat

      ip = GSL::Interp.alloc(type, slon, da0) # 初期化
      d0 = ip.eval(slon, da0, plon)
      d1 = ip.eval(slon, da1, plon)
      d2 = ip.eval(slon, da2, plon)
      d3 = ip.eval(slon, da3, plon)
      pval = interp1D(slat, NArray[d0,d1,d2,d3], plat, type)

    when "nearest-neighbor"
      dlon0 = (plon - glon[i]).abs; dlon1 = (glon[i+1] - plon).abs
      dlat0 = (plat - glat[j]).abs; dlat1 = (glat[j+1] - plat).abs
      if (dlon0 < dlon1) then it = i else it = i + 1 end
      if (dlat0 < dlat1) then jt = j else jt = j + 1 end
      pval = data[it,jt]

    else
      raise "#{type} interpolation is not supported."
    end

    return pval
  end


  def interp2D_2(u, v, plon, plat, i, j, type="linear") # 2次元補完（１点のみ, 2変数同時）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"
    if (1 <= i && i <= @gslv_lon_ext.length-3) then
      slon = @gslv_lon_ext[(i-1)..(i+2)]
      ua0 =  GSL::Vector[ u[(i-1)..(i+2), j-1].to_a]
      ua1 =  GSL::Vector[ u[(i-1)..(i+2), j  ].to_a]
      ua2 =  GSL::Vector[ u[(i-1)..(i+2), j+1].to_a]
      ua3 =  GSL::Vector[ u[(i-1)..(i+2), j+2].to_a]
      va0 =  GSL::Vector[ v[(i-1)..(i+2), j-1].to_a]
      va1 =  GSL::Vector[ v[(i-1)..(i+2), j  ].to_a]
      va2 =  GSL::Vector[ v[(i-1)..(i+2), j+1].to_a]
      va3 =  GSL::Vector[ v[(i-1)..(i+2), j+2].to_a]

    else
      slon = GSL::Vector[ @gslv_lon_ext[i-1], @gslv_lon_ext[i], @gslv_lon_ext[i+1], @gslv_lon_ext[i+2] ]
      ua0 = GSL::Vector[ u[i-1,j-1], u[i,j-1], u[i+1,j-1], u[i+2,j-1] ]
      ua1 = GSL::Vector[ u[i-1,j  ], u[i,j  ], u[i+1,j  ], u[i+2,j  ] ]
      ua2 = GSL::Vector[ u[i-1,j+1], u[i,j+1], u[i+1,j+1], u[i+2,j+1] ]
      ua3 = GSL::Vector[ u[i-1,j+2], u[i,j+2], u[i+1,j+2], u[i+2,j+2] ]
      va0 = GSL::Vector[ v[i-1,j-1], v[i,j-1], v[i+1,j-1], v[i+2,j-1] ]
      va1 = GSL::Vector[ v[i-1,j  ], v[i,j  ], v[i+1,j  ], v[i+2,j  ] ]
      va2 = GSL::Vector[ v[i-1,j+1], v[i,j+1], v[i+1,j+1], v[i+2,j+1] ]
      va3 = GSL::Vector[ v[i-1,j+2], v[i,j+2], v[i+1,j+2], v[i+2,j+2] ]

    end
    if (1 <= j && j <= @gslv_lat_ext.length-3) then
     slat = @gslv_lat_ext[(j-1)..(j+2)]
    else
     slat = GSL::Vector[ @gslv_lat_ext[j-1], @gslv_lat_ext[j], @gslv_lat_ext[j+1], @gslv_lat_ext[j+2] ]
    end

    #puts slon.to_a, slat.to_a, plon, plat

    ip = GSL::Interp.alloc(type, slon, ua0) # 初期化
    u0 = ip.eval(slon, ua0, plon)
    u1 = ip.eval(slon, ua1, plon)
    u2 = ip.eval(slon, ua2, plon)
    u3 = ip.eval(slon, ua3, plon)
    v0 = ip.eval(slon, va0, plon)
    v1 = ip.eval(slon, va1, plon)
    v2 = ip.eval(slon, va2, plon)
    v3 = ip.eval(slon, va3, plon)

    uaa = GSL::Vector[u0,u1,u2,u3]
    vaa = GSL::Vector[v0,v1,v2,v3]
    ip = GSL::Interp.alloc(type, slat, uaa) # 初期化
    uval = ip.eval(slat, uaa, plat)
    vval = ip.eval(slat, vaa, plat)
    # else
    # end
    return [uval, vval]
  end


  def interp2D_2all(glon, glat, u, v, plon, plat, type="linear") # 2次元補完（2変数同時）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"

    # 準備
    rlonlen = plon.shape[0]
    rlatlen = plon.shape[1]

    uval = NArray.sfloat(rlonlen, rlatlen)
    vval = NArray.sfloat(rlonlen, rlatlen)

    ip = GSL::Interp.alloc(type, glon, u[true,0]) # 初期化
    u_i_lon = NArray.sfloat(rlonlen, rlatlen, glat.length)
    v_i_lon = NArray.sfloat(rlonlen, rlatlen, glat.length)
    glat.length.times{|j|
      u_i_lon[true,true, j] = ip.eval(glon,u[true,j],GSL::Vector.to_gv(plon)).to_na.reshape(rlonlen, rlatlen)
      v_i_lon[true,true, j] = ip.eval(glon,v[true,j],GSL::Vector.to_gv(plon)).to_na.reshape(rlonlen, rlatlen)
    }

    ip = GSL::Interp.alloc(type, glat, u[0,true]) # 初期化
    rlatlen.times{|n|
      rlonlen.times{|m|
       uval[m,n] = ip.eval(glat,u_i_lon[m,n,true],plat[m,n])
       vval[m,n] = ip.eval(glat,v_i_lon[m,n,true],plat[m,n])
      }
    }
    return uval, vval
  end

  def interp3D(u, v, plon, plat, pz, i, j, k, type="linear") # 3次元補完（１点のみ, 2変数同時）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"
    sz = @gslv_z[k..(k+1)]
    # 鉛直内挿の係数
    c0 = (sz[1]-pz); c1 = (pz-sz[0]); c2 = 1.0/(sz[1] - sz[0])
    if (1 <= i && i <= @gslv_lon_ext.length-3) then
      slon = @gslv_lon_ext[(i-1)..(i+2)]
      # 鉛直内挿して、GSL::Vectorにする。
      ua0 =  GSL::Vector[ ((u[(i-1)..(i+2),j-1,k]*c0+u[(i-1)..(i+2),j-1,k+1]*c1)*c2).to_a]
      ua1 =  GSL::Vector[ ((u[(i-1)..(i+2),j  ,k]*c0+u[(i-1)..(i+2),j  ,k+1]*c1)*c2).to_a]
      ua2 =  GSL::Vector[ ((u[(i-1)..(i+2),j+1,k]*c0+u[(i-1)..(i+2),j+1,k+1]*c1)*c2).to_a]
      ua3 =  GSL::Vector[ ((u[(i-1)..(i+2),j+2,k]*c0+u[(i-1)..(i+2),j+2,k+1]*c1)*c2).to_a]
      va0 =  GSL::Vector[ ((v[(i-1)..(i+2),j-1,k]*c0+v[(i-1)..(i+2),j-1,k+1]*c1)*c2).to_a]
      va1 =  GSL::Vector[ ((v[(i-1)..(i+2),j  ,k]*c0+v[(i-1)..(i+2),j  ,k+1]*c1)*c2).to_a]
      va2 =  GSL::Vector[ ((v[(i-1)..(i+2),j+1,k]*c0+v[(i-1)..(i+2),j+1,k+1]*c1)*c2).to_a]
      va3 =  GSL::Vector[ ((v[(i-1)..(i+2),j+2,k]*c0+v[(i-1)..(i+2),j+2,k+1]*c1)*c2).to_a]
    else
      slon = GSL::Vector[ @gslv_lon_ext[i-1], @gslv_lon_ext[i], @gslv_lon_ext[i+1], @gslv_lon_ext[i+2] ]
      ua0 = (GSL::Vector[u[i-1,j-1,k  ], u[i,j-1,k  ], u[i+1,j-1,k  ], u[i+2,j-1,k  ]]*c0 +
             GSL::Vector[u[i-1,j-1,k+1], u[i,j-1,k+1], u[i+1,j-1,k+1], u[i+2,j-1,k+1]]*c1   )*c2
      ua1 = (GSL::Vector[u[i-1,j  ,k  ], u[i,j  ,k  ], u[i+1,j  ,k  ], u[i+2,j  ,k  ]]*c0 +
             GSL::Vector[u[i-1,j  ,k+1], u[i,j  ,k+1], u[i+1,j  ,k+1], u[i+2,j  ,k+1]]*c1   )*c2
      ua2 = (GSL::Vector[u[i-1,j+1,k  ], u[i,j+1,k  ], u[i+1,j+1,k  ], u[i+2,j+1,k  ]]*c0 +
             GSL::Vector[u[i-1,j+1,k+1], u[i,j+1,k+1], u[i+1,j+1,k+1], u[i+2,j+1,k+1]]*c1   )*c2
      ua3 = (GSL::Vector[u[i-1,j+2,k  ], u[i,j+2,k  ], u[i+1,j+2,k  ], u[i+2,j+2,k  ]]*c0 +
             GSL::Vector[u[i-1,j+2,k+1], u[i,j+2,k+1], u[i+1,j+2,k+1], u[i+2,j+2,k+1]]*c1   )*c2
      va0 = (GSL::Vector[v[i-1,j-1,k  ], v[i,j-1,k  ], v[i+1,j-1,k  ], v[i+2,j-1,k  ]]*c0 +
             GSL::Vector[v[i-1,j-1,k+1], v[i,j-1,k+1], v[i+1,j-1,k+1], v[i+2,j-1,k+1]]*c1   )*c2
      va1 = (GSL::Vector[v[i-1,j  ,k  ], v[i,j  ,k  ], v[i+1,j  ,k  ], v[i+2,j  ,k  ]]*c0 +
             GSL::Vector[v[i-1,j  ,k+1], v[i,j  ,k+1], v[i+1,j  ,k+1], v[i+2,j  ,k+1]]*c1   )*c2
      va2 = (GSL::Vector[v[i-1,j+1,k  ], v[i,j+1,k  ], v[i+1,j+1,k  ], v[i+2,j+1,k  ]]*c0 +
             GSL::Vector[v[i-1,j+1,k+1], v[i,j+1,k+1], v[i+1,j+1,k+1], v[i+2,j+1,k+1]]*c1   )*c2
      va3 = (GSL::Vector[v[i-1,j+2,k  ], v[i,j+2,k  ], v[i+1,j+2,k  ], v[i+2,j+2,k  ]]*c0 +
             GSL::Vector[v[i-1,j+2,k+1], v[i,j+2,k+1], v[i+1,j+2,k+1], v[i+2,j+2,k+1]]*c1   )*c2
    end
    if (1 <= j && j <= @gslv_lat_ext.length-3) then
     slat = @gslv_lat_ext[(j-1)..(j+2)]
    else
     slat = GSL::Vector[ @gslv_lat_ext[j-1], @gslv_lat_ext[j], @gslv_lat_ext[j+1], @gslv_lat_ext[j+2] ]
    end

    #水平補完
    ip = GSL::Interp.alloc(type, slon, ua0) # 初期化
    u0 = ip.eval(slon, ua0, plon)
    u1 = ip.eval(slon, ua1, plon)
    u2 = ip.eval(slon, ua2, plon)
    u3 = ip.eval(slon, ua3, plon)
    v0 = ip.eval(slon, va0, plon)
    v1 = ip.eval(slon, va1, plon)
    v2 = ip.eval(slon, va2, plon)
    v3 = ip.eval(slon, va3, plon)

    uaa = GSL::Vector[u0,u1,u2,u3]
    vaa = GSL::Vector[v0,v1,v2,v3]
    ip = GSL::Interp.alloc(type, slat, uaa) # 初期化
    uval = ip.eval(slat, uaa, plat)
    vval = ip.eval(slat, vaa, plat)
    # else
    # end
    return [uval, vval]
  end

  def interp4D_gp(gu, gv, plon, plat, pz, i, j, k, t, t_div, n, type="linear") # 3次元補完（１点のみ, 2変数同時）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"
    sz = @gslv_z[k..(k+1)]
    # 鉛直内挿の係数
    c0 = (sz[1]-pz); c1 = (pz-sz[0]); c2 = 1.0/(sz[1] - sz[0])

    if (1 <= i && i <= gu.shape[0]-3 && 1 <= j && j <= gu.shape[1]-3) then # indexがはみ出さない場合
      # 必要なデータの切り出し、NArray化
      u = gu[(i-1)..(i+2),(j-1)..(j+2),k..(k+1),t..(t+1)].val
      v = gv[(i-1)..(i+2),(j-1)..(j+2),k..(k+1),t..(t+1)].val
      # 時間内挿
      u = (u[true,true,true,0]*(t_div-n)+u[true,true,true,1]*n)/t_div
      v = (v[true,true,true,0]*(t_div-n)+v[true,true,true,1]*n)/t_div
      # 鉛直内挿
      u = (u[true,true,0]*c0+u[true,true,1]*c1)*c2
      v = (v[true,true,0]*c0+v[true,true,1]*c1)*c2
      # 水平補完の準備　
      slon = @gslv_lon_ext[(i-1)..(i+2)]
      ua0 =  GSL::Vector[ u[true, 0].to_a]; ua1 =  GSL::Vector[ u[true, 1].to_a]; ua2 =  GSL::Vector[ u[true, 2].to_a]; ua3 =  GSL::Vector[ u[true, 3].to_a]
      va0 =  GSL::Vector[ v[true, 0].to_a]; va1 =  GSL::Vector[ v[true, 1].to_a]; va2 =  GSL::Vector[ v[true, 2].to_a]; va3 =  GSL::Vector[ v[true, 3].to_a]
    else # indexがはみ出す場合
      @u_halo = data_ext_halo_only(gu, 4, true) unless @u_halo
      @v_halo = data_ext_halo_only(gv, 4, true) unless @v_halo
      if (i < 1 && 1 <= j && j <= gu.shape[1]-3) then # 左側
        u = @u_halo[0][(4+i-1)..(4+i+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
        v = @v_halo[0][(4+i-1)..(4+i+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
      elsif (gu.shape[0]-3 < i && 1 <= j && j <= gu.shape[1]-3) then # 右側
        u = @u_halo[1][(i-gu.shape[0]+3-1)..(i-gu.shape[0]+3+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
        v = @v_halo[1][(i-gu.shape[0]+3-1)..(i-gu.shape[0]+3+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
      elsif (1 <= i && i <= gu.shape[0]-3 && j < 1 ) then # 下側
        u = @u_halo[2][(i-1)..(i+2), (4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        v = @v_halo[2][(i-1)..(i+2), (4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
      elsif (1 <= i && i <= gu.shape[0]-3 && gu.shape[1]-3 < j) then # 上側
        u = @u_halo[3][(i-1)..(i+2),(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2),k..(k+1),t..(t+1)]
        v = @v_halo[3][(i-1)..(i+2),(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2),k..(k+1),t..(t+1)]
      elsif (j < 1 ) then# 下の隅
        u = NArray.sfloat(4,4,2,2)
        u[0,false] = @u_halo[2][(i-1)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        u[1,false] = @u_halo[2][(i  )%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        u[2,false] = @u_halo[2][(i+1)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        u[3,false] = @u_halo[2][(i+2)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        v = NArray.sfloat(4,4,2,2)
        v[0,false] = @v_halo[2][(i-1)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        v[1,false] = @v_halo[2][(i  )%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        v[2,false] = @v_halo[2][(i+1)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        v[3,false] = @v_halo[2][(i+2)%gu.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
      elsif (gu.shape[1]-3 < j) then # 上の隅
        u = NArray.sfloat(4,4,2,2)
        u[0,false] = @u_halo[3][(i-1)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        u[1,false] = @u_halo[3][(i  )%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        u[2,false] = @u_halo[3][(i+1)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        u[3,false] = @u_halo[3][(i+2)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        v = NArray.sfloat(4,4,2,2)
        v[0,false] = @v_halo[3][(i-1)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        v[1,false] = @v_halo[3][(i  )%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        v[2,false] = @v_halo[3][(i+1)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
        v[3,false] = @v_halo[3][(i+2)%gu.shape[0],(j-gu.shape[1]+3-1)..(j-gu.shape[1]+3+2), k..(k+1),t..(t+1)]
      else
        binding.pry # 当てはまらないはず
      end

      # 時間内挿
      u = (u[true,true,true,0]*(t_div-n)+u[true,true,true,1]*n)/t_div
      v = (v[true,true,true,0]*(t_div-n)+v[true,true,true,1]*n)/t_div
      # 鉛直内挿
      u = (u[true,true,0]*c0+u[true,true,1]*c1)*c2
      v = (v[true,true,0]*c0+v[true,true,1]*c1)*c2
      # 水平補完の準備
      slon = GSL::Vector[ @gslv_lon_ext[i-1], @gslv_lon_ext[i], @gslv_lon_ext[i+1], @gslv_lon_ext[i+2] ]
      ua0 =  GSL::Vector[ u[true, 0].to_a]; ua1 =  GSL::Vector[ u[true, 1].to_a]; ua2 =  GSL::Vector[ u[true, 2].to_a]; ua3 =  GSL::Vector[ u[true, 3].to_a]
      va0 =  GSL::Vector[ v[true, 0].to_a]; va1 =  GSL::Vector[ v[true, 1].to_a]; va2 =  GSL::Vector[ v[true, 2].to_a]; va3 =  GSL::Vector[ v[true, 3].to_a]


      #
      #
      # p i, j
      # p "I0"
      # imax = gu.shape[0] ; ia = [(i-1)%imax, i%imax, (i+1)%imax, (i+2)%imax] # 東西は周期境界
      # jmax = gu.shape[1]-1; ja = [j-1, j, j+1, j+2]
      # a = NArray.sfloat(4,4).fill(1.0) # ベクトル成分を反転するかどうかの係数
      # 4.times{|m|
      #   if (ja[m] > jmax) then
      #     ja[m] = 2*jmax-ja[m]; a[true,m] = -1.0; ia[m] = (ia[m]+imax/2)%imax # 東西半球を反転させる
      #   elsif (ja[m] < 0) then
      #     ja[m] = -ja[m]      ; a[true,m] = -1.0; ia[m] = (ia[m]+imax/2)%imax # 東西半球を反転させる
      #   end
      # }
      # p "I1"
      # # 必要なデータの切り出し、NArray化
      # gu = gu.val
      # uk0t0 = NArray.to_na(
      #       [ [gu[ia[0],ja[0],k,t], gu[ia[1],ja[0],k,t], gu[ia[2],ja[0],k,t], gu[ia[3],ja[0],k,t]],
      #         [gu[ia[0],ja[1],k,t], gu[ia[1],ja[1],k,t], gu[ia[2],ja[1],k,t], gu[ia[3],ja[1],k,t]],
      #         [gu[ia[0],ja[2],k,t], gu[ia[1],ja[2],k,t], gu[ia[2],ja[2],k,t], gu[ia[3],ja[2],k,t]],
      #         [gu[ia[0],ja[3],k,t], gu[ia[1],ja[3],k,t], gu[ia[2],ja[3],k,t], gu[ia[3],ja[3],k,t]] ] )
      # uk1t0 = NArray.to_na(
      #       [ [gu[ia[0],ja[0],k+1,t], gu[ia[1],ja[0],k+1,t], gu[ia[2],ja[0],k+1,t], gu[ia[3],ja[0],k+1,t]],
      #         [gu[ia[0],ja[1],k+1,t], gu[ia[1],ja[1],k+1,t], gu[ia[2],ja[1],k+1,t], gu[ia[3],ja[1],k+1,t]],
      #         [gu[ia[0],ja[2],k+1,t], gu[ia[1],ja[2],k+1,t], gu[ia[2],ja[2],k+1,t], gu[ia[3],ja[2],k+1,t]],
      #         [gu[ia[0],ja[3],k+1,t], gu[ia[1],ja[3],k+1,t], gu[ia[2],ja[3],k+1,t], gu[ia[3],ja[3],k+1,t]] ] )
      # uk0t1 = NArray.to_na(
      #       [ [gu[ia[0],ja[0],k,t+1], gu[ia[1],ja[0],k,t+1], gu[ia[2],ja[0],k,t+1], gu[ia[3],ja[0],k,t+1]],
      #         [gu[ia[0],ja[1],k,t+1], gu[ia[1],ja[1],k,t+1], gu[ia[2],ja[1],k,t+1], gu[ia[3],ja[1],k,t+1]],
      #         [gu[ia[0],ja[2],k,t+1], gu[ia[1],ja[2],k,t+1], gu[ia[2],ja[2],k,t+1], gu[ia[3],ja[2],k,t+1]],
      #         [gu[ia[0],ja[3],k,t+1], gu[ia[1],ja[3],k,t+1], gu[ia[2],ja[3],k,t+1], gu[ia[3],ja[3],k,t+1]] ] )
      # uk1t1 = NArray.to_na(
      #       [ [gu[ia[0],ja[0],k+1,t+1], gu[ia[1],ja[0],k+1,t+1], gu[ia[2],ja[0],k+1,t+1], gu[ia[3],ja[0],k+1,t+1]],
      #         [gu[ia[0],ja[1],k+1,t+1], gu[ia[1],ja[1],k+1,t+1], gu[ia[2],ja[1],k+1,t+1], gu[ia[3],ja[1],k+1,t+1]],
      #         [gu[ia[0],ja[2],k+1,t+1], gu[ia[1],ja[2],k+1,t+1], gu[ia[2],ja[2],k+1,t+1], gu[ia[3],ja[2],k+1,t+1]],
      #         [gu[ia[0],ja[3],k+1,t+1], gu[ia[1],ja[3],k+1,t+1], gu[ia[2],ja[3],k+1,t+1], gu[ia[3],ja[3],k+1,t+1]] ] )
      # uk0t0.map!{|i| i}; uk0t0 = uk0t0.to_type(4)*a
      # uk1t0.map!{|i| i}; uk1t0 = uk1t0.to_type(4)*a
      # uk0t1.map!{|i| i}; uk0t1 = uk0t1.to_type(4)*a
      # uk1t1.map!{|i| i}; uk1t1 = uk1t1.to_type(4)*a
      # p "I2"
      # vk0t0 = NArray.to_na(
      #       [ [gv[ia[0],ja[0],k,t], gv[ia[1],ja[0],k,t], gv[ia[2],ja[0],k,t], gv[ia[3],ja[0],k,t]],
      #         [gv[ia[0],ja[1],k,t], gv[ia[1],ja[1],k,t], gv[ia[2],ja[1],k,t], gv[ia[3],ja[1],k,t]],
      #         [gv[ia[0],ja[2],k,t], gv[ia[1],ja[2],k,t], gv[ia[2],ja[2],k,t], gv[ia[3],ja[2],k,t]],
      #         [gv[ia[0],ja[3],k,t], gv[ia[1],ja[3],k,t], gv[ia[2],ja[3],k,t], gv[ia[3],ja[3],k,t]] ] )
      # vk1t0 = NArray.to_na(
      #       [ [gv[ia[0],ja[0],k+1,t], gv[ia[1],ja[0],k+1,t], gv[ia[2],ja[0],k+1,t], gv[ia[3],ja[0],k+1,t]],
      #         [gv[ia[0],ja[1],k+1,t], gv[ia[1],ja[1],k+1,t], gv[ia[2],ja[1],k+1,t], gv[ia[3],ja[1],k+1,t]],
      #         [gv[ia[0],ja[2],k+1,t], gv[ia[1],ja[2],k+1,t], gv[ia[2],ja[2],k+1,t], gv[ia[3],ja[2],k+1,t]],
      #         [gv[ia[0],ja[3],k+1,t], gv[ia[1],ja[3],k+1,t], gv[ia[2],ja[3],k+1,t], gv[ia[3],ja[3],k+1,t]] ] )
      # vk0t1 = NArray.to_na(
      #       [ [gv[ia[0],ja[0],k,t+1], gv[ia[1],ja[0],k,t+1], gv[ia[2],ja[0],k,t+1], gv[ia[3],ja[0],k,t+1]],
      #         [gv[ia[0],ja[1],k,t+1], gv[ia[1],ja[1],k,t+1], gv[ia[2],ja[1],k,t+1], gv[ia[3],ja[1],k,t+1]],
      #         [gv[ia[0],ja[2],k,t+1], gv[ia[1],ja[2],k,t+1], gv[ia[2],ja[2],k,t+1], gv[ia[3],ja[2],k,t+1]],
      #         [gv[ia[0],ja[3],k,t+1], gv[ia[1],ja[3],k,t+1], gv[ia[2],ja[3],k,t+1], gv[ia[3],ja[3],k,t+1]] ] )
      # vk1t1 = NArray.to_na(
      #       [ [gv[ia[0],ja[0],k+1,t+1], gv[ia[1],ja[0],k+1,t+1], gv[ia[2],ja[0],k+1,t+1], gv[ia[3],ja[0],k+1,t+1]],
      #         [gv[ia[0],ja[1],k+1,t+1], gv[ia[1],ja[1],k+1,t+1], gv[ia[2],ja[1],k+1,t+1], gv[ia[3],ja[1],k+1,t+1]],
      #         [gv[ia[0],ja[2],k+1,t+1], gv[ia[1],ja[2],k+1,t+1], gv[ia[2],ja[2],k+1,t+1], gv[ia[3],ja[2],k+1,t+1]],
      #         [gv[ia[0],ja[3],k+1,t+1], gv[ia[1],ja[3],k+1,t+1], gv[ia[2],ja[3],k+1,t+1], gv[ia[3],ja[3],k+1,t+1]] ] )
      # vk0t0.map!{|i| i.val[0]}; vk0t0 = vk0t0.to_type(4)*a
      # vk1t0.map!{|i| i.val[0]}; vk1t0 = vk1t0.to_type(4)*a
      # vk0t1.map!{|i| i.val[0]}; vk0t1 = vk0t1.to_type(4)*a
      # vk1t1.map!{|i| i.val[0]}; vk1t1 = vk1t1.to_type(4)*a
      # p "I3"
      #
      # # 時間内挿
      # uk0 = (uk0t0*(t_div-n)+uk0t1*n)/t_div; uk1 = (uk1t0*(t_div-n)+uk1t1*n)/t_div
      # vk0 = (vk0t0*(t_div-n)+vk0t1*n)/t_div; vk1 = (vk1t0*(t_div-n)+vk1t1*n)/t_div
      # # 鉛直内挿
      # u = (uk0*c0+uk1*c1)*c2
      # v = (vk0*c0+vk1*c1)*c2
    end

    if (1 <= j && j <= @gslv_lat_ext.length-3) then
     slat = @gslv_lat_ext[(j-1)..(j+2)]
    else
     slat = GSL::Vector[ @gslv_lat_ext[j-1], @gslv_lat_ext[j], @gslv_lat_ext[j+1], @gslv_lat_ext[j+2] ]
    end

    #水平補完
    ip = GSL::Interp.alloc(type, slon, ua0) # 初期化
    u0 = ip.eval(slon, ua0, plon)
    u1 = ip.eval(slon, ua1, plon)
    u2 = ip.eval(slon, ua2, plon)
    u3 = ip.eval(slon, ua3, plon)
    v0 = ip.eval(slon, va0, plon)
    v1 = ip.eval(slon, va1, plon)
    v2 = ip.eval(slon, va2, plon)
    v3 = ip.eval(slon, va3, plon)

    uaa = GSL::Vector[u0,u1,u2,u3]
    vaa = GSL::Vector[v0,v1,v2,v3]
    ip = GSL::Interp.alloc(type, slat, uaa) # 初期化
    uval = ip.eval(slat, uaa, plat)
    vval = ip.eval(slat, vaa, plat)
    # else
    # end
    return [uval, vval]
  end

  def interp3D_1(u, plon, plat, pz, i, j, k, type="linear") # 3次元補完（１点のみ, 1変数）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"
    sz = @gslv_z[k..(k+1)]
    # 鉛直内挿の係数
    c0 = (sz[1]-pz); c1 = (pz-sz[0]); c2 = 1.0/(sz[1] - sz[0])
    if (1 <= i && i <= @gslv_lon_ext.length-3) then
      slon = @gslv_lon_ext[(i-1)..(i+2)]
      # 鉛直内挿して、GSL::Vectorにする。
      ua0 =  GSL::Vector[ ((u[(i-1)..(i+2),j-1,k]*c0+u[(i-1)..(i+2),j-1,k+1]*c1)*c2).to_a]
      ua1 =  GSL::Vector[ ((u[(i-1)..(i+2),j  ,k]*c0+u[(i-1)..(i+2),j  ,k+1]*c1)*c2).to_a]
      ua2 =  GSL::Vector[ ((u[(i-1)..(i+2),j+1,k]*c0+u[(i-1)..(i+2),j+1,k+1]*c1)*c2).to_a]
      ua3 =  GSL::Vector[ ((u[(i-1)..(i+2),j+2,k]*c0+u[(i-1)..(i+2),j+2,k+1]*c1)*c2).to_a]
    else
      slon = GSL::Vector[ @gslv_lon_ext[i-1], @gslv_lon_ext[i], @gslv_lon_ext[i+1], @gslv_lon_ext[i+2] ]
      ua0 = (GSL::Vector[u[i-1,j-1,k  ], u[i,j-1,k  ], u[i+1,j-1,k  ], u[i+2,j-1,k  ]]*c0 +
             GSL::Vector[u[i-1,j-1,k+1], u[i,j-1,k+1], u[i+1,j-1,k+1], u[i+2,j-1,k+1]]*c1   )*c2
      ua1 = (GSL::Vector[u[i-1,j  ,k  ], u[i,j  ,k  ], u[i+1,j  ,k  ], u[i+2,j  ,k  ]]*c0 +
             GSL::Vector[u[i-1,j  ,k+1], u[i,j  ,k+1], u[i+1,j  ,k+1], u[i+2,j  ,k+1]]*c1   )*c2
      ua2 = (GSL::Vector[u[i-1,j+1,k  ], u[i,j+1,k  ], u[i+1,j+1,k  ], u[i+2,j+1,k  ]]*c0 +
             GSL::Vector[u[i-1,j+1,k+1], u[i,j+1,k+1], u[i+1,j+1,k+1], u[i+2,j+1,k+1]]*c1   )*c2
      ua3 = (GSL::Vector[u[i-1,j+2,k  ], u[i,j+2,k  ], u[i+1,j+2,k  ], u[i+2,j+2,k  ]]*c0 +
             GSL::Vector[u[i-1,j+2,k+1], u[i,j+2,k+1], u[i+1,j+2,k+1], u[i+2,j+2,k+1]]*c1   )*c2
    end
    if (1 <= j && j <= @gslv_lat_ext.length-3) then
     slat = @gslv_lat_ext[(j-1)..(j+2)]
    else
     slat = GSL::Vector[ @gslv_lat_ext[j-1], @gslv_lat_ext[j], @gslv_lat_ext[j+1], @gslv_lat_ext[j+2] ]
    end

    #水平補完
    ip = GSL::Interp.alloc(type, slon, ua0) # 初期化
    begin
      u0 = ip.eval(slon, ua0, plon)
      u1 = ip.eval(slon, ua1, plon)
      u2 = ip.eval(slon, ua2, plon)
      u3 = ip.eval(slon, ua3, plon)
    rescue
      binding.pry
    end

    uaa = GSL::Vector[u0,u1,u2,u3]
    ip = GSL::Interp.alloc(type, slat, uaa) # 初期化
    uval = ip.eval(slat, uaa, plat)
    # else
    # end
    return uval
  end

  def interp4D_1_gp(gw, plon, plat, pz, i, j, k, t, t_div, n, type="linear") # 3次元補完（１点のみ, 1変数）
    # グリッド、データ、補完座標、補完座標の左と下のインデックス
    # 内部で方向分離して、1次元補完を呼ぶ。
    # case type
    # when "linear", "cspline"
    sz = @gslv_z[k..(k+1)]
    # 鉛直内挿の係数
    c0 = (sz[1]-pz); c1 = (pz-sz[0]); c2 = 1.0/(sz[1] - sz[0])
    if (1 <= i && i <= gw.shape[0]-3 && 1 <= j && j <= gw.shape[1]-3) then # indexがはみ出さない場合
      # 必要なデータの切り出し、NArray化
      w = gw[(i-1)..(i+2),(j-1)..(j+2),k..(k+1),t..(t+1)].val
      # 時間内挿
      w = (w[true,true,true,0]*(t_div-n)+w[true,true,true,1]*n)/t_div
      # 鉛直内挿
      w = (w[true,true,0]*c0+w[true,true,1]*c1)*c2
      # 水平補完の準備　
      slon = @gslv_lon_ext[(i-1)..(i+2)]
      wa0 =  GSL::Vector[ w[true, 0].to_a]; wa1 =  GSL::Vector[ w[true, 1].to_a]; wa2 =  GSL::Vector[ w[true, 2].to_a]; wa3 =  GSL::Vector[ w[true, 3].to_a]
    else # indexがはみ出す場合
      @w_halo = data_ext_halo_only(gw,4,false) unless @w_halo
      if (i < 1 && 1 <= j && j <= gw.shape[1]-3) then # 左側
        w = @w_halo[0][(4+i-1)..(4+i+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
      elsif (gw.shape[0]-3 < i && 1 <= j && j <= gw.shape[1]-3) then # 右側
        w = @w_halo[1][(i-gw.shape[0]+3-1)..(i-gw.shape[0]+3+2),(j-1)..(j+2),k..(k+1),t..(t+1)]
      elsif (1 <= i && i <= gw.shape[0]-3 && j < 1 ) then # 下側
        w = @w_halo[2][(i-1)..(i+2), (4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
      elsif (1 <= i && i <= gw.shape[0]-3 && gw.shape[1]-3 < j) then # 上側
        w = @w_halo[3][(i-1)..(i+2),(j-gw.shape[1]+3-1)..(j-gw.shape[1]+3+2),k..(k+1),t..(t+1)]
      elsif (j < 1 ) then# 下の隅
        w = NArray.sfloat(4,4,2,2)
        w[0,false] = @w_halo[2][(i-1)%gw.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        w[1,false] = @w_halo[2][(i  )%gw.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        w[2,false] = @w_halo[2][(i+1)%gw.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
        w[3,false] = @w_halo[2][(i+2)%gw.shape[0],(4+j-1)..(4+j+2), k..(k+1),t..(t+1)]
      elsif (gw.shape[1]-3 < j) then # 上の隅
        w = NArray.sfloat(4,4,2,2)
        w[0,false] = @w_halo[3][(i-1)%gw.shape[0],(j-gw.shape[1]+3-1)..(j-gw.shape[1]+3+2), k..(k+1),t..(t+1)]
        w[1,false] = @w_halo[3][(i  )%gw.shape[0],(j-gw.shape[1]+3-1)..(j-gw.shape[1]+3+2), k..(k+1),t..(t+1)]
        w[2,false] = @w_halo[3][(i+1)%gw.shape[0],(j-gw.shape[1]+3-1)..(j-gw.shape[1]+3+2), k..(k+1),t..(t+1)]
        w[3,false] = @w_halo[3][(i+2)%gw.shape[0],(j-gw.shape[1]+3-1)..(j-gw.shape[1]+3+2), k..(k+1),t..(t+1)]
      else
        binding.pry # 当てはまらないはず
      end

      # 時間内挿
      w = (w[true,true,true,0]*(t_div-n)+w[true,true,true,1]*n)/t_div
      # 鉛直内挿
      w = (w[true,true,0]*c0+w[true,true,1]*c1)*c2
      # 水平補完の準備
      slon = GSL::Vector[ @gslv_lon_ext[i-1], @gslv_lon_ext[i], @gslv_lon_ext[i+1], @gslv_lon_ext[i+2] ]
      wa0 =  GSL::Vector[ w[true, 0].to_a]; wa1 =  GSL::Vector[ w[true, 1].to_a]; wa2 =  GSL::Vector[ w[true, 2].to_a]; wa3 =  GSL::Vector[ w[true, 3].to_a]
    end
    if (1 <= j && j <= @gslv_lat_ext.length-3) then
      slat = @gslv_lat_ext[(j-1)..(j+2)]
    else
      slat = GSL::Vector[ @gslv_lat_ext[j-1], @gslv_lat_ext[j], @gslv_lat_ext[j+1], @gslv_lat_ext[j+2] ]
    end

    #水平補完
    ip = GSL::Interp.alloc(type, slon, wa0) # 初期化
    begin
      w0 = ip.eval(slon, wa0, plon)
      w1 = ip.eval(slon, wa1, plon)
      w2 = ip.eval(slon, wa2, plon)
      w3 = ip.eval(slon, wa3, plon)
    rescue
      binding.pry
    end

    waa = GSL::Vector[w0,w1,w2,w3]
    ip = GSL::Interp.alloc(type, slat, waa) # 初期化
    wval = ip.eval(slat, waa, plat)
    # else
    # end
    return wval
  end

  def find_midpoint(lon, lat, u, v, dt, r=1.0) #中間点を推定する
    x, y, z = lonlat2xyz(lon, lat, r)
    xd, yd, zd = uv2xyz(u,v,lon,lat)
    xm, ym, zm = find_midpoint_xyz(x,y,z,xd,yd,zd,dt,r)
    mlon, mlat = xyz2lonlat(xm,ym,zm, r)
    return mlon, mlat
  end


  def find_departurepoint(lon, lat, mlon, mlat, r=1.0) #始点を求める
    x, y, z = lonlat2xyz(lon, lat, r)
    xm, ym, zm = lonlat2xyz(mlon, mlat, r)
    x0, y0, z0 = find_departurepoint_xyz(x, y, z, xm, ym, zm, r)
    lon0, lat0 = xyz2lonlat(x0,y0,z0, r)
    return lon0, lat0
  end


  def find_departurepoint_Hortal2002(lon, lat, u, v, slon, slat, su, sv, dt, r=1.0) # 始点を求める（Hortal 2002）
    x, y, z = lonlat2xyz(lon, lat, r)
    xd, yd, zd = uv2xyz(u,v,lon,lat)
    sxd, syd, szd = uv2xyz(su,sv,slon,slat)
    # 球面上に束縛するために係数
    tmp = ((xd+sxd)**2+(yd+syd)**2+(zd+szd)**2)*(dt**2)*0.25/(r**2) -
          ((xd+sxd)*x+(yd+syd)*y+(zd+szd)*z)*dt/(r**2) + 1.0
    b = 1.0/(tmp**0.5)
    sx = (x - 0.5*dt*(xd + sxd))*b
    sy = (y - 0.5*dt*(yd + syd))*b
    sz = (z - 0.5*dt*(zd + szd))*b
    slon, slat = xyz2lonlat(sx,sy,sz,r)
    return slon, slat
  end

  # 初期化ルーチン
  def init(gp)
    # gslとNArrayの連携が切れているので、配列は GSL::Vector オブジェクトにしておく。
    # グリッドの取り出しと拡張
    lon = gp.coordinate(0).val*D2R
    lat = gp.coordinate(1).val*D2R
    lon_ext, lat_ext = grid_ext(lon, lat, 4)

    @gslv_lon = GSL::Vector[lon.to_a]
    @gslv_lat = GSL::Vector[lat.to_a]
    @gslv_lon_ext = GSL::Vector[lon_ext.to_a]
    @gslv_lat_ext = GSL::Vector[lat_ext.to_a]
  end

  # 初期化ルーチン
  def init_3D(gp)
    # gslとNArrayの連携が切れているので、配列は GSL::Vector オブジェクトにしておく。
    # グリッドの取り出しと拡張
    lon = gp.coordinate(0).val*D2R
    lat = gp.coordinate(1).val*D2R
    z = gp.coordinate(2).val
    lon_ext, lat_ext = grid_ext(lon, lat, 4)

    @gslv_lon = GSL::Vector[lon.to_a]
    @gslv_lat = GSL::Vector[lat.to_a]
    @gslv_lon_ext = GSL::Vector[lon_ext.to_a]
    @gslv_lat_ext = GSL::Vector[lat_ext.to_a]
    @gslv_z = GSL::Vector[z.to_a]
  end


  # 始点グリッド(slon, slat, sz)と水平流速場(u,v)を与えて、
  # 終点グリッド(elon, elat, ez)を返す。u, v は GPhysオブジェクト
  # dt もGPhysオブジェクト（単位つき）で渡す
  def particle_advection_2D(p_lon, p_lat, gu, gv, dt)
    vTYPE = "cspline"
    if dt.class == VArray then# dtを負で与えることで、始点探索を終点探索に変える。
      if dt.units.to_s.include?("since") then
        dt.units=dt.units.to_s.split("since")[0].strip
      end
      dt_in_sec = -(UNumeric[0.0,"s"] + dt).val[0]
    else
      dt_in_sec = -(UNumeric[0.0,"s"] + dt).to_f
    end
    # 始点位置のGPhysをNArray化。単位変換（Deg→Rad）も実施。
    if (p_lon.units.to_s.downcase.include?("deg")) then
      slon = p_lon.val*D2R; slat = p_lat.val*D2R
    else
      slon = p_lon.val    ; slat = p_lat.val
    end
    # 流速のNArray化と拡張
    if (gu.val.class == NArrayMiss) then
      u_ext = data_ext(gu.val.to_na, 4)
      v_ext = data_ext(gv.val.to_na, 4)
    else
      u_ext = data_ext(gu.val, 4)
      v_ext = data_ext(gv.val, 4)
    end

    # 初期推定値
    mlon = slon; mlat = slat
    if (mlon.class == NArray) then
      i_st, j_st = find_stencil(slon, slat)
      mu = NArray.sfloat(slon.length); mv = NArray.sfloat(slon.length)
      slon.length.times{|n|
        mu[n], mv[n] = interp2D_2(u_ext,v_ext,slon[n],slat[n],i_st[n],j_st[n], vTYPE)
      }
    # elsif (true) then # 以下も遅い
    #   mu = NArray.sfloat(slon.length); mv = NArray.sfloat(slon.length)
    #   slon.length.times{|n|
    #     mu[n] = gu.cut(slon[n],slat[n]).val
    #     mv[n] = gv.cut(slon[n],slat[n]).val
    #   }
    # elsif (true) then # 以下は遅い
    #   tmp = NMatrix.to_na(gu.cut(slon,slat).val.to_a) # NMatrix 生成
    #   mu = NArray.to_na( (tmp - tmp.diagonal(0)).sum(1).flatten.to_a ) # 対角成分だけ取り出したい
    #   tmp = NMatrix.to_na(gv.cut(slon,slat).val.to_a) # NMatrix 生成
    #   mv = NArray.to_na( (tmp - tmp.diagonal(0)).sum(1).flatten.to_a ) # 対角成分だけ取り出したい
    else
      mu = gu.cut(slon,slat).val[0]; mv = gv.cut(slon,slat).val[0]
    end

    # Ritchie (1987)「大円上の等角速度運動を仮定」の方法
    10.times{|it|
      # 前回の推定値を保存
      pmlon = mlon; pmlat = mlat
      # 中間点の推定
      mlon, mlat = find_midpoint(slon, slat, mu, mv, dt_in_sec, GAnalysis::Planet.radius.to_f)
      if (it > 0) then
        if (mlon.class == NArray) then
          dellon = NMath.sin(mlon) - NMath.sin(pmlon)
          dellat = NMath.sin(mlat) - NMath.sin(pmlat)
          delta = [dellon.abs.max, dellat.abs.max].max
          delmean = [dellon.abs.mean, dellat.abs.mean]
        else
          dellon = Math.sin(mlon) - Math.sin(pmlon)
          dellat = Math.sin(mlat) - Math.sin(pmlat)
          delta = [dellon.abs, dellat.abs].max
        end
        if (delta < 0.01*D2R) then
          break
        end
      end
      if (it == 9) then
        break
      end
      # ステンシルを求める
      if (true) then
        i_st, j_st = find_stencil(mlon, mlat)
        # 中間点の u, v を求める
        if (mlon.class == NArray) then
          mlon.length.times{|n|
            mu[n], mv[n] = interp2D_2(u_ext,v_ext,mlon[n],mlat[n],i_st[n],j_st[n], vTYPE)
          }
        else
          mu, mv = interp2D_2(u_ext,v_ext,mlon,mlat,i_st,j_st,vTYPE)
        end
      else
        # 中間点の u, v を求める
        # rlatlen.times{|j|
        #   rlonlen.times{|i|
        #     next if (rad_org[i,j] == miss_value) # 欠損値グリッドはスキップ
        #     next if ([dellon[i,j].abs, dellat[i,j].abs].max < EPS)
        #     i_st[i,j] = GSL::Interp.bsearch(gslv_gvlon, mlon[i,j], 0, gvlon_ext.length-1)
        #     j_st[i,j] = GSL::Interp.bsearch(gslv_gvlat, mlat[i,j], 0, gvlat_ext.length-1)
        #     mu[i,j], mv[i,j] = interp2D_2(lon_ext,lat_ext,u_ext,v_ext,mlon[i,j],mlat[i,j],i_st[i,j],j_st[i,j], vTYPE)
        #   }
        # }
      end
    }# イタレーション ここまで

    # 始点を求める
    elon, elat = find_departurepoint(slon, slat, mlon, mlat, GAnalysis::Planet.radius.to_f)
    ######### Ritchie (1987) ここまで ############
    if (p_lon.units.to_s.downcase.include?("deg")) then
      p_lon.replace_val(elon*R2D)
      p_lat.replace_val(elat*R2D)
    else
      p_lon.replace_val(elon)
      p_lat.replace_val(elat)
    end
    t = p_lon.coord(-1)
    p_lon.axis(-1).set_pos(t+dt.val)
    p_lat.axis(-1).set_pos(t+dt.val)
    return p_lon, p_lat
  end

  # 始点グリッド(slon, slat, sz)と水平流速場(u,v)を与えて、
  # 終点グリッド(elon, elat, ez)を返す。u, v は GPhysオブジェクト
  # dt もGPhysオブジェクト（単位つき）で渡す
  # 計算手順：水平と鉛直は方向分離で計算する
  #         sz の高度で中間点(mlon, mlat)、終点(elon,elat)をもとめ、
  #         中間点の鉛直流の値で sz → ezの鉛直変位を計算する。

  def particle_advection_3D(p_lon, p_lat, p_z, gu, gv, gw, dt, t, t_div, n)
    vTYPE = "cspline"
    if dt.class == VArray then# dtを負で与えることで、始点探索を終点探索に変える。
      if dt.units.to_s.include?("since") then
        dt.units=dt.units.to_s.split("since")[0].strip
      end
      dt_in_sec = -(UNumeric[0.0,"s"] + dt).val[0]
    else
      dt_in_sec = -(UNumeric[0.0,"s"] + dt).to_f
    end
    # 始点位置のGPhysをNArray化。単位変換（Deg→Rad）も実施。
    if (p_lon.units.to_s.downcase.include?("deg")) then
      slon = p_lon.val*D2R; slat = p_lat.val*D2R
    else
      slon = p_lon.val    ; slat = p_lat.val
    end
    # 鉛直位置のNArray化
    sz = p_z.val
    # 流速のNArray化と拡張
    # if (gu.val.class == NArrayMiss) then
      # u_ext = data_ext(gu.val.to_na, 4)
      # v_ext = data_ext(gv.val.to_na, 4)
      # w_ext = data_ext(gw.val.to_na, 4)
    # else
      # u_ext = data_ext(gu.val, 4)
      # v_ext = data_ext(gv.val, 4)
      # w_ext = data_ext(gw.val, 4)
    # end
    # @u_halo = data_ext_halo_only(gu, 4, true)
    # @v_halo = data_ext_halo_only(gv, 4, true)
    # @w_halo = data_ext_halo_only(gw, 4, false)
    @u_halo = nil; @v_halo = nil; @w_halo = nil

    # 初期推定値
    mlon = slon; mlat = slat
    if (mlon.class == NArray) then
      i_st, j_st = find_stencil(slon, slat)
      k_st = find_stencil_1D(sz)
      mu = NArray.sfloat(slon.length); mv = NArray.sfloat(slon.length); mw = NArray.sfloat(slon.length)
      slon.length.times{|n|
#        mu[n], mv[n] = interp3D(u_ext,v_ext,slon[n],slat[n],sz[n],i_st[n],j_st[n],k_st[n],vTYPE)
        mu[n], mv[n] = interp4D_gp(gu,gv,slon[n],slat[n],sz[n],i_st[n],j_st[n],k_st[n], t, t_div, n, vTYPE)
      }
    else
      mu = gu.cut(slon,slat).val[0]; mv = gv.cut(slon,slat).val[0]
    end

    # Ritchie (1987)「大円上の等角速度運動を仮定」の方法
    10.times{|it|
      # 前回の推定値を保存
      pmlon = mlon; pmlat = mlat
      # 中間点の推定
      mlon, mlat = find_midpoint(slon, slat, mu, mv, dt_in_sec, GAnalysis::Planet.radius.to_f)
      if (it > 0) then
        if (mlon.class == NArray) then
          dellon = NMath.sin(mlon) - NMath.sin(pmlon)
          dellat = NMath.sin(mlat) - NMath.sin(pmlat)
          delta = [dellon.abs.max, dellat.abs.max].max
          delmean = [dellon.abs.mean, dellat.abs.mean]
        else
          dellon = Math.sin(mlon) - Math.sin(pmlon)
          dellat = Math.sin(mlat) - Math.sin(pmlat)
          delta = [dellon.abs, dellat.abs].max
        end
        if (delta < 0.01*D2R) then
          break
        end
      end
      if (it == 9) then
        break
      end
      # ステンシルを求める
      if (true) then
        i_st, j_st = find_stencil(mlon, mlat)
        # 中間点の u, v を求める
        if (mlon.class == NArray) then
          mlon.length.times{|n|
#            mu[n], mv[n] = interp3D(u_ext,v_ext,mlon[n],mlat[n],sz[n],i_st[n],j_st[n],k_st[n], vTYPE)
            mu[n], mv[n] = interp4D_gp(gu,gv,mlon[n],mlat[n],sz[n],i_st[n],j_st[n],k_st[n], t, t_div, n, vTYPE)
          }
        else
          mu, mv = interp3D(u_ext,v_ext,mlon,mlat,sz[n],i_st,j_st,k_st,vTYPE)
        end
      else
        # 中間点の u, v を求める
        # rlatlen.times{|j|
        #   rlonlen.times{|i|
        #     next if (rad_org[i,j] == miss_value) # 欠損値グリッドはスキップ
        #     next if ([dellon[i,j].abs, dellat[i,j].abs].max < EPS)
        #     i_st[i,j] = GSL::Interp.bsearch(gslv_gvlon, mlon[i,j], 0, gvlon_ext.length-1)
        #     j_st[i,j] = GSL::Interp.bsearch(gslv_gvlat, mlat[i,j], 0, gvlat_ext.length-1)
        #     mu[i,j], mv[i,j] = interp2D_2(lon_ext,lat_ext,u_ext,v_ext,mlon[i,j],mlat[i,j],i_st[i,j],j_st[i,j], vTYPE)
        #   }
        # }
      end
    }# イタレーション ここまで
    # 始点を求める
    elon, elat = find_departurepoint(slon, slat, mlon, mlat, GAnalysis::Planet.radius.to_f)

    # 中間点での w を求める
    i_st, j_st = find_stencil(mlon, mlat)
    if (mlon.class == NArray) then
      mlon.length.times{|n|
#        mw[n] = interp3D_1(w_ext,mlon[n],mlat[n],sz[n],i_st[n],j_st[n],k_st[n], vTYPE)
        mw[n] = interp4D_1_gp(gw,mlon[n],mlat[n],sz[n],i_st[n],j_st[n],k_st[n], t, t_div, n, vTYPE)
      }
    else
      mw = interp3D_1(w_ext,mlon,mlat,sz[n],i_st,j_st,k_st,vTYPE)
    end
    # 中間点の w の値で鉛直移流させる
    sz.add!(mw*(-dt_in_sec))
    sz.mul!(sz.gt(0)) # 高度が負になった場合に ゼロ に補正する。


    ######### Ritchie (1987) ここまで ############
    if (p_lon.units.to_s.downcase.include?("deg")) then
      p_lon.replace_val(elon*R2D)
      p_lat.replace_val(elat*R2D)
    else
      p_lon.replace_val(elon)
      p_lat.replace_val(elat)
    end
    p_z.replace_val(sz)
    t = p_lon.coord(-1)
    p_lon.axis(-1).set_pos(t+dt.val)
    p_lat.axis(-1).set_pos(t+dt.val)
    p_z.axis(-1).set_pos(t+dt.val)
    return p_lon, p_lat, p_z




  end





end
