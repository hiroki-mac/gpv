#--
# visulization methods for gpv.
#++
class GPV
  def tone_full(gphys, newframe=true, options=nil)
    gropn_1_if_not_yet
    if newframe!=true && newframe!=false
      raise ArgumentError, "2nd arg (newframe) must be true or false"
    end
    opts = @@tone_options.interpret(options)

    if (gphys.class == GPhys) then
      xax, yax, closer, gp = data_prep_2D([gphys], newframe, opts)
      fig(xax, yax) if newframe
      DCL.uwsgxa(xax.val) if xax.rank==1
      DCL.uwsgya(yax.val) if yax.rank==1
      val = gp.data.val
      # set colored value range
      zmin = opts['min'] || val.min
      zmax = opts['max'] || val.max
      DCL.uiscrg(zmin,zmax)
      DCL.uiscsq([zmin,zmax],[0,DCL.uifpac(255,255,255)]) if @OPT_nocolor
      # draw full color fig.
      DCL.uipdat(val)
       if newframe
         axes_or_map_and_ttl(gp,opts, xax, yax)
       end
       closer.call if closer
      return [zmin, zmax]

    elsif (gphys.class == Array) then
      gp1 = gphys[0]; gp2 = gphys[1]
      xax, yax, closer, gp2 = data_prep_2D([gp2], newframe, opts)
      xax, yax, closer, gp1 = data_prep_2D([gp1], newframe, opts)
      fig(xax, yax) if newframe
      DCL.uwsgxa(xax.val) if xax.rank==1
      DCL.uwsgya(yax.val) if yax.rank==1
      val1 = gp1.data.val; val2 = gp2.data.val
      # set colored value range
      zmin1 = opts['min'][0] || val1.min
      zmax1 = opts['max'][0] || val1.max
      zmin2 = opts['min'][1] || val2.min
      zmax2 = opts['max'][1] || val2.max

      if (gphys.length >= 3) then
        gp3 = gphys[2]
        xax, yax, closer, gp3 = data_prep_2D([gp3], newframe, opts)
        val3 = gp3.data.val
        zmin3 = opts['min'][2] || val3.min
        zmax3 = opts['max'][2] || val3.max
      end


      if (false) then
        DCL.uipset("lcycle",true)
        # HSV to RGB translation: c1 左上、c2 右上、c3 左下、c4 右下
        c1 = [180,100,65]
        c2 = [60,100,65]
        c3 = [240,100,35]
        c4 = [0,100,35]
        colors = [c1, c2, c3, c4]
        colors.map!{|hsl|
          r = Color::HSL.new(hsl[0],hsl[1],hsl[2]).to_rgb.red
          b = Color::HSL.new(hsl[0],hsl[1],hsl[2]).to_rgb.blue
          g = Color::HSL.new(hsl[0],hsl[1],hsl[2]).to_rgb.green
          [r,g,b]
        }
        colors.map!{|rgb| DCL.uifpac(rgb[0],rgb[1],rgb[2]) }
        DCL.uiscmp(colors[0],colors[1],colors[2],colors[3])
        DCL.uiscr2(zmin1,zmax1,zmin2,zmax2)
        DCL.uipda2(val1, val2)
      end

      # set hue and values
      sang = 320; wang = 280 # start and width of hue circle [sang = 320, wang = 350]
      valwidth = 220.0; valpower = 1.5
      # 0 <= h < 360.0, 0.0 <= s <= 1.0, 0 <= v <= 255
      nval1= (val1 - zmin1)/(zmax1-zmin1)
        lm = nval1.lt(0); rm = nval1.ge(1.0)
        nval1 = nval1 - lm*nval1 - rm*nval1 + rm*1.0
        nval1 = ((1.0-nval1)*wang + sang).to_i%360
  #    nval1= (sang + 360*(1 - (val1 - zmin1)/(zmax1-zmin1))).to_i%360
      nval2 = (((val2-zmin2)/(zmax2-zmin2))**valpower  )*valwidth + (255.0-valwidth)
        lm = nval2.lt(255.0-valwidth); rm = nval2.gt(255.0)
        nval2 = nval2 - lm*nval2 - rm*nval2 + lm*(255.0-valwidth) + rm*255.0
      if (gphys.length >= 3) then
        nval3 = ((val3 - zmin3)/(zmax3-zmin3))
          lm = nval3.lt(0); rm = nval3.ge(1.0)
          nval3 = nval3 - lm*nval3 - rm*nval3 + rm*1.0
      end
      c = NArray.sfloat(nval1.shape[0],nval1.shape[1]).fill(0.0)
      if (gphys.length == 2 ) then
        r,g,b = DCL.clsvrg(nval1, c+1.0, nval2) # このr,g,bは0-255で出てくる。
      else
        r,g,b = DCL.clsvrg(nval1, nval3, nval2) # このr,g,bは0-255で出てくる。
      end
      r = (c + r)/255.0; g = (c + g)/255.0;  b = (c + b)/255.0
      # draw full color fig.
      DCL.uipda3(r, g, b) # ここのr,g,bは0.0-1.0に正規化したもの。

      axes_or_map_and_ttl(gp1,opts, xax, yax) if newframe
      closer.call if closer

      # colorbar-----------
      unless (@OPT_nocolorbar) then
        vx0, vx1, vy0, vy1 = DCL.sgqvpt # viewportの記憶
        DCL.uzfact(0.6)
        DCL.uzpset("LABELXB",false);  DCL.uzpset("LABELYL",false)
        DCL.uiscr2(zmin1,zmax1,zmin2,zmax2)
        aspect = (vx1 - vx0) / (vy1 - vy0)
        if (aspect > 1.5)
          DCL.uipcmp(vx0,vx0+0.25,vy0-0.15,vy0-0.05,"B")
        else
          DCL.uipcmp(vx1+0.01,vx1+0.21,vy0,vy0+0.2,"B")
        end

        nm = 20
        hue = NArray.sfloat(nm,nm); val = NArray.sfloat(nm,nm)
        sat = NArray.sfloat(nm,nm).fill(0.0)
        nm.times{|n| hue[n,true] = (sang+wang*(1-(n.to_f/nm))).to_i%360
                     val[true,n] = ((n.to_f/nm)**valpower)*valwidth + (255.0-valwidth) }
        r,g,b = DCL.clsvrg(hue, sat+1.0, val)
        r = (sat + r)/255.0; g = (sat + g)/255.0;  b = (sat + b)/255.0
        DCL.uwsgxb(zmin1,zmax1,nm)
        DCL.uwsgyb(zmin2,zmax2,nm)
        DCL.uipda3(r, g, b)

    #    DCL.uxptmk("B",2, [240,260,280,300])
        DCL.uzpset("LABELXB",true)
        if (aspect >= 1.5) then
          DCL.uzpset("LABELYL",true)
        else
          DCL.uzpset("LABELYR",true)
        end

        DCL.uzinit
        DCL.usdaxs
        DCL.uxsttl("B", gp1.name+" ["+gp1.units.to_s+"]", 0)
        if (aspect > 2.0) then
          DCL.uysttl("L", gp2.name+" ["+gp2.units.to_s+"]", 0)
        else
          DCL.uysttl("R", gp2.name+" ["+gp2.units.to_s+"]", 0)
        end
        DCL.sgsvpt(vx0, vx1, vy0, vy1) # viewportを戻す
        DCL.uzfact(1.0/0.6)
      end
      #-----------------
    end
  end

  def draw_setup

    # set missing value
    DCLExt.gl_set_params('lmiss'=>true)#, 'rmiss'=>-999)

    # fontsize
    DCL.sgpset('lcntl', true)
    DCL.sgpset('lfull', true)               # use full area in the window
    DCL.sgpset('lfprop',true)               # use proportional font
    DCL.uscset('cyspos', 'B' )              # move unit y axis
    DCL.sgiset('ifont', 2)
    DCL.uzfact(0.6* (@OPT_fact|| 1).to_f ) unless $flag_uzfact
    $flag_uzfact = true unless $flag_uzfact

  #  DCL.udpset('icycle', 5)
   #  DCL.udpset('label', false) #ラベルの有無
   #  DCL.udpset('lmsg', false) #等値線間隔表示の有無
   # DCL.udpset('indxmj', 1) #計曲線
   # DCL.udpset('indxmn', 1) #主曲線


    # 線の色と太さの設定
    if (@OPT_miscindex) then
      # マージン
      DCL.sgpset('index', @OPT_miscindex.to_i)
      # 等値線
      DCL.udpset('indxmj', @OPT_miscindex.to_i) #計曲線
      DCL.udpset('indxmn', @OPT_miscindex.to_i) #主曲線
      DCL.udpset('label', false) #ラベルの有無
      DCL.udpset('ldash', false) #負を破線にするかどうか

      # ベクトル図
      DCL.ugpset('index', @OPT_miscindex.to_i) #ベクトル
      DCL.ugpset('iuindx', @OPT_miscindex.to_i) #単位ベクトルのタイトル
      # 折れ線図
      DCL.sgpset('indexc', @OPT_miscindex.to_i) #ラベル
      # 直交直線軸
      DCL.uzpset('indext1',@OPT_miscindex.to_i) # 小さい座標軸の目盛
      DCL.uzpset('indext2',@OPT_miscindex.to_i) # 大きい座標軸の目盛
      DCL.uzpset('indexl1',@OPT_miscindex.to_i) # 小さいラベル、タイトル
      DCL.uzpset('indexl2',@OPT_miscindex.to_i) # 大きいラベル、タイトル
    #  DCL.uzpset('rsizec2', 0.01)
      # 地図
      DCL.umpset('indexmj', @OPT_miscindex.to_i) #計曲線
      DCL.umpset('indexmn', @OPT_miscindex.to_i) #主曲線
      DCL.umpset('indexbnd', @OPT_miscindex.to_i) #地図のふち
      DCL.umpset('indexout', @OPT_miscindex.to_i) #海岸線
      DCL.umpset('lgridmj',true) # 緯経線(major)
      DCL.umpset('lgridmn',true) # 緯経線(minor)

    end


    # viewport size
    GGraph.set_fig('viewport'=>@VIEWPORT)

    GGraph.set_fig( 'itr'=>(@OPT_itr == nil) ? 1 : @OPT_itr.to_i )
    GGraph.set_fig("xrev"=>"units:mb,units:hPa,units:millibar,positive:down",
                   "yrev"=>"units:mb,units:hPa,units:millibar,positive:down")

    # viewport layout
    if (@OPT_sldiv)
      sldiv_ary = @OPT_sldiv.split(",")
      DCL.sldiv(sldiv_ary[0],sldiv_ary[1].to_i,sldiv_ary[2].to_i)
    end

    # margin size
    DCL.slmgn(0,0,0.05,0) unless @OPT_panelfit


    # set options
    min_range,  max_range  = __split_range(@OPT_range)
    min_crange, max_crange = __split_range(@OPT_crange)
    min_srange, max_srange = __split_range(@OPT_srange)
    min_irange, max_irange = __split_range(@OPT_irange)

    min_clr, max_clr = __split_range(@OPT_clr_range)

    if (@OPT_clr_range) # set color range
      DCL.uepset("icolor1", min_clr); DCL.uepset("icolor2",max_clr)
    end


    GGraph.set_linear_contour_options(
                                      'int' => ( @OPT_cint   || @OPT_interval || @OPT_int ).to_f,
                                      'min' => ( min_crange  || min_range ),
                                      'max' => ( max_crange  || max_range )
                                      )
    GGraph.set_linear_tone_options(
                                      'int' => ( @OPT_sint   || @OPT_interval || @OPT_int ).to_f,
                                      'min' => ( min_srange  || min_irange || min_range),
                                      'max' => ( max_srange  || max_irange || max_range),
                                      'inf_max' => @OPT_irange, 'inf_min'=> @OPT_irange,
                                      'lbound'=>(@OPT_srange || @OPT_irange || @OPT_range) # if false, min and/or max may be shifted automatically.
                                   )
    if ( @OPT_clevels || @OPT_levels )
      @OPT_label=(@OPT_clevels || @OPT_levels).split(',')
  #    @OPT_label=["10|-10\"","","10|-8\"","","10|-6\"","","10|-4\"","","10|-2\"","","10|0\"","","10|2\"","","10|4\""]
      @OPT_clevels=(@OPT_clevels || @OPT_levels).split(',').map!{|v| v.to_f }
    end

    if ( @OPT_slevels || @OPT_levels )
      @OPT_slevels=(@OPT_slevels||@OPT_levels).split(',').map!{|v| v.to_f }
    end

    if ( @OPT_patterns )
      @OPT_patterns=@OPT_patterns.split(',').map!{|v| v.to_f }
    end


    # set tone levels and pattern with +/- infinity => # not used since gphys 1.5.0.1
    # if (@OPT_irange) then
    #   min_range,  max_range  = __split_range(@OPT_irange)
    #   if (@OPT_int) then
    #     if (@OPT_int.to_f > 0) then
    #       dx = @OPT_int.to_f
    #       nlev = ((max_range - min_range)/dx).to_i
    #     elsif (@OPT_int.to_i < 0) then
    #       nlev = -@OPT_int.to_i
    #       dx = (max_range - min_range)/nlev
    #     end
    #   else
    #     nlev = 10 # default number of patterns within given range
    #     dx = (max_range - min_range)/nlev
    #   end
    #
    #   levels = Array.new(nlev+1)
    #   (nlev+1).times{|i| levels[i] = min_range + dx*i }
    #
    #   dx = (DCL.uepget("icolor2").to_f-DCL.uepget("icolor1").to_f)/(nlev+1)
    #   patterns = Array.new(nlev+2)
    #   (nlev+2).times{|i| patterns[i] = (DCL.uepget("icolor1").to_i + (dx*i).to_i)*1000 + 999 }
    #   @OPT_slevels = levels
    #   @OPT_patterns = patterns
    #   @OPT_clevels = levels unless (@OPT_clevels or @OPT_cint)
    # end




    # similar projection
    if (@OPT_similar)
      if /([\d\-.]*),([\d\-.]*),([\d\-.]*)/ =~ @OPT_similar
        similar=[$1.to_f,$2.to_f,$3.to_f]
      elsif /([\d\-.]*),([\d\-.]*)/ =~ @OPT_similar
        similar=[$1.to_f,$2.to_f,0]
      elsif /([\d\-.]*)/ =~ @OPT_similar
        similar=[$1.to_f,0,0]
      end
      GGraph.set_fig('similar'=>similar)
    end

    # similar projection
    if (@OPT_map_axis)
      if /([\d\-.]*),([\d\-.]*),([\d\-.]*)/ =~ @OPT_map_axis
        @map_axis=[$1.to_f,$2.to_f,$3.to_f]
      elsif /([\d\-.]*),([\d\-.]*)/ =~ @OPT_map_axis
        @map_axis=[$1.to_f,$2.to_f,0]
      elsif /([\d\-.]*)/ =~ @OPT_similar
        @map_axis=[$1.to_f,0,0]
      end
      GGraph.set_fig('map_axis'=>@map_axis)
    end

    # clipping parameter
    if ( @OPT_map_radius ) then
      map_radius=@OPT_map_radius.to_f
    else
      map_radius=90.0
    end
    GGraph.set_fig('map_radius'=>map_radius)

    # simulate satelite view
    if (@OPT_sateliteview) then
      distance = @OPT_sateliteview.to_f
      clip_deg = acos(1.0/distance)*R2D
      GGraph.set_fig('map_rsat'=>distance,'map_radius'=>clip_deg)
    end



    # map
    if ( @OPT_m || @OPT_map)
      map_type = "coast_world"     if @OPT_m
      map_type = @OPT_map          if @OPT_map
      GGraph::set_map(map_type=>true)
      DCL.umpset("lgridmj",false)
      DCL.umpset("lgridmn",false)
      DCL.umpset("indexmj",1)
      DCL.umpset("itypemn",1)

      DCL.umpmap(map_type)
    end


    # time axis
    if (@OPT_time_ax)
      @OPT_time_ax = false if @OPT_time_ax == "false"
      GGraph.set_axes('time_ax'=>@OPT_time_ax)
    end

    if (@OPT_xint)
      xint, xint_div = @OPT_xint.split(":"); xint_div = 1 if xint_div == nil
      GGraph.set_axes('xtickint'=>(xint.to_f / xint_div.to_i )); GGraph.set_axes('xlabelint'=>xint.to_f)
    end

    if (@OPT_yint)
      yint, yint_div = @OPT_yint.split(":"); yint_div = 1 if yint_div == nil
      GGraph.set_axes('ytickint'=>(yint.to_f / yint_div.to_i )); GGraph.set_axes('ylabelint'=>yint.to_f)
    end


    if (@OPT_nolabels) then
      GGraph.set_axes("yunits"=>"")
      DCL::uzlset('LABELYL',false); DCL::uzlset('LABELYR',false)
      DCL::uzlset('LABELXT',false);DCL::uzlset('LABELXB',false)
    end



    # DCL::sgiset('NBITS', 32) # ラインタイプの指定を32bitに
  end

  # setup the window and draw_setup
  def window_setup
    if (DCLVERNUM > 710) then
      DCL.swlset("lsysfnt",true)
    else
      DCL.swlset("sysfont",true)
    end
    DCL.swcset("fontname","Hiragino Sans W5") unless ENV["SW_FONTNAME"] # font name
  #  DCL.swcset("fontname", "Helvectica")
    if (@OPT_lwf) then
      DCL.swpset("rlwfact", (@OPT_lwf.to_i)*DCL.swpget("rlwfact")) # line width factor for pdf.
    else
      DCL.swpset("rlwfact", 999)
    end

    DCL.swpset('lalt',true)
    if(@OPT_size) then
      if (@OPT_size.include?('auto'))
        sizefact = (@OPT_size.split('x')[1] || 1).to_f
        aspect = (@OPT_aspect || 2).to_f
        xn = @OPT_sldiv ? @OPT_sldiv.split(',')[1].to_i : 1
        yn = @OPT_sldiv ? @OPT_sldiv.split(',')[2].to_i : 1
        xlen = aspect*xn; ylen = 1.0*yn
        xsizefact = xlen/[xlen,ylen].max*sizefact
        ysizefact = ylen/[xlen,ylen].max*sizefact
      else
        xsizefact, ysizefact = @OPT_size.split("x")
        if (ysizefact == nil)
          xsizefact = xsizefact.to_f; ysizefact = xsizefact*0.75
        else
          xsizefact = xsizefact.to_f; ysizefact = ysizefact.to_f
        end
      end
    else
      xsizefact = 1.0; ysizefact = 0.75
    end

    if (@OPT_file) then
      @OPT_file, fname = @OPT_file.split(",")
      case @OPT_file
      when "png"; DCL.swiset('ifl',1); fname.sub!(".png","") if fname
        xsizefact = 2.0*xsizefact unless @OPT_size; ysizefact = 2.0*ysizefact unless @OPT_size
      when "eps"; DCL.swiset('ifl',2); fname.sub!(".eps","") if fname
      when "svg"; DCL.swiset('ifl',3); fname.sub!(".svg","") if fname
      when "pdf"; DCL.swiset('ifl',4); fname.sub!(".pdf","") if fname
      end
      DCL.swpset("iwidth",900*xsizefact)
      DCL.swpset("iheight",900*ysizefact)
      DCL.swpset("fname",fname) if fname

      DCL.gropn(2)

    else
      # 以下は set_vpsize に移動
      # Xwindowの大きさを@OPT_size倍にする。（Defalut: x 1.25）
      xsizefact = 1.25*xsizefact; ysizefact = 1.25*ysizefact
      iw = 900.0; ih = 650.0
      DCL.swpset("iwidth",900*xsizefact)
      DCL.swpset("iheight",900*ysizefact)

      DCL.gropn(@OPT_wsn||1)

    end

    @VIEWPORT = set_vpsize( VIEWPORT, (@OPT_aspect||2.0) ) # ここに移すのは良くないかも...

    draw_setup
    @flag_window_open += 1
    DCL.sgpset("LSOFTF","true") if (@OPT_fill and (@OPT_wsn.to_i == 2))
  end

  def set_draw_range(gp)
    xmin = nil; xmax = nil; ymin = nil; ymax = nil;
    if (@OPT_xrange) then
      xrange = ymd2date(@OPT_xrange)
      xmin, xmax = xrange.split(":")
      xmin = (eval xmin)
      xmax = (eval xmax)
      if (xmax.class == Date or xmin.class == Date)
        axis_unit = gp.axis(0).to_gphys.units unless @OPT_exch
        axis_unit = gp.axis(1).to_gphys.units if @OPT_exch
        since_date = DateTime.parse(axis_unit.to_s)
        cal_unit = axis_unit.to_s.split("since")[0].strip
        xmin = date_diff(xmin, since_date, cal_unit) if xmin.class == Date
        xmax = date_diff(xmax, since_date, cal_unit) if xmax.class == Date
        xmax = xmin + (xmax-xmin)/360*365.25 if @OPT_calendar360
        xmax = xmin + (xmax-xmin)/365*365.25 if @OPT_calendar365
      end
    end
    if (@OPT_yrange) then
      yrange = ymd2date(@OPT_yrange)
      ymin, ymax = yrange.split(":")
      ymin = (eval ymin)
      ymax = (eval ymax)
      if (ymax.class == Date or ymin.class == Date)
        axis_unit = gp.axis(1).to_gphys.units unless @OPT_exch
        axis_unit = gp.axis(0).to_gphys.units if @OPT_exch
        since_date = DateTime.parse(axis_unit.to_s)
        cal_unit = axis_unit.to_s.split("since")[0].strip
        ymin = date_diff(ymin, since_date, cal_unit) if ymin.class == Date
        ymax = date_diff(ymax, since_date, cal_unit) if ymax.class == Date
        ymax = ymin + (ymax-ymin)/360*365.25 if @OPT_calendar360
        ymax = ymin + (ymax-ymin)/365*365.25 if @OPT_calendar365
      end
    end
    GGraph.next_fig("window"=>[xmin,xmax,ymin,ymax])
  end

  def func(x)
    #y = sin(x*PI/180.0)
    y = (eval @OPT_top_axis.split(",")[0])
    return y
  end

  def ifunc(y)
    #x = asin(y)*180.0/PI
    x = (eval @OPT_top_axis.split(",")[1])
    return x
  end


  def draw(gp, draw_flag)
    if (draw_flag == "nofig") then
      val = gp.val
      if (val.class == NArray or val.class == NArrayMiss)
        p val[0]
      else
        p val
      end
      return
    end
    ## open work station and setup draw parameters
    window_setup if (@flag_window_open == 0 or @flag_window_open == 2)

    if (@OPT_itr.to_i >= 5) then  # itr >= 5 (地図投影)
      title = @OPT_title || gp.long_name || gp.name # タイトルの自動選択はしない
      @OPT_title = ""
    end
    # draw hontai
    case draw_flag
  #-------------------------------------------------------------------------------------------------
    when "line", "mark"
      if ( @Overplot == 1 )
        if (@OPT_right_axis or @OPT_top_axis)
          rsizet1 = DCL.uzpget("RSIZET1"); rsizet2 = DCL.uzpget("RSIZET2")
          DCL.uzpset("RSIZET1",0); DCL.uzpset("RSIZET2",0)
          DCL::uzlset('LABELYL',false); DCL::uzlset('LABELYR',false)
          DCL::uzlset('LABELXT',false);DCL::uzlset('LABELXB',false) #; DCL.uscget("cyfmt")
          if (@OPT_range) then
            @OPT_right_axis = (@OPT_range).split(/\s*,\s*/)[1] ? (@OPT_range).split(/\s*,\s*/)[1] : ""
            @OPT_range      = (@OPT_range).split(/\s*,\s*/)[0]
          end
        end

        set_draw_range(gp) if @OPT_xrange  #line plot で横軸の範囲を設定する場合

        if (@OPT_zerocenter == "") then
          max = [gp.min.val.abs, gp.max.val.abs].max;  @OPT_range = "#{-max}:#{max}"
        elsif (@OPT_zerocenter) then
          max = @OPT_zerocenter.to_f;  @OPT_range = "#{-max}:#{max}"
        end


        GGraph.next_axes("yside"=>"", "xside"=>"") if @OPT_axis_options

        if (draw_flag == "line") then
          GGraph.line((@OPT_step ? step_shape(gp) : gp),true,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||1),
                      "type" =>(@OPT_type ||1),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1],
                      "legend" =>!@OPT_nolegend)

        elsif (draw_flag == "mark") then
          GGraph.mark(gp, true,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||1),
                      "type" =>(@OPT_type ||1),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1])

          if (@OPT_line) then
            GGraph.line(gp, false,
                        "index"=>(@OPT_index||1),
                        "type" =>(@OPT_type ||1),
                        "exchange"=>@OPT_exch,
                        "annotate"=>@annotate)
          end

        end

        if (@OPT_axis_options) then
          @OPT_axis_options.split(",").each{|i|
            var, val = i.split("="); DCL.uspset(var, val)
          }
          GGraph.axes(gp.axis(0).to_gphys,gp,"yside"=>"LR","xside"=>"BT") unless @OPT_exch
          GGraph.axes(gp,gp.axis(0).to_gphys,"yside"=>"LR","xside"=>"BT") if @OPT_exch
        end

        if (@Overplot)
          @Overplot_ymin = gp.max.val;  @Overplot_ymax = gp.min.val
        end

        # x軸が時間でカレンダー属性を持つ場合、時間軸の単位(since)を合わせる。
        if (gp.axis(0).to_gphys.get_att("calendar")) then
            get_calendar(gp,0)
        end


        if (@OPT_right_axis or @OPT_top_axis) then
          DCL.uzpset("ROFFXT",0.0)
          DCL.uzpset("RSIZET1",rsizet1); DCL.uzpset("RSIZET2",rsizet2)
          DCL::uzlset('LABELYL', true); DCL::uzlset('LABELXB', true)
          if (@OPT_right_axis and !@OPT_top_axis) then
            GGraph.axes(gp.axis(0).to_gphys,gp,"yside"=>"L","xside"=>"BT")
          elsif (!@OPT_right_axis and @OPT_top_axis) then
            GGraph.axes(gp.axis(0).to_gphys,gp,"yside"=>"LR","xside"=>"B")
          elsif (@OPT_right_axis and @OPT_top_axis)
            GGraph.axes(gp.axis(0).to_gphys,gp,"yside"=>"L","xside"=>"B")
          end
          if (@OPT_top_axis) then # これ何する機能だっけ？
            #xmin = DCL.sgpget("UXMIN")  ; xmax = DCL.sgpget("UXMAX"); old_axis = gp.axis(0).pos.val
            # gp.axis(0).set_pos(VArray.new(NMath::sin(old_axis*PI/180.0), nil, "")); new_axis = gp.axis(0).pos.val
            # #DCL.grstrn(1)
            # DCL.uscset("cxfmt",""); DCL.usrset("xfac","-999") ;DCL.usrset("xoff","-999")
            # DCL.sgpset("UXMIN",new_axis.min)  ; DCL.sgpset("UXMAX",new_axis.max)
            DCL::uzlset('LABELXT', true)
            # GGraph.axes(gp.axis(0).to_gphys,gp,"yside"=>"n","xside"=>"T","xtickint"=>-999,"xlabelint"=>-999)
            xmin = DCL.sgpget("UXMIN"); xmax = DCL.sgpget("UXMAX");
            dxt = DCL.uspget("dxt"); dxl = DCL.uspget("dxl"); bx = DCL.rgnlt(xmin)
            ux1 = []; ux2 = []; ch = []
            if (true) then
              if (xmin*xmax > 0) then #ゼロを含まない場合
                ((xmax-xmin)/dxt+2).abs.to_i.times{|i| ux1 << bx+dxt*i if (bx+dxt*i>xmin && bx+dxt*i<xmax) }
                ((xmax-xmin)/dxl+2).abs.to_i.times{|i| ux2 << bx+dxl*i if (bx+dxl*i>xmin && bx+dxl*i<xmax) }
              else #ゼロを含む場合
                ux1 << 0; ux2 << 0; i = 1
                while (dxt*i <= [xmin.abs, xmax.abs].max)
                  ux1 << dxt*i; ux1 << -dxt*i; ux2 << dxl*i; ux2 << -dxl*i; i=i+1
                end
                ux1.select!{|a| a <= [xmin,xmax].max && a >= [xmin,xmax].min  }
                ux2.select!{|a| a <= [xmin,xmax].max && a >= [xmin,xmax].min  }
              end
            elsif (true) then
              nxmin = func(xmin); nxmax = func(xmax); i = 0
              bx1 = DCL.rgnlt([nxmin,nxmax].min); bx2 = DCL.rgngt([nxmin,nxmax].min); dxt = (bx1-bx2).abs*0.5
              r = [((nxmax-nxmin)/dxt/10).round.to_i, 2].max # 大きい目盛の割合
              # p bx1, dxt
              while (bx1 + dxt*i <= [nxmax,nxmin].max)
                ux2 << ifunc(bx1+dxt*i) if i%r == 0
                ux1 << ifunc(bx1+dxt*i)
                #ii = log10(((bx1+dxt*i)/bx1).abs).to_i  + 1#桁数-1
                #p bx1+dxt*i
                i = i + 1 #10**(ii-1)
              end
              ux1.select!{|a| a <= [xmin,xmax].max && a >= [xmin,xmax].min  }
              ux2.select!{|a| a <= [xmin,xmax].max && a >= [xmin,xmax].min  }
            end
            ch  = ux2.map{|i| func(i).round(1).to_s}
            DCL.uxaxlb("T",ux1,ux2,ch,1)
            # DCL.usrset("xfac","-999");DCL.usrset("xoff","-999");DCL.uscset("cxfmt","")
            # DCL.sgpset("UXMIN",xmin)  ; DCL.sgpset("UXMAX",xmax)
          end
        end

        if (@OPT_overplot_rm) then
          GGraph.line(gp.running_mean(0,@OPT_overplot_rm.to_i,10,@OPT_overplot_rm.to_i),false,
                    "index"=>(@OPT_index||1)+3,
                    "type" =>(@OPT_type ||1),"exchange"=>@OPT_exch)
        end

        if (@OPT_overplot_stddev) then
          GGraph.line(gp+@stddev,false,
                    "index"=>(@OPT_index||1),
                    "type" =>(@OPT_type ||3),"exchange"=>@OPT_exch)
          GGraph.line(gp-@stddev,false,
                    "index"=>(@OPT_index||1),
                    "type" =>(@OPT_type ||3),"exchange"=>@OPT_exch)
        end

        line_fill((@OPT_step ? step_shape(gp) : gp)) if (@OPT_fill == "")

  #-----OVERPLOT----------------------------
  #    elsif ( @OPT_overplot_color || @OPT_Opc )
      elsif ( !@OPT_nocolor )
        @Overplot_ymin = [@Overplot_ymin, gp.min.val].min
        @Overplot_ymax = [@Overplot_ymax, gp.max.val].max
        print "MIN: #{@Overplot_ymin}, MAX: #{@Overplot_ymax} \n" if @Overplot == @Overplot_max

        if (@OPT_right_axis)
          ymin = DCL.sgpget("UYMIN"); ymax = DCL.sgpget("UYMAX")
          gpmin = (@OPT_right_axis == "") ? gp.min.val : __split_range(@OPT_right_axis)[0]
          gpmax = (@OPT_right_axis == "") ? gp.max.val : __split_range(@OPT_right_axis)[1]
          # scale = (ymax - ymin)/(gpmax-gpmin); offset = gpmin - ymin
          # gp = (gp - gpmin)*scale + ymin
          DCL.sgpset("UYMIN",gpmin); DCL.sgpset("UYMAX",gpmax)
          GGraph.axes(nil,gp,"yside"=>"R","xside"=>"n","ytickint"=>-999,"ylabelint"=>-999)
        end

        gp = set_calendar(gp,0) if @Calendar

        if (draw_flag == "line") then
          GGraph.line((@OPT_step ? step_shape(gp) : gp), false,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index),
                      "type" =>(@OPT_type ||1),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1],
                      "legend" => !@OPT_nolegend
                      )
        elsif (draw_flag == "mark") then
          GGraph.mark(gp, false,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index),
                      "type" =>(@OPT_type ||1),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1])

          if (@OPT_line) then
            GGraph.line(gp, false,
                        "index"=>(@OPT_index||1),
                        "type" =>(@OPT_type ||1),
                        "exchange"=>@OPT_exch,
                        "annotate"=>@annotate)
          end

        end

        if (@OPT_right_axis)
          DCL::uzlset('LABELYR', true)
          #  DCL::uzrset('YFACT', 1.0/scale);  DCL::uzrset('YOFFSET', gpmin-ymin/scale)
          # -999は未定義を表し、内部で決める。前の軸の値が残っているので、改めて未定義 -999 を入れる必要がある。
          DCL.usrset("yfac","-999");DCL.usrset("yoff","-999");DCL.uscset("cyfmt","")
          #  DCL.uzpset('indext1', DCL.uzpget('indext1')+@OPT_index) # 小さい座標軸の目盛
          #  DCL.uzpset('indext2', DCL.uzpget('indext2')+@OPT_index) # 大きい座標軸の目盛
          DCL.uzpset('indexl1', DCL.uzpget('indexl1')+((@OPT_index || 1.0)/10)*10) # 小さいラベル、タイトル
          GGraph.axes(nil,gp,"yside"=>"R","xside"=>"n","ytickint"=>-999,"ylabelint"=>-999)
          DCL.uzpset('indexl1', DCL.uzpget('indexl1')-((@OPT_index || 1.0)/10)*10) # 小さいラベル、タイトル
        end
        if (@OPT_overplot_rm) then
          GGraph.line(gp.running_mean(0,@OPT_overplot_rm.to_i,10,@OPT_overplot_rm.to_i),false,
                  "index"=>((@OPT_index || 1.0)+3),
                  "type" =>(@OPT_type ||1),"exchange"=>@OPT_exch)
        end
        if (@OPT_overplot_stddev) then
          GGraph.line(gp+@stddev,false,
                  "index"=>(@OPT_index),
                  "type" =>(@OPT_type ||3),"exchange"=>@OPT_exch)
          GGraph.line(gp-@stddev,false,
                  "index"=>((@OPT_index || 1.0)),
                  "type" =>(@OPT_type ||3),"exchange"=>@OPT_exch)
        end
        if (@OPT_fill && @Overplot == 1) then
  #      if (@OPT_fill) then
          line_fill((@OPT_step ? step_shape(gp) : gp)) if (@OPT_fill == "")
          line_fill([@prev_gp,gp]) if (@OPT_fill.split(",").include?(@Overplot.to_s))
        end



      else #OVERPLOT WITH NOCOLOR
        if (@OPT_right_axis)
          ymin = DCL.sgpget("UYMIN"); ymax = DCL.sgpget("UYMAX")
          gpmin = (@OPT_right_axis == "") ? gp.min.val : __split_range(@OPT_right_axis)[0]
          gpmax = (@OPT_right_axis == "") ? gp.max.val : __split_range(@OPT_right_axis)[1]
          DCL.sgpset("UYMIN",gpmin); DCL.sgpset("UYMAX",gpmax)
          GGraph.axes(nil,gp,"yside"=>"R","xside"=>"n","ytickint"=>-999,"ylabelint"=>-999)
        end
        if (draw_flag == "line") then
          GGraph.line((@OPT_step ? step_shape(gp) : gp),
                      false,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||1),
                      "type" =>(@OPT_type ||@Overplot),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1],
                      "legend" => !@OPT_nolegend
                      )
        elsif (draw_flag == "mark") then
          GGraph.mark(gp, false,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||1),
                      "type" =>(@OPT_type ||@Overplot),
                      "exchange"=>@OPT_exch,
                      "annotate"=>@annotate,
                      "min" => __split_range(@OPT_range)[0],
                      "max" => __split_range(@OPT_range)[1]
                      )
        end
        if (@OPT_right_axis)
          DCL::uzlset('LABELYR', true)
          # -999は未定義を表し、内部で決める。前の軸の値が残っているので、改めて未定義 -999 を入れる必要がある。
          DCL.usrset("yfac","-999");DCL.usrset("yoff","-999");DCL.uscset("cyfmt","")
          GGraph.axes(nil,gp,"yside"=>"R","xside"=>"n","ytickint"=>-999,"ylabelint"=>-999)
        end
        if (@OPT_overplot_rm) then
          GGraph.line(gp.running_mean(0,@OPT_overplot_rm.to_i,10,@OPT_overplot_rm.to_i),false,
                    "index"=>(@OPT_index||1)+3,
                    "type" =>(@OPT_type ||@Overplot),"exchange"=>@OPT_exch)
        end

      end

      if ( @Overplot < @Overplot_max )
        @Overplot += 1
      else
        @Overplot = 1
      end

  #-------------------------------------------------------------------------------------------------
    when "full", "nocont", "noshade"
      if @Overplot == 1
        new_page = true
      else
        new_page = false
      end

      if (@OPT_xmean or @OPT_ymean) then
        vxmin = @VIEWPORT[0]; vxmax = @VIEWPORT[1]; vymin = @VIEWPORT[2]; vymax = @VIEWPORT[3]
        vxlen = vxmax - vxmin; vylen = vymax - vymin
        GGraph.set_fig('viewport'=>[vxmin, vxmin+vxlen*0.8 , vymax-vylen*0.8 , vymax])
        if (@OPT_ymean) then
          DCL.uzpset("LABELXB",false)
          DCL.udpset("lmsg",false)
          GGraph.next_axes("yunits"=>"")
        end
      end

      if (draw_flag == "noshade") then
        cont_tf = true
        mj = DCL.udpget('indxmj'); mn = DCL.udpget('indxmn')
        title = @OPT_title       ; annotate = @annotate
      else
        cont_tf = false
      end

      set_draw_range(gp) if (@OPT_xrange or @OPT_yrange)

      gp2D = gp.first2D.copy

  #    gp2D = cyclic_extension(gp2D, @OPT_cyclic_extension,4) if @OPT_cyclic_extension

      if (@OPT_bothsides) then
        vxmin = @VIEWPORT[0]; vxmax = @VIEWPORT[1]; vymin = @VIEWPORT[2]; vymax = @VIEWPORT[3]
        vxlen = vxmax - vxmin; vylen = vymax - vymin
        GGraph.next_fig('viewport'=>[vxmin, vxmin+vxlen*0.5 , vymin, vymax])
      end

      if (draw_flag == "full" or draw_flag == "nocont") then
        GGraph.tone(gp2D,new_page,
                    "title"=>@OPT_title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "auto"=>@auto,
                    "tonf"=>@tonf,
                    "tonb"=>@tonb,
                    "tonc"=>@tonc,
                    "fullcolor"=>@tone_fullcolor,
                    "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord,
                    "log"=>@log_int)
      end
      if (draw_flag == "full" or draw_flag == "noshade") then
        GGraph.contour(gp2D, cont_tf,
               "title"=>title,
               "annotate"=>annotate,
               "transpose"=>@OPT_exch,
               "levels"=>@OPT_clevels,
               "nozero"=>@OPT_nozero,
               "label"=>@OPT_label,
               "xcoord"=>@OPT_xcoord,
               "ycoord"=>@OPT_ycoord,
               "color"=>false,
               "log"=>@log_int )
      end

      if (@OPT_bothsides) then
        other_side =  [DCL.sgpget('plx')+180, -DCL.sgpget('ply'), DCL.sgpget('plrot')] # 惑星の裏側
        GGraph.next_fig('map_axis'=>other_side,'viewport'=>[vxmin+vxlen*0.5, vxmax, vymin, vymax], "new_frame"=>false)
        GGraph.tone(gp2D,new_page,
                    "title"=>@OPT_title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "auto"=>@auto,
                    "tonf"=>@tonf,
                    "tonb"=>@tonb,
                    "tonc"=>@tonc,
                    "fullcolor"=>@tone_fullcolor,
                    "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord,
                    "log"=>@log_int)
        DCL.sgsvpt(vxmin, vxmax, vymin, vymax)
      end

      # テスト
      if (false) then
        gp2D = gp.first2D.val
        max = gp.first2D.max
        min = gp.first2D.min
        xm, ym = gp2D.shape
        data = (gp2D.to_na - min) / (max - min) * 255
        xcoord = gp.first2D.coord(0).val
        ycoord = gp.first2D.coord(1).val
        xm.times{|i|
          ym.times{|j|
            vx0, vy0  = DCL.stftrf(xcoord[[0, i-1].max],ycoord[j])
            vx1, vy1  = DCL.stftrf(xcoord[i],ycoord[j])
            vx2, vy2  = DCL.stftrf(xcoord[[i+1, xm-1].min],ycoord[j])
            vr1 = ((vx1 - vx0)**2 + (vy1 - vy0)**2)**0.5
            vr2 = ((vx2 - vx1)**2 + (vy2 - vy1)**2)**0.5
            vr = [vr1, vr2].max*1.1
            rgb = DCL.isgrgb(data[i,j],data[i,j],data[i,j])
            DCL.sgpmxu([xcoord[i]],[ycoord[j]],11,7,rgb,vr)
          }
        }
      end



      if ( @Overplot < @Overplot_max )
        @Overplot += 1
      else
        @Overplot = 1
      end



  #-------------------------------------------------------------------------------------------------
    when "fullcolor"
      if (@OPT_xmean or @OPT_ymean) then
        vxmin = @VIEWPORT[0]; vxmax = @VIEWPORT[1]; vymin = @VIEWPORT[2]; vymax = @VIEWPORT[3]
        vxlen = vxmax - vxmin; vylen = vymax - vymin
        GGraph.set_fig('viewport'=>[vxmin, vxmin+vxlen*0.8 , vymax-vylen*0.8 , vymax])
        if (@OPT_ymean) then
          DCL.uzpset("LABELXB",false)
          DCL.udpset("lmsg",false)
          GGraph.next_axes("yunits"=>"")
        end
      end

      set_draw_range(gp) if (@OPT_xrange or @OPT_yrange)

      zmin, zmax = tone_full(gp,
                  true,
                  "title"=>@OPT_title,
                  "annotate"=>@annotate,
                  "transpose"=>@OPT_exch,
                  "levels"=>@OPT_slevels,
                  "patterns"=>@OPT_patterns,
                  "xcoord"=>@OPT_xcoord,
                  "ycoord"=>@OPT_ycoord,
                  "log"=>@log_int,
                  "min" => __split_range(@OPT_range)[0],
                  "max" => __split_range(@OPT_range)[1]
                  )
      unless (@OPT_nocont) then
        GGraph.contour(gp,
                       false,
                       "transpose"=>@OPT_exch,
                       "levels"=>@OPT_clevels,
                       "nozero"=>@OPT_nozero,
                       "label"=>@OPT_label,
                       "xcoord"=>@OPT_xcoord,
                       "ycoord"=>@OPT_ycoord,
                       "log"=>@log_int)
      end

  #-------------------------------------------------------------------------------------------------
    when "histogram1D"
      if (@Overplot == 1) then
        if (@OPT_range)
          ranges = (@OPT_range).split(/\s*,\s*/)
          xmin = __split_range(ranges[0])[0]
          xmax = __split_range(ranges[0])[1]
          ymin = __split_range(ranges[1])[0]
          ymax = __split_range(ranges[1])[1]
        end
        nb = -(@OPT_int).to_i if ((@OPT_int).to_i < 0)
        gph = gp.histogram("nbins"=>nb,"min"=>xmin,"max"=>xmax)
        if (@OPT_histogram != "")
          gph = gph/gph.sum*100
          gph.set_att("long_name","ratio (%) of bins")
          gph.units="%"
        end
        DCL.uusfri(@OPT_index||1) #index of box
        DCL.uusfrt(@OPT_type||1) #type of box
        ci = (@OPT_index||1.0)/10; wi = (@OPT_index||1.0)%10; pi = (@Overplot+1)%7
        fi = ci*1000+pi*100+wi*10+2
        GGraph.histogram(gph,
                         true,
                         "title"=>@OPT_title,
                         "exchange"=>@OPT_exch,
                         "window"=>[xmin,xmax,ymin,ymax],
                         "fill"=>@OPT_fill,
                         "fill_pattern"=>fi)
      else
        if (@OPT_range)
          ranges = (@OPT_range).split(/\s*,\s*/)
          xmin = __split_range(ranges[0])[0]
          xmax = __split_range(ranges[0])[1]
          ymin = __split_range(ranges[1])[0]
          ymax = __split_range(ranges[1])[1]
        end

        nb = -(@OPT_int).to_i if ((@OPT_int).to_i < 0)
        gph = gp.histogram("nbins"=>nb,"min"=>xmin,"max"=>xmax)
        if (@OPT_histogram != "")
          gph = gph/gph.sum*100
          gph.set_att("long_name","ratio (%) of bins")
          gph.units="%"
        end
        DCL.uusfri((@OPT_index || 1.0)) #index of box
        DCL.uusfrt(@OPT_type||1) #type of box
        ci = (@OPT_index || 1.0)/10; wi = (@OPT_index || 1.0)%10; pi = (@Overplot+1)%7
        fi = ci*1000+pi*100+wi*10+2
        GGraph.histogram(gph,false,
                         "exchange"=>@OPT_exch,
                         "fill"=>@OPT_fill,
                         "fill_pattern"=>fi
                         )
      end

      if ( @Overplot < @Overplot_max )
        @Overplot += 1
      else
        @Overplot = 1
      end

    end #end case

  #-------------------------------------------------------------------------------------------------
    # plot linear line
    draw_linearline(gp) if @OPT_linearline && @Overplot == 1
  #-------------------------------------------------------------------------------------------------
    if (@OPT_itr.to_i >= 5) then
      title(title)
      @OPT_title = title
      @OPT_title = nil if (title == gp.long_name || title == gp.name)
    end

    # lon-lat axis for itr=10
      if (@OPT_itr.to_i == 10) then
        DCL.uxaxdv("T",10,30); DCL.uxaxdv("B",10,30)
        DCL.uyaxdv("L",10,20); DCL.uyaxdv("R",10,20)
      end



    title(@OPT_subtitle, 1) if @OPT_subtitle && @Overplot == 1

    # test for counter
    # counter("t = ", 0.0, 1.00, " Earth days")
    # DCL::sgtxzr(0.135 , 0.3725, "90", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.865 , 0.3725, "90", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.3255 , 0.185, "0", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.676 , 0.185, "0", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.3255 , 0.56, "180", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.676 , 0.56, "180", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
    # DCL::sgtxzr(0.5 , 0.375, "270", DCL.uzpget('rsizec1')*0.75, 0, 0, 993)

  #---small window below the figure for y-mean line plot
    if (@OPT_ymean) then
      itr =  (@OPT_itr.to_i == 3) ? 3 : 1
      mdim = @OPT_exch ? 0 : 1
      if (@OPT_itr.to_i < 5) then
        GGraph.next_fig('viewport'=>[vxmin, vxmin+vxlen*0.8 , vymin, vymax-vylen*0.8-0.01],
                      "new_frame"=>false, "itr"=>itr)
        GGraph.next_axes("xtitle"=>nil, "ytitle"=>gp.axis(mdim).to_gphys.name+" mean", "xunits"=>nil, "yunits"=>gp.units.to_s)
      else
        GGraph.next_fig('viewport'=>[vxmin, vxmin+vxlen*0.8 , vymin, vymax-vylen*0.8-0.01],
                      "new_frame"=>false, "itr"=>itr)
        GGraph.next_axes("xtitle"=>"", "ytitle"=>gp.axis(mdim).to_gphys.name+" mean", "xunits"=>nil, "yunits"=>gp.units.to_s,
                         "xtickint"=>10.0, "xlabelint"=>90.0)
      end

      DCL::uzlset('LABELXB', true)
      GGraph.line(gp.mean(mdim),true, "title"=>"", "annotate"=>false,
                  "min" => __split_range(@OPT_range)[0],
                  "max" => __split_range(@OPT_range)[1],
                  "index"=>(@OPT_index||1), "type" =>(@OPT_type ||1)
                  )

      GGraph.set_fig("viewport"=>@VIEWPORT)
    end

    #-------------------------------------------------------------------------------------------------
      # color bar
      if ( ( draw_flag == "full") || ( draw_flag == "nocont") ) && @colorbar then
        if (@OPT_aspect.to_f <= 2.0 && @OPT_itr.to_i <= 4 && @OPT_aspect) then
          GGraph::color_bar("left"=> false,"landscape" => false, "units" => "["+gp.units.to_s+"]")
        elsif (@OPT_itr.to_i <= 4) then
          GGraph::color_bar("left"=> true,"landscape" => true, "units" => "["+gp.units.to_s+"]")
        else
          GGraph::color_bar("left"=> true,"landscape" => true, "units" => "["+gp.units.to_s+"]")
        end
    #                      "vlength" => 0.55
      elsif (draw_flag == "fullcolor" && @colorbar) then
        vx0, vx1, vy0, vy1 = DCL.sgqvpt # viewportの記憶
        aspect = (vx1 - vx0) / (vy1 - vy0)
        if (aspect >= 2.0)
          offset = vy0 - 0.1875
          DCL.uixbar(0.35,0.65,0.035+offset,0.0575+offset,zmin,zmax,"B")
        else
          DCL.uzpset("LABELXB",false);  DCL.uzpset("LABELYL",false)
          DCL.uiybar(vx1+0.02,vx1+0.045,vy0,vy0+0.3,zmin,zmax,"R")
          DCL.uxsttl("T", "["+gp.units.to_s+"]", -1)
          DCL.uzpset("LABELYR",true)
          DCL.uzinit; DCL.usdaxs
          #DCL.uzpset("ROFFYR",0); DCL.uzpset("IROTCYR",4)
          #DCL.uysttl("R", "["+gp.units.to_s+"]", 1)
          DCL.sgsvpt(vx0, vx1, vy0, vy1) # viewportを戻す
        end
        DCL.uzpset("LABELXB",true);  DCL.uzpset("LABELYL",true)
        DCL.uzpset("LABELXT",false);  DCL.uzpset("LABELYR",false)
      end

  #-------------------------------------------------------------------------------------------------
      # mask shading
      maskshading(gp) if @OPT_maskshading
  #-------------------------------------------------------------------------------------------------






  # ---- small window on right side for x-mean line plot
      if (@OPT_xmean) then
        itr =  (@OPT_itr.to_i == 2) ? 2 : 1
        if (@OPT_itr.to_i < 5) then
          GGraph.next_fig('viewport'=>[vxmin+vxlen*0.8+0.01, vxmax,   vymax-vylen*0.8 , vymax],
                        "new_frame"=>false, "itr"=>itr)
          GGraph.next_axes("xtitle"=>"", "ytitle"=>"", "xunits"=>gp.units.to_s, "yunits"=>"")
        else
          GGraph.next_fig('viewport'=>[vxmin+vxlen*0.8+0.01, vxmax,   vymax-vylen*0.8+0.00325 , vymax-0.00325],
                        "new_frame"=>false, "itr"=>itr)
          GGraph.next_axes("xtitle"=>"", "ytitle"=>"", "xunits"=>gp.units.to_s, "yunits"=>"",
                           "ytickint"=>10.0, "ylabelint"=>20.0)
        end
        DCL::uzlset('LABELXB', true) unless @OPT_nolabels
        DCL::uzlset('LABELYL', false)
        mdim = @OPT_exch ? 1 : 0
        GGraph.line(gp.mean(mdim),true, "exch"=>true, "title"=>"", "annotate"=>false,
                    "min" => __split_range(@OPT_range)[0],
                    "max" => __split_range(@OPT_range)[1],
                    "index"=>(@OPT_index||1), "type" =>(@OPT_type ||1)
                    )
        if ( DCL.sgpget("UXMIN")*DCL.sgpget("UXMAX") < 0) then
          DCL::sgplzu([0,0],[DCL.sgpget("UYMIN"),DCL.sgpget("UYMAX")],3,1) # ゼロ線を描く
        end

        title(gp.axis(mdim).to_gphys.name+" mean", 1) unless @OPT_title == ""
        DCL::uzlset('LABELYL', true)
        GGraph.set_fig("viewport"=>@VIEWPORT)
      end

  end

  def draw_multi_vars(gp, draw_flag)

    ## open work station and setup draw parameters
    window_setup if (@flag_window_open == 0 or @flag_window_open == 2)

    if (@OPT_range and !@OPT_rmap)
      ranges = (@OPT_range).split(/\s*,\s*/)
      xmin = __split_range(ranges[0])[0]
      xmax = __split_range(ranges[0])[1]
      ymin = __split_range(ranges[1])[0]
      ymax = __split_range(ranges[1])[1]
      zmin = __split_range(ranges[2])[0]
      zmax = __split_range(ranges[2])[1]
      GGraph.next_fig("window"=>[xmin,xmax,ymin,ymax]) unless @OPT_fullcolor
    end

    if (@OPT_itr.to_i >= 5) then  # itr >= 5 (地図投影)
      title = @OPT_title || gp[0].long_name || gp[0].name # タイトルの自動選択はしない
      @OPT_title = ""
    end
    if (@OPT_anim) then
      adim = gp[0].axis(@OPT_anim).to_gphys.val.to_a
      adim = adim.reverse if @OPT_reverse || @OPT_Gr
      gp_all = gp.clone
    else
      adim = [1]
    end

    anim_counter = 0

    while (adim[0])
      if (@OPT_anim) then
        gp = gp_all.map{|g| g.cut(@OPT_anim=>adim[0]) }
        anim_counter += 1
      end
      adim.shift

      case draw_flag
    #-------------------------------------------------------------------------------------------------
      when "scatter"
        xtitle, ytitle = (@OPT_scatter).split(/\s*,\s*/)
        GGraph.next_axes("xtitle"=>xtitle, "ytitle"=>ytitle)
        GGraph.scatter(gp[0],gp[1],
                       true,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||@index_array[0]),
                      "type" =>(@OPT_type ||@type_array[0]),
                      "annotate"=>@annotate,
                      "correlation"=>@annotate)
        GGraph.regression_line(gp[0],gp[1],"annot_intercept"=>@annotate, "annot_slope"=> @annotate, "index"=>(@OPT_index||@index_array[0])) # 回帰直線を引く

        if (gp.length/2 > 1) then
          (gp.length/2 - 1).times{|i|
            GGraph.scatter(gp[2*i+2],gp[2*i+3],
                           true,
                          "index"=>(@index_array[i+1]||1),
                          "type" =>(@type_array[i+1]||2),
                          "annotate"=>@annotate,
                          "correlation"=>@annotate)
            GGraph.regression_line(gp[2*i+2],gp[2*i+3],"annot_intercept"=>@annotate,"annot_slope"=> @annotate,  "index"=>(@index_array[i+1])) # 回帰直線を引く

          }
        end

    #-------------------------------------------------------------------------------------------------
      when "color_scatter"
        xtitle, ytitle = (@OPT_color_scatter).split(/\s*,\s*/)
        GGraph.next_axes("xtitle"=>xtitle, "ytitle"=>ytitle)

        # tmp = NArray.float(gp[2].length).indgen!+2020
        # gp[2].replace_val(tmp)
        # gp[5].replace_val(tmp)
        # gp[8].replace_val(tmp)


        # --color_scatterのときに、2つのgturlしかないときは、最初のgturlの
        # ２つ目の軸(緯度を想定), なければ１つめの軸（時間を想定）の値で色付けする。
        if (true) then
          gp << (gp[0]*0.0 + NArray.float(gp[0].length).indgen! + 2021)
        elsif (gp[2]==nil && gp[0].rank > 1) then
          gp << gp[0]*0.0 + gp[0].axis(1).to_gphys
        elsif (gp[2]==nil)
          gp << gp[0]*0.0 + gp[0].axis(0).to_gphys
        end
        GGraph.color_scatter(gp[0],gp[1],gp[2],
                       true,
                      "title"=>@OPT_title,
                      "index"=>(@OPT_index||1),
                      "type" =>(@OPT_type ||2),
                      "annotate"=>@annotate,
                      "min" => zmin,
                      "max" => zmax,
                      "nlev" => ( @OPT_sint   || @OPT_interval || @OPT_int ),
                      "correlation"=>true
                       )
        GGraph.regression_line(gp[0],gp[1],"annot_intercept"=>true) # 回帰直線を引く

        GGraph.color_bar( "left"=> false, "landscape" => false)

        if (gp.length/3 > 1) then
          (gp.length/3 - 1).times{|i|
            GGraph.color_scatter(gp[3*i+3],gp[3*i+4],gp[3*i+5],
                                false,
                                "index"=>(@index_array[i+1]||1),
                                "type" =>(@type_array[i+1]||2),
                                "annotate"=>@annotate,
                                "min" => zmin,
                                "max" => zmax,
                                "correlation"=>true)
          }
        end

    #-------------------------------------------------------------------------------------------------
      when "contour_over_tone"
        if (@OPT_bothsides) then
          vxmin = @VIEWPORT[0]; vxmax = @VIEWPORT[1]; vymin = @VIEWPORT[2]; vymax = @VIEWPORT[3]
          vxlen = vxmax - vxmin; vylen = vymax - vymin
          GGraph.next_fig('viewport'=>[vxmin, vxmin+vxlen*0.5 , vymin, vymax])
        end
        title = @OPT_notitle ? "" : gp[0].name + " & "+gp[1].name
        GGraph.tone(gp[0],
                    true,
                    "title"=> title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "auto"=>@auto,
                    "tonf"=>@tonf,
                    "tonb"=>@tonb,
                    "tonc"=>@tonc,
                    "fullcolor"=>@tone_fullcolor,
                    "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord,
                    "log"=>@log_int
                    )
          # GGraph.contour(gp[0],
      		#             true, "title"=>"",
       		#             "transpose"=>@OPT_exch,
      		#             "levels"=>@OPT_slevels,
      		#             "nozero"=>@OPT_nozero,
          #             "label"=>@OPT_label,
          #             "xcoord"=>@OPT_xcoord,
          #             "ycoord"=>@OPT_ycoord,
          #             "log"=>@log_int,
          #             "coloring"=>true
          #             )


  #      DCL.sgscmn(5) # 色コンターの色番号
  #       DCL.udpset('indxmj', 25) #計曲線
  #       DCL.udpset('indxmn', 27) #主曲線

        min_crange, max_crange = __split_range(@OPT_crange)
        GGraph.set_linear_contour_options('int' => @OPT_cint.to_f, 'min' =>min_crange,'max' =>max_crange)
        GGraph.contour(gp[1], false, "title"=>"", "transpose"=>@OPT_exch, "levels"=>@OPT_clevels,
    		            "nozero"=>@OPT_nozero, "label"=>@OPT_label, "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord, "log"=>@log_int, "coloring"=>false)#,
  #                  "clr_min"=>50, "clr_max"=>60, "line_type"=>3 )

  #      GGraph.contour(gp[2], false, "title"=>"", "transpose"=>@OPT_exch, "levels"=>@OPT_clevels,
  #                   "nozero"=>@OPT_nozero, "label"=>@OPT_label, "xcoord"=>@OPT_xcoord,
  #                   "ycoord"=>@OPT_ycoord, "log"=>@log_int, "coloring"=>true)#,
  #                   "clr_min"=>50, "clr_max"=>60, "line_type"=>3 )


        DCL.sgscmn(@OPT_clrmap||63) # カラーマップ戻す

        if (@OPT_bothsides) then
          other_side =  [@map_axis[0]+180, -@map_axis[1], @map_axis[2]] # 惑星の裏側
          GGraph.next_fig('map_axis'=>other_side,'viewport'=>[vxmin+vxlen*0.5, vxmax, vymin, vymax], "new_frame"=>false)
          GGraph.tone(gp[0],
                      true,
                      "title"=>"",
                      "annotate"=>@annotate,
                      "transpose"=>@OPT_exch,
                      "levels"=>@OPT_slevels,
                      "patterns"=>@OPT_patterns,
                      "auto"=>@auto,
                      "tonf"=>@tonf,
                      "tonb"=>@tonb,
                      "tonc"=>@tonc,
                      "fullcolor"=>@tone_fullcolor,
                      "xcoord"=>@OPT_xcoord,
                      "ycoord"=>@OPT_ycoord,
                      "log"=>@log_int
                      )
          DCL.sgscmn(5) # 白黒コンター
          GGraph.contour(gp[1],
      		            false, "title"=>"",
       		            "transpose"=>@OPT_exch,
      		            "levels"=>@OPT_clevels,
      		            "nozero"=>@OPT_nozero,
                      "label"=>@OPT_label,
                      "xcoord"=>@OPT_xcoord,
                      "ycoord"=>@OPT_ycoord,
                      "log"=>@log_int,
                      "coloring"=>false
                      )
          DCL.sgscmn(@OPT_clrmap||63) # カラーマップ戻す
          DCL.sgsvpt(vxmin, vxmax, vymin, vymax)
        end

        maskshading(gp) if @OPT_maskshading

        GGraph::color_bar("left"=> true,"landscape" => true) unless @OPT_nocolorbar
  #      GGraph::color_bar("left"=> false,"landscape" => false)
        # test for counter
  #      counter("t = ", 0.0, 1.00, " Earth days")
  #      DCL::sgtxzr(0.135 , 0.3725, "90", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.865 , 0.3725, "90", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.3255 , 0.185, "0", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.676 , 0.185, "0", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.3255 , 0.56, "180", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.676 , 0.56, "180", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)
  #      DCL::sgtxzr(0.5 , 0.375, "270", DCL.uzpget('rsizec1')*0.75, 0, 0, 993)
  #      DCL::sgtxzr(0.72 , 0.1475, "[m s-1]", DCL.uzpget('rsizec1')*0.75, 0, 0, 3)


    #-------------------------------------------------------------------------------------------------
      when "histogram2D"
        nb = -(@OPT_int).to_i if ((@OPT_int).to_i < 0)
        gph = GAnalysis.histogram2D(gp[0],gp[1],"nbins0"=>nb,"nbins1"=>nb,
                                               "min0"=>xmin,"max0"=>xmax,
                                               "min1"=>ymin,"max1"=>ymax)
        if (@OPT_histogram2D != "")
          gph = gph/gph.sum*100
          gph.set_att("long_name","ratio (%) of bins")
          gph.units="%"
        end

        GGraph.tone(gph,true,
                    "title"=>@OPT_title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "min" => __split_range(@OPT_srange)[0],
                    "max" => __split_range(@OPT_srange)[1],
                    "nlev"=> (@OPT_sint.to_i < 0 ? -@OPT_sint.to_i : nil),
                    "interval"=>(@OPT_sint.to_i > 0 ? @OPT_sint.to_i : nil),
                    "tonb"=>true,
                    "log"=>@log_int
                    )
        GGraph::color_bar("left"=> true,"landscape" => true)

      when "rmap"
        if (@OPT_rmap.to_i.to_s == @OPT_rmap ) then
          dim = @OPT_rmap.to_i
        else
          dim = gp[0].axnames.index(@OPT_rmap)
        end

        gprmap, ndiv = corelation(gp[0],gp[1],dim)
    #    gprmap.set_att("long_name","correlation along "+gp[0].axnames[dim])
        GGraph.tone(gprmap,true,
                    "title"=>@OPT_title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "auto"=>@auto,
                    "tonf"=>@tonf,
                    "tonb"=>@tonb,
                    "tonc"=>@tonc,
                    "fullcolor"=>@tone_fullcolor,
                    "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord,
                    "log"=>@log_int
                    )
        GGraph::color_bar("left"=> true,"landscape" => true)

      when "fullcolor2"


        tone_full(gp,
                    true,
                    "title"=>@OPT_title,
                    "annotate"=>@annotate,
                    "transpose"=>@OPT_exch,
                    "levels"=>@OPT_slevels,
                    "patterns"=>@OPT_patterns,
                    "xcoord"=>@OPT_xcoord,
                    "ycoord"=>@OPT_ycoord,
                    "log"=>@log_int,
                    "min" => [xmin,ymin],
                    "max" => [xmax,ymax]
                    )
        when "vector"
          londim, latdim = GAnalysis::Planet.find_lon_lat_dims(gp[0])
          @xintv, @yintv = @OPT_vint.split(",") if @OPT_vint
          @xintv = (@xintv || (gp[0].shape[0]/30).to_i + 1).to_i
          @yintv = (@yintv || gp[0].shape[1].to_f/gp[0].shape[0]*@xintv*(@OPT_aspect||2).to_f).to_i
          @vfact = (@OPT_vfact || 1).to_f
          @vkeep = true if (anim_counter != 1 && @OPT_vkeep )
          if (@OPT_itr == '5') then
             GGraph.vector( gp[0], gp[1], true,
                           "title"=>@OPT_title, "annotate"=>@annotate, "exchange"=>@OPT_exch,
                           "flow_vect"=>false,  "flow_itr5"=>true,     "xintv"=>@xintv,
                           "yintv"=>@yintv,     "factor"=>@vfact,      "unit_vect"=>true,
                           "max_unit_vect"=>true, "ux_unit"=>@ux_unit, "uy_unit"=>@uy_unit,
                           "keep" =>@vkeep)

          elsif ( (londim == nil && latdim == 0) || (londim == 0 && latdim == nil) ) then # meridonal vector plot
            require 'numru/ggraph_on_merdional_section'
            if (gp.size == 2) then # vectors only
              GGraph.vector_on_merdional_section(gp[0], gp[1], true,
                                'fact'=>@vfact, 'xintv'=>6, 'yintv'=>1,'unit'=>true, 'annot'=>false, "use_before_scale" => @vkeep)
            else
              GGraph.tone(gp[0],true,
                                "title"=>@OPT_title, "annotate"=>@annotate,
                                "transpose"=>@OPT_exch, "levels"=>@OPT_slevels,
                                "patterns"=>@OPT_patterns, "auto"=>@auto,
                                "tonf"=>@tonf, "tonb"=>@tonb, "tonc"=>@tonc,
                                "fullcolor"=>@tone_fullcolor,
                                "xcoord"=>@OPT_xcoord, "ycoord"=>@OPT_ycoord,
                                "log"=>@log_int
                                )
              GGraph::color_bar("left"=> true,"landscape" => true) unless @OPT_nocolorbar
              if (gp[3]) then
                GGraph.set_linear_contour_options('int' => @OPT_cint.to_f, 'min' =>min_crange,'max' =>max_crange)
                GGraph.contour(gp[3], false, "title"=>"", "transpose"=>@OPT_exch, "levels"=>@OPT_clevels,
                              "nozero"=>@OPT_nozero, "label"=>@OPT_label, "xcoord"=>@OPT_xcoord,
                              "ycoord"=>@OPT_ycoord, "log"=>@log_int, "coloring"=>false )
              end
              GGraph.vector_on_merdional_section(gp[1], gp[2], false,
                                'fact'=>@vfact, 'xintv'=>@xintv, 'yintv'=>@yintv,'unit'=>true, 'annot'=>false, "use_before_scale" => @vkeep)

              # DCL.ugpset('index',13) # ベクトルのラインインデックス
              # GGraph.vector_on_merdional_section(gp[1]*gp[2].le(0), gp[2]*gp[2].le(0), false,
              #                   'fact'=>@vfact, 'xintv'=>6, 'yintv'=>1,'unit'=>true, 'annot'=>false, "use_before_scale" => (anim_counter != 1))
              # DCL.ugpset('index',23) # ベクトルのラインインデックス
              # GGraph.vector_on_merdional_section(gp[1]*gp[2].ge(0), gp[2]*gp[2].ge(0), false,
              #                   'fact'=>@vfact, 'xintv'=>6, 'yintv'=>1,'unit'=>true, 'annot'=>false, "use_before_scale" => true)

                                min_crange, max_crange = __split_range(@OPT_crange)
            end

          else # horizontal vector plots
            if (gp.size == 2) then # vectors only
              GGraph.vector(gp[0], gp[1], true,
                               "title"=>@OPT_title, "annotate"=>@annotate,
                               "exchange"=>@OPT_exch, "flow_vect"=>true,
                               "xintv"=>@xintv, "yintv"=>@yintv,
                               "factor"=>@vfact, "unit_vect"=>true,
                               "max_unit_vect"=>false, "ux_unit"=>@ux_unit,
                               "uy_unit"=>@uy_unit, "keep"=>@vkeep
                               )
            else
              GGraph.tone(gp[0],true,
                                "title"=>@OPT_title, "annotate"=>@annotate,
                                "transpose"=>@OPT_exch, "levels"=>@OPT_slevels,
                                "patterns"=>@OPT_patterns, "auto"=>@auto,
                                "tonf"=>@tonf, "tonb"=>@tonb, "tonc"=>@tonc,
                                "fullcolor"=>@tone_fullcolor,
                                "xcoord"=>@OPT_xcoord, "ycoord"=>@OPT_ycoord,
                                "log"=>@log_int
                                )
              GGraph::color_bar("left"=> true,"landscape" => true) unless @OPT_nocolorbar
              if (gp[3]) then
                GGraph.set_linear_contour_options('int' => @OPT_cint.to_f, 'min' =>min_crange,'max' =>max_crange)
                GGraph.contour(gp[3], false, "title"=>"", "transpose"=>@OPT_exch, "levels"=>@OPT_clevels,
                              "nozero"=>@OPT_nozero, "label"=>@OPT_label, "xcoord"=>@OPT_xcoord,
                              "ycoord"=>@OPT_ycoord, "log"=>@log_int, "coloring"=>true )
              end
              # DCL.ugpset("index",995)
              GGraph.vector( gp[1], gp[2], false,
                              "title"=>@OPT_title, "annotate"=>@annotate,
                              "exchange"=>@OPT_exch, "flow_vect"=>true,
                              "xintv"=>@xintv, "yintv"=>@yintv,
                              "factor"=>@vfact, "unit_vect"=>true,
                              "max_unit_vect"=>false, "ux_unit"=>@ux_unit,
                              "uy_unit"=>@uy_unit, "keep"=>@vkeep
                              )
            end
          end
          maskshading(gp) if @OPT_maskshading

      end # end case
    #-------------------------------------------------------------------------------------------------
      draw_linearline(gp[0]) if @OPT_linearline
    #-------------------------------------------------------------------------------------------------
        if (@OPT_itr.to_i >= 5) then
          title(title)
          @OPT_title = title
          @OPT_title = nil if (title == gp[0].long_name || title == gp[0].name)
        end
        title(@OPT_subtitle, 1) if @OPT_subtitle
      end # end while (for anim loop)
  end

  def draw_linearline(gp=nil)
    if (gp) then
      unit_0 = gp.axis(0).to_gphys.units.to_s
      flag_0 = unit_0.include?("since")
      unit_1 = gp.axis(1).to_gphys.units.to_s if gp.rank >= 2
      flag_1 = unit_1.include?("since") if gp.rank >= 2
    end
    index = 1; type  = 1
    line_paras = (@OPT_linearline).split(/\s*,\s*/)
    line_paras.each{|lp|
      axis, value = lp.split(/\s*=\s*/)
      if (axis=="i")
        index = value.to_i
      elsif (axis=="t")
        type  = value.to_i
      elsif (lp=="x=y" or lp == "y=x")
        minval = [DCL.sgpget("UXMIN"), DCL.sgpget("UYMIN")].max
        maxval = [DCL.sgpget("UXMAX"), DCL.sgpget("UYMAX")].min
        DCL::sgplzu([minval,maxval],[minval,maxval],type,index)
      elsif (lp.include?("x^"))
        raise "linear c*x^p line can be drawn only when itr = 4 (log-log plot)" unless @OPT_itr.to_i == 4
        cn, pn = lp.split("x^");pn = (eval pn).to_f; cn = ((eval cn) || 1).to_f
        spx = DCL.sgpget("UXMIN"); epx = DCL.sgpget("UXMAX"); spy = DCL.sgpget("UYMAX")
        spy = cn*spx**pn; epy = cn*epx**pn
        ymax = DCL.sgpget("UYMAX"); ymin = DCL.sgpget("UYMIN")
        if (spy < ymin) then spy = ymin; spx = (ymin/cn)**(1.0/pn) end
        if (epy < ymin) then epy = ymin; epx = (ymin/cn)**(1.0/pn) end
        if (spy > ymax) then spy = ymax; spx = (ymax/cn)**(1.0/pn) end
        if (epy > ymax) then epy = ymax; epx = (ymax/cn)**(1.0/pn) end
        DCL::sgplzu([spx,epx],[spy,epy],type,index)
  #    elsif (axis.include?("x") && value.include?("y"))
  #    elsif (axis.include?("y") && value.include?("x"))
      elsif (axis=="x")
        if (gp.rank == 1) then
          value = ymd2time(unit_0, value)
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/360*365.25 if @OPT_calendar360 && flag_0 && !@OPT_exch
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/365*365.25 if @OPT_calendar365 && flag_0 && !@OPT_exch
          # ↑ DCL::sgplzuだと、UXMINから1年365日のスケールで座標が取られているみたい。
        elsif (@OPT_exch) then
          value = ymd2time(unit_1, value)
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/360*365.25 if @OPT_calendar360 && flag_1
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/365*365.25 if @OPT_calendar365 && flag_1
        else
          value = ymd2time(unit_0, value)
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/360*365.25 if @OPT_calendar360 && flag_0
          value = DCL.sgpget("UXMIN") + (value - DCL.sgpget("UXMIN"))/365*365.25 if @OPT_calendar365 && flag_0
        end
        DCL::sgplzu([value,value],[DCL.sgpget("UYMIN"),DCL.sgpget("UYMAX")],type,index)
      elsif (axis=="y")
        if (gp.rank == 1) then
          value = ymd2time(unit_0, value)
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/360*365.25 if @OPT_calendar360 && flag_0 && @OPT_exch
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/365*365.25 if @OPT_calendar365 && flag_0 && @OPT_exch
        elsif (@OPT_exch) then
          value = ymd2time(unit_0, value)
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/360*365.25 if @OPT_calendar360 && flag_0
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/365*365.25 if @OPT_calendar365 && flag_0
        else
          value = ymd2time(unit_1, value)
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/360*365.25 if @OPT_calendar360 && flag_1
          value = DCL.sgpget("UYMIN") + (value - DCL.sgpget("UYMIN"))/365*365.25 if @OPT_calendar365 && flag_1
        end
        DCL::sgplzu([DCL.sgpget("UXMIN"),DCL.sgpget("UXMAX")],[value,value],type,index)
      elsif (axis=="dx")
        x = DCL.sgpget("UXMIN"); xmax = DCL.sgpget("UXMAX"); dx = value.to_f
        x = x + dx
        while (x < xmax)
          DCL::sgplzu([x,x],[DCL.sgpget("UYMIN"),DCL.sgpget("UYMAX")],type,index)
          x = x + dx
        end
      elsif (axis=="dy")
        y = DCL.sgpget("UYMIN"); ymax = DCL.sgpget("UYMAX"); dy = value.to_f
        y = y + dy
        while (y < ymax)
          DCL::sgplzu([DCL.sgpget("UXMIN"),DCL.sgpget("UXMAX")],[y,y],type,index)
          y = y + dy
        end
      end
      }
  end

  def maskshading(gp)
    ary = @OPT_maskshading.split(","); newary = []; newary << ary[0]
    if (ary[1]) then
      (ary.length-1).times{|i|
        if (/^[^\d-]/ =~ ary[i+1] )
          newary[0] = newary[0]+","+ary[i+1]
        else
          newary << ary[i+1]
        end
      }
    end
    if (/\D/ =~ newary[0]) then
      mask = open_gturl_wildcard(newary[0])[0]
    else
      mask = gp
    end
    rlev = (newary[1]||0); ipat = (newary[2]||1515)
    DCL.uepset("rlev",rlev); DCL.uepset("ipat",ipat)
    GGraph.tone(mask,false,"ltone"=>false, "transpose"=>@OPT_exch)
  end

  def line_fill(gp)
    if (gp.class == GPhys) then
      ci = (@OPT_index || 1.0)/10; wi = (@OPT_index || 1.0)%10; pi = (@Overplot+1)%7
      wi = 1
      fi = ci*1000+pi*100+wi*10+4; #fi = ci*1000+999 #べた塗り
      DCL.uusarp(fi, fi)
      xmin = DCL.sgpget("UXMIN"); ymin = DCL.sgpget("UYMIN")
      xmax = DCL.sgpget("UXMAX"); ymax = DCL.sgpget("UYMAX")
      if (xmin*xmax < 0) then xline = 0.0
      elsif (xmax < 0) then xline = xmax
      else xline = xmin  end
      if (ymin*ymax < 0) then yline = 0.0
      elsif (ymax < 0) then yline = ymax
      else yline = ymin end
      DCL.uvdif(gp.axis(0).to_gphys.val,gp.val,gp.val*0+yline) unless @OPT_exch
      DCL.uhdif(gp.val,gp.val*0+xline,gp.axis(0).to_gphys.val) if @OPT_exch
    elsif (gp.class == Array) then
      ci = (@OPT_index || 1.0)/10; wi = (@OPT_index || 1.0)%10; pi = (@Overplot+1)%7
      wi = 1
      fi = ci*1000+999   #pi*100+wi*10+4; #fi = ci*1000+999 #べた塗り
      DCL.uusarp(fi, fi)
      gp0 = gp[0]; gp1 = gp[1]
      DCL.uvdif(gp0.axis(0).to_gphys.val,gp0.val,gp1.val) unless @OPT_exch
      DCL.uhdif(gp0.val,gp1.val,gp0.axis(0).to_gphys.val) if @OPT_exch
    end
  end


  def set_vpsize( default_vp, aspect=2.0 )

    if (@OPT_panelfit) then
      vxmin = DCL.sgpget("vxmin");    vxmax = DCL.sgpget("vxmax")
      vymin = DCL.sgpget("vymin");    vymax = DCL.sgpget("vymax")
      if (@OPT_sldiv) then
        vaspect = @OPT_sldiv.split(",")[1].to_f/@OPT_sldiv.split(",")[2].to_f
        vymax = vymax*vaspect if vaspect < 1
        vxmax = vxmax/vaspect if vaspect > 1
        vxymax = [vxmax,vymax].max
        vxmax = vxmax/vxymax
        vymax = vymax/vxymax
      end
      unless @OPT_panelfit == "" then
        fact = @OPT_panelfit.to_f
        vxmargin = (vxmax - vxmin)*(1 - fact)*0.5
        vymargin = (vymax - vymin)*(1 - fact)*0.5
        vxmin = vxmin+vxmargin; vxmax = vxmax-vxmargin
        vymin = vymin+vymargin; vymax = vymax-vymargin
      end
      return [vxmin, vxmax, vymin, vymax]
    end

    aspect = aspect.to_f
    iw = DCL.swpget('iwidth'); ih = DCL.swpget('iheight')
    margin = 0.15; yoffset = -0.00
    if (@OPT_sldiv)
     xn = @OPT_sldiv.split(",")[1].to_f; yn = @OPT_sldiv.split(",")[2].to_f
     if (@OPT_sldiv.split(",")[3] == 'full')
       margin = 0.025; yoffset = -0.025
     end
    else
     xn = 1.0; yn = 1.0
    end
    if (iw/xn >= ih/yn) then # 1パネルの描画領域のviewportの最大値
      vxlimit = 1.0; vylimit = (ih/yn)/(iw/xn)
      if (aspect >= 1 && vylimit >= 1.0/aspect) then
        vxmin = 0.0 + margin; vxmax = 1.0 - margin
        vymin = (vylimit/2) - 0.5*(vxmax-vxmin)/aspect + yoffset
        vymax = (vylimit/2) + 0.5*(vxmax-vxmin)/aspect + yoffset
      else
        vymin = 0.0 + (margin + yoffset)*vylimit ; vymax = vylimit + (-margin + yoffset)*vylimit
        vxmin = 0.5 - 0.5*(vymax-vymin)*aspect
        vxmax = 0.5 + 0.5*(vymax-vymin)*aspect
      end
    else
      vxlimit = (iw/xn)/(ih/yn); vylimit = 1.0
      if (aspect < 1.0 && vxlimit >= 1.0*aspect) then
        vymin = 0.0 + margin + yoffset; vymax = 1.0 - margin + yoffset
        vxmin = (vxlimit/2) - 0.5*(vymax-vymin)*aspect
        vxmax = (vxlimit/2) + 0.5*(vymax-vymin)*aspect
      else
        vxmin = 0.0 + margin*vxlimit ; vxmax = vxlimit - margin*vxlimit
        vymin = 0.5 - 0.5*(vxmax-vxmin)/aspect + yoffset*vxlimit
        vymax = 0.5 + 0.5*(vxmax-vxmin)/aspect + yoffset*vxlimit
      end
    end
    return [vxmin,vxmax,vymin,vymax]
  end

  def title(string, size=2)
    v = DCL::sgqvpt
    vx = (v[0]+v[1])/2
    if (size == 2) then
      txhgt = DCL.uzpget('rsizec2')     # rsizec2: larger text
      vy = v[3] + (0.5 + DCL.uzpget('pad1')) * txhgt * size
    else
      txhgt = DCL.uzpget('rsizec1')     # rsizec1: small text
      vy = v[3] + (0.5 + DCL.uzpget('pad1')) * txhgt * size
    end

    lw = 3 unless @OPT_miscindex
    lw = @OPT_miscindex.to_i if @OPT_miscindex
    DCL::sgtxzr(vx , vy, string, txhgt, 0, 0, lw)
  #      DCL.uxmttl('t',string,0.0)

    # top-left-corner
    if (@OPT_sldiv) then
      @prev_panelname = nil unless @prev_panelname
      panelnum = (DCL.sgpget("nframe")-1) % (@OPT_sldiv.split(",")[1].to_i * @OPT_sldiv.split(",")[2].to_i)
      panelname = "("+@alphabet[panelnum]+")"
      DCL::sgtxzr(v[0]*0.5 , v[3]*1.075, panelname, txhgt*1.35, 0, -1, lw) unless @prev_panelname == panelname
      @prev_panelname = panelname
    end
    #p DCL.sgpget("nframe"), DCL.sgpget("npage"), DCL.sgpget("nlevel")
  end

  # 1次元のgpをline plotしたときにstep状になるように変更する。
  # 元々の格子点を x[n]、データを v[n] としたとき、いかのような配列にする。長さは2倍になる。
  # 軸：x[0], x[1], x[1], x[2], x[2], ..., x[N], x[N], x[N]+dx
  # 値：v[0], v[0], v[1], v[1], v[2], ..., v[N-1], v[N], v[N]
  def step_shape(gp)
    axis = gp.axis(0).to_gphys; aval = axis.val; gpval = gp.val; alen = aval.length
    axna = NArray.sfloat(2*alen)
    (alen-1).times{|i| axna[2*i]=aval[i]; axna[2*i+1]=aval[i+1]}
    axna[2*(alen-1)]=aval[(alen-1)]; axna[2*(alen-1)+1]=aval[(alen-1)]+(aval[(alen-1)]-aval[(alen-2)])
    vana = NArray.sfloat(2*alen)
    gpval.length.times{|i| vana[2*i]=gpval[i]; vana[2*i+1]=gpval[i]}
    xaxis = Axis.new
    xaxis.pos = VArray.new(axna,{"long_name"=>axis.get_att("long_name"),"units"=>axis.units.to_s},axis.name)
    grid = Grid.new(xaxis)
    varray = VArray.new(vana,{"long_name"=>gp.get_att("long_name"),"units"=>gp.units.to_s},gp.name)
    gp = GPhys.new(grid,varray)
    return gp
  end

  def draw_text(text, gary)
    vxmin = @VIEWPORT[0]; vxmax = @VIEWPORT[1]; vymin = @VIEWPORT[2]; vymax = @VIEWPORT[3]
    DCL.grfrm ; DCL.grstrn(1); DCL.grsvpt(vxmin,vxmax,vymin,vymax)
    DCL.grswnd(0,1,0,1); DCL.grstrf
    rsizel1 = DCL.uzpget("rsizel1"); charsize = 0.025; charindex = 3
    charfact = charsize/rsizel1; offset = 0.0

    # large horizontal color bar
    if (@textbox_lc) then
      if (gary[-1].rank >= 2) then
        GGraph::color_bar("top"=>true,"landscape"=>true,"title" => "["+gary[-1].units.to_s+"]",
                          "voff"=>-0.01, "vlength"=>(vxmax-vxmin)*0.80, "vwidth"=>(vymax-vymin)*0.1,
                          "charfact"=>charfact)
        offset = (vymax-vymin)*0.1 + charsize*2
      end

      # large legend of line plot
      if (gary[-1].rank == 1) then
        com_block = common_blocks(@sources[-@Overplot_max..-1])
        @Overplot_max.times{|i|
          if (i==0) then
            first = true; vy=vymax; vx=vxmin*1.1
          else
            first = false; vy=nil; vx=vxmin*1.1
          end
          if (@OPT_nocolor) then
            index = 3; type =@type_array[i]
          else
            index=@index_array[i]+2; type = 1
          end
          line = (!@OPT_mark) ? true : false
          size=charsize; dx=nil; mark_size=nil
          legend = @sources[-@Overplot_max+i]
          com_block.each{|s| legend.sub!(s,"").sub!("__","_") }
          legend = legend[1..-1] if legend[0] == "_"
          legend = legend[0..-2] if legend[-1] == "_"
          DCLExt::legend(legend,type,index,line,size,vx,dx,vy,first)
      }
      offset = (@Overplot_max)*charsize*1.5 + charsize
      end
    else offset = 0
    end
    DCL.sgtxzv(vxmin,vymax-offset,text,charsize,0,0,charindex) # draw text

  end

  def visualize_and_output(g)
    ## set title from @OPT_title_array
    if (@OPT_title_array && @Overplot == 1)
      @OPT_title = @title_array[0]
      @title_array.shift
    end

    ## set index from @index_array
    @OPT_index = @index_array[@Overplot-1] unless @OPT_nocolor

    ## set type from @type_array
    @OPT_type = @type_array[@Overplot-1] if @OPT_type_array

    # judge draw kind
    if (g.class == GPhys)
      unless @kind_of_fig
        gp_rank = (GPhys::VERSION >= "1.5.4") ? g.rank_len_ne_1 : g.rank
        if (@OPT_mark)
          @kind_of_fig = "mark"
        elsif (@OPT_histogram)
          @kind_of_fig = "histogram1D"
        elsif (@OPT_line || gp_rank == 1)
          @kind_of_fig = "line"
        elsif (@OPT_fullcolor)
          @kind_of_fig = "fullcolor"
        elsif (!@OPT_line && gp_rank >= 2) && !@OPT_noshade && @OPT_nocont
          @kind_of_fig = "nocont"
        elsif (!@OPT_line && gp_rank >= 2) && @OPT_noshade && !@OPT_nocont
          @kind_of_fig = "noshade"
        elsif (!@OPT_line && gp_rank >= 2) && !@OPT_noshade && !@OPT_nocont
          @kind_of_fig = "full"
        elsif (gp_rank == 0) then
          @kind_of_fig = "nofig"
        end
      end

      #g.axis(1).set_pos(VArray.new(NArray.to_na([1,2,3,4,5,6,7,8,9,10,11,12]), {"units"=>"month"}, "time") )

      # change the netcdf var name and attributes.
      if (@OPT_rename || @OPT_long_name || @OPT_unit) then
        g = g.copy if (g.data.file != nil || g.data.file.class == NArray)
        g.rename(@OPT_rename) if @OPT_rename
        g.set_att("long_name", @OPT_long_name) if @OPT_long_name
        g.set_att("long_name", @OPT_rename) if (!@OPT_long_name && @OPT_rename)
        g.units=@OPT_unit if @OPT_unit
      end

      ### output operated gphys object as NetCDF file
      if (@OPT_nc)
        if (@OPT_output_assoc_coords) then
          g = g.assoc_coord_gphys(g.assoccoordnames[0])
        end


        g = g.to_type(@OPT_ntype) if @OPT_ntype
        g = g.copy if (g.data.file != nil || g.data.file.class == NArray)
        if (@OPT_nc == ""); outfilename = "out.nc"; else; outfilename = @OPT_nc; end
        unless @outncfile then
          #NetCDF.creation_format=(NetCDF::NC_NETCDF4 | NetCDF::NC_CLASSIC_MODEL) # for netcdf ver 4 compression
          @outncfile=NetCDF.create(outfilename)
        end
        if (g.val.class == NArrayMiss) then # to set missing values to output netcdf files.
          if (g.val.get_mask.eq(0).sum > 0) then # skip if there are no missing value points.
            val = g.val
            rmiss = (val.to_na*val.get_mask.eq(0).to_f).max
            rmiss = (val.to_na*val.get_mask.eq(0).to_f).min if rmiss == 0.0 # get missing value
            g.set_att("missing_value", [rmiss])
          end
        end
        g.del_att("positive") if g.get_att("positive")
        GPhys::IO.write( @outncfile, g )
        # @outncfile.var(g.name).deflate(1,true) # netcdf ver 4 compression
        print "Operated GPhys object is written as #{outfilename}.\n" unless @OPT_silent
      end

      ### output operated gphys object as CSV file
      output_csv(g) if (@OPT_csv)
      ###


      if (@OPT_scatter or @OPT_color_scatter or @OPT_histogram2D or @OPT_cot or @OPT_rmap or @OPT_fullcolor == "2" or @OPT_vector)
        return g # do not draw single var plot; store vars in array.
      else
        draw(g, @kind_of_fig) unless (@OPT_scatter or @OPT_color_scatter or @OPT_nodraw)
        @prev_gp = g
        return g #if (@OPT_pry or @OPT_auto_pry or @OPT_nc4)
      end
    else
      unless @OPT_silent
        p g # if g is scalar, just print the value
      else
        # print g.val, ","
        p g.val
      end
      return g
    end
  end

  def check_dclopts
    indices = []
      ARGV.each_index{|i|
      if /^-(.)(.):(.*)=(.*)/ =~ ARGV[i]
        indices << i
      end
      }

    dclopts = []
      indices.reverse_each{|i|
      dclopts << ARGV[i]
      ARGV.slice!(i)
      }

    dclopts.each{|opt|
      pkg = opt[1..2]
      name = opt[4..-1].split('=')[0]
      value = opt[4..-1].split('=')[1]
      dcl_set_params(pkg,name,value)
    }

  end

  def dcl_set_params(pkg,name,value)
    set = 'stx'
    case name
    when /^c/i
      eval( "DCL.#{pkg}c#{set}(name,value.to_s)" )
    when /^l/i
      if /(.*)(T|t)(.*)/ =~ value
        eval( "DCL.#{pkg}l#{set}(name,true)" )
      else
      if /(.*)(F|f)(.*)/ =~ value
        eval( "DCL.#{pkg}l#{set}(name,false)" )
      else
        raise "value of logical parameter must include 't' or 'f'"
      end
      end
    when /^[i-n]/i
      eval( "DCL.#{pkg}i#{set}(name,value.to_i)" )
    else
      eval( "DCL.#{pkg}r#{set}(name,value.to_f)" )
    end

  end


end
