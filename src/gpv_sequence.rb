#--
# =DESCRIPTION
#  Scripts to execute gpv sequentially.
# =AUTHOR
#  Hiroki Kashimura
#++

module GPVSequence
  module_function
  def p_e(command)
    print "  "+command+"\n"
    print "    "+`#{command}`
  end 


  # SCALE-GMのnetcdf出力の単位を補完する 
  def scale_units 
    p_e "gpv prs.nc --edit_ncatt prs:units:Pa"
    p_e "gpv ps.nc --edit_ncatt ps:units:Pa"    
    p_e "gpv u.nc --edit_ncatt u:units:m.s-1"
    p_e "gpv v.nc --edit_ncatt v:units:m.s-1"
    p_e "gpv w.nc --edit_ncatt w:units:m.s-1"
    p_e "gpv t.nc --edit_ncatt t:units:K"
    p_e "gpv rho.nc --edit_ncatt rho:units:kg.m-3"
  end

end
