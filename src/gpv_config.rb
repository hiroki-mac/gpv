#--
# =DESCRIPTION
#  Config part of gpv.
# =AUTHOR
#  Hiroki Kashimura
#++

HOMEDIR = File.expand_path('~')

require "getoptlong"        # for option_parse
begin
  require "numru/narray"
rescue LoadError
  require "narray"
end

require "numru/fftw3"
require "numru/gphys"

require "numru/ggraph"      # ggraph library
require "numru/gphys/gpcommon"
require "numru/ganalysis"
require "csv"
require "pry"
require "ruby-progressbar"
require "parallel"

include NumRu
include NMath

begin
  # require "RubyFits"
  # include Fits
rescue LoadError
end

# numru-narrayかどうか判定する
begin
  NumRu::NArray
  NArrayType = "bigmem"
rescue NameError
  NArrayType = "standard"
end
#####################################################
## Default param.
VIEWPORT = [0.15,0.85,0.125,0.6] #gpview default
D2R = PI/180.0
R2D = 180.0/PI
URLfmt = "path[@|/]varname[,dimname=pos1[:pos2[:thinning_intv]][,dimname=...]]"
AX = "ax"

## get DCL version number
begin
  DCLVERNUM = `cdclconfig --dclvernum`.to_i
rescue
  DCLVERNUM = `dclconfig --dclvernum`.to_i
end

## judge OS type
host_os = RbConfig::CONFIG['host_os']
case host_os
when /darwin|mac os/
  OS_NAME = "macOS"
  EXT_PATH_DIR = "/opt/local/bin/"
when /linux/
  OS_NAME = "linux"
  EXT_PATH_DIR = "/usr/bin/"
else
  OS_NAME = "unknown"
  EXT_PATH_DIR = "/usr/bin/"
end

## set PATH of external commands
## set nil if the coomand is unavailable.
CMD_SSED    = EXT_PATH_DIR + "ssed"
CMD_CONVERT = EXT_PATH_DIR + "convert"
#CMD_FFMPEG  = EXT_PATH_DIR + "ffmpeg"
CMD_FFMPEG  = "ffmpeg"
CMD_MOGRIFY = EXT_PATH_DIR + "mogrify"
CMD_EXIFTOOL= EXT_PATH_DIR + "exiftool"
CMD_CPDF    = "cpdf"
CMD_PDFCROP = "/Library/TeX/texbin/pdfcrop"
CMD_RM      = "/bin/rm"
CMD_MV      = "/bin/mv"
