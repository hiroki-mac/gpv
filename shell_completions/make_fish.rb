# Run: ruby make_fish.rb > gpv.fish 
# Then: mv gpv.fish ~/.config/fish/completions/.
require "getoptlong"        # for option_parse
require "../src/gpv_options.rb"
require "pry"

OPTIONS.each{|opt|
  opt[0].slice!("--")
  opt[-1].tr!("\n", " ")
  opt[-1].squeeze!(" ")
  print "complete -c gpv -l #{opt[0]} -d \"#{opt[-2]} | #{opt[-1]} \"\n"
  }
