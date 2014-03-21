#!/hgsc_software/ruby/latest/bin/ruby 

#name: make_p3cmd.rb
#author: Oliver Alexander Hampton
#date: Mar. 13, 2014
#
#input:  BreakPointRegions.txt
#
#output: Primer3 cmd file
#
#description: makes the Primer3 command file assuming primer3 (libprimer3 release 2.3.5)
#             "primer3_core" is an executable command, full path is set in #CONSTANTS
#             as is the full path to primer3_config for the varialbe THERMODYNAMIC_PARAMETERS_PATH   
#             Tier1 and Tier2 primer design parameters are optimized for Sanger Sequencing. 

require 'getoptlong'
require 'fileutils'
require 'rubygems'

#  METHODS
###############################################################################

def processArguments()
  optsArray =  [  ['--BreakPointRegion',  '-b', GetoptLong::REQUIRED_ARGUMENT],
                  ['--help', '-h', GetoptLong::NO_ARGUMENT]  ]

  progOpts = GetoptLong.new(*optsArray)
  missingOpts = Array.new
  optsArray.each{ |aa|
    if( aa[2] == 1 )
      missingOpts.push( aa[0] )
    end
  }

  optsHash = Hash.new{|hh,kk| hh[kk]=""}
  begin
    progOpts.each do |opt,arg|
      optsHash[opt] = arg
      missingOpts.delete( opt )
    end
  end

  if(missingOpts.length != 0 || optsHash.key?("--help") || optsHash.empty?)
    puts "make_p3cmd.rb -b BreakpointRegion.txt "
    puts "optional arguments:"
    puts "help                  --help          [-h]"
  end
  return optsHash
end


#  CONSTANTS
################################################################################

PRIMER3 = "/hgsc_software/primer3/latest/src/primer3_core"
THERMODYNAMIC_PARAMETERS_PATH = "/hgsc_software/primer3/latest/src/primer3_config/"

PARAMS_TIER1 = "PRIMER_MIN_SIZE=20
PRIMER_OPT_SIZE=25
PRIMER_MAX_SIZE=27
PRIMER_MIN_TM=60.0
PRIMER_OPT_TM=60.0
PRIMER_MAX_TM=65.0
PRIMER_PRODUCT_SIZE_RANGE=200-500 200-750 200-1000 100-1000
PRIMER_MAX_POLY_X=3
PRIMER_MAX_SELF_END=3.00
PRIMER_MAX_SELF_ANY=4.00
PRIMER_GC_CLAMP=1
PRIMER_MIN_GC=45.0
PRIMER_MAX_GC=60.0
PRIMER_NUM_RETURN=5
PRIMER_MAX_NUM_NS_ACCEPTED=0
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_EXPLAIN_FLAG=1"

PARAMS_TIER2 = "PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=22
PRIMER_MAX_SIZE=26
PRIMER_MIN_TM=60.0
PRIMER_OPT_TM=64.0
PRIMER_MAX_TM=68.0
PRIMER_PRODUCT_SIZE_RANGE=200-500 200-750 200-1000 100-1000
PRIMER_MAX_POLY_X=3
PRIMER_MAX_SELF_END=3.00
PRIMER_MAX_SELF_ANY=4.00
PRIMER_GC_CLAMP=0
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=70.0
PRIMER_NUM_RETURN=5
PRIMER_MAX_NUM_NS_ACCEPTED=0
PRIMER_PAIR_MAX_DIFF_TM=10.0
PRIMER_EXPLAIN_FLAG=1"

################################################################################

optHash = processArguments()

definition = ""
sequence = ""
target = ""

input_label = "#{optHash["--BreakPointRegion"]}_p3cmd"
input_label.gsub!(".txt", "")
breakpoint_primer3_file = input_label + ".txt"
 
breakpoint_out = File.open( breakpoint_primer3_file, "w" )

reader = File.open( optHash["--BreakPointRegion"] )
reader.each("=\n"){ |rec|
  rec.each_line{ |line|
    if( line =~ /^SEQUENCE_ID=/ )
      definition = line.chomp!
    elsif( line =~ /^SEQUENCE_TEMPLATE=/ )
      sequence = line.chomp!
    elsif( line =~ /^SEQUENCE_TARGET=/ )
      target = line.chomp!
    end
  }
  breakpoint_out.puts definition + "|Tier1"
  breakpoint_out.puts sequence
  breakpoint_out.puts PARAMS_TIER1
  breakpoint_out.puts target
  breakpoint_out.puts "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=#{THERMODYNAMIC_PARAMETERS_PATH}"
  breakpoint_out.puts "="
  breakpoint_out.puts definition + "|Tier2"
  breakpoint_out.puts sequence
  breakpoint_out.puts PARAMS_TIER2
  breakpoint_out.puts target
  breakpoint_out.puts "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=#{THERMODYNAMIC_PARAMETERS_PATH}"
  breakpoint_out.puts "="
}
breakpoint_out.close
reader.close

# Execute primer3 by cat'ing into primer3_core

cleanOutput = "#{input_label}_primer3"

p3cmdOK = system("cat '#{breakpoint_primer3_file}' | #{PRIMER3} > #{cleanOutput}.raw 2> #{cleanOutput}.err")

# Check command succeeded:
if(p3cmdOK)
	
else
  puts "PRIMER3 EXECUTION ERROR : cat '#{breakpoint_primer3_file}' | #{PRIMER3} > #{cleanOutput}.raw 2> #{cleanOutput}.primer3.err"
end


exit 0
