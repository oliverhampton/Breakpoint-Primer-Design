#!/usr/bin/env ruby

#name: make_p3cmd.rb
#author: Oliver Alexander Hampton
#date: Sept 18, 2009
#
#input:  BreakPointRegions.txt
#
#output: Primer3 cmd file
#
#description: makes the Primer3 command file for primer3 execution
#             program assumes "primer3_core" is an executable command
#             for Primer3 on current system

require 'ftools'
require 'brl/util/util'
require 'brl/util/textFileUtil'


#  METHODS
###############################################################################

def processArguments()
  # We want to add all the prop_keys as potential command line options
  optsArray =  [  ['--BreakPointRegion',  '-b', GetoptLong::OPTIONAL_ARGUMENT],
                   ['--help', '-h', GetoptLong::NO_ARGUMENT]
               ]
  progOpts = GetoptLong.new(*optsArray)
	optsHash = progOpts.to_hash
	missingOpts = progOpts.getMissingOptions()
	if(missingOpts.length != 0 || optsHash.key?("--help") || optsHash.empty?)
		puts "make_p3cmd.rb"
		puts "optional arguments:"
		puts "BreakPointRegions.txt --BreakPointRegion   [-b]"
		puts "help                  --help               [-h]"
  end
	return optsHash
end


#  CONSTANTS
################################################################################


PARAMS_TIER1 = "PRIMER_MIN_SIZE=20
PRIMER_OPT_SIZE=25
PRIMER_MAX_SIZE=27
PRIMER_MIN_TM=60.0
PRIMER_OPT_TM=60.0
PRIMER_MAX_TM=65.0
PRIMER_PRODUCT_SIZE_RANGE=50-200 50-400 50-600 50-800 50-1000
PRIMER_MAX_POLY_X=3
PRIMER_SELF_END=3.00
PRIMER_SELF_ANY=4.00
PRIMER_GC_CLAMP=1
PRIMER_MIN_GC=45.0
PRIMER_MAX_GC=60.0
PRIMER_NUM_RETURN=5
PRIMER_NUM_NS_ACCEPTED=0
PRIMER_MAX_DIFF_TM=5.0
PRIMER_EXPLAIN_FLAG=1"

PARAMS_TIER2 = "PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=22
PRIMER_MAX_SIZE=26
PRIMER_MIN_TM=60.0
PRIMER_OPT_TM=64.0
PRIMER_MAX_TM=68.0
PRIMER_PRODUCT_SIZE_RANGE=50-200 50-400 50-600 50-800 50-1000
PRIMER_MAX_POLY_X=3
PRIMER_SELF_END=3.00
PRIMER_SELF_ANY=4.00
PRIMER_GC_CLAMP=0
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=70.0
PRIMER_NUM_RETURN=5
PRIMER_NUM_NS_ACCEPTED=0
PRIMER_MAX_DIFF_TM=10.0
PRIMER_EXPLAIN_FLAG=1"


################################################################################

optHash = processArguments()

definition = ""
sequence = ""
target = ""

input_label = "#{optHash["--BreakPointRegion"]}_p3cmd"
input_label.gsub!("_p3cmd.txt", "")
breakpoint_primer3_file = input_label + ".txt"
 
breakpoint_out = File.open( breakpoint_primer3_file, "w" )

reader = BRL::Util::TextReader.new(optHash["--BreakPointRegion"])
reader.each("=\n"){ |rec|
	rec.each{ |line|
		if( line =~ /^PRIMER_SEQUENCE_ID=/ )
			definition = line.chomp!
		elsif( line =~ /^SEQUENCE=/ )
			sequence = line.chomp!
		elsif( line =~ /^TARGET=/ )
			target = line.chomp!
		end
	}
	breakpoint_out.puts definition + "|Tier1"
	breakpoint_out.puts sequence
	breakpoint_out.puts PARAMS_TIER1
	breakpoint_out.puts target
	breakpoint_out.puts "="
	breakpoint_out.puts definition + "|Tier2"
	breakpoint_out.puts sequence
	breakpoint_out.puts PARAMS_TIER2
	breakpoint_out.puts target
	breakpoint_out.puts "="
}
breakpoint_out.close

# Execute primer3 by cat'ing into primer3_core

cleanOutput = "#{input_label}_primer3"

p3cmdOK = system("cat '#{breakpoint_primer3_file}' | primer3_core > #{cleanOutput}.raw 2> #{cleanOutput}.err")

# Check command succeeded:
if(p3cmdOK)
	
else
	puts "PRIMER3 EXECUTION ERROR : cat '#{breakpoint_primer3_file}' | primer3_core > #{cleanOutput}.raw 2> #{cleanOutput}.primer3.err"
end


exit 0
