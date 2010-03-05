#!/usr/bin/env ruby

# name: make_isPCR.rb
# author: Oliver Alexander Hampton
# date: Sept 18, 2007
#
# input:  breakpoint primer3 output file
#         HG18 one-line DNA files [not command line input]
#
# output: lff annotations of Primer3 sequences
#         one file for breakpoints
#         outputs the PRIMER_SEQUENCE_ID for those breakpoints without 
#         passing primers
#
# description: checks the primer3 primers for making spurious amplicons
#
# NOTE: This program must be executed on the LINUX(woodstock) cluster or NED

require 'ftools'
require 'brl/util/util'
require 'brl/util/textFileUtil'


#  METHODS
###############################################################################

def processArguments()
  # We want to add all the prop_keys as potential command line options
  optsArray =  [  ['--BreakPointPrimer3',  '-b', GetoptLong::REQUIRED_ARGUMENT],
                  ['--help', '-h', GetoptLong::NO_ARGUMENT]
               ]
  progOpts = GetoptLong.new(*optsArray)
	optsHash = progOpts.to_hash
	missingOpts = progOpts.getMissingOptions()
	if(missingOpts.length != 0 || optsHash.key?("--help") || optsHash.empty?)
		puts "make_isPCR.rb -b BreakPointPrimer3.raw"
		puts "optional arguments:"
		puts "help                  --help          [-h]"
  end
	return optsHash
end


#  CLASSES
################################################################################

class LFF
	attr_accessor :lffclass, :name, :type, :subtype, :chrm, :start, :stop, :strand, :phase, :score, :qStart, :qStop, :attributes, :sequence, :comments, :amplicon_size
	def initialize(lffclass, name, type, subtype, chrm, start, stop, strand, attributes, sequence, amplicon_size)
		@lffclass = lffclass
		@name = name
		@type = type
		@subtype = subtype
		@chrm = chrm
		@start = start
		@stop = stop
		@strand = strand
		@phase = "0"
		@score = "0"
		@qStart = "."
		@qStop = "."
		@attributes = attributes
		@sequence = sequence
		@comments = ""
		@amplicon_size = amplicon_size
	end

	def write()
		return "#{@lffclass}\t#{@name}\t#{@type}\t#{@subtype}\t#{@chrm}\t#{@start.to_s}\t#{@stop.to_s}\t#{@strand}\t#{@phase}\t#{@score}\t#{@qStart}\t#{@qStop}\t#{@attributes}\t#{@sequence}\t#{@comments}"
	end
end


#  CONSTANTS
################################################################################

##PATH TO CHROMOSOME .nib FILES FOR THE GENOME USED IN PRIMER DESIGN
DIR = "/data/brl_workspace/brl/blastdb/Hs.GoldenPath/hg18/chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}.fa.nib"

##SAME AS make_breakpoint_regions.rb "SeqRead" VARIABLE
SeqRead = 1000

################################################################################

optHash = processArguments()
breakpoints = Hash.new{|hh,kk| hh[kk]=Hash.new{|mm,nn| mm[nn]=Hash.new{|jj,vv| jj[vv]=Hash.new{|pp,qq| pp[qq]=Hash.new} } } }
controls = Hash.new{|hh,kk| hh[kk]=Hash.new{|mm,nn| mm[nn]=Hash.new{|jj,vv| jj[vv]=Hash.new{|pp,qq| pp[qq]=Hash.new} } } }
breakpointPrimerNames = Hash.new{ |hh,kk| hh[kk]=0 } 
str_breakpointPrimerNames = "" 
controlPrimerNames = Hash.new{ |hh,kk| hh[kk]=0 } 

local_attributes = []
target_attributes = [] 
local_primer_seq = []  ## LEFT  PRIMER
target_primer_seq = [] ## RIGHT PRIMER
local_primer_start = []
local_primer_stop = []
target_primer_start = []
target_primer_stop = []
amplicon_size = []

# BREAK POINT PRIMERS
################################################################################

base_name = "#{optHash["--BreakPointPrimer3"]}"
base_name.gsub!(".txt_p3cmd_primer3.raw", "")

isPCR_input_file = "#{base_name}_isPcr_input.txt"
isPCR_input = File.open( isPCR_input_file, "w")

reader = BRL::Util::TextReader.new(optHash["--BreakPointPrimer3"])
reader.each("=\n"){ |rec|
	local_attributes.clear
	target_attributes.clear
	local_primer_seq.clear
	target_primer_seq.clear
	local_primer_start.clear
	local_primer_stop.clear
	target_primer_start.clear
	target_primer_stop.clear
	amplicon_size.clear

	name = nil
	type = nil
	#attr = nil
	local_chrm = nil
	local_start = nil
	local_stop = nil
	local_strand = nil
	target_chrm = nil
	target_start = nil
	target_stop = nil
	target_strand = nil
	mask = nil
	tier = nil
	num = nil

	rec.each{ |line|
		line.strip!
		if( line =~ /^PRIMER_SEQUENCE_ID=(.+)/ )
			arrID = $1.split(/\|/)
			name = arrID[0]
			type = arrID[1]
			local_chrm = arrID[2]
			local_start = arrID[3].to_i
			local_stop = arrID[4].to_i
			local_strand = arrID[5]
			target_chrm = arrID[6]
			target_start = arrID[7].to_i
			target_stop = arrID[8]
			target_strand = arrID[9]
			#attr = arrID[10]
			mask = arrID[10]
			tier = arrID[11]
			breakpointPrimerNames[name]
			str_breakpointPrimerNames << name
		elsif( line =~ /^PRIMER_PAIR_EXPLAIN=.+(\d+)$/ )
			num = $1.to_i
		elsif( line =~ /^PRIMER_LEFT(_\d)?_SEQUENCE=(\w+)$/ )
			reg_index = $1
			reg_seq = $2
			if( reg_index.nil? )
				index = 0
			else
				index = reg_index.gsub("_", "").to_i
			end
			local_primer_seq[index] = reg_seq
		elsif( line =~ /^PRIMER_RIGHT(_\d)?_SEQUENCE=(\w+)$/ )
			reg_index = $1
			reg_seq = $2
			if(  reg_index.nil? )
				index = 0
			else
				index = reg_index.gsub("_","").to_i
			end
			target_primer_seq[index] = reg_seq
		elsif( line =~ /^PRIMER_LEFT(_\d)?=(\d+),(\d+)/ )
			reg_index = $1
			reg_start = $2
			reg_length = $3
      if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end		
			if( local_strand == "+" )
				local_primer_start[index] = local_start + reg_start.to_i + 1
				local_primer_stop[index] = local_start + reg_start.to_i + reg_length.to_i
			elsif( local_strand == "-" )
				local_primer_start[index] = ( ( local_stop - reg_start.to_i ) - reg_length.to_i ) 
				local_primer_stop[index] =  ( local_stop - reg_start.to_i ) - 1
			else
				$stderr.puts local_strand
				$stderr.puts "FAILING LOCAL PRIMER START AND STOP CALCULATIONS"
			end
			### CREATE EDGE URL LINK IN TARGET ATTRIBUTES
			local_edge_url_start = local_primer_start[index] - ((local_primer_stop[index] - local_primer_start[index]) * 0.2).to_i
			local_edge_url_stop = local_primer_stop[index] + ((local_primer_stop[index] - local_primer_start[index]) * 0.2).to_i
			target_attributes[index] = "edge=http://www.genboree.org/java-bin/gbrowser.jsp?refSeqId=???&entryPointId=#{local_chrm}&from=#{local_edge_url_start}&to=#{local_edge_url_stop}; "
		elsif( line =~ /^PRIMER_RIGHT(_\d)?=(\d+),(\d+)/ )
			reg_index = $1
			reg_start = $2
			reg_length = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
			if( target_strand == "-" )
				target_primer_start[index] = ( target_start + ( reg_start.to_i - SeqRead ) ) - ( reg_length.to_i - 1 )
				target_primer_stop[index] = ( target_start + ( reg_start.to_i - SeqRead ) )
			elsif( target_strand == "+" )
        target_primer_start[index] = target_start + ( SeqRead - ( reg_start.to_i - SeqRead ) )
				target_primer_stop[index] = target_start + ( SeqRead - ( reg_start.to_i - SeqRead ) ) + ( reg_length.to_i - 1 )
			else
				$stderr.puts target_strand
				$stderr.puts "FAILING TARGET PRIMER START AND STOP CALCULATIONS"
			end
			### CREATE EDGE URL LINK IN LOCAL ATTRIBUTES
			target_edge_url_start = target_primer_start[index] - ((target_primer_stop[index] - target_primer_start[index]) * 0.2).to_i
			target_edge_url_stop = target_primer_stop[index] + ((target_primer_stop[index] - target_primer_start[index]) * 0.2).to_i
			local_attributes[index] = "edge=http://www.genboree.org/java-bin/gbrowser.jsp?refSeqId=???&entryPointId=#{target_chrm}&from=#{target_edge_url_start}&to=#{target_edge_url_stop}; "
		elsif( line =~ /^(PRIMER_LEFT)(_\d)?(_TM=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
			local_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_RIGHT)(_\d)?(_TM=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_LEFT)(_\d)?(_GC_PERCENT=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      local_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_RIGHT)(_\d)?(_GC_PERCENT=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_LEFT)(_\d)?(_SELF_ANY=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      local_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_RIGHT)(_\d)?(_SELF_ANY=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_LEFT)(_\d)?(_SELF_END=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      local_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_RIGHT)(_\d)?(_SELF_END=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_LEFT)(_\d)?(_END_STABILITY=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      local_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_RIGHT)(_\d)?(_END_STABILITY=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_PAIR)(_\d)?(_COMPL_ANY=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
			local_attributes[index] << "#{reg_label}#{reg_term}; "
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_PAIR)(_\d)?(_COMPL_END=[0-9.]+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
			local_attributes[index] << "#{reg_label}#{reg_term}; "
      target_attributes[index] << "#{reg_label}#{reg_term}; "
		elsif( line =~ /(PRIMER_PRODUCT_SIZE)(_\d)?(=\d+)/ )
			reg_label = $1
			reg_index = $2
			reg_term = $3
			if( reg_index.nil? )
        index = 0
      else
        index = reg_index.gsub("_","").to_i
      end
			amplicon_size[index] = reg_term.gsub("=", "").to_i
			local_attributes[index] << "#{reg_label}#{reg_term}; "
      target_attributes[index] << "#{reg_label}#{reg_term}; "
			local_attributes[index] << "NUM_PRIMERS_PASSED=#{num}; "
      target_attributes[index] << "NUM_PRIMERS_PASSED=#{num}; "
		end
	}	
	# Insert Primer set records into primer hash
	# Write Primer set sequences to in-silico pcr input file
	unless( num == 0 )
		if( num > 5 )
			num = 5
		end
		0.upto(num-1){ |ii|
			if( ii == 0 )
				lff_name = name
			else
				lff_name = name + "_#{ii}"
			end
			##local_attributes[ii] << "#{attr}; "
			##target_attributes[ii] << "#{attr}; "
			lff_left_primer = LFF.new("Primer", lff_name, type, "Primer-Span", local_chrm, local_primer_start[ii], local_primer_stop[ii], local_strand, local_attributes[ii], local_primer_seq[ii], amplicon_size[ii])
			lff_right_primer = LFF.new("Primer", lff_name, type, "Primer-Span", target_chrm, target_primer_start[ii], target_primer_stop[ii], target_strand, target_attributes[ii], target_primer_seq[ii], amplicon_size[ii])
			breakpoints[name][mask][tier]["left"][ii.to_s] = lff_left_primer
			breakpoints[name][mask][tier]["right"][ii.to_s] = lff_right_primer
      if( local_primer_seq[ii] != nil && target_primer_seq[ii] != nil )
        isPCR_input.puts "#{lff_name}|#{mask}|#{tier}\t#{local_primer_seq[ii]}\t#{target_primer_seq[ii]}"
      end
		}		
	end
}
isPCR_input.close

# Execute In-Silico PCR on BreakPoint Primers
isPCR_file = "#{base_name}_isPCR.out"
if( File.exist?(isPCR_file) )
	File.truncate(isPCR_file, 0)
end

Dir::glob(DIR).each{ |file_name|
	system( "isPcr #{file_name} #{isPCR_input_file} stdout -out=bed -maxSize=6000 >> #{isPCR_file}" ) 
}

# PARSE In-Silico BreakPoint Primers Output & Generate LFF Annotations for PASSING Primers
reader = BRL::Util::TextReader.new(isPCR_file)
reader.each{ |line|
	arrSplit = line.split(/\t/)
	arrID = arrSplit[3].split(/\|/)
	if( arrID[0] =~ /(_\d)?$/ )
		reg_index = $1
		if( reg_index.nil? )
			index = 0
		else
			index = reg_index.gsub("_", "").to_i
			arrID[0].chop!
			arrID[0].chop!
		end
	end
	breakpoints[arrID[0]][arrID[1]][arrID[2]]["left"].delete(index.to_s)
	breakpoints[arrID[0]][arrID[1]][arrID[2]]["right"].delete(index.to_s)
}

# CLEAN UP BREAK POINT HASH
breakpoints.each_key{ |name_key|
  breakpoints[name_key].each_key{ |mask_key|
    breakpoints[name_key][mask_key].each_key{ |tier_key|
      breakpoints[name_key][mask_key][tier_key].each_key{ |dir_key|
        breakpoints[name_key][mask_key][tier_key][dir_key].each_key{ |ii|
          if( breakpoints[name_key][mask_key][tier_key][dir_key][ii].amplicon_size == nil )
            breakpoints[name_key][mask_key][tier_key]["left"].delete(ii)
            breakpoints[name_key][mask_key][tier_key]["right"].delete(ii)
          end
        }
      }
    }
  }
}
breakpoints.each_key{ |name_key|
	breakpoints[name_key].each_key{ |mask_key|
		breakpoints[name_key][mask_key].each_key{ |tier_key|
			breakpoints[name_key][mask_key][tier_key].each_key{ |dir_key|
				if( breakpoints[name_key][mask_key][tier_key][dir_key].empty? )
					breakpoints[name_key][mask_key][tier_key].delete("left")
          breakpoints[name_key][mask_key][tier_key].delete("right")
				end
			}
		}
	}
}
breakpoints.each_key{ |name_key|
	breakpoints[name_key].each_key{ |mask_key|
		breakpoints[name_key][mask_key].each_key{ |tier_key|
			if( breakpoints[name_key][mask_key][tier_key].empty? )
				breakpoints[name_key][mask_key].delete(tier_key)
			end
		}
	}
}
breakpoints.each_key{ |name_key|
	breakpoints[name_key].each_key{ |mask_key|
	if( breakpoints[name_key][mask_key].empty? )
		breakpoints[name_key].delete(mask_key)
	end
	}
}
breakpoints.each_key{ |name_key|
	if( breakpoints[name_key].empty? )
		breakpoints.delete(name_key)
	end
}
# PRINT OUT LFF PRIMER ANNOTATIONS
isPCR_lff_file = "#{base_name}_isPCR_primers.lff"
isPCR_lff = File.open( isPCR_lff_file, "w")
isPCR_fail_file = "#{base_name}_isPCR_primer_FAILURE.txt"
isPCR_fail = File.open( isPCR_fail_file, "w")
mask1 = "RepeatMasked"
mask2 = "NonRepeatMasked"
tier1 = "Tier1"
tier2 = "Tier2"

breakpointPrimerNames.each_key{ |name_key|
	if( breakpoints.key?(name_key) )
		if( breakpoints[name_key].key?(mask1) )
			if( breakpoints[name_key][mask1].key?(tier1) )
        arrLFF = breakpoints[name_key][mask1][tier1]["left"].values.sort{ |aa,bb|
              retVal = aa.amplicon_size <=> bb.amplicon_size
              retVal = breakpoints[name_key][mask1][tier1]["left"].index(aa) <=>  breakpoints[name_key][mask1][tier1]["left"].index(bb)
              retVal
        }
        largest_key = breakpoints[name_key][mask1][tier1]["left"].index(arrLFF.first)
        breakpoints[name_key][mask1][tier1]["left"][largest_key].attributes << "PRIMER_MASK=#{mask1}; PRIMER_TIER=#{tier1}; "
        isPCR_lff.puts breakpoints[name_key][mask1][tier1]["left"][largest_key].write()
        breakpoints[name_key][mask1][tier1]["right"][largest_key].attributes << "PRIMER_MASK=#{mask1}; PRIMER_TIER=#{tier1}; "
        isPCR_lff.puts breakpoints[name_key][mask1][tier1]["right"][largest_key].write()
      elsif( breakpoints[name_key][mask1].key?(tier2) )
        arrLFF = breakpoints[name_key][mask1][tier2]["left"].values.sort{ |aa,bb|
              retVal = aa.amplicon_size <=> bb.amplicon_size
              retVal = breakpoints[name_key][mask1][tier2]["left"].index(aa) <=>  breakpoints[name_key][mask1][tier2]["left"].index(bb)
              retVal
        }
        largest_key = breakpoints[name_key][mask1][tier2]["left"].index(arrLFF.first)
        breakpoints[name_key][mask1][tier2]["left"][largest_key].attributes << "PRIMER_MASK=#{mask1}; PRIMER_TIER=#{tier2}; "
        isPCR_lff.puts breakpoints[name_key][mask1][tier2]["left"][largest_key].write()
        breakpoints[name_key][mask1][tier2]["right"][largest_key].attributes << "PRIMER_MASK=#{mask1}; PRIMER_TIER=#{tier2}; "
        isPCR_lff.puts breakpoints[name_key][mask1][tier2]["right"][largest_key].write()
			end
		elsif( breakpoints[name_key].key?(mask2) )
			if( breakpoints[name_key][mask2].key?(tier1) )
        arrLFF = breakpoints[name_key][mask2][tier1]["left"].values.sort{ |aa,bb|
              retVal = aa.amplicon_size <=> bb.amplicon_size
              retVal = breakpoints[name_key][mask2][tier1]["left"].index(aa) <=>  breakpoints[name_key][mask2][tier1]["left"].index(bb)
              retVal
        }
        largest_key = breakpoints[name_key][mask2][tier1]["left"].index(arrLFF.first)
        breakpoints[name_key][mask2][tier1]["left"][largest_key].attributes << "PRIMER_MASK=#{mask2}; PRIMER_TIER=#{tier1}; "
        isPCR_lff.puts breakpoints[name_key][mask2][tier1]["left"][largest_key].write()
        breakpoints[name_key][mask2][tier1]["right"][largest_key].attributes << "PRIMER_MASK=#{mask2}; PRIMER_TIER=#{tier1}; "
        isPCR_lff.puts breakpoints[name_key][mask2][tier1]["right"][largest_key].write()
			elsif( breakpoints[name_key][mask2].key?(tier2) )
        arrLFF = breakpoints[name_key][mask2][tier2]["left"].values.sort{ |aa,bb|
              retVal = aa.amplicon_size <=> bb.amplicon_size
              retVal = breakpoints[name_key][mask2][tier2]["left"].index(aa) <=>  breakpoints[name_key][mask2][tier2]["left"].index(bb)
              retVal
        }
        largest_key = breakpoints[name_key][mask2][tier2]["left"].index(arrLFF.first)
        breakpoints[name_key][mask2][tier2]["left"][largest_key].attributes << "PRIMER_MASK=#{mask2}; PRIMER_TIER=#{tier2}; "
        isPCR_lff.puts breakpoints[name_key][mask2][tier2]["left"][largest_key].write()
        breakpoints[name_key][mask2][tier2]["right"][largest_key].attributes << "PRIMER_MASK=#{mask2}; PRIMER_TIER=#{tier2}; "
        isPCR_lff.puts breakpoints[name_key][mask2][tier2]["right"][largest_key].write()
			end
		end
	else
		isPCR_fail.puts name_key
	end
}
isPCR_lff.close
isPCR_fail.close

breakpointPrimerNames.clear
breakpoints.clear


exit 0
