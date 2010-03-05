#!/usr/bin/env ruby

#name: make_mcf7_breakpoint_regions.rb
#author: Oliver Alexander Hampton
#date: Sept 16, 2009
#
#input:  breakpoint or structural variation lff file
#        HG18 one-line DNA files [not command line input]
#output: Primer3 input file
#
#description: makes the Primer3 input file for MCF7 breakpoint primer
#             designs.  Take in an .lff file of breakpoint annotations.
#             NOTE: program check .lff format: class must be "StructVar"
#             else program fails. Two identically named breakpoint annotations
#             are required per breakpoint, one for each side of the breakpoint.


require 'ftools'
require 'brl/util/util'
require 'brl/util/textFileUtil'


#  METHODS
###############################################################################

def processArguments()
  # We want to add all the prop_keys as potential command line options
  optsArray =  [  ['--lffFile',  '-f', GetoptLong::REQUIRED_ARGUMENT],
                  ['--help', '-h', GetoptLong::NO_ARGUMENT]
               ]
  progOpts = GetoptLong.new(*optsArray)
	optsHash = progOpts.to_hash
	missingOpts = progOpts.getMissingOptions()
	if(missingOpts.length != 0 || optsHash.key?("--help") || optsHash.empty?)
		puts "make_breakpoint_regions.rb -f Breakpoint.lff"
		puts "optional arguments:"
		puts "help                  --help          [-h]"
  end
	return optsHash
end


#  CONSTANTS
################################################################################

##Manually choose the Genome (Repeat Masked or Non Repeat Masked) in single line form
##Oliver Hampton's path is given as example
# REPEAT MASKED (TIER 1)
DIR = "/users/oh147993/brl/blastdb/Hs.GoldenPath/hg18/repeatMasked/oneLineSeq/"
OneLine = ".oneLineSeq.masked.seq"
Label = "RepeatMasked"


# NON-REPEAT MASKED (TIER 2)
#DIR = "/users/oh147993/brl/blastdb/Hs.GoldenPath/hg18/oneLineSeq/"
#OneLine = ".oneLineSeq.seq"
#Label = "NonRepeatMasked"


##Original Parameters
##Annotation Buffer lengths
## to pad breakpoint annotations so primers are not designed on breakpoint junction
## and sets contatenated annotation lengths to 1000bp each
TargetLength = 20
SeqRead = 1000

##CHM LIMIT Hash for HG18
$limitHash =  { "chr1" => 247249719,
  "chr10" => 135374737,
  "chr11" => 134452384,
  "chr12" => 132349534,
  "chr13" => 114142980,
  "chr14" => 106368585,
  "chr15" => 100338915,
  "chr16" => 88827254,
  "chr17" => 78774742,
  "chr18" => 76117153,
  "chr19" => 63811651,
  "chr2" => 242951149,
  "chr20" => 62435964,
  "chr21" => 46944323,
  "chr22" => 49691432,
  "chr3" => 199501827,
  "chr4" => 191273063,
  "chr5" => 180857866,
  "chr6" => 170899992,
  "chr7" => 158821424,
  "chr8" => 146274826,
  "chr9" => 140273252,
  "chrX" => 154913754,
  "chrY" => 57772954    }

################################################################################


#  CLASSES
################################################################################

class Edge
	attr_accessor :name, :type, :phase, :plate, :chrm, :start, :stop
	def initialize(theLine)
		parseEdgeLine(theLine)
	end

	def parseEdgeLine(theLine)
		arrSplit = theLine.split(/\t/)
		@name = arrSplit[1]
    @type = arrSplit[2]
		@chrm = arrSplit[4]
		@start = arrSplit[5].to_i
		@stop = arrSplit[6].to_i
    @phase = arrSplit[7]
	end

  def valid?()
    if( (@start - SeqRead) > 0 && (@stop + SeqRead) < $limitHash[@chrm] )
      return true
    else
      return false
    end
  end

end

################################################################################

class Join
	attr_accessor :name, :def_line, :sequence, :targets
	def initialize(name, def_line, sequence, targets)
		@name = name
		@def_line = def_line
		@sequence = sequence
		@targets = targets
	end
end

################################################################################


optHash = processArguments()
edges = Hash.new{|hh,kk| hh[kk]=Array.new}
joins = {}

primer3_file = "#{optHash["--lffFile"]}_#{Label}_regions.txt" 
primer3_file.gsub!(".lff", "")
breakpoint_out = File.open( primer3_file, "w" )


reader = BRL::Util::TextReader.new(optHash["--lffFile"])
reader.each{ |line|
	if( line =~ /^StructVar/ )
		line.chomp!
		arrSplit = line.split(/\t/)
		edge = Edge.new(line)
    if( edge.valid? )
      edges[arrSplit[1]].push(edge)
    else
      puts "StructVar NOT VALID: #{line}"
    end
  end
}


edges.each_key{ |kk|
if( edges[kk].size == 2 )
	if( edges[kk][0].phase != edges[kk][1].phase )
		if( edges[kk][0].phase == "+" )
			### SEQUENCE CONSTRUCTION
			def_line = "#{kk}|#{edges[kk][0].type}|#{edges[kk][0].chrm}|#{(edges[kk][0].stop - SeqRead).to_s}|#{edges[kk][0].stop.to_s}|#{edges[kk][0].phase}|#{edges[kk][1].chrm}|#{edges[kk][1].start.to_s}|#{(edges[kk][1].start + SeqRead).to_s}|#{edges[kk][1].phase}"
			local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
			local_file.seek(edges[kk][0].stop - SeqRead)
			local_seq = local_file.read(SeqRead)
			local_file.close
			target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
			target_file.seek(edges[kk][1].start - 1)
			target_seq = target_file.read(SeqRead)
			target_file.close
			sequence = local_seq + target_seq
			### TARGETS PARAM CONSTRUCTION
			targets_start = local_seq.size - (TargetLength/2)   
			targets = "#{targets_start.to_s},#{TargetLength.to_s}"
			join = Join.new(kk, def_line, sequence, targets)
			joins[kk] = join
		else
			### SEQUENCE CONSTRUCTION
			def_line = "#{kk}|#{edges[kk][1].type}|#{edges[kk][1].chrm}|#{(edges[kk][1].stop - SeqRead).to_s}|#{edges[kk][1].stop.to_s}|#{edges[kk][1].phase}|#{edges[kk][0].chrm}|#{edges[kk][0].start.to_s}|#{(edges[kk][0].start + SeqRead).to_s}|#{edges[kk][0].phase}"
			local_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
			local_file.seek(edges[kk][1].stop - SeqRead)
			local_seq = local_file.read(SeqRead)
			local_file.close
			target_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
			target_file.seek(edges[kk][0].start - 1)
			target_seq = target_file.read(SeqRead)
			target_file.close
			sequence = local_seq + target_seq
			### TARGETS PARAM CONSTRUCTION
			targets_start = local_seq.size - (TargetLength/2)   
			targets = "#{targets_start.to_s},#{TargetLength.to_s}"
			join = Join.new(kk, def_line, sequence, targets)
			joins[kk] = join
		end
	else
    if( edges[kk][0].phase == "+" )
			### SEQUENCE CONSTRUCTION
			def_line = "#{kk}|#{edges[kk][0].type}|#{edges[kk][0].chrm}|#{(edges[kk][0].stop - SeqRead).to_s}|#{edges[kk][0].stop.to_s}|#{edges[kk][0].phase}|#{edges[kk][1].chrm}|#{(edges[kk][1].stop - SeqRead).to_s}|#{edges[kk][1].stop.to_s}|#{edges[kk][1].phase}"
      local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
      local_file.seek(edges[kk][0].stop - SeqRead)
      local_seq = local_file.read(SeqRead)
      local_file.close
			target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
			target_file.seek(edges[kk][1].stop - SeqRead)
			target_seq = target_file.read(SeqRead)
			target_seq_RV = target_seq.tr("ACGTacgt","TGCAtgca").reverse
			target_file.close
			sequence = local_seq + target_seq_RV
			### TARGETS PARAM CONSTRUCTION
			targets_start = local_seq.size - (TargetLength/2)   
			targets = "#{targets_start.to_s},#{TargetLength.to_s}"
			join = Join.new(kk, def_line, sequence, targets)
			joins[kk] = join
		else
      ### SEQUENCE CONSTRUCTION
      def_line = "#{kk}|#{edges[kk][1].type}|#{edges[kk][1].chrm}|#{edges[kk][1].start.to_s}|#{(edges[kk][1].start + SeqRead).to_s}|#{edges[kk][1].phase}|#{edges[kk][0].chrm}|#{edges[kk][0].start.to_s}|#{(edges[kk][0].start + SeqRead).to_s}|#{edges[kk][0].phase}"
			local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
			local_file.seek(edges[kk][0].start - 1)
			local_seq = local_file.read(SeqRead)
			local_seq_RV = local_seq.tr("ACGTacgt","TGCAtgca").reverse
			local_file.close
			target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
			target_file.seek(edges[kk][1].start - 1)
			target_seq = target_file.read(SeqRead)
			target_file.close
			sequence = local_seq_RV + target_seq
			### TARGETS PARAM CONSTRUCTION
			targets_start = local_seq.size - (TargetLength/2)
			targets = "#{targets_start.to_s},#{TargetLength.to_s}"
			join = Join.new(kk, def_line, sequence, targets)
			joins[kk] = join
		end
	end
end
}
edges.clear

joins.each_key{ |kk|
	breakpoint_out.puts "PRIMER_SEQUENCE_ID=#{joins[kk].def_line}|#{Label}"
	breakpoint_out.puts "SEQUENCE=#{joins[kk].sequence}"
	breakpoint_out.puts "TARGET=#{joins[kk].targets}"
	breakpoint_out.puts"="
}
breakpoint_out.close


exit 0
