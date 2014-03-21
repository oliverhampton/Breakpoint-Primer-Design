#!/hgsc_software/ruby/latest/bin/ruby

#name: make_breakpoint_regions.rb
#author: Oliver Alexander Hampton
#date: Feb. 26, 2014
#
#input:       breakpoint or structural variation lff file
#             one-line DNA reference file path
#
#output:      Primer3 input file
#
#description: makes the Primer3 input file for breakpoint primer
#             designs.  Take in an .lff file of breakpoint annotations.
#             NOTE: Two identically named breakpoint annotations
#             are required per breakpoint, one for each side of the breakpoint.

require 'getoptlong'
require 'fileutils'
require 'rubygems'

#  METHODS
###############################################################################

def processArguments()
  optsArray =  [  ['--lffFile',  '-f', GetoptLong::REQUIRED_ARGUMENT],
                  ['--referenceDir', '-r', GetoptLong::REQUIRED_ARGUMENT],
                  ['--target_length', '-t', GetoptLong::REQUIRED_ARGUMENT],
                  ['--seq_read_length', '-s', GetoptLong::REQUIRED_ARGUMENT],
                  ['--reference_label', '-l', GetoptLong::REQUIRED_ARGUMENT],
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
    puts "make_breakpoint_regions.rb -f Breakpoint.lff -r Reference_Dir (path) -l Reference_Label -t Target_Length (int) -s Sequence_Read_Length (int)"
    puts "optional arguments:"
    puts "help                  --help          [-h]"
  end
  return optsHash
end


#  CONSTANTS
################################################################################

OneLine = ".fa.oneLine"

# SET PARAMETERS FROM INPUT
##TargetLength = 20
##SeqRead = 1000


#  CLASSES
################################################################################

class Edge
  attr_accessor :name, :type, :ori, :plate, :chrm, :start, :stop
  def initialize(theLine)
    parseEdgeLine(theLine)
  end
  
  def parseEdgeLine(theLine)
    arrSplit = theLine.split(/\t/)
    @name = arrSplit[1]
    @type = arrSplit[3]
    @chrm = arrSplit[4]
    @start = arrSplit[5].to_i
    @stop = arrSplit[6].to_i
    @ori = arrSplit[7]
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

LFF = optHash["--lffFile"]
DIR = optHash["--referenceDir"]
SeqRead = optHash["--seq_read_length"].to_i
TargetLength = optHash["--target_length"].to_i
Label = optHash["--reference_label"]

edges = Hash.new{|hh,kk| hh[kk]=Array.new}
joins = {}

primer3_file = "#{LFF}_#{Label}_regions.txt" 
breakpoint_out = File.open( primer3_file, "w" )

reader = File.open( optHash["--lffFile"] )
reader.each{ |line|
  line.chomp!
  arrSplit = line.split(/\t/)
  edge = Edge.new(line)
  edges[arrSplit[1]].push(edge)
}


flag = 0
edges.each_key{ |kk|
  if( edges[kk].size == 2 )
    flag = 0
    if( edges[kk][0].ori != edges[kk][1].ori )
      if( edges[kk][0].ori == "+" )
        ### SEQUENCE CONSTRUCTION
        def_line = "#{kk}|#{edges[kk][0].type}|#{edges[kk][0].chrm}|#{(edges[kk][0].stop - SeqRead).to_s}|#{edges[kk][0].stop.to_s}|#{edges[kk][0].ori}|#{edges[kk][1].chrm}|#{edges[kk][1].start.to_s}|#{(edges[kk][1].start + SeqRead).to_s}|#{edges[kk][1].ori}"
        local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")

        ###WORKING HERE
        if( (edges[kk][0].stop - SeqRead) >= 0 )
          local_file.seek(edges[kk][0].stop - SeqRead)
          local_seq = local_file.read(SeqRead)
          local_file.close
          if( local_seq.length != SeqRead )
            flag = 1
          end
        else
          flag = 1
        end

        target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
        target_file.seek(edges[kk][1].start - 1)
        target_seq = target_file.read(SeqRead)
        target_file.close
        if( target_seq.length != SeqRead )
          flag = 1
        end
        sequence = local_seq + target_seq
        ### TARGETS PARAM CONSTRUCTION
        targets_start = local_seq.size - (TargetLength/2)   
        targets = "#{targets_start.to_s},#{TargetLength.to_s}"
        join = Join.new(kk, def_line, sequence, targets)
        if( flag == 0 )
          joins[kk] = join
        end
      else
        ### SEQUENCE CONSTRUCTION
        def_line = "#{kk}|#{edges[kk][1].type}|#{edges[kk][1].chrm}|#{(edges[kk][1].stop - SeqRead).to_s}|#{edges[kk][1].stop.to_s}|#{edges[kk][1].ori}|#{edges[kk][0].chrm}|#{edges[kk][0].start.to_s}|#{(edges[kk][0].start + SeqRead).to_s}|#{edges[kk][0].ori}"
        local_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
        if( (edges[kk][1].stop - SeqRead) >= 0 )
          local_file.seek(edges[kk][1].stop - SeqRead)        
          local_seq = local_file.read(SeqRead)
          local_file.close
          if( local_seq.length != SeqRead )
            flag = 1
          end
        else
          flag = 1
        end
        
        target_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
        target_file.seek(edges[kk][0].start - 1)
        target_seq = target_file.read(SeqRead)
        target_file.close
        if( target_seq.length != SeqRead )
          flag = 1
        end
        sequence = local_seq + target_seq
        ### TARGETS PARAM CONSTRUCTION
        targets_start = local_seq.size - (TargetLength/2)   
        targets = "#{targets_start.to_s},#{TargetLength.to_s}"
        join = Join.new(kk, def_line, sequence, targets)
        if( flag == 0 )
          joins[kk] = join
        end
      end
    else
      if( edges[kk][0].ori == "+" )
        ### SEQUENCE CONSTRUCTION
        def_line = "#{kk}|#{edges[kk][0].type}|#{edges[kk][0].chrm}|#{(edges[kk][0].stop - SeqRead).to_s}|#{edges[kk][0].stop.to_s}|#{edges[kk][0].ori}|#{edges[kk][1].chrm}|#{(edges[kk][1].stop - SeqRead).to_s}|#{edges[kk][1].stop.to_s}|#{edges[kk][1].ori}"
        local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")
        if( (edges[kk][0].stop - SeqRead) >= 0 )
          local_file.seek(edges[kk][0].stop - SeqRead)
          local_seq = local_file.read(SeqRead)
          local_file.close
          if( local_seq.length != SeqRead )
            flag = 1
          end
        else
          flag = 1
        end

        target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
        target_file.seek(edges[kk][1].stop - SeqRead)
        target_seq = target_file.read(SeqRead)
        target_seq_RV = target_seq.tr("ACGTacgt","TGCAtgca").reverse
        target_file.close
        if( target_seq.length != SeqRead )
          flag = 1
        end
        sequence = local_seq + target_seq_RV
        ### TARGETS PARAM CONSTRUCTION
        targets_start = local_seq.size - (TargetLength/2)   
        targets = "#{targets_start.to_s},#{TargetLength.to_s}"
        join = Join.new(kk, def_line, sequence, targets)
        if( flag == 0 )	
          joins[kk] = join
        end
      else
      ### SEQUENCE CONSTRUCTION
        def_line = "#{kk}|#{edges[kk][1].type}|#{edges[kk][1].chrm}|#{edges[kk][1].start.to_s}|#{(edges[kk][1].start + SeqRead).to_s}|#{edges[kk][1].ori}|#{edges[kk][0].chrm}|#{edges[kk][0].start.to_s}|#{(edges[kk][0].start + SeqRead).to_s}|#{edges[kk][0].ori}"
        local_file = File.open("#{DIR}#{edges[kk][0].chrm}#{OneLine}")	
        local_file.seek(edges[kk][0].start - 1)
        local_seq = local_file.read(SeqRead)
        local_seq_RV = local_seq.tr("ACGTacgt","TGCAtgca").reverse
        local_file.close
        if( local_seq.length != SeqRead )
          flag = 1
        end

        target_file = File.open("#{DIR}#{edges[kk][1].chrm}#{OneLine}")
        target_file.seek(edges[kk][1].start - 1)
        target_seq = target_file.read(SeqRead)
        target_file.close
        if( target_seq.length != SeqRead )
          flag = 1
        end
        sequence = local_seq_RV + target_seq
        ### TARGETS PARAM CONSTRUCTION
        targets_start = local_seq.size - (TargetLength/2)
        targets = "#{targets_start.to_s},#{TargetLength.to_s}"
        join = Join.new(kk, def_line, sequence, targets)
        if( flag == 0 )	
          joins[kk] = join
        end

      end
    end
  end
}
edges.clear

joins.each_key{ |kk|
	breakpoint_out.puts "SEQUENCE_ID=#{joins[kk].def_line}|#{Label}"
	breakpoint_out.puts "SEQUENCE_TEMPLATE=#{joins[kk].sequence}"
	breakpoint_out.puts "SEQUENCE_TARGET=#{joins[kk].targets}"
	breakpoint_out.puts"="
}
breakpoint_out.close


exit 0
