#REQUIRED LIBRARIES
require 'open-uri'

#Given a file path, reads a local fasta file and returns a hash map with sequence ID and sequence
#as key-value pairs
def read_fasta(fasta_file)
    fasta_name = ""
    fasta_seqs = {}
    seq = ""
    File.open(fasta_file).each do |line|
        if !(line.nil?)
            if line.start_with? ">"
                seq = ""
                fasta_name = line.chomp.split(">")[1]
            elsif
                seq = seq + line.chomp
                fasta_seqs[fasta_name] = seq
            end
        end
    end
    return fasta_seqs
end

#Given a url pointing to a online fasta file, reads a fasta and returns a hash map with sequence ID and sequence
#as key-value pairs
def read_online_fasta(url)
    fasta_name = ""
    fasta_seqs = {}
    seq = ""
    open(url) do |file|
        file.read.split("\n").each do |line|
            if !(line.nil?)
                if line.start_with? ">"
                    seq = ""
                    fasta_name = line.split(">")[1]
                elsif
                    seq = seq + line
                    fasta_seqs[fasta_name] = seq
                end
            end
        end
    end
    return fasta_seqs
end

#Reads a Fast-Q file
def read_fastq(filename)
    sequences = []
    qualities = []
    counter = 1     #Support variable to know what lines are needed to save
    current = "seq"
    File.foreach(filename) do |line|
        if counter % 2 == 0
            if current == "seq"
                sequences = sequences.concat([line.chomp])
                current = "qual"
            else
                qualities = qualities.concat([line.chomp])
                current = "seq"
            end
        end
        counter += 1
    end

    return [sequences, qualities]
end

#Given an array of Phred33 qualities, it creates a histogram showing the overall quality reads 
def create_Qreads_hist_data(qualities)
    hist = [0] * 50
    qualities.each do |qual|
        qual.each_char do |phred|
            q = convert_phred33_to_q(phred)
            hist[q] += 1
        end
    end
    return hist
end


def create_GC_hist_data(reads)
    gc = [0] * 100
    totals = [0] * 100
    reads.each do |xread|
        for i in (0...xread.length)
            if xread[i] == "C" || xread[i] == "G"
                gc[i] += 1
            end
            totals[i] += 1
        end
    end

    for i in (0...gc.length)
        if totals[i] > 0
            gc[i] /= totals[i].to_f
        end
    end

    return gc
end

#Turn Q into Phred+33 ASCII-encoded quality
def convert_q_to_phred33(q_value)
    return (q_value.round() + 33).chr
end

#Turn Phred+33 ASCII encoded quality into Q
def convert_phred33_to_q(phred)
    return (phred.ord - 33)
end


secs, quals = read_fastq("test.fastq")

require 'matplotlib'
Matplotlib.use('Agg')
require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot

nts = Hash.new(0)
secs.each do |sec|
    sec.each_char do |nt|
        nts[nt] += 1
end
end

print(nts)

gc = create_GC_hist_data(secs)
plt.plot( (0...gc.length).to_a , gc  )
plt.show
plt.savefig("test2.png")


#Read a text file, then store and return all the numbers in a numeric array
def read_numbers(numeric_file)
    numbers = []
    File.open(numeric_file).each do |line|
        if !(line.nil?)
            numbers.concat(line.chomp.split(" ").map{|item| item.to_i})
        end
    end
    return numbers
end
