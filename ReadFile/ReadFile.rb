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
