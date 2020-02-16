#Given a DNA/RNA sequence and a contiguous DNA/RNA motif, returns the positions where the given motif is found 
#inside the sequence
def get_motif (dna_seq, substring)
    pos = []
    for index in (0..(dna_seq.length - substring.length))
        if dna_seq[index...(index + substring.length)] == substring
            pos.push(index + 1)
        end
    end   
    return pos
end

#Given a DNA sequence, it returns the GC percent 
def get_gc_percent(dna_seq)
    gc_counter = 0
    dna_seq.upcase.each_char do |nt|
        if nt == "C" || nt == "G"
            gc_counter += 1
        end
    end
    gc_perc = (((gc_counter.to_f / dna_seq.length)) * 100).round(6)
    return gc_perc
end 

#Given two DNA sequences of equal length, it returns the hamming distance, that is the number of differences
#between the two sequences, as a method to meassure the number of point mutations
def count_point_mutations(dna_seq1, dna_seq2)
    n_muts = 0
    for nt_pos in (0...(dna_seq1.length))
        if dna_seq1[nt_pos] != dna_seq2[nt_pos]
            n_muts += 1
        end
    end
    return n_muts
end

#Returns the number of each nucleotid of a given DNA sequence 
def count_nt(dna_seq)
    nts_array = dna_seq.upcase.split("")
    nts = {"A" => 0 , "G" => 0, "C" => 0, "T" => 0}
    
    nts_array.each do |nt|
        if nt == "A" 
            nts["A"] += 1
        end
        if nt == "C" 
            nts["C"] += 1
        end
        if nt == "T" 
            nts["T"] += 1
        end
        if nt == "G" 
            nts["G"] += 1
        end
    end
    return nts
end

#Returns the RNA transcript of a DNA template (but not his complement)
def get_rna_transcript(dna_seq)
    nts_array = dna_seq.upcase.split("")
    transcribed = []
    nts_array.each do |nt|
        if nt == "T" 
            transcribed.concat(["U"])
        else
            transcribed.concat([nt])
        end
    end
    transcribed = transcribed.join("")
    return transcribed
end

#Returns the reverse complement chain of a DNA sequence
def get_reverse_complement(dna_seq)
    nts_array = dna_seq.upcase.split("")
    reverse = nts_array.reverse
    complement = []
    reverse.each do |nt|
        if nt == "T"
            complement.concat(["A"])
        elsif nt == "A"
            complement.concat(["T"])
        elsif nt == "C"
            complement.concat(["G"])
        elsif nt == "G"
            complement.concat(["C"])
        end
    end
    complement = complement.join("")
    return complement            
end

#Given an array of DNA sequences it return the consensus sequence and the nucleotides statistics at each position,
# and the profile matrix. Sequences need to be the same length in order to work
def get_consensus_sequence(seqs_array)
    a_array = Array.new(seqs_array[0].length).fill(0)
    t_array = Array.new(seqs_array[0].length).fill(0)
    c_array = Array.new(seqs_array[0].length).fill(0)
    g_array = Array.new(seqs_array[0].length).fill(0)
    consensus = []
    for pos in (0...seqs_array[0].length)
        for seq in (0...seqs_array.length)
            if seqs_array[seq][pos] == "A" || seqs_array[seq][pos] == "a"
                a_array[pos] += 1
            elsif seqs_array[seq][pos] == "T" || seqs_array[seq][pos] == "t"
                t_array[pos] += 1
            elsif seqs_array[seq][pos] == "C" || seqs_array[seq][pos] == "c"
                c_array[pos] += 1
            elsif seqs_array[seq][pos] == "G" || seqs_array[seq][pos] == "g"
                g_array[pos] += 1
            end
        end
    end
    for pos in (0...a_array.length)
        consensus_nt = [a_array[pos], t_array[pos], g_array[pos], c_array[pos]].max()
        if consensus_nt == a_array[pos]
            consensus[pos] = "A"
        elsif consensus_nt == t_array[pos]
            consensus[pos] = "T"
        elsif consensus_nt == c_array[pos]
            consensus[pos] = "C"
        elsif consensus_nt == g_array[pos]
            consensus[pos] = "G"
        end
    end
    consensus = consensus.join("")
    a_array = a_array.join(" ")
    t_array = t_array.join(" ")
    g_array = g_array.join(" ")
    c_array = c_array.join(" ")
    profile_matrix = {"A" => a_array, "T" => t_array, "C" => c_array, "G" => g_array}
    return [consensus, profile_matrix]
end

#Given an array of sequences it returns the biggest contiguous shared motif in all the sequences
def get_shared_motif (seqs_array)
    template = seqs_array.min()
    seqs_array.delete(template)
    temp_len = template.length
    consensus = ""
    for stop in (0..temp_len).to_a.reverse
        start = temp_len - stop
        for start_iter in (0..start)
            if start_iter == (start_iter + stop)
                break
            end
            #To see how it works
            #print("range => #{start_iter}-#{start_iter+stop}, seq =>  ")
            #puts("#{template[start_iter...(start_iter + stop)]} ")
            if (seqs_array.all? { |seq| seq.include?(template[start_iter...(stop + start_iter)]) }) && consensus.length < ((template[start_iter...(stop + start_iter)]).length)
                #If you want to see other hits, place this return at the end of
                #the function and activate de line below
                #puts("hit with =>  #{template[start_iter...(start_iter + stop)]}")
                consensus = template[start_iter...(stop + start_iter)]
                return consensus
            end
        end
    end
end

#Given a DNA sequence it return position and length of every palindromic sequence inside
def get_palindromic_sequences(dna_seq)
    seq_len = dna_seq.length
    min_len = 3
    max_len = 12
    pos_and_len = []
    for palin_len in (min_len..max_len)
        rango = seq_len - palin_len
        for nt in (0...rango)
            temp = dna_seq[nt..(nt + palin_len)]
            if temp == get_reverse_complement(temp)
                puts("#{nt + 1} #{palin_len + 1}")
                puts("#{temp} - #{get_reverse_complement(temp)}")
                puts()
                pos_and_len.push([nt + 1, palin_len + 1])
            end
        end
    end
    return pos_and_len
end

#Translate a mRNA sequence into a protein sequence. If stop_required is true, then it returns nil if 
# no stop codon was found
def get_prot_seq(mrna, stop_required = false)
    prot = ""
    codons = {"GCU" => "A","GCC" => "A","GCA" => "A","GCG" => "A","GUU" => "V","GUC" => "V",
             "GUA" => "V","GUG" => "V","GAU" => "D","GAC" => "D","GAA" => "E","GAG" => "E",
             "GGU" => "G","GGC" => "G","GGA" => "G","GGG" => "G","AUU" => "I","AUC" => "I",
             "AUA" => "I","AUG" => "M","ACU" => "T","ACC" => "T","ACA" => "T","ACG" => "T",
             "AAU" => "N","AAC" => "N","AAA" => "K","AAG" => "K","AGU" => "S","AGC" => "S",
             "AGA" => "R","AGG" => "R","CUU" => "L","CUC" => "L","CUA" => "L","CUG" => "L",
             "CCU" => "P","CCC" => "P","CCA" => "P","CCG" => "P","CAU" => "H","CAC" => "H",
             "CAA" => "Q","CAG" => "Q","CGU" => "R","CGC" => "R","CGA" => "R","CGG" => "R",
             "UUU" => "F","UUC" => "F","UUA" => "L","UUG" => "L","UCU" => "S","UCC" => "S",
             "UCA" => "S","UCG" => "S","UAU" => "Y","UAC" => "Y","UAA" => "STOP",
             "UAG" => "STOP","UGU" => "C","UGC" => "C","UGA" => "STOP","UGG" => "W"}
    for i in ((0...mrna.length).step(3))
        if codons[mrna[i...(i+3)]] == "STOP"
            return prot
        elsif
            prot += codons[mrna[i...(i+3)]]
        end
    end
    if stop_required
        return nil
    end
    return prot
end

#Given a DNA sequence and a array with of introns, it returns the ORF of the DNA sequence (only with the exons)
def eliminate_introns(dna_seq, introns)
    template = dna_seq
    introns.each do |intron|
        template = template.split(intron)
        template = template.join("")
    end
    return template
end

#Given a DNA/RNA/protein sequence and a DNA/RNA/protein motif, it return the coordinates where the 
#the spliced motif can be found inside the DNA sequence, or returns Spliced motif
#not found if no full asociation could be done
def get_spliced_motif (dna_seq, spliced_motif)
    dna_seq = dna_seq.split("")
    spliced_motif = spliced_motif.split("")
    act_pos = 0
    motif_coordinates = []
    spliced_motif.each do |motif_nt|
        for seq_nt in (0...dna_seq.length)
            if (motif_nt == dna_seq[seq_nt]) && (seq_nt > act_pos)
                motif_coordinates.push(seq_nt + 1)
                act_pos = seq_nt
                break
            end
        end
    end
    
    if motif_coordinates.length == spliced_motif.length
        return motif_coordinates
    else
        return "Spliced motif not found"
    end
end

#Given a protein sequence it return every location where the N-Glycosylation motif is found
def find_nglycosylation_motif(prot_seq)
    pos = []
    for ac_pos in (0...(prot_seq.length - 4))
        if prot_seq[ac_pos] == "N"
            if prot_seq[ac_pos + 1] != "P"
                if prot_seq[ac_pos + 2] == "S" || prot_seq[ac_pos + 2] == "T"
                    if prot_seq[ac_pos + 3] != "P"
                        pos = pos.concat([(ac_pos + 1)])
                    end
                end
            end
        end
    end
    if !(pos.nil?)
        return pos
    end
end

#Given a protein IDs array it searches the aminoacidic sequence at uniprot web and give a hash with all the 
#protein ids who has N-Glycosylation motifs in his sequence and the locations of the motifs. This function needs
#to use read_online_fasta function from ReadFile.rb in order to work
def get_prots_with_nglyc_motif(prots_id_array)
    prots_hash = {}
    motifs = {}
    prots_id_array.each do |id|
        prots_hash[id] = (read_online_fasta("http://www.uniprot.org/uniprot/#{id}.fasta")).values[0]
    end
    prots_hash.each do |id, seq|
        pos = find_nglycosylation_motif(seq)
        if !(pos.empty?)
            motifs[id] = pos
        end
    end
    return motifs
end

#Given a protein sequence it returns the monoisotopic mass in dalton
def get_prot_mass(prot_seq)
    weights = {"A" =>71.03711,"C"  => 103.00919,"D"  => 115.02694,"E"  => 129.04259,"F"  => 147.06841,
        "G"  => 57.02146,"H"  => 137.05891,"I"  => 113.08406,"K"  => 128.09496,"L"  => 113.08406,
        "M"  => 131.04049,"N"  => 114.04293,"P"  => 97.05276,"Q"  => 128.05858,"R"  => 156.10111,
        "S"  => 87.03203,"T"  => 101.04768,"V"  => 99.06841,"W"  => 186.07931,"Y" =>  163.06333,
        "H2O" => 18.01056 }
    total_w = 0
    prot_seq.each_char do |ac|
        total_w += weights[ac]
    end
    return (total_w.round(3))
end

#Given an array of DNA IDs and an array of sequences, it returns the directed graph associations
#(as an adjacency list) obtained by verifying that the tail of sequence X overlap with the head of 
#sequence Y with length given by window's variable
def get_direct_graph_associations(dna_ids_array, dna_seqs_array, window )
    assoc = []
    for tail in (0...dna_seqs_array.length)
        for head in (0...dna_seqs_array.length)
            if tail != head
                if dna_seqs_array[tail][((dna_seqs_array[tail].length - window)...(dna_seqs_array[tail].length))] == dna_seqs_array[head][(0...window)]
                    assoc.push([dna_ids_array[tail],dna_ids_array[head]])
                end
            end
        end
    end
    return assoc
end

#Given a certain DNA sequence, return all the theorical ORF (from ATG to STOP codon)
#found in that sequence and his reverse complement, in the 3 different offsets
def search_orf(dna_seq)
    rna_seq = get_rna_transcript(dna_seq)
    rna_revcom = get_rna_transcript(get_reverse_complement(dna_seq))
    orf1 = []
    orf2 = []
    orf3 = []
    orf4 = []
    orf5 = []
    orf6 = []
    
    for i in (0...(dna_seq.length - 2))
        if ((i+3)%3) == 0     
            if (rna_seq[(i...i+3)]).length == 3
                orf1.push(rna_seq[(i...i+3)])
            end
            if (rna_seq[(i+1...i+4)]).length == 3
                orf2.push(rna_seq[(i+1...i+4)])
            end
            if (rna_seq[(i+2...i+5)]).length == 3
                orf3.push(rna_seq[(i+2...i+5)])
            end
        end
    end
    for i in (0...(dna_seq.length - 2))
        if ((i+3)%3) == 0
            if (rna_revcom[(i...i+3)]).length == 3
                orf4.push(rna_revcom[(i...i+3)])
            end
            if (rna_revcom[(i+1...i+4)]).length == 3
                orf5.push(rna_revcom[(i+1...i+4)])
            end
            if (rna_revcom[(i+2...i+5)]).length == 3
                orf6.push(rna_revcom[(i+2...i+5)])
            end
        end
    end

    orf1 = orf1.join("")
    orf2 = orf2.join("")
    orf3 = orf3.join("")
    orf4 = orf4.join("")
    orf5 = orf5.join("")
    orf6 = orf6.join("")
    orfs = [orf1, orf2, orf3, orf4, orf5, orf6]
    prots = []
    orfs.each do |orf|
        for j in (0...(orf.length))
            if ((j+3)%3) == 0
                if orf[(j...j+3)] == "AUG" 
                    prot = get_prot_seq(orf[(j...(orf.length))], stop_required = true)
                    if !(prot.nil?) && !(prot.empty?)
                        prots.push(prot)
                    end
                end
            end
        end
    end
    return prots.uniq
end

#Given an array of sequences(reads) of the same length and a threshold of nucleotides length for the overlap
#(a percentage from 1 to 100), i.e. the amount of nucleotides length to validate the overlaps, returns 
#a superstring containing all the overlapping sequences, that is, the genome assembly of the reads, with a 
#maximum error(percentage from 1 to 100)to validate each pair of overlaps
def assamble_reads(dna_seqs_array, min_overlap = 50, max_error = 40)
    similarity_index = (100-max_error).to_f / 100
    total_read_length = dna_seqs_array[0].length
    matches = {}
    superstring = ""
    for head in (0...dna_seqs_array.length)
        skip = false
        for tail in (0...dna_seqs_array.length)
            if skip
                break
            end
            if head != tail && !(matches.values.include?(tail))
                if check_overlap(dna_seqs_array[head], dna_seqs_array[tail], min_overlap, similarity_index)
                    matches[head] = tail
                    skip = true
                end
            end
        end
    end

    start = nil
    stop = nil
    for pos in (0...dna_seqs_array.length)
        if !(start.nil?()) && !(stop.nil?())
            break
        end
        if matches.keys.include?(pos) && !(matches.values.include?(pos))
            start = pos
        elsif !(matches.keys.include?(pos)) && (matches.values.include?(pos))
            stop = pos
        end
    end

    while start != stop
        indexx = check_overlap(dna_seqs_array[start], dna_seqs_array[matches[start]], min_overlap, similarity_index, return_index = true)
        superstring += dna_seqs_array[start][(0...indexx)]
        start = matches[start]
    end
    superstring += dna_seqs_array[stop]
    return superstring
end

#Support function for assamble_reads
def check_overlap(head, tail, min_overlap, sim_index, return_index = false)
    min_overlap_length = ((([head, tail].min()).length) * ((min_overlap.to_f)/100)).to_i
    for window in (0...min_overlap_length)
        if get_similarity(head[(window...head.length)], tail) >= sim_index
            if return_index
                return window
            end
            return true
        end
    end
    return false
end

#Given two DNA sequences it returns a percentage of similarity between the two sequences. It can be used alone
#or as a support function for assamble_reads
def get_similarity(seq1, seq2)
    total = 0
    for nt in (0...(([seq1,seq2].min()).length))
        if seq1[nt] == seq2[nt]
            total += 1
        end
    end
    return (total.to_f / (([seq1,seq2].min()).length))
end

#Given two DNA sequences, it calculates the transition/transversion ratio, that is, the number of mutations of the
# type purine to purine (ex: G -> A) or pyrimidine to pyrimidine (T -> C) versus the number of transversion,
# that is, from purine to pyrimidine or viceversa (ex: G -> C) 
def calculate_transition_transversion_ratio(seq1, seq2)
    transi = 0
    transve = 0
    for nt in (0...(([seq1,seq2].max()).length))
        if seq1[nt] != seq2[nt]
            if (seq1[nt] == "G" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "G")
                transi += 1
            elsif (seq1[nt] == "G" && seq2[nt] == "T") || (seq1[nt] == "T" && seq2[nt] == "G")
                transve += 1
            elsif (seq1[nt] == "G" && seq2[nt] == "C") || (seq1[nt] == "C" && seq2[nt] == "G")
                transve += 1
            elsif (seq1[nt] == "T" && seq2[nt] == "C") || (seq1[nt] == "C" && seq2[nt] == "T")
                transi += 1
            elsif (seq1[nt] == "T" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "T")
                transve += 1
            elsif (seq1[nt] == "C" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "C")
                transve += 1
            end
        end
    end
    return (((transi.to_f) / transve).round(11))
end

#Given a protein sequence it return the number of theorical mRNAs that may
#produce the same protein
def calculate_number_of_possible_RNAs_from_prot(protseq)
    result = 1  #Variable for returning product of each aac number of codons
    stop = 3 #Variable for remembering that we have to take account of the stop codons too
    codons ={"GCU"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A","GUU"=>"V","GUC"=>"V",
             "GUA"=>"V","GUG"=>"V","GAU"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
             "GGU"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G","AUU"=>"I","AUC"=>"I",
             "AUA"=>"I","AUG"=>"M","ACU"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
             "AAU"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K","AGU"=>"S","AGC"=>"S",
             "AGA"=>"R","AGG"=>"R","CUU"=>"L","CUC"=>"L","CUA"=>"L","CUG"=>"L",
             "CCU"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P","CAU"=>"H","CAC"=>"H",
             "CAA"=>"Q","CAG"=>"Q","CGU"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
             "UUU"=>"F","UUC"=>"F","UUA"=>"L","UUG"=>"L","UCU"=>"S","UCC"=>"S",
             "UCA"=>"S","UCG"=>"S","UAU"=>"Y","UAC"=>"Y","UAA"=>"STOP",
             "UAG"=>"STOP","UGU"=>"C","UGC"=>"C","UGA"=>"STOP","UGG"=>"W"}
    
    protseq.each_char do |aac|
        count = 0 #it saves the number of codons that produce a certain aminoacid
        codons.each do |codon, translated_codon|
            if translated_codon == aac
                count += 1
            end
        end
        result *= count
    end
    return (result * stop) % 1000000
end