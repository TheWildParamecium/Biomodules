#Fibonacci sequence. Months represent the number of generations after we want to retrieve a results,
#and offspring the offspring produced by each pair at each generation
def fibonacci_rabbits(months, offspring, sequence = false)
    pairs = [1, 1]
    x = 2
    while x < months
        pairs.push( (pairs[x-1]) + (pairs[(x-2)]*offspring)  )
        x += 1
    end
    if sequence
        return pairs
    end
    return pairs[months - 1]
end

#Mortal rabbits Fibonacci sequence. Months represent the number of generations after we want to retrieve a results,
#and lived months the number of months that each pair of rabbits produced at a given generation live
def mortal_fibonacci_rabbits(months, lived_months)
    pairs = [1, 1]
    x = 2
    while x < months
        total = 0
        if lived_months > x
            pairs.push( (pairs[x-1]) + (pairs[(x-2)])  )
        else
            for n in (2..lived_months)
                total += pairs[x - n]
            end
            pairs.push(total)
        end
        x += 1
    end
    return pairs[months - 1]
end

#Mendel firt Law. Given k,m,n as the number of homozygous dominant, heterozygous and homozygous recessive indiviuals,
#it gives the probability that two random organisms selected will produce a dominant allele (either homozygous or 
# heterozygous)
def get_dominant_allele_probability(k, m, n)
    poblation = k + m + n
    probs = {}
    probs["mn"] = ( (m.to_f/poblation)*(n.to_f/(poblation - 1)) + (n.to_f/poblation)*(m.to_f/(poblation -1)) ) * 0.5
    probs["kn"] = ( (k.to_f/poblation)*(n.to_f/(poblation - 1)) + (n.to_f/poblation)*(k.to_f/(poblation - 1)) ) * 1
    probs["km"] = ( (k.to_f/poblation)*(m.to_f/(poblation - 1)) + (m.to_f/poblation)*(k.to_f/(poblation - 1)) ) * 1
    probs["kk"] = ( (k.to_f/poblation)*((k.to_f - 1)/(poblation - 1))) * 1
    probs["mm"] = ( (m.to_f/poblation)*((m.to_f - 1)/(poblation - 1))) * 0.75
    return (probs.values.sum()).round(5)
end

#Given the number of couples with a certain genetic background at a given allele and the number of childs for
#each couple, it returns the expected number organisms with the dominant phenotype 
def calculate_expected_dominant_offspring(xAA_AA, xAA_Aa, xAA_aa, xAa_Aa, xAa_aa, xaa_aa, offspring_each = 2)
    xAA_AA_prob = offspring_each.to_f * xAA_AA * 1 
    xAA_Aa_prob = offspring_each.to_f * xAA_Aa * 1 
    xAA_aa_prob = offspring_each.to_f * xAA_aa * 1 
    xAa_Aa_prob = offspring_each.to_f * xAa_Aa * 0.75 
    xAa_aa_prob = offspring_each.to_f * xAa_aa * 0.5 
    xaa_aa_prob = offspring_each.to_f * xaa_aa * 0 
    return xAA_AA_prob + xAA_Aa_prob + xAA_aa_prob + xAa_Aa_prob + xAA_aa_prob + xAa_aa_prob + xaa_aa_prob
end


#Given a n value it returns all the possible permutations of length n with that n numbers
def get_permutations(n)
    return (1..n).to_a.permutation.to_a
end

#Given numbers N and K, counts the total number of partial permutations, also known as variations
#of k objects that can be formed from a collection of n objects
def get_partial_permutations(n,k)
    result = 1
    while (k > 0)
        result *= n
        n -= 1
        k -= 1
    end
    return (result % 1000000)
end

#Given a number N it returns every permutations of length N of the N numbers between -N and N, without repeating
#the same number inside every permutation
def get_signed_permutations(n)
    numbers = (-n..n).to_a
    numbers.delete(0)
    perms = numbers.permutation(n).to_a
    perms = perms.select{|item| item.map{|subitem|subitem.abs()}.uniq.length == item.length}
    return perms.uniq
end

#Returns the factorial of a given number
def get_factorial(n)
    total = 1
    while n > 1
        total *= n
        n -= 1
    end
    return total
end
