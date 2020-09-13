#K-Nearest Neighbor Algorithm. It classifies new data according to the k_number of nearest neighbors points.
#Parameters are the following:
#-training is an array of arrays with each subarray being a data (with all of his numeric variables) and the last
#item in the array being the response/category;
#-prediction is the data (with all of his numeric variables) that we want to predict the response/category; 
#-k_number is the number of nearest point that are took into account in order to stablish the category of the prediction.    
def use_KNN(training, prediction, k_number = 5)
    y = []
    distances = []
    training.each do |vector|
        y.push(vector.pop)
        distances.push(calculate_distance(prediction, vector))
    end
    k_first = []

    k_number.times do
        index = distances.index(distances.min)
        k_first.push(y[index])

        distances.delete_at(index || distances.length)
        y.delete_at(index || distances.length)
    end

    result = k_first.max_by { |i| k_first.count(i) }
    return result
end

#Support function for use_KNN. It measures the distance between two points in a p1.length-dimensional space.
def calculate_distance(p1, p2)
    sum = 0
    for i in (0...p1.length)
        sum += (p1[i] - p2[i])**2
    end
    sum = sum.to_f
    return Math.sqrt(sum)
end

#Used to read Iris dataset. It can read every CSV file whose last column contains the response variable.
def read_csv_numeric_dataset(file)
    dataset = []
    File.readlines(file).each do |line|
        if !(line.empty?)        
            vector = []
            fields = line.split(",")
            response = fields.pop
            fields.each do |value|
                vector.push(value.to_f)
            end
            vector.push(response.strip)
            dataset.push(vector)
        end
    end
    return dataset
end


=begin
#If you want to test it, remove the comment block (=begin and =end)
dataset = read_csv_numeric_dataset("iris.data")
result = use_KNN(dataset, [6.5, 3.0, 5.2, 2.0])
puts(result)
=end