#Binary search function. Given a sorted array and a item, it returns the position where the item is in the array
#(if it exist), otherwise it returns nil(null). It runs at O(log n) speed.
def use_binary_search(list, item)
    low = 0
    high = list.length - 1
    while low <= high
        mid = (low + high)
        guess = list[mid]
        if guess == item
            return mid
        end
        if guess > item
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return nil
end

#Algoritmic function to find the smallest item in a given array. It runs at O(n) speed as it assumes 
# the array is unsorted. It works as a helper function for use_selection_sort
def find_smallest(arr)
    smallest = arr[0] #Stores the smallest value
    smallest_index = 0 #Stores the index of the smallest value
    for i in (1...arr.length)
        if arr[i] < smallest
            smallest = arr[i]
            smallest_index = i
        end
    end
    return smallest_index
end

#Algoritmic function to sort an unsorted array. It runs at O(n x n) speed.
def use_selection_sort(arr)
    newArr = []
    for i in (0...arr.length)
        smallest = find_smallest(arr)
        puts(arr[smallest])
        newArr.push(arr.delete_at(smallest))
    end
    return newArr
end

#Algoritmic function to sort an unsorted array. It runs at O(n x n) speed in the worst case, O(n x log n) on 
#the average case.
def use_quick_sort(arr)
    if arr.length < 2
        return arr #Base case: arrays with 0 or 1 element are already “sorted.”
    else
    pivot = arr[0] #Recursive case
    less =  arr[(1...arr.length)].select{|x| x <= pivot} #Sub-array of all the elements less than the pivot
    greater = arr[(1...arr.length)].select{|x| x > pivot} #Sub-array of all the elements greater than the pivot
    return use_quick_sort(less) + [pivot] + use_quick_sort(greater)
    end
end

#Breadth first search. An Algoritmic function to be used with unweigthed graphs. It searches from "item" until it
# finds another item that fulfills some lambda logic_function (be a certain item, be greater than, etc) and 
#returns some of the following options:
#1)true if searched item is found, (then use returned = "found");
#2)the number of steps required, (then use returned = "steps")
def use_breadth_first(item, graph, logic_function = ->(x){graph[x].empty?} , returned = "steps" )
    search_queue = []
    steps = {}
    search_queue = search_queue.concat(graph[item])
    searched = []
    #Setting up initial steps 
    if !search_queue.empty?
        search_queue.each do |term|
            steps[term] = 1
        end
    end
    #Here goes the graph algorithm
    while !search_queue.empty?
        person = search_queue.shift()
        if !( searched.include?(person) )

            if logic_function.call(person)
                if returned == "steps"
                    return steps[person]
                end
                if returned == "found"
                    return true
                end
            else
                if !(graph[person].nil?) 
                    graph[person].each do |related|
                        steps[related] = steps[person] + 1  #Setting up the steps of parents of the current element in the queue
                    end
                    search_queue = search_queue.concat(graph[person])
                end
            end

        end
    end
    return false
end