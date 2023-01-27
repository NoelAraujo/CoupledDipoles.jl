function save_1D_simulation(f::Function, xi::Tuple{AbstractArray, AbstractArray}, file_name::String, configs)
    indices_to_execute, index_span = xi ## ex: 1:12, 1:14 -> access the first 12 elements of the vector of size 14
    x_length = length(indices_to_execute)

    x_length = length(indices_to_execute)
    all_results = [Any[] for i in indices_to_execute]

    channelSize = x_length
    channel_data = RemoteChannel(() -> Channel(channelSize));

    effetice_name = file_name*".jld2"
    @async fill_results(channel_data, all_results)
    @async save_periodically(effetice_name, all_results, configs)

    p = Progress(x_length, showspeed=true)
    progress_pmap(indices_to_execute, progress=p) do i
        x = index_span[i]
        single_data_point = f(x)

        payload = i, single_data_point
        put!(channel_data, payload)
    end

    ## saving one last time before finishing
    FileIO.save(effetice_name, Dict("payload"=>all_results , "configs"=>configs))

    println("Results (payload) saved on file '$(effetice_name)'")
    true
end

function save_2D_simulation(f::Function, xi::Tuple{AbstractArray, AbstractArray}, yi::Tuple{AbstractArray, AbstractArray}, file_name::String, configs)
    indices_to_execute_x, x_index_span = xi ## ex: 1:12, ones(14) -> access the first 12 elements of the vector of size 14
    indices_to_execute_y, y_index_span = yi

    x_length = length(x_index_span)
    y_length = length(y_index_span)
    all_results = [Any[] for i in indices_to_execute_x, j in indices_to_execute_y]

    channelSize = 20*x_length*y_length
    channel_data = RemoteChannel(() -> Channel(channelSize));

    effetice_name = file_name*".jld2"
    @async fill_results(channel_data, all_results)
    @async save_periodically(effetice_name, all_results, configs)


    all_pair_combinations = CartesianIndices((x_length,y_length)) # I need to know all combinations
    all_parts_numbers = LinearIndices(all_pair_combinations) # to figure out the (i,j) pair number in column-major structure
    local_pair_combinations = CartesianIndices((indices_to_execute_x, indices_to_execute_y)) # to make parallel computing, i avoid nested loops, and use a map

    p = Progress(length(local_pair_combinations), showspeed=true)
    progress_pmap(local_pair_combinations, progress=p) do ij_index
        i,j = ij_index[1], ij_index[2]
        x, y = x_index_span[i], y_index_span[j]
        single_data_point = f((x,y))

        payload = ij_index, single_data_point
        put!(channel_data, payload)

    end
    ## saving one last time before finishing
    FileIO.save(effetice_name, Dict("payload"=>all_results , "configs"=>configs))

    println("Results (payload) saved on file '$(effetice_name)'")
    true
end


function fill_results(getData_channel, all_results)
    while true
        try
            payload = take!(getData_channel)
            ij_index, single_data_point = payload
            push!(all_results[ij_index], single_data_point)
        catch
            println("ERROR TO SAVE ON MATRIX RESULT")
        end
    end
    nothing
end

using FileIO
function save_periodically(file_name, all_results, configs)
    while true
        sleep(15) # arbitrary time period
        try
            FileIO.save(file_name, Dict("payload"=>all_results, "configs"=>configs))
        catch
            println("ERROR TO SAVE ON DISK")
        end
    end
    nothing
end

function load_simulation(file_name::String)
    results = load(file_name*".jld2")
    payload = results["payload"]
    configs = results["configs"]
    return payload, configs
end

