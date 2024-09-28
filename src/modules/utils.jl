module Utils

using DelimitedFiles
using Dates
using JSON


"""
    saveData(data [; path::String = "./data/", name::Union{Missing, String} = missing])

Save `data` to JSON file `name`.json in `path` directory. If `name` not provided, current time in format "Y-mm-dd_HH:MM:SS.s" is used as file name. Requires `typeof(data) <: AbstractDict`.
"""
function saveData(data; path::String = "./data/", name::Union{Missing, String} = missing)
    if name === missing
        name = Dates.format(now(), "Y-mm-dd_HH:MM:SS.s")
    end

    file = open(string(path, name, ".json"), "w")
    JSON.print(file, data, 2)
    close(file)

    return nothing
end


end
