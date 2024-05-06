# Define your array of strings
A = ["Hello_World", "Hello_Carl", "Carl_Dieter", "HelloWorld", "Hello", "World"]
check_str = ["Hello", "World"]
# Function that checks a string to only contain defined substrings
function check_str_only_def_substr(s::String, strings_to_check::Vector{String})
    parts = split(s, "_")

    check_vec = falses(length(parts))
    i = 1

    for part in parts
        for stri in strings_to_check
            if part == stri
                check_vec[i] = true
                break
            end
        end
        i += 1
    end

    return all(check_vec)
end

# Function that iterates over an array of strings and checks if it only contains defined substrings (or "_")
function check_array_string_only_substrings(s_vec::Vector{String}, strings_to_check::Vector{String})
    check_res = falses(length(s_vec))
    j = 1
    for stri in s_vec
        res = check_str_only_def_substr(stri, [strings_to_check..., "_"])
        check_res[j] = res
        j += 1
    end
    return hcat(s_vec,check_res)
end