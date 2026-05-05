using Nemo

function tuple_to_sbox_row(tuple)
    A, B = tuple
    return vcat(matrix_to_sbox(A), matrix_to_sbox(B))
end

function tuples_to_sbox_rows(tuples)
    return [tuple_to_sbox_row(tuple) for tuple in tuples]
end

function write_c_array(file, rows)
    if isempty(rows)
        write(file, "int AllTuples[0][0] = {};\n")
        return
    end

    row_size = length(rows[1])

    write(file, "#define N_TUPLES $(length(rows))\n\n")
    write(file, "int AllTuples[N_TUPLES][$row_size] = {\n")

    for i in 1:length(rows)
        row = rows[i]
        row_text = join(row, ", ")

        if i < length(rows)
            write(file, "    {$row_text},\n")
        else
            write(file, "    {$row_text}\n")
        end
    end

    write(file, "};\n")
end

function generate_tuples_file(tuples, filename::String = "tuples/AllTuples.h")
    rows = tuples_to_sbox_rows(tuples)

    open(filename, "w") do file
        write_c_array(file, rows)
    end
end