module APNLib

using Nemo

include("core/finite_field.jl")
include("core/matrix.jl")
include("core/nemo_compat.jl")
include("core/rcf.jl")
include("apn/tuples.jl")
include("apn/export.jl")
include("apn/fileUtils.jl")
include("apn/search.jl")

export gf2

export int_to_bits,
       int_to_column_vector,
       column_vector_to_int,
       each_column_vector

export check_square,
       check_same_size,
       check_same_field,
       check_compatible_pair

export matrix_multiplicative_order,
       matrix_is_similar

export check_order_space,
       filtro_proposicao_4,
       filtro_proposicao_5,
       addDDTInformation,
       removeDDTInformation,
       isComplete,
       nextFreePosition,
       nextVal,
       APNSearchContext,
       APNSearch,
       APNsearch,
       matrix_to_sbox,
       is_permutation_tuple,
       is_power_similar,
       is_extended_power_similar,
       gen_permutation_tuples

export blocks_for_rcf,
       get_rcfs

export generate_tuples_file,
       generate_tuple_constants_file,
       generate_tuple_constants_files,
       generate_tuple_matrix_constants_file,
       generate_tuple_matrix_constants_files,
       extrair_matrizes,
       extract_matrices_from_tuple_lut,
       load_precomputed_tuple_constants,
       load_precomputed_tuple_matrix_constants,
       precomputed_tuple_matrices,
       precomputed_tuple_row,
       precomputed_tuple_sboxes,
       tuple_to_sbox_row,
       tuples_to_sbox_rows

end
