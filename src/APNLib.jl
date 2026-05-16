module APNLib

using Nemo

function include_if_exists(relative_path::String)::Bool
    path = joinpath(@__DIR__, relative_path)
    isfile(path) || return false
    include(relative_path)
    return true
end

include("core/finite_field.jl")
include("core/matrix.jl")
include("core/nemo_compat.jl")
include("core/rcf.jl")
include("apn/tuples.jl")
include("apn/fileUtils.jl")
include("apn/apn.jl")
include("apn/search.jl")
include("core/trace.jl")
include("core/multiplicities.jl")
include("core/representations.jl")
include("apn/families.jl")
include_if_exists("apn/polynomial_json.jl")
include_if_exists("apn/database.jl")
include_if_exists("apn/database_import.jl")

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
       proposition4_filter,
       proposition5_filter,
       standard_visit_order,
       offset_visit_order,
       c_reference_visit_order,
       addDDTInformation,
       removeDDTInformation,
       isComplete,
       nextFreePosition,
       nextVal,
       APNSearchContext,
       APNSearch,
       APNSearchClasses,
       int_to_field_element,
       interpolate_sbox_polynomial,
       format_sbox_polynomial,
       absolute_trace_to_field,
       trace_to_field,
       absolute_trace_bit,
       trace_sign,
       field_elements,
       walsh_coefficient,
       walsh_spectrum,
       walsh_coefficient_table,
       multiplicities_sigma,
       ANFCoordinate,
       ANFVector,
       field_element_to_int,
       univariate_to_lut,
       lut_to_univariate,
       lut_to_anf,
       anf_to_lut,
       lut_to_graph,
       graph_to_lut,
       univariate_to_graph,
       anf_to_graph,
       graph_to_univariate,
       graph_to_anf,
       univariate_to_anf,
       anf_to_univariate,
       is_apn,
       matrix_to_sbox,
       permutation_cycle_structure,
       matrix_cycle_structure,
       same_cycle_structure,
       is_permutation_tuple,
       is_power_similar,
       is_extended_power_similar,
       gen_permutation_tuples,
       APNCoefficient,
       OneCoefficient,
       PowerCoefficient,
       APNTerm,
       APNTraceTerm,
       APNFunction,
       APNFamilyMatch,
       x,
       pterm,
       Tr,
       parse_apn_function,
       classify_family,
       belongs_to_family,
       parse_polynomial_expression_terms,
       polynomial_expression_json,
       build_polynomial_from_json,
       open_apn_database,
       init_apn_database!,
       with_apn_database,
       apn_function_id,
       apn_public_function_id,
       upsert_table_values!,
       update_table_values!,
       json_table,
       json_row_dict,
       split_int_values,
       require_json_columns,
       insert_apn_function_json!,
       insert_apn_functions!,
       insert_apn_function_table_json!,
       family_c3,
       family_c4,
       family_c5,
       family_c6

export blocks_for_rcf,
       get_rcfs

export generate_tuples_file,
       generate_tuple_constants_file,
       generate_tuple_constants_files,
       generate_tuple_matrix_constants_file,
       generate_tuple_matrix_constants_files,
       extract_matrices_from_tuple_lut,
       load_precomputed_tuple_constants,
       load_precomputed_tuple_matrix_constants,
       all_precomputed_tuple_class_indices,
       parse_class_indices_argument,
       normalize_precomputed_tuple_classes,
       precomputed_tuple_matrices,
       precomputed_tuple_row,
       precomputed_tuple_sboxes,
       save_search_result_constant,
       search_result_constant_name,
       tuple_to_sbox_row,
       tuples_to_sbox_rows

end
