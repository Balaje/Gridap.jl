module Mappings

#@fverdugo
# For me it is very confusing to import so many parts of the old Gridap.
# It is difficult to understand what is working thanks to the old implementation
# and what thanks to the new one. In particular, importing Gridap.Arrays which has a lot of implementation
# that collides with what is done in this new module.
using Gridap.Helpers
using Gridap.Inference
using Gridap.Arrays
using Gridap.Algebra: mul!
using FillArrays
using Test

import Gridap.Inference: return_type

# Mapping interface

export Mapping
export return_cache
export evaluate!
export evaluate
export return_type
export test_mapping

import Gridap.Arrays: testitem
import Gridap.Arrays: getindex!
import Gridap.Arrays: uses_hash
import Gridap.Arrays: IndexStyle

# MappedArray

export MappedArray

export apply
export test_mapped_array
import Gridap.Arrays: array_cache
# import Gridap.Arrays: apply

export BroadcastMapping
export OperationMapping
export Operation
export operation

# Field

using Gridap.TensorValues
using ForwardDiff
using LinearAlgebra: mul!, Transpose

import LinearAlgebra: det, inv, transpose
import Base: +, -, *, /
import LinearAlgebra: ⋅

export Field
export GenericField
export FieldGradient
export FieldHessian
export BroadcastField
export ZeroField
export MockField
export Point

export evaluate_gradient!
export return_gradient_type
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache
export return_hessian_type

export gradient
export ∇
export hessian

export test_field
export test_field_array
export test_operation_field_array
export test_broadcast_field_array

export mock_field

export TransposeFieldVector
export TransposeFieldIndices
export BroadcastOpFieldArray
export DotOpFieldVectors
export LinearCombinationField
export CompositionFieldArrayField
export FieldGradientArray
export FieldHessianArray

export linear_combination
export integrate

export MatMul
export LinCombVal
export Integrate

include("MappingInterfaces.jl")

include("MappedArrays.jl")

include("FieldsInterfaces.jl")

include("FieldArrays.jl")

include("ApplyOptimizations.jl")

include("MockFields.jl")

include("AutoDiff.jl")

include("ArraysMappings.jl")

end
