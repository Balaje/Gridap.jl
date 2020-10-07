module Mappings

using Gridap.Arrays: CachedArray
using Gridap.Arrays: setsize!

# using Gridap.Arrays: getitems!
# using Gridap.Arrays: get_array

using Gridap.TensorValues: get_array

import Gridap.Arrays: array_cache
import Gridap.Arrays: testitem
import Gridap.Arrays: getindex!
# import Gridap.Arrays: uses_hash
import Gridap.Arrays: IndexStyle

using Gridap.Helpers: GridapType, @abstractmethod, @notimplemented
using Gridap.Helpers: @notimplementedif, @unreachable

using Gridap.Inference: testvalue
import Gridap.Inference: return_type

using Gridap.Algebra: mul!

using FillArrays
using Test

# Mapping interface

# export Mapping

# export return_cache
# export evaluate!
# export evaluate
# export return_type
# export test_mapping

# LazyArray

# export LazyArray

# import Gridap.Arrays: lazy_map
# import Gridap.Arrays: array_cache

# export lazy_map
# export array_cache
# export test_lazy_array

export Broadcasting
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

# Mapping

using Gridap.Arrays: Mapping
using Gridap.Arrays: Operation
using Gridap.Arrays: Broadcasting
using Gridap.Arrays: evaluate
using Gridap.Arrays: test_mapping
using Gridap.Arrays: test_lazy_array

import Gridap.Arrays: evaluate!
import Gridap.Arrays: return_cache
import Gridap.Arrays: return_type

export evaluate!
export return_type
export return_cache

# Lazy Arrays

using Gridap.Arrays: LazyArray
using Gridap.Arrays: getindex!

import Gridap.Arrays: array_cache
import Gridap.Arrays: lazy_map

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

include("FieldsInterfaces.jl")

include("MockFields.jl")

include("FieldArrays.jl")

include("ApplyOptimizations.jl")

include("AutoDiff.jl")

include("ArraysMappings.jl")

end
