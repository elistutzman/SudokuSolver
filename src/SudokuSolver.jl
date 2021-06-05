module SudokuSolver

using JuMP
using GLPK

export ClassicSudoku, solve!

struct ClassicSudoku{T, N} <: AbstractArray{T, N}
	data::Dict{NTuple{N, Int}, T}
	dims::NTuple{N, Int}
	vals::Vector{T}
end

ClassicSudoku(dims::Int) = _ClassicSudoku(Int, (dims, dims), collect(1:dims), true)

ClassicSudoku(dims::Int, vals::Vector{T}) where {T} = _ClassicSudoku(T, (dims, dims), vals, true)

ClassicSudoku(dims::NTuple{2, Int}) = _ClassicSudoku(Int, dims, collect(1:dims[1]), true)

ClassicSudoku(dims::NTuple{2, Int}, vals::Vector{T}) where {T} = _ClassicSudoku(T, dims, vals, true)

ClassicSudoku(vals::Vector{T}) where {T} = _ClassicSudoku(T, (length(vals), length(vals)), vals, true)

ClassicSudoku(::Type{T}, dims::Int, vals::Vector{T}) where {T} = _ClassicSudoku(T, (dims, dims), vals, true)

function _ClassicSudoku(::Type{T}, dims::NTuple{2, Int}, vals::Vector{T}, init::Bool) where {T}
	if init
		if dims[1] != dims[2]
			error("Dimensions are not equal")
		end

		try
			Int(sqrt(dims[1]))
		catch
			error("Dimensions are not square")
		end

		if length(vals) != dims[1]
			error("Length of vals does not match dimensions")
		end

		if length(unique(vals)) != dims[1]
			error("Number of unique vals does not match dimensions")
		end
	end

    ClassicSudoku{T, 2}(Dict{NTuple{2, Int}, T}(), dims, vals)
end

Base.size(s::ClassicSudoku) = s.dims

Base.similar(s::ClassicSudoku, ::Type{T}, dims::Dims) where {T} = _ClassicSudoku(T, dims, s.vals, false)

function Base.getindex(s::ClassicSudoku{T, N}, I::Vararg{Int, N}) where {T, N}
	if hasmethod(zero, Tuple{T})
		try
			get(s.data, I, zero(T))
		catch
			Missing
		end
	else
		get(s.data, I, Missing)
	end
end

function Base.setindex!(s::ClassicSudoku{T, N}, v, I::Vararg{Int, N}) where {T, N}
    if !(v in s.vals)
		error("Added value is not in set of values")
	end
	
	s.data[I] = v
end

function nonzero_constraint(m::Model, dims::Int)
	X = m[:X]
	for i=1:dims, j=1:dims
		@constraint(m, sum(X[i, j, k] for k=1:dims) == 1)
	end
end

function column_constraint(m::Model, dims::Int)
	X = m[:X]
	for i=1:dims, k=1:dims
		@constraint(m, sum(X[i, j, k] for j=1:dims) == 1)
	end
end

function row_constraint(m::Model, dims::Int)
	X = m[:X]
	for j=1:dims, k=1:dims
		@constraint(m, sum(X[i, j, k] for i=1:dims) == 1)
	end
end

function cage_unique_constraint(m::Model, dims::Int, cage_dims::NTuple{2, Int})
	X = m[:X]
	n_cage_row = Int(dims / cage_dims[1])
	n_cage_col = Int(dims / cage_dims[2])
	for a=1:cage_dims[1], b=1:cage_dims[2], k=1:dims
		@constraint(m, sum(X[i, j, k] for i=(n_cage_row * a - (cage_dims[1] - 1)):(n_cage_row * a) for j=(n_cage_col * b - (cage_dims[2] - 1)):(n_cage_col * b)) == 1)
	end
end

function initial_constraint(m::Model, s::ClassicSudoku, dims::Int)
	X = m[:X]
	for i=1:dims, j=1:dims
		if s[i, j] in s.vals
			@constraint(m, X[i, j, findfirst(isequal(s[i, j]), s.vals)] == 1)
		end
	end
end

function solve!(s::ClassicSudoku)
	n = size(s)[1]
	m = Int(sqrt(n))
	model = Model(GLPK.Optimizer)

	@variable(model, X[1:n, 1:n, 1:n], Bin)

	# Each cell has a non-zero value
	nonzero_constraint(model, n)

	# Each column has one of each value
	column_constraint(model, n)

	# Each row has one of each value
	row_constraint(model, n)
	
	# Each cage has one of each value
	cage_unique_constraint(model, n, (m, m))

	# Initialize starting conditions
	initial_constraint(model, s, n)

	optimize!(model)
	status = Int(termination_status(model))
	if status != 1
		error("No solution for this Classic Sudoku puzzle")
	end

	solved_X = value.(X)

	for i=1:n, j=1:n, k=1:n
		if solved_X[i, j, k] > 0
			s[i, j] = s.vals[k]
		end
	end

end

end # module
