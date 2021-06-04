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
	try
		get(s.data, I, zero(T))
	catch
		Missing
	end
end

function Base.setindex!(s::ClassicSudoku{T, N}, v, I::Vararg{Int, N}) where {T, N}
    if !(v in s.vals)
		error("Added value is not in set of values")
	end
	
	s.data[I] = v
end

function solve!(s::ClassicSudoku)
	n = size(s)[1]
	m = Int(sqrt(n))
	model = Model(GLPK.Optimizer)

	@variable(model, X[1:n, 1:n, 1:n], Bin)

	# Each cell has a non-zero value
	for i=1:n, j=1:n
		@constraint(model, sum(X[i, j, k] for k=1:n) == 1)
	end

	# Each column has one of each value
	for i=1:n, k=1:n
		@constraint(model, sum(X[i, j, k] for j=1:n) == 1)
	end

	# Each row has one of each value
	for j=1:n, k=1:n
		@constraint(model, sum(X[i, j, k] for i=1:n) == 1)
	end
	
	# Each cage has one of each value
	for a=1:m, b=1:m, k=1:n
		@constraint(model, sum(X[i, j, k] for i=3a-2:3a for j=3b-2:3b) == 1)
	end

	# Initialize starting conditions
	for i=1:n, j=1:n
		if s[i, j] ∈ s.vals
			@constraint(model, X[i, j, findfirst(isequal(s[i, j]), s.vals)] == 1)
		end
	end

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
