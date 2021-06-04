module SudokuSolver

using JuMP
using GLPK

export ClassicSudoku, solve!

struct ClassicSudoku{T, N} <: AbstractArray{T, N}
	data::Dict{NTuple{N, Int}, T}
	dims::NTuple{N, Int}
end

ClassicSudoku(dims::Int) = ClassicSudoku((dims, dims))

ClassicSudoku(dims::NTuple{2, Int}) = ClassicSudoku{Int, 2}(Dict{NTuple{2, Int}, Int}(), dims)

Base.size(s::ClassicSudoku) = s.dims

Base.similar(s::ClassicSudoku, ::Type{T}, dims::Dims) where {T} = ClassicSudoku(dims)

Base.getindex(s::ClassicSudoku{T, N}, I::Vararg{Int, N}) where {T, N} = get(s.data, I, zero(T))

Base.setindex!(s::ClassicSudoku{T, N}, v, I::Vararg{Int, N}) where {T, N} = (s.data[I] = v)

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
		if s[i, j] != 0
			@constraint(model, X[i, j, s[i, j]] == 1)
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
			s[i, j] = k
		end
	end

end

end # module
