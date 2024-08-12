# import Pkg
using JuMP, HiGHS
# Pkg.add("DataFrames")
using DataFrames

@enum(Direction, NORTH, EAST, SOUTH, WEST)

using Random
Random.seed!(0)

function main(Φ, d, w1, w2, t0, v_avg, v_max, avg_a, max_a, grid_dim, L)
    n = sum(length, d)
    model = Model(HiGHS.Optimizer)

    # Use the intersection topology together with the initial distribution of vehicles on the grid to decide how many t_access variables are needed
    # The counting convention is as follows:
    #               I_1                   |               I_2                   | ... | I_{grid_dim[2]}
    #        I_{grid_dim[2] + 1}          |        I_{grid_dim[2] + 2}          | ... | I_{2*grid_dim[2]}
    #               ...                   |               ...                   | ... |        ...
    # I_{(grid_dim[1]-1)*grid_dim[2] + 1} | I_{(grid_dim[1]-1)*grid_dim[2] + 2} | ... | I_{grid_dim[1]*grid_dim[2]}
    custom_mod(x, n) = (r = x % n; r == 0 ? n : r) # Help function
    indices = []
    for i in 1:size(Φ)[1]
        for j in 1:size(Φ[i])[1]
            if i==1
                vehicle_number = j
            elseif i>1
                vehicle_number = sum(length(Φ[k]) for k in 1:i-1) + j
            end
            if Φ[i][j] == NORTH
                number_of_intersections_north = ceil(i/grid_dim[2])  #i%grid_dim[1] (this commented code was wrong)
                indices = vcat(indices, [(vehicle_number, j) for j in i:-1*grid_dim[2]:custom_mod(i, grid_dim[2])])
            elseif Φ[i][j] == EAST
                number_of_intersections_east =  (grid_dim[2] - i%grid_dim[2])%grid_dim[2] + 1
                indices = vcat(indices, [(vehicle_number, j) for j in i:(i+number_of_intersections_east-1)])
            elseif Φ[i][j] == SOUTH
                number_of_intersections_south = grid_dim[1] - ceil(i/grid_dim[2]) + 1 # grid_dim[1] - i%grid_dim[1] (this commented code was wrong)
                indices = vcat(indices, [(vehicle_number, j) for j in i:1*grid_dim[2]:((grid_dim[1]-1)*grid_dim[2] + custom_mod(i, grid_dim[2]))])
            elseif  Φ[i][j] == WEST
                number_of_intersections_west = custom_mod(i,grid_dim[2])
                indices = vcat(indices, [(vehicle_number, j) for j in i:-1:(i-number_of_intersections_west+1)])
            end
        end
    end
    indices = sort(indices)
    @variable(model,t_access[indices]>=0)

    # Define the desired access times t_access_des
    Φ_combined = vcat(Φ...)
    d_combined = hcat(d...)
    @expression(model, t_access_des[i in indices], t_access[i] .+ L/v_avg)
    for i in 1:n
        # Extract the indices from vh1
        indices_vh1 = [(k, l) for (k, l) in indices if k == i]
        # Determine all the intersections vh1 is passing
        intersections_vh1 = last.(indices_vh1)
        if length(intersections_vh1) > 1
            if Φ_combined[i] == NORTH
                for j in 1:(length(intersections_vh1)-1)
                    t_access_des[(i, intersections_vh1[j])] = t_access[(i, intersections_vh1[j+1])] .+ L/v_avg
                end
                t_access_des[(i, intersections_vh1[end])] = t0 .+ d_combined[i]/v_avg
            elseif Φ_combined[i] == EAST
                for j in 1:(length(intersections_vh1)-1)
                    t_access_des[(i, intersections_vh1[j+1])] = t_access[(i, intersections_vh1[j])] .+ L/v_avg
                end
                t_access_des[(i, intersections_vh1[1])] = t0 .+ d_combined[i]/v_avg
            elseif Φ_combined[i] == SOUTH
                for j in 1:(length(intersections_vh1)-1)
                t_access_des[(i, intersections_vh1[j+1])] = t_access[(i, intersections_vh1[j])] .+ L/v_avg
                end
                t_access_des[(i, intersections_vh1[1])] = t0 .+ d_combined[i]/v_avg
            elseif  Φ_combined[i] == WEST
                for j in 1:(length(intersections_vh1)-1)
                    t_access_des[(i, intersections_vh1[j])] = t_access[(i, intersections_vh1[j+1])] .+ L/v_avg
                end
                t_access_des[(i, intersections_vh1[end])] = t0 .+ d_combined[i]/v_avg
            end
        elseif length(intersections_vh1) == 1
            t_access_des[(i, intersections_vh1[1])] = t0 + d_combined[i]/v_avg
        end
    end

    # Define the min access time t_access_min
    # The following three lines merely serve as an initialization of the t_access_min expression. We will change its values accordingly in the for loop below
    @expression(model, Δt1[i in indices], min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* L) .- v_avg) / max_a))
    @expression(model, Δt2[i in indices], max.(L ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0))
    @expression(model, t_access_min, Δt1 .+ Δt2 .+ 0*t_access) # The 0*t_access is there to ensure t_access_min is also of type AffEpr instead of Float64
    for i in 1:n
        # Extract the indices from vh1
        indices_vh1 = [(k, l) for (k, l) in indices if k == i]
        # Determine all the intersections vh1 is passing
        intersections_vh1 = last.(indices_vh1)
        if Φ_combined[i] == NORTH
            if length(intersections_vh1)>1
                for j in 1:(length(intersections_vh1)-1)
                    Δt1[(i,intersections_vh1[j])] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* L) .- v_avg) / max_a)
                    Δt2[(i,intersections_vh1[j])] = max.(L ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
                    t_access_min[(i,intersections_vh1[j])] = Δt1[(i,intersections_vh1[j])] + Δt2[(i,intersections_vh1[j])] + t_access[(i,intersections_vh1[j+1])]
                end
            end

            j = intersections_vh1[end]
            Δt1[(i,j)] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* d_combined[i]) .- v_avg) / max_a)
            Δt2[(i,j)] = max.(d_combined[i] ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
            t_access_min[(i,j)] = Δt1[(i,j)] + Δt2[(i,j)]

        elseif Φ_combined[i] == EAST
            if length(intersections_vh1)>1
                for j in 1:(length(intersections_vh1)-1)
                    Δt1[(i,intersections_vh1[j+1])] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* L) .- v_avg) / max_a)
                    Δt2[(i,intersections_vh1[j+1])] = max.(L ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
                    t_access_min[(i,intersections_vh1[j+1])] = Δt1[(i,intersections_vh1[j+1])] + Δt2[(i,intersections_vh1[j+1])] + t_access[(i,intersections_vh1[j])]
                end
            end

            j = intersections_vh1[1]
            Δt1[(i,j)] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* d_combined[i]) .- v_avg) / max_a)
            Δt2[(i,j)] = max.(d_combined[i] ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
            t_access_min[(i,j)] = Δt1[(i,j)] + Δt2[(i,j)]

        elseif Φ_combined[i] == SOUTH
            if length(intersections_vh1)>1
                for j in 1:(length(intersections_vh1)-1)
                    Δt1[(i,intersections_vh1[j+1])] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* L) .- v_avg) / max_a)
                    Δt2[(i,intersections_vh1[j+1])] = max.(L ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
                    t_access_min[(i,intersections_vh1[j+1])] = Δt1[(i,intersections_vh1[j+1])] + Δt2[(i,intersections_vh1[j+1])] + t_access[(i,intersections_vh1[j])]
                end
            end

            j = intersections_vh1[1]
            Δt1[(i,j)] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* d_combined[i]) .- v_avg) / max_a)
            Δt2[(i,j)] = max.(d_combined[i] ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
            t_access_min[(i,j)] = Δt1[(i,j)] + Δt2[(i,j)]

        elseif  Φ_combined[i] == WEST
            if length(intersections_vh1)>1
                for j in 1:(length(intersections_vh1)-1)
                    Δt1[(i,intersections_vh1[j])] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* L) .- v_avg) / max_a)
                    Δt2[(i,intersections_vh1[j])] = max.(L ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
                    t_access_min[(i,intersections_vh1[j])] = Δt1[(i,intersections_vh1[j])] + Δt2[(i,intersections_vh1[j])] + t_access[(i,intersections_vh1[j+1])]
                end
            end

            j = intersections_vh1[end]
            Δt1[(i,j)] = min.(((v_max .- v_avg) / max_a), (sqrt.(v_avg.^2 .+ 2 .* max_a .* d_combined[i]) .- v_avg) / max_a)
            Δt2[(i,j)] = max.(d_combined[i] ./ v_max .- (v_max^2 .- v_avg^2) / (2 * max_a * v_max), 0)
            t_access_min[(i,j)] = Δt1[(i,j)] + Δt2[(i,j)]
        end
    end

    # Speed limit and maximum acceleration constraints
    @constraint(model, [i in indices], t_access[i] .>= t_access_min[i])

    # Safety gap on the same movement constraint:
    for i in 1:n
        for j in 1:n
            if ((i != j) && (Φ_combined[i] == Φ_combined[j]))
                # Extract the indices from vh1 and vh2
                indices_vh1 = [(k, l) for (k, l) in indices if k == i]
                indices_vh2 = [(k, l) for (k, l) in indices if k == j]
                # Determine all the intersections vh1 and vh2 are passing
                intersections_vh1 = Set(last.(indices_vh1))
                intersections_vh2 = Set(last.(indices_vh2))
                intersection_overlap = intersect(intersections_vh1, intersections_vh2)
                if length(intersections_vh1) != length(intersections_vh2)
                    if Φ_combined[i] == NORTH
                        if maximum(intersections_vh1) >= maximum(intersections_vh2)
                            @constraint(model,[k in intersection_overlap], t_access[(i,k)] - t_access[(j,k)] >= tgap1)
                        end
                    elseif Φ_combined[i] == EAST
                        if minimum(intersections_vh1) >= minimum(intersections_vh2)
                            @constraint(model,[k in intersection_overlap], t_access[(i,k)] - t_access[(j,k)] >= tgap1)
                        end
                    elseif Φ_combined[i] == SOUTH
                        if minimum(intersections_vh1) >= minimum(intersections_vh2)
                            @constraint(model,[k in intersection_overlap], t_access[(i,k)] - t_access[(j,k)] >= tgap1)
                        end
                    elseif  Φ_combined[i] == WEST
                        if maximum(intersections_vh1) >= maximum(intersections_vh2)
                            @constraint(model,[k in intersection_overlap], t_access[(i,k)] - t_access[(j,k)] >= tgap1)
                        end
                    end
                else
                    if d_combined[i] >= d_combined[j]
                        @constraint(model,[k in intersection_overlap], t_access[(i,k)] - t_access[(j,k)] >= tgap1)
                    end
                end
            end
        end
    end

    # Safety gap on different movements constraint:
    for i in 1:n
        for j in 1:n
            if (((Int(Φ_combined[i]) - Int(Φ_combined[j])) % 2 != 0) && (j>i))
                # Extract the indices from vh1 and vh2
                indices_vh1 = [(k, l) for (k, l) in indices if k == i]
                indices_vh2 = [(k, l) for (k, l) in indices if k == j]
                # Determine all the intersections vh1 and vh2 are passing
                intersections_vh1 = Set(last.(indices_vh1))
                intersections_vh2 = Set(last.(indices_vh2))
                intersection_overlap = intersect(intersections_vh1, intersections_vh2)
                if length(intersection_overlap)>0
                    b1 = @variable(model, binary = true)
                    @constraint(model,[k in intersection_overlap], (t_access[(i,k)] - t_access[(j,k)] + m * b1) >= tgap2)
                    @constraint(model,[k in intersection_overlap], (t_access[(j,k)] - t_access[(i,k)] + m * (1-b1)) >= tgap2 )
                end
            end
        end
    end

    # Discontinuity in Cost Function J1
    @variable(model, t_access_max >= t0)
    @constraint(model, t_access_max .>= t_access)

    # Discontinuity in Cost Function J2
    @variable(model, Δt_access_abs[i in indices] >= t0)
    @constraint(model, [i in indices], Δt_access_abs[i] ≥ t_access[i] - t_access_des[i])
    @constraint(model, [i in indices], Δt_access_abs[i] ≥ -(t_access[i] - t_access_des[i]))

    # Linear combination of objective functions J1 and J2
    @objective(model, Min, w1 * t_access_max + w2 * sum(Δt_access_abs))

    optimize!(model)
    println(solution_summary(model))
    return model
end

grid_rowsize = 2 # Number of rows in the intersection grid
grid_columnsize = 2 # Number of columns in the intersection grid
grid_dim = [grid_rowsize, grid_columnsize]
d1 = [100 200 300] 
d2 = [100 200 300] 
d3 = [100 200 300] 
d4 = [100 200 300]
d = [d1, d2, d3, d4] # Combine the distance vectors of all intersections into a matrix

n = sum(length, d)

ϕ1 = [WEST, SOUTH, EAST]
ϕ2 = [SOUTH, EAST, NORTH]
ϕ3 = [EAST, NORTH, WEST]
ϕ4 = [NORTH, WEST, SOUTH]
Φ = [ϕ1, ϕ2, ϕ3, ϕ4] # Combine the phases of all intersections into a matrix
L = 500 # The distance in meters between the different intersection points

"""2x2 grid:
Scenario 1:

d1 = [100 200 300] 
d2 = [100 200 300] 
d3 = [100 200 300] 
d4 = [100 200 300] 

ϕ1 = [WEST, SOUTH, EAST]
ϕ2 = [SOUTH, EAST, NORTH]
ϕ3 = [EAST, NORTH, WEST]
ϕ4 = [NORTH, WEST, SOUTH]

Scenario 2:

d1 = [600 700 800 900 1000 1100 1200 1300 1400]
d2 = [600 700 800 900 1000 1100 1200 1300 1400]
d3 = [600 700 800 900 1000 1100 1200 1300 1400]
d4 = [600 700 800 900 1000 1100 1200 1300 1400]

ϕ1 = [WEST, SOUTH, EAST, NORTH, EAST, SOUTH, NORTH, WEST, SOUTH]
ϕ2 = [SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH]
ϕ3 = [EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST]

3x2 grid:

Scenario 1:

d1 = [600 700 800 900 1000 1100 1200 1300 1400]
d2 = [600 700 800 900 1000 1100 1200 1300 1400]
d3 = [600 700 800 900 1000 1100 1200 1300 1400]
d4 = [600 700 800 900 1000 1100 1200 1300 1400]
d5 = [600 700 800 900 1000 1100 1200 1300 1400]
d6 = [600 700 800 900 1000 1100 1200 1300 1400]

ϕ1 = [WEST, SOUTH, EAST, NORTH, EAST, SOUTH, NORTH, WEST, SOUTH]
ϕ2 = [SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH]
ϕ3 = [EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST]
ϕ4 = [NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH]
ϕ5 = [NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST, NORTH]
ϕ6 = [WEST, NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST]

Scenario 2:

d1 = [100 200 300 400] 
d2 = [100 200 300 400] 
d3 = [100 200 300 400] 
d4 = [100 200 300 400] 
d5 = [100 200 300 400] 
d6 = [100 200 300 400] 

ϕ1 = [WEST, SOUTH, EAST, NORTH]
ϕ2 = [SOUTH, EAST, NORTH, WEST]
ϕ3 = [EAST, NORTH, WEST, SOUTH]
ϕ4 = [NORTH, WEST, SOUTH, EAST]
ϕ5 = [NORTH, EAST, SOUTH, WEST]
ϕ6 = [WEST, NORTH, EAST, SOUTH]

3x3 grid:

Scenario 1:

d1 = [100 200 300 400 500 600 700 800 900] 
d2 = [100 200 300 400 500 600 700 800 900] 
d3 = [100 200 300 400 500 600 700 800 900] 
d4 = [100 200 300 400 500 600 700 800 900] 
d5 = [100 200 300 400 500 600 700 800 900] 
d6 = [100 200 300 400 500 600 700 800 900]
d7 = [100 200 300 400 500 600 700 800 900]
d8 = [100 200 300 400 500 600 700 800 900]
d9 = [100 200 300 400 500 600 700 800 900]

ϕ1 = [WEST, SOUTH, EAST, NORTH, EAST, SOUTH, NORTH, WEST, SOUTH]
ϕ2 = [SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH]
ϕ3 = [EAST, NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST]
ϕ4 = [NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH]
ϕ5 = [NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST, NORTH]
ϕ6 = [WEST, NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST]
ϕ7 = [NORTH, WEST, SOUTH, EAST, NORTH, WEST, SOUTH, EAST, NORTH]
ϕ8 = [NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST, NORTH]
ϕ9 = [WEST, NORTH, EAST, SOUTH, WEST, NORTH, EAST, SOUTH, WEST]"""

# 0,1 and 1,0 and 0.5,0.5 (3 different optimizations)
w1 = 0.5 # Weight for objective function 1
w2 = 0.5 # weight for objective function 2

t0 = 0.0
v_avg = 50000 / 3600
v_max = 80000 / 3600
max_a = 3
avg_a = 1

tgap1 = 2
tgap2 = 6.5

m = 2500 # This is the m-value used in the big-M method when dealing with OR constraints

model = main(Φ, d, w1, w2, t0, v_avg, v_max, avg_a, max_a, grid_dim, L)

# Visualize the solutions using a dataframe
result = [(model[:t_access].axes[1][i], 
value.(model[:t_access]).data[i], 
value.(model[:t_access_des]).data[i], 
value.(model[:t_access_min]).data[i]) 
for i in 1:length(value.(model[:t_access]).data)]
header = ["(Vehicle, Intersection)", "t_access", "t_access_des", "t_access_min"]
df = DataFrame(result,header);
println(df)

println("\nThe maximum t_access value is: $(maximum(value.(model[:t_access])))")
