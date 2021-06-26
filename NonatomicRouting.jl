using Plots
using JuMP
using Ipopt
using Polynomials


"""
Contains static information about a transportation network. 
"""
struct RoutingGame{T<:AbstractFloat, U<:Int}
    # Edges of the digraph
    edges::Vector{Tuple{U, U}}
    # Marginal cost along each edge as a function of its flow
    marginalcosts::Vector{Polynomial{T}}
    # Total cost =  flow[edge] * costs[edge]
    totalcosts::Vector{Polynomial{T}}
    # Demand source and sink nodes
    stpairs::Vector{Tuple{U, U}}
    # Quantity of demand for each s-t pair
    traffic::Vector{T}
    # Paths (lists of lists of edge indices) that each category of demand can use
    paths::Vector{Vector{Vector{U}}}
end


"""
    RoutingGame(edges, marginalcosts, stpairs, traffic)

Construct a `RoutingGame` instance using minimal input.
"""
function RoutingGame(edges::Vector{Tuple{U, U}},
                     marginalcosts::Vector{Polynomial{T}},
                     stpairs::Vector{Tuple{U, U}},
                     traffic::Vector{T}) where {T<:AbstractFloat, U<:Int}
    paths = [findpaths(edges, s, t) for (s, t) in stpairs]
    totalcosts = marginalcosts .* Polynomial([0, 1])
    return RoutingGame{T, U}(edges, marginalcosts, totalcosts, stpairs, traffic, paths)
end


"""
    findstpath(edges, s, t)

Find an `s`-`t` path in the graph `edges`. 
"""
function findstpath(edges::Vector{Tuple{Int, Int}}, s::Int, t::Int)::Vector{Int}
    n_nodes = mapreduce(maximum, max, edges)
    visited = falses(n_nodes)
    
    if s==t
        return Int[]
    end
    
    nextedge = findfirst(isequal((s, t)), edges)
    
    if nextedge===nothing
        nextedge = findfirst([a==s && !visited[b] for (a, b) in edges])
        @show nextedge
        return vcat(nextedge, findstpath(edges, edges[nextedge][2], t))
    else
        @show nextedge
        return [nextedge]
    end
end


"""
    findpaths(edges, s, t)

Find all the `s`-`t` paths in the graph defined by `edges`.
"""
function findpaths(edges::Vector{Tuple{Int, Int}}, s::Int, t::Int)::Vector{Vector{Int}}
    n_nodes = mapreduce(maximum, max, edges)
    
    if s==t
        return [Int[]]
    else
        edgesout = findall(x->x[1]==s, edges)
        return [(vcat(j, p)) for j in edgesout for p in findpaths(edges, edges[j][2], t)]
    end
end


"""
    solveroutinggame(game)

Solve a `RoutingGame` instance. Compute the socially optimal flow and 
equilibrium flow and print the price of anarchy. Returns two `NamedTuple`s,
each one containing `.pathflows` and `.edgeflows`.
"""
function solveroutinggame(game::RoutingGame)
    n_edges = length(game.edges)
    n_commodities = length(instance.paths)
    
    model = Model(Ipopt.Optimizer)
    @variable(model, pathwiseflow[c=1:n_commodities, p=1:length(game.paths[c])] ≥ 0)
    
    # Flow in each edge equals sum of pathwise flows containing that edge
    @expression(model, flow[e = 1:n_edges],
                            sum(e in path ? pathwiseflow[c, p] : 0
                                for c in 1:n_commodities
                                for (p, path) in enumerate(game.paths[c])))
    
    # Flow in each edge is nonnegative
    @constraint(model, FlowNonneg[e = 1:n_edges], flow[e] ≥ 0)
    
    # Sum of pathwise flows assoc with each commodity equals total demand
    @constraint(model, TrafficVolume[c = 1:n_commodities],
                sum(pathwiseflow[c, p] for p in 1:length(game.paths[c])) == game.traffic[c])


    # Total cost. Need to "manually" unpack the polynomial for JuMP to parse.
    @expression(model, TotalCost,
        sum(b * flow[e]^(k-1) for e in 1:n_edges
                              for (k, b) in enumerate(game.totalcosts[e])))
        # = sum(game.totalcosts[e](flow[e]) for e in 1:n_edges)
    
    @objective(model, Min, TotalCost)
    
    optimize!(model)
    systemopt = (pathflows = value.(pathwiseflow), edgeflows = value.(flow))
    poa = objective_value(model)
    
    
    # Optimize again using the potential function, which gives the equilibrium flow
    @expression(model, Potential,
        sum(b * flow[e]^(k-1) for e in 1:n_edges
                              for (k, b) in enumerate(integrate(game.marginalcosts[e]))))
    @objective(model, Min, Potential)
    
    optimize!(model)
    equilibrium = (pathflows = value.(pathwiseflow), edgeflows = value.(flow))
    poa = value.(TotalCost) / poa
    println("Price of anarchy: $poa")
        
    return systemopt, equilibrium
end


"""
    shownet(game, node_coords, flow=nothing)

Display a nonatomic routing game with nodes at the given coordinates
and, optionally, a flow, which is a `NamedTuple` containing `.pathflows` and `.edgeflows`.
"""
function shownet(game::RoutingGame, node_coords, flow=nothing)
    local pl
    pl = plot(ticks=nothing, border=nothing)
    
    for (i, e) in enumerate(game.edges)
        x = node_coords[e[1], 1]
        y = node_coords[e[1], 2]
        u = node_coords[e[2], 1] - x
        v = node_coords[e[2], 2] - y
        quiver!(pl, [x], [y],
                quiver=([u], [v]),
                arrow=Plots.Arrow(:open, :none, 344, 284),
                color=:gray90, lw=5, label="", ms=8)
        quiver!(pl, [x], [y],
                quiver=(0.5*[u], 0.5*[v],),
                arrow=Plots.Arrow(:open, :head, 344, 284),
                color=:gray90, lw=6, label="", ms=8)

        if !(flow===nothing)
            MC = round(game.marginalcosts[i](flow.edgeflows[i]), digits=4)
            SC = round(game.totalcosts[i](flow.edgeflows[i]), digits=4)
            fl = round(flow.edgeflows[i], digits=4)
            annotate!(pl, [(x+0.5*u, y+0.5*v, Plots.text("MC: $MC\nSC: $SC\nFlow: $fl", 7, :left))])
        end
    end
    
    
    colors = [:cornflowerblue, :crimson, :olivedrab, :gold]
    for (i, (s, t)) in enumerate(game.stpairs)
        annotate!(node_coords[s, 1]-0.05, node_coords[s, 2]+0.05, Plots.text("(-$(traffic[i]))", colors[i], 9))
        annotate!(node_coords[t, 1]-0.05, node_coords[t, 2]+0.05, Plots.text("($(traffic[i]))", colors[i], 9))
    end
    
    if !(flow===nothing)
        for c in 1:length(game.paths)
            for (p, path) in enumerate(game.paths[c])
                for edge_idx in path
                    e = game.edges[edge_idx]
                    x = node_coords[e[1], 1] + .01c
                    y = node_coords[e[1], 2] + .01c
                    u = node_coords[e[2], 1] - x + .01c
                    v = node_coords[e[2], 2] - y + .01c
                    quiver!(pl, [x], [y],
                            quiver=([u], [v],),
                            arrow=Plots.Arrow(:open, :none, 344, 284),
                            color=colors[c], lw=10*flow.pathflows[c, p], label="", ms=8)
                end
            end
        end
    end
        
    
    scatter!(pl, node_coords[:, 1], node_coords[:, 2], mc=:black, msw=0, ms=10, label="")
    annotate!(pl, [(x, y, Plots.text(string(i), :white, 9)) for (i, (x, y)) in enumerate(eachrow(node_coords))])
    
    
    return pl
end


"""
    pathcosts(game, flow)

Gives the total marginal cost along each path, and prints in a format that can be inspected
to determined if `flow` is an equilibrium. Flow should be a `NamedTuple` or `struct` such that
`flow.edgeflows` gives the flow in each edge.
"""
function pathcosts(game::RoutingGame, flow)
    res = []
    for c in 1:length(game.paths)
        @show c
        thiscommoditycosts = []
        for (p, path) in enumerate(game.paths[c])
            pathmarginalcost = sum(game.marginalcosts[edge_idx](flow.edgeflows[edge_idx]) for edge_idx in path)
            pathfl = equilibrium.pathflows[c, p]
            pathfl = pathfl > 1e-5 ? pathfl : "0"
            println("  Path: $(game.edges[game.paths[c][p]])")
            println("    Flow      : $pathfl")
            println("    Marg. cost: $pathmarginalcost")
            push!(thiscommoditycosts, pathmarginalcost)
        end
        push!(res, thiscommoditycosts)
    end

    return res
end