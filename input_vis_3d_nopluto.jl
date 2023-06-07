using GLMakie
using TetGen
using StructArrays
using StaticArrays
#using PlutoUI
using GeometryBasics
using GridVisualize
using SimplexGridFactory
using ExtendableGrids
using Printf

#=
"""
### Functions that I stole from the source code
#### https://github.com/j-fu/GridVisualize.jl/blob/main/src/makie.jl
#### (they are not exported from the module)
"""
=#

function makescene_grid(ctx)
	#Kommentiere den Block ein für Farbbalken
    XMakie = ctx[:Plotter]
    GL = XMakie.GridLayout(ctx[:figure])
    GL[1, 1] = ctx[:scene]

	
    ncol = length(ctx[:cmap])
    nbcol = length(ctx[:cmap])

	# fontsize=0.5*ctx[:fontsize],ticklabelsize=0.5*ctx[:fontsize]
    if ctx[:colorbar] == :vertical
        GL[1, 2] = XMakie.Colorbar(
            ctx[:figure];
            colormap = XMakie.cgrad(ctx[:cmap]; categorical = true),
            limits = (1, ncol),
            width = 15,
        )
        GL[1, 3] = XMakie.Colorbar(
            ctx[:figure];
            colormap = XMakie.cgrad(ctx[:bcmap]; categorical = true),
            limits = (1, nbcol),
            width = 15,
        )
    elseif ctx[:colorbar] == :horizontal
        GL[2, 1] = XMakie.Colorbar(
            ctx[:figure];
            colormap = XMakie.cgrad(ctx[:cmap]; categorical = true),
            limits = (1, ncol),
            heigth = 15,
            vertical = false,
        )
        GL[3, 1] = XMakie.Colorbar(
            ctx[:figure];
            colormap = XMakie.cgrad(ctx[:bcmap]; categorical = true),
            limits = (1, nbcol),
            heigth = 15,
            vertical = false,
        )
    end
	
    GL
end


function scene_interaction(
    update_scene,
    scene,
    XMakie,
    switchkeys::Vector{Symbol} = [:nothing],
)
# Check if pixel position pos sits within the scene
    function _inscene(scene, pos)
        area = scene.px_area[]
        pos[1] > area.origin[1] &&
            pos[1] < area.origin[1] + area.widths[1] &&
            pos[2] > area.origin[2] &&
            pos[2] < area.origin[2] + area.widths[2]
    end

    # Initial active switch key is the first in the vector passed
    activeswitch = Observable(switchkeys[1])

    # Handle mouse position within scene
    mouseposition = Observable((0.0, 0.0))

    XMakie.on(scene.events.mouseposition) do m
        mouseposition[] = m
        false
    end

    # Set keyboard event callback
    XMakie.on(scene.events.keyboardbutton) do buttons
        if _inscene(scene, mouseposition[])
            # On pressing a switch key, pass control
            for i = 1:length(switchkeys)
                if switchkeys[i] != :nothing &&
                   XMakie.ispressed(scene, getproperty(XMakie.Keyboard, switchkeys[i]))
                    activeswitch[] = switchkeys[i]
                    update_scene(0, switchkeys[i])
                    return true
                end
            end

            # Handle change values via up/down control
            if XMakie.ispressed(scene, XMakie.Keyboard.up)
                update_scene(1, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.down)
                update_scene(-1, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.page_up)
                update_scene(10, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.page_down)
                update_scene(-10, activeswitch[])
                return true
            end
        end
        return false
    end
end

scenekwargs(ctx) = Dict(
    #:xticklabelsize => 0.5*ctx[:fontsize],
    #:yticklabelsize => 0.5*ctx[:fontsize],
    #:zticklabelsize => 0.5*ctx[:fontsize],
    #:xlabelsize => 0.5*ctx[:fontsize],
    #:ylabelsize => 0.5*ctx[:fontsize],
    #:zlabelsize => 0.5*ctx[:fontsize],
    #:xlabeloffset => 20,
    #:ylabeloffset => 20,
    #:zlabeloffset => 20,
    :titlesize => ctx[:fontsize],
)

function makeaxis3d(ctx)
    XMakie = ctx[:Plotter]
    if ctx[:scene3d] == :LScene
        # "Old" LScene with zoom-in functionality
        XMakie.LScene(ctx[:figure])
		println("Option 1")
    else
        # "New" Axis3 with prospective new stuff by Julius.
        XMakie.Axis3(
            ctx[:figure];
            aspect = (1,1,1), #:data,
            viewmode = :fit,
            elevation = ctx[:elev] * π / 180,
            azimuth   = ctx[:azim] * π / 180,
            perspectiveness = ctx[:perspectiveness],
            #title = map(data -> data.t, ctx[:data]),
            scenekwargs(ctx)...,
        )
    end
end



function xyzminmax(grid::ExtendableGrid)
    coord = grid[Coordinates]
    ndim = size(coord, 1)
    xyzmin = zeros(ndim)
    xyzmax = ones(ndim)
    for idim = 1:ndim
        @views mn, mx = extrema(coord[idim, :])
        xyzmin[idim] = mn
        xyzmax[idim] = mx
    end
    xyzmin, xyzmax
end

#="""
#### unimport support functions
#### e.g. the builder is created
"""
=#

function triangulation_of_cube_with_subregion()
    
    builder=SimplexGridBuilder(Generator=TetGen)
    
    p1=point!(builder,0,0,0)
    p2=point!(builder,0,0,1)
    p3=point!(builder,0,1,0)
    p4=point!(builder,0,1,1)
    p5=point!(builder,1,0,0)
    p6=point!(builder,1,0,1)
    p7=point!(builder,1,1,0)
    p8=point!(builder,1,1,1)
    #p5=point!(builder,-1,2)

	facetregion!(builder,1)
    facet!(builder,p1,p2,p3)
	facet!(builder,p2,p3,p4)
    facet!(builder,p1,p2,p5)
	facet!(builder,p2,p5,p6)
    facet!(builder,p1,p3,p5)
    facet!(builder,p3,p5,p7)
	
	facetregion!(builder,2)	
    facet!(builder,p2,p4,p6)
    facet!(builder,p4,p6,p8)
    facet!(builder,p3,p4,p7)
    facetregion!(builder,2)	
    facet!(builder,p4,p7,p8)
    facet!(builder,p5,p6,p7)
    facet!(builder,p6,p7,p8)

	facetregion!(builder,3)
	facet!(builder,p2,p3,p5) #p6

	cellregion!(builder,1)
    maxvolume!(builder, 1)
    regionpoint!(builder, 0.2,0.2,0.1)

	cellregion!(builder,2)
    maxvolume!(builder, 1)
    regionpoint!(builder, 0.8,0.7,0.01)
	
    options!(builder,maxvolume=1.0)
	builder
end

# ╔═╡ fccb08c7-6707-478c-bc92-8de2a47ebce0
function triangulation_of_cube_with_innerboundary()
    
    builder=SimplexGridBuilder(Generator=TetGen)

	
    p1=point!(builder,0,0,0)
    p2=point!(builder,0,0,1)
    p3=point!(builder,0,1,0)
    p4=point!(builder,0,1,1)
    p5=point!(builder,1,0,0)
    p6=point!(builder,1,0,1)
    p7=point!(builder,1,1,0)
    p8=point!(builder,1,1,1)
	
    p9 =point!(builder,0.25,0.25,0.25)
    p10=point!(builder,0.25,0.25,0.50)
    p11=point!(builder,0.25,0.50,0.25)
    p12=point!(builder,0.25,0.50,0.50)
    p13=point!(builder,0.50,0.25,0.25)
    p14=point!(builder,0.50,0.25,0.50)
    p15=point!(builder,0.50,0.50,0.25)
    p16=point!(builder,0.50,0.50,0.50)
	

	
	#outer bnd
	facetregion!(builder,1)
    facet!(builder,p1,p2,p3)
	facet!(builder,p2,p3,p4)
    facet!(builder,p1,p2,p5)
	
	facet!(builder,p2,p5,p6)
    facet!(builder,p1,p3,p5)
    facet!(builder,p3,p5,p7)
		
    facet!(builder,p2,p4,p6)
	
	facetregion!(builder,2)
    facet!(builder,p4,p6,p8)
    facet!(builder,p3,p4,p7)
    #facetregion!(builder,3)	
    facet!(builder,p4,p7,p8)
    facet!(builder,p5,p6,p7)
    facet!(builder,p6,p7,p8)
	

	
	#inner bnd	
    facetregion!(builder,3)
	facet!(builder,p9,p10,p11)
	facet!(builder,p10,p11,p12)
    facet!(builder,p9,p10,p13)
	facet!(builder,p10,p13,p14)
    facet!(builder,p9,p11,p13)
    facet!(builder,p11,p13,p15)
	
	facetregion!(builder,4)	
    facet!(builder,p10,p12,p14)
    facet!(builder,p12,p14,p16)
    #facet!(builder,p11,p12,p15)	
    facet!(builder,p12,p15,p16)
    facet!(builder,p13,p14,p15)
    facet!(builder,p14,p15,p16)

	facetregion!(builder,3)	
    facet!(builder,p11,p12,p15)	
    
	
    #options!(builder,maxvolume=maxvol)
	builder
end

#=
"""
### (important support) functions I wrote
#### culminating in the boundaryplot3d function
"""
=#

function xyzminmax(builder::SimplexGridBuilder)
	
	return ([minimum(builder.pointlist.points[1,:]),
			 minimum(builder.pointlist.points[2,:]), 	                               minimum(builder.pointlist.points[3,:])], 
			[maximum(builder.pointlist.points[1,:]),
			 maximum(builder.pointlist.points[2,:]), 	                               maximum(builder.pointlist.points[3,:])])
end

function gridplot33d(ctx, TP::Type{MakieType}, ::Type{Val{3}}, grid)
    #this function is essentially = GridVisualize.gridplot!, but it is also influenced by angle changes etc
	
	make_mesh(pts, fcs) = GeometryBasics.Mesh(GridVisualize.meta(pts; normals = GridVisualize.normals(pts, fcs)), fcs)
	

    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)

    XMakie = ctx[:Plotter]
    xyzmin, xyzmax = xyzminmax(grid)
    xyzstep = (xyzmax - xyzmin) / 25

    function adjust_planes()
        ctx[:xplanes][1] = max(xyzmin[1], min(xyzmax[1], ctx[:xplanes][1]))
        ctx[:yplanes][1] = max(xyzmin[2], min(xyzmax[2], ctx[:yplanes][1]))
        ctx[:zplanes][1] = max(xyzmin[3], min(xyzmax[3], ctx[:zplanes][1]))
    end

    adjust_planes()

    if !haskey(ctx, :scene)
        ctx[:data] = Observable((
            g = grid,
            x = ctx[:xplanes][1],
            y = ctx[:yplanes][1],
            z = ctx[:zplanes][1],
            t = ctx[:title],
        ))

        ctx[:scene] = makeaxis3d(ctx)
        cmap = GridVisualize.region_cmap(nregions)
        ctx[:cmap] = cmap
        bcmap = GridVisualize.bregion_cmap(nbregions)
        ctx[:bcmap] = bcmap

        ############# Interior cuts
        # We draw a mesh for each color.
        if ctx[:interior]
            ctx[:celldata] = map(
                d -> GridVisualize.extract_visible_cells3D(
                    d.g,
                    (d.x, d.y, d.z);
                    primepoints = hcat(xyzmin, xyzmax),
                    Tp = Point3f,
                    Tf = GridVisualize.GLTriangleFace,
                ),
                ctx[:data],
            )

            ctx[:cellmeshes] =
                map(d -> [make_mesh(d[1][i], d[2][i]) for i = 1:nregions], ctx[:celldata])

            for i = 1:nregions
                XMakie.mesh!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:cellmeshes]);
                    color = cmap[i],
                    backlight = 1.0f0,
                )

                if ctx[:linewidth] > 0
                    XMakie.wireframe!(
                        ctx[:scene],
                        map(d -> d[i], ctx[:cellmeshes]);
                        color = :black,
                        strokecolor = :black,
                        strokewidth = ctx[:linewidth],
                        linewidth = ctx[:linewidth],
                    )
                end
            end
        end

        ############# Visible boundary faces
        ctx[:facedata] = map(
            d -> GridVisualize.extract_visible_bfaces3D(
                d.g,
                (d.x, d.y, d.z);
                primepoints = hcat(xyzmin, xyzmax),
                Tp = Point3f,
                Tf = GridVisualize.GLTriangleFace,
            ),
            ctx[:data],
        )

        ctx[:facemeshes] =
            map(d -> [make_mesh(d[1][i], d[2][i]) for i = 1:nbregions], ctx[:facedata])

        for i = 1:nbregions
            XMakie.mesh!(
                ctx[:scene],
                map(d -> d[i], ctx[:facemeshes]);
                color = bcmap[i],
                backlight = 1.0f0,
            )
            if ctx[:linewidth] > 0
                XMakie.wireframe!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:facemeshes]);
                    color = :black,
                    strokecolor = :black,
                    linewidth = ctx[:linewidth],
                )
            end
        end

        ############# Transparent outline

        if ctx[:outlinealpha] > 0.0
            ctx[:outlinedata] = map(
                d -> GridVisualize.extract_visible_bfaces3D(
                    d.g,
                    xyzmax;
                    primepoints = hcat(xyzmin, xyzmax),
                    Tp = Point3f,
                    Tf = GridVisualize.GLTriangleFace,
                ),
                ctx[:data],
            )
            ctx[:outlinemeshes] = map(
                d -> [make_mesh(d[1][i], d[2][i]) for i = 1:nbregions],
                ctx[:outlinedata],
            )

            for i = 1:nbregions
                XMakie.mesh!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:outlinemeshes]);
                    color = (bcmap[i], ctx[:outlinealpha]),
                    transparency = true,
                    backlight = 1.0f0,
                )
            end
        end

        ##### Interaction
        scene_interaction(ctx[:scene].scene, XMakie, [:z, :y, :x, :q]) do delta, key
            if key == :x
                ctx[:xplanes][1] += delta * xyzstep[1]
                ctx[:status][] = @sprintf("x=%.3g", ctx[:xplanes][1])
            elseif key == :y
                ctx[:yplanes][1] += delta * xyzstep[2]
                ctx[:status][] = @sprintf("y=%.3g", ctx[:yplanes][1])
            elseif key == :z
                ctx[:zplanes][1] += delta * xyzstep[3]
                ctx[:status][] = @sprintf("z=%.3g", ctx[:zplanes][1])
            elseif key == :q
                ctx[:status][] = " "
            end
            adjust_planes()

            ctx[:data][] = (
                g = grid,
                x = ctx[:xplanes][1],
                y = ctx[:yplanes][1],
                z = ctx[:zplanes][1],
                t = ctx[:title],
            )
        end

        ctx[:status] = Observable(" ")

        GridVisualize.add_scene!(ctx, makescene_grid(ctx))

    else
        ctx[:data][] = (
            g = grid,
            x = ctx[:xplanes][1],
            y = ctx[:yplanes][1],
            z = ctx[:zplanes][1],
            t = ctx[:title],
        )
    end
    reveal(ctx, TP)
end

function give_builder_vertsfaces(builder; Tp = SVector{3, Float32}, Tf = SVector{3, Int32}) #return the points and faces of the builder in the required format, works, but only for one region
	#builder = builder3
	nbregions = maximum(builder.facetregions)
	xyz = builder.pointlist.points
	tri = builder.facets	
	num_of_tris_in_region = [length(findall(x->x==i, builder.facetregions)) for i=1:length(tri)]
	faces = [Vector{Tf}(undef, num_of_tris_in_region[i]) for i=1:nbregions]
	verts = [Vector{Tp}(undef, size(xyz,2)) for i=1:nbregions]
	count = Int.(ones(nbregions))
	for j=1:nbregions
		for i in 1:size(xyz,2)
			verts[j][i] = Point(xyz[1,i], xyz[2,i], xyz[3,i])
			#push!(verts[j], Point(xyz[1,i], xyz[2,i], xyz[3,i]))
		end
	end
	
	for i in 1:length(tri)
		faces[builder.facetregions[i]][count[builder.facetregions[i]]] = TriangleFace(tri[i][1], tri[i][2], tri[i][3])
		count[builder.facetregions[i]] += 1
		#push!(faces[builder.facetregions[i]], TriangleFace(tri[i][1], tri[i][2], tri[i][3]))
	end
	verts, faces
end

function boundaryplot3d(builder)
	make_mesh(pts, fcs) = 
	GeometryBasics.Mesh(
		GridVisualize.meta(pts; normals = GridVisualize.normals(pts, fcs)), fcs
	)

	visualizer = GridVisualizer(; Plotter=GLMakie, show=true)
	ctx = visualizer[1,1]
	
    XMakie = ctx[:Plotter]
    xyzmin, xyzmax = xyzminmax(builder) #([-0.1, -0.1, -0.1], [1.1, 1.1, 1.1]) #xyzminmax(grid)
    #xyzmin, xyzmax = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]) #xyzminmax(grid)
    xyzstep = (xyzmax - xyzmin) / 25

    function adjust_planes()
        ctx[:xplanes][1] = max(xyzmin[1], min(xyzmax[1], ctx[:xplanes][1]))
        ctx[:yplanes][1] = max(xyzmin[2], min(xyzmax[2], ctx[:yplanes][1]))
        ctx[:zplanes][1] = max(xyzmin[3], min(xyzmax[3], ctx[:zplanes][1]))
    end

    adjust_planes()

	### colors and transparencies
	if(size(builder.regionpoints,2) > 0) 
		nregions  = maximum(builder.regionnumbers)
	else
		nregions = 0 #	minimum(maximum, ) #num_cellregions(grid)
	end
    nbregions = maximum(builder.facetregions) #num_bfaceregions(grid)
	ctx[:cmap] = GridVisualize.region_cmap(nregions) # = cmap
	ctx[:bcmap] = GridVisualize.bregion_cmap(nbregions) # = bcmap
	alpha = 0.2*ones(nbregions) #reverse([0.1, 0.1, 0.1, 0.1])
	ctx[:scene] = makeaxis3d(ctx)

	
	### create the boundary mesh, and then its observable
	ctx[:celldata] = map(
		d -> give_builder_vertsfaces(builder; Tp = Point3f,
			Tf = GridVisualize.GLTriangleFace),
		Observable([1]),
	)

	ctx[:cellmeshes] = map(d -> [make_mesh(d[1][i], d[2][i]) for i = 1:nbregions], ctx[:celldata])
	
	
	
	### plot the boundary (wireframe and colors)
	for i = 1:nbregions
		XMakie.mesh!(
			ctx[:scene], map(d -> d[i], ctx[:cellmeshes]);
			color=(ctx[:bcmap][i], alpha[i]), transparency = true)
	
		XMakie.wireframe!(
			ctx[:scene], map(d -> d[i], ctx[:cellmeshes]);
			color = :black, strokecolor = :black, strokewidth = ctx[:linewidth],
			linewidth = ctx[:linewidth],
		)
	end
	
	
	### scatter of regionpoints
	if(size(builder.regionpoints,2) > 0)
		XMakie.scatter!(ctx[:scene], Matrix(builder.regionpoints), color=ctx[:cmap][1:nregions])
	end
			
	##### Interaction
	
	scene_interaction(ctx[:scene].scene, XMakie, [:z, :y, :x, :q]) do delta, key
		#println("key = ", key)
		if key == :x
			ctx[:xplanes][1] += delta * xyzstep[1]
			ctx[:status][] = @sprintf("x=%.3g", ctx[:xplanes][1])
		elseif key == :y
			ctx[:yplanes][1] += delta * xyzstep[2]
			ctx[:status][] = @sprintf("y=%.3g", ctx[:yplanes][1])
		elseif key == :z
			ctx[:zplanes][1] += delta * xyzstep[3]
			ctx[:status][] = @sprintf("z=%.3g", ctx[:zplanes][1])
		elseif key == :q
			ctx[:status][] = " "
		end
		adjust_planes()
	end
    ctx[:status] = Observable(" ")
    GridVisualize.add_scene!(ctx, makescene_grid(ctx))
		
    reveal(ctx, plottertype(ctx[:Plotter]))
	
end

### playground


builder_SR = triangulation_of_cube_with_subregion()
builder_IB = triangulation_of_cube_with_innerboundary()
grid_SR = simplexgrid(builder_SR)
grid_IB = simplexgrid(builder_IB)

function execute_gridplot(grid)
	visualizer = GridVisualizer(; Plotter=GLMakie, show=true)
	ctx = visualizer[1,1]
	gridplot33d(ctx, plottertype(ctx[:Plotter]), Val{dim_space(grid)}, grid)
end



