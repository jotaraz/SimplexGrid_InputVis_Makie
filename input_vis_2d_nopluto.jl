using GLMakie
using CairoMakie
using Colors


function plot_input_builder(builder; XMakie=GLMakie)
	fig = XMakie.Figure()
	Axis(fig[1, 1])
	
	r = builder.pointlist.points
	f = builder.facets
	n = size(f,1)  # number of faces

	#different colors we'll use for the faces, each boundaryregion gets its own color
	cols = distinguishable_colors(maximum(builder.facetregions), [RGB(1.0, 0.0, 0.0), RGB(0.0, 1.0, 0.0), RGB(0.0, 0.0, 1.0)]) 
	for j in 1:n
		linesegments!([r[1,f[j][1]], r[1,f[j][2]]], [r[2,f[j][1]], r[2,f[j][2]]], color=cols[builder.facetregions[j]])
		
	end

	XMakie.Colorbar(
			fig[1, 2];
			label="boundary regions",
			colormap = XMakie.cgrad(cols; categorical=true),
			limits = (1,n),
	) #colorbar for the boundary

	#only consider regionpoints, if there are regionpoints
	if(size(builder.regionpoints, 2) > 0) 
		#different colors we'll use for the cellregion, each cellregion gets its own color
		cols = distinguishable_colors(maximum(builder.regionnumbers), [RGB(1.0, 0.0, 0.0), RGB(0.0, 1.0, 0.0), RGB(0.0, 0.0, 1.0)])[builder.regionnumbers]  
		scatter!(builder.regionpoints[1,:], builder.regionpoints[2,:], color=cols)
		
		XMakie.Colorbar(
				fig[1, 3];
				label="cell regions",
				colormap = XMakie.cgrad(cols; categorical=true),
				limits = (1,maximum(builder.regionnumbers)),
		) #colorbar for the cellregions
	end

	fig
end
