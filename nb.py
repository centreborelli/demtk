# # DEMTK : A Toolkit for Digital Elevation Models
#
# A DEM is a two-dimensional array of numbers.  Each number
# represents a height in meters.  The positions inside the array are mapped into
# geographical coordinates (typically UTM).  The mapping between positions in
# the array and geographical coordinates is called georeferencing.
# Thus, a geo-referenced DEM can be interpreted as a cloud of points in 3D
# space.  For this document we forget about georeferencing and consider DEM
# just as raster images.
#

# ## 0. Sources of DEM
#
# In this notebook we will work with two concrete DEM:
#
# * DEM **x**, an SRTM raster of the area around Fuji at a resolution of 90m
# per pixel
# * DEM **y**, the projection of a high-resolution Lidar point cloud of an
# urban region, at 2m per pixel
#
# Since they are raster images, we can display them directly, just by
# rescaling their values to [0,255]:

import demtk  # DEM toolkit
import iio    # image input/output

x = iio.read("i/fuji.tif")[:,:,0]
y = iio.read("i/terrassa.tif")[:,:,0]
# NOTE: The rest of the notebook uses the global variables **x** and **y**,
# you can replace them here to test all the methods with other data.

# Show "x" and "y with manually chosen value scalings
iio.gallery([
	x,           # displaying "x" alone clips heights to [0,255]
	255*x/4000,  # rescale [0,4000] to [0,255] for a complete display
	y-250        # to display "y" we just the values to a visible range
])




# ## 1. DEM rendering
#
# How to produce an image that visually displays the contents of a DEM

iio.gallery([
	demtk.render(x),
   	demtk.render(x,p=0.1),
	demtk.render(y),
    demtk.render(y,p=5)
	])

# ### 1.1. Color palette

# ### 1.2. Hillshading

# ### 1.3. Shadows

# ### 1.4. Curvature map

# ### 1.5. Ambient occlusion

# ### 1.6. Combined rendering


# ## 2. DEM denoising
#
# How to clean-out and smooth the values of a DEM

# ### 2.1. morphological operations

# ### 2.2. cc filter

# ### 2.3. bilateral filter



# ## 3. DEM interpolation
#
# How to fill-in missing data in a DEM

# ### 3.1. cc border, avg, min, min5pc

# ### 3.2. poisson, biharmonic, TV, with Dirichlet boundary conditions

# ### 3.3. neumann boundaries at selected points


# ## 4. DEM registration
#
# How to correct slightly the position of a DEM so that it fits well another

# ### 4.1. apply shift (manually selected)

# ### 4.2. multiscale correlation

# ### 4.3. phase correlation in the fourier domain


# ## 5. DEM comparison
#
# Criteria to compare two different DEM (e.g., a ground truth and a computed
# DEM)

# ### 5.1. difference of registered images, signed palette

# ### 5.2. filtering of difference

# ### 5.3. metrics, statistics and histograms


# ## 6. DEM fusion
#
# How to combine several DEM into one

# ### 6.1. nanavg

# ### 6.2. nanmed

# ### 6.3. xmedians

# ### 6.4. cnt

# ### 6.5. filtered fusion (by iqd, cnt, ...)


# ## 7. DEM flattening
#
# Or how to transform a point cloud, or a mesh, into a DEM

# ### 7.1. plyflatten

# ### 7.2. multiscale plyflatten

# ### 7.3. splatting


# ## 8. DEM elevation
#
# How to convert a DEM into a 3D mesh

# ### 8.1. isolated points

# ### 8.2. local mesh

# ### 8.3. refined mesh heuristics


# ## 9. DEM colorization
#
# How to color the points of a DEM (or of the associated mesh) according to
# an image or a set of images

# ### 9.1. naive colorization

# ### 9.2. shadowed colorization

# ### 9.3. orthophoto merging (color median)

# ### 9.4. orthophoto merging by Poisson

# ### 9.5. orthophoto merging by osmosis

# ### 9.6. pansharpening and other multispectral issues






import demtk
from numpy import fmin, newaxis


s = demtk.qauto(demtk.render_shading(x),p=1)
a = demtk.qauto(fmin(0,demtk.filter_riesz(x,-0.5)),p=1)
sa = demtk.qauto((s*1.0*a)**0.5, p=1)
p = demtk.render_palette_dem(x)
iio.gallery([
    demtk.qauto(x),
    a,
    s,
    sa,
    p,
    demtk.qauto(sa[:,:,newaxis]*1.0*p,p=0.1),
    demtk.render(x,p=1),
    demtk.render(x,p=0),
    demtk.render(x,p=0.1),
    demtk.renderclean(x)
])







iio.gallery([
    demtk.renderclean(x),
    demtk.render(x),
    demtk.render_palette_dem(x),
    demtk.render_shading(x),
    demtk.render_lssao(x),
    fmin(0,demtk.render_lssao(x))
], qauto=True)

iio.display(demtk.colorize(x))







# vim:set tw=77 filetype=python spell spelllang=en:
