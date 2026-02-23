# a toolkit for processing digital elevation maps (DEM)

#def grid_adjacency(w, h):
#	from scipy.sparse import eye, kron
#	x = eye(w - 1, w, 1) - eye(w  -1, w)  # path of length W
#	y = eye(h - 1, h, 1) - eye(h - 1, h)  # path of length H
#	B = kron(x, y)                        # kronecker product
#	return B

def grid_incidence(h, w):
	from scipy.sparse import eye, kron, vstack
	x = eye(w - 1, w, 1) - eye(w  -1, w)             # path of length W
	y = eye(h - 1, h, 1) - eye(h - 1, h)             # path of length H
	B = vstack([ kron(eye(h),x) , kron(y,eye(w)) ])  # kronecker sum
	return B

def grid_laplacian(h, w):
	B = grid_incidence(h, w)
	L = - B.T @ B
	return L

def grid_structuring(h, w):
	L = grid_laplacian(h, w)
	E = (L != 0).astype(float)
	return E

def graph_dilation(E, x):
	from scipy.sparse import diags
	y = (diags(x.flatten())@E).max(axis=0).toarray().squeeze()
	return y

def graph_median(E, x):
	from sklearn.utils.sparsefuncs import csc_median_axis_0 as med
	from scipy.sparse import diags
	y = med((diags(x.flatten())@E).tocsc()).squeeze()
	print("ERROR: this code does not work")
	return y

def graph_erosion(E, x):
	m = 1 + x.max()
	t = m - x
	y = m - graph_dilation(E, t)
	return y

def cross_dilation(d):
	E = grid_structuring(*d.shape)
	r = graph_dilation(E, d.flatten())
	return r.reshape(*d.shape)

def cross_erosion(d):
	E = grid_structuring(*d.shape)
	r = graph_erosion(E, d.flatten())
	return r.reshape(*d.shape)

def cross_median(d,p=3):
	from scipy.signal import medfilt
	return medfilt(d,p).squeeze()
	#E = grid_structuring(*d.shape)
	#r = graph_median(E, d.flatten())
	#return r.reshape(*d.shape)

#def fill_nans_by_laplace_equation(xx):
#	from numpy import isfinite
#	from scipy.sparse import eye, spsolve, diags
#
#	# extract data
#	w,h = xx.shape
#	x = x.flatten()
#	m = isnan(x).astype(float)
#
#	# build matrices
#	B = grid_incidence(w, h)
#	L = -B.T @ B            # laplacian operator
#	M = diags(0 + (m > 0))  # R.O.I. mask operator
#	I = eye(w*h)            # identity matrix
#
#	# build the linear system
#	A = M @ L + I - M
#	b = (I - M) @ x
#
#	# solve the symmetrized linear system
#	u = spsolve(A.T @ A, A.T @ b)
#
#	return u.reshape(w, h)
#
#def fill_nans_by_laplace_equation(xx):
#	from numpy import isnan
#	from scipy.sparse import eye, spsolve, diags
#	w = xx.shape[0]                # width of the image domain
#	h = xx.shape[1]                # height of the image domain
#	x = x.flatten()                # flattened data
#	B = grid_incidence(w, h)       # signed incidence matrix
#	L = -B.T @ B                   # laplacian operator
#	M = diags(0.0 + isnan(x))      # mask operator
#	I = eye(w*h)                   # identity matrix of the appropriate size
#	A = M @ L + I - M              # left hand-side matrix of the system Au=b
#	b = (I - M) @ x                # right hand side
#	u = spsolve(A.T @ A, A.T @ b)  # solve the symmetrized system
#	return u.reshape(w, h)


def fill_nans_by_laplace_equation(xx):
	from numpy import isnan, isinf, nan_to_num
	from scipy.sparse import eye, diags
	from scipy.sparse.linalg import spsolve
	s = xx.shape                    # shape of the image domain
	x = xx.flatten()                # flattened data
	x[isinf(x)] = float("nan")      # normalize inf to nan
	L = grid_laplacian(s[0], s[1])  # laplacian matrix
	M = diags(0.0 + isnan(x))       # complementary mask operator
	I = eye(s[0] * s[1])            # identity matrix
	A = M @ L + I - M               # left hand-side matrix of the system
	b = (I - M) @ nan_to_num(x)     # right hand side vector
	u = spsolve(A.T @ A, A.T @ b)   # solve the symmetrized system
	return u.reshape(s[0], s[1])    # reshape and return

def descending_neumann_interpolation(xx):
	from numpy import isnan, isinf, nan_to_num
	from scipy.sparse import eye, diags
	from scipy.sparse.linalg import spsolve
	s = xx.shape                    # shape of the image domain
	x = xx.flatten()                # flattened data
	m = 0.0 + isnan(x)              # mask indicator function
	x[isinf(x)] = float("nan")      # normalize inf to nan
	x = nan_to_num(x)               # sanitize data (fix numpy regression)
	M = diags(m)                    # mask operator
	I = eye(s[0] * s[1])            # identity matrix
	B = grid_incidence(s[0], s[1])  # signed incidence matrix
	L = -B.T @ B                    # laplacian matrix

	A = M @ L + I - M               # left hand-side matrix of the system
	b = (I - M) @ x                 # right hand side vector
	u = spsolve(A.T @ A, A.T @ b)   # solve the symmetrized system

	n = 0 <= (B@m) * (B@u)          # descending edge complementary mask
	N = diags(0.0 + n)              # descending mask operator
	B = N@B                         # new graph without descending edges
	L = -B.T @ B                    # laplacian matrix

	A = M @ L + I - M               # left hand-side matrix of the system
	b = (I - M) @ x                 # right hand side vector
	u = spsolve(A.T @ A, A.T @ b)   # solve the symmetrized system
	return u.reshape(s[0], s[1])    # reshape and return



# API: d = demtk.fuse(D)
def fuse(D):
	"""
	Merge several registered DEM into a single one

	Input  D : stack of DEM (a three-dimensional array)
	Output d : DEM of the same size
	"""
	from numpy import nanmedian
	d = nanmedian(D, axis=0)
	return d

# API: e = demtk.fuse(d)
def fill(d):
	"""
	Fill-in the missing data of a DEM by smooth interpolation

	Input  d : DEM with NANs
	Output e : DEM without NANs
	"""
	e = fill_nans_by_laplace_equation(d)
	return e


# API: x = demtk.register(d, e)
def register(d, e):
	"""
	Find the 3D translation that best registers two DEM

	Input  d : DEM
	Input  e : DEM
	Output x : a 3D vector mapping the domain of e into d
	"""
	from numpy import array
	return array([0,0,0])

# API: e = demtk.shift(d, x)
def shift(d):
	"""
	Apply a 3D translation to the domain of a DEM

	Input  d : a DEM
	Input  x : a 3D vector
	Output e : a DEM with the same data as d, shifted by x
	"""
	e = d + x[2]
	return e

# API: x = demtk.render(d)
def render(d, p=1, s=-1):
	"""
	Render a DEM using simple hillshading and color palette

	Input  d : a DEM
	Output x : a 8-bit color image
	"""
	from numpy import newaxis, sqrt
	d_pal = render_palette_dem(d).astype(float)
	d_lam = qauto(render_shading(d)[:,:,newaxis], p).astype(float)
	d_sao = qauto(render_lao(d, s)[:,:,newaxis], p).astype(float)

	return qauto(sqrt(d_lam*d_sao)*d_pal, p)

	# does not work due to brain-damaged broadcasting rules
	#d_pal = render_palette_dem(d)
	#d_lam = qauto(render_shading(d), 5)
	#d_sao = qauto(render_lao(d), 5)
	#return qauto(sqrt(d_lam*d_sao)*d_pal, 2)

# API: x = demtk.render(d)
def renderclean(d, p=1):
	"""
	Render a DEM using simple hillshading and color palette

	Input  d : a DEM
	Output x : a 8-bit color image
	"""
	from numpy import newaxis, sqrt
	d_pal = render_palette_dem(d).astype(float)
	d_lam = qauto(render_shading(d)[:,:,newaxis], p).astype(float)
	#d_sao = qauto(render_lao(d,p=1)[:,:,newaxis], 5).astype(float)

	return qauto(sqrt(d_lam)*d_pal, p)

# project a 3D point cloud into a DEM
def project(xyz, xmin, xmax, ymin, ymax, r):
	"""
	Project a 3D point cloud into a DEM

	Input xyz : a matrix of size Nx3 representing the point cloud
	Input xmin, xmax, ymin, ymax : the bounding box to project into
	Input r  : the desired resolution

	Output d : a DEM
	"""
	return d

# elevate a DEM into a 3D point cloud
def elevate(d, xmin, xmax, ymin, ymax):
	"""
	Elevate a DEM into a 3D point cloud

	Input d : a DEM
	Input xmin, xmax, ymin, ymax : georeferencing data of the DEM

	Output xyz : a 3D point cloud
	"""
	return xyz

# colorize a DEM using a satellite image
def colorize(d, img, rpc):
	"""
	Obtain an orthophoto from a DEM and a satellite image

	Input d   : a DEM of size WxH
	Input img : a satellite image with D-dimensional pixels
	Input rpc : the projection model of the iamge

	Output x  : an image of size WxH with D-dimensional pixels
	"""
	return x



# directional derivative o a along s
def render_shading(a, s=(1,1)):
	from numpy import pad
	x = pad(a, ((0,0),(0,1)), 'edge')[:,1:] - a     # da / dx
	y = pad(a, ((0,1),(0,0)), 'edge')[1:,:] - a     # da / dy
	z = s[0] * x + s[1] * y                         # da / ds
	return z

# cast shadows from a given sun position
def render_shadows(d, s=(1,1,1)):
	import tempfile, iio, os
	fi = tempfile.NamedTemporaryFile(suffix=".tif").name
	fo = tempfile.NamedTemporaryFile(suffix=".png").name
	iio.write(fi, d)
	os.system(f"shadowcast -M {s[0]} {s[1]} {s[2]} {fi} {fo}")
	z = iio.read(fo).squeeze()
	os.system(f"rm -f {fi} {fo}")
	return z

# Riesz scale-space (periodic boundary conditions)
def filter_riesz0(x, s):
	from numpy.fft import fft2, ifft2, fftfreq
	from numpy import repeat, hypot
	h, w = x.shape
	p = repeat(w*fftfreq(w).reshape(1,w), h, axis=0)  # x-frequencies
	q = repeat(h*fftfreq(h).reshape(h,1), w, axis=1)  # y-frequencies
	r = hypot(p, q)           # image of spectral radius
	r[0,0] = 1                # avoid warnings when 1/0
	X = fft2(x) / r**s        # apply the filter in the frequency domain
	X[0,0] = 0        # for negative s, set the mean to zero
	return ifft2(X).real

# Riesz scale-space (symmetric boundary conditions)
def filter_riesz(x, s):
	from numpy import pad
	h, w = x.shape
	y = pad(x, ((0,h),(0,w)), 'symmetric')
	z = filter_riesz0(y, s)
	return z[0:h,0:w]

# linear (screen space) ambient occlusion
def render_lao(d, s=-1):
	from numpy import fmin
	z = fmin(0, filter_riesz(d, s))
	return z

# quantize a floating-point image into 8 bits
def qauto(x, p=0):
	from numpy import uint8
	if p > 0:
		from numpy import percentile
		m = percentile(x, p)
		M = percentile(x, 100-p)
		X = (((x.astype(float) - m)/(M-m)).clip(0,1)*255.0).astype(uint8)
	else:
		m = x.min()
		M = x.max()
		X  = (255 * (x.astype(float) - m) / (M - m) ).astype(uint8)
	return X

def render_palette_dem(x):
	import iio
	img_terrain = iio.read("i/DEM_poster.png")
	pal_terrain = img_terrain[0][0:256]
	return pal_terrain[qauto(x,p=0.1)]

version = 4

# vim:set tw=80 filetype=python ts=8 sw=8 sts=0 noexpandtab:
