import sys, piio, demtk
x = piio.read(sys.argv[1]).squeeze()
piio.write(sys.argv[2], demtk.render(r))
