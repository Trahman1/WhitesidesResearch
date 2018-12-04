# interpreted by ovitos
# written in Python 3.4

import ovito
import sys

def log(s):
	print(s)


print("Initiating OVITO %i.%i.%i" %ovito.version)

def main():
	script, data, viewfile, renderer, isPrint, outfile = sys.argv
	print([data, viewfile, renderer, isPrint, outfile])
	rs = ovito.vis.RenderSettings(size=(800,600), filename=outfile)
	if renderer == "t":
		rs.renderer = ovito.vis.TachyonRenderer()
	if isPrint == "v":
		print("Is a video!")
		rs.range = ovito.vis.RenderSettings.Range.ANIMATION
		node = ovito.io.import_file(data, multiple_frames = True)
	elif isPrint == "p":
		node = ovito.io.import_file(data)
	
	node.add_to_scene()	
	node.source.cell.display.line_width = 0.

	print("Saving to: " + outfile)
	
	if viewfile == "video1":
		vp = ovito.vis.Viewport()
		vp.type = ovito.vis.Viewport.Type.ORTHO
		vp.camera_pos = (33.0949, 41.088, 18.7402)
		vp.camera_dir = (.756395,-.577565, -.30706)
		vp.fov = 45.8893
		vp.render(rs)

	else:
		ovito.dataset.viewports.active_vp.render(rs)
	

if __name__=='__main__': main()

