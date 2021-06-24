import gl
import sys
print(sys.version)
print(gl.version())
# gl.resetdefaults()

help(gl)

gl.linewidth(3)

gl.opacity(1, 40)
gl.opacity(2, 0)
gl.colorname(1, 'jet')
gl.colorname(2, 'jet')


# gl.viewcoronal(True)
#gl.view(2)
#gl.orthoviewmm(0.5,0.5,0.5)
#gl.savebmp('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Software/mricrogl/doesthisworknow')