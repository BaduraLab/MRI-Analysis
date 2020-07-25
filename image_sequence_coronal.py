import gl
import sys

print(sys.version)
print(gl.version())
# gl.resetdefaults()

gl.linewidth(1)
# gl.viewcoronal(True)
#help(gl)
gl.view(2)
gl.orthoviewmm(0.5,0.5,0.5)
gl.savebmp('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Software/mricrogl/doesthisworknow')