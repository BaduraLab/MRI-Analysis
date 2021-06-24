#!/usr/bin/env python

import slicer
import sys

if __name__ == '__main__':
    # Parse arguments
    image_input = sys.argv[1]
    image_output = sys.argv[2]

    [success, loadedVolumeNode] = slicer.util.loadVolume(image_input, returnNode=True)
    slicer.util.saveNode(loadedVolumeNode, image_output)

    exit()